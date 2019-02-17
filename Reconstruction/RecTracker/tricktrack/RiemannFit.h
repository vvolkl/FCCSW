#ifndef TRICKTRACK_RIEMANNFIT_H
#define TRICKTRACK_RIEMANNFIT_H

#include <Eigen/Core>
#include <Eigen/Eigenvalues>

using namespace std;

namespace tricktrack {

using namespace Eigen;

constexpr double d = 1.e-4;         //!< used in numerical derivative (J2 in Circle_fit())
constexpr unsigned int max_nop = 8;  //!< In order to avoid use of dynamic memory

using MatrixNd = Eigen::Matrix<double, Dynamic, Dynamic, 0, max_nop, max_nop>;
using ArrayNd = Eigen::Array<double, Dynamic, Dynamic, 0, max_nop, max_nop>;
using Matrix2Nd = Eigen::Matrix<double, Dynamic, Dynamic, 0, 2 * max_nop, 2 * max_nop>;
using Matrix3Nd = Eigen::Matrix<double, Dynamic, Dynamic, 0, 3 * max_nop, 3 * max_nop>;
using Matrix2xNd = Eigen::Matrix<double, 2, Dynamic, 0, 2, max_nop>;
using Array2xNd = Eigen::Array<double, 2, Dynamic, 0, 2, max_nop>;
using Matrix3xNd = Eigen::Matrix<double, 3, Dynamic, 0, 3, max_nop>;
using MatrixNx3d = Eigen::Matrix<double, Dynamic, 3, 0, max_nop, 3>;
using MatrixNx5d = Eigen::Matrix<double, Dynamic, 5, 0, max_nop, 5>;
using VectorNd = Eigen::Matrix<double, Dynamic, 1, 0, max_nop, 1>;
using Vector2Nd = Eigen::Matrix<double, Dynamic, 1, 0, 2 * max_nop, 1>;
using Vector3Nd = Eigen::Matrix<double, Dynamic, 1, 0, 3 * max_nop, 1>;
using RowVectorNd = Eigen::Matrix<double, 1, Dynamic, 1, 1, max_nop>;
using RowVector2Nd = Eigen::Matrix<double, 1, Dynamic, 1, 1, 2 * max_nop>;
using Matrix5d = Eigen::Matrix<double, 5, 5>;
using Matrix6d = Eigen::Matrix<double, 6, 6>;
using Vector5d = Eigen::Matrix<double, 5, 1>;
using u_int    = unsigned int;

struct circle_fit {
  Vector3d par;  //!< parameter: (X0,Y0,R)
  Matrix3d cov;
  /*!< covariance matrix: \n
      |cov(X0,X0)|cov(Y0,X0)|cov( R,X0)| \n
      |cov(X0,Y0)|cov(Y0,Y0)|cov( R,Y0)| \n
      |cov(X0, R)|cov(Y0, R)|cov( R, R)|
  */
  int q;  //!< particle charge
  double chi2;
};

struct line_fit {
  Vector2d par;  //!<(cotan(theta),Zip)
  Matrix2d cov;
  /*!<
      |cov(c_t,c_t)|cov(Zip,c_t)| \n
      |cov(c_t,Zip)|cov(Zip,Zip)|
  */
  double chi2;
};

struct helix_fit {
  Vector5d par;  //!<(phi,Tip,pt,cotan(theta)),Zip)
  Matrix5d cov;
  /*!< ()->cov() \n
      |(phi,phi)|(Tip,phi)|(p_t,phi)|(c_t,phi)|(Zip,phi)| \n
      |(phi,Tip)|(Tip,Tip)|(p_t,Tip)|(c_t,Tip)|(Zip,Tip)| \n
      |(phi,p_t)|(Tip,p_t)|(p_t,p_t)|(c_t,p_t)|(Zip,p_t)| \n
      |(phi,c_t)|(Tip,c_t)|(p_t,c_t)|(c_t,c_t)|(Zip,c_t)| \n
      |(phi,Zip)|(Tip,Zip)|(p_t,Zip)|(c_t,Zip)|(Zip,Zip)|
  */
  int q;  //!< particle charge
  double chi2_circle;
  double chi2_line;
  Vector4d fast_fit;
  VectorXd time;  // TO FIX just for profiling
};

/*!
    \brief raise to square.
*/

inline double sqr(const double& a) { return a * a; }

/*!
    \brief Compute cross product of two 2D vector (assuming z component 0),
    returning z component of the result.

    \param a first 2D vector in the product.
    \param b second 2D vector in the product.

    \return z component of the cross product.
*/

inline double cross2D(const Vector2d& a, const Vector2d& b) {
  return a.x() * b.y() - a.y() * b.x();
}

/*!
    \brief Compute the covariance matrix (in radial coordinates) of points in
    the transverse plane due to multiple Coulomb scattering.

    \param p2D 2D points in the transverse plane.
    \param fast_fit fast_fit Vector4d result of the previous pre-fit
    structured in this form:(X0, Y0, R, Theta)).

    \return scatter_cov_rad errors due to multiple scattering.

    \warning input points must be ordered radially from the detector center
    (from inner layer to outer ones; points on the same layer must ordered too).
    \bug currently works only for points in the barrel.

    \details Only the tangential component is computed (the radial one is
    negligible).

 */
// X in input TO FIX
inline MatrixNd Scatter_cov_rad(const Matrix2xNd& p2D, const Vector4d& fast_fit, VectorNd const & rad) {
  u_int n = p2D.cols();
  double theta = fast_fit(3);
  double p = fast_fit(2) / cos(fast_fit(3));
  double X = 0.04;

  MatrixNd scatter_cov_rad = MatrixXd::Zero(n, n);
  const float cos_f_theta = cos(theta);
  const double sig2 = sqr(0.015 / std::abs(p) * (1 + 0.038 * log(X / std::abs(cos_f_theta)))) * (X / cos_f_theta) ;
  for (u_int k = 0; k < n; ++k) {
    for (u_int l = k; l < n; ++l) {
      for (u_int i = 0; i < std::min(k, l); ++i) {
        scatter_cov_rad(k, l) += (rad(k) - rad(i)) * (rad(l) - rad(i)) * sig2 / sqr(sin(theta));
        scatter_cov_rad(l, k) = scatter_cov_rad(k, l);
      }
    }
  }
  return scatter_cov_rad;
}

/*!
    \brief Transform covariance matrix from radial (only tangential component)
    to Cartesian coordinates (only transverse plane component).

    \param p2D 2D points in the transverse plane.
    \param cov_rad covariance matrix in radial coordinate.

    \return cov_cart covariance matrix in Cartesian coordinates.
*/

inline Matrix2Nd cov_radtocart(const Matrix2xNd& p2D,
                               const MatrixNd& cov_rad,
                               const VectorNd &rad) {
  u_int n = p2D.cols();
  Matrix2Nd cov_cart = MatrixXd::Zero(2 * n, 2 * n);
  VectorNd rad_inv = rad.cwiseInverse();
  for (u_int i = 0; i < n; ++i) {
    for (u_int j = i; j < n; ++j) {
      cov_cart(i, j) = cov_rad(i, j) * p2D(1, i) * rad_inv(i) * p2D(1, j) * rad_inv(j);
      cov_cart(i + n, j + n) = cov_rad(i, j) * p2D(0, i) * rad_inv(i) * p2D(0, j) * rad_inv(j);
      cov_cart(i, j + n) = -cov_rad(i, j) * p2D(1, i) * rad_inv(i) * p2D(0, j) * rad_inv(j);
      cov_cart(i + n, j) = -cov_rad(i, j) * p2D(0, i) * rad_inv(i) * p2D(1, j) * rad_inv(j);

      cov_cart(j, i) = cov_cart(i, j);
      cov_cart(j + n, i + n) = cov_cart(i + n, j + n);
      cov_cart(j + n, i) = cov_cart(i, j + n);
      cov_cart(j, i + n) = cov_cart(i + n, j);
    }
  }
  return cov_cart;
}

/*!
    \brief Transform covariance matrix from Cartesian coordinates (only
    transverse plane component) to radial coordinates (both radial and
    tangential component but only diagonal terms, correlation between different
    point are not managed).

    \param p2D 2D points in transverse plane.
    \param cov_cart covariance matrix in Cartesian coordinates.

    \return cov_rad covariance matrix in raidal coordinate.

    \warning correlation between different point are not computed.
*/
inline MatrixNd cov_carttorad(const Matrix2xNd& p2D,
                       const Matrix2Nd& cov_cart,
                       const VectorNd& rad) {
  u_int n = p2D.cols();
  MatrixNd cov_rad = MatrixXd::Zero(n, n);
  const VectorNd rad_inv2 = rad.cwiseInverse().array().square();
  for (u_int i = 0; i < n; ++i) {
    //!< in case you have (0,0) to avoid dividing by 0 radius
    if (rad(i) < 1.e-4)
      cov_rad(i, i) = cov_cart(i, i);
    else {
      cov_rad(i, i) =
          rad_inv2(i) * (cov_cart(i, i) * sqr(p2D(1, i)) + cov_cart(i + n, i + n) * sqr(p2D(0, i)) -
                         2. * cov_cart(i, i + n) * p2D(0, i) * p2D(1, i));
    }
  }
  return cov_rad;
}

/*!
    \brief Transform covariance matrix from Cartesian coordinates (only
    transverse plane component) to coordinates system orthogonal to the
    pre-fitted circle in each point.
    Further information in attached documentation.

    \param p2D 2D points in transverse plane.
    \param cov_cart covariance matrix in Cartesian coordinates.
    \param fast_fit fast_fit Vector4d result of the previous pre-fit
    structured in this form:(X0, Y0, R, tan(theta))).

    \return cov_rad covariance matrix in the pre-fitted circle's
    orthogonal system.

*/

inline MatrixNd cov_carttorad_prefit(const Matrix2xNd& p2D, const Matrix2Nd& cov_cart,
                              const Vector4d& fast_fit,
                              const VectorNd& rad) {
  u_int n = p2D.cols();
  MatrixNd cov_rad = MatrixXd::Zero(n, n);
  for (u_int i = 0; i < n; ++i) {
    //!< in case you have (0,0) to avoid dividing by 0 radius
    if (rad(i) < 1.e-4)
      cov_rad(i, i) = cov_cart(i, i);  // TO FIX
    else {
      Vector2d a = p2D.col(i);
      Vector2d b = p2D.col(i) - fast_fit.head(2);
      const double x2 = a.dot(b);
      const double y2 = cross2D(a, b);
      const double tan_c = - y2/x2;
      const double tan_c2 = sqr(tan_c);
      cov_rad(i, i) =
          1. / (1. + tan_c2) *
          (cov_cart(i, i) + cov_cart(i + n, i + n) * tan_c2 + 2 * cov_cart(i, i + n) * tan_c);
    }
  }
  return cov_rad;
}

/*!
    \brief Compute the points' weights' vector for the circle fit when multiple
    scattering is managed.
    Further information in attached documentation.

    \param p2D 2D points in transverse plane.
    \param cov_rad_inv covariance matrix inverse in radial coordinated
    (or, beter, pre-fitted circle's orthogonal system).

    \return weight VectorNd points' weights' vector.

    \bug I'm not sure this is the right way to compute the weights for non
    diagonal cov matrix. Further investigation needed.
*/

inline VectorNd Weight_circle(const Matrix2xNd& /*p2D*/, const MatrixNd& cov_rad_inv) {
  return cov_rad_inv.colwise().sum().transpose();
}

/*!
    \brief Compute the points' weights' vector for the line fit (ODR).
    Results from a pre-fit is needed in order to take the orthogonal (to the
    line) component of the errors.

    \param x_err2 squared errors in the x axis.
    \param y_err2 squared errors in the y axis.
    \param tan_theta tangent of theta (angle between y axis and line).

    \return weight points' weights' vector for the line fit (ODR).
*/

inline VectorNd Weight_line(const ArrayNd& x_err2, const ArrayNd& y_err2, const double& tan_theta) {
  return (1. + sqr(tan_theta)) * 1. / (x_err2 + y_err2 * sqr(tan_theta));
}

/*!
    \brief Find particle q considering the  sign of cross product between
    particles velocity (estimated by the first 2 hits) and the vector radius
    between the first hit and the center of the fitted circle.

    \param p2D 2D points in transverse plane.
    \param par_uvr result of the circle fit in this form: (X0,Y0,R).

    \return q int 1 or -1.
*/

inline int Charge(const Matrix2xNd& p2D, const Vector3d& par_uvr) {
  return ((p2D(0, 1) - p2D(0, 0)) * (par_uvr.y() - p2D(1, 0)) -
              (p2D(1, 1) - p2D(1, 0)) * (par_uvr.x() - p2D(0, 0)) >
          0)
             ? -1
             : 1;
}

/*!
    \brief Transform circle parameter from (X0,Y0,R) to (phi,Tip,p_t) and
    consequently covariance matrix.

    \param circle_uvr parameter (X0,Y0,R), covariance matrix to
    be transformed and particle charge.
    \param B magnetic field in Gev/cm/c unit.
    \param error flag for errors computation.
*/

inline void par_uvrtopak(circle_fit& circle, const double& B, const bool& error) {
  Vector3d par_pak;
  const double temp0 = circle.par.head(2).squaredNorm();
  const double temp1 = sqrt(temp0);
  par_pak << atan2(circle.q * circle.par(0), -circle.q * circle.par(1)),
      circle.q * (temp1 - circle.par(2)), circle.par(2) * B;
  if (error) {
    const double temp2 = sqr(circle.par(0)) * 1. / temp0;
    const double temp3 = 1. / temp1 * circle.q;
    Matrix3d J4;
    J4 << -circle.par(1) * temp2 * 1. / sqr(circle.par(0)), temp2 * 1. / circle.par(0), 0.,
        circle.par(0) * temp3, circle.par(1) * temp3, -circle.q, 0., 0., B;
    circle.cov = J4 * circle.cov * J4.transpose();
  }
  circle.par = par_pak;
}

/*!
    \brief Compute the error propagation to obtain the square errors in the
    x axis for the line fit. If errors have not been computed in the circle fit
    than an'approximation is made.
    Further information in attached documentation.

    \param V hits' covariance matrix.
    \param circle result of the previous circle fit (only the covariance matrix
    is needed) TO FIX
    \param J Jacobian of the transformation producing x values.
    \param error flag for error computation.

    \return x_err2 squared errors in the x axis.
*/

inline VectorNd X_err2(const Matrix3Nd& V, const circle_fit& circle, const MatrixNx5d& J,
                const bool& error, u_int n) {
  VectorNd x_err2(n);
  for (u_int i = 0; i < n; ++i) {
    Matrix5d Cov = MatrixXd::Zero(5, 5);
    if (error) Cov.block(0, 0, 3, 3) = circle.cov;
    Cov(3, 3) = V(i, i);
    Cov(4, 4) = V(i + n, i + n);
    Cov(3, 4) = Cov(4, 3) = V(i, i + n);
    x_err2(i) = J.row(i) * Cov * J.row(i).transpose();
  }
  return x_err2;
}

/*!
    \brief Compute the eigenvector associated to the minimum eigenvalue.

    \param A the Matrix you want to know eigenvector and eigenvalue.
    \param chi2 the double were the chi2-related quantity will be stored.

    \return the eigenvector associated to the minimum eigenvalue.

    \warning double precision is needed for a correct assessment of chi2.

    \details The minimus eigenvalue is related to chi2.
    We exploit the fact that the matrix is symmetrical and small (2x2 for line
    fit and 3x3 for circle fit), so the SelfAdjointEigenSolver from Eigen
    library is used, with the computedDirect  method (available only for 2x2
    and 3x3 Matrix) wich computes eigendecomposition of given matrix using a
    fast closed-form algorithm.
    For this optimization the matrix type must be known at compiling time.

*/

inline Vector3d min_eigen3D(const Matrix3d& A, double& chi2) {
  SelfAdjointEigenSolver<Matrix3d> solver(3);
  solver.computeDirect(A);
  int min_index;
  chi2 = solver.eigenvalues().minCoeff(&min_index);
  return solver.eigenvectors().col(min_index);
}

/*!
    \brief A faster version of min_eigen3D() where double precision is not
    needed.

    \param A the Matrix you want to know eigenvector and eigenvalue.
    \param chi2 the double were the chi2-related quantity will be stored

    \return the eigenvector associated to the minimum eigenvalue.

    \detail The computedDirect() method of SelfAdjointEigenSolver for 3x3 Matrix
    indeed, use trigonometry function (it solves a third degree equation) which
    speed up in  single precision.
*/

inline Vector3d min_eigen3D_fast(const Matrix3d& A) {
  SelfAdjointEigenSolver<Matrix3f> solver(3);
  solver.computeDirect(A.cast<float>());
  int min_index;
  solver.eigenvalues().minCoeff(&min_index);
  return solver.eigenvectors().col(min_index).cast<double>();
}

/*!
    \brief 2D version of min_eigen3D().

    \param A the Matrix you want to know eigenvector and eigenvalue.
    \param chi2 the double were the chi2-related quantity will be stored

    \return the eigenvector associated to the minimum eigenvalue.

    \detail The computedDirect() method of SelfAdjointEigenSolver for 2x2 Matrix
    do not use special math function (just sqrt) therefore it doesn't speed up
    significantly in single precision.
*/

inline Vector2d min_eigen2D(const Matrix2d& A, double& chi2) {
  SelfAdjointEigenSolver<Matrix2d> solver(2);
  solver.computeDirect(A);
  int min_index;
  chi2 = solver.eigenvalues().minCoeff(&min_index);
  return solver.eigenvectors().col(min_index);
}

/*!
    \brief A very fast helix fit: it fits a circle by three points (first, middle
    and last point) and a line by two points (first and last).

    \param hits points to be fitted

    \return result in this form: (X0,Y0,R,tan(theta)).

    \warning points must be passed ordered (from internal layer to external) in
    order to maximize accuracy and do not mistake tan(theta) sign.

    \details This fast fit is used as pre-fit which is needed for:
    - weights estimation and chi2 computation in line fit (fundamental);
    - weights estimation and chi2 computation in circle fit (useful);
    - computation of error due to multiple scattering.
*/

inline Vector4d Fast_fit(const Matrix3xNd& hits) {
  Vector4d result;
  u_int n = hits.cols(); // get the number of hits

  // CIRCLE FIT
  // Make segments between middle-to-first(b) and last-to-first(c) hits
  const Vector2d b = hits.block(0, n / 2, 2, 1) - hits.block(0, 0, 2, 1);
  const Vector2d c = hits.block(0, n - 1, 2, 1) - hits.block(0, 0, 2, 1);
  // Compute their lengths
  const double b2 = b.squaredNorm();
  const double c2 = c.squaredNorm();
  double X0;
  double Y0;
  // The algebra has been verified (MR). The usual approach has been followed:
  // * use an orthogonal reference frame passing from the first point.
  // * build the segments (chords)
  // * build orthogonal lines through mid points
  // * make a system and solve for X0 and Y0.
  // * add the initial point
  if (abs(b.x()) > abs(b.y())) {  //!< in case b.x is 0 (2 hits with same x)
    const double k = c.x() / b.x();
    const double div = 2. * (k * b.y() - c.y());
    // if aligned TO FIX
    Y0 = (k * b2 - c2) / div;
    X0 = b2 / (2 * b.x()) - b.y() / b.x() * Y0;
  } else {
    const double k = c.y() / b.y();
    const double div = 2. * (k * b.x() - c.x());
    // if aligned TO FIX
    X0 = (k * b2 - c2) / div;
    Y0 = b2 / (2 * b.y()) - b.x() / b.y() * X0;
  }

  result(0) = X0 + hits(0, 0);
  result(1) = Y0 + hits(1, 0);
  result(2) = sqrt(sqr(X0) + sqr(Y0));

  // LINE FIT
  const Vector2d d = hits.block(0, 0, 2, 1) - result.head(2);
  const Vector2d e = hits.block(0, n - 1, 2, 1) - result.head(2);
  // Compute the arc-length between first and last point: L = R * theta = R *  atan (tan (Theta) )
  const double dr = result(2) * atan2(cross2D(d, e), d.dot(e));
  // Simple difference in Z between last and first hit
  const double dz = hits(2, n - 1) - hits(2, 0);

  result(3) = (dr / dz);

  return result;
}

/*!
    \brief Fit a generic number of 2D points with a circle using Riemann-Chernov
    algorithm. Covariance matrix of fitted parameter is optionally computed.
    Multiple scattering (currently only in barrel layer) is optionally handled.

    \param hits2D 2D points to be fitted.
    \param hits_cov2D covariance matrix of 2D points.
    \param fast_fit pre-fit result in this form: (X0,Y0,R,tan(theta)).
    (tan(theta) is not used).
    \param error flag for error computation.
    \param scattering flag for multiple scattering

    \return circle circle_fit:
    -par parameter of the fitted circle in this form (X0,Y0,R); \n
    -cov covariance matrix of the fitted parameter (not initialized if
    error = false); \n
    -q charge of the particle; \n
    -chi2.

    \warning hits must be passed ordered from inner to outer layer (double hits
    on the same layer must be ordered too) so that multiple scattering is
    treated properly.
    \warning Multiple scattering for barrel is still not tested.
    \warning Multiple scattering for endcap hits is not handled (yet). Do not
    fit endcap hits with scattering = true !

    \bug for small pt (<0.3 Gev/c) chi2 could be slightly underestimated.
    \bug further investigation needed for error propagation with multiple
    scattering.
*/

inline circle_fit Circle_fit(const Matrix2xNd& hits2D, const Matrix2Nd& hits_cov2D,
                      const Vector4d& fast_fit, VectorNd const & rad,
                      const bool& error = true,
                      const bool& scattering = false) {
  // INITIALIZATION
  Matrix2Nd V = hits_cov2D;
  u_int n = hits2D.cols();

  // WEIGHT COMPUTATION
  VectorNd weight;
  MatrixNd G;
  double renorm;
  {
    MatrixNd cov_rad;
    cov_rad = cov_carttorad_prefit(hits2D, V, fast_fit, rad);
    // cov_rad = cov_carttorad(hits2D, V);

    if (scattering) {
      MatrixNd scatter_cov_rad = Scatter_cov_rad(hits2D, fast_fit, rad);
      V += cov_radtocart(hits2D, scatter_cov_rad, rad);
      cov_rad += scatter_cov_rad;
      G = cov_rad.inverse();
      renorm = G.sum();
      G *= 1. / renorm;
      weight = Weight_circle(hits2D, G);
    } else {
      weight = cov_rad.diagonal().cwiseInverse();
      renorm = weight.sum();
      weight *= 1. / renorm;
    }
  }

  // SPACE TRANSFORMATION

  // center
  const Vector2d h_ = hits2D.rowwise().mean();  // centroid
  Matrix3xNd p3D(3, n);
  p3D.block(0, 0, 2, n) = hits2D.colwise() - h_;
  Vector2Nd mc(2 * n);  // centered hits, used in error computation
  mc << p3D.row(0).transpose(), p3D.row(1).transpose();

  // scale
  const double q = mc.squaredNorm();
  const double s = sqrt(n * 1. / q);  // scaling factor
  p3D *= s;

  // project on paraboloid
  p3D.row(2) = p3D.block(0, 0, 2, n).colwise().squaredNorm();

  // COST FUNCTION

  // compute
  Matrix3d A = Matrix3d::Zero();
  const Vector3d r0 = p3D * weight;  // center of gravity
  const Matrix3xNd X = p3D.colwise() - r0;
  if (scattering)
    A = X * G * X.transpose();
  else {
    for (u_int i = 0; i < n; ++i) A += weight(i) * (X.col(i) * X.col(i).transpose());
  }


  // minimize
  double chi2;
  Vector3d v = min_eigen3D(A, chi2);
  v *= (v(2) > 0) ? 1 : -1;  // TO FIX dovrebbe essere N(3)>0
  const double c = -v.transpose() * r0;

  // COMPUTE CIRCLE PARAMETER

  // auxiliary quantities
  const double h = sqrt(1. - sqr(v(2)) - 4. * c * v(2));
  const double v2x2_inv = 1. / (2. * v(2));
  const double s_inv = 1. / s;
  Vector3d par_uvr_;  // used in error propagation
  par_uvr_ << -v(0) * v2x2_inv, -v(1) * v2x2_inv, h * v2x2_inv;

  circle_fit circle;
  circle.par << par_uvr_(0) * s_inv + h_(0), par_uvr_(1) * s_inv + h_(1), par_uvr_(2) * s_inv;
  circle.q = Charge(hits2D, circle.par);
  circle.chi2 = abs(chi2) * renorm * 1. / sqr(2 * v(2) * par_uvr_(2) * s);

  // ERROR PROPAGATION
  if (error) {
    ArrayNd Vcs_[2][2];  // cov matrix of center & scaled points
    {
      const Matrix2Nd Vcs = sqr(s) * V + sqr(sqr(s)) * 1. / (4. * q * n) *
                                             (2. * V.squaredNorm() + 4. * mc.transpose() * V * mc) *
                                             mc * mc.transpose();
      Vcs_[0][0] = Vcs.block(0, 0, n, n);
      Vcs_[0][1] = Vcs.block(0, n, n, n);
      Vcs_[1][1] = Vcs.block(n, n, n, n);
      Vcs_[1][0] = Vcs_[0][1].transpose();
    }

    MatrixNd C[3][3];  // cov matrix of 3D transformed points
    {
      const ArrayNd t0 = (VectorXd::Constant(n, 1.) * p3D.row(0));
      const ArrayNd t1 = (VectorXd::Constant(n, 1.) * p3D.row(1));
      const ArrayNd t00 = p3D.row(0).transpose() * p3D.row(0);
      const ArrayNd t01 = p3D.row(0).transpose() * p3D.row(1);
      const ArrayNd t11 = p3D.row(1).transpose() * p3D.row(1);
      const ArrayNd t10 = t01.transpose();
      C[0][0] = Vcs_[0][0];
      C[0][1] = Vcs_[0][1];
      C[0][2] = 2. * (Vcs_[0][0] * t0 + Vcs_[0][1] * t1);
      C[1][1] = Vcs_[1][1];
      C[1][2] = 2. * (Vcs_[1][0] * t0 + Vcs_[1][1] * t1);
      C[2][2] = 2. * (Vcs_[0][0] * Vcs_[0][0] + Vcs_[0][0] * Vcs_[0][1] + Vcs_[1][1] * Vcs_[1][0] +
                      Vcs_[1][1] * Vcs_[1][1]) +
                4. * (Vcs_[0][0] * t00 + Vcs_[0][1] * t01 + Vcs_[1][0] * t10 + Vcs_[1][1] * t11);
    }

    Matrix3d C0;  // cov matrix of center of gravity (r0.x,r0.y,r0.z)
    for (u_int i = 0; i < 3; ++i) {
      for (u_int j = i; j < 3; ++j) {
        C0(i, j) = weight.transpose() * C[i][j] * weight;
        C0(j, i) = C0(i, j);
      }
    }

    const MatrixNd W = weight * weight.transpose();
    const MatrixNd H = MatrixXd::Identity(n, n).rowwise() - weight.transpose();
    const MatrixNx3d s_v = H * p3D.transpose();

    MatrixNd D_[3][3];  // cov(s_v)
    {
      D_[0][0] = (H * C[0][0] * H.transpose()).cwiseProduct(W);
      D_[0][1] = (H * C[0][1] * H.transpose()).cwiseProduct(W);
      D_[0][2] = (H * C[0][2] * H.transpose()).cwiseProduct(W);
      D_[1][1] = (H * C[1][1] * H.transpose()).cwiseProduct(W);
      D_[1][2] = (H * C[1][2] * H.transpose()).cwiseProduct(W);
      D_[2][2] = (H * C[2][2] * H.transpose()).cwiseProduct(W);
      D_[1][0] = D_[0][1].transpose();
      D_[2][0] = D_[0][2].transpose();
      D_[2][1] = D_[1][2].transpose();
    }

    constexpr u_int nu[6][2] = {{0, 0}, {0, 1}, {0, 2}, {1, 1}, {1, 2}, {2, 2}};

    Matrix6d E;  // cov matrix of the 6 independent elements of A
    for (u_int a = 0; a < 6; ++a) {
      const u_int i = nu[a][0], j = nu[a][1];
      for (u_int b = a; b < 6; ++b) {
        const u_int k = nu[b][0], l = nu[b][1];
        VectorNd t0(n);
        VectorNd t1(n);
        if (l == k) {
          t0 =  (D_[j][l] * s_v.col(l)) * 2.;
          if (i == j)
            t1 = t0;
          else
            t1 = (D_[i][l] * s_v.col(l)) * 2.;
        } else {
          t0 = D_[j][l] * s_v.col(k) + D_[j][k] * s_v.col(l);
          if (i == j)
            t1 = t0;
          else
            t1 = D_[i][l] * s_v.col(k) + D_[i][k] * s_v.col(l);
        }

        if (i == j)
          E(a, b) = 0. + s_v.col(i).transpose() * (t0 + t1);
        else
          E(a, b) = 0. + (s_v.col(i).transpose() * t0) + (s_v.col(j).transpose() * t1);
        if (b != a) E(b, a) = E(a, b);
      }
    }

    Matrix<double, 3, 6> J2;  // Jacobian of min_eigen() (numerically computed)
    for (u_int a = 0; a < 6; ++a) {
      const u_int i = nu[a][0], j = nu[a][1];
      Matrix3d Delta = Matrix3d::Zero();
      Delta(i, j) = Delta(j, i) = abs(A(i, j) * d);
      J2.col(a) = min_eigen3D_fast(A + Delta);
      const int sign = (J2.col(a)(2) > 0) ? 1 : -1;
      J2.col(a) = (J2.col(a) * sign - v) / Delta(i, j);
    }

    Matrix4d Cvc;  // joint cov matrix of (v0,v1,v2,c)
    {
      Matrix3d t0 = J2 * E * J2.transpose();
      Vector3d t1 = -t0 * r0;
      Cvc.block(0, 0, 3, 3) = t0;
      Cvc.block(0, 3, 3, 1) = t1;
      Cvc.block(3, 0, 1, 3) = t1.transpose();
      Cvc(3, 3) =
          (v.transpose() * C0 * v) + (C0.cwiseProduct(t0)).sum() + (r0.transpose() * t0 * r0);
    }

    Matrix<double, 3, 4> J3;  // Jacobian (v0,v1,v2,c)->(X0,Y0,R)
    {
      const double t = 1. / h;
      J3 << -v2x2_inv, 0, v(0) * sqr(v2x2_inv) * 2., 0, 0, -v2x2_inv, v(1) * sqr(v2x2_inv) * 2., 0,
          0, 0, -h * sqr(v2x2_inv) * 2. - (2. * c + v(2)) * v2x2_inv * t, -t;
    }

    /// @todo (vvolkl): investigate calculation of errors here and resulting compilation error
    const RowVector2Nd Jq = mc.transpose() * s * 1. / n;  // var(q)
    double JqProd = Jq * V * Jq.transpose();


    
    Matrix3d cov_uvr = J3 * Cvc * J3.transpose() * sqr(s_inv)  // cov(X0,Y0,R)
                       + (par_uvr_ * par_uvr_.transpose()) * JqProd;

    circle.cov = cov_uvr;
  }


  return circle;
}

/*!
    \brief Fit of helix parameter cotan(theta)) and Zip by projection on the
    pre-fitted cylinder  and line fit on its surface.

    \param hits hits coordinates.
    \param hits_cov covariance matrix of the hits.
    \param circle cylinder parameter, their covariance (if computed, otherwise
    uninitialized) and particle charge.
    \param fast_fit result of the previous fast fit in this form:
    (X0,Y0,R,cotan(theta))).
    \param error flag for error computation.

    \return line line_fit:
    -par parameter of the line in this form: (cotan(theta)), Zip); \n
    -cov covariance matrix of the fitted parameter; \n
    -chi2.

    \warning correlation between R and z are neglected, this could be relevant
    if geometry detector provides sloped modules in the R/z plane.

    \bug chi2 and errors could be slightly underestimated for small eta (<0.2)
    when pt is small (<0.3 Gev/c).

    \todo multiple scattering treatment.

    \details Line fit is made by orthogonal distance regression where
    correlation between coordinates in the transverse plane (x,y) and z are
    neglected (for a barrel + endcap geometry this is a very good
    approximation).
    Covariance matrix of the fitted parameter is optionally computed.
    Multiple scattering is not handled (yet).
    A fast pre-fit is performed in order to evaluate weights and to compute
    errors.
*/

inline line_fit Line_fit(const Matrix3xNd& hits, const Matrix3Nd& hits_cov, const circle_fit& circle,
                  const Vector4d& fast_fit, const bool& error = true) {
  u_int n = hits.cols();
  // PROJECTION ON THE CILINDER
  Matrix2xNd p2D(2, n);
  MatrixNx5d Jx(n, 5);

  // x & associated Jacobian
  // cfr https://indico.cern.ch/event/663159/contributions/2707659/attachments/1517175/2368189/Riemann_fit.pdf
  // Slide 11
  // a ==> -o i.e. the origin of the circle in XY plane, negative
  // b ==> p i.e. distances of the points wrt the origin of the circle.
  const Vector2d o(circle.par(0), circle.par(1));
  for (u_int i = 0; i < n; ++i) {  // x
    Vector2d p = hits.block(0, i, 2, 1) - o;
    const double cross = cross2D(-o, p);
    const double dot = (-o).dot(p);
    // atan2(cross, dot) give back the angle in the transverse plane so tha the final equation reads:
    // x_i = -q*R*theta (theta = angle returned by atan2)
    const double atan2_ = -circle.q * atan2(cross, dot);
    p2D(0, i) = atan2_ * circle.par(2);

    // associated Jacobian, used in weights and errors computation
    const double temp0 = -circle.q * circle.par(2) * 1. / (sqr(dot) + sqr(cross));
    double d_X0 = 0, d_Y0 = 0, d_R = 0.;  // good approximation for big pt and eta
    if (error) {
      d_X0 = - temp0 * ((p(1) + o(1)) * dot - (p(0) - o(0)) * cross);
      d_Y0 = temp0 * ((p(0) + o(0)) * dot - (o(1) - p(1)) * cross);
      d_R = atan2_;
    }
    const double d_x = temp0 * (o(1) * dot + o(0) * cross);
    const double d_y = temp0 * (-o(0) * dot + o(1) * cross);
    Jx.row(i) << d_X0, d_Y0, d_R, d_x, d_y;
  }
  // Math of d_{X0,Y0,R,x,y} all verified by hand

  // y
  p2D.row(1) = hits.row(2);


  // WEIGHT COMPUTATION
  VectorNd x_err2 = X_err2(hits_cov, circle, Jx, error, n);
  VectorNd y_err2 = hits_cov.block(2 * n, 2 * n, n, n).diagonal();

  const VectorNd err2_inv = Weight_line(x_err2, y_err2, fast_fit(3));
  const VectorNd weight = err2_inv * 1. / err2_inv.sum();
  // COST FUNCTION

  // compute
  // r0 represents the weighted mean of "x" and "y".
  const Vector2d r0 = p2D * weight;
  // This is the X  vector that will be used to build the
  // scatter matrix S = X^T * X
  const Matrix2xNd X = p2D.colwise() - r0;
  Matrix2d A = Matrix2d::Zero();
  for (u_int i = 0; i < n; ++i) {
    A += err2_inv(i) * (X.col(i) * X.col(i).transpose());
  }
  // minimize
  double chi2;
  Vector2d v = min_eigen2D(A, chi2);
  // n *= (chi2>0) ? 1 : -1; //TO FIX
  const double c = -v.transpose() * r0;

  // COMPUTE LINE PARAMETER
  line_fit line;
  line.par << -v(0) / v(1),                          // cotan(theta))
      -c * sqrt(sqr(v(0)) + sqr(v(1))) * 1. / v(1);  // Zip
  line.chi2 = abs(chi2);

  // ERROR PROPAGATION
  if (error) {
    const double v0_2 = sqr(v(0));
    const double v1_2 = sqr(v(1));

    Matrix3d C;  // cov(v,c)
    {
      double norm_chernov = 0.;
      for (u_int i = 0; i < n; ++i)
        norm_chernov += err2_inv(i) * (v(0) * p2D(0, i) + v(1) * p2D(1, i) + c)
          * (v(0) * p2D(0, i) + v(1) * p2D(1, i) + c);
      norm_chernov /= float(n);
      // Indeed it should read:
      // * compute the average error in the orthogonal direction: err2_inv.cwiseInverse().sum()/sqr(n)
      // * normalize the A(0,0)+A(1,1) dividing by err2_inv.sum(), since those have been weighted
      const double norm = (err2_inv.cwiseInverse().sum())*err2_inv.sum()*1./sqr(n);
      const double sig2 = 1./(A(0,0) + A(1,1))*norm;
//      const double sig2 = 1. / (A(0, 0) + A(1, 1));
      C(0, 0) = sig2 * v1_2;
      C(1, 1) = sig2 * v0_2;
      C(0, 1) = C(1, 0) = -sig2 * v(0) * v(1);
      const VectorNd weight_2 = (weight).array().square();
      const Vector2d C0(weight_2.dot(x_err2), weight_2.dot(y_err2));
      C.block(0, 2, 2, 1) = C.block(2, 0, 1, 2).transpose() = -C.block(0, 0, 2, 2) * r0;
      C(2, 2) = v0_2 * C0(0) + v1_2 * C0(1) + C0(0) * C(0, 0) + C0(1) * C(1, 1) +
                (r0.transpose() * C.block(0, 0, 2, 2) * r0);
    }

    Matrix<double, 2, 3> J;  // Jacobian of (v,c) -> (cotan(theta)),Zip)
    {
      const double t0 = 1. / v(1);
      const double t1 = sqr(t0);
      const double sqrt_ = sqrt(v1_2 + v0_2);
      const double t2 = 1. / sqrt_;
      J << -t0, v(0) * t1, 0, -c * v(0) * t0 * t2, v0_2 * c * t1 * t2, -sqrt_ * t0;
    }

    line.cov = J * C * J.transpose();
  }

  return line;
}

/*!
    \brief Helix fit by three step:
    -fast pre-fit (see Fast_fit() for further info); \n
    -circle fit of hits projected in the transverse plane by Riemann-Chernov
        algorithm (see Circle_fit() for further info); \n
    -line fit of hits projected on cylinder surface by orthogonal distance
        regression (see Line_fit for further info). \n
    Points must be passed ordered (from inner to outer layer).

    \param hits Matrix3xNd hits coordinates in this form: \n
        |x0|x1|x2|...|xn| \n
        |y0|y1|y2|...|yn| \n
        |z0|z1|z2|...|zn|

    \param hits_cov Matrix3Nd covariance matrix in this form (()->cov()): \n

   |(x0,x0)|(x1,x0)|(x2,x0)|.|(y0,x0)|(y1,x0)|(y2,x0)|.|(z0,x0)|(z1,x0)|(z2,x0)| \n
   |(x0,x1)|(x1,x1)|(x2,x1)|.|(y0,x1)|(y1,x1)|(y2,x1)|.|(z0,x1)|(z1,x1)|(z2,x1)| \n
   |(x0,x2)|(x1,x2)|(x2,x2)|.|(y0,x2)|(y1,x2)|(y2,x2)|.|(z0,x2)|(z1,x2)|(z2,x2)| \n
       .       .       .    .    .       .       .    .    .       .       .     \n
   |(x0,y0)|(x1,y0)|(x2,y0)|.|(y0,y0)|(y1,y0)|(y2,x0)|.|(z0,y0)|(z1,y0)|(z2,y0)| \n
   |(x0,y1)|(x1,y1)|(x2,y1)|.|(y0,y1)|(y1,y1)|(y2,x1)|.|(z0,y1)|(z1,y1)|(z2,y1)| \n
   |(x0,y2)|(x1,y2)|(x2,y2)|.|(y0,y2)|(y1,y2)|(y2,x2)|.|(z0,y2)|(z1,y2)|(z2,y2)| \n
       .       .       .    .    .       .       .    .    .       .       .     \n
   |(x0,z0)|(x1,z0)|(x2,z0)|.|(y0,z0)|(y1,z0)|(y2,z0)|.|(z0,z0)|(z1,z0)|(z2,z0)| \n
   |(x0,z1)|(x1,z1)|(x2,z1)|.|(y0,z1)|(y1,z1)|(y2,z1)|.|(z0,z1)|(z1,z1)|(z2,z1)| \n
   |(x0,z2)|(x1,z2)|(x2,z2)|.|(y0,z2)|(y1,z2)|(y2,z2)|.|(z0,z2)|(z1,z2)|(z2,z2)|

   \param B magnetic field in the center of the detector in Gev/cm/c
   unit, in order to perform pt calculation.
   \param error flag for error computation.
   \param scattering flag for multiple scattering treatment.
   (see Circle_fit() documentation for further info).

   \warning see Circle_fit(), Line_fit() and Fast_fit() warnings.

   \bug see Circle_fit(), Line_fit() and Fast_fit() bugs.
*/

inline helix_fit Helix_fit(const Matrix3xNd& hits, const Matrix3Nd& hits_cov, const double& B,
                    const bool& error = true, const bool& scattering = false) {
  u_int n = hits.cols();
  VectorNd rad = (hits.block(0, 0, 2, n).colwise().norm());

  // Fast_fit gives back (X0, Y0, R, theta) w/o errors, using only 3 points.
  const Vector4d fast_fit = Fast_fit(hits);

  circle_fit circle = Circle_fit(hits.block(0, 0, 2, n), hits_cov.block(0, 0, 2 * n, 2 * n),
                                 fast_fit, rad, error, scattering);

  const line_fit line = Line_fit(hits, hits_cov, circle, fast_fit, error);

  par_uvrtopak(circle, B, error);

  helix_fit helix;
  helix.par << circle.par, line.par;
  if (error) {
    helix.cov = MatrixXd::Zero(5, 5);
    helix.cov.block(0, 0, 3, 3) = circle.cov;
    helix.cov.block(3, 3, 2, 2) = line.cov;
  }
  helix.q = circle.q;
  helix.chi2_circle = circle.chi2;
  helix.chi2_line = line.chi2;

  return helix;
}





}
#endif /* TRICKTRACK_RIEMANNFIT_H */
