
#include "DetInterface/IGeoSvc.h"
#include "DetInterface/ITrackingGeoSvc.h"
#include "RecInterface/ITrackSeedingTool.h"

#include "GaudiKernel/SystemOfUnits.h"



#include "datamodel/TrackStateCollection.h"
#include "datamodel/PositionedTrackHitCollection.h"


#include <cmath>
#include <random>

#include <Eigen/Dense>

#include "KalmanFilter.h"


DECLARE_ALGORITHM_FACTORY(KalmanFilter)

KalmanFilter::KalmanFilter(const std::string& name, ISvcLocator* svcLoc) : GaudiAlgorithm(name, svcLoc) 
{
  declareProperty("FittedTracks", m_fittedTracks, "tracks/FittedTracks");
  declareProperty("TrackSeeds", m_trackSeeds, "tracks/TrackSeeds");
  declareProperty("TrackerHits", m_hits, "TrackerHits");
}

StatusCode KalmanFilter::initialize() {

  return StatusCode::SUCCESS;
}

StatusCode KalmanFilter::execute() {


    /* TODO: Input, initial state */
    const fcc::TrackStateCollection* trackSeeds = m_trackSeeds.get();
    for (const auto& trackSeed: *trackSeeds) {
      std::cout << "trackSeed:\t" << trackSeed.d0() << "\t" << trackSeed.z0() << "\t" << trackSeed.theta() << "\t" << trackSeed.phi() << "\t" << trackSeed.qOverP() << std::endl; 
      std::cout << "\t refPoint: " << trackSeed.referencePoint().x << std::endl;
    }


    const fcc::PositionedTrackHitCollection* trackHits = m_hits.get();
    for (const auto& hit: *trackHits) {
      std::cout << "measurement: " <<  hit.position().x << std::endl;
 
    }
    // output
    auto fittedTrackCollection = m_fittedTracks.createAndPut();




    /// p .... track state
    Eigen::Matrix<double, 5, 1> p;
    p << 0.1, 0.2, 0.3, 0.4, 0.5;
    /// C ... track state covariance 
    Eigen::Matrix<double, 5, 5> C = Eigen::Matrix<double, 5, 5>::Random();



    std::list<std::pair<Eigen::Matrix<double, 5, 1>, Eigen::Matrix<double, 5, 5>>> filtered_states;
    std::list<std::pair<Eigen::Matrix<double, 5, 1>, Eigen::Matrix<double, 5, 5>>> smoothed_states;
    std::list<std::pair<Eigen::Matrix<double, 5, 1>, Eigen::Matrix<double, 5, 5>>> predicted_states;

    // dummy
     auto predict = [](Eigen::Matrix<double, 5, 1> theState) { return theState; };
    
    // A) Forward Filter

    // outsourcing the selection of new measurements 
    // the idea is that the next measurement can be selected in a number of ways
    // so the fitter needs an interface that returns measurements as long as there are more steps to take
    // also effectively decouples the fitter logic from the extrapolation and geometry 


    // prepare some dummy measurements and covariances 
    std::vector<Eigen::Matrix<double, 2, 1>> measurement_candidates;
    std::vector<Eigen::Matrix<double, 2, 2>> measurement_candidates_cov;
    for (int i = 0; i < 10; ++i) {
      measurement_candidates.push_back(Eigen::Matrix<double, 2 , 1>::Random());
      measurement_candidates_cov.push_back(Eigen::Matrix<double, 2 , 2>::Random());
    }


    // dummy next_measurement function
    int currentStep = 0;
    auto get_next_measurement = [&measurement_candidates, &measurement_candidates_cov, &currentStep]() {
      if (currentStep < measurement_candidates.size()) {
         auto _m = measurement_candidates[currentStep];
         auto _m_cov = measurement_candidates_cov[currentStep];
         currentStep++;
         return std::make_pair(_m, _m_cov);
      } else {
        Eigen::Matrix<double, 2, 1> m;
        m << std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN();
        return std::make_pair(m, Eigen::Matrix<double, 2, 2>());
      }
      };


    // 2) Find next measurement
    //   * this could be from a given pattern recognition
    //   * or be used iteratively with a number of candidates

    // 3) Gain Matrix Update

    // get matrix H (projector)
    // in the default fcc hh usecase, H just selects the two local coordinates ( eloc0 and eloc1 )
    // H *  covariance * H.transpose()  selects just the submatrix of the  

    Eigen::Matrix<double, 2, 5> H;

    H.row(0) << 1,0,0,0,0;
    H.row(1) << 0,1,0,0,0;
    std::cout << "H" << std::endl;
    std::cout << H << std::endl;

    // helper function
    auto  is_nan = [](const Eigen::Matrix<double, 2, 1>& x)
    {
      return !((x.array() == x.array())).all();
    };



    std::pair<Eigen::Matrix<double, 2, 1>, Eigen::Matrix<double, 2, 2>> meas_pair;
    while (1) {

        auto predicted_state = predict(p);
        predicted_states.push_back(std::make_pair(p, C));


        meas_pair = get_next_measurement();
        if ( is_nan(meas_pair.first)) {
          break;
        }
        std::cout << "next measurement: " << std::endl;
        std::cout << meas_pair.first << std::endl;
        
        /// m ... measurement in local coordinates 
        auto m = meas_pair.first;
        /// V ... measurement covariance
        auto V = meas_pair.second;



        /// K ... Gain matrix
        Eigen::Matrix<double, 5, 2> K = C * H.transpose() * (V + H * C * H.transpose()).inverse();

        static const Eigen::Matrix<double, 5, 5>  unit
                         = Eigen::Matrix<double, 5, 5>::Identity();

        Eigen::Matrix<double, 5, 1> new_p
                  = p + K * (m - H * p);

        std::cout << "updated parameters" << std::endl;

        std::cout << new_p << std::endl;
        Eigen::Matrix<double, 5, 5> new_C = (unit - K * H) * C;

        std::cout << "updated covariance" << std::endl;

        std::cout << new_C << std::endl;


        filtered_states.push_back(std::make_pair(new_p, new_C));



      }





  std::cout << "number of filtered states: " << filtered_states.size() << std::endl;

  //B Smoothing
  auto it = filtered_states.rbegin();

     // for the last measurement the filtered state and the smoothed state are
         // equal


  decltype(it) pLast = it++;
  unsigned int backCounter = 0;
  for (; it != filtered_states.rend(); ++it, ++pLast) {

    /// F ... Transport Jacobian
    Eigen::Matrix<double, 5, 5> F =  Eigen::Matrix<double, 5, 5>::Random();

    std::cout << "filtered states" << std::endl;
    std::cout << (*it).first  << std::endl;
    std::cout << "filtered cov" << std::endl;
    std::cout << (*it).second  << std::endl;
    auto filtered_C = (*it).second;
    /// A: todo
    auto A = filtered_C * F * filtered_C.inverse();

  backCounter ++;
  }








  return StatusCode::SUCCESS;
}

StatusCode KalmanFilter::finalize() {
  StatusCode sc = GaudiAlgorithm::finalize();
  return sc;
}
