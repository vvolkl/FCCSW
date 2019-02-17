#ifndef RECINTERFACE_IKALMANMEASUREMENTTOOL_H
#define RECINTERFACE_IKALMANMEASUREMENTTOOL_H

// Gaudi
#include "GaudiKernel/IAlgTool.h"

#include "datamodel/TrackStateCollection.h"

#include <Eigen/Core>



class IKalmanMeasurementTool : virtual public IAlgTool {
public:
  DeclareInterfaceID(IKalmanMeasurementTool, 1, 0);

  virtual std::pair<Eigen::Matrix<double, 2, 1>, Eigen::Matrix<double, 2, 2>> getMeasurement(const fcc::TrackState* theState) = 0;
};

#endif /* RECINTERFACE_IKALMANMEASUREMENTTOOL_H */
