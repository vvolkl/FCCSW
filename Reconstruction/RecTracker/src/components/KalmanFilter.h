#ifndef RECTRACKER_KALMANFILTER_H
#define RECTRACKER_KALMANFILTER_H

// GAUDI
#include "GaudiAlg/GaudiAlgorithm.h"
#include "GaudiKernel/ToolHandle.h"
#include "GaudiKernel/RndmGenerators.h"

// FCCSW
#include "FWCore/DataHandle.h"



namespace fcc {
class TrackStateCollection;
}

class KalmanFilter : public GaudiAlgorithm {
public:
  KalmanFilter(const std::string& name, ISvcLocator* svcLoc);

  ~KalmanFilter() = default;

  StatusCode initialize() override final;

  StatusCode execute() override final;

  StatusCode finalize() override final;

private:

  /// input for the kalman fit
  DataHandle<fcc::TrackStateCollection> m_trackStates{"TrackStates", Gaudi::DataHandle::Writer,
                                                                      this};
};



#endif /* RECTRACKER_KALMANFILTER_H */
