
#include "DetInterface/IGeoSvc.h"
#include "DetInterface/ITrackingGeoSvc.h"
#include "RecInterface/ITrackSeedingTool.h"

#include "GaudiKernel/SystemOfUnits.h"



#include "datamodel/TrackStateCollection.h"


#include <cmath>
#include <random>

#include "KalmanFilter.h"


DECLARE_ALGORITHM_FACTORY(KalmanFilter)

KalmanFilter::KalmanFilter(const std::string& name, ISvcLocator* svcLoc) : GaudiAlgorithm(name, svcLoc) 
{
  declareProperty("TrackStates", m_trackStates, "tracks/TrackStates");
}

StatusCode KalmanFilter::initialize() {

  return StatusCode::SUCCESS;
}

StatusCode KalmanFilter::execute() {

  fcc::TrackStateCollection* trackStateCollection = new fcc::TrackStateCollection();


  m_trackStates.put(trackStateCollection);
  return StatusCode::SUCCESS;
}

StatusCode KalmanFilter::finalize() {
  StatusCode sc = GaudiAlgorithm::finalize();
  return sc;
}
