
#include "DetInterface/IGeoSvc.h"

#include "GaudiKernel/PhysicalConstants.h"

#include "datamodel/PositionedTrackHitCollection.h"
#include "datamodel/TrackCollection.h"
#include "datamodel/TrackHitCollection.h"
#include "datamodel/TrackStateCollection.h"

#include "DD4hep/Detector.h"
#include "DD4hep/Volumes.h"
#include "DDRec/API/IDDecoder.h"
#include "DDSegmentation/BitField64.h"
#include "DDSegmentation/CartesianGridXZ.h"

#include <cmath>
#include <random>

#include "RecInterface/ITrackSeedingTool.h"
#include "RecInterface/ITrackFittingTool.h"
#include "RecTrackAlg.h"
#include "RecTracker/TrackingUtils.h"



DECLARE_ALGORITHM_FACTORY(RecTrackAlg)

RecTrackAlg::RecTrackAlg(const std::string& name, ISvcLocator* svcLoc) : GaudiAlgorithm(name, svcLoc) {

  declareProperty("TrackerPositionedHits", m_positionedTrackHits, "Tracker hits (Input)");
  declareProperty("Tracks", m_tracks, "Tracks (Output)");
  declareProperty("TrackStates", m_trackStates, "TrackStates (Output)");
  declareProperty("TrackSeedingTool", m_trackSeedingTool);
  declareProperty("TrackFittingTool", m_trackFittingTool);
}

StatusCode RecTrackAlg::initialize() {
  debug() << "initialize" << endmsg;
  return StatusCode::SUCCESS;
}

StatusCode RecTrackAlg::execute() {

  // get hits from event store
  const fcc::PositionedTrackHitCollection* hits = m_positionedTrackHits.get();
  debug() << "hit collection size: " << hits->size() << endmsg;
  // delegate track seeding to tool
  auto seedmap = m_trackSeedingTool->findSeeds(hits);
  // delegate track fitting to tool
  auto tracksAndTrackstates = m_trackFittingTool->fitTracks(hits, seedmap);

  m_tracks.put(tracksAndTrackstates.first);
  m_trackStates.put(tracksAndTrackstates.second);

  return StatusCode::SUCCESS;
}

StatusCode RecTrackAlg::finalize() {
  StatusCode sc = GaudiAlgorithm::finalize();
  return sc;
}
