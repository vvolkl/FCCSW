#include "DetInterface/IGeoSvc.h"
#include "DetInterface/ITrackingGeoSvc.h"
#include "RecInterface/ITrackSeedingTool.h"

#include "GaudiKernel/IRndmGenSvc.h"
#include "GaudiKernel/SystemOfUnits.h"


#include "datamodel/PositionedTrackHitCollection.h"
#include "datamodel/TrackStateCollection.h"

#include "datamodel/GenVertexCollection.h"
#include "datamodel/MCParticleCollection.h"

#include "ACTS/Detector/TrackingGeometry.hpp"
#include "ACTS/EventData/Measurement.hpp"
#include "ACTS/Extrapolation/ExtrapolationCell.hpp"
#include "ACTS/Extrapolation/ExtrapolationEngine.hpp"
#include "ACTS/Extrapolation/IExtrapolationEngine.hpp"
#include "ACTS/Extrapolation/MaterialEffectsEngine.hpp"
#include "ACTS/Extrapolation/RungeKuttaEngine.hpp"
#include "ACTS/Extrapolation/StaticEngine.hpp"
#include "ACTS/Extrapolation/StaticNavigationEngine.hpp"
#include "ACTS/Fitter/KalmanFitter.hpp"
#include "ACTS/Fitter/KalmanUpdator.hpp"
#include "ACTS/MagneticField/ConstantBField.hpp"
#include "ACTS/Surfaces/PerigeeSurface.hpp"
#include "ACTS/Utilities/Definitions.hpp"
#include "ACTS/Utilities/Identifier.hpp"
#include "ACTS/Utilities/Logger.hpp"



#include <cmath>
#include <random>

#include "ExtrapolationTest.h"

#include "TLorentzVector.h"
#include "TVector3.h"

using namespace Acts;
using DefaultCovMatrix = ActsSymMatrix<ParValue_t, NGlobalPars>;

DECLARE_ALGORITHM_FACTORY(ExtrapolationTest)

ExtrapolationTest::ExtrapolationTest(const std::string& name, ISvcLocator* svcLoc)
    : GaudiAlgorithm(name, svcLoc), m_extrapolationTool(nullptr) {
  declareProperty("extrapolationTool", m_extrapolationTool,
                  "Pointer to extrapolation tool, needed to extrapolate through the tracker.");
  declareProperty("ExtrapolatedTrackStates", m_extrapolatedTrackStates, "ExtrapolatedTrackStates");
  declareProperty("genParticles", m_genParticles, "Handle for the EDM MC particles to be read");
}

StatusCode ExtrapolationTest::initialize() {

  StatusCode sc = GaudiAlgorithm::initialize();
  if (sc.isFailure()) return sc;
  // retrieve the extrapolation tool
  if (!m_extrapolationTool.retrieve()) {
    error() << "Extrapolation tool cannot be retrieved" << endmsg;
    return StatusCode::FAILURE;
  }

  return sc;
}

StatusCode ExtrapolationTest::execute() {

  // get the input mc particles
  const fcc::MCParticleCollection* mcparticles = m_genParticles.get();
  // create the TrackStateCollection to be written out
  auto exTrackStateCollection = m_extrapolatedTrackStates.createAndPut();
  // go through all particles to be extrapolated for this event
  for (const auto& mcparticle : *mcparticles) {

    auto vertex = mcparticle.startVertex();
    auto p4 = mcparticle.core().p4;
    double phi = std::atan2(p4.py, p4.px);
    double p3Mag = std::sqrt(std::pow(p4.px,2) + std::pow(p4.py,2) + std::pow(p4.pz,2));

    double theta = std::acos(p4.pz / p3Mag);
    double qOverP = 1. / p3Mag * mcparticle.charge() * -1;
    double d0 = 0;
    double z0 = 0;

    auto theTrackState = fcc::TrackState(phi, theta, qOverP, d0, z0, vertex.position(), std::array<float, 15ul>());
    debug() << "start extrapolation ..." << endmsg;
    auto stateVector = m_extrapolationTool->extrapolate(theTrackState);
    for (auto t: stateVector) {
      auto theTrackstate = exTrackStateCollection->create();
      theTrackstate.referencePoint(t.referencePoint());

      //debug
      fcc::Point _rp = t.referencePoint();
      std::cout << "Track position, r: \t" << std::sqrt(std::pow(_rp.x,2) + std::pow(_rp.y,2))  << " phi: \t" << std::atan2(_rp.y, _rp.x) << " z: \t" << _rp.z << std::endl;
    }

  }

  debug() << "got " << exTrackStateCollection->size() << " extrapolation steps" << endmsg;

  return StatusCode::SUCCESS;
}

StatusCode ExtrapolationTest::finalize() {
  StatusCode sc = GaudiAlgorithm::finalize();
  return sc;
}
