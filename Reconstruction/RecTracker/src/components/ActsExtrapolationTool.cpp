#include "ActsExtrapolationTool.h"
#include "ACTS/Detector/TrackingGeometry.hpp"
#include "ACTS/Extrapolation/ExtrapolationCell.hpp"
#include "ACTS/Extrapolation/ExtrapolationEngine.hpp"
#include "ACTS/Extrapolation/IExtrapolationEngine.hpp"
#include "ACTS/Extrapolation/MaterialEffectsEngine.hpp"
#include "ACTS/Extrapolation/RungeKuttaEngine.hpp"
#include "ACTS/Extrapolation/StaticEngine.hpp"
#include "ACTS/Extrapolation/StaticNavigationEngine.hpp"
#include "ACTS/MagneticField/ConstantBField.hpp"
#include "ACTS/Utilities/Logger.hpp"

#include "RecTracker/ACTSLogger.h"

#include "GaudiKernel/SystemOfUnits.h"


DECLARE_TOOL_FACTORY(ActsExtrapolationTool)

ActsExtrapolationTool::ActsExtrapolationTool(const std::string& type, const std::string& name, const IInterface* parent)
    : GaudiTool(type, name, parent),
      m_trkGeoSvc("TrackingGeoSvc", "ExtrapolationTool"),
      m_collectSensitive(true),
      m_collectPassive(false),
      m_collectBoundary(false),
      m_collectMaterial(false),
      m_sensitiveCurvilinear(false),
      m_searchMode(1),
      m_pathLimit(-1) {
  declareInterface<ITrackExtrapolationTool>(this);
  declareProperty("trackingGeometrySvc", m_trkGeoSvc, "The geometry service providing the tracking geometry");
  declareProperty("collectSensitive", m_collectSensitive, "Switch if sensitive hits should be collected");
  declareProperty("collectPassive", m_collectPassive, "Switch if hits on passive material should be collected");
  declareProperty("collectBoundary", m_collectBoundary, "Switch if hits on boundaries should be collected");
  declareProperty("collectMaterial", m_collectMaterial, "Switch if material should be collected along the way");
  declareProperty("sensitiveCurvilinear", m_sensitiveCurvilinear,
                  "Stay with curvilinear parameters for sensitive mode");
  declareProperty("searchMode", m_searchMode, "Depth of search applied");
  declareProperty("pathLimit", m_pathLimit, "The given path limit (-1 if no limit)");
}

StatusCode ActsExtrapolationTool::initialize() {
  StatusCode sc = GaudiTool::initialize();
  if (sc.isFailure()) {
    return sc;
  }
  if (!m_trkGeoSvc) {
    error() << "Did not retrieve tracking geometry service!" << endmsg;
    return StatusCode::FAILURE;
  }

  auto trackingGeometry = m_trkGeoSvc->trackingGeometry();
  if (!trackingGeometry) {
    error() << "Could not retrieve tracking geometry!" << endmsg;
    return StatusCode::FAILURE;
  }

  info() << "Set up the Extrapolator" << endmsg;
  // EXTRAPOLATOR - set up the extrapolator

  /// @todo hardcode engines, possibly create interfaces for future
  // (a) RungeKuttaPropagtator
  using RKEngine = Acts::RungeKuttaEngine<Acts::ConstantBField>;
  RKEngine::Config propConfig;
  /// @todo check units
  propConfig.fieldService = std::make_shared<Acts::ConstantBField>(0., 0., 0.);
  auto propEngine = std::make_shared<RKEngine>(propConfig);
  // (b) MaterialEffectsEngine
  auto matConfig = Acts::MaterialEffectsEngine::Config();
  auto materialEngine = std::make_shared<Acts::MaterialEffectsEngine>(matConfig);
  // (c) StaticNavigationEngine
  auto navConfig = Acts::StaticNavigationEngine::Config();
  navConfig.propagationEngine = propEngine;
  navConfig.materialEffectsEngine = materialEngine;
  navConfig.trackingGeometry = trackingGeometry;
  auto navEngine = std::make_shared<Acts::StaticNavigationEngine>(navConfig);
  // (d) the StaticEngine
  auto statConfig = Acts::StaticEngine::Config();
  statConfig.propagationEngine = propEngine;
  statConfig.navigationEngine = navEngine;
  statConfig.materialEffectsEngine = materialEngine;
  auto statEngine = std::make_shared<Acts::StaticEngine>(statConfig);
  // (e) the material engine
  auto exEngineConfig = Acts::ExtrapolationEngine::Config();
  exEngineConfig.trackingGeometry = trackingGeometry;
  exEngineConfig.propagationEngine = propEngine;
  exEngineConfig.navigationEngine = navEngine;
  exEngineConfig.extrapolationEngines = {statEngine};
  m_extrapolationEngine = std::make_unique<Acts::ExtrapolationEngine>(exEngineConfig);

  return sc;
}

Acts::ExtrapolationCell<Acts::TrackParameters>
ActsExtrapolationTool::extrapolate(fcc::TrackState theTrackState) {
  // create the start parameters
  auto refPoint = theTrackState.referencePoint();
  Acts::Vector3D perigee(refPoint.x, refPoint.y, refPoint.z);
  Acts::PerigeeSurface surface(perigee);
  double d0 = theTrackState.d0();
  double z0 = theTrackState.z0();
  double phi = theTrackState.phi();
  double theta = theTrackState.theta();
  double qop = theTrackState.qOverP();
  // parameters
  Acts::ActsVectorD<5> pars;
  pars << d0, z0, phi, theta, qop;
  // some screen output
  std::unique_ptr<Acts::ActsSymMatrixD<5>> cov = nullptr;
  // create the bound parameters
  Acts::BoundParameters startParameters(std::move(cov), std::move(pars), surface);
  // create the extrapolation cell & configure it
  Acts::ExtrapolationCell<Acts::TrackParameters> ecc(startParameters);
  ecc.addConfigurationMode(Acts::ExtrapolationMode::StopAtBoundary);
  ecc.addConfigurationMode(Acts::ExtrapolationMode::FATRAS);
  ecc.searchMode = m_searchMode;
  // now set the behavioral bits
  debug() << "Extrapolation Modes turned on are: " << endmsg;
  if (m_collectSensitive) {
    debug() << "- collect sensitive" << endmsg;
    ecc.addConfigurationMode(Acts::ExtrapolationMode::CollectSensitive);
  }
  if (m_collectPassive) {
    debug() << "- collect passive" << endmsg;
    ecc.addConfigurationMode(Acts::ExtrapolationMode::CollectPassive);
  }
  if (m_collectBoundary) {
    debug() << "- collect boundary" << endmsg;
    ecc.addConfigurationMode(Acts::ExtrapolationMode::CollectBoundary);
  }
  if (m_collectMaterial) {
    debug() << "- collect material" << endmsg;
    ecc.addConfigurationMode(Acts::ExtrapolationMode::CollectMaterial);
  }
  if (m_sensitiveCurvilinear) {
    debug() << "sensitive curvilinear set to true" << endmsg;
    ecc.sensitiveCurvilinear = true;
  }

  // force a stop in the extrapoaltion mode
  if (m_pathLimit > 0.) {
    debug() << "path limit set to: " << m_pathLimit << endmsg;
    ecc.pathLimit = m_pathLimit;
    ecc.addConfigurationMode(Acts::ExtrapolationMode::StopWithPathLimit);
  }
  debug() << "search mode set to: " << m_searchMode << endmsg;
  // screen output
  debug() << "===> forward extrapolation - collecting information <<===" << endmsg;
  Acts::ExtrapolationCode eCode = m_extrapolationEngine->extrapolate(ecc);
  if (eCode.isFailure()) error() << ("Extrapolation failed.") << endmsg;
  if (eCode.isSuccess()) info() << ("Extrapolation finished successfully") << endmsg;

  return ecc;
}

StatusCode ActsExtrapolationTool::finalize() { return GaudiTool::finalize(); }
