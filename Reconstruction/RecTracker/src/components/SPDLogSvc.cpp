#include "SPDLogSvc.h"

// Gaudi
#include "GaudiKernel/IToolSvc.h"


DECLARE_SERVICE_FACTORY(SPDLogSvc)

SPDLogSvc::SPDLogSvc(const std::string& aName, ISvcLocator* aSL):base_class(aName, aSL) {
}

SPDLogSvc::~SPDLogSvc() {}

StatusCode SPDLogSvc::initialize() {
  // Initialize necessary Gaudi components
  if (Service::initialize().isFailure()) {
    error() << "Unable to initialize Service()" << endmsg;
    return StatusCode::FAILURE;
  }

  // create a spd log sink that forwards msg to gaudi msg service
  auto l_gaudi_sink = std::make_shared<spdlog::sinks::gaudi_sink_st>(info());

  auto fcc_spd_logger = std::make_shared<spdlog::logger>("ttlog", l_gaudi_sink);
  spdlog::drop_all();
  spdlog::register_logger(fcc_spd_logger);
  return StatusCode::SUCCESS;
}



StatusCode SPDLogSvc::finalize() {
  return Service::finalize();
}
