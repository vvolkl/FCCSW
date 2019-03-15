
//#include "DreamActionInitialization.hh"

#include "SimG4DreamActions.h"


DECLARE_COMPONENT(SimG4DreamActions)

SimG4DreamActions::SimG4DreamActions(const std::string& type, const std::string& name, const IInterface* parent)
    : AlgTool(type, name, parent) {
  declareInterface<ISimG4ActionTool>(this);
}

SimG4DreamActions::~SimG4DreamActions() {}

StatusCode SimG4DreamActions::initialize() {
  if (AlgTool::initialize().isFailure()) {
    return StatusCode::FAILURE;
  }
  return StatusCode::SUCCESS;
}

StatusCode SimG4DreamActions::finalize() { return AlgTool::finalize(); }

G4VUserActionInitialization* SimG4DreamActions::userActionInitialization() {
  //return new DreamActionInitialization();
  return nullptr;
}
