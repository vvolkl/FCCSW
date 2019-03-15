
#include "SimG4Dream/DreamDetectorConstruction.hh"

#include "SimG4DreamDetector.h"

// Geant4
#include "G4VUserDetectorConstruction.hh"

DECLARE_TOOL_FACTORY(SimG4DreamDetector)

SimG4DreamDetector::SimG4DreamDetector(const std::string& aType, const std::string& aName, const IInterface* aParent)
    : GaudiTool(aType, aName, aParent) {
  declareInterface<ISimG4DetectorConstruction>(this);
}

SimG4DreamDetector::~SimG4DreamDetector() {}

StatusCode SimG4DreamDetector::initialize() {
  if (GaudiTool::initialize().isFailure()) {
    return StatusCode::FAILURE;
  }
  return StatusCode::SUCCESS;
}

StatusCode SimG4DreamDetector::finalize() { return GaudiTool::finalize(); }

G4VUserDetectorConstruction* SimG4DreamDetector::detectorConstruction() { 
  return new DreamDetectorConstruction();
}
