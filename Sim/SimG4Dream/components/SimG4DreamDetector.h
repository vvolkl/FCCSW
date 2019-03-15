
#ifndef SIMG4_DREAMDETECTOR
#define SIMG4_DREAMDETECTOR

// Gaudi
#include "GaudiAlg/GaudiTool.h"

// FCCSW
#include "SimG4Interface/ISimG4DetectorConstruction.h"

class SimG4DreamDetector : public GaudiTool, virtual public ISimG4DetectorConstruction {
public:
  explicit SimG4DreamDetector(const std::string& aType, const std::string& aName, const IInterface* aParent);
  virtual ~SimG4DreamDetector();
  /**  Initialize.
   *   @return status code
   */
  virtual StatusCode initialize();
  /**  Finalize.
   *   @return status code
   */
  virtual StatusCode finalize();
  /** Get the initilization of the geometry.
   *  @return pointer to G4VUserDetectorConstruction (ownership is transferred to the caller)
   */
  virtual G4VUserDetectorConstruction* detectorConstruction();
};

#endif /* SIMG4_DREAMDETECTOR */
