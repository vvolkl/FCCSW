#include "SimG4SaveParticleHistory.h"

// FCCSW
#include "SimG4Common/ParticleInformation.h"
#include "SimG4Common/Units.h"

// Geant4
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4PrimaryParticle.hh"
#include "G4PrimaryVertex.hh"

// datamodel
#include "datamodel/GenVertexCollection.h"
#include "datamodel/MCParticleCollection.h"

DECLARE_COMPONENT(SimG4SaveParticleHistory)

SimG4SaveParticleHistory::SimG4SaveParticleHistory(const std::string& aType, const std::string& aName,
                                                   const IInterface* aParent)
    : GaudiTool(aType, aName, aParent) {
  declareInterface<ISimG4SaveOutputTool>(this);
  declareProperty("mcParticles", m_mcParticles, "Handle to the secondary particles");
  declareProperty("genVertices", m_genVertices, "Handle to the decay vertices");
}


StatusCode SimG4SaveParticleHistory::saveOutput(const G4Event& aEvent) {
  auto evtinfo = dynamic_cast<sim::EventInformation*>(aEvent.GetUserInformation());
  // take over ownership of particle and vertex collections
  auto genVertexColl = evtinfo->releaseVertexCollection();
  auto mcParticleColl = evtinfo->releaseParticleCollection(); 
  info() << "Saved " << mcParticleColl->size() << " particles from Geant4 history." << endmsg;
  m_mcParticles.put(mcParticleColl.release());
  m_genVertices.put(genVertexColl.release());

  return StatusCode::SUCCESS;
}
