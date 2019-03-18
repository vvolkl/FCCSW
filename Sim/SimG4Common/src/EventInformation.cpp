#include "SimG4Common/EventInformation.h"

#include "G4Track.hh"
#include "G4LorentzVector.hh"

#include "SimG4Common/Units.h"

#include "datamodel/GenVertexCollection.h"
#include "datamodel/MCParticleCollection.h"

namespace sim {
EventInformation::EventInformation() { 
  m_particleData = new std::vector<fcc::G4ParticleData>();
}

void EventInformation::setCollections(std::vector<fcc::G4ParticleData>*& aParticleDataVector ) {
  // ownership is transferred here - to SaveTool which is supposed to put it in the event store
  aParticleDataVector = m_particleData;
}

void EventInformation::addParticle(const G4Track* aSecondary, G4LorentzVector initialPos,
G4LorentzVector initialEnergy) {
  fcc::G4ParticleData p;
  p.motherId = aSecondary->GetParentID();
  p.trackId = aSecondary->GetTrackID();
  p.pdgId = aSecondary->GetDynamicParticle()->GetDefinition()->GetPDGEncoding();

  auto g4mom = aSecondary->GetMomentum();
  auto g4energy = aSecondary->GetTotalEnergy();
  float mass = g4energy * g4energy - g4mom.mag2();
  mass = sqrt(fabs(mass));

  p.p4Final.px = g4mom.x() * sim::g42edm::energy;
  p.p4Final.py = g4mom.y() * sim::g42edm::energy;
  p.p4Final.pz = g4mom.z() * sim::g42edm::energy;
  p.p4Final.mass = mass * sim::g42edm::energy;

  p.p4Initial.px = initialEnergy.px() * sim::g42edm::energy;
  p.p4Initial.py = initialEnergy.py() * sim::g42edm::energy;
  p.p4Initial.pz = initialEnergy.pz() * sim::g42edm::energy;
  p.p4Initial.mass = initialEnergy.m() * sim::g42edm::energy;



  auto g4EndPos = aSecondary->GetPosition();
  p.endVertex.position.x = g4EndPos.x() * sim::g42edm::length;
  p.endVertex.position.y = g4EndPos.y() * sim::g42edm::length;
  p.endVertex.position.z = g4EndPos.z() * sim::g42edm::length;
  p.endVertex.ctau =aSecondary->GetGlobalTime() * sim::g42edm::length;


  auto g4StartPos = initialPos;
  p.startVertex.position.x = g4StartPos.x() * sim::g42edm::length;
  p.startVertex.position.y = g4StartPos.y() * sim::g42edm::length;
  p.startVertex.position.z = g4StartPos.z() * sim::g42edm::length;
  p.startVertex.ctau = (aSecondary->GetGlobalTime() - aSecondary->GetLocalTime()) * sim::g42edm::length;
  m_particleData->emplace_back(p);
}
}
