#include "SimG4Full/ParticleHistoryAction.h"

#include "SimG4Common/EventInformation.h"

#include "G4EventManager.hh"
#include "G4LorentzVector.hh"

namespace sim {

ParticleHistoryAction::ParticleHistoryAction(double aEnergyCut): m_energyCut(aEnergyCut) {}

void ParticleHistoryAction::PreUserTrackingAction(const G4Track* aTrack) {
  auto g4EvtMgr = G4EventManager::GetEventManager();
  auto evtinfo = dynamic_cast<sim::EventInformation*>(g4EvtMgr->GetUserInformation());
  m_initialPos = G4LorentzVector(aTrack->GetGlobalTime() - aTrack->GetLocalTime(), aTrack->GetVertexPosition());

  m_initialEnergy = G4LorentzVector(aTrack->GetMomentum(), aTrack->GetTotalEnergy());
}

void ParticleHistoryAction::PostUserTrackingAction(const G4Track* aTrack) {
  
  auto g4EvtMgr = G4EventManager::GetEventManager();
  auto evtinfo = dynamic_cast<sim::EventInformation*>(g4EvtMgr->GetUserInformation());
  if (selectSecondary(*aTrack, m_energyCut)) {
    evtinfo->addParticle(aTrack, m_initialPos, m_initialEnergy);
  }
  
  }

bool ParticleHistoryAction::selectSecondary(const G4Track& aTrack, double aEnergyCut) {
  G4LorentzVector p4(aTrack.GetMomentum(), aTrack.GetTotalEnergy());
  if (p4.e() < aEnergyCut) { 
    return false;
  }
  return true;
}
}
