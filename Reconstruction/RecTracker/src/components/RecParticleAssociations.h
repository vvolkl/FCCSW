#ifndef RECTRACKER_RECPARTICLEASSOCIATIONS_H
#define RECTRACKER_RECPARTICLEASSOCITAIONS_H

// GAUDI
#include "GaudiAlg/GaudiAlgorithm.h"
#include "GaudiKernel/RndmGenerators.h"
#include "GaudiKernel/ToolHandle.h"

// FCCSW
#include "FWCore/DataHandle.h"


namespace fcc {
class TrackHitCollection;
class TrackCollection;
class TrackStateCollection;
class PositionedTrackHitCollection;
class ParticleCollection;
class MCParticleCollection;
class ParticleMCParticleAssociationCollection;
}

class RecParticleAssociations : public GaudiAlgorithm {
public:
  RecParticleAssociations(const std::string& name, ISvcLocator* svcLoc);

  ~RecParticleAssociations() = default;

  StatusCode initialize() override final;

  StatusCode execute() override final;

  StatusCode finalize() override final;

private:
  /// TrackHits as Input to the track seeding
  DataHandle<fcc::ParticleCollection> m_recParticles{"TrackRecoParticles", Gaudi::DataHandle::Reader, this};
  DataHandle<fcc::MCParticleCollection> m_simParticles{"SimParticles", Gaudi::DataHandle::Reader, this};
  DataHandle<fcc::ParticleMCParticleAssociationCollection> m_associations{"SimRecoParticleAssociation", Gaudi::DataHandle::Writer, this};
};

#endif /* RECTRACKER_RECPARTICLEASSOCIATIONS_H */
