#ifndef SIMG4COMMON_EVENTINFORMATION_H
#define SIMG4COMMON_EVENTINFORMATION_H

#include "G4VUserEventInformation.hh"
#include "G4LorentzVector.hh"

#include "SimG4Common/G4ParticleData.h"

#include <iostream>
#include <map>
#include <vector>

class G4Track;
namespace fcc {
class GenVertexCollection;
class MCParticleCollection;
}

/** @class sim::EventInformation SimG4Common/SimG4Common/EventInformation.h EventInformation.h
 *
 * Additional event information.
 *
 * Currently holds the particle history in form of edm particles and vertices
 *
 * @author J. Lingemann
 */

namespace sim {
class EventInformation : public G4VUserEventInformation {
public:
  /// Default constructor
  EventInformation();
  /// Destructor
  virtual ~EventInformation() = default;
  /** Set external pointers to point at the particle and vertex collections.
   * @param[in] aGenVertexCollection  pointer to a collection that should take ownership of the particles saved here
   * @param[in] aMCParticleCollection  pointer to a collection that should take ownership of the particles saved here
   */
  void setCollections(std::vector<fcc::G4ParticleData>*& aParticleDataVector);
  /// Add a particle to be tracked in the EDM collections
  void addParticle(const G4Track* aSecondary, G4LorentzVector initialPos, G4LorentzVector
  initialEnergy);

  void Print() const {};

private:
  std::vector<fcc::G4ParticleData>* m_particleData;
};
}
#endif /* define SIMG4COMMON_EVENTINFORMATION_H */
