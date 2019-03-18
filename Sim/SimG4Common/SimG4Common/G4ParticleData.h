#ifndef SIMG4COMMON_G4PARTICLEDATA_H
#define SIMG4COMMON_G4PARTICLEDATA_H

#include "datamodel/LorentzVector.h"
#include "datamodel/GenVertexData.h"

namespace fcc {

struct G4ParticleData {

  LorentzVector p4Initial;
  LorentzVector p4Final;
  int pdgId;
  unsigned int trackId;
  unsigned int motherId;
  GenVertexData startVertex;
  GenVertexData endVertex;
};

}


#endif
