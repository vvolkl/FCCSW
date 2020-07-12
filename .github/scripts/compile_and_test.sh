#!/bin/bash

cd /Package
source init.sh
ls /cvmfs/geant4.cern.ch/share/data/G4ENSDFSTATE2.2
export G4ENSFGSTATEDATA=/cvmfs/geant4.cern.ch/share/data/G4ENSDFSTATE2.2
export G4PARTICLEXSDATA=/cvmfs/geant4.cern.ch/share/data/G4PARTICLEXS2.1
export G4NEUTRONHPDATA=/cvmfs/geant4.cern.ch/share/data/G4NDL4.6
export G4LEDATA=/cvmfs/geant4.cern.ch/share/data/G4EMLOW7.9.1
export G4LEVELGAMMADATA=/cvmfs/geant4.cern.ch/share/data/PhotonEvaporation5.5
export G4RADIOACTIVEDATA=/cvmfs/geant4.cern.ch/share/data/RadioactiveDecay5.4
export G4PIIDATA=/cvmfs/geant4.cern.ch/share/data/G4PII1.3
export G4REALSURFACEDATA=/cvmfs/geant4.cern.ch/share/data/RealSurface2.1.1
export G4SAIDXSDATA=/cvmfs/geant4.cern.ch/share/data/G4SAIDDATA2.0
export G4ABLADATA=/cvmfs/geant4.cern.ch/share/data/G4ABLA3.1
export G4INCLDATA=/cvmfs/geant4.cern.ch/share/data/G4INCL1.0

mkdir build install
cd build
cmake -DCMAKE_INSTALL_PREFIX=../install -DCMAKE_CXX_STANDARD=17  .. && \
make  -j `getconf _NPROCESSORS_ONLN` && \
make install && \
ls -alhR $G4LEDATA && \
ctest -j `getconf _NPROCESSORS_ONLN` --output-on-failure
