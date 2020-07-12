#!/bin/bash

cd /Package
source init.sh
ls /cvmfs/geant4.cern.ch/share/data/G4ENSDFSTATE2.2
ls $G4ENSDFSTATEDATA
mkdir geant4data
cp -r $G4ENSDFSTATEDATA/* ./geant4data; export G4ENSDFSTATEDATA=$PWD/geant4data; 
cp -r $G4PARTICLEXSDATA/* ./geant4data; export G4PARTICLEXSDATA=$PWD/geant4data; 
cp -r $G4LEDATA/* ./geant4data; export G4LEDATA=$PWD/geant4data; 

mkdir build install
cd build
cmake -DCMAKE_INSTALL_PREFIX=../install -DCMAKE_CXX_STANDARD=17  .. && \
make  -j `getconf _NPROCESSORS_ONLN` && \
make install && \
ctest -j `getconf _NPROCESSORS_ONLN` --output-on-failure
