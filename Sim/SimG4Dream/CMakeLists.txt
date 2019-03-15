################################################################################
# Package: SimG4Dream
################################################################################
gaudi_subdir(SimG4Dream v1r0)

# this declaration will not be needed in the future
gaudi_depends_on_subdirs(GaudiAlg FWCore Sim/SimG4Common Detector/DetInterface Sim/SimG4Interface Detector/DetCommon)

find_package(Geant4)

find_package(FCCEDM)
find_package(PODIO)

include(${Geant4_USE_FILE})
find_package(ROOT COMPONENTS MathCore GenVector Geom REQUIRED)

gaudi_install_headers(SimG4Dream)
gaudi_install_python_modules()


gaudi_add_library(SimG4Dream
                 src/*.cc
                 INCLUDE_DIRS FWCore Geant4 ROOT FCCEDM PODIO
                 LINK_LIBRARIES GaudiAlgLib Geant4 ROOT  FCCEDM PODIO
                 PUBLIC_HEADERS SimG4Dream)

gaudi_add_module(SimG4DreamComponents
                 components/*.cpp
                 INCLUDE_DIRS Geant4 FWCore SimG4Common SimG4Dream SimG4Interface DetCommon  ROOT
                 LINK_LIBRARIES GaudiAlgLib SimG4Dream Geant4 FWCore SimG4Common DetCommon  ROOT)


include(CTest)

gaudi_add_test(DreamCaloGeantinoScan
               WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
               FRAMEWORK tests/options/geant_fullsim_gdml.py)
