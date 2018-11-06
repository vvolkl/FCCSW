################################################################################
# Package: TestCaloHitMatching
################################################################################
gaudi_subdir(TestCaloHitMatching v1r0)

gaudi_depends_on_subdirs(GaudiAlg GaudiKernel FWCore Detector/DetInterface Detector/DetCommon Detector/DetSegmentation)
find_package(Geant4)
find_package(ACTS COMPONENTS Core)
find_package(DD4hep)

include(${Geant4_USE_FILE})

#gaudi_add_module(TestTestPlugins
#                 src/*.cpp
#                 INCLUDE_DIRS FWCore DetInterface TrackingUtils Geant4 DetCommon DetSegmentation GaudiKernel ACTS::ACTSCore
#                 LINK_LIBRARIES GaudiAlgLib FWCore Geant4 TrackingUtils DetCommon DetSegmentation ACTS::ACTSCore)

#include(CTest)
#gaudi_add_test(LookForNeighbours
#               WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/Test/TestCaloHitMatching/tests/
#               FRAMEWORK tests/options/testneighbours.py)

