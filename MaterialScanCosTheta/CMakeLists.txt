gaudi_subdir(MaterialScanCosTheta v1r0)

gaudi_depends_on_subdirs(GaudiKernel
                         FWCore
                         Detector/DetInterface
                         Detector/DetSegmentation
                         Detector/DetCommon)

find_package(ROOT COMPONENTS MathCore GenVector Geom Tree Physics RIO Hist)

find_package(DD4hep COMPONENTS DDRec REQUIRED)

find_package(Geant4)

gaudi_add_module(MaterialScanCosTheta
                 src/*.cpp
                 INCLUDE_DIRS GaudiKernel ROOT DD4hep  FWCore DetInterface DetSegmentation DetCommon  Geant4
                 LINK_LIBRARIES GaudiKernel ROOT DD4hep ${DD4hep_COMPONENT_LIBRARIES} DetSegmentation DetCommon FWCore Geant4 )
