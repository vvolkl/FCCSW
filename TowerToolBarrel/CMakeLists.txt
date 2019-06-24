gaudi_subdir(TowerToolBarrel v1r0)

gaudi_depends_on_subdirs(GaudiAlg FWCore Detector/DetInterface Detector/DetSegmentation Detector/DetCommon Reconstruction/RecInterface)



find_package(FCCEDM)
find_package(PODIO)
find_package(ROOT COMPONENTS MathCore GenVector Geom Tree Physics RIO Hist)

find_package(DD4hep COMPONENTS DDRec REQUIRED)

find_package(Geant4)


gaudi_add_module(TowerToolBarrel
                 src/*.cpp
                 INCLUDE_DIRS GaudiKernel ROOT DD4hep  FWCore DetInterface DetSegmentation DetCommon FCCEDM RecInterface PODIO Geant4
                 LINK_LIBRARIES GaudiKernel ROOT DD4hep ${DD4hep_COMPONENT_LIBRARIES} DetSegmentation DetCommon FWCore FCCEDM GaudiAlgLib PODIO Geant4)
