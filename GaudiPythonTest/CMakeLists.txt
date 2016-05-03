gaudi_subdir(GaudiPythonTest v1r0)

# this declaration will not be needed in the future
gaudi_depends_on_subdirs(GaudiAlg FWCore)

find_package(FCCEDM)
find_package(PODIO)

gaudi_add_module(GaudiPythonTest
                 *.cpp
                 INCLUDE_DIRS FCCEDM PODIO
                 LINK_LIBRARIES FCCEDM GaudiAlgLib PODIO)

