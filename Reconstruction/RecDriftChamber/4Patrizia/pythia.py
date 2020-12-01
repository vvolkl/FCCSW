

from Gaudi.Configuration import *
import GaudiKernel.SystemOfUnits as units

from Configurables import ApplicationMgr
app = ApplicationMgr()
app.EvtSel = "NONE"
app.EvtMax = 10
app.ExtSvc = []
app.TopAlg = []
app.OutputLevel = DEBUG


from Configurables import FCCDataSvc
podioevent = FCCDataSvc("EventDataSvc")
app.ExtSvc += [podioevent]

from Configurables import PythiaInterface
pythia8gentool = PythiaInterface()
### Example of pythia configuration file to generate events
pythiafilename = "Pythia_LHEinput.cmd"
#path_to_pythiafile = os.environ.get("FCCSW_SHARE_DIR", "")
#pythiafile = os.path.join(path_to_pythiafile, pythiafilename)
# Example of pythia configuration file to read LH event file
#pythiafile="options/Pythia_LHEinput.cmd"
pythia8gentool.Filename = pythiafilename
pythia8gentool.doEvtGenDecays = False
pythia8gentool.printPythiaStatistics = False

from Configurables import GenAlg
pythia8gen = GenAlg("Pythia8")
pythia8gen.SignalProvider = pythia8gentool
pythia8gen.hepmc.Path = "hepmc"
ApplicationMgr().TopAlg += [pythia8gen]

# reads an HepMC::GenEvent from the data service and writes a collection of EDM Particles
from Configurables import HepMCToEDMConverter
hepmc_converter = HepMCToEDMConverter("Converter")
hepmc_converter.hepmc.Path="hepmc"
hepmc_converter.genparticles.Path="GenParticlesUnfiltered"
hepmc_converter.genvertices.Path="GenVerticesUnfiltered"
ApplicationMgr().TopAlg += [hepmc_converter]

from Configurables import GenParticleFilter
genfilter = GenParticleFilter("StableParticles")
genfilter.accept = [1]
genfilter.allGenParticles.Path = "GenParticlesUnfiltered"
genfilter.filteredGenParticles.Path = "allGenParticles"
ApplicationMgr().TopAlg += [genfilter]
