

from Gaudi.Configuration import *
import GaudiKernel.SystemOfUnits as units

from Configurables import ApplicationMgr
app = ApplicationMgr()
app.EvtSel = "NONE"
app.EvtMax = 100
app.ExtSvc = []
app.TopAlg = []
app.OutputLevel = DEBUG


from Configurables import FCCDataSvc
podioevent = FCCDataSvc("EventDataSvc")
app.ExtSvc += [podioevent]

from Configurables import MomentumRangeParticleGun
from GaudiKernel import PhysicalConstants as constants
guntool = MomentumRangeParticleGun()
guntool.ThetaMin = 0 
guntool.ThetaMax = 2 * constants.pi 
guntool.PdgCodes = [11]
from Configurables import GenAlg
gen = GenAlg()
gen.SignalProvider=guntool
gen.hepmc.Path = "hepmc"
app.TopAlg += [gen]

# reads an HepMC::GenEvent from the data service and writes a collection of EDM Particles
from Configurables import HepMCToEDMConverter
hepmc_converter = HepMCToEDMConverter("Converter")
hepmc_converter.hepmc.Path="hepmc"
hepmc_converter.genparticles.Path="allGenParticles"
hepmc_converter.genvertices.Path="allGenVertices"
app.TopAlg += [hepmc_converter]
