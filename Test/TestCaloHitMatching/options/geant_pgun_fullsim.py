from Gaudi.Configuration import *

from Configurables import FCCDataSvc
## Data service
podioevent = FCCDataSvc("EventDataSvc")

from Configurables import GenAlg, MomentumRangeParticleGun
guntool = MomentumRangeParticleGun(PdgCodes=[13], ThetaMin=1.5 , ThetaMax=1.59, PhiMin=0.3,
PhiMax=.7,MomentumMin=100000)
gen = GenAlg("ParticleGun", SignalProvider=guntool, VertexSmearingTool="FlatSmearVertex")
gen.hepmc.Path = "hepmc"

from Configurables import HepMCToEDMConverter
hepmc_converter = HepMCToEDMConverter("Converter")
hepmc_converter.hepmc.Path="hepmc"
hepmc_converter.genparticles.Path="allGenParticles"
hepmc_converter.genvertices.Path="allGenVertices"

from Configurables import GeoSvc
geoservice = GeoSvc("GeoSvc", detectors=['file:Detector/DetFCChhBaseline1/compact/FCChh_DectEmptyMaster.xml',
                                         'file:Detector/DetFCChhTrackerSimple/compact/Tracker.xml',
                                         'file:Detector/DetFCChhECalInclined/compact/FCChh_ECalBarrel_withCryostat.xml',
                                         #'file:Detector/DetFCChhCalDiscs/compact/Endcaps_coneCryo.xml',
                                         #'file:Detector/DetFCChhCalDiscs/compact/Forward_coneCryo.xml',
                                         ],
                    OutputLevel = INFO)

from Configurables import SimG4Svc, SimG4FullSimActions
actions = SimG4FullSimActions()
actions.enableHistory=True

from Configurables import SimG4ConstantMagneticFieldTool
field = SimG4ConstantMagneticFieldTool(
    "SimG4ConstantMagneticFieldTool", FieldOn=True, IntegratorStepper="ClassicalRK4")

from Configurables import SimG4Svc
## Geant4 service
# Configures the Geant simulation: geometry, physics list and user actions
geantservice = SimG4Svc("SimG4Svc", detector='SimG4DD4hepDetector', physicslist="SimG4FtfpBert",
                        actions=actions, magneticField=field)


from Configurables import SimG4Alg, SimG4SaveTrackerHits, SimG4SaveCalHits, SimG4PrimariesFromEdmTool
## Geant4 algorithm
# Translates EDM to G4Event, passes the event to G4, writes out outputs via tools
# first, create a tool that saves the tracker hits
# Name of that tool in GAUDI is "XX/YY" where XX is the tool class name ("SimG4SaveTrackerHits")
# and YY is the given name ("saveTrackerHits")

from Configurables import SimG4SaveParticleHistory
savehisttool = SimG4SaveParticleHistory("saveHistory")
savehisttool.mcParticles.Path = "simParticles"
savehisttool.genVertices.Path = "simVertices"

savetrackertool = SimG4SaveTrackerHits("saveTrackerHits", readoutNames = ["TrackerBarrelReadout", "TrackerEndcapReadout"])
savetrackertool.positionedTrackHits.Path = "positionedHits"
savetrackertool.trackHits.Path = "hits"
# and a tool that saves the calorimeter hits with a name "SimG4SaveCalHits/saveCalHits"
saveecaltool = SimG4SaveCalHits("saveECalBarrelHits", readoutNames = ["ECalBarrelEta"])
saveecaltool.positionedCaloHits.Path = "ECalBarrelPositionedHits"
saveecaltool.caloHits.Path = "ECalBarrelHits"
#saveendcaptool = SimG4SaveCalHits("saveECalEndcapHits", readoutNames = ["EMECPhiEta"])
#saveendcaptool.positionedCaloHits.Path = "ECalEndcapPositionedHits"
#saveendcaptool.caloHits.Path = "ECalEndcapHits"
#savefwdtool = SimG4SaveCalHits("saveECalFwdHits", readoutNames = ["EMFwdPhiEta"])
#savefwdtool.positionedCaloHits.Path = "ECalFwdPositionedHits"
#savefwdtool.caloHits.Path = "ECalFwdHits"
# next, create the G4 algorithm, giving the list of names of tools ("XX/YY")

particle_converter = SimG4PrimariesFromEdmTool("EdmConverter")
particle_converter.genParticles.Path = "allGenParticles"
geantsim = SimG4Alg("SimG4Alg",
                    outputs = [ "SimG4SaveTrackerHits/saveTrackerHits", 
                                "SimG4SaveCalHits/saveECalBarrelHits", 
                                #"SimG4SaveCalHits/saveECalEndcapHits", 
                                #"SimG4SaveCalHits/saveECalFwdHits", 
                                "SimG4SaveParticleHistory/saveHistory"
                              ],
                    eventProvider=particle_converter)

from Configurables import PodioOutput
out = PodioOutput("out",
                   OutputLevel=DEBUG)
out.outputCommands = ["keep *"]
out.filename = "output_calo_hit_matching.root"

from Configurables import ApplicationMgr
ApplicationMgr( TopAlg=[gen, hepmc_converter, geantsim, out],
                EvtSel='NONE',
                EvtMax=3,
                ## order is important, as GeoSvc is needed by SimG4Svc
                ExtSvc=[podioevent, geoservice, geantservice],
                OutputLevel=INFO
 )
