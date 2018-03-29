import os
import argparse
from GaudiKernel import SystemOfUnits as units
from GaudiKernel import PhysicalConstants as constants
from Gaudi.Configuration import *

parser = argparse.ArgumentParser()
parser.add_argument('--inputfile', type=str, default='muons_for_seeding.root', help='specify an input file')
parser.add_argument('--pileupfiles',nargs="+", type=str, default='', help='specify an input file')
parser.add_argument('--outputfile', type=str, default='muons_for_seeding_tracks.root', help='specify an output file')
parser.add_argument('--nevents', type=int, default=10, help='specify number of events to process')
parser.add_argument('--npileup', type=int, default=0, help='specify number of events to process')
parser.add_argument('--noSignal', action="store_true",  help='use old collection names')
parser.add_argument('--cleanHits', action="store_true",  help='use old collection names')
parser.add_argument('--legacyCollectionNames', action="store_true",  help='use old collection names')
parser.add_argument('--overlayCollectionNames', action="store_true",  help='use old collection names')
parser.add_argument('--trajectories', action="store_true",  help='read trajectories')
parser.add_argument('--truthseeding', action="store_true",  help='use truth seeding toolG')
parser.add_argument('--recoHelix', action="store_true",  help='run recoHelix')

args, _ = parser.parse_known_args()

inputfile = args.inputfile

algList = []

from Configurables import AuditorSvc, ChronoAuditor
chra = ChronoAuditor()
audsvc = AuditorSvc()
audsvc.Auditors = [chra]

# easily switch between different branch naming conventions
collectionNames = {}
if args.legacyCollectionNames:
  collectionNames["GenParticles"] = "allGenParticles"
  collectionNames["GenVertices"] = "allGenVertices"
  collectionNames["SimParticles"] = "simParticles"
  collectionNames["SimVertices"] = "simVertices"
  collectionNames["TrackerPositionedHits"] = "positionedHits"
  collectionNames["TrackerHits"] = "hits"
elif args.overlayCollectionNames:
  collectionNames["GenParticles"] = "mergedGenParticles"
  collectionNames["GenVertices"] = "mergedGenVertices`"
  collectionNames["SimVertices"] = "mergedSimVertices"
  collectionNames["SimParticles"] = "mergedSimParticles"
  collectionNames["TrackerPositionedHits"] = "mergedTrackerPositionedHits"
  collectionNames["TrackerHits"] = "mergedTrackerHits"
else:
  collectionNames["GenParticles"] = "GenParticles"
  collectionNames["GenVertices"] = "GenVertices"
  collectionNames["SimParticles"] = "SimParticles"
  collectionNames["SimVertices"] = "SimVertices"
  collectionNames["TrackerPositionedHits"] = "TrackerPositionedHits"
  collectionNames["TrackerHits"] = "TrackerHits"
if args.trajectories:
  collectionNames["trajectory"] = "trajectory"
  collectionNames["trajectoryPoints"] = "trajectoryPoints"

if not args.noSignal:
  from Configurables import PodioInput
  podioinput = PodioInput("PodioReader", 
                          collections=collectionNames.values(),
                            OutputLevel=DEBUG,
                            )
  podioinput.AuditExecute = True
  algList += [podioinput]
print "# pileup: ", args.npileup
if args.npileup > 0:
  pileupFilenames = args.pileupfiles

  from Configurables import PileupParticlesMergeTool
  particlemergetool = PileupParticlesMergeTool("PileupGenParticlesMergeTool")
  # branchnames for the pileup
  particlemergetool.genParticlesBranch = "GenParticles"
  particlemergetool.genVerticesBranch = "GenVertices"
  # branchnames for the signal
  particlemergetool.signalGenVertices.Path = "GenVertices"
  particlemergetool.signalGenParticles.Path = "GenParticles"
  # branchnames for the output
  particlemergetool.mergedGenParticles.Path = "mergedGenParticles"
  particlemergetool.mergedGenVertices.Path = "mergedGenVertices"

  simparticlemergetool = PileupParticlesMergeTool("PileupSimParticlesMergeTool")
  # branchnames for the pileup
  simparticlemergetool.genParticlesBranch = "SimParticles"
  simparticlemergetool.genVerticesBranch = "SimVertices"
  # branchnames for the signal
  simparticlemergetool.signalGenVertices.Path = "SimVertices"
  simparticlemergetool.signalGenParticles.Path = "SimParticles"
  # branchnames for the output
  simparticlemergetool.mergedGenParticles.Path = "mergedSimParticles"
  simparticlemergetool.mergedGenVertices.Path = "mergedSimVertices"

  from Configurables import PileupTrackHitMergeTool
  trackhitsmergetool = PileupTrackHitMergeTool("PileupTrackHitMergeTool")
  # branchnames for the pileup
  trackhitsmergetool.pileupHitsBranch = "TrackerHits"
  trackhitsmergetool.pileupPositionedHitsBranch = "TrackerPositionedHits"
  # branchnames for the signal
  trackhitsmergetool.signalHits = "TrackerHits"
  trackhitsmergetool.signalPositionedHits = "TrackerPositionedHits"
  # branchnames for the output
  trackhitsmergetool.mergedHits = "mergedTrackerHits"
  trackhitsmergetool.mergedPositionedHits = "mergedTrackerPositionedHits"

  # use the pileuptool to specify the number of pileup
  from Configurables import ConstPileUp
  pileuptool = ConstPileUp("ConstPileUp", numPileUpEvents=args.npileup)

  # algorithm for the overlay
  from Configurables import PileupOverlayAlg
  overlay = PileupOverlayAlg()
  overlay.pileupFilenames = pileupFilenames
  overlay.randomizePileup = False
  overlay.noSignal = args.noSignal
  overlay.mergeTools = [particlemergetool, simparticlemergetool, trackhitsmergetool]
  overlay.PileUpTool = pileuptool
  overlay.AuditExecute = True
  algList += [overlay]

from Configurables import FCCDataSvc
podioevent   = FCCDataSvc("EventDataSvc", input=args.inputfile)




# TrickTrack Seeding Configuration
from Configurables import FastHitFilterTool
hitfiltertool = FastHitFilterTool("FastHitFilterTool")

from Configurables import FullInnerLayerGraphTool
layergraphtool = FullInnerLayerGraphTool()


from Configurables import TrickTrackSeedingTool
tricktrack_seed_tool = TrickTrackSeedingTool()
tricktrack_seed_tool.LayerGraphTool = layergraphtool
tricktrack_seed_tool.deltaZ=860
tricktrack_seed_tool.deltaT=3000000.0
tricktrack_seed_tool.deltaPhi=1.4
tricktrack_seed_tool.seedingLayerIndices = [[0, 0], [0,1], [0,2], [0,3], [2,0], [2,2], [2,4], [2,6], [2,1],[2,3],[2,5],[2,7],[2,9],[2,8]]

tricktrack_seed_tool.ptMin = 0.0
tricktrack_seed_tool.phiCut = 12.80
#tricktrack_seed_tool.thetaCut = 2.
#tricktrack_seed_tool.hardPtCut= 100.0
tricktrack_seed_tool.regionOriginRadius= 120.5

tricktrack_seed_tool.cleanHits = args.cleanHits


# Alternative: TruthSeeding
from Configurables import TruthSeedingTool
truth_seeds = TruthSeedingTool()

from Configurables import RecTrackAlg
RecTrackAlg = RecTrackAlg()
if args.truthseeding:
  RecTrackAlg.TrackSeedingTool = truth_seeds
else:
  RecTrackAlg.TrackSeedingTool = tricktrack_seed_tool
RecTrackAlg.TrackerPositionedHits.Path = collectionNames["TrackerPositionedHits"]
RecTrackAlg.OutputLevel = DEBUG




from Configurables import PodioOutput
out = PodioOutput("out",
                   filename=args.outputfile,
                   OutputLevel=DEBUG)
out.outputCommands = ["keep *"]


RecTrackAlg.AuditExecute = True
out.AuditExecute = True
algList += [RecTrackAlg]
if args.recoHelix:
  from Configurables import RecHelixTrajectory
  rht = RecHelixTrajectory()
  rht.AuditExecute = True
  algList += [rht]
algList += [out]

from Configurables import ApplicationMgr
ApplicationMgr( TopAlg =  algList,
                EvtSel = 'NONE',
                EvtMax   = args.nevents,
                ExtSvc = [podioevent],
                OutputLevel=DEBUG,
 )
