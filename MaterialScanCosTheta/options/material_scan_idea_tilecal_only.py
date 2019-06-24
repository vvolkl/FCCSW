import os
from Gaudi.Configuration import *


from Configurables import FCCDataSvc
podioevent   = FCCDataSvc("EventDataSvc")

from Configurables import GeoSvc
geoservice = GeoSvc("GeoSvc")
geoservice.detectors = [
  'file:/home/vali/repo/FCCSW/Detector/DetFCCeeHCalTile/compact/standalone.xml',
  ]

from Configurables import MaterialScanCosTheta
materialservice = MaterialScanCosTheta()
materialservice.filename="fccee-lar-materialscan_idea-tilecal.root"
materialservice.cosThetaBinning=0.05
materialservice.cosThetaMax=0.99
materialservice.nPhiTrials=1

from Configurables import PodioOutput
## PODIO algorithm
out = PodioOutput("out", OutputLevel=DEBUG)

from Configurables import ApplicationMgr
ApplicationMgr( TopAlg = [out],
                EvtSel = 'NONE',
                EvtMax   = 1,
                # order is important, as GeoSvc is needed by SimG4Svc
                ExtSvc = [podioevent, geoservice, materialservice],
                OutputLevel=DEBUG
 )
