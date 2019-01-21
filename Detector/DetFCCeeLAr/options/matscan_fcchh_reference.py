import os
from Gaudi.Configuration import *

from Configurables import GeoSvc

from Configurables import FCCDataSvc
podioevent   = FCCDataSvc("EventDataSvc")
geoservice = GeoSvc("GeoSvc")
geoservice.detectors=[
                      'file:Detector/DetFCCeeLAr/compact/hh_reference_ecal_master.xml',
                      ]
geoservice.OutputLevel = DEBUG

from Configurables import MaterialScan
# Material scan is done from the interaction point to the end of world volume.
# In order to use other end boundary, please provide the name of a thin, e.g. cylindrical volume.
# For instance adding envelopeName="BoundaryPostCalorimetry" will perform the scan only till the end of calorimetry.
# BoundaryPostCalorimetry is defined in Detector/DetFCChhECalInclined/compact/envelopePreCalo.xml
materialservice = MaterialScan()
materialservice.filename="matscan_fcchh_ecal.root"
materialservice.etaBinning=0.05
materialservice.etaMax=6
materialservice.nPhiTrials=10

from Configurables import PodioOutput
## PODIO algorithm
out = PodioOutput("out")
out.OutputLevel = DEBUG

from Configurables import ApplicationMgr
ApplicationMgr( TopAlg = [out],
                EvtSel = 'NONE',
                EvtMax   = 1,
                # order is important, as GeoSvc is needed by SimG4Svc
                ExtSvc = [podioevent, geoservice, materialservice],
                OutputLevel=DEBUG
 )
