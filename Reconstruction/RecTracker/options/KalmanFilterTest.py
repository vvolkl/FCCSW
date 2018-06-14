from Gaudi.Configuration import *

from Configurables import FCCDataSvc
podioevent   = FCCDataSvc("EventDataSvc", input="output_trkExtrapolationTest.root")

from Configurables import PodioInput
podioinput = PodioInput("PodioReader", 
                        collections= [
                                      "ExtrapolatedTrackstates",
                                      "SimTrackerPositionedHits",
                                      ],
                          OutputLevel=DEBUG,
                          )



from Configurables import KalmanFilter
kalman = KalmanFilter()
kalman.FittedTracks.Path = "FittedTracks"
kalman.TrackSeeds.Path = "ExtrapolatedTrackstates"
kalman.TrackerHits.Path = "TrackerHits"


from Configurables import ApplicationMgr
ApplicationMgr( TopAlg = [podioinput, kalman],
                EvtSel = 'NONE',
                EvtMax   = 10,
                # order is important, as GeoSvc is needed by SimG4Svc
                ExtSvc = [podioevent],
                OutputLevel=DEBUG,
 )
