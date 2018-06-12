from Gaudi.Configuration import *

from Configurables import FCCDataSvc
podioevent   = FCCDataSvc("EventDataSvc")




from Configurables import KalmanFilter
kalman = KalmanFilter()


from Configurables import ApplicationMgr
ApplicationMgr( TopAlg = [kalman],
                EvtSel = 'NONE',
                EvtMax   = 1,
                # order is important, as GeoSvc is needed by SimG4Svc
                ExtSvc = [],
                OutputLevel=DEBUG,
 )
