from Gaudi.Configuration import *
from Configurables import PodioTest
from Configurables import PodioOutput

podiotest = PodioTest()
podiotest.DataOutputs.trackClusters.Path = "clusters"


from Configurables import ApplicationMgr
from Configurables import FCCDataSvc
podioevent   = FCCDataSvc("EventDataSvc")
## PODIO algorithm
out = PodioOutput("out", OutputLevel=DEBUG)
out.outputCommands = ["keep *"]

ApplicationMgr(
    # all algorithms should be put here
    TopAlg=[podiotest, out],
    EvtSel='NONE',
    # number of events
    EvtMax=1000,
    # all services should be put here
    ExtSvc=[podioevent],
    # possible: DEBUG, INFO, WARNING, ERROR, FATAL
    OutputLevel=INFO
)

from GaudiPython import AppMgr, gbl, Helper, InterfaceCast
from cppyy import addressof, bind_object
import numpy as np
gaudi = AppMgr()
gaudi.run(1)


# access transient event store
evt = gaudi.evtSvc()
unreadable_dataobject = evt['/Event/clusters']
_datasvc = Helper.service(gbl.Gaudi.svcLocator(), "EventDataSvc")
podiosvc = InterfaceCast(gbl.PodioDataSvc)(_datasvc)
coll_vec = podiosvc.getCollections()
assert coll_vec.size() > 0
coll = coll_vec.at(0) 
coll_name, coll_base = coll.first, coll.second
print 'collection name: ', coll.first
coll_buffptr = coll_base.getBufferAddress()
trackclusterdata_vec = bind_object(coll_buffptr[0], gbl.std.vector(gbl.fcc.TrackClusterData))
trackclusterdata = trackclusterdata_vec.at(100)
print 'some cluster coordinate: ', trackclusterdata.Core.position.X


# use tool
toolsvc = gaudi.toolsvc()
pileuptool = toolsvc.create("PoissonPileUp", interface=gbl.IPileUpTool)
print 'tool created ...'
print 'some tool method: ', pileuptool.numberOfPileUp()
