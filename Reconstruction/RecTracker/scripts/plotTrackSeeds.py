
import os
import sys
import matplotlib
matplotlib.use("Agg")
import ROOT
import argparse
from EventStore import EventStore
import numpy as np

import matplotlib.pyplot as plt


# command line arguments
parser = argparse.ArgumentParser()
parser.add_argument("filename", help="edm file to read")
parser.add_argument("--no_process_hits", action="store_true", help="process hits from edm file")
parser.add_argument("--no_process_tracks", action="store_true", help="process tracks from edm file")
parser.add_argument("--no_process_genparticles", action="store_true", help="process genparticles from edm file")
parser.add_argument("--no_process_trajectories", action="store_true", help="process trajectores from edm file")
parser.add_argument("--no_process_rechelix", action="store_true", help="process trajectores from edm file")
parser.add_argument("--plot_event", help="event for which to create plot", nargs="+", type=int, default=[0])
parser.add_argument("--rCut", help="rcut determining whether hits are drawn", type=float, default=800)
parser.add_argument("--zCut", help="zcut determining whether hits are drawn", type=float, default=2000)
parser.add_argument("--rlimlow", help="event for which to create plot", type=float, default=0)
parser.add_argument("--rlimhigh", help="event for which to create plot", type=float, default=450)
parser.add_argument("--zlimlow", help="event for which to create plot", type=float, default=-2000)
parser.add_argument("--zlimhigh", help="event for which to create plot", type=float, default=2000)
parser.add_argument("--nevents", help="max number of events to process", type=int, default=1000)
parser.add_argument("--legacyCollectionNames", action="store_true",  help="create fake rate histograms", )
parser.add_argument("--overlayCollectionNames", action="store_true",  help="create fake rate histograms", )
parser.add_argument("--plotprefix", help="where to store the plots", type=str, default=os.environ.get("FCCPLOTS", ""))
args = parser.parse_args()

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
  collectionNames["SimParticles"] = "mergedSimParticles"
  collectionNames["TrackerPositionedHits"] = "mergedTrackerPositionedHits"
else:
  collectionNames["GenParticles"] = "GenParticles"
  collectionNames["SimParticles"] = "SimParticles"
  collectionNames["TrackerPositionedHits"] = "TrackerPositionedHits"


#TODO
def Particle2Track(p, vertex, charge=1., B=4.):
    pT = np.sqrt(p[0]**2 + p[1]**2)
    pX = p[0]
    pY = p[1]
    pZ = p[2]
    q_pT = charge / pT
    R = pT / (0.003 * charge * B) 
    print R, np.sqrt(vertex[0]**2 + vertex[1]**2)
    d0 = -charge*(np.sqrt(vertex[0]**2 + vertex[1]**2))
    phi0 = np.arctan(pY / pX)
    cotTheta = pZ / pT
    z0 = vertex[2] 
    return [d0, z0, phi0, cotTheta, q_pT]

#TODO
def Track2Particle(params):
    d0, z0, phi0, cotTheta, q_pT = params
    charge = np.sign(q_pT)
    pT = charge / q_pT
    pZ  =  cotTheta * pT
    vertex = [d0*np.cos(phi0), d0*np.sin(phi0), z0]
    p = [pT*np.cos(phi0), pT*np.sin(phi0), pZ]
    return [p, vertex]

#TODO
def helix(params, time_points=(0,2,2000), B=4.):
  path_points = np.linspace(*time_points)
  d0, z0, phi0, cotTheta, q_pT = params
  charge = np.sign(q_pT)
  # motion in z linear, unchanged by magnetic field in z
  # cartesian coordinates of point of closest approach
  x0 =  d0 * np.cos(phi0)
  y0 = d0 * np.sin(phi0)
  rho = charge / np.abs(q_pT) / 0.0003 / B
  z = z0 +  cotTheta * np.abs(rho) * path_points 
  # cartesian coordinates of helix center
  xc = x0 + rho * np.sin(phi0)
  yc = y0 + rho * np.cos(phi0)
  # generate real space points from track parameters
  x = xc + rho * np.cos(charge*path_points + phi0 - np.pi*0.5)
  y = yc + rho * np.sin(charge*path_points + phi0 - np.pi*0.5)
  return np.array([x, y, z]).T


# plot detector cartoon 
plt.figure("xy")
ax = plt.gca()
for r in [25, 60, 100, 150]:
  ax.add_artist(plt.Circle((0,0), r,fill=False, color="blue", alpha=0.2, lw=6 ))
plt.figure("rz")
ax = plt.gca()
for z in [25, 60, 100, 150]:
  ax.add_artist(plt.Rectangle((-850,z - 3), 1700, 5 ,fill=True, color="red", alpha=0.4, lw=0 ))
for z in [950, 1178, 1462, 1813, -950, -1178, -1462, -1813]: 
  ax.add_artist(plt.Rectangle((z - 3, 25), 20, 400 ,fill=True, color="red", alpha=0.6, lw=1, ))
events = EventStore([args.filename])
print len(events),  " events in rootfile ", args.filename

filter_trajectories = []
# main event loop
#for i, store in enumerate([events[args.plot_event]]):
for i in args.plot_event:
      store = events[i]
      print "event ", i
      if not args.no_process_genparticles:
        genparticles = store.get(collectionNames["GenParticles"])
        print "processing allGenParticles ..."
        for p in genparticles:
            vertex = [p.startVertex().position().x, p.startVertex().position().y, p.startVertex().position().z] 

            if p.core().status != 201:
              momentum = [p.core().p4.px, p.core().p4.py, p.core().p4.pz, p.core().p4.mass]
              pm = ROOT.TLorentzVector(*momentum)
              print "\t sim Eta:", pm.Eta() 
              print "\tsim trackID: ", p.core().bits, "sim pdgId: ", p.core().pdgId, "momentum: ", [p.core().p4.px, p.core().p4.py, p.core().p4.pz]
              print "\tsim vertex: ", vertex
              print "\tsim phi: ", np.arctan(p.core().p4.py / p.core().p4.px)
              print "\tsim cottheta: ", p.core().p4.pz / np.sqrt(p.core().p4.px**2 +  p.core().p4.py**2)
              print "\tsim q_pT: ", 1. / np.sqrt(p.core().p4.px**2 +  p.core().p4.py**2)
              if pm.Pt() > 0 and np.abs(pm.Eta()) < 6 :
                #filter_trajectories.append(p.core().bits)
                filter_trajectories.append(1)
              charge = 1
              #vertex = [0,0,0]
              cpar = Particle2Track(momentum, vertex, charge=charge)
              print "\tsim calculated params:", cpar
              cmom = Track2Particle(cpar)
              print "\tsim calculated particle: ", cmom

      if not args.no_process_tracks:
        print "processing trackStates ..."
        trackStates = store.get("TrackStates")
        print len(trackStates)
        for ts in trackStates:
            pass
        print "processing tracks ..."
        tracks = store.get('Tracks')
        trackCounter = 0
        for t in tracks:
              if t.bits() == 1 or True:
                ts = t.states(0)
                trackparams = [
                  ts.d0(),
                  ts.z0(),
                  ts.phi(),
                  ts.theta(),
                  ts.qOverP() * -10,
                  ]

                
                pt = 1./trackparams[-1]

                print "\trec pt: ", pt
                print "\trec params: ", trackparams
                pos = []
                for j in range(t.hits_size()):
                  cor = t.hits(j).position()
                  pos.append([cor.x, cor.y, cor.z])
                pos = np.array(pos)
                plt.figure("xy")
                plt.plot(pos[:,0],pos[:,1], '--', color="black", label="Track Seed" if trackCounter == 0 else None)
                plt.figure("rz")
                plt.plot(pos[:,2], np.sqrt(pos[:,0]**2 + pos[:,1]**2), '--', color="black", label="Track Seed" if trackCounter == 0 else None)
                trackCounter +=1
                """

                helix_points = helix(trackparams)
                pos = np.array(helix_points)
                plt.figure("xy")
                plt.plot(pos[:,0],pos[:,1], '--', color="magenta", label="Track Helix" if trackCounter == 0 else None)
                plt.figure("rz")
                plt.plot(pos[:,2], np.sqrt(pos[:,0]**2 + pos[:,1]**2), '--', color="magenta", label="Track Seed" if trackCounter == 0 else None)
                """

      if not args.no_process_trajectories:
        plt.figure("xy")
        print "processing trajectories ..."
        clusters = store.get('trajectoryPoints')
        pos = []
        for c in clusters:
              #print "trajectory ID: ", c.bits()
              cor = c.position()
              if c.bits() in filter_trajectories or True:
                pos.append([cor.x, cor.y, cor.z, c.bits()])
        p = np.array(pos)
        if p.size > 0:
          p =  p[(np.linalg.norm(p[:,:2], axis=1) < args.rCut)]
          p = p[np.abs(p[:,2]) < args.zCut] 
          for _id in np.unique(p[:,3]):
            maxpoints = -1
            pf = p[p[:,3] == _id]
            plt.figure("xy")
            plt.plot(pf[:,0],pf[:,1], '-', color="green", alpha=0.3)#, label="MCTruth trajectory")
            plt.figure("rz")
            plt.plot(pf[:,2], np.sqrt(pf[:,0]**2 + pf[:,1]**2), '-', color="green", alpha=0.3)#, label="MCTruth trajectory")

      if not args.no_process_rechelix:
        plt.figure("xy")
        print "processing recHelix ..."
        clusters = store.get('RecHelixPoints')
        pos = []
        for c in clusters:
              #print "trajectory ID: ", c.bits()
              cor = c.position()
              if c.bits() in filter_trajectories or True:
                pos.append([cor.x, cor.y, cor.z, c.bits()])
        p = np.array(pos)
        if p.size > 0:
          p =  p[(np.linalg.norm(p[:,:2], axis=1) < args.rCut)]
          p = p[np.abs(p[:,2]) < args.zCut] 
          for _id in np.unique(p[:,3]):
            maxpoints = -1
            pf = p[p[:,3] == _id]
            plt.figure("xy")
            plt.plot(pf[:,0],pf[:,1], '--', color="purple", alpha=0.3)#, label="MCTruth trajectory")
            plt.figure("rz")
            plt.plot(pf[:,2], np.sqrt(pf[:,0]**2 + pf[:,1]**2), '--', color="purple", alpha=0.3)#, label="MCTruth trajectory")

      if not args.no_process_hits:
        print "processing hits ..."
        hits = store.get(collectionNames['TrackerPositionedHits'])
        pos = []
        for c in hits:
            cor = c.position()
            if np.abs(cor.z) < args.zCut:
              #if c.bits() in filter_trajectories:
              pos.append([cor.x, cor.y, cor.z])
        if pos:
          p = np.array(pos)
          p =  p[(np.linalg.norm(p[:,:2], axis=1) < args.rCut)]
          plt.figure("xy")
          plt.plot(p[:,0],p[:,1], '.', color="midnightblue", label="Tracker Hit")
          plt.figure("rz")
          plt.plot(p[:,2], np.sqrt(p[:,0]**2 + p[:,1]**2), '.', color="midnightblue", label="Tracker Hit")


f = args.filename.replace(".root", "") + "_event_" + str(args.plot_event) + "_trackDisplay"
plt.figure("xy")
plt.xlim(-1* args.rlimhigh, args.rlimhigh)
plt.ylim(-1* args.rlimhigh, args.rlimhigh)
plt.xlabel("x [mm]")
plt.ylabel("y [mm]")
#plt.legend(loc="best")
#plt.savefig(os.path.join(args.plotprefix, f + "_xy.png"))
#plt.savefig(os.path.join(args.plotprefix, f + "_xy.pdf"))
plt.figure("rz")
plt.xlim(-1* args.zlimhigh, args.zlimhigh)
plt.ylim(args.rlimlow, args.rlimhigh)
plt.xlabel("z [mm]")
plt.ylabel("r [mm]")
#plt.legend(loc="best")
#plt.savefig(os.path.join(args.plotprefix, f + "_rz.png"))
#plt.savefig(os.path.join(args.plotprefix, f + "_rz.pdf"))

figs = [plt.figure(n) for n in plt.get_fignums()]
for fig in figs:
  fig.savefig(os.path.join(args.plotprefix, f + "_" + fig._label + ".png"))

