"""
Example reconstruction analysis script for FCCSW track reconstruction.

Performs calculation of efficiencies and fake rates for FCCSW simulated data
based on the MC Truth and the reconstucted tracks.
"""

import ROOT
import argparse
# PODIO python interface: requires the fcc-edm library (libdatamodel.so in LD_LIBRARY_PATH)
# and PODIO in the python path 
from EventStore import EventStore
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats
import os



# command line arguments
parser = argparse.ArgumentParser()
parser.add_argument("filename", help="edm file to read")
parser.add_argument("--nevents", help="max number of events to process", type=int, default=1000)
parser.add_argument("--process", help="name of the process for the textbox", type=str, default="MinBias")
parser.add_argument("--ptCut", help="min pt for particles to be included in efficiency definition [GeV]", type=float, default=1.)
parser.add_argument("--etaCut", help="max eta for particles to be included in efficiency definiton ", type=float, default=5.)
parser.add_argument("--textBoxPositionX", help="textBoxPositionX", type=float, default=0.5)
parser.add_argument("--textBoxPositionY", help="textBoxPositionY", type=float, default=0.7)
parser.add_argument("--histos", action="store_true",  help="plot efficiency histograms", )
parser.add_argument("--ignore_secondaries", action="store_true",  help="don't include secondary particles in efficiency definitions", )
parser.add_argument("--process_genparticles", action="store_true",  help="use genparticles in the event loop", )
parser.add_argument("--process_trackerhits", action="store_true",  help="use trackerhits in the event loop", )
parser.add_argument("--fakeRateHistos", action="store_true",  help="create fake rate histograms", )
parser.add_argument("--legacyCollectionNames", action="store_true",  help="create fake rate histograms", )
parser.add_argument("--overlayCollectionNames", action="store_true",  help="create fake rate histograms", )
parser.add_argument("--plotprefix", help="where to store the plots", type=str, default=os.environ["FCCPLOTS"])
parser.add_argument("--checkOnlyPrimary", help="where to store the plots", action="store_true")
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

process_string = ""
if "MinBias" in args.filename:
  process_string = "MinBias, $\sqrt s = 100$ TeV"
if "Single" in args.filename:
  process_string = "Single $\mu^-$" 




boxtext = r"\textsf{\textit{\textbf{FCC}-hh Sim}, Det \textbf{v3.03}} \\  $\qquad$ " + process_string + ", n = " + str(args.nevents) + \
      r" \newline Cuts: $ |\eta| < " + str(args.etaCut) + r",  \quad p_T >  " + str(args.ptCut) + r"$ GeV"

events = EventStore([args.filename])
print len(events),  " events in rootfile ", args.filename
args.nevents = min(len(events), args.nevents)


def DelphesTrackingEfficiency(p):
  if np.abs(p.Eta()) > 6:
    return 0
  if np.abs(p.Pt()) < 0.5:
    return 0 
  if np.abs(p.Pt()) < 1:
    return 0.9
  else:
    return 0.99

# arrays to hold ntuples
delphes_tracking_efficiencies = {}
etas = {}
cotThetas = {}
pts = {}
efficiencies = []
fake_rates = []
successfully_found_track = {}
reco_track_is_fake = {}
reconstructed_track_ids = {}
reco_track_pt = {}
reco_track_eta = {}
bad_event = []
# main event loop
for i, store in enumerate(events):
      if i > args.nevents:
        break  
      ### print "Processing generated Particles from branch 'allGenParticles'" # 
      if args.process_genparticles: # currently unused
        for p in store.get(collectionNames["GenParticles"]):
            pm = [ p.core().p4.px, 
                   p.core().p4.py, 
                   p.core().p4.pz, 
                   p.core().p4.mass, ]
            pm = ROOT.TLorentzVector(*pm)
            pass

      ### print "processing trackerhits from branch 'positionedhits' " #########
      if args.process_trackerhits: # currently unused
        hits = store.get(collectionNames["TrackerPositionedHits"])
        ids = []
        for h in hits:
              ids.append(h.core().bits)
              if h.core().bits == 1:
                pass
                #print "hit r, z:", np.sqrt(h.position().x**2 + h.position().y**2), h.position().z
        id_unique, id_counts = np.unique(ids, return_counts=True)
        print dict(zip(id_unique, id_counts))
        unreconstructable = id_unique[id_counts < 4]
        print "unreconstructable: ", unreconstructable

      true_track_ids = []
      numFakes = 0
      num_ignored_secondaries = 0
      currentEta = 0
      currentPt = 0
      for p in store.get(collectionNames["SimParticles"]):
            pm = [ p.core().p4.px, 
                   p.core().p4.py, 
                   p.core().p4.pz, 
                   p.core().p4.mass, ]
            pm = ROOT.TLorentzVector(*pm)
            currentEta = pm.Eta()
            currentPt = pm.Pt()

            if np.abs(pm.Pt()) > args.ptCut and np.abs(pm.Eta()) < args.etaCut:
              if not args.ignore_secondaries or p.core().status != 201: # fcc convention for secondaries
                if not args.checkOnlyPrimary or p.core().bits == 1:
                  true_track_ids.append(p.core().bits)
                  #true_track_ids.append(1)
                  etas[(i, p.core().bits)] = pm.Eta()
                  pts[(i, p.core().bits)] = pm.Pt()
                  cotThetas[(i, p.core().bits)] = 1. / np.tan(pm.Theta())
                  successfully_found_track[(i, p.core().bits)] = 0 # initialization, per default not found
                  delphes_tracking_efficiencies[(i, p.core().bits)] = DelphesTrackingEfficiency(pm)
              else:
                  num_ignored_secondaries +=1
          
      ### print "processing reco tracks from branch 'tracks' " #########
      correctly_identified_tracks_ids = []
      tracks = store.get('tracks')
      numTracks = len(tracks)
      
      for trackCounter, t in enumerate(tracks):
              reco_track_is_fake[(i, trackCounter)] = 1
              reco_track_ids = []
              # go through the hits comprising the track
              for j in range(t.hits_size()):
                # access mc truth for the track Id
                reco_track_ids.append(t.hits(j).core().bits) 
              # if all the hits have the same track id the track was identified correctly
              if len(set(reco_track_ids)) == 1:
                # we might choose to ignore this track for the efficiency nonetheless,
                # (based on pt, eta ... )
                if reco_track_ids[0] in true_track_ids:
                  correctly_identified_tracks_ids.append(reco_track_ids[0])
                  successfully_found_track[(i, reco_track_ids[0])] = 1
                  reconstructed_track_ids[(i, reco_track_ids[0])] = reco_track_ids[0]
                  reco_track_is_fake[(i, trackCounter)] = 0
              else:
                numFakes += 1
              if args.fakeRateHistos:
                ts = t.states(0)
                trackparams = [
                  ts.d0(),
                  ts.z0(),
                  ts.phi(),
                  ts.theta(), #cotan theta
                  ts.qOverP() * -10,
                  ]
                actual_pt = np.abs(1./trackparams[-1]) 
                #if actual_pt > args.ptCut:
                if t.bits() == 1:
                  if np.abs(currentPt) > args.ptCut and np.abs(currentEta) < args.etaCut:
                    reco_track_pt[(i, trackCounter)] = actual_pt
                    #actual_theta = np.arctan( 1. / ts.theta())
                    reco_track_eta[(i, trackCounter)] = ts.theta() #-1* np.log(np.tan(actual_theta * 0.5)) 
                
      # safeguard against empty lists
      if numTracks == 0:
        print "!!!!!!!!!! no tracks in event! eta, pt: ", i, currentEta, currentPt 
        #efficiencies.append(0)
        #fake_rates.append(0)
        #bad_event.append(i)
      else:
        try:
          # efficiency: proportion of true tracks that were successfully reconstructe
          efficiency = len(set(correctly_identified_tracks_ids)) / float(len(set(true_track_ids)))
        except ZeroDivisionError:
          print "!!!!!!!!!! no true_track_ids in event!"
          efficiency = 0 
        # fake rate: proportion of incorrectly reconstructed tracks
        fake_rate = numFakes / float(numTracks)
        print  " %i tracks in event %i (eta %.2f), eff.: %.3f, fake rate:  %.3f | %i secondary particles ignored" \
             % (len(true_track_ids), i, currentEta,  efficiency, fake_rate, num_ignored_secondaries)
        efficiencies.append(efficiency)
        fake_rates.append(fake_rate)
        if efficiency < 1:
          bad_event.append(i)


sk = sorted(etas.iterkeys())
etas = np.array([etas[k] for k in sk])
cotThetas2 = np.array([cotThetas[k] for k in sk])
pts = np.array([pts[k] for k in sk])
skr_pt = sorted(reco_track_pt.iterkeys())
skr_eta = sorted(reco_track_eta.iterkeys())
reco_etas = np.array([reco_track_eta[k] for k in skr_eta])
reco_pts = np.array([reco_track_pt[k] for k in skr_pt])
reco_fr = np.array([reco_track_is_fake[k] for k in skr_pt])
efficiencies = np.array([successfully_found_track.get(k, 0) for k in successfully_found_track.keys()])
delphes_tracking_efficiencies = np.array([delphes_tracking_efficiencies.get(k, 0) for k in sk])

print "Delphes parametrized efficiency: ", np.mean(delphes_tracking_efficiencies)

plt.figure("fit_validation_mctruth")
plt.plot(cotThetas2, np.abs(pts), 'o')
plt.yscale("log")
plt.ylim(1e-1, 1e10)
plt.figure("fit_validation_reco")
plt.yscale("log")
plt.ylim(1e-1, 1e10)
plt.plot(reco_etas, np.abs(reco_pts), 'o')




def efficiency_plot(x, y, conf, ax=plt.gca()):
  means_result = scipy.stats.binned_statistic(x, [y, y**2], 
      bins=conf["numBins"], 
      range=conf["histo_range"], 
      statistic='mean')
  count_result = scipy.stats.binned_statistic(x, [y, y**2], 
      bins=conf["numBins"], 
      range=conf["histo_range"], 
      statistic='count')
  counts, counts2 = count_result.statistic
  means, means2 = means_result.statistic
  standard_deviations = np.sqrt(means2 - means**2)
  bin_edges = means_result.bin_edges
  bin_centers = (bin_edges[:-1] + bin_edges[1:])/2.
  bin_width = 0.45 * (bin_edges[1:] - bin_edges[0:-1])
  means = np.array(means)
  counts[counts == 0] = 0.1
  # clip errorbars
  standard_deviations_high = 0.5 *   np.minimum(np.sqrt(counts) / counts, 1. - means )
  standard_deviations_low =  0.5 * np.minimum(np.sqrt(counts) / counts ,  means)
  if not conf.get("doYErr", 1):
    standard_deviations_high = [0] * len(standard_deviations_high)
    standard_deviations_low = [0] * len(standard_deviations_low)

  # do the plot
  (_, caps, _) = plt.errorbar(x=bin_centers, y=means, 
        xerr=bin_width, yerr=[standard_deviations_low, standard_deviations_high], 
        linestyle = '', marker='o',capsize=2, capthick=1)

  for cap in caps:
    cap.set_markeredgewidth(1)
  plt.xlabel(conf["xlabel"], horizontalalignment='right', x=1.0)
  plt.ylabel(conf["ylabel"], horizontalalignment='right', y=.95, #rotation="horizontal", 
    labelpad=3)
  plt.text(args.textBoxPositionX, args.textBoxPositionY, conf["boxtext"], {'size': 24}, transform=plt.gca().transAxes, bbox={"boxstyle": "round", "alpha": 0.9, "fc": "white"})
  plt.gca().fill_between(bin_centers, 0, means, color=conf.get("fillbetweencolor", "blue"), alpha=0.2, step="mid")
  if conf["setlogx"]:
    plt.xscale("log")

  plt.xlim(conf["xlimlow"], conf["xlimhigh"])
  plt.ylim(-0.05,1.05)


  #for f in conf["filenames"]: # allow different formats to be saved
  #  print "Saving ", f
  #  plt.savefig(f)

# plot efficiencies and fake rates split in pt / eta regimes
if args.histos:
  plt.figure("eff_eta")
  # do some preselection
  eta_selection = pts > args.ptCut
  x = etas[eta_selection]
  y = efficiencies[eta_selection]
  eta_histo_conf = {}
  eta_histo_conf["numBins"] = 70
  eta_histo_conf["histo_range"] = (-6, 6)
  eta_histo_conf["xlabel"] = r"$\eta$"
  eta_histo_conf['ylabel'] = "seeding efficiency"
  eta_histo_conf["filenames"] = [args.filename.replace(".root", "") + "_eff_eta.png"]
  eta_histo_conf["setlogx"] = False
  eta_histo_conf["xlimlow"] = -6
  eta_histo_conf["xlimhigh"] = 6
  eta_histo_conf["boxtext"] = boxtext
  efficiency_plot(x, y, eta_histo_conf)

  plt.figure("eff_pt")
  # do some preselection
  pt_selection = np.abs(etas) < args.etaCut
  x = pts[pt_selection]
  y = efficiencies[pt_selection]
  pt_histo_conf = {}
  pt_histo_conf["numBins"] = 200 #np.logspace(-1,6,35)
  pt_histo_conf["histo_range"] = (0, 20)
  pt_histo_conf["xlabel"] = r"$p_T$ [GeV]"
  pt_histo_conf['ylabel'] = "seeding efficiency"
  pt_histo_conf["filenames"] = [args.filename.replace(".root", "") + "_eff_pt.png"]
  pt_histo_conf["setlogx"] = False
  pt_histo_conf["xlimlow"] = 0.05
  pt_histo_conf["xlimhigh"] = 20#10**5
  pt_histo_conf["boxtext"] = boxtext
  efficiency_plot(x, y, pt_histo_conf)

  plt.figure("delphes eff_pt")
  # do some preselection
  pt_selection = np.abs(etas) < args.etaCut
  x = pts[pt_selection]
  y = delphes_tracking_efficiencies[pt_selection]
  pt_histo_conf = {}
  pt_histo_conf["numBins"] = np.logspace(-1,6,35)
  pt_histo_conf["histo_range"] = None
  pt_histo_conf["xlabel"] = r"$p_T$ [GeV]"
  pt_histo_conf['ylabel'] = "delphes parametrized efficiency"
  pt_histo_conf["filenames"] = [args.filename.replace(".root", "") + "_delpheseff_pt.png"]
  pt_histo_conf["setlogx"] = True
  pt_histo_conf["doYErr"] = False
  pt_histo_conf["fillbetweencolor"] = "green"
  pt_histo_conf["xlimlow"] = 0.05
  pt_histo_conf["xlimhigh"] = 10**5
  pt_histo_conf["boxtext"] = boxtext
  efficiency_plot(x, y, pt_histo_conf)

if args.fakeRateHistos:
  plt.figure("fakeRate_pt")
  pt_selection = np.abs(reco_pts) > args.ptCut
  x = reco_pts[pt_selection]
  y = reco_fr[pt_selection]
  pt_histo_conf = {}
  pt_histo_conf["numBins"] = np.logspace(-1,6,35)
  pt_histo_conf["histo_range"] = None
  pt_histo_conf["xlabel"] = r"$p_T$ [GeV]"
  pt_histo_conf['ylabel'] = "fake rate"
  pt_histo_conf["setlogx"] = True
  pt_histo_conf["xlimlow"] = 0.05
  pt_histo_conf["fillbetweencolor"] = "red"
  pt_histo_conf["xlimhigh"] = 10**5
  pt_histo_conf["boxtext"] = boxtext
  efficiency_plot(x, y, pt_histo_conf)

  plt.figure("fakeRate_eta")
  pt_selection = np.abs(reco_pts) > 1
  x = reco_etas[pt_selection]
  y = reco_fr[pt_selection]
  eta_histo_conf = {}
  eta_histo_conf["numBins"] = 70
  eta_histo_conf["histo_range"] = (-6, 6)
  eta_histo_conf["xlabel"] = r"$\eta$"
  eta_histo_conf['ylabel'] = "fake rate"
  eta_histo_conf["setlogx"] = False
  eta_histo_conf["xlimlow"] = 0.05
  eta_histo_conf["fillbetweencolor"] = "red"
  eta_histo_conf["xlimlow"] = -6
  eta_histo_conf["xlimhigh"] = 6
  eta_histo_conf["boxtext"] = boxtext 
  efficiency_plot(x, y, eta_histo_conf)

print " ####################################### "
print "mean efficiency \t mean fakerate "
print "\t%.3f \t\t\t %.3f" % (np.mean(efficiencies), np.mean(fake_rates))
print " ####################################### "



f = "SeedingEff_" + args.filename.replace(".root", "") + "_etaCut" +str(args.etaCut) + "_ptCut" + str(args.ptCut) + "_events%i" % args.nevents

figs = [plt.figure(n) for n in plt.get_fignums()]
for fig in figs:
  fig.savefig(os.path.join(args.plotprefix, f + "_" + fig._label + ".png"))
print "bad events: ", bad_event
