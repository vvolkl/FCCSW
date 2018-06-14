

from EventStore import EventStore
import argparse
import ROOT
from ROOT import fcc, std
import numpy as np


# command line arguments
parser = argparse.ArgumentParser()
parser.add_argument("filename", help="edm file to read")
parser.add_argument("--nevents", help="max events to process (takes length of input root file by default)", type=int, default=-1)
parser.add_argument('--output', type=str, help="name of rootfile to write", default="ExtrapolationCheckPlots.root")
args = parser.parse_args()


print "creating root file and trees ..."
f = ROOT.TFile(args.filename, "open")
outfile = ROOT.TFile(args.output, "recreate")
t = f.Get("events")

ROOT.gROOT.SetBatch()
outfile.cd()
    

c1 = ROOT.TCanvas()


histos_for_ratio = []
titles  = { 
            "ExtrapolatedTrackstates.referencePoint": "ACTS Extrapolation",
            "SimTrackerPositionedHits.position": "Geant4 Geantino Hits",
          }
for _htype in ["ExtrapolatedTrackstates.referencePoint", "SimTrackerPositionedHits.position"]:
  t.Draw(_htype+".x:"+_htype+".y")
  graph_xy = ROOT.TGraph(t.GetSelectedRows(), t.GetV2(), t.GetV1());
  graph_xy.SetMarkerStyle(4)
  graph_xy.SetTitle(titles[_htype] + " xy")
  graph_xy.GetXaxis().SetTitle("X [mm]")
  graph_xy.GetYaxis().SetTitle("Y [mm]")


  graph_xy.SetLineColorAlpha(ROOT.kWhite, 1.0);
  graph_xy.SetMarkerColor(ROOT.kBlue);
  graph_xy.Write(graph_xy.GetTitle())

  t.Draw("sqrt(pow("+_htype+".x,2)+pow("+_htype+".y,2)):"+_htype+".z")
  graph_rz = ROOT.TGraph(t.GetSelectedRows(), t.GetV2(), t.GetV1());
  graph_rz.SetMarkerStyle(4)
  graph_rz.SetTitle(titles[_htype] + " rz")
  graph_rz.GetXaxis().SetTitle("Z [mm]")
  graph_rz.GetYaxis().SetTitle("R [mm]")
  graph_rz.SetLineColorAlpha(ROOT.kWhite, 1.0);
  graph_rz.SetMarkerColor(ROOT.kBlue);
  graph_rz.Write(graph_rz.GetTitle())



  
  h2 = ROOT.TH2F(_htype+"numSize2D", _htype+"numSize2D", 100, 0., 6., 100, 0., 29.5)
  h2.GetXaxis().SetTitle("#eta")
  h2.GetYaxis().SetTitle("numHits")
  t.Project(_htype+"numSize2D", "@"+_htype.split(".")[0]+".size():atanh(allGenParticles.core.p4.pz / sqrt(pow(allGenParticles.core.p4.px,2) + pow(allGenParticles.core.p4.py,2) + pow(allGenParticles.core.p4.pz,2)))")

  pf1 = h2.ProfileX(_htype+"Profile")
  pf1.GetXaxis().SetTitle("#eta")
  pf1.GetYaxis().SetTitle("numHits")

  histos_for_ratio.append(pf1)
  pf1.Write()
  h2.Write()

if len(histos_for_ratio) > 1:
  hr = histos_for_ratio[0].Clone("ACTS / Geant4 extrapolation ratio")
  hr.Divide(histos_for_ratio[1])
  hr.Draw("HIST P0")


  hr.Write()


print "... writing root file " + args.output
