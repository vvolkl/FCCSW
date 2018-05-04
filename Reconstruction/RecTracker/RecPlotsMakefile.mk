###########################
# FCC Track Reco Plots Makefile
# run from FCCSW directory
# with 
#
# make -f Reconstruction/RecTracker/RecPlotsMakeFile.mk plots
#
#
#
###########################
.PHONY: all
all: plots

plots: tricktrack_seeding_example.root
	./run python Reconstruction/RecTracker/scripts/plotTrackSeeds.py tricktrack_tracks.root
tracks: muons_for_seeding.root
	./run fccrun.py Reconstruction/RecTracker/options/TrickTrackReco.py 
muons: muons_for_seeding.log

muons_for_seeding.log:
	./run fccrun.py  Reconstruction/RecTracker/options/geantSim_TrackerPerformance.py -N 1000 -s 0123456 --outName muons_for_seeding.root --singlePart --particle 13 --etaMin 0 --etaMax 6 --discretePt --pt 1000 2000 5000 10000 100000 1000000 --trajectories | tee muons_for_seeding.log


single_particle_resolutions: single_particle_resolutions.log

single_particle_resolutions.log:
	./run fccrun.py Reconstruction/RecTracker/options/single_particle_trackFits.py
	
