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

all: seeding_plots 
seeding_plots: tracks
	./run python Reconstruction/RecTracker/scripts/plotTrackSeeds.py muons_for_seeding_tracks.root
# run the track seeding 
tracks: muons_for_seeding_tracks.root
muons_for_seeding_tracks.root: muons
	./run fccrun.py Reconstruction/RecTracker/options/TrickTrackReco.py  --trajectories --recoHelix

# simulate single muons to have a sample for track seeding
muons: muons_for_seeding.log
muons_for_seeding.log:
	./run fccrun.py  Reconstruction/RecTracker/options/geantSim_TrackerPerformance.py -N 1000 -s 0123456 --outName muons_for_seeding.root --singlePart --particle 13 --etaMin 0 --etaMax 6 --discretePt --pt 1000 2000 5000 10000 100000 1000000 --trajectories | tee muons_for_seeding.log


# use truth seeding + RiemannFit to calculate single particle resolutions from the muon sample
single_particle_resolutions: single_particle_resolutions.log
single_particle_resolutions.log:
	./run fccrun.py Reconstruction/RecTracker/options/single_particle_trackFits.py | tee single_particle_resolutions.log
	
