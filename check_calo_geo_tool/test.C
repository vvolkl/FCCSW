
auto magic = fcc::MCParticleData();


auto _get_layerid = [](fcc::BareHit in) {
  int result;
  int systemId = in.cellId % 16;
  int layerId = (in.cellId >> 4) % 32;
  result = systemId * 100 + layerId;
  return result;


};

auto get_layerid = [](std::vector<fcc::PositionedTrackHitData> in) {

  std::vector<int> result;
  for (int i = 0; i < in.size(); ++i) {
  result.push_back(_get_layerid(in[i].core));
  }

return result;
};


int test(std::string filename) {
  gROOT->SetBatch();
  ROOT::RDataFrame df("events", //"geantinos_for_seeding_single_direction_tracks.root");
    filename);
  auto dff = df
               .Define("layerid", get_layerid, {"SortedQuadruplets"})
              ;


  std::string basefilename(filename);
  basefilename.replace(basefilename.end() - 5, basefilename.end(), "");

  std::string snapshotname(filename);
  snapshotname.replace(snapshotname.end() - 5, snapshotname.end(), "_snapshot.root");
  std::cout << "Write snapshot to " <<  snapshotname << " ... " <<  std::endl;
  dff.Snapshot("events", snapshotname, 
    {
      "SimParticles",
      "Gen_Eta",
      "Gen_Pt",
      "GenParticles_Pt", 
      "GenParticles_Eta",
      "TrackRecoParticles_Eta",
      "TrackRecoParticles_Pt",
      "SimParticles_Eta",
      "SimParticles_Pt",
      "CorrectTrackRecoParticles",
      "FakeTrackRecoParticles",
      "PrimaryFound",
      "SimParticleFound",
      "SortedQuadruplets",
      "PassCurvatureCheck",
      "DistanceBeamspot",
      "VertexR",
      "PassAlignmentCheck",
      "RZAlignment",
      "CurvatureRadius",
      "TimeDifferences",
      "ZDifferences",
      "PhiDifferences",
      "PassAllChecks",
      });


  auto sum_curvature_check = dff.Sum("PassCurvatureCheck");
  int sum_curvature_check2 = *sum_curvature_check;
  std::cout << "Events that pass curvature check: " << sum_curvature_check2 << std::endl;

  auto sum_align_check = dff.Sum("PassAlignmentCheck");
  int sum_align_check2 = *sum_align_check;
  std::cout << "Events that pass alignment check: " << sum_align_check2 << std::endl;





  gApplication->Terminate();
  return 0;


}
