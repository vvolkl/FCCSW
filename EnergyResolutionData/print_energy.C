void print_energy(std::string filename) {

  TFile f = TFile(filename.c_str());
  TTree* events = (TTree*) f.Get("events");
  std::vector<fcc::MCParticleData>* p = new std::vector<fcc::MCParticleData>();
  events->SetBranchAddress("GenParticles", &p);
  events->GetEvent(1);
  fcc::MCParticleData part = (*p)[0];
  std::cout << events->GetEntries() << "\t" << std::sqrt(pow(part.core.p4.px,2) + pow(part.core.p4.py,2)+pow(part.core.p4.pz,2)) << std::endl;
  gApplication->Terminate();


}
