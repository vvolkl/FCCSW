
void drawtext() {

   TPaveText *pt = new TPaveText(0.1592498,0.8973723,0.4499682,0.8913649,"NDC");
   pt->SetName("title______dfsdf");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->SetFillStyle(0);
   pt->SetTextSize(0.04);
   pt->SetTextAlign(10);
   TText *pt_LaTex = pt->AddText("#it{FCC-ee} Simulation");
   //TText *pt_LaTex3 = pt->AddText("#bf{Tracker layout v3.03}");
   //pt_LaTex3->SetTextSize(0.03);
   //TText *pt_LaTex4 = pt->AddText("#bf{Single muons}");
   //pt_LaTex4->SetTextSize(0.03);
   pt->Draw();

}


void plot() {
  {
  TFile f= TFile("fccee-lar-materialscan_lar.root", "OPEN");
  f.Print();
  f.ls();
  TH1F* h = (TH1F*)  f.Get("nX0");
  std::cout << h->GetName() << std::endl;
  TCanvas* C = new TCanvas;
  h->SetLineColor(kBlue );
  h->SetLineColor(kBlack );
  h->SetLineWidth(1);
  h->SetFillColor(kBlue );
  h->Draw();
  drawtext();
  C->Print("fccee-lar-materialscan_lar_nX0.png");
  C->Print("fccee-lar-materialscan_lar_nX0.pdf");
  }

  {
  TFile f= TFile("fccee-lar-materialscan_idea-tilecal.root", "OPEN");
  f.Print();
  f.ls();
  TH1F* h = (TH1F*)  f.Get("nLambda");
  std::cout << h->GetName() << std::endl;
  TCanvas* C = new TCanvas;
  h->SetLineColor(kOrange );
  h->SetLineColor(kBlack );
  h->SetLineWidth(1 );
  h->SetFillColor(kOrange );
  h->Draw();
  drawtext();
  C->Print("fccee-lar-materialscan_tilecal_nLambda.png");
  C->Print("fccee-lar-materialscan_tilecal_nLambda.pdf");
  }

  {
  TFile f= TFile("fccee-lar-materialscan_idea.root", "OPEN");
  f.Print();
  f.ls();
  TH1F* h = (TH1F*)  f.Get("nX0");
  std::cout << h->GetName() << std::endl;
  TCanvas* C = new TCanvas;
  h->SetLineColor(kBlack );
  h->SetLineWidth(1 );
  h->SetFillColor(kGray );
  h->Draw();

  drawtext();
  C->Print("fccee-lar-materialscan_idea_nX0.png");
  C->Print("fccee-lar-materialscan_idea_nX0.pdf");



  }




gApplication->Terminate();
}
