#include <TH2Poly.h>
#include <TH2D.h>
#include <vector>
#include <TMath.h>
#include <TFile.h>
#include <TTree.h>
#include <THStack.h>
#include <iostream> 

void plotXSCorr(std::string file1){
   
  TFile* f1 = TFile::Open(file1.c_str());
  std::string outname = file1+"_cov";
  TFile * fileOut = new TFile((outname+".root").c_str(),"RECREATE");
  
  TCanvas *c = new TCanvas("canv", "canv", 1280, 1080);
  c->SetTopMargin(0.07);
  c->SetBottomMargin(0.27);
  c->SetLeftMargin(0.23);
  c->SetRightMargin(0.18);

  TH2D* h1 = f1->Get("postfit_corr_plot");
  
  h1->GetXaxis()->SetRangeUser(100,144);
  h1->GetYaxis()->SetRangeUser(100,144);

  h1->GetXaxis()->SetLabelSize(0.03);
  h1->GetYaxis()->SetLabelSize(0.03);
  h1->GetXaxis()->LabelsOption("v");
  h1->SetTitle("");

  TPaletteAxis* palette = (TPaletteAxis*) h1->GetListOfFunctions()->FindObject("palette");
  palette->SetX1NDC(0.83);
  palette->SetX2NDC(0.88);
  palette->SetY1NDC(0.27);
  palette->SetY2NDC(0.93);

  h1->Draw("colz");
  c->Print("5q2xs.pdf");
}

