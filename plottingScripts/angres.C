#include "TH1D.h"
#include "TH1F.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TFile.h"
#include "TKey.h"
#include "THStack.h"
#include "TROOT.h"
#include "TGaxis.h"
#include <iostream>
#include <algorithm>


void angres(std::string file0) {

  TCanvas *c = new TCanvas("canv", "canv", 1580, 1080);
  TFile *file_ = TFile::Open(file0.c_str());

  std::string kin = "angres";
 
  c->SetTopMargin(0.05);
  c->SetBottomMargin(0.12);
  c->SetLeftMargin(0.09);
  c->SetRightMargin(0.16);
  
  TCanvas *c1 = file_->Get("Canvas_1")->Clone();
  TH1* test;
  TH2* h2;
  h2 = (TH2*)c1->GetPrimitive("h1");
  c1->Clear();
  c1->Draw();
  gStyle->SetPalette(51);
  c1->SetRightMargin(0.18);
  c1->SetLeftMargin(0.16);
  c1->SetBottomMargin(0.15);
  h2->GetXaxis()->SetTitle("Reconstructed cos #theta_{#mu}");
  h2->GetXaxis()->SetTitleSize(0.04);
  h2->GetYaxis()->SetTitleSize(0.04);
  h2->GetZaxis()->SetTitleSize(0.04);
  h2->GetXaxis()->SetTitleOffset(1.3);
  h2->GetYaxis()->SetTitleOffset(1.7);
  h2->GetZaxis()->SetTitleOffset(1.5);
  h2->GetXaxis()->SetLabelSize(0.03);
  h2->GetYaxis()->SetLabelSize(0.03);
  h2->GetZaxis()->SetTitle("Events");
  h2->SetTitle("");
  h2->GetYaxis()->SetTitle("True cos #theta_{#mu}");
  h2->GetXaxis()->SetNdivisions(1005);
  palette = (TPaletteAxis*) h2->GetListOfFunctions()->FindObject("palette");
  palette->SetX1NDC(0.83);
  palette->SetX2NDC(0.88);
  h2->GetXaxis()->SetLabelSize(0.03);
  h2->GetZaxis()->SetLabelSize(0.03);
  h2->Draw("colz");
  int nbins = h2->GetXaxis()->GetNbins();
  c1->Print((kin+std::string("2d.pdf")).c_str());
  
  TH1 *hbins[500];
  for (int i=0;i<500;i++){
    hbins[i] = h2->ProjectionX(Form("bin%d",i+1),i+1,i+2); 
  }

  TH1 *hmean;
  hmean = new TH1F("","",500,0.95,1.0);
 
  for(int i=0; i<500; i++){
    TString name = Form("binY%d",i+1);
    hbins[i]->Fit("gaus");
    hmean->SetBinContent(i+1,hbins[i]->GetRMS());//GetFunction("gaus")->GetParameter(2));
  }

  //hbins[450]->Draw();
  c1->SetCanvasSize(1580,1080);
  c1->SetLeftMargin(0.14);
  c1->SetRightMargin(0.1);
  //  hmean->GetXaxis()->SetRangeUser(0.96,0.99);
  hmean->GetXaxis()->SetTitle("True cos #theta_{#mu}");
  ///hmean->GetYaxis()->SetRangeUser(0,160);
  hmean->GetYaxis()->SetTitleOffset(1.3);
  hmean->GetXaxis()->SetTitleOffset(1.2);
  hmean->GetYaxis()->SetTitle("RMS");
  hmean->Fit("pol1","","",0.96,0.99);
  hmean->GetFunction("pol1")->SetLineColor(kRed);
  hmean->Draw("colz");
  c1->Print((kin+std::string("1d.pdf")).c_str());
  std::cout << nbins << std::endl;
}
