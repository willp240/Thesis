#include "TH1D.h"
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


void RatioHistos(std::string file1, std::string file2, std::string file3) {

  TCanvas *c = new TCanvas("canv", "canv", 1080, 1080);
  TFile *onefile = TFile::Open(file1.c_str());
  TFile *twofile = TFile::Open(file2.c_str());
  TFile *threefile = TFile::Open(file3.c_str());
  TFile *outfile = new TFile("25OctEbPeakRatios.root", "RECREATE");
 
  c->SetTopMargin(0.07);
  c->SetBottomMargin(0.11);
  c->SetLeftMargin(0.10);
  c->SetRightMargin(0.18);
  
  c->Print((std::string("25OctEbPeakRatios")+std::string(".pdf[")).c_str());

  std::string FGD1Loc,FGD2Loc,h1Name,h2Name;  
  FGD1Loc = "FGD1 numuCC 0pi/FGD1_numuCC_0pi_mean";
  FGD2Loc = "FGD2 numuCC 0pi/FGD2_numuCC_0pi_mean";

  ///////////////////////////////////// Peak 1
  h1Name = "Peak1FGD1";
  h2Name = "Peak1FGD2";
  TH2D *h1FGD1 = (TH2D*)(onefile->Get((FGD1Loc).c_str()));
  if(!h1FGD1)
    std::cout << "No h1FGD1!" << std::endl;
  h1FGD1->SetName(h1Name.c_str());
  h1FGD1->SetTitle("CC0Pi FGD1: 25 MeV < Eb_C < 30 MeV");
  c->Clear();
  h1FGD1->GetXaxis()->SetRangeUser(0,5000);
  h1FGD1->GetXaxis()->SetTitleOffset(1.1);
  h1FGD1->GetYaxis()->SetTitleOffset(1.1);
  h1FGD1->GetZaxis()->SetTitleOffset(1.5);
  h1FGD1->Draw("colz");
  gPad->Update();
  TPaletteAxis* palette = (TPaletteAxis*) h1FGD1->GetListOfFunctions()->FindObject("palette");
  palette->SetX1NDC(0.83);
  palette->SetX2NDC(0.88);
  h1FGD1->GetXaxis()->SetLabelSize(0.03);
  h1FGD1->GetZaxis()->SetLabelSize(0.03);
  c->Print((std::string("25OctEbPeakRatios")+std::string(".pdf")).c_str());
  h1FGD1->Write();

  TH2D *h1FGD2 = (TH2D*)(onefile->Get((FGD2Loc).c_str()));
  if(!h1FGD2)
    std::cout << "No h1FGD2!" << std::endl;
  h1FGD2->SetName(h2Name.c_str());
  h1FGD2->SetTitle("CC0Pi FGD2: 25 MeV < Eb_C < 30 MeV");
  c->Clear();
  h1FGD2->GetXaxis()->SetRangeUser(0,5000);
  h1FGD2->GetXaxis()->SetTitleOffset(1.1);
  h1FGD2->GetYaxis()->SetTitleOffset(1.1);
  h1FGD2->GetZaxis()->SetTitleOffset(1.5);
  h1FGD2->Draw("colz");
  gPad->Update();
  palette = (TPaletteAxis*) h1FGD2->GetListOfFunctions()->FindObject("palette");
  palette->SetX1NDC(0.83);
  palette->SetX2NDC(0.88);
  h1FGD2->GetXaxis()->SetLabelSize(0.03);
  h1FGD2->GetZaxis()->SetLabelSize(0.03);
  c->Print((std::string("25OctEbPeakRatios")+std::string(".pdf")).c_str());
  h1FGD2->Write();

  ///////////////////////////////////// Peak 2
  h1Name = "Peak2FGD1";
  h2Name = "Peak2FGD2";
  TH2D *h2FGD1 = (TH2D*)(twofile->Get((FGD1Loc).c_str()));
  if(!h2FGD1)
    std::cout << "No h2FGD1!" << std::endl;
  h2FGD1->SetName(h2Name.c_str());
  h2FGD1->SetTitle("CC0Pi FGD1: 30 MeV < Eb_C < 33 MeV");
  c->Clear();
  h2FGD1->GetXaxis()->SetRangeUser(0,5000);
  h2FGD1->GetXaxis()->SetTitleOffset(1.1);
  h2FGD1->GetYaxis()->SetTitleOffset(1.1);
  h2FGD1->GetZaxis()->SetTitleOffset(1.5);
  h2FGD1->Draw("colz");
  gPad->Update();
  palette = (TPaletteAxis*) h2FGD1->GetListOfFunctions()->FindObject("palette");
  palette->SetX1NDC(0.83);
  palette->SetX2NDC(0.88);
  h2FGD1->GetXaxis()->SetLabelSize(0.03);
  h2FGD1->GetZaxis()->SetLabelSize(0.03);
  c->Print((std::string("25OctEbPeakRatios")+std::string(".pdf")).c_str());
  h2FGD1->Write();

  TH2D *h2FGD2 = (TH2D*)(twofile->Get((FGD2Loc).c_str()));
  if(!h2FGD2)
    std::cout << "No h2FGD2!" << std::endl;
  h2FGD2->SetName(h2Name.c_str());
  h2FGD2->SetTitle("CC0Pi FGD2: 30 MeV < Eb_C < 33 MeV");
  c->Clear();
  h2FGD2->GetXaxis()->SetRangeUser(0,5000);
  h2FGD2->GetXaxis()->SetTitleOffset(1.1);
  h2FGD2->GetYaxis()->SetTitleOffset(1.1);
  h2FGD2->GetZaxis()->SetTitleOffset(1.5);
  h2FGD2->Draw("colz");
  gPad->Update();
  palette = (TPaletteAxis*) h2FGD2->GetListOfFunctions()->FindObject("palette");
  palette->SetX1NDC(0.83);
  palette->SetX2NDC(0.88);
  h2FGD2->GetXaxis()->SetLabelSize(0.03);
  h2FGD2->GetZaxis()->SetLabelSize(0.03);
  c->Print((std::string("25OctEbPeakRatios")+std::string(".pdf")).c_str());
  h2FGD2->Write();
	
  h1Name = "Peak2FGD1Ratio";
  h2Name = "Peak2FGD2Ratio";
  TH2D *h2FGD1Ratio = (TH2D*)(twofile->Get((FGD1Loc).c_str()));
  if(!h2FGD1Ratio)
    std::cout << "No h2FGD1Ratio!" << std::endl;
  h2FGD1Ratio->SetName(h1Name.c_str());
  h2FGD1Ratio->SetTitle("CC0Pi FGD1 Ratio (25 MeV < Eb_C < 30 MeV):(30 MeV < Eb_C < 33 MeV)");
  h2FGD1Ratio->Divide(h1FGD1);
  c->Clear();
  h2FGD1Ratio->GetXaxis()->SetRangeUser(0,5000);
  h2FGD1Ratio->GetXaxis()->SetTitleOffset(1.1);
  h2FGD1Ratio->GetYaxis()->SetTitleOffset(1.1);
  h2FGD1Ratio->GetZaxis()->SetTitleOffset(1.5);
  h2FGD1Ratio->GetZaxis()->SetRangeUser(0.9,1.1);
  h2FGD1Ratio->Draw("colz");
  gPad->Update();
  palette = (TPaletteAxis*) h2FGD1Ratio->GetListOfFunctions()->FindObject("palette");
  palette->SetX1NDC(0.83);
  palette->SetX2NDC(0.88);
  h2FGD1Ratio->GetXaxis()->SetLabelSize(0.03);
  h2FGD1Ratio->GetZaxis()->SetLabelSize(0.03);
  c->Print((std::string("25OctEbPeakRatios")+std::string(".pdf")).c_str());
  h2FGD1Ratio->Write();

  TH2D *h2FGD2Ratio = (TH2D*)(twofile->Get((FGD2Loc).c_str()));
  if(!h2FGD2Ratio)
    std::cout << "No h2FGD2Ratio!" << std::endl;
  h2FGD2Ratio->SetName(h2Name.c_str());
  h2FGD2Ratio->SetTitle("CC0Pi FGD2 Ratio (25 MeV < Eb_C < 30 MeV):(30 MeV < Eb_C < 33 MeV)");
  h2FGD2Ratio->Divide(h1FGD2);
  c->Clear();
  h2FGD2Ratio->GetXaxis()->SetRangeUser(0,5000);
  h2FGD2Ratio->GetXaxis()->SetTitleOffset(1.1);
  h2FGD2Ratio->GetYaxis()->SetTitleOffset(1.1);
  h2FGD2Ratio->GetZaxis()->SetTitleOffset(1.5);
  h2FGD2Ratio->GetZaxis()->SetRangeUser(0.9,1.1);
  h2FGD2Ratio->Draw("colz");
  gPad->Update();
  palette = (TPaletteAxis*) h2FGD2Ratio->GetListOfFunctions()->FindObject("palette");
  palette->SetX1NDC(0.83);
  palette->SetX2NDC(0.88);
  h2FGD2Ratio->GetXaxis()->SetLabelSize(0.03);
  h2FGD2Ratio->GetZaxis()->SetLabelSize(0.03);
  c->Print((std::string("25OctEbPeakRatios")+std::string(".pdf")).c_str());
  h2FGD2Ratio->Write();

  ///////////////////////////////////// Peak 3
  h1Name = "Peak3FGD1";
  h2Name = "Peak3FGD2";
  TH2D *h3FGD1 = (TH2D*)(threefile->Get((FGD1Loc).c_str()));
  if(!h3FGD1)
    std::cout << "No h3FGD1!" << std::endl;
  h3FGD1->SetName(h1Name.c_str());
  h3FGD1->SetTitle("CC0Pi FGD1: 30 MeV < Eb_C < 33 MeV");
  c->Clear();
  h3FGD1->GetXaxis()->SetRangeUser(0,5000);
  h3FGD1->GetXaxis()->SetTitleOffset(1.1);
  h3FGD1->GetYaxis()->SetTitleOffset(1.1);
  h3FGD1->GetZaxis()->SetTitleOffset(1.5);
  h3FGD1->Draw("colz");
  gPad->Update();
  palette = (TPaletteAxis*) h3FGD1->GetListOfFunctions()->FindObject("palette");
  palette->SetX1NDC(0.83);
  palette->SetX2NDC(0.88);
  h3FGD1->GetXaxis()->SetLabelSize(0.03);
  h3FGD1->GetZaxis()->SetLabelSize(0.03);
  c->Print((std::string("25OctEbPeakRatios")+std::string(".pdf")).c_str());
  h3FGD1->Write();

  TH2D *h3FGD2 = (TH2D*)(threefile->Get((FGD2Loc).c_str()));
  if(!h3FGD2)
    std::cout << "No h3FGD2!" << std::endl;
  h3FGD2->SetName(h2Name.c_str());
  h3FGD2->SetTitle("CC0Pi FGD2: 30 MeV < Eb_C < 33 MeV");
  c->Clear();
  h3FGD2->GetXaxis()->SetRangeUser(0,5000);
  h3FGD2->GetXaxis()->SetTitleOffset(1.1);
  h3FGD2->GetYaxis()->SetTitleOffset(1.1);
  h3FGD2->GetZaxis()->SetTitleOffset(1.5);
  h3FGD2->Draw("colz");
  gPad->Update();
  palette = (TPaletteAxis*) h3FGD2->GetListOfFunctions()->FindObject("palette");
  palette->SetX1NDC(0.83);
  palette->SetX2NDC(0.88);
  h3FGD2->GetXaxis()->SetLabelSize(0.03);
  h3FGD2->GetZaxis()->SetLabelSize(0.03);
  c->Print((std::string("25OctEbPeakRatios")+std::string(".pdf")).c_str());
  h3FGD2->Write();
	
  h1Name = "Peak3FGD1Ratio";
  h2Name = "Peak3FGD2Ratio";
  TH2D *h3FGD1Ratio = (TH2D*)(threefile->Get((FGD1Loc).c_str()));
  if(!h3FGD1Ratio)
    std::cout << "No h3FGD1Ratio!" << std::endl;
  h3FGD1Ratio->SetName(h1Name.c_str());
  h3FGD1Ratio->SetTitle("CC0Pi FGD1 Ratio (25 MeV < Eb_C < 30 MeV):(33 MeV < Eb_C < 37 MeV)");
  h3FGD1Ratio->Divide(h1FGD1);
  c->Clear();
  h3FGD1Ratio->GetXaxis()->SetRangeUser(0,5000);
  h3FGD1Ratio->GetXaxis()->SetTitleOffset(1.1);
  h3FGD1Ratio->GetYaxis()->SetTitleOffset(1.1);
  h3FGD1Ratio->GetZaxis()->SetTitleOffset(1.5);
  h3FGD1Ratio->GetZaxis()->SetRangeUser(0.9,1.1);
  h3FGD1Ratio->Draw("colz");
  gPad->Update();
  palette = (TPaletteAxis*) h3FGD1Ratio->GetListOfFunctions()->FindObject("palette");
  palette->SetX1NDC(0.83);
  palette->SetX2NDC(0.88);
  h3FGD1Ratio->GetXaxis()->SetLabelSize(0.03);
  h3FGD1Ratio->GetZaxis()->SetLabelSize(0.03);
  c->Print((std::string("25OctEbPeakRatios")+std::string(".pdf")).c_str());
  h3FGD1Ratio->Write();

  TH2D *h3FGD2Ratio = (TH2D*)(threefile->Get((FGD2Loc).c_str()));
  if(!h3FGD2Ratio)
    std::cout << "No h3FGD2Ratio!" << std::endl;
  h3FGD2Ratio->SetName(h2Name.c_str());
  h3FGD2Ratio->SetTitle("CC0Pi FGD2 Ratio (25 MeV < Eb_C < 30 MeV):(33 MeV < Eb_C < 37 MeV)");
  h3FGD2Ratio->Divide(h1FGD2);
  c->Clear();
  h3FGD2Ratio->GetXaxis()->SetRangeUser(0,5000);
  h3FGD2Ratio->GetXaxis()->SetTitleOffset(1.1);
  h3FGD2Ratio->GetYaxis()->SetTitleOffset(1.1);
  h3FGD2Ratio->GetZaxis()->SetTitleOffset(1.5);
  h3FGD2Ratio->GetZaxis()->SetRangeUser(0.9,1.1);
  h3FGD2Ratio->Draw("colz");
  gPad->Update();
  palette = (TPaletteAxis*) h3FGD2Ratio->GetListOfFunctions()->FindObject("palette");
  palette->SetX1NDC(0.83);
  palette->SetX2NDC(0.88);
  h3FGD2Ratio->GetXaxis()->SetLabelSize(0.03);
  h3FGD2Ratio->GetZaxis()->SetLabelSize(0.03);
  c->Print((std::string("25OctEbPeakRatios")+std::string(".pdf")).c_str());
  h3FGD2Ratio->Write();
	  
  std::cout << "**********"<< std::endl;
   
  c->Print((std::string("25OctEbPeakRatios")+std::string(".pdf]")).c_str());
}
