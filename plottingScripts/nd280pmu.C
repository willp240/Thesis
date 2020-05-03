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


void nd280pmu(std::string file0) {

  gStyle->SetPalette(51);

  TCanvas *c = new TCanvas("canv", "canv", 1580, 1080);
  TFile *file_ = TFile::Open(file0.c_str());
 
  c->SetTopMargin(0.05);
  c->SetBottomMargin(0.14);
  c->SetLeftMargin(0.1);
  c->SetRightMargin(0.16);
  
  std::string samp[7] = {"cc0pi0p", "cc0pi1p", "cc0piNp", "cc1pi0p", "cc1pi1p", "cc1piNp", "ccOther"};

  for(int s=0; s<7; s++){

    std::string name = samp[s];

    TH2D *hist = (TH2D*) file_->Get(name.c_str())->Clone();

    hist->SetTitle("");
    hist->GetYaxis()->SetTitleOffset(0.85);
    hist->GetXaxis()->SetTitleOffset(1.0);
    hist->GetZaxis()->SetTitleOffset(0.9);
    hist->GetXaxis()->SetTitleSize(0.05);
    hist->GetXaxis()->SetRangeUser(0,5000);
    hist->GetZaxis()->SetRangeUser(-0.1,hist->GetMaximum());
    hist->GetYaxis()->SetTitleSize(0.05);
    hist->GetZaxis()->SetTitleSize(0.05);
    hist->GetXaxis()->SetTitle("p_{#mu} (MeV)");
    hist->GetYaxis()->SetTitle("cos#theta_{#mu}");
    hist->GetZaxis()->SetTitle("Events");
    palette = (TPaletteAxis*) hist->GetListOfFunctions()->FindObject("palette");
    palette->SetX1NDC(0.85);
    palette->SetX2NDC(0.90);
    palette->SetY1NDC(0.14);
    palette->SetY2NDC(0.95);
    hist->Draw("colz");
    c->Print((std::string("nd280_pmtmuu_")+name+std::string(".pdf")).c_str());
    std::cout << "got " << name << std::endl;
  }
}
