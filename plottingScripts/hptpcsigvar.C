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


void hptpcsigvar() {

  gStyle->SetPalette(51);

  TCanvas *c = new TCanvas("canv", "canv", 1580, 1080);
  TFile *filemu = TFile::Open("inputFiles/STV/RFG_PmuTmu.Spectra.root");
  TFile *filep = TFile::Open("inputFiles/STV/RFGMAQE_PpTmu.Spectra.root"); 
  TFile *filepi = TFile::Open("inputFiles/STV/HPTPC_Ppion_costheta_NormFakeData.root");
  TFile *filept = TFile::Open("inputFiles/STV/ND280_MixedSelections.Spectra.root");

  c->SetTopMargin(0.05);
  c->SetBottomMargin(0.14);
  c->SetLeftMargin(0.1);
  c->SetRightMargin(0.16);
  
  std::string samp[7] = {"cc0pi0p", "cc0pi1p", "cc0piNp", "nu_cc1pi0p_1", "cc1pi1p", "cc1piNp", "ccOther"};

  std::string xtitle[7] = {"p_{#mu} (MeV)", "p_{p} (MeV)", "p_{p} (MeV)", "p_{#pi} (MeV)", "#delta#it{p}_{T} (MeV)","p_{p} (MeV)", "p_{#mu} (MeV)"};
  std::string ytitle[7] = {"cos#theta_{#mu}", "cos#theta_{#mu}", "cos#theta_{#mu}", "cos#theta_{#mu}", "cos#theta_{#mu}", "cos#theta_{#mu}", "cos#theta_{#mu}"};

  for(int s=0; s<7; s++){

    std::string name = samp[s];

    if(s==0 || s==6)
      TH2D *hist = (TH2D*) filemu->Get(name.c_str())->Clone();
    if(s==1 || s== 2 || s==5)
      TH2D *hist = (TH2D*) filep->Get(name.c_str())->Clone();
    if(s==3)
      TH2D *hist = (TH2D*) filepi->Get(name.c_str())->Clone();
    if(s==4)
      TH2D *hist = (TH2D*) filept->Get(name.c_str())->Clone();
    
    hist->SetTitle("");
    hist->GetYaxis()->SetTitleOffset(0.85);
    hist->GetXaxis()->SetTitleOffset(1.0);
    hist->GetZaxis()->SetTitleOffset(0.9);
    hist->GetXaxis()->SetTitleSize(0.05);
    hist->GetXaxis()->SetRangeUser(0,5000);
    hist->GetZaxis()->SetRangeUser(-0.1,hist->GetMaximum());
    hist->GetYaxis()->SetTitleSize(0.05);
    hist->GetZaxis()->SetTitleSize(0.05);
    hist->GetXaxis()->SetTitle(xtitle[s].c_str());
    hist->GetYaxis()->SetTitle(ytitle[s].c_str());
    hist->GetZaxis()->SetTitle("Events");
    if(s!=3){
      palette = (TPaletteAxis*) hist->GetListOfFunctions()->FindObject("palette");
      palette->SetX1NDC(0.85);
      palette->SetX2NDC(0.90);
      palette->SetY1NDC(0.14);
      palette->SetY2NDC(0.95);
    }
    hist->Draw("colz");
    //    std::cout << "print to " << "hptpc_sigvar_" << std::endl;
    c->Print((std::string("hptpc_sigvar_")+name+std::string(".pdf")).c_str());
    std::cout << "got " << name << std::endl;
  }
}
