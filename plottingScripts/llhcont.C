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
#include "TH2Poly.h"

void llhcont(std::string file0) {

  //  gStyle->SetPalette(56);

  TCanvas *c = new TCanvas("canv", "canv", 1580, 1080);
  TFile *file_ = TFile::Open(file0.c_str());
 
  c->SetTopMargin(0.05);
  c->SetBottomMargin(0.14);
  c->SetLeftMargin(0.1);
  c->SetRightMargin(0.16);
  
  std::string samp[18] = {"FGD1 numuCC 0pi", "FGD1 numuCC 1pi", "FGD1 numuCC other", "FGD1 anti-numuCC 0pi", "FGD1 anti-numuCC 1pi", "FGD1 anti-numuCC other", "FGD1 NuMuBkg CC0pi in AntiNu Mode", "FGD1 NuMuBkg CC1pi in AntiNu Mode", "FGD1 NuMuBkg CCother in AntiNu Mode", "FGD2 numuCC 0pi", "FGD2 numuCC 1pi", "FGD2 numuCC other", "FGD2 anti-numuCC 0pi","FGD2 anti-numuCC 1pi", "FGD2 anti-numuCC other", "FGD2 NuMuBkg CC0pi in AntiNu Mode", "FGD2 NuMuBkg CC1pi in AntiNu Mode", "FGD2 NuMuBkg CCother in AntiNu Mode"};  

  std::string samp_us[18] = {"FGD1_numuCC_0pi", "FGD1_numuCC_1pi", "FGD1_numuCC_other", "FGD1_anti-numuCC_0pi", "FGD1_anti-numuCC_1pi", "FGD1_anti-numuCC_other", "FGD1_NuMuBkg_CC0pi_in_AntiNu_Mode", "FGD1_NuMuBkg_CC1pi_in_AntiNu_Mode", "FGD1_NuMuBkg_CCother_in_AntiNu_Mode", "FGD2_numuCC_0pi", "FGD2_numuCC_1pi", "FGD2_numuCC_other", "FGD2_anti-numuCC_0pi","FGD2_anti-numuCC_1pi", "FGD2_anti-numuCC_other", "FGD2_NuMuBkg_CC0pi_in_AntiNu_Mode", "FGD2_NuMuBkg_CC1pi_in_AntiNu_Mode", "FGD2_NuMuBkg_CCother_in_AntiNu_Mode"};

  for(int s=0; s<18; s++){

    std::string name = samp[s] + "/" + samp_us[s] + "_nom_MeanlnL";
    //    name = "FGD1 numuCC 0pi/FGD1_numuCC_0pi_nom_mean_ratio";
    TH2Poly* hist = (TH2Poly*)file_->Get(name.c_str())->Clone();

    hist->SetTitle("");
    hist->GetYaxis()->SetTitleOffset(0.85);
    hist->GetXaxis()->SetTitleOffset(1.1);
    hist->GetZaxis()->SetTitleOffset(1.0);
    hist->GetXaxis()->SetTitleSize(0.05);
    hist->GetXaxis()->SetRangeUser(0,5000);
    hist->GetZaxis()->SetRangeUser(-1,15);
    hist->GetYaxis()->SetRangeUser(0.7,1.0);
    hist->GetYaxis()->SetTitleSize(0.05);
    hist->GetZaxis()->SetTitleSize(0.05);
    hist->GetXaxis()->SetTitle("p_{#mu} (MeV)");
    hist->GetZaxis()->SetTitle("-2LLH_{Sample} x sign(MC-Data)");
    hist->Draw("colz");
    c->Print((std::string("llhcont_")+samp_us[s]+std::string(".pdf")).c_str());
  }
}
