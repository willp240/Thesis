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

void th2polynom(std::string file0) {

  gStyle->SetPalette(51);

  TCanvas *c = new TCanvas("canv", "canv", 1580, 1080);
  TFile *file_ = TFile::Open(file0.c_str());
 
  c->SetTopMargin(0.05);
  c->SetBottomMargin(0.14);
  c->SetLeftMargin(0.1);
  c->SetRightMargin(0.16);
  
  std::string samp[18] = {"FGD1_numuCC_0pi", "FGD1_numuCC_1pi", "FGD1_numuCC_other", "FGD1_anti-numuCC_0pi", "FGD1_anti-numuCC_1pi", "FGD1_anti-numuCC_other", "FGD1_NuMuBkg_CC0pi_in_AntiNu_Mode", "FGD1_NuMuBkg_CC1pi_in_AntiNu_Mode", "FGD1_NuMuBkg_CCother_in_AntiNu_Mode", "FGD2_numuCC_0pi", "FGD2_numuCC_1pi", "FGD2_numuCC_other", "FGD2_anti-numuCC_0pi","FGD2_anti-numuCC_1pi", "FGD2_anti-numuCC_other", "FGD2_NuMuBkg_CC0pi_in_AntiNu_Mode", "FGD2_NuMuBkg_CC1pi_in_AntiNu_Mode", "FGD2_NuMuBkg_CCother_in_AntiNu_Mode"};  

  for(int s=0; s<18; s++){

    std::cout << "getting " << samp[s] << std::endl;
    std::string name = "MC_" + samp[s];

    TH2Poly *hist = (TH2Poly*) file_->Get(name.c_str())->Clone();
    

    hist->SetTitle("");
    hist->GetYaxis()->SetTitleOffset(0.85);
    hist->GetXaxis()->SetTitleOffset(1.1);
    hist->GetZaxis()->SetTitleOffset(1.0);
    hist->GetXaxis()->SetTitleSize(0.05);
    hist->GetXaxis()->SetRangeUser(0,5000);
    hist->GetYaxis()->SetTitleSize(0.05);
    hist->GetZaxis()->SetTitleSize(0.05);
    hist->GetXaxis()->SetTitle("p_{#mu} (MeV)");
    hist->Draw("colz");
    c->Print((std::string("TH2PolyNom_")+name+std::string(".pdf")).c_str());
    std::cout << "got " << name << std::endl;
  }
}
