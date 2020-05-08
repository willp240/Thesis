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


void divide(std::string file0) {

  TCanvas *c = new TCanvas("canv", "canv", 1080, 1080);
  TFile *infile = TFile::Open(file0.c_str());
 
  c->SetTopMargin(0.07);
  c->SetBottomMargin(0.13);
  c->SetLeftMargin(0.13);
  c->SetRightMargin(0.18);

  int nsamples=18;

  std::string sampleName[18] = {"FGD1_numuCC_0pi", "FGD1_numuCC_1pi", "FGD1_numuCC_other", "FGD2_numuCC_0pi", "FGD2_numuCC_1pi", "FGD2_numuCC_other", "FGD1_anti-numuCC_0pi", "FGD1_anti-numuCC_1pi", "FGD1_anti-numuCC_other", "FGD2_anti-numuCC_0pi", "FGD2_anti-numuCC_1pi", "FGD2_anti-numuCC_other", "FGD1_NuMuBkg_CC0pi_in_AntiNu_Mode", "FGD1_NuMuBkg_CC1pi_in_AntiNu_Mode", "FGD1_NuMuBkg_CCother_in_AntiNu_Mode", "FGD2_NuMuBkg_CC0pi_in_AntiNu_Mode", "FGD2_NuMuBkg_CC1pi_in_AntiNu_Mode", "FGD2_NuMuBkg_CCother_in_AntiNu_Mode"};
  
  for (int s=0; s<nsamples; s++){

    std::string mcname = "MC_"+sampleName[s];
    std::string dataname = "DATA_"+sampleName[s];
    std::string ratname = "RATIO_"+sampleName[s];
    
    TH2Poly *hMC = (TH2Poly*)(infile->Get((mcname).c_str()));
    if(!hMC)
      std::cout << "No hMC!" << std::endl;
    
    TH2Poly *hDAT = (TH2Poly*)(infile->Get((dataname).c_str()));
    if(!hDAT)
      std::cout << "No hDAT!" << std::endl;
    
    TH2Poly *hRatio = (TH2Poly*)(infile->Get((mcname).c_str()));
    hRatio->SetTitle(ratname.c_str());
    hRatio->SetName(ratname.c_str());
    
    for(int i=1; i<hRatio->GetNumberOfBins(); i++){
      std::cout << "Bin " << i << " MC: " << hMC->GetBinContent(i) << " Dat: " << hDAT->GetBinContent(i) << " Ratio: " << hMC->GetBinContent(i)/hDAT->GetBinContent(i) << std::endl;
      if(hDAT->GetBinContent(i)>0)
	hRatio->SetBinContent(i,hMC->GetBinContent(i)/hDAT->GetBinContent(i));
      // else
      //std::cout << "Bin " << i << " Dat: " << hDAT->GetBinContent(i) << ", MC: " << hMC->GetBinContent(i) << std::endl;
    }
    hRatio->GetZaxis()->SetRangeUser(0,2);
    hRatio->GetXaxis()->SetRangeUser(0,5000);
    hRatio->Draw("colz");
    c->Print((sampleName[s]+"_ratioflucpoly.pdf").c_str());
  }
}
