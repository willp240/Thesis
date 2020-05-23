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

void ndpred(std::string prior, std::string pred) {

  TFile *fileprior = TFile::Open(prior.c_str()); 
  TFile *filepred = TFile::Open(pred.c_str());


  //  TString spec[5] = {"numu", "numubar", "nue", "nuebar", "nue1pi"};
  TString samp[18] = {"FGD1_numuCC_0pi", "FGD1_numuCC_1pi", "FGD1_numuCC_other", "FGD1_anti-numuCC_0pi", "FGD1_anti-numuCC_1pi", "FGD1_anti-numuCC_other", "FGD1_NuMuBkg_CC0pi_in_AntiNu_Mode", "FGD1_NuMuBkg_CC1pi_in_AntiNu_Mode", "FGD1_NuMuBkg_CCother_in_AntiNu_Mode", "FGD2_numuCC_0pi", "FGD2_numuCC_1pi", "FGD2_numuCC_other", "FGD2_anti-numuCC_0pi","FGD2_anti-numuCC_1pi", "FGD2_anti-numuCC_other", "FGD2_NuMuBkg_CC0pi_in_AntiNu_Mode", "FGD2_NuMuBkg_CC1pi_in_AntiNu_Mode", "FGD2_NuMuBkg_CCother_in_AntiNu_Mode"};

  TH1D *total = new TH1D("","", 1000,0,1000);
  double totvec[2000];

  for(int i=0; i<2000; i++)
    totvec[i]=0.0;

  for(int s=0; s<1; s++){

    TH1D *ratesprior = new TH1D("","", 100,0,100); 
    TH1D *ratespred = new TH1D("","", 100,0,100);

    for(int i=0; i<2000; i++){
      
      TString name = Form("%s_%d",samp[s].Data(),i);
      //      std::cout << name << std::endl;
      TH2Poly* histprior = (TH2Poly*)fileprior->Get(name)->Clone();
      histprior->SetName("Prior");
      histprior->SetTitle("Prior");
      TH2Poly* histpred = (TH2Poly*)filepred->Get(name)->Clone();
      histpred->SetName("Post");
      histpred->SetTitle("Post");

      ratesprior->Fill(histprior->GetBinContent(10));
      ratespred->Fill(histpred->GetBinContent(10));
      //      totvec[i]+=hist->Integral();
    }
    ratesprior->Fit("gaus");
    ratespred->Fit("gaus");
    std::cout << "prior " << samp[s] << " Mean: " << ratesprior->GetMean() << ", RMS: " << ratesprior->GetRMS() << std::endl;
    std::cout << "pred " << samp[s] << " Mean: " << ratespred->GetMean() << ", RMS: " << ratespred->GetRMS() << std::endl;
    ratesprior->Draw();
    ratespred->SetLineColor(kRed);
    //    ratespred->Draw("");
    // ratesprior->Draw("same");
    ratesprior->SetAxisRange(ratesprior->GetMean()-16*ratesprior->GetRMS(), ratesprior->GetMean()+16*ratesprior->GetRMS());
    ratespred->SetAxisRange(ratespred->GetMean()-16*ratespred->GetRMS(), ratespred->GetMean()+16*ratespred->GetRMS());
    ratesprior->Fit("gaus");
    ratespred->Fit("gaus");
    std::cout << "prior " << samp[s] << " Mean: " << ratesprior->GetMean() << ", RMS: " << ratesprior->GetRMS() << std::endl;
    std::cout << "pred " << samp[s] << " Mean: " << ratespred->GetMean() << ", RMS: " << ratespred->GetRMS() << std::endl;
  }

  //  for(int i=0; i<2000; i++){
  //total->Fill(totvec[i]);
  // }

  // total->Fit("gaus");
  // std::cout < <"Total Mean: " << total->GetMean() << ", RMS: " << total->GetRMS() << std::endl;
  
}
