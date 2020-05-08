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


void RatioPolys(){

  //open up poly template file
  TFile *file1 = TFile::Open("5MayFlucAsimovFDFit_TH2DBinned_0_ND280_nom.root");
  TFile *file2 = TFile::Open("th2dFlucAsimov_binnedTH2D2.root");
  TFile *outfile = new TFile("ratioth2d.root","RECREATE");

  int nsamples=18;
  //define sample names
  std::string sampleName[18] = {"FGD1_numuCC_0pi", "FGD1_numuCC_1pi", "FGD1_numuCC_other", "FGD1_anti-numuCC_0pi", "FGD1_anti-numuCC_1pi", "FGD1_anti-numuCC_other", "FGD1_NuMuBkg_CC0pi_in_AntiNu_Mode", "FGD1_NuMuBkg_CC1pi_in_AntiNu_Mode", "FGD1_NuMuBkg_CCother_in_AntiNu_Mode", "FGD2_numuCC_0pi", "FGD2_numuCC_1pi", "FGD2_numuCC_other", "FGD2_anti-numuCC_0pi", "FGD2_anti-numuCC_1pi", "FGD2_anti-numuCC_other", "FGD2_NuMuBkg_CC0pi_in_AntiNu_Mode", "FGD2_NuMuBkg_CC1pi_in_AntiNu_Mode", "FGD2_NuMuBkg_CCother_in_AntiNu_Mode"};

  TH2Poly* polysample1[18];
  TH2Poly* polysample2[18];
  TH2Poly* ratiosample[18];

  for(int s=0; s<18; s++){

    polysample1[s] = (TH2Poly*)file1->Get(("MC_"+sampleName[s]).c_str());
    polysample2[s] = (TH2Poly*)file2->Get(("MC_"+sampleName[s]).c_str());
    ratiosample[s] = (TH2Poly*)file2->Get(("MC_"+sampleName[s]).c_str());

    //loop over bins
    for(int i=1; i<=polysample1[s]->GetNumberOfBins(); i++){

      double p1 = polysample1[s]->GetBinContent(i);
      double p2 = polysample2[s]->GetBinContent(i);

      double ratio;
      
      if(p2>0)
	ratio=p1/p2;
      else
	ratio=p1;

      if(ratio!=1.0)
	std::cout << "Sample " << s << ", Bin " << i << ", " << p1 << " " << p2 << " " << ratio << std::endl;

      ratiosample[s]->SetBinContent(i, ratio);	
    }
  }
  for(int s=0; s<18; s++){
    outfile->cd();
    ratiosample[s]->Write(("Ratio_"+sampleName[s]).c_str());
  }
}
