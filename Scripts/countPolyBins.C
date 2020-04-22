#include <TH2Poly.h>
#include <TH2D.h>
#include <vector>
#include <TMath.h>
#include <TFile.h>
#include <TTree.h>
#include <THStack.h>
#include <iostream> 

void countPolyBins(std::string PolyFileName){
   
  TFile* filePoly = TFile::Open(PolyFileName.c_str());
  
  const double nsamples = 18;
  std::string sampleName[nsamples] = {"FGD1_numuCC_0pi", "FGD1_numuCC_1pi", "FGD1_numuCC_other", "FGD1_anti-numuCC_0pi", "FGD1_anti-numuCC_1pi", "FGD1_anti-numuCC_other", "FGD1_NuMuBkg_CC0pi_in_AntiNu_Mode", "FGD1_NuMuBkg_CC1pi_in_AntiNu_Mode", "FGD1_NuMuBkg_CCother_in_AntiNu_Mode", "FGD2_numuCC_0pi", "FGD2_numuCC_1pi", "FGD2_numuCC_other", "FGD2_anti-numuCC_0pi", "FGD2_anti-numuCC_1pi", "FGD2_anti-numuCC_other", "FGD2_NuMuBkg_CC0pi_in_AntiNu_Mode", "FGD2_NuMuBkg_CC1pi_in_AntiNu_Mode", "FGD2_NuMuBkg_CCother_in_AntiNu_Mode"};
  int count=0;
  //loop over all samples
  for(int s = 0; s<nsamples; s++){//nsamples; s++){

    std::cout << "Doing for " << sampleName[s] << std::endl;
    //get the poly
    std::string getname = "" + sampleName[s];

    TH2Poly *poly=(TH2Poly*)(filePoly->Get(getname.c_str()));
    if(!poly)
      std::cout << "No hist " << getname << std::endl;
    poly->SetName(("p"+getname).c_str());
    poly->SetTitle(("p"+getname).c_str());
    count+=poly->GetNumberOfBins();
    std::cout << sampleName[s] << " " << poly->GetNumberOfBins() << std::endl;
  }
  std::cout << count << std::endl;
}//end convertPolyToTH2D

