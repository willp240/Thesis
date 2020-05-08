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


void fillbins(std::string file0){

  //open up poly template file
  TFile *polyfile = TFile::Open("../inputs/polybinning_minwidth100MeV_template.root");
  TFile *th2dfile = TFile::Open("../inputs/polybinning_BANFFTH2Ds_template.root");
  //open th2d file with events
  TFile *flucfile = TFile::Open(file0.c_str());
  TFile *outfileNew = new TFile("binnedTH2D2_No0Bins.root", "RECREATE");
  TFile *outfileOrig = new TFile("binnedPoly2_No0Bins.root","RECREATE");

  double p,t;
  int sam;
  TTree *t1 = flucfile->Get("EventTree");
  t1->SetBranchAddress("momentum",&p);
  t1->SetBranchAddress("costheta",&t);
  t1->SetBranchAddress("sample",&sam);
  
  bool debug =false;

  int nsamples=18;
  //define sample names
  std::string sampleName[18] = {"FGD1_numuCC_0pi", "FGD1_numuCC_1pi", "FGD1_numuCC_other", "FGD1_anti-numuCC_0pi", "FGD1_anti-numuCC_1pi", "FGD1_anti-numuCC_other", "FGD1_NuMuBkg_CC0pi_in_AntiNu_Mode", "FGD1_NuMuBkg_CC1pi_in_AntiNu_Mode", "FGD1_NuMuBkg_CCother_in_AntiNu_Mode", "FGD2_numuCC_0pi", "FGD2_numuCC_1pi", "FGD2_numuCC_other", "FGD2_anti-numuCC_0pi", "FGD2_anti-numuCC_1pi", "FGD2_anti-numuCC_other", "FGD2_NuMuBkg_CC0pi_in_AntiNu_Mode", "FGD2_NuMuBkg_CC1pi_in_AntiNu_Mode", "FGD2_NuMuBkg_CCother_in_AntiNu_Mode"};

  TH2Poly* polysample[18];
  TH2Poly* th2dsample[18];

  for(int s=0; s<18; s++){
    polysample[s] = (TH2Poly*)polyfile->Get(sampleName[s].c_str())->Clone();
    th2dsample[s] = (TH2Poly*)th2dfile->Get(sampleName[s].c_str())->Clone();
    polysample[s]->Reset("");
    th2dsample[s]->Reset("");
  }

  //loop over entries
  for(int i=0; i<t1->GetEntries(); i++){
    t1->GetEntry(i);

    int polybin = polysample[sam]->FindBin(p,t);
    int th2dbin = th2dsample[sam]->FindBin(p,t);

    if(polybin<=0 || polybin>polysample[sam]->GetNumberOfBins()){
      std::cout << "ERROR: PolyBin " << polybin << ". Eve " << i << ", Sample:" << sam << p << ", " << t << std::endl;
      return 0;
    }
    if(th2dbin<=0 || th2dbin>th2dsample[sam]->GetNumberOfBins()){
      std::cout << "ERROR: Th2dBin " << th2dbin << ". Eve " << i << ", Sample:" << sam << p << ", " << t << std::endl;
      return 0;
    }

    // fill the polys
    polysample[sam]->SetBinContent(polybin,polysample[sam]->GetBinContent(polybin)+1);
    th2dsample[sam]->SetBinContent(th2dbin,th2dsample[sam]->GetBinContent(th2dbin)+1);
  }
  
  for(int s=0; s<18; s++){

    for(int i=1; i<th2dsample[s]->GetNumberOfBins(); i++){
      if(th2dsample[s]->GetBinContent(i)<1){
	std::cout << "Sample " << s << " Bin " << i << " Content " << th2dsample[s]->GetBinContent(i) << std::endl;
	th2dsample[s]->SetBinContent(i,th2dsample[s]->GetBinContent(i)+1);
      }
    }

    outfileNew->cd();
    th2dsample[s]->Write(("MC_"+sampleName[s]).c_str());
    outfileOrig->cd();
    polysample[s]->Write(("MC_"+sampleName[s]).c_str());
  }
}
