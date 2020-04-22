#include "TH1D.h"
#include "TH2D.h"
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

void plotPOTInfo(std::string fileMC)  {

  //TCanvas *c = new TCanvas();
  TFile *mcfile = TFile::Open(fileMC.c_str());
  TBranch** ptr;
  void* add;
  double POT[13];


  //  c->SetTopMargin(0.07);
  //c->SetBottomMargin(011);
  //c->SetLeftMargin(0.10);
  //c->SetRightMargin(0.18);

  TTree *TreeMC = (TTree*)mcfile->Get("header");
  // TreeMC->SetBranchAddress("POTInfo_v2", add, ptr);
  
  for(int i = 0; i < 60; i++) {
    TreeMC->GetEntry(0);
    TreeMC->SetBranchAddress("POTInfo_v2_POTInfo",      &POT);
    std::cout << POT[5] << std::endl;
  }

  //  TreeMC->Draw("_POTInfo>>hPOTMC","","");
  
  //TH1D *hMC = (TH1D*)gDirectory->Get("hPOTMC");


  
    //  hMC->SetLineColor(kRed);
    //hMC->GetXaxis()->SetRangeUser(100E15,2000E15);
    //hMC->GetYaxis()->SetRangeUser(0,10000);
    //hMC->Draw("");

}
