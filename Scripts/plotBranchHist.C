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

void plotBranchHist(std::string fileMC, std::string fileData)  {

  TCanvas *c = new TCanvas();
  c->Divide(1,2);
  TFile *mcfile = TFile::Open(fileMC.c_str());
  TFile *datafile = TFile::Open(fileData.c_str());

  c->SetTopMargin(0.07);
  c->SetBottomMargin(011);
  c->SetLeftMargin(0.10);
  c->SetRightMargin(0.18);

  TTree *TreeMC = (TTree*)mcfile->Get("header");
  TTree *TreeData = (TTree*)datafile->Get("header");

  TreeMC->Draw("_POTInfo>>hPOTMC","","");
  TreeData->Draw("_POTInfo>>hPOTData","","");
  
  TH1D *hMC = (TH1D*)gDirectory->Get("hPOTMC");
  TH1D *hData = (TH1D*)gDirectory->Get("hPOTData");

  std::cout << "MC Integral " << hMC->Integral() << std::endl;
  std::cout << "Data Integral " << hData->Integral() << std::endl;
  
  hMC->SetLineColor(kRed);
  hData->SetLineColor(kBlue);
  hMC->GetXaxis()->SetRangeUser(100E15,2000E15);
  hMC->GetYaxis()->SetRangeUser(0,10000);
  hData->GetXaxis()->SetRangeUser(100E15,2000E15);
  hData->GetYaxis()->SetRangeUser(0,10000);
  hMC->Draw("");
  hData->Draw("same");

}
