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

void pval2(std::string file0) {

  gStyle->SetPalette(51);

  TCanvas *c = new TCanvas("canv", "canv", 1580, 1080);
  TFile *file_ = TFile::Open(file0.c_str());
 
  c->SetTopMargin(0.05);
  c->SetBottomMargin(0.14);
  c->SetLeftMargin(0.1);
  c->SetRightMargin(0.16);
  
  std::string name = "lnLDrawHist";
  //    name = "FGD1 numuCC 0pi/FGD1_numuCC_0pi_nom_mean_ratio";
  TH2D* hist = (TH2D*)file_->Get(name.c_str())->Clone();

  hist->GetXaxis()->SetRangeUser(4200,5400);

  // Now add the TLine going diagonally
  double minimum = hist->GetXaxis()->GetBinLowEdge(1);
  minimum=4200;
  if (hist->GetYaxis()->GetBinLowEdge(1) > minimum) {
    minimum = hist->GetYaxis()->GetBinLowEdge(1);
  }
  double maximum = hist->GetXaxis()->GetBinLowEdge(hist->GetXaxis()->GetNbins());
  maximum=5400;
  if (hist->GetYaxis()->GetBinLowEdge(hist->GetYaxis()->GetNbins()) < maximum) {
    maximum = hist->GetYaxis()->GetBinLowEdge(hist->GetYaxis()->GetNbins());
  }
  TLine *TempLine = new TLine(minimum, minimum, maximum, maximum);
  TempLine->SetLineColor(kRed);
  TempLine->SetLineWidth(2);

  hist->SetTitle("");
  hist->GetYaxis()->SetTitleOffset(1);
  hist->GetXaxis()->SetTitleOffset(1.2);
  //  hist->GetZaxis()->SetTitleOffset(1.0);
  hist->GetXaxis()->SetTitleSize(0.05);
  hist->GetYaxis()->SetTitleSize(0.05);
  hist->GetZaxis()->SetTitleSize(0.05);
  hist->GetZaxis()->SetTitle("Steps");
  hist->Draw("colz");
  TempLine->Draw("same");
  c->Print((std::string("pval2_")+std::string(".pdf")).c_str());

}
