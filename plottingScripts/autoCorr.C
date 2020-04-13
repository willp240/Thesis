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


void autoCorr(std::string file0) {

  TCanvas *c = new TCanvas("canv", "canv", 1080, 1080);
  TFile *file_ = TFile::Open(file0.c_str());
 
  c->SetTopMargin(0.07);
  c->SetBottomMargin(0.11);
  c->SetLeftMargin(0.10);
  c->SetRightMargin(0.18);
  
  c->Print((std::string("autoCorr")+std::string(".pdf[")).c_str());

  file_->cd("Auto_corr");
  
  xsec_0->SetLineColor(kRed+2);
  xsec_23->SetLineColor(kBlue+2);
  xsec_37->SetLineColor(kYellow+2);
  xsec_39->SetLineColor(kGreen+2);
  xsec_0->SetTitle("");
  xsec_0->GetYaxis()->SetTitleOffset(1.2);

  xsec_0->Draw();
  xsec_23->Draw("same");
  xsec_37->Draw("same");
  xsec_39->Draw("same");

  auto legend = new TLegend(0.49,0.73,0.79,0.89);
  legend->AddEntry(xsec_23,"Step Size Scale 0.04","l");
  legend->AddEntry(xsec_0,"Step Size Scale 0.05","l");
  legend->AddEntry(xsec_37,"Step Size Scale 0.06","l");
  legend->AddEntry(xsec_39,"Step Size Scale 0.07","l");
  legend->Draw();
  c->Print((std::string("autoCorr")+std::string(".pdf")).c_str());
  c->Print((std::string("autoCorr")+std::string(".pdf]")).c_str());
}
