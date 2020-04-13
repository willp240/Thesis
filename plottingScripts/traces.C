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


void traces(std::string file0) {

  TCanvas *c = new TCanvas("canv", "canv", 1580, 1080);
  TFile *file_ = TFile::Open(file0.c_str());
 
  c->SetTopMargin(0.07);
  c->SetBottomMargin(0.11);
  c->SetLeftMargin(0.10);
  c->SetRightMargin(0.08);
  
  c->Print((std::string("burnin")+std::string(".pdf[")).c_str());

  file_->cd("Trace");

  b_0->GetYaxis()->SetTitleOffset(1.4);
  b_0->GetXaxis()->SetRangeUser(0,50000);
  b_0->SetTitle("");
  b_0->Draw();
  c->Print("burnin.png");
  c->Print((std::string("burnin")+std::string(".pdf")).c_str());
  c->Print((std::string("burnin")+std::string(".pdf]")).c_str());
}
