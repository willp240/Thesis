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


void batch(std::string file0) {

  TCanvas *c = new TCanvas("canv", "canv", 1580, 1080);
  TFile *file_ = TFile::Open(file0.c_str());
 
  c->SetTopMargin(0.05);
  c->SetBottomMargin(0.14);
  c->SetLeftMargin(0.09);
  c->SetRightMargin(0.09);
  
  c->Print((std::string("batch")+std::string(".pdf[")).c_str());

  file_->cd("Batched_means");

  b_4_batch->GetYaxis()->SetTitleOffset(1.3);
  b_4_batch->GetXaxis()->SetTitleOffset(1.7);
  b_4_batch->GetYaxis()->SetTitle("Mean Parameter Value");
  b_4_batch->GetXaxis()->SetTitle("Steps");
  b_4_batch->SetTitle("");
  b_4_batch->Draw();
  c->Print((std::string("batch")+std::string(".pdf")).c_str());
  c->Print((std::string("batch")+std::string(".pdf]")).c_str());
}
