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


void par1d(std::string file0) {

  gStyle->SetPalette(51);

  TCanvas *c = new TCanvas("canv", "canv", 1580, 1080);
  TFile *file_ = TFile::Open(file0.c_str());

  gDirectory->cd("params");

  TCanvas *c1 = (TCanvas*) file_->Get("params/xsec_18");

  c1->Draw();

  c1->Print("asmvEbCnu.pdf");
 
  c->SetTopMargin(0.07);
  c->SetBottomMargin(0.11);
  c->SetLeftMargin(0.11);
  c->SetRightMargin(0.14);
  
  gDirectory->cd();
  

}
