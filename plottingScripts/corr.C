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


void corr(std::string file0) {

  gStyle->SetPalette(51);

  TCanvas *c = new TCanvas("canv", "canv", 1580, 1080);
  TFile *file_ = TFile::Open(file0.c_str());
 
  c->SetTopMargin(0.07);
  c->SetBottomMargin(0.11);
  c->SetLeftMargin(0.11);
  c->SetRightMargin(0.14);
  
  TTree *pos = (TTree*)file_->Get("posteriors");
  file_->cd();

  /*  pos->Draw("b_0:b_1>>h0(200,1.05,1.3,200,1.05,1.35)","step>200000");
  TH1D *h0 = (TH1D*)gDirectory->Get("h0");
  h0->Draw("colz");
  h0->GetYaxis()->SetTitle("Beam Pararmeter 1");
  h0->GetYaxis()->SetTitleOffset(1.0);
  h0->GetXaxis()->SetTitle("Beam Parameter 0");
  h0->GetZaxis()->SetRangeUser(-1,h0->GetMaximum());
  h0->SetTitle("");
  h0->GetXaxis()->SetNoExponent(1);
  c->Print("corr1.png");*/

  pos->Draw("xsec_0:xsec_20>>h1(200,0.75,1.05,200,0.94,1.1)","step>200000");
  TH1D *h1 = (TH1D*)gDirectory->Get("h1");
  h1->Draw("colz");
  h1->GetYaxis()->SetTitle("#it{M_{A}^{QE}}");
  h1->GetYaxis()->SetTitleOffset(1.0);
  h1->GetXaxis()->SetTitle("#it{M_{A}^{RES}}");
  h1->GetZaxis()->SetRangeUser(-1,h1->GetMaximum());
  h1->SetTitle("");
  h1->GetXaxis()->SetNoExponent(1);
  c->Print("corr2.png");
  /*
  pos->Draw("LogL_systematic_xsec_cov:step>>h2(1000,0,800000,300,0,30)","LogL_systematic_xsec_cov<1000");
  TH1D *h2 = (TH1D*)gDirectory->Get("h2");
  h2->GetXaxis()->SetRangeUser(0,600000);
  h2->Draw("colz");
  h2->GetYaxis()->SetTitle("-2LogL_{Xsec}");
  h2->GetYaxis()->SetTitleOffset(1.0);
  h2->GetXaxis()->SetTitle("Step");
  h2->GetZaxis()->SetRangeUser(-1,150);
  h2->SetTitle("");
  h2->GetXaxis()->SetNoExponent(1);
  c->Print("llh_xsec.png");

  pos->Draw("LogL_systematic_nddet_cov:step>>h3(1000,0,800000,300,0,400)","LogL_systematic_nddet_cov<1000");
  TH1D *h3 = (TH1D*)gDirectory->Get("h3");
  h3->GetXaxis()->SetRangeUser(0,600000);
  h3->Draw("colz");
  h3->GetYaxis()->SetTitle("-2LogL_{NDDet}");
  h3->GetYaxis()->SetTitleOffset(1.0);
  h3->GetXaxis()->SetTitle("Step");
  h3->GetZaxis()->SetRangeUser(-1,400);
  h3->SetTitle("");
  h3->GetXaxis()->SetNoExponent(1);
  c->Print("llh_nddet.png");*/
}
