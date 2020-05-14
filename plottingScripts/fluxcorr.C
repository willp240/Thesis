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
#include "TColor.h"

void fluxcorr(std::string file0) {

  TCanvas *c = new TCanvas("canv", "canv", 1580, 1380);
  TFile *file_ = TFile::Open(file0.c_str());
  TFile *outfile = new TFile("test.root", "RECREATE");
  
  c->SetTopMargin(0.05);
  c->SetBottomMargin(0.18);
  c->SetLeftMargin(0.18);
  c->SetRightMargin(0.22);

  TMatrixT<double>* cov;
  cov = (TMatrixT<double>*)file_->Get("total_flux_cov");

  TH2D* cov2 = new TH2D("","",100,0,100,100,0,100);
  TH2D* corr = new TH2D("","",100,0,100,100,0,100);
  
  for (int i=1;i<101;i++)
    {
      //      cov2->GetXaxis()->SetBinLabel(i,"");
      //cov2->GetYaxis()->SetBinLabel(i,"");
      for( int j=1; j<101; j++){
	cov2->SetBinContent(i,j,cov(i-1,j-1));
      }
    }

  for(int i=1; i<= cov2->GetXaxis()->GetNbins(); i++){
    for(int j=1; j<= cov2->GetYaxis()->GetNbins();j++){
      corr->SetBinContent(i,j, cov2->GetBinContent(i,j)/sqrt(cov2->GetBinContent(i,i)*cov2->GetBinContent(j,j)));
    }
  }
  
  corr->GetXaxis()->SetBinLabel(1,"ND280 FHC #nu_{#mu}");//11
  corr->GetXaxis()->SetBinLabel(12,"ND280 FHC #bar{#nu_{#mu}}");//5
  corr->GetXaxis()->SetBinLabel(17,"ND280 FHC #nu_{e}");//7
  corr->GetXaxis()->SetBinLabel(24,"ND280 FHC #bar{#nu_{e}}");//2
  corr->GetXaxis()->SetBinLabel(26,"ND280 RHC #nu_{#mu}");//5
  corr->GetXaxis()->SetBinLabel(31,"ND280 RHC #bar{#nu_{#mu}}");//11
  corr->GetXaxis()->SetBinLabel(42,"ND280 RHC #nu_{e}");//2
  corr->GetXaxis()->SetBinLabel(44,"ND280 RHC #bar{#nu_{e}}");//7
  corr->GetXaxis()->SetBinLabel(51,"SK FHC #nu_{#mu}");//11
  corr->GetXaxis()->SetBinLabel(62,"SK FHC #bar{#nu_{#mu}}");//5
  corr->GetXaxis()->SetBinLabel(67,"SK FHC #nu_{e}");//7
  corr->GetXaxis()->SetBinLabel(74,"SK FHC #bar{#nu_{e}}");//2
  corr->GetXaxis()->SetBinLabel(76,"SK RHC #nu_{#mu}");//5
  corr->GetXaxis()->SetBinLabel(81,"SK RHC #bar{#nu_{#mu}}");//11
  corr->GetXaxis()->SetBinLabel(92,"SK RHC #nu_{e}");//2
  corr->GetXaxis()->SetBinLabel(94,"SK RHC #bar{#nu_{e}}");//7

  corr->GetYaxis()->SetBinLabel(1,"ND280 FHC #nu_{#mu}");//11
  corr->GetYaxis()->SetBinLabel(12,"ND280 FHC #bar{#nu_{#mu}}");//5
  corr->GetYaxis()->SetBinLabel(17,"ND280 FHC #nu_{e}");//7
  corr->GetYaxis()->SetBinLabel(24,"ND280 FHC #bar{#nu_{e}}");//2
  corr->GetYaxis()->SetBinLabel(26,"ND280 RHC #nu_{#mu}");//5
  corr->GetYaxis()->SetBinLabel(31,"ND280 RHC #bar{#nu_{#mu}}");//11
  corr->GetYaxis()->SetBinLabel(42,"ND280 RHC #nu_{e}");//2
  corr->GetYaxis()->SetBinLabel(44,"ND280 RHC #bar{#nu_{e}}");//7
  corr->GetYaxis()->SetBinLabel(51,"SK FHC #nu_{#mu}");//11
  corr->GetYaxis()->SetBinLabel(62,"SK FHC #bar{#nu_{#mu}}");//5
  corr->GetYaxis()->SetBinLabel(67,"SK FHC #nu_{e}");//7
  corr->GetYaxis()->SetBinLabel(74,"SK FHC #bar{#nu_{e}}");//2
  corr->GetYaxis()->SetBinLabel(76,"SK RHC #nu_{#mu}");//5
  corr->GetYaxis()->SetBinLabel(81,"SK RHC #bar{#nu_{#mu}}");//11
  corr->GetYaxis()->SetBinLabel(92,"SK RHC #nu_{e}");//2
  corr->GetYaxis()->SetBinLabel(94,"SK RHC #bar{#nu_{e}}");//7

  corr->GetXaxis()->SetLabelSize(0.03);
  corr->GetYaxis()->SetLabelSize(0.03);
  corr->GetXaxis()->SetNdivisions(1010);
  corr->GetZaxis()->SetRangeUser(-1,1);
  corr->GetZaxis()->SetTitle("Correlation");
  corr->GetZaxis()->SetTitleOffset(1.5);
  corr->Draw("colz");
  c->Print((std::string("fluxcorr")+std::string(".pdf")).c_str());
  
  outfile->cd();
  corr->Write();
}
