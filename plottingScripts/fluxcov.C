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

void fluxcov(std::string file0) {

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
  
  for (int i=1;i<101;i++)
    {
      //      cov2->GetXaxis()->SetBinLabel(i,"");
      //cov2->GetYaxis()->SetBinLabel(i,"");
      for( int j=1; j<101; j++){
	cov2->SetBinContent(i,j,cov(i-1,j-1));
      }
    }
  
  cov2->GetXaxis()->SetBinLabel(1,"ND280 FHC #nu_{#mu}");//11
  cov2->GetXaxis()->SetBinLabel(12,"ND280 FHC #bar{#nu_{#mu}}");//5
  cov2->GetXaxis()->SetBinLabel(17,"ND280 FHC #nu_{e}");//7
  cov2->GetXaxis()->SetBinLabel(24,"ND280 FHC #bar{#nu_{e}}");//2
  cov2->GetXaxis()->SetBinLabel(26,"ND280 RHC #nu_{#mu}");//5
  cov2->GetXaxis()->SetBinLabel(31,"ND280 RHC #bar{#nu_{#mu}}");//11
  cov2->GetXaxis()->SetBinLabel(42,"ND280 RHC #nu_{e}");//2
  cov2->GetXaxis()->SetBinLabel(44,"ND280 RHC #bar{#nu_{e}}");//7
  cov2->GetXaxis()->SetBinLabel(51,"SK FHC #nu_{#mu}");//11
  cov2->GetXaxis()->SetBinLabel(62,"SK FHC #bar{#nu_{#mu}}");//5
  cov2->GetXaxis()->SetBinLabel(67,"SK FHC #nu_{e}");//7
  cov2->GetXaxis()->SetBinLabel(74,"SK FHC #bar{#nu_{e}}");//2
  cov2->GetXaxis()->SetBinLabel(76,"SK RHC #nu_{#mu}");//5
  cov2->GetXaxis()->SetBinLabel(81,"SK RHC #bar{#nu_{#mu}}");//11
  cov2->GetXaxis()->SetBinLabel(92,"SK RHC #nu_{e}");//2
  cov2->GetXaxis()->SetBinLabel(94,"SK RHC #bar{#nu_{e}}");//7

  cov2->GetYaxis()->SetBinLabel(1,"ND280 FHC #nu_{#mu}");//11
  cov2->GetYaxis()->SetBinLabel(12,"ND280 FHC #bar{#nu_{#mu}}");//5
  cov2->GetYaxis()->SetBinLabel(17,"ND280 FHC #nu_{e}");//7
  cov2->GetYaxis()->SetBinLabel(24,"ND280 FHC #bar{#nu_{e}}");//2
  cov2->GetYaxis()->SetBinLabel(26,"ND280 RHC #nu_{#mu}");//5
  cov2->GetYaxis()->SetBinLabel(31,"ND280 RHC #bar{#nu_{#mu}}");//11
  cov2->GetYaxis()->SetBinLabel(42,"ND280 RHC #nu_{e}");//2
  cov2->GetYaxis()->SetBinLabel(44,"ND280 RHC #bar{#nu_{e}}");//7
  cov2->GetYaxis()->SetBinLabel(51,"SK FHC #nu_{#mu}");//11
  cov2->GetYaxis()->SetBinLabel(62,"SK FHC #bar{#nu_{#mu}}");//5
  cov2->GetYaxis()->SetBinLabel(67,"SK FHC #nu_{e}");//7
  cov2->GetYaxis()->SetBinLabel(74,"SK FHC #bar{#nu_{e}}");//2
  cov2->GetYaxis()->SetBinLabel(76,"SK RHC #nu_{#mu}");//5
  cov2->GetYaxis()->SetBinLabel(81,"SK RHC #bar{#nu_{#mu}}");//11
  cov2->GetYaxis()->SetBinLabel(92,"SK RHC #nu_{e}");//2
  cov2->GetYaxis()->SetBinLabel(94,"SK RHC #bar{#nu_{e}}");//7
  
  cov2->GetXaxis()->SetLabelSize(0.03);
  cov2->GetYaxis()->SetLabelSize(0.03);
  cov2->GetXaxis()->SetNdivisions(1010);
  cov2->GetZaxis()->SetRangeUser(-0.0125,0.0125);
  cov2->GetZaxis()->SetTitle("Covariance");
  cov2->GetZaxis()->SetTitleOffset(1.5);
  cov2->Draw("colz");
  c->Print((std::string("fluxcov")+std::string(".pdf")).c_str());

  cov2->GetXaxis()->SetRangeUser(0,11);
  cov2->GetYaxis()->SetRangeUser(0,11);

  cov2->GetXaxis()->SetBinLabel(1,"0-0.4 GeV");
  cov2->GetXaxis()->SetBinLabel(2,"0.4-0.5 GeV");
  cov2->GetXaxis()->SetBinLabel(3,"0.5-0.6 GeV");
  cov2->GetXaxis()->SetBinLabel(4,"0.6-0.7 GeV");
  cov2->GetXaxis()->SetBinLabel(5,"0.7-1.0 GeV");
  cov2->GetXaxis()->SetBinLabel(6,"1.0-1.5 GeV");
  cov2->GetXaxis()->SetBinLabel(7,"1.5-2.5 GeV");
  cov2->GetXaxis()->SetBinLabel(8,"2.5-3.5 GeV");
  cov2->GetXaxis()->SetBinLabel(9,"3.5-5.0 GeV");
  cov2->GetXaxis()->SetBinLabel(10,"5.0-7.0 GeV");
  cov2->GetXaxis()->SetBinLabel(11,"7.0-30.0 GeV");
  cov2->GetYaxis()->SetBinLabel(1,"0-0.4 GeV");
  cov2->GetYaxis()->SetBinLabel(2,"0.4-0.5 GeV");
  cov2->GetYaxis()->SetBinLabel(3,"0.5-0.6 GeV");
  cov2->GetYaxis()->SetBinLabel(4,"0.6-0.7 GeV");
  cov2->GetYaxis()->SetBinLabel(5,"0.7-1.0 GeV");
  cov2->GetYaxis()->SetBinLabel(6,"1.0-1.5 GeV");
  cov2->GetYaxis()->SetBinLabel(7,"1.5-2.5 GeV");
  cov2->GetYaxis()->SetBinLabel(8,"2.5-3.5 GeV");
  cov2->GetYaxis()->SetBinLabel(9,"3.5-5.0 GeV");
  cov2->GetYaxis()->SetBinLabel(10,"5.0-7.0 GeV");
  cov2->GetYaxis()->SetBinLabel(11,"7.0-30.0 GeV");

  cov2->Draw("colz");
  c->Print((std::string("fluxcov1sample")+std::string(".pdf")).c_str());  

  outfile->cd();
  cov2->Write();
}
