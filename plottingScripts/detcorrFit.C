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

void detcorrFit(std::string file0) {

  TCanvas *c = new TCanvas("canv", "canv", 1580, 1380);
  TFile *file_ = TFile::Open(file0.c_str());
  
  c->SetTopMargin(0.05);
  c->SetBottomMargin(0.22);
  c->SetLeftMargin(0.21);
  c->SetRightMargin(0.2);

  TH2D* cov;
  cov = (TH2D*)file_->Get("Total_Covariance_Matrix");
  TH2D* corr = (TH2D*)file_->Get("Total_Covariance_Matrix");

  for(int i=1; i<= cov->GetXaxis()->GetNbins(); i++){
    for(int j=1; j<= cov->GetYaxis()->GetNbins();j++){
      corr->SetBinContent(i,j, cov->GetBinContent(i,j)/sqrt(cov->GetBinContent(i,i)*cov->GetBinContent(j,j)));
    }
  }

  corr->SetTitle("");
  corr->GetXaxis()->SetLabelSize(0.03);
  corr->GetYaxis()->SetLabelSize(0.03);
  corr->GetXaxis()->SetNdivisions(1010);
  corr->GetZaxis()->SetRangeUser(-1,1);
  corr->GetZaxis()->SetTitle("Correlation");
  corr->GetZaxis()->SetTitleOffset(1.5);
  corr->GetXaxis()->SetRangeUser(0,3070);
  corr->GetYaxis()->SetRangeUser(0,3070);

  for (int i=1;i<3071;i++)
    {
      corr->GetXaxis()->SetBinLabel(i,"");
      corr->GetYaxis()->SetBinLabel(i,"");
    }
  
  corr->GetXaxis()->SetBinLabel(1,"FGD1 FHC #nu_{#mu} CC0#pi");
  corr->GetXaxis()->SetBinLabel(395,"FGD1 FHC #nu_{#mu} CC1#pi");
  corr->GetXaxis()->SetBinLabel(703,"FGD1 FHC #nu_{#mu} CCOther");
  corr->GetXaxis()->SetBinLabel(1009,"FGD2 FHC #nu_{#mu} CC0#pi");
  corr->GetXaxis()->SetBinLabel(1391,"FGD2 FHC #nu_{#mu} CC1#pi");
  corr->GetXaxis()->SetBinLabel(1669,"FGD2 FHC #nu_{#mu} CCOther");
  corr->GetXaxis()->SetBinLabel(1990,"FGD1 RHC #bar{#nu_{#mu} CC0#pi}");
  corr->GetXaxis()->SetBinLabel(2165,"FGD1 RHC #bar{#nu_{#mu} CC1#pi}");
  corr->GetXaxis()->SetBinLabel(2207,"FGD1 RHC #bar{#nu_{#mu} CC1Other}");
  corr->GetXaxis()->SetBinLabel(2278,"FGD2 RHC #bar{#nu_{#mu} CC0#pi}");
  corr->GetXaxis()->SetBinLabel(2450,"FGD2 RHC #bar{#nu_{#mu} CC1#pi}");
  corr->GetXaxis()->SetBinLabel(2489,"FGD2 RHC #bar{#nu_{#mu} CC1Other}");
  corr->GetXaxis()->SetBinLabel(2558,"FGD1 RHC #nu_{#mu} CC0#pi");
  corr->GetXaxis()->SetBinLabel(2697,"FGD1 RHC #nu_{#mu} CC1#pi");
  corr->GetXaxis()->SetBinLabel(2759,"FGD1 RHC #nu_{#mu} CCOther");
  corr->GetXaxis()->SetBinLabel(2822,"FGD2 RHC #nu_{#mu} CC0#pi");
  corr->GetXaxis()->SetBinLabel(2958,"FGD2 RHC #nu_{#mu} CC1#pi");
  corr->GetXaxis()->SetBinLabel(3010,"FGD2 RHC #nu_{#mu} CCOther");
    
  corr->GetYaxis()->SetBinLabel(1,"FGD1 FHC #nu_{#mu} CC0#pi");
  corr->GetYaxis()->SetBinLabel(395,"FGD1 FHC #nu_{#mu} CC1#pi");
  corr->GetYaxis()->SetBinLabel(703,"FGD1 FHC #nu_{#mu} CCOther");
  corr->GetYaxis()->SetBinLabel(1009,"FGD2 FHC #nu_{#mu} CC0#pi");
  corr->GetYaxis()->SetBinLabel(1391,"FGD2 FHC #nu_{#mu} CC1#pi");
  corr->GetYaxis()->SetBinLabel(1669,"FGD2 FHC #nu_{#mu} CCOther");
  corr->GetYaxis()->SetBinLabel(1990,"FGD1 RHC #bar{#nu_{#mu} CC0#pi}");
  corr->GetYaxis()->SetBinLabel(2165,"FGD1 RHC #bar{#nu_{#mu} CC1#pi}");
  corr->GetYaxis()->SetBinLabel(2207,"FGD1 RHC #bar{#nu_{#mu} CC1Other}");
  corr->GetYaxis()->SetBinLabel(2278,"FGD2 RHC #bar{#nu_{#mu} CC0#pi}");
  corr->GetYaxis()->SetBinLabel(2450,"FGD2 RHC #bar{#nu_{#mu} CC1#pi}");
  corr->GetYaxis()->SetBinLabel(2489,"FGD2 RHC #bar{#nu_{#mu} CC1Other}");
  corr->GetYaxis()->SetBinLabel(2558,"FGD1 RHC #nu_{#mu} CC0#pi");
  corr->GetYaxis()->SetBinLabel(2697,"FGD1 RHC #nu_{#mu} CC1#pi");
  corr->GetYaxis()->SetBinLabel(2759,"FGD1 RHC #nu_{#mu} CCOther");
  corr->GetYaxis()->SetBinLabel(2822,"FGD2 RHC #nu_{#mu} CC0#pi");
  corr->GetYaxis()->SetBinLabel(2958,"FGD2 RHC #nu_{#mu} CC1#pi");
  corr->GetYaxis()->SetBinLabel(3010,"FGD2 RHC #nu_{#mu} CCOther");

  corr->Draw("colz");
  //  c->Print((std::string("detcorrFit")+std::string(".pdf")).c_str());
  c->SaveAs("detcorrFit.png");
 }
