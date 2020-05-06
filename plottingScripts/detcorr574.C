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

void detcorr574(std::string file0) {

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
  corr->GetXaxis()->SetRangeUser(0,573);
  corr->GetYaxis()->SetRangeUser(0,573);

  for (int i=1;i<574;i++)
    {
      corr->GetXaxis()->SetBinLabel(i,"");
      corr->GetYaxis()->SetBinLabel(i,"");
    }
  
  corr->GetXaxis()->SetBinLabel(1,"FGD1 FHC #nu_{#mu} CC0#pi");
  corr->GetXaxis()->SetBinLabel(65,"FGD1 FHC #nu_{#mu} CC1#pi");
  corr->GetXaxis()->SetBinLabel(146,"FGD1 FHC #nu_{#mu} CCOther");
  corr->GetXaxis()->SetBinLabel(236,"FGD2 FHC #nu_{#mu} CC0#pi");
  corr->GetXaxis()->SetBinLabel(300,"FGD2 FHC #nu_{#mu} CC1#pi");
  corr->GetXaxis()->SetBinLabel(381,"FGD2 FHC #nu_{#mu} CCOther");
  corr->GetXaxis()->SetBinLabel(471,"FGD1 RHC #bar{#nu_{#mu} CC0#pi}");
  corr->GetXaxis()->SetBinLabel(491,"FGD1 RHC #bar{#nu_{#mu} CC1#pi}");
  corr->GetXaxis()->SetBinLabel(495,"FGD1 RHC #bar{#nu_{#mu} CC1Other}");
  corr->GetXaxis()->SetBinLabel(507,"FGD2 RHC #bar{#nu_{#mu} CC0#pi}");
  corr->GetXaxis()->SetBinLabel(527,"FGD2 RHC #bar{#nu_{#mu} CC1#pi}");
  corr->GetXaxis()->SetBinLabel(531,"FGD2 RHC #bar{#nu_{#mu} CC1Other}");
  corr->GetXaxis()->SetBinLabel(543,"FGD1 RHC #nu_{#mu} CC0#pi");
  corr->GetXaxis()->SetBinLabel(549,"FGD1 RHC #nu_{#mu} CC1#pi");
  corr->GetXaxis()->SetBinLabel(555,"FGD1 RHC #nu_{#mu} CCOther");
  corr->GetXaxis()->SetBinLabel(559,"FGD2 RHC #nu_{#mu} CC0#pi");
  corr->GetXaxis()->SetBinLabel(565,"FGD2 RHC #nu_{#mu} CC1#pi");
  corr->GetXaxis()->SetBinLabel(571,"FGD2 RHC #nu_{#mu} CCOther");
    
  corr->GetYaxis()->SetBinLabel(1,"FGD1 FHC #nu_{#mu} CC0#pi");
  corr->GetYaxis()->SetBinLabel(65,"FGD1 FHC #nu_{#mu} CC1#pi");
  corr->GetYaxis()->SetBinLabel(146,"FGD1 FHC #nu_{#mu} CCOther");
  corr->GetYaxis()->SetBinLabel(236,"FGD2 FHC #nu_{#mu} CC0#pi");
  corr->GetYaxis()->SetBinLabel(300,"FGD2 FHC #nu_{#mu} CC1#pi");
  corr->GetYaxis()->SetBinLabel(381,"FGD2 FHC #nu_{#mu} CCOther");
  corr->GetYaxis()->SetBinLabel(471,"FGD1 RHC #bar{#nu_{#mu} CC0#pi}");
  corr->GetYaxis()->SetBinLabel(491,"FGD1 RHC #bar{#nu_{#mu} CC1#pi}");
  corr->GetYaxis()->SetBinLabel(495,"FGD1 RHC #bar{#nu_{#mu} CC1Other}");
  corr->GetYaxis()->SetBinLabel(507,"FGD2 RHC #bar{#nu_{#mu} CC0#pi}");
  corr->GetYaxis()->SetBinLabel(527,"FGD2 RHC #bar{#nu_{#mu} CC1#pi}");
  corr->GetYaxis()->SetBinLabel(531,"FGD2 RHC #bar{#nu_{#mu} CC1Other}");
  corr->GetYaxis()->SetBinLabel(543,"FGD1 RHC #nu_{#mu} CC0#pi");
  corr->GetYaxis()->SetBinLabel(549,"FGD1 RHC #nu_{#mu} CC1#pi");
  corr->GetYaxis()->SetBinLabel(555,"FGD1 RHC #nu_{#mu} CCOther");
  corr->GetYaxis()->SetBinLabel(559,"FGD2 RHC #nu_{#mu} CC0#pi");
  corr->GetYaxis()->SetBinLabel(565,"FGD2 RHC #nu_{#mu} CC1#pi");
  corr->GetYaxis()->SetBinLabel(571,"FGD2 RHC #nu_{#mu} CCOther");

  corr->Draw("colz");
  //  c->Print((std::string("detcov574")+std::string(".pdf")).c_str());
  c->SaveAs("detcorr574.png");  
 }
