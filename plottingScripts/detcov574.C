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

void detcov574(std::string file0) {

  TCanvas *c = new TCanvas("canv", "canv", 1580, 1380);
  TFile *file_ = TFile::Open(file0.c_str());
  
  c->SetTopMargin(0.05);
  c->SetBottomMargin(0.22);
  c->SetLeftMargin(0.21);
  c->SetRightMargin(0.2);

  TH2D* cov;
  cov = (TH2D*)file_->Get("Total_Covariance_Matrix");

  cov->SetTitle("");
  cov->GetXaxis()->SetLabelSize(0.03);
  cov->GetYaxis()->SetLabelSize(0.03);
  cov->GetXaxis()->SetNdivisions(1010);
  cov->GetZaxis()->SetRangeUser(-0.0125,0.0125);
  cov->GetZaxis()->SetTitle("Covariance");
  cov->GetZaxis()->SetTitleOffset(1.5);
  cov->GetXaxis()->SetRangeUser(0,573);
  cov->GetYaxis()->SetRangeUser(0,573);

  for (int i=1;i<574;i++)
    {
      cov->GetXaxis()->SetBinLabel(i,"");
      cov->GetYaxis()->SetBinLabel(i,"");
    }
  
  cov->GetXaxis()->SetBinLabel(1,"FGD1 FHC #nu_{#mu} CC0#pi");
  cov->GetXaxis()->SetBinLabel(65,"FGD1 FHC #nu_{#mu} CC1#pi");
  cov->GetXaxis()->SetBinLabel(146,"FGD1 FHC #nu_{#mu} CCOther");
  cov->GetXaxis()->SetBinLabel(236,"FGD2 FHC #nu_{#mu} CC0#pi");
  cov->GetXaxis()->SetBinLabel(300,"FGD2 FHC #nu_{#mu} CC1#pi");
  cov->GetXaxis()->SetBinLabel(381,"FGD2 FHC #nu_{#mu} CCOther");
  cov->GetXaxis()->SetBinLabel(471,"FGD1 RHC #bar{#nu_{#mu} CC0#pi}");
  cov->GetXaxis()->SetBinLabel(491,"FGD1 RHC #bar{#nu_{#mu} CC1#pi}");
  cov->GetXaxis()->SetBinLabel(495,"FGD1 RHC #bar{#nu_{#mu} CC1Other}");
  cov->GetXaxis()->SetBinLabel(507,"FGD2 RHC #bar{#nu_{#mu} CC0#pi}");
  cov->GetXaxis()->SetBinLabel(527,"FGD2 RHC #bar{#nu_{#mu} CC1#pi}");
  cov->GetXaxis()->SetBinLabel(531,"FGD2 RHC #bar{#nu_{#mu} CC1Other}");
  cov->GetXaxis()->SetBinLabel(543,"FGD1 RHC #nu_{#mu} CC0#pi");
  cov->GetXaxis()->SetBinLabel(549,"FGD1 RHC #nu_{#mu} CC1#pi");
  cov->GetXaxis()->SetBinLabel(555,"FGD1 RHC #nu_{#mu} CCOther");
  cov->GetXaxis()->SetBinLabel(559,"FGD2 RHC #nu_{#mu} CC0#pi");
  cov->GetXaxis()->SetBinLabel(565,"FGD2 RHC #nu_{#mu} CC1#pi");
  cov->GetXaxis()->SetBinLabel(571,"FGD2 RHC #nu_{#mu} CCOther");
    
  cov->GetYaxis()->SetBinLabel(1,"FGD1 FHC #nu_{#mu} CC0#pi");
  cov->GetYaxis()->SetBinLabel(65,"FGD1 FHC #nu_{#mu} CC1#pi");
  cov->GetYaxis()->SetBinLabel(146,"FGD1 FHC #nu_{#mu} CCOther");
  cov->GetYaxis()->SetBinLabel(236,"FGD2 FHC #nu_{#mu} CC0#pi");
  cov->GetYaxis()->SetBinLabel(300,"FGD2 FHC #nu_{#mu} CC1#pi");
  cov->GetYaxis()->SetBinLabel(381,"FGD2 FHC #nu_{#mu} CCOther");
  cov->GetYaxis()->SetBinLabel(471,"FGD1 RHC #bar{#nu_{#mu} CC0#pi}");
  cov->GetYaxis()->SetBinLabel(491,"FGD1 RHC #bar{#nu_{#mu} CC1#pi}");
  cov->GetYaxis()->SetBinLabel(495,"FGD1 RHC #bar{#nu_{#mu} CC1Other}");
  cov->GetYaxis()->SetBinLabel(507,"FGD2 RHC #bar{#nu_{#mu} CC0#pi}");
  cov->GetYaxis()->SetBinLabel(527,"FGD2 RHC #bar{#nu_{#mu} CC1#pi}");
  cov->GetYaxis()->SetBinLabel(531,"FGD2 RHC #bar{#nu_{#mu} CC1Other}");
  cov->GetYaxis()->SetBinLabel(543,"FGD1 RHC #nu_{#mu} CC0#pi");
  cov->GetYaxis()->SetBinLabel(549,"FGD1 RHC #nu_{#mu} CC1#pi");
  cov->GetYaxis()->SetBinLabel(555,"FGD1 RHC #nu_{#mu} CCOther");
  cov->GetYaxis()->SetBinLabel(559,"FGD2 RHC #nu_{#mu} CC0#pi");
  cov->GetYaxis()->SetBinLabel(565,"FGD2 RHC #nu_{#mu} CC1#pi");
  cov->GetYaxis()->SetBinLabel(571,"FGD2 RHC #nu_{#mu} CCOther");

  cov->Draw("colz");
  //  c->Print((std::string("detcov574")+std::string(".pdf")).c_str());
  c->SaveAs("detcov574.png");  

  cov->GetXaxis()->SetRangeUser(0,63);
  cov->GetYaxis()->SetRangeUser(0,63);
  cov->GetXaxis()->SetBinLabel(1,"0-300 MeV");
  cov->GetXaxis()->SetBinLabel(9,"300-1000 MeV");
  cov->GetXaxis()->SetBinLabel(17,"1000-1250 MeV");
  cov->GetXaxis()->SetBinLabel(25,"1250-1500 MeV");
  cov->GetXaxis()->SetBinLabel(33,"1500-2000 MeV");
  cov->GetXaxis()->SetBinLabel(41,"2000-3000 MeV");
  cov->GetXaxis()->SetBinLabel(49,"3000-5000 MeV");
  cov->GetXaxis()->SetBinLabel(57,"5000-30000 MeV");
  cov->GetYaxis()->SetBinLabel(1,"0-300 MeV");
  cov->GetYaxis()->SetBinLabel(9,"300-1000 MeV");
  cov->GetYaxis()->SetBinLabel(17,"1000-1250 MeV");
  cov->GetYaxis()->SetBinLabel(25,"1250-1500 MeV");
  cov->GetYaxis()->SetBinLabel(33,"1500-2000 MeV");
  cov->GetYaxis()->SetBinLabel(41,"2000-3000 MeV");
  cov->GetYaxis()->SetBinLabel(49,"3000-5000 MeV");
  cov->GetYaxis()->SetBinLabel(57,"5000-30000 MeV");

  cov->Draw("colz");
  c->SaveAs("detcov574_1sample.png");
 }
