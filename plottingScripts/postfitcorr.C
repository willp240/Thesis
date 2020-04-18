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

void postfitcorr(std::string file0) {

  TCanvas *c = new TCanvas("canv", "canv", 3000, 2550);
  TFile *file_ = TFile::Open(file0.c_str());
  
  c->SetTopMargin(0.02);
  c->SetBottomMargin(0.1);
  c->SetLeftMargin(0.08);
  c->SetRightMargin(0.16);

  TH2D* cov2 = (TH2D*)file_->Get("postfit_cov_plot");//new TH2D("","",100,0,100,100,0,100);
  
  for (int i=1;i<cov2->GetXaxis()->GetNbins();i++)
    {
      for(int j=1; j<cov2->GetXaxis()->GetNbins(); j++){
	if(cov2->GetBinContent(i,j)>0)
	  cov2->SetBinContent(i,j,sqrt(cov2->GetBinContent(i,j)));
	else
	  cov2->SetBinContent(i,j,-sqrt(-cov2->GetBinContent(i,j)));
      }
    }
  
  cov2->GetXaxis()->SetBinLabel(1,"ND280 FHC #nu_{#mu}");//11
  cov2->GetXaxis()->SetBinLabel(12,"ND280 FHC #nu_{e}");//7
  cov2->GetXaxis()->SetBinLabel(19,"ND280 FHC #bar{#nu_{#mu}}");//5
  cov2->GetXaxis()->SetBinLabel(24,"ND280 FHC #bar{#nu_{e}}");//2
  cov2->GetXaxis()->SetBinLabel(26,"ND280 RHC #bar{#nu_{#mu}}");//11
  cov2->GetXaxis()->SetBinLabel(37,"ND280 RHC #bar{#nu_{e}}");//7
  cov2->GetXaxis()->SetBinLabel(44,"ND280 RHC #nu_{#mu}");//5
  cov2->GetXaxis()->SetBinLabel(49,"ND280 RHC #nu_{e}");//2
  cov2->GetXaxis()->SetBinLabel(51,"SK FHC #nu_{#mu}");//11
  cov2->GetXaxis()->SetBinLabel(62,"SK FHC #nu_{e}");//7
  cov2->GetXaxis()->SetBinLabel(69,"SK FHC #bar{#nu_{#mu}}");//5
  cov2->GetXaxis()->SetBinLabel(74,"SK FHC #bar{#nu_{e}}");//2
  cov2->GetXaxis()->SetBinLabel(76,"SK RHC #bar{#nu_{#mu}}");//11
  cov2->GetXaxis()->SetBinLabel(87,"SK RHC #bar{#nu_{e}}");//7
  cov2->GetXaxis()->SetBinLabel(94,"SK RHC #nu_{#mu}");//5
  cov2->GetXaxis()->SetBinLabel(99,"SK RHC #nu_{e}");//2

  cov2->GetYaxis()->SetBinLabel(1,"ND280 FHC #nu_{#mu}");//11
  cov2->GetYaxis()->SetBinLabel(12,"ND280 FHC #nu_{e}");//7
  cov2->GetYaxis()->SetBinLabel(19,"ND280 FHC #bar{#nu_{#mu}}");//5
  cov2->GetYaxis()->SetBinLabel(24,"ND280 FHC #bar{#nu_{e}}");//2
  cov2->GetYaxis()->SetBinLabel(26,"ND280 RHC #bar{#nu_{#mu}}");//11
  cov2->GetYaxis()->SetBinLabel(37,"ND280 RHC #bar{#nu_{e}}");//7
  cov2->GetYaxis()->SetBinLabel(44,"ND280 RHC #nu_{#mu}");//5
  cov2->GetYaxis()->SetBinLabel(49,"ND280 RHC #nu_{e}");//2
  cov2->GetYaxis()->SetBinLabel(51,"SK FHC #nu_{#mu}");//11
  cov2->GetYaxis()->SetBinLabel(62,"SK FHC #nu_{e}");//7
  cov2->GetYaxis()->SetBinLabel(69,"SK FHC #bar{#nu_{#mu}}");//5
  cov2->GetYaxis()->SetBinLabel(74,"SK FHC #bar{#nu_{e}}");//2
  cov2->GetYaxis()->SetBinLabel(76,"SK RHC #bar{#nu_{#mu}}");//11
  cov2->GetYaxis()->SetBinLabel(87,"SK RHC #bar{#nu_{e}}");//7
  cov2->GetYaxis()->SetBinLabel(94,"SK RHC #nu_{#mu}");//5
  cov2->GetYaxis()->SetBinLabel(99,"SK RHC #nu_{e}");//2

  PrettifyTitles(cov2);

  cov2->SetTitle("");
  cov2->GetXaxis()->SetLabelSize(0.01);
  cov2->GetYaxis()->SetLabelSize(0.01);
  cov2->GetXaxis()->SetNdivisions(1010);
  cov2->GetZaxis()->SetRangeUser(-0.2,0.2);
  cov2->GetZaxis()->SetTitle("Covariance");
  cov2->GetZaxis()->SetTitleOffset(1.5);
  cov2->Draw("colz");
  c->Print((std::string("postfitcov")+std::string(".pdf")).c_str());  

}

void PrettifyTitles(TH2D *Hist) {

  for (int i = 0; i < Hist->GetXaxis()->GetNbins(); ++i) {
    std::string title = Hist->GetXaxis()->GetBinLabel(i+1);
    if (title == "MAQE") title = "M_{A}^{QE}";

    else if (title == "2p2h_norm_nu") title = "2p2h norm #nu";
    else if (title == "2p2h_norm_nubar") title = "2p2h norm #bar{#nu}";
    else if (title == "2p2h_normCtoO") title = "2p2h norm ^{12}C/^{16}O";

    else if (title == "2p2h_shape_C") title = "2p2h shape ^{12}C";
    else if (title == "2p2h_shape_O") title = "2p2h shape ^{16}O";

    else if (title == "2p2h_Edep_lowEnu") title = "2p2h Edep Low #nu";
    else if (title == "2p2h_Edep_highEnu") title = "2p2h Edep High #nu";
    else if (title == "2p2h_Edep_lowEnubar") title = "2p2h Edep Low #bar{#nu}";
    else if (title == "2p2h_Edep_highEnubar") title = "2p2h Edep High #bar{#nu}";

    else if (title == "Q2_norm_0") title = "Q^{2} Norm 0";
    else if (title == "Q2_norm_1") title = "Q^{2} Norm 1";
    else if (title == "Q2_norm_2") title = "Q^{2} Norm 2";
    else if (title == "Q2_norm_3") title = "Q^{2} Norm 3";
    else if (title == "Q2_norm_4") title = "Q^{2} Norm 4";
    else if (title == "Q2_norm_5") title = "Q^{2} Norm 5";
    else if (title == "Q2_norm_6") title = "Q^{2} Norm 6";
    else if (title == "Q2_norm_7") title = "Q^{2} Norm 7";

    else if (title == "CA5") title = "C_{5}^{A}";
    else if (title == "MARES") title = "M_{A}^{RES}";
    else if (title == "ISO_BKG") title = "Non-res I_{1/2}";
    else if (title == "ISO_BKG_LowPPi") title = "Non-res I_{1/2} Low p_{#pi}";

    else if (title == "CC_norm_nu") title = "CC Norm #nu";
    else if (title == "CC_norm_nubar") title = "CC Norm #bar{nu}";
    else if (title == "nue_numu") title = "#nu_{e}/#nu_{#mu}";
    else if (title == "nuebar_numubar") title = "#bar{#nu}_{e}/#bar{#nu}_{#mu}";

    else if (title == "CC_BY_DIS") title = "CC BY DIS";
    else if (title == "CC_BY_MPi") title = "CC BY Multi #pi";
    else if (title == "CC_AGKY_Mult") title = "CC AGKY Multi #pi";
    else if (title == "CC_Misc") title = "CC Misc";
    else if (title == "CC_BY_MPi") title = "CC BY Multi #pi";
    else if (title == "CC_DIS_MultPi_Norm_Nu") title = "CC DIS/M#pi norm #nu";
    else if (title == "CC_DIS_MultPi_Norm_Nubar") title = "CC DIS/M#pi norm #bar{#nu}";

    else if (title == "CC_Coh_C") title = "CC coh. ^{12}C";
    else if (title == "CC_Coh_O") title = "CC coh. ^{16}O";
    else if (title == "NC_Coh") title = "NC coh.";
    else if (title == "NC_1gamma") title = "NC 1#gamma";
    else if (title == "NC_other_near") title = "NC oth. ND280";
    else if (title == "NC_other_far") title = "NC oth. SK";

    else if (title == "FSI_INEL_LO") title = "FSI Inel lo";
    else if (title == "FSI_INEL_HI") title = "FSI Inel hi";
    else if (title == "FSI_PI_PROD") title = "FSI Pi Prod.";
    else if (title == "FSI_PI_ABS") title = "FSI Pi Abs.";
    else if (title == "FSI_CEX_LO") title = "FSI Cex lo";
    else if (title == "FSI_CEX_HI") title = "FSI Cex hi";

    else if (title == "EB_dial_C_nu") title = "EB Dial C Nu";
    else if (title == "EB_dial_C_nubar") title = "EB Dial C Nubar";
    else if (title == "EB_dial_O_nu") title = "EB Dial O Nu";
    else if (title == "EB_dial_O_nubar") title = "EB Dial O Nubar";


    Hist->GetXaxis()->SetBinLabel(i+1, title.c_str());
    Hist->GetYaxis()->SetBinLabel(i+1, title.c_str());
  }

}
