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

void PrettifyTitles(TH2D*);

void xseccorr(std::string file0) {

  TCanvas *c = new TCanvas("canv", "canv", 1580, 1380);
  TFile *file_ = TFile::Open(file0.c_str());
  
  c->SetTopMargin(0.05);
  c->SetBottomMargin(0.18);
  c->SetLeftMargin(0.18);
  c->SetRightMargin(0.18);

  TH2D* cov = (TH2D*) file_->Get("hcov")->Clone();
 
  TH2D* corr = (TH2D*) file_->Get("hcov")->Clone();
  corr->GetZaxis()->SetTitle("Correlation");
  corr->GetZaxis()->SetTitleOffset(1.1);
  corr->GetZaxis()->SetRangeUser(-1,1);

  for(int i=1; i<= cov->GetXaxis()->GetNbins(); i++){
    for(int j=1; j<= cov->GetYaxis()->GetNbins();j++){
      corr->SetBinContent(i,j, cov->GetBinContent(i,j)/sqrt(cov->GetBinContent(i,i)*cov->GetBinContent(j,j)));
    }
  }

  PrettifyTitles(corr);
  corr->Draw("colz");
  c->Print((std::string("xseccorr")+std::string(".pdf")).c_str());
  
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

    else if (title == "FEFQE") title = "FSI Inel lo";
    else if (title == "FEFQEH") title = "FSI Inel hi";
    else if (title == "FEFINEL") title = "FSI Pi Prod.";
    else if (title == "FEFABS") title = "FSI Pi Abs.";
    else if (title == "FEFCX") title = "FSI Cex lo";
    else if (title == "FEFCXH") title = "FSI Cex hi";

    else if (title == "EB_dial_C_nu") title = "EB Dial C Nu";
    else if (title == "EB_dial_C_nubar") title = "EB Dial C Nubar";
    else if (title == "EB_dial_O_nu") title = "EB Dial O Nu";
    else if (title == "EB_dial_O_nubar") title = "EB Dial O Nubar";


    Hist->GetXaxis()->SetBinLabel(i+1, title.c_str());
    Hist->GetYaxis()->SetBinLabel(i+1, title.c_str());
  }

}
