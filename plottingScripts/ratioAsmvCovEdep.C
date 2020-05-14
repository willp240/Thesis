#include "TH1D.h"
#include "TH2D.h"
#include "TMatrixT.h"
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

int MaCh3ToBanff(int mach3index);
void PrettifyTitles(TH2D*);

void ratioAsmvCovEdep() {

  TCanvas *c = new TCanvas("canv", "canv", 3000, 2550);
  c->SetTopMargin(0.02);
  c->SetBottomMargin(0.1);
  c->SetLeftMargin(0.08);
  c->SetRightMargin(0.16);

  TCanvas *c1 = new TCanvas("canv1", "canv1", 3000, 2550);
  c1->SetTopMargin(0.02);
  c1->SetBottomMargin(0.1);
  c1->SetLeftMargin(0.08);
  c1->SetRightMargin(0.16);

  TFile *m3file = TFile::Open("/vols/t2k/users/wparker2/OA2019/Cedar/27AprTH2D574Data/27AprTH2D574Data_drawCorr.root");
  TFile *asmvfile = TFile::Open("inputFiles/27AprTH2D574Asimov_drawCorr.root");  

  TObjArray* param_list = (TObjArray*)asmvfile->Get("param_list");

  TH2D* asmvcov = (TH2D*)asmvfile->Get("postfit_cov_plot")->Clone();
  TH2D* asmvcov2 = (TH2D*)asmvfile->Get("postfit_cov_plot")->Clone();
  TH2D* mach3cov = (TH2D*)m3file->Get("postfit_cov_plot")->Clone();
  TH2D* mach3cov2 = (TH2D*)m3file->Get("postfit_cov_plot")->Clone();
  TH2D* ratio = (TH2D*) m3file->Get("postfit_cov_plot")->Clone();

  for(int i=1; i<=mach3cov->GetXaxis()->GetNbins(); i++){
    for(int j=1; j<=mach3cov->GetYaxis()->GetNbins(); j++){

      double m3 = mach3cov->GetBinContent(i,j)/(sqrt(mach3cov->GetBinContent(i,i)*mach3cov->GetBinContent(j,j)));
      double asmv = asmvcov->GetBinContent(i,j)/(sqrt(asmvcov->GetBinContent(i,i)*asmvcov->GetBinContent(j,j)));

      mach3cov2->SetBinContent(i,j,m3);
      asmvcov2->SetBinContent(i,j,asmv);
      ratio->SetBinContent(i,j,(m3-asmv));

      /*      if(i==125||i==142||j==125||j==142){
	mach3cov2->SetBinContent(i,j,0);
	asmvcov2->SetBinContent(i,j,0);
	ratio->SetBinContent(i,j,0);

	if(i==j){
	  mach3cov2->SetBinContent(i,j,1);
	  asmvcov2->SetBinContent(i,j,1);
	  ratio->SetBinContent(i,j,0);
	}
	}*/
    }
  }
  std::cout << "set values" << std::endl;
  ratio->SetTitle("");
  ratio->GetZaxis()->SetTitle("Difference in Correlations");
  ratio->GetZaxis()->SetTitleOffset(1.1);
  ratio->GetZaxis()->SetRangeUser(-1,1);
  PrettifyTitles(ratio);
  ratio->Draw("colz");
  c1->Print("MaCh3AsmvDataCorrDifference.png");
  ratio->GetXaxis()->SetRangeUser(100,147);
  ratio->GetYaxis()->SetRangeUser(100,147);
  ratio->Draw("colz");
  c1->Print("MaCh3AsmvDataCorrDifferenceXsec.png");

  asmvcov2->SetTitle("");
  asmvcov2->GetZaxis()->SetTitle("Correlation");
  asmvcov2->GetZaxis()->SetRangeUser(-1.0,1.0);
  PrettifyTitles(asmvcov2);
  asmvcov2->Draw("colz");
  c1->Print("Mach3AsmvCorr.png");
  asmvcov2->GetXaxis()->SetRangeUser(100,147);
  asmvcov2->GetYaxis()->SetRangeUser(100,147);
  asmvcov2->Draw("colz");
  c1->Print("Mach3AsmvCorrXsec.png");
  c->cd();/*
  mach3cov2->SetTitle("");
  mach3cov2->GetZaxis()->SetTitle("Correlation");
  mach3cov2->GetZaxis()->SetRangeUser(-1.0,1.0);
  PrettifyTitles(mach3cov2);
  mach3cov2->Draw("colz");
  c->Print("mach3.png");
  mach3cov2->GetXaxis()->SetRangeUser(100,147);
  mach3cov2->GetYaxis()->SetRangeUser(100,147);
  mach3cov2->Draw("colz");
  c->Print("mach3xsec.png");
  c->cd();*/
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

  Hist->GetXaxis()->SetBinLabel(1,"ND280 FHC #nu_{#mu}");//11
  Hist->GetXaxis()->SetBinLabel(12,"ND280 FHC #bar{#nu_{#mu}}");//5
  Hist->GetXaxis()->SetBinLabel(17,"ND280 FHC #nu_{e}");//7
  Hist->GetXaxis()->SetBinLabel(24,"ND280 FHC #bar{#nu_{e}}");//2
  Hist->GetXaxis()->SetBinLabel(26,"ND280 RHC #nu_{#mu}");//5
  Hist->GetXaxis()->SetBinLabel(31,"ND280 RHC #bar{#nu_{#mu}}");//11
  Hist->GetXaxis()->SetBinLabel(42,"ND280 RHC #nu_{e}");//2
  Hist->GetXaxis()->SetBinLabel(44,"ND280 RHC #bar{#nu_{e}}");//7
  Hist->GetXaxis()->SetBinLabel(51,"SK FHC #nu_{#mu}");//11
  Hist->GetXaxis()->SetBinLabel(62,"SK FHC #bar{#nu_{#mu}}");//5
  Hist->GetXaxis()->SetBinLabel(67,"SK FHC #nu_{e}");//7
  Hist->GetXaxis()->SetBinLabel(74,"SK FHC #bar{#nu_{e}}");//2
  Hist->GetXaxis()->SetBinLabel(76,"SK RHC #nu_{#mu}");//5
  Hist->GetXaxis()->SetBinLabel(81,"SK RHC #bar{#nu_{#mu}}");//11
  Hist->GetXaxis()->SetBinLabel(92,"SK RHC #nu_{e}");//2
  Hist->GetXaxis()->SetBinLabel(94,"SK RHC #bar{#nu_{e}}");//7

  Hist->GetYaxis()->SetBinLabel(1,"ND280 FHC #nu_{#mu}");//11
  Hist->GetYaxis()->SetBinLabel(12,"ND280 FHC #bar{#nu_{#mu}}");//5
  Hist->GetYaxis()->SetBinLabel(17,"ND280 FHC #nu_{e}");//7
  Hist->GetYaxis()->SetBinLabel(24,"ND280 FHC #bar{#nu_{e}}");//2
  Hist->GetYaxis()->SetBinLabel(26,"ND280 RHC #nu_{#mu}");//5
  Hist->GetYaxis()->SetBinLabel(31,"ND280 RHC #bar{#nu_{#mu}}");//11
  Hist->GetYaxis()->SetBinLabel(42,"ND280 RHC #nu_{e}");//2
  Hist->GetYaxis()->SetBinLabel(44,"ND280 RHC #bar{#nu_{e}}");//7
  Hist->GetYaxis()->SetBinLabel(51,"SK FHC #nu_{#mu}");//11
  Hist->GetYaxis()->SetBinLabel(62,"SK FHC #bar{#nu_{#mu}}");//5
  Hist->GetYaxis()->SetBinLabel(67,"SK FHC #nu_{e}");//7
  Hist->GetYaxis()->SetBinLabel(74,"SK FHC #bar{#nu_{e}}");//2
  Hist->GetYaxis()->SetBinLabel(76,"SK RHC #nu_{#mu}");//5
  Hist->GetYaxis()->SetBinLabel(81,"SK RHC #bar{#nu_{#mu}}");//11
  Hist->GetYaxis()->SetBinLabel(92,"SK RHC #nu_{e}");//2
  Hist->GetYaxis()->SetBinLabel(94,"SK RHC #bar{#nu_{e}}");//7
  
}

int MaCh3ToBanff(int mach3index){

  int banffindex;

  if(mach3index<100)
    banffindex=mach3index;

  else if(mach3index>=100 && mach3index<118)
    banffindex=mach3index+574+5;

  else if(mach3index>=118 && mach3index<122)
    banffindex=mach3index+574+5+20;

  else if(mach3index>=122 &&mach3index<142)
    banffindex=mach3index+574+5-4;

  else if(mach3index>=142)
    banffindex=mach3index+574-42;

  else{
    std::cout << "Hang on fella, " << mach3index << "not in expected range" << std::endl;
    return 0;
  }

  return banffindex;

}
