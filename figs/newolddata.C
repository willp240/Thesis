//Use on output of DrawComp to get different plotting style of postfit parameters
//Run with something like root -l -b -q 'GetAsimovPlots_HPDOnly.C("outputFromDrawComp.root")'
//Originally written by Clarence, with some changes by Will
//Central postfit value taken is the Highest Posterior Density. Watch out for parameter names and number of parameters per plot being quite hardcoded

#include <algorithm>
// Rebin the plots
const int nBinsSameSignMu = 11;
const int nBinsSameSignE = 7;
const int nBinsWrongSignMu = 5;
const int nBinsWrongSignE = 2;

double SameSignMu[nBinsSameSignMu+1]    = {0.3, 0.4, 0.5, 0.6, 0.7, 1.0, 1.5, 2.5, 3.5, 5.0, 7.0, 30.0};
double SameSignE[nBinsSameSignE+1]      = {0.3, 0.5, 0.7, 0.8, 1.5, 2.5, 4.0, 30.0};
double WrongSignMu[nBinsWrongSignMu+1]  = {0.3, 0.7, 1.0, 1.5, 2.5, 30.0};
double WrongSignE[nBinsWrongSignE+1]    = {0.3, 2.5, 30.0};

TH1D* PrettifyMe(TH1D *One, int which) {

  TH1D* OneCopy = NULL;
  std::string Name = One->GetName();

  int offset = -1;
  switch (which) {
    case 0:
      OneCopy = new TH1D((Name+"_ND280_FHC_numu").c_str(), "ND280 FHC #nu_{#mu}", nBinsSameSignMu, SameSignMu);
      offset = 0;
      break;
    case 1:
      OneCopy = new TH1D((Name+"_ND280_FHC_nue").c_str(), "ND280 FHC #nu_{e}", nBinsSameSignE, SameSignE);
      offset = nBinsSameSignMu;
      break;
    case 2:
      OneCopy = new TH1D((Name+"_ND280_FHC_numubar").c_str(), "ND280 FHC #bar{#nu}_{#mu}", nBinsWrongSignMu, WrongSignMu);
      offset = nBinsSameSignE+nBinsSameSignMu;
      break;
    case 3:
      OneCopy = new TH1D((Name+"_ND280_FHC_nuebar").c_str(), "ND280 FHC #bar{#nu}_{e}", nBinsWrongSignE, WrongSignE);
      offset = nBinsWrongSignMu+nBinsSameSignE+nBinsSameSignMu;
      break;

    case 4:
      OneCopy = new TH1D((Name+"_ND280_RHC_numubar").c_str(), "ND280 RHC #bar{#nu}_{#mu}", nBinsSameSignMu, SameSignMu);
      offset = nBinsWrongSignE+nBinsWrongSignMu+nBinsSameSignE+nBinsSameSignMu;
      break;
    case 5:
      OneCopy = new TH1D((Name+"_ND280_RHC_nuebar").c_str(), "ND280 RHC #bar{#nu}_{e}", nBinsSameSignE, SameSignE);
      offset = nBinsWrongSignE+nBinsWrongSignMu+nBinsSameSignE+nBinsSameSignMu*2;
      break;
    case 6:
      OneCopy = new TH1D((Name+"_ND280_RHC_numu").c_str(), "ND280 RHC #nu_{#mu}", nBinsWrongSignMu, WrongSignMu);
      offset = nBinsWrongSignE+nBinsWrongSignMu+nBinsSameSignE*2+nBinsSameSignMu*2;
      break;
    case 7:
      OneCopy = new TH1D((Name+"_ND280_RHC_nue").c_str(), "ND280 RHC #nu_{e}", nBinsWrongSignE, WrongSignE);
      offset = nBinsWrongSignE+nBinsWrongSignMu*2+nBinsSameSignE*2+nBinsSameSignMu*2;
      break;


    case 8:
      OneCopy = new TH1D((Name+"_SK_FHC_numu").c_str(), "SK FHC #nu_{#mu}", nBinsSameSignMu, SameSignMu);
      offset = nBinsWrongSignE*2+nBinsWrongSignMu*2+nBinsSameSignE*2+nBinsSameSignMu*2;
      break;
    case 9:
      OneCopy = new TH1D((Name+"_SK_FHC_nue").c_str(), "SK FHC #nu_{e}", nBinsSameSignE, SameSignE);
      offset = nBinsWrongSignE*2+nBinsWrongSignMu*2+nBinsSameSignE*2+nBinsSameSignMu*3;
      break;
    case 10:
      OneCopy = new TH1D((Name+"_SK_FHC_numubar").c_str(), "SK FHC #bar{#nu}_{#mu}", nBinsWrongSignMu, WrongSignMu);
      offset = nBinsWrongSignE*2+nBinsWrongSignMu*2+nBinsSameSignE*3+nBinsSameSignMu*3;
      break;
    case 11:
      OneCopy = new TH1D((Name+"_SK_FHC_nuebar").c_str(), "SK FHC #bar{#nu}_{e}", nBinsWrongSignE, WrongSignE);
      offset = nBinsWrongSignE*2+nBinsWrongSignMu*3+nBinsSameSignE*3+nBinsSameSignMu*3;
      break;

    case 12:
      OneCopy = new TH1D((Name+"_SK_RHC_numubar").c_str(), "SK RHC #bar{#nu}_{#mu}", nBinsSameSignMu, SameSignMu);
      offset = nBinsWrongSignE*3+nBinsWrongSignMu*3+nBinsSameSignE*3+nBinsSameSignMu*3;
      break;
    case 13:
      OneCopy = new TH1D((Name+"_SK_RHC_nuebar").c_str(), "SK RHC #bar{#nu}_{e}", nBinsSameSignE, SameSignE);
      offset = nBinsWrongSignE*3+nBinsWrongSignMu*3+nBinsSameSignE*3+nBinsSameSignMu*4;
      break;
    case 14:
      OneCopy = new TH1D((Name+"_SK_RHC_numu").c_str(), "SK RHC #nu_{#mu}", nBinsWrongSignMu, WrongSignMu);
      offset = nBinsWrongSignE*3+nBinsWrongSignMu*3+nBinsSameSignE*4+nBinsSameSignMu*4;
      break;
    case 15:
      OneCopy = new TH1D((Name+"_SK_RHC_nue").c_str(), "SK RHC #nu_{e}", nBinsWrongSignE, WrongSignE);
      offset = nBinsWrongSignE*3+nBinsWrongSignMu*4+nBinsSameSignE*4+nBinsSameSignMu*4;
      break;
  }
  OneCopy->GetYaxis()->SetTitle("Variation rel. nom.");
  OneCopy->GetXaxis()->SetTitle("E_{#nu} (GeV)");
  OneCopy->GetXaxis()->SetTitleOffset(OneCopy->GetXaxis()->GetTitleOffset()*1.2);

  OneCopy->SetFillStyle(One->GetFillStyle());
  OneCopy->SetFillColor(One->GetFillColor());
  OneCopy->SetMarkerStyle(One->GetMarkerStyle());
  OneCopy->SetMarkerSize(One->GetMarkerSize());
  OneCopy->SetMarkerColor(One->GetMarkerColor());
  OneCopy->SetLineColor(One->GetLineColor());
  OneCopy->GetYaxis()->SetRangeUser(0.7, 1.3);
  //OneCopy->GetXaxis()->SetMoreLogLabels();

  for (int j = 0; j < OneCopy->GetXaxis()->GetNbins(); ++j) {
    OneCopy->SetBinContent(j+1, One->GetBinContent(j+1+offset));
    OneCopy->SetBinError(j+1, One->GetBinError(j+1+offset));
  }

  return OneCopy;
}

void PrettifyTitles(TH1D *Hist) {
  
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
  }

}

TH1D* MakeRatioPlot(TH1D* OneCopy, TH1D *TwoCopy) {

  TH1D *Ratio = (TH1D*)OneCopy->Clone();
  Ratio->GetYaxis()->SetTitle("(x_{fit}-#mu_{Prior})/#sigma_{Prior}");
  Ratio->SetMinimum(-3.7);
  Ratio->SetMaximum(3.7);
  for (int j = 0; j < Ratio->GetXaxis()->GetNbins(); ++j) {
    if (TwoCopy->GetBinError(j+1) != 0) {
      Ratio->SetBinContent(j+1, (TwoCopy->GetBinContent(j+1)-OneCopy->GetBinContent(j+1))/OneCopy->GetBinError(j+1));
    } else {
      Ratio->SetBinContent(j+1, 0);
    }
    double up = (TwoCopy->GetBinContent(j+1)+TwoCopy->GetBinError(j+1)-OneCopy->GetBinContent(j+1))/OneCopy->GetBinError(j+1);
    double down = (TwoCopy->GetBinContent(j+1)-TwoCopy->GetBinError(j+1)-OneCopy->GetBinContent(j+1))/OneCopy->GetBinError(j+1);

    double maximum = up-Ratio->GetBinContent(j+1);
    double minimum = Ratio->GetBinContent(j+1)-down;

    Ratio->SetBinError(j+1, std::max(maximum, minimum));
  }

  Ratio->SetFillStyle(0);
  Ratio->SetFillColor(0);

  Ratio->SetLineColor(TwoCopy->GetFillColor());
  if (Ratio->GetLineColor() == 0) Ratio->SetLineColor(kBlack);
  Ratio->SetLineWidth(2);
  Ratio->SetTitle("");

  Ratio->SetMarkerSize(2);
  Ratio->SetMarkerStyle(20);

  Ratio->GetYaxis()->SetTitleSize(30);
  Ratio->GetYaxis()->SetTitleFont(43);
  Ratio->GetYaxis()->SetTitleOffset(2.0);
  Ratio->GetYaxis()->SetLabelFont(43);
  Ratio->GetYaxis()->SetLabelSize(25);
  Ratio->GetYaxis()->CenterTitle();
  Ratio->GetYaxis()->SetNdivisions(5,2,0);

  Ratio->GetXaxis()->SetTitleSize(30);
  Ratio->GetXaxis()->SetTitleFont(43);
  Ratio->GetXaxis()->SetTitleOffset(4.0);
  Ratio->GetXaxis()->SetLabelFont(43);
  Ratio->GetXaxis()->SetLabelSize(25);

  return Ratio;
}

void newolddata(std::string FileName1, std::string FileName2, std::string FileName3) {
  
  TFile *File1 = new TFile(FileName1.c_str());
  TFile *File2 = new TFile(FileName2.c_str());
  TFile *File3 = new TFile(FileName3.c_str());

  //gStyle->SetPalette(51);
  TCanvas *canv = new TCanvas("canv", "canv", 1024, 1024);
  canv->SetLeftMargin(0.12);
  canv->SetBottomMargin(0.12);
  canv->SetTopMargin(0.08);
  canv->SetRightMargin(0.04);
  //canv->Print((FileName1+"_"+FileName2+"_"+FileName3+".pdf[").c_str());
  TPad *p1 = new TPad("p1", "p1", 0.0, 0.3, 1.0, 1.0);
  TPad *p2 = new TPad("p2", "p2", 0.0, 0.0, 1.0, 0.3);
  p1->SetLeftMargin(canv->GetLeftMargin());
  p1->SetRightMargin(canv->GetRightMargin());
  p1->SetTopMargin(canv->GetTopMargin());
  p1->SetBottomMargin(0);
  p2->SetLeftMargin(canv->GetLeftMargin());
  p2->SetRightMargin(canv->GetRightMargin());
  p2->SetTopMargin(0);
  p2->SetBottomMargin(0.25);


  TH1D *One = (TH1D*)File1->Get("param_flux_prefit_0");
  TH1D *Two = (TH1D*)File1->Get("param_flux_HPD_0");
  TH1D *Three = (TH1D*)File2->Get("param_flux_HPD_0");
  TH1D *Four = (TH1D*)File3->Get("param_flux_HPD_0");

  One->SetFillColor(kRed);
  One->SetFillStyle(3003);
  One->SetMarkerStyle(7);
  One->SetMarkerColor(kRed);
  One->SetLineColor(kRed);
  Two->SetFillColor(kBlue);
  Two->SetFillStyle(3003);
  Two->SetMarkerColor(kBlue);
  Two->SetMarkerStyle(7);
  Two->SetLineColor(kBlue);
  Three->SetMarkerColor(kBlack);
  Three->SetLineColor(kBlack);
  Three->SetMarkerStyle(7);
  Four->SetMarkerColor(kGreen+3);
  Four->SetLineColor(kGreen+3);
  Four->SetFillColor(kGreen+3);
  Four->SetMarkerStyle(7);

  // Make a Legend page
  TLegend *leg = new TLegend(0.0, 0.0, 1.0, 1.0);
  if (One != NULL) leg->AddEntry(One, "Prior", "lpf");
  if (Two != NULL) leg->AddEntry(Two, "Runs 2-9", "lpf");
  if (Three != NULL) leg->AddEntry(Three, "Runs 2-6", "lpf");
  if (Four != NULL) leg->AddEntry(Four, "Runs 7-9", "lp");

  canv->cd();
  canv->Clear();
  leg->Draw();
  canv->Print("newolddatafits_leg.pdf");
  delete leg;

  int nBins = One->GetXaxis()->GetNbins();

  // Make eight flux plots
  const int nFluxPlots = 16;
  TH1D* OneCopy[nFluxPlots];
  TH1D* TwoCopy[nFluxPlots];
  TH1D* ThreeCopy[nFluxPlots];
  TH1D* FourCopy[nFluxPlots];

  int nFluxParams = 0;
  for (int i = 0; i < nFluxPlots; ++i) {

    OneCopy[i]  = PrettifyMe(One, i);
    if (Two != NULL) TwoCopy[i]  = PrettifyMe(Two, i);
    if (Three != NULL) ThreeCopy[i]  = PrettifyMe(Three, i);    
    if (Four != NULL) FourCopy[i] = PrettifyMe(Four, i);

    canv->cd();
    p1->Draw();
    p1->cd();
    p1->SetLogx(true);
    OneCopy[i]->Draw("e2");
    if (Two != NULL) TwoCopy[i]->Draw("e2,same");
    if (Three != NULL) ThreeCopy[i]->Draw("e1,same");
    if (Four != NULL) FourCopy[i]->Draw("e1,same");

    //OneCopy[i]->GetXaxis()->SetMoreLogLabels();
    OneCopy[i]->GetYaxis()->SetLabelSize(0.);
    OneCopy[i]->GetYaxis()->SetTitleSize(0.05);
    OneCopy[i]->GetYaxis()->SetTitleOffset(1.3);
    canv->Update();
    TGaxis *axis = new TGaxis(OneCopy[i]->GetXaxis()->GetBinLowEdge(OneCopy[i]->GetXaxis()->GetFirst()), gPad->GetUymin()*1.05,
        OneCopy[i]->GetXaxis()->GetBinLowEdge(OneCopy[i]->GetXaxis()->GetFirst()), gPad->GetUymax(),
        gPad->GetUymin()*1.05, gPad->GetUymax(),
        510, "");
    axis->SetLabelFont(43);
    axis->SetLabelSize(25);
    axis->Draw();

    canv->cd();
    p2->SetLogx(true);
    p2->Draw();
    p2->cd();
    TH1D *Ratio1 = MakeRatioPlot(OneCopy[i], TwoCopy[i]);
    Ratio1->SetMarkerColor(Ratio1->GetLineColor());
    Ratio1->Draw("p");
    TH1D *Ratio2 = MakeRatioPlot(OneCopy[i], ThreeCopy[i]);
    Ratio2->SetMarkerColor(Ratio2->GetLineColor());
    Ratio2->Draw("p same");
    TH1D *Ratio3 = MakeRatioPlot(OneCopy[i], FourCopy[i]);
    Ratio3->SetMarkerColor(Ratio3->GetLineColor());
    Ratio3->Draw("p same");

    TLine *line = new TLine(OneCopy[i]->GetXaxis()->GetXmin(), 0.0, OneCopy[i]->GetXaxis()->GetXmax(), 0.0);
    TLine *line2 = new TLine(OneCopy[i]->GetXaxis()->GetXmin(), -1.0, OneCopy[i]->GetXaxis()->GetXmax(), -1.0);
    TLine *line3 = new TLine(OneCopy[i]->GetXaxis()->GetXmin(), 1.0, OneCopy[i]->GetXaxis()->GetXmax(), 1.0);

    line->SetLineColor(kRed);
    line->SetLineStyle(kDashed);
    line->SetLineWidth(2);
    line2->SetLineColor(kRed);
    line2->SetLineStyle(kDashed);
    line2->SetLineWidth(2);
    line3->SetLineColor(kRed);
    line3->SetLineStyle(kDashed);
    line3->SetLineWidth(2);

    line->Draw("same");
    line2->Draw("same");
    line3->Draw("same");
    
    std::stringstream ss;
    ss << i;
    std::string fluxnum = ss.str();

    canv->Print(("newolddatafitsflux_"+fluxnum+".pdf").c_str());

    nFluxParams += OneCopy[i]->GetXaxis()->GetNbins();

    delete axis;
    delete line;
    delete line2;
    delete line3;
    delete Ratio1;
    delete Ratio2;
    delete Ratio3;
  }

  canv->cd();
  canv->SetLogx(false);
  canv->Clear();
  canv->SetBottomMargin(canv->GetBottomMargin()*1.7);


  ///////////////////////////////////////////////////////

  /*Hard code Eb to be relative to 'neut nominal', which is 25(27) MeV for targetC(O). It's a bit messy because we already shift by one as nominal is 0 in the fit. But basic idea is for the priors (Histo named One) for the nu parameters we undo that shift in central value, scale the central value and error by 1/25(27), then reapply the +1 shift.
For the nubar parameters we don't need to undo/redo the shift to the prior as it's already just 1 anyway. So we just scale the errors by 25(27).

Then for the postfits (Histo named Two) they haven't been shifted yet. So we do the scaling then shift the central value. 
  */

  //
  // C nu
  One->SetBinContent(nFluxParams+19,1+(One->GetBinContent(nFluxParams+19)-1)/25);
  One->SetBinError(nFluxParams+19,One->GetBinError(nFluxParams+19)/25);
  Two->SetBinContent(nFluxParams+19,1+Two->GetBinContent(nFluxParams+19)/25);
  Two->SetBinError(nFluxParams+19,Two->GetBinError(nFluxParams+19)/25);
  Three->SetBinContent(nFluxParams+19,1+Three->GetBinContent(nFluxParams+19)/25);
  Three->SetBinError(nFluxParams+19,Three->GetBinError(nFluxParams+19)/25);
  Four->SetBinContent(nFluxParams+19,1+Four->GetBinContent(nFluxParams+19)/25);
  Four->SetBinError(nFluxParams+19,Four->GetBinError(nFluxParams+19)/25);
  // C nubar
  One->SetBinError(nFluxParams+20,One->GetBinError(nFluxParams+20)/25);
  Two->SetBinContent(nFluxParams+20,1+Two->GetBinContent(nFluxParams+20)/25);
  Two->SetBinError(nFluxParams+20,Two->GetBinError(nFluxParams+20)/25);
  Three->SetBinContent(nFluxParams+20,1+Three->GetBinContent(nFluxParams+20)/25);
  Three->SetBinError(nFluxParams+20,Three->GetBinError(nFluxParams+20)/25);
  Four->SetBinContent(nFluxParams+20,1+Four->GetBinContent(nFluxParams+20)/25);
  Four->SetBinError(nFluxParams+20,Four->GetBinError(nFluxParams+20)/25);
  // O nu
  One->SetBinContent(nFluxParams+21,1+(One->GetBinContent(nFluxParams+21)-1)/27);
  One->SetBinError(nFluxParams+21,One->GetBinError(nFluxParams+21)/27);
  Two->SetBinContent(nFluxParams+21,1+Two->GetBinContent(nFluxParams+21)/27);
  Two->SetBinError(nFluxParams+21,Two->GetBinError(nFluxParams+21)/27);
  Three->SetBinContent(nFluxParams+21,1+Three->GetBinContent(nFluxParams+21)/27);
  Three->SetBinError(nFluxParams+21,Three->GetBinError(nFluxParams+21)/27);
  Four->SetBinContent(nFluxParams+21,1+Four->GetBinContent(nFluxParams+21)/27);
  Four->SetBinError(nFluxParams+21,Four->GetBinError(nFluxParams+21)/27);
  // O nubar
  One->SetBinError(nFluxParams+22,One->GetBinError(nFluxParams+22)/27);
  Two->SetBinContent(nFluxParams+22,1+Two->GetBinContent(nFluxParams+22)/27);
  Two->SetBinError(nFluxParams+22,Two->GetBinError(nFluxParams+22)/27);
  Three->SetBinContent(nFluxParams+22,1+Three->GetBinContent(nFluxParams+22)/27);
  Three->SetBinError(nFluxParams+22,Three->GetBinError(nFluxParams+22)/27);
  Four->SetBinContent(nFluxParams+22,1+Four->GetBinContent(nFluxParams+22)/27);
  Four->SetBinError(nFluxParams+22,Four->GetBinError(nFluxParams+22)/27);

  //and hard code fixed 2p2hEdep + I1/2 Low Pi Mom
  
  for (int i=7; i<11; i++) { 
    Two->SetBinContent(nFluxParams+i,-999);
    Two->SetBinError(nFluxParams+i,0);
    Three->SetBinContent(nFluxParams+i,-999);
    Three->SetBinError(nFluxParams+i,0);
    Four->SetBinContent(nFluxParams+i,-999);
    Four->SetBinError(nFluxParams+i,0);
  }

  Two->SetBinContent(nFluxParams+25,-999);
  Two->SetBinError(nFluxParams+25,0);
  Three->SetBinContent(nFluxParams+25,-999);
  Three->SetBinError(nFluxParams+25,0);
  Four->SetBinContent(nFluxParams+25,-999);
  Four->SetBinError(nFluxParams+25,0);

  //Flat prior for low Q2s
  for (int i=11; i<16; i++) {
    One->SetBinError(nFluxParams+i,2.0);
  }


  /////////////////////////////////////////////////// Done hard coding

  // Do some fancy replacements
  PrettifyTitles(One);
  const int XsecPlots = 6;
  int XsecOffset[XsecPlots] = {nFluxParams, nFluxParams+10, nFluxParams+10+12,nFluxParams+10+12+16, nFluxParams+10+12+16+4, nBins};
  std::string titles[XsecPlots-1] = {"CC0#pi", "Q^{2} and E_{b}", "CC1#pi, #nu_{e}, CC DIS, CC Multi #pi, CC coh", "NC", "Pion FSI"};
  One->GetYaxis()->SetTitle("Variation rel. nom.");
  One->GetYaxis()->SetTitleOffset(One->GetYaxis()->GetTitleOffset()*1.2);

  TPad *p3 = new TPad("p3", "p3", 0.0, 0.4, 1.0, 1.0);
  TPad *p4 = new TPad("p4", "p4", 0.0, 0.0, 1.0, 0.4);
  p3->SetLeftMargin(canv->GetLeftMargin());
  p3->SetRightMargin(canv->GetRightMargin());
  p3->SetTopMargin(canv->GetTopMargin());
  p3->SetBottomMargin(0);
  p4->SetLeftMargin(canv->GetLeftMargin());
  p4->SetRightMargin(canv->GetRightMargin());
  p4->SetTopMargin(0);
  p4->SetBottomMargin(0.50);

  //This was from a time I set the nominal xsec values of a fit to be the best fit from the last analysis (this was for det cov binning anyway studies, which we didn't want to tune on data). Anway what I wanted was the prior in these plots to be hardcoded to show any deviation in the fit so you can prolly just ignore this. Cheers thanks, Will x
  //    double bestfit[34] = {0.94, 1.07, 1.0, 1.65, 0.8, 0.95, 1.85, 1.67, 1.25, 1.55, 0.9, 1.0, 1.0, 1.0, 0.9, 1.0, 1.0, 1.0, 1.0, 1.0, 1.3, 0.85, 0.85, 0.9, 1.2, 1.3, 1.0, 0.72, 1.0, 1.0, 0.7, 1.1};

  //for(int i=0; i<32; i++){
  // One->SetBinContent(nFluxParams+i+1,bestfit[i]);
  //}

  for (int i = 1; i < XsecPlots; ++i) {
    
    canv->cd();
    p3->Draw();
    p3->cd();
    One->SetTitle(titles[i-1].c_str());
    One->GetXaxis()->SetRangeUser(XsecOffset[i-1], XsecOffset[i]);
    One->GetYaxis()->SetRangeUser(0, 2.0);
    if(i==3)
      One->GetYaxis()->SetRangeUser(0, 2.3);
    One->GetYaxis()->SetLabelSize(0.);
    One->GetYaxis()->SetTitleSize(0.05);
    One->GetYaxis()->SetTitleOffset(1.3);
    One->Draw("e2");
    if (Two != NULL) {
      Two->GetXaxis()->SetRangeUser(XsecOffset[i-1], XsecOffset[i]);
      Two->Draw("e2,same");
    }
    if (Three != NULL) {
      Three->GetXaxis()->SetRangeUser(XsecOffset[i-1], XsecOffset[i]);
      Three->Draw("e1,same");
    }
    if (Four != NULL) {
      Four->GetXaxis()->SetRangeUser(XsecOffset[i-1], XsecOffset[i]);
      Four->Draw("e1,same");
    }

    canv->Update();
    TGaxis *axis = new TGaxis(One->GetXaxis()->GetBinLowEdge(One->GetXaxis()->GetFirst()), gPad->GetUymin()+0.01,
        One->GetXaxis()->GetBinLowEdge(One->GetXaxis()->GetFirst()), gPad->GetUymax(),
        gPad->GetUymin()+0.01, gPad->GetUymax(),
        510, "");
    axis->SetLabelFont(43);
    axis->SetLabelSize(25);
    axis->Draw();

    canv->cd();
    p4->Draw();
    p4->cd();
    TH1D *RatioOne = MakeRatioPlot(One, Two);
    TH1D *RatioTwo = MakeRatioPlot(One, Three);
    TH1D *RatioThree = MakeRatioPlot(One, Four);

    RatioOne->SetMarkerColor(RatioOne->GetLineColor());
    RatioOne->Draw("p");
    RatioTwo->SetMarkerColor(RatioTwo->GetLineColor());
    RatioTwo->Draw("p same");
    RatioThree->SetMarkerColor(RatioThree->GetLineColor());
    RatioThree->Draw("p same");

    TLine *line = new TLine(RatioOne->GetXaxis()->GetBinLowEdge(RatioOne->GetXaxis()->GetFirst()), 0.0, RatioOne->GetXaxis()->GetBinLowEdge(RatioOne->GetXaxis()->GetLast()+1), 0.0);
    TLine *line2 = new TLine(RatioOne->GetXaxis()->GetBinLowEdge(RatioOne->GetXaxis()->GetFirst()), 1.0, RatioOne->GetXaxis()->GetBinLowEdge(RatioOne->GetXaxis()->GetLast()+1), 1.0);
    TLine *line3 = new TLine(RatioOne->GetXaxis()->GetBinLowEdge(RatioOne->GetXaxis()->GetFirst()), -1.0, RatioOne->GetXaxis()->GetBinLowEdge(RatioOne->GetXaxis()->GetLast()+1), -1.0);

    line->SetLineColor(kRed);
    line->SetLineStyle(kDashed);
    line->SetLineWidth(2);
    line2->SetLineColor(kRed);
    line2->SetLineStyle(kDashed);
    line2->SetLineWidth(2);
    line3->SetLineColor(kRed);
    line3->SetLineStyle(kDashed);
    line3->SetLineWidth(2);

    line->Draw("same");
    line2->Draw("same");
    line3->Draw("same");

    std::stringstream ss;
    ss << i;
    std::string xsecnum = ss.str();

    canv->Print(("newolddatafitsxsec_"+xsecnum+".pdf").c_str());
    
    delete axis;
    delete line;
    delete line2;
    delete line3;
   
    delete RatioOne;
    delete RatioTwo;
    delete RatioThree;
   }

  //  canv->Print((FileName1+"_"+FileName2+"_"+FileName3+".pdf]").c_str());
}
