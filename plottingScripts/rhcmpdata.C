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
    else if (title == "pF_C") title = "p_{F}^{C}";
    else if (title == "pF_O") title = "p_{F}^{O}";

    else if (title == "2p2h_norm_nu") title = "2p2h norm #nu";
    else if (title == "2p2h_norm_nubar") title = "2p2h norm #bar{#nu}";
    else if (title == "2p2h_normCtoO") title = "2p2h norm ^{12}C/^{16}O";

    else if (title == "2p2h_shape_C") title = "2p2h shape ^{12}C";
    else if (title == "2p2h_shape_O") title = "2p2h shape ^{16}O";

    else if (title == "BeRPA_A") title = "BeRPA A";
    else if (title == "BeRPA_B") title = "BeRPA B";
    else if (title == "BeRPA_D") title = "BeRPA D";
    else if (title == "BeRPA_E") title = "BeRPA E";
    else if (title == "BeRPA_U") title = "BeRPA U";

    else if (title == "CA5") title = "C_{5}^{A}";
    else if (title == "MARES") title = "M_{A}^{RES}";
    else if (title == "ISO_BKG") title = "Non-res I_{1/2}";

    else if (title == "CC_norm_nu") title = "CC Norm #nu";
    else if (title == "CC_norm_nubar") title = "CC Norm #bar{nu}";
    else if (title == "nue_numu") title = "#nu_{e}/#nu_{#mu}";
    else if (title == "nuebar_numubar") title = "#bar{#nu}_{e}/#bar{#nu}_{#mu}";

    else if (title == "CC_DIS") title = "CC DIS";
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

    else if (title == "EB_dial_C") title = "EB Dial C";
    else if (title == "EB_dial_O") title = "EB Dial O";


    Hist->GetXaxis()->SetBinLabel(i+1, title.c_str());
  }

}

TH1D* MakeRatioPlot(TH1D* OneCopy, TH1D *TwoCopy) {

  TH1D *Ratio = (TH1D*)OneCopy->Clone();
  Ratio->GetYaxis()->SetTitle("(x_{fit}-#mu_{Prior})/#sigma_{Prior}");
  Ratio->SetMinimum(-2.3);
  Ratio->SetMaximum(2.3);
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

void rhcmpdata(std::string FileName1) {

  TFile *File1 = new TFile(FileName1.c_str());

  //gStyle->SetPalette(51);
  TCanvas *canv = new TCanvas("canv", "canv", 1024, 1024);
  canv->SetLeftMargin(0.12);
  canv->SetBottomMargin(0.12);
  canv->SetTopMargin(0.08);
  canv->SetRightMargin(0.04);
  //canv->Print((FileName1+FileName2+".pdf[").c_str());
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
  TH1D *Three = (TH1D*)File1->Get("param_flux_HPD_1");

  /*  if (Two != NULL)
    {
      Two->SetFillStyle(3001);
      Two->SetFillColor(kBlue-1);
      Two->SetMarkerColor(Two->GetFillColor());
      Two->SetMarkerStyle(20);
      Two->SetLineColor(Two->GetFillColor());
      Two->SetMarkerSize(Two->GetMarkerSize());
    }

  if (Three != NULL)
    {
      Three->SetFillStyle(3244);
      Three->SetFillColor(kOrange-5);
      Three->SetMarkerColor(Three->GetFillColor());
      Three->SetMarkerStyle(23);
      Three->SetLineWidth(2);
      Three->SetMarkerSize((One->GetMarkerSize())*0.75);
      Three->SetLineColor(Three->GetFillColor());
      }*/

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

  // Make a Legend page
  TLegend *leg = new TLegend(0.0, 0.0, 1.0, 1.0);
  if (One != NULL) leg->AddEntry(One, "Prior", "lpf");
  if (Two != NULL) leg->AddEntry(Two, "RHC multi-track", "lpf");
  if (Three != NULL) leg->AddEntry(Three, "RHC multi-#pi", "lp");
  canv->cd();
  canv->Clear();
  leg->Draw();
  canv->Print("rhcmpdat248_leg.pdf");
  //  canv->Print((FileName1+FileName2+".pdf").c_str());
  delete leg;

  int nBins = One->GetXaxis()->GetNbins();

  // Make eight flux plots
  const int nFluxPlots = 16;
  TH1D* OneCopy[nFluxPlots];
  TH1D* TwoCopy[nFluxPlots];
  TH1D* ThreeCopy[nFluxPlots];

  int nFluxParams = 0;
  for (int i = 0; i < nFluxPlots; ++i) {

    OneCopy[i]  = PrettifyMe(One, i);
    if (Two != NULL) TwoCopy[i]  = PrettifyMe(Two, i);
    if (Three != NULL) ThreeCopy[i] = PrettifyMe(Three, i);

    canv->cd();
    p1->Draw();
    p1->cd();
    p1->SetLogx(true);
    OneCopy[i]->Draw("e2");
    if (Two != NULL) TwoCopy[i]->Draw("e2,same");
    if (Three != NULL) ThreeCopy[i]->Draw("e1,same");

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

    TLine *line = new TLine(OneCopy[i]->GetXaxis()->GetXmin(), 0.0, OneCopy[i]->GetXaxis()->GetXmax(), 0.0);
    TLine *line2 = new TLine(OneCopy[i]->GetXaxis()->GetXmin(), -1.0, OneCopy[i]->GetXaxis()->GetXmax(), -1.0);
    line->SetLineColor(kRed);
    line->SetLineStyle(kDashed);
    line->SetLineWidth(2);
    line2->SetLineColor(kRed);
    line2->SetLineStyle(kDashed);
    line2->SetLineWidth(2);

    line->Draw("same");
    line2->Draw("same");
    Ratio2->Draw("p,same");

    std::stringstream ss;
    ss << i;
    std::string fluxnum = ss.str();

    //canv->Print((FileName1+FileName2+".pdf").c_str());
    canv->Print(("rhcmpdat248flux_"+fluxnum+".pdf").c_str());
    nFluxParams += OneCopy[i]->GetXaxis()->GetNbins();

    delete axis;
    delete line;
    delete line2;
    delete Ratio1;
    delete Ratio2;
  }

  canv->cd();
  canv->SetLogx(false);
  canv->Clear();
  canv->SetBottomMargin(canv->GetBottomMargin()*1.7);

  // Do some fancy replacements
  PrettifyTitles(One);
  const int XsecPlots = 5;
  int XsecOffset[XsecPlots] = {nFluxParams, nFluxParams+13, nFluxParams+13+10, nFluxParams+13+10+4, nBins};
  std::string titles[XsecPlots-1] = {"CC0#pi", "CC1#pi, #nu_{e}, CC DIS, CC coh", "NC", "Pion FSI + EB"};
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

  for (int i = 1; i < XsecPlots; ++i) {
    canv->cd();
    p3->Draw();
    p3->cd();
    One->SetTitle(titles[i-1].c_str());
    One->GetXaxis()->SetRangeUser(XsecOffset[i-1], XsecOffset[i]);
    One->GetYaxis()->SetRangeUser(0, 2.0);
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

    RatioOne->SetMarkerColor(RatioOne->GetLineColor());
    RatioOne->Draw("p");

    TH1D *RatioTwo = MakeRatioPlot(One, Three);
    RatioTwo->SetMarkerColor(RatioTwo->GetLineColor());

    TLine *line = new TLine(RatioOne->GetXaxis()->GetBinLowEdge(RatioOne->GetXaxis()->GetFirst()), 0.0, RatioOne->GetXaxis()->GetBinLowEdge(RatioOne->GetXaxis()->GetLast()+1), 0.0);
    TLine *line2 = new TLine(RatioOne->GetXaxis()->GetBinLowEdge(RatioOne->GetXaxis()->GetFirst()), 1.0, RatioOne->GetXaxis()->GetBinLowEdge(RatioOne->GetXaxis()->GetLast()+1), 1.0);
    line->SetLineColor(kRed);
    line->SetLineStyle(kDashed);
    line->SetLineWidth(2);
    line2->SetLineColor(kRed);
    line2->SetLineStyle(kDashed);
    line2->SetLineWidth(2);
    line->Draw("same");
    line2->Draw("same");

    RatioTwo->Draw("p,same");

    //    canv->Print((FileName1+FileName2+".pdf").c_str());
    
    std::stringstream ss;
    ss << i;
    std::string xsecnum = ss.str();

    //canv->Print((FileName1+FileName2+".pdf").c_str());
    canv->Print(("rhcmpdatxsec248_"+xsecnum+".pdf").c_str());

    delete axis;
    delete line;
    delete line2;

    delete RatioOne;
    delete RatioTwo;
  }

  //  canv->Print((FileName1+FileName2+".pdf]").c_str());
}
