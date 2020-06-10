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
#include "TH2Poly.h"

void prioronly1D_v2(std::string pfile) {

  //  gStyle->SetPalette(51);

  TCanvas *c = new TCanvas("canv", "canv", 1580, 1080);
  TFile *filep = TFile::Open(pfile.c_str());//d
 
  c->SetTopMargin(0.1);
  c->SetBottomMargin(0.16);
  c->SetLeftMargin(0.12);
  c->SetRightMargin(0.1);

  TPad *p1 = new TPad("p1", "p1", 0.0, 0.3, 1.0, 1.0);
  TPad *p2 = new TPad("p2", "p2", 0.0, 0.0, 1.0, 0.3);
  p1->SetLeftMargin(c->GetLeftMargin());
  p1->SetRightMargin(c->GetRightMargin());
  p1->SetTopMargin(c->GetTopMargin());
  p1->SetBottomMargin(0);
  p2->SetLeftMargin(c->GetLeftMargin());
  p2->SetRightMargin(c->GetRightMargin());
  p2->SetTopMargin(0);
  p2->SetBottomMargin(0.35);
  
  std::string samp[18] = {"FGD1 numuCC 0pi", "FGD1 numuCC 1pi", "FGD1 numuCC other", "FGD1 anti-numuCC 0pi", "FGD1 anti-numuCC 1pi", "FGD1 anti-numuCC other", "FGD1 NuMuBkg CC0pi in AntiNu Mode", "FGD1 NuMuBkg CC1pi in AntiNu Mode", "FGD1 NuMuBkg CCother in AntiNu Mode", "FGD2 numuCC 0pi", "FGD2 numuCC 1pi", "FGD2 numuCC other", "FGD2 anti-numuCC 0pi","FGD2 anti-numuCC 1pi", "FGD2 anti-numuCC other", "FGD2 NuMuBkg CC0pi in AntiNu Mode", "FGD2 NuMuBkg CC1pi in AntiNu Mode", "FGD2 NuMuBkg CCother in AntiNu Mode"};  

  std::string samp_us[18] = {"FGD1_numuCC_0pi", "FGD1_numuCC_1pi", "FGD1_numuCC_other", "FGD1_anti-numuCC_0pi", "FGD1_anti-numuCC_1pi", "FGD1_anti-numuCC_other", "FGD1_NuMuBkg_CC0pi_in_AntiNu_Mode", "FGD1_NuMuBkg_CC1pi_in_AntiNu_Mode", "FGD1_NuMuBkg_CCother_in_AntiNu_Mode", "FGD2_numuCC_0pi", "FGD2_numuCC_1pi", "FGD2_numuCC_other", "FGD2_anti-numuCC_0pi","FGD2_anti-numuCC_1pi", "FGD2_anti-numuCC_other", "FGD2_NuMuBkg_CC0pi_in_AntiNu_Mode", "FGD2_NuMuBkg_CC1pi_in_AntiNu_Mode", "FGD2_NuMuBkg_CCother_in_AntiNu_Mode"};

  TLegend *leg = new TLegend(0.0, 0.0, 1.0, 1.0);

  for(int s=0; s<18; s++){

    std::string name = samp[s] + "/" + samp_us[s] + "_nom_mean_x_x";
    std::cout << name << std::endl;
    TH1D* priorX = (TH1D*)filep->Get(name.c_str())->Clone();
    //    name = samp_us[s] + "_post_x";
    //TH1D* postX = (TH1D*)filep->Get(name.c_str())->Clone();
    name = samp[s] + "/" + samp_us[s] + "_nom_mean_y_y";
    TH1D* priorY = (TH1D*)filep->Get(name.c_str())->Clone();
    //name = samp_us[s] + "_post_y";
    //TH1D* postY = (TH1D*)filep->Get(name.c_str())->Clone();

    std::string datname = samp[s] + "/" + samp_us[s] + "_data_x_x";
    TH1D* dathistX = (TH1D*)filep->Get(datname.c_str())->Clone();
    datname = samp[s] + "/" +samp_us[s] + "_data_y_y";
    TH1D* dathistY = (TH1D*)filep->Get(datname.c_str())->Clone();

    priorX->Sumw2();
    //postX->Sumw2();
    dathistX->Sumw2();
    priorY->Sumw2();
    //postY->Sumw2();
    dathistY->Sumw2();

    priorX->Scale(1, "width");
    // postX->Scale(1, "width");
    dathistX->Scale(1, "width");
    priorY->Scale(1, "width");
    //postY->Scale(1, "width");
    dathistY->Scale(1, "width");
    
    dathistX->SetLineColor(kBlack);
    //dathistX->SetLineWidth(2);
    priorX->SetLineColor(kRed);
    //priorX->SetLineWidth(2);
    //postX->SetLineColor(kBlue);
    //postX->SetLineWidth(2);
    dathistY->SetLineColor(kBlack);
    //dathistY->SetLineWidth(2);
    priorY->SetLineColor(kRed);
    //priorY->SetLineWidth(2);
    //postY->SetLineColor(kBlue);
    //    postY->SetLineWidth(2);    

    if(s==0){
      leg->AddEntry(dathistX, "Data", "lp");
      leg->AddEntry(priorX, "Prior Prediction", "lp");
      //leg->AddEntry(postX, "Posterior Prediction", "lp");
      leg->Draw();
      c->Print((std::string("prior1dleg.pdf")).c_str());
      //      for (int b=1; b<=postX->GetNbinsX(); b++){
      //std::cout << b << " " << postX->GetBinError(b) << std::endl;
      //      }
    }

    TH1D* priorratioX = (TH1D*)priorX->Clone();
    //TH1D* postratioX = (TH1D*)postX->Clone();
    TH1D* priorratioY = (TH1D*)priorY->Clone();
    // TH1D* postratioY = (TH1D*)postY->Clone();
    TH1D* dataratioX = (TH1D*)dathistX->Clone();
    TH1D* dataratioY = (TH1D*)dathistY->Clone();

    priorratioX->Sumw2();
    // postratioX->Sumw2();
    priorratioX->Divide(dathistX);
    //postratioX->Divide(dathistX);
    priorratioY->Sumw2();
    //postratioY->Sumw2();
    priorratioY->Divide(dathistY);
    //postratioY->Divide(dathistY);
    dataratioY->Sumw2();
    dataratioX->Sumw2();
    dataratioY->Divide(dathistY);
    dataratioX->Divide(dathistX);

    int maxX;
    if(s<2 || (s>8 && s<12))
      maxX=5000;
    else if(s==3||s==12||s==5||s==14||s==6||s==15||s==8||s==17)
      maxX=4000;
    else if(s==4||s==13)
      maxX=2500;
    else if(s==7||s==16)
      maxX=1500;

    gPad->Clear();

    c->cd();
    p1->Draw();
    p1->cd();
    dathistX->SetTitle("");
    dathistX->Draw("e");
    priorX->Draw("e same");
    //    postX->Draw("e same");
    dathistX->SetTitle("");
    dathistX->GetXaxis()->SetTitle("p_{#mu} (MeV)");
    dathistX->GetYaxis()->SetTitle("Events/MeV");
    dathistX->GetYaxis()->SetTitleOffset(0.77);
    dathistX->GetXaxis()->SetRangeUser(0,maxX);
    dathistX->GetYaxis()->SetRangeUser(0.1,1.25*dathistX->GetMaximum());
    dathistX->GetYaxis()->SetLabelSize(0.05);
    dathistX->GetYaxis()->SetTitleSize(0.07);
    c->Update();

    c->cd();
    p2->Draw();
    p2->cd();
    priorratioX->SetTitle("");
    priorratioX->GetXaxis()->SetTitle("p_{#mu} (MeV)");
    priorratioX->GetYaxis()->SetTitle("MC/Data");
    priorratioX->GetYaxis()->SetTitleOffset(0.36);
    priorratioX->GetXaxis()->SetRangeUser(0,maxX);
    priorratioX->GetYaxis()->SetRangeUser(0.51,1.49);
    priorratioX->GetYaxis()->SetLabelSize(0.12);
    priorratioX->GetYaxis()->SetNdivisions(305);
    priorratioX->GetYaxis()->SetTitleSize(0.15);
    priorratioX->GetXaxis()->SetLabelSize(0.12);
    priorratioX->GetXaxis()->SetTitleSize(0.15);
    priorratioX->Draw("e");
    //postratioX->Draw("e same");
    dataratioX->SetFillStyle(3144);
    dataratioX->SetFillColor(kBlack);
    dataratioX->SetFillStyle(3003);
    dataratioX->SetMarkerStyle(7);
    dataratioX->SetMarkerColor(kBlack);
    dataratioX->SetLineColor(kBlack);
    dataratioX->Draw("e2 same");
    TLine *line = new TLine(0, 1.0, maxX, 1.0);
    line->SetLineColor(kBlack);
    // line->Draw("same");
    
    c->Print((std::string("prioronly1D_p_")+samp_us[s]+std::string(".pdf")).c_str());

    p1->cd();
    dathistY->SetTitle("");
    dathistY->Draw("e");
    priorY->Draw("esame");
    //postY->Draw("esame");
    dathistY->SetTitle("");
    dathistY->GetXaxis()->SetTitle("cos #theta_{#mu}");
    dathistY->GetYaxis()->SetTitle("Events/MeV");
    dathistY->GetYaxis()->SetTitleOffset(0.77);
    dathistY->GetXaxis()->SetRangeUser(0.7,1.0);
    dathistY->GetYaxis()->SetRangeUser(0.1,1.25*dathistY->GetMaximum());
    dathistY->GetYaxis()->SetLabelSize(0.05);
    dathistY->GetYaxis()->SetTitleSize(0.07);
    c->Update();

    c->cd();
    p2->Draw();
    p2->cd();
    priorratioY->SetTitle("");
    priorratioY->GetXaxis()->SetTitle("cos #theta_{#mu}");
    priorratioY->GetYaxis()->SetTitle("MC/Data");
    priorratioY->GetYaxis()->SetTitleOffset(0.36);
    priorratioY->GetXaxis()->SetRangeUser(0.7,1.0);
    priorratioY->GetYaxis()->SetRangeUser(0.51,1.49);
    priorratioY->GetYaxis()->SetLabelSize(0.12);
    priorratioY->GetYaxis()->SetNdivisions(305);
    priorratioY->GetYaxis()->SetTitleSize(0.15);
    priorratioY->GetXaxis()->SetLabelSize(0.12);
    priorratioY->GetXaxis()->SetTitleSize(0.15);
    priorratioY->Draw("e");
    //postratioY->Draw("e same");
    dataratioY->Draw("e2 same");
    dataratioY->SetFillStyle(3144);
    dataratioY->SetFillColor(kBlack);
    dataratioY->SetFillStyle(3003);
    dataratioY->SetMarkerStyle(7);
    dataratioY->SetMarkerColor(kBlack);
    dataratioY->SetLineColor(kBlack);
    TLine *line2 = new TLine(0.7, 1.0, 1.0, 1.0);
    line2->SetLineColor(kBlack);
    //    line2->Draw("same");

    c->Print((std::string("prioronly1D_t_")+samp_us[s]+std::string(".pdf")).c_str());
  }


}
