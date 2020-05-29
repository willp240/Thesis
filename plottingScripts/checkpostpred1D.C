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

void checkpostpred1D(std::string file0) {

  //  gStyle->SetPalette(51);

  TCanvas *c = new TCanvas("canv", "canv", 1580, 1080);
  TFile *file_ = TFile::Open(file0.c_str());
  TFile *outfile = new TFile("outtest1D.root","RECREATE");
 
  c->SetTopMargin(0.1);
  c->SetBottomMargin(0.16);
  c->SetLeftMargin(0.12);
  c->SetRightMargin(0.1);
  TPad *p1 = new TPad("p1", "p1", 0.0, 0.3, 1.0, 1.0);
  TPad *p2 = new TPad("p2", "p2", 0.0, 0.0, 1.0, 0.3);
  p1->SetLeftMargin(canv->GetLeftMargin());
  p1->SetRightMargin(canv->GetRightMargin());
  p1->SetTopMargin(canv->GetTopMargin());
  p1->SetBottomMargin(0);
  p2->SetLeftMargin(canv->GetLeftMargin());
  p2->SetRightMargin(canv->GetRightMargin());
  p2->SetTopMargin(0);
  p2->SetBottomMargin(0.35);
  
  std::string samp[18] = {"FGD1 numuCC 0pi", "FGD1 numuCC 1pi", "FGD1 numuCC other", "FGD1 anti-numuCC 0pi", "FGD1 anti-numuCC 1pi", "FGD1 anti-numuCC other", "FGD1 NuMuBkg CC0pi in AntiNu Mode", "FGD1 NuMuBkg CC1pi in AntiNu Mode", "FGD1 NuMuBkg CCother in AntiNu Mode", "FGD2 numuCC 0pi", "FGD2 numuCC 1pi", "FGD2 numuCC other", "FGD2 anti-numuCC 0pi","FGD2 anti-numuCC 1pi", "FGD2 anti-numuCC other", "FGD2 NuMuBkg CC0pi in AntiNu Mode", "FGD2 NuMuBkg CC1pi in AntiNu Mode", "FGD2 NuMuBkg CCother in AntiNu Mode"};  

  std::string samp_us[18] = {"FGD1_numuCC_0pi", "FGD1_numuCC_1pi", "FGD1_numuCC_other", "FGD1_anti-numuCC_0pi", "FGD1_anti-numuCC_1pi", "FGD1_anti-numuCC_other", "FGD1_NuMuBkg_CC0pi_in_AntiNu_Mode", "FGD1_NuMuBkg_CC1pi_in_AntiNu_Mode", "FGD1_NuMuBkg_CCother_in_AntiNu_Mode", "FGD2_numuCC_0pi", "FGD2_numuCC_1pi", "FGD2_numuCC_other", "FGD2_anti-numuCC_0pi","FGD2_anti-numuCC_1pi", "FGD2_anti-numuCC_other", "FGD2_NuMuBkg_CC0pi_in_AntiNu_Mode", "FGD2_NuMuBkg_CC1pi_in_AntiNu_Mode", "FGD2_NuMuBkg_CCother_in_AntiNu_Mode"};

  for(int s=0; s<1; s++){
    
    std::stringstream ss;
    ss << samp_us[s] << "_" << 0;
    std::string name = ss.str();
    std::cout << "Getting " << name << std::endl;
    
    TH1D* projX = (TH1D*)file_->Get(name.c_str())->ProjectionX()->Clone();
    TH1D* projY = (TH1D*)file_->Get(name.c_str())->ProjectionY()->Clone();

    double rms_sum=0;
    
    for (int b=1; b<projX->GetNumberOfBins(); b++){
      TString name2 = samp_us[s] + "_bin_";
      name2+=b;
      TH1D* events = new TH1D("h1", "h1", 2000,0,300);
      //      events->SetName(name2);
      // events->SetTitle(name2);

      for(int t=0; t<2000; t++){
	std::stringstream ss3;
	ss3 << samp_us[s] << "_" << t;
	std::string name3 = ss3.str();
	TH1D* projX2 = (TH1D*)file_->Get(name.c_str())->ProjectionX()->Clone();
	events->Fill(projX2->GetBinContent(b));
	delete projX2;
      }
      std::cout << "writing " << name2 << std::endl;
      outfile->cd();
      events->Write(name2);
      rms_sum+=events->GetRMS()*events->GetRMS();
      std::cout << "rms=" << events->GetRMS() << " rms^2=" << events->GetRMS() << " sum=" << rms_sum << std::endl;
      delete events;
    }
    delete projX;
    delete projY;
  }
  
    /*
    std::string name = samp[s] + "/" + samp_us[s] + "_nom_mean_y_y";
    //    name = "FGD1 numuCC 0pi/FGD1_numuCC_0pi_nom_mean_ratio";
    TH1D* hist = (TH1D*)file_->Get(name.c_str())->Clone();
    std::string datname = samp[s] + "/" + samp_us[s] + "_data_y_y";
    TH1D* dathistX = (TH1D*)file_->Get(datname.c_str())->Clone();
    TH1D* datratioX = (TH1D*)file_->Get(datname.c_str())->Clone();
    datratioX->Sumw2();
    hist->Sumw2();
    datratioX->Sumw2();

    //    for(int i=0; i<hist->GetXaxis()->GetNbins(); i++){
    //hist->SetBinContent(i, hist->GetBinContent(i)/hist->GetBinWidth(i));
    //dathistX->SetBinContent(i, dathistX->GetBinContent(i)/dathistX->GetBinWidth(i));
    //datratioX->SetBinContent(i, datratioX->GetBinContent(i)/datratioX->GetBinWidth(i));
    //}

    //    hist->Scale(1, "width");
    //dathistX->Scale(1, "width");
    //datratioX->Scale(1, "width");

    datratioX->Divide(hist);

    dathistX->SetLineColor(kBlack);
    datratioX->SetLineColor(kBlack);
    hist->SetLineColor(kRed);
    hist->SetLineWidth(2);
    dathistX->SetLineWidth(2);
    dathistX->SetLineColor(kBlack);

    int maxX;
    if(s<2 || (s>8 && s<12))
      maxX=5000;
    else if(s==3||s==12||s==5||s==14||s==6||s==15||s==8||s==17)
      maxX=4000;
    else if(s==4||s==13)
      maxX=2500;
    else if(s==7||s==16)
      maxX=1500;

    maxX=1.0;

    gPad->Clear();

    c->cd();
    p1->Draw();
    p1->cd();
    hist->GetXaxis()->SetRangeUser(0.6,1.0);
    dathistX->GetXaxis()->SetRangeUser(0.6,1.0);
    hist->GetXaxis()->SetRangeUser(0.6,1.0);
    dathistX->SetTitle("");
    dathistX->Draw("e");
    hist->Draw("same");
    dathistX->SetTitle("");
    dathistX->GetXaxis()->SetTitle("cos #theta_{#mu}");
    dathistX->GetYaxis()->SetTitle("Events");
    dathistX->GetYaxis()->SetTitleOffset(0.85);
    dathistX->GetXaxis()->SetRangeUser(0.6,1.0);
    hist->GetXaxis()->SetRangeUser(0.6,1.0);
    dathistX->GetYaxis()->SetRangeUser(0.1,1.25*dathistX->GetMaximum());
    dathistX->GetYaxis()->SetLabelSize(0.05);
    dathistX->GetYaxis()->SetTitleSize(0.07);
    c->Update();

    c->cd();
    p2->Draw();
    p2->cd();
    datratioX->SetTitle("");
    datratioX->GetXaxis()->SetTitle("cos #theta_{#mu}");
    datratioX->GetYaxis()->SetTitle("Data/MC ");
    datratioX->GetYaxis()->SetTitleOffset(0.36);
    datratioX->GetXaxis()->SetRangeUser(0.6,1.0);
    datratioX->GetYaxis()->SetRangeUser(0.01,1.99);
    datratioX->GetYaxis()->SetNdivisions(305);
    datratioX->GetYaxis()->SetLabelSize(0.1);
    datratioX->GetYaxis()->SetTitleSize(0.15);
    datratioX->GetXaxis()->SetLabelSize(0.12);
    datratioX->GetXaxis()->SetTitleSize(0.15);
    datratioX->Draw("e");
    TLine *line = new TLine(datratioX->GetXaxis()->GetXmin(), 1.0, maxX, 1.0);
    line->SetLineColor(kRed);
    line->Draw("same");

    c->Print((std::string("postpred1D_")+samp_us[s]+std::string(".pdf")).c_str());
    }*/
}
