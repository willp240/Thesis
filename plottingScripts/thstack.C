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
#include <vector>
#include "TString.h"

void thstack(std::string filename) {

    TCanvas *c = new TCanvas("canv", "canv", 1580, 1080);
    TFile *file_ = TFile::Open(filename.c_str());

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

    
    std::string samp[18] = {"FGD1_numuCC_0pi", "FGD1_numuCC_1pi", "FGD1_numuCC_other", "FGD1_anti-numuCC_0pi", "FGD1_anti-numuCC_1pi", "FGD1_anti-numuCC_other", "FGD1_NuMuBkg_CC0pi_in_AntiNu_Mode", "FGD1_NuMuBkg_CC1pi_in_AntiNu_Mode", "FGD1_NuMuBkg_CCother_in_AntiNu_Mode", "FGD2_numuCC_0pi", "FGD2_numuCC_1pi", "FGD2_numuCC_other", "FGD2_anti-numuCC_0pi","FGD2_anti-numuCC_1pi", "FGD2_anti-numuCC_other", "FGD2_NuMuBkg_CC0pi_in_AntiNu_Mode", "FGD2_NuMuBkg_CC1pi_in_AntiNu_Mode", "FGD2_NuMuBkg_CCother_in_AntiNu_Mode"};  
    
    std::string mode[13] = {"", "CCQE", "2p2h", "CC1pi", "CCcoh", "CCMpi", "CCDIS", "CCMisc", "NC1pi0", "NC1pipm", "NCcoh", "NCoth", "NC1gam"};

    std::string title[13] = {"", "CCQE", "2p2h", "CC 1#pi", "CC coherent", "CC mult-#pi", "CC DIS", "CC miscellaneous", "NC 1#pi^{0}", "NC 1#pi^{#pm}", "NC coherent", "NC other", "NC 1#gamma"}; 

    TLegend *leg = new TLegend(0.0, 0.0, 1.0, 1.0);

    for(int s=0; s<2; s++){
      
      // Stacks for the projected samples
      std::string momname = samp[s]+"_x";
      std::string thname = samp[s]+"_y";
      THStack *MomStack = new THStack((momname).c_str(), (momname).c_str());
      THStack *ThStack = new THStack((thname).c_str(), (thname).c_str());
      std::cout << "start " << std::endl;
      for(int m=0; m<13; m++){
	
	if(m==0){
	  
	  std::string datname = "DATA_" + samp[s];
	  std::cout << "getting " << datname << std::endl;
	  momname = datname+"_x";
	  thname = datname+"_y";
	  TH1D *dathistX = (TH1D*) file_->Get((momname).c_str())->Clone();
	  TH1D *dathistY = (TH1D*) file_->Get((thname).c_str())->Clone();

	  for (int i=1; i<dathistX->GetXaxis()->GetNbins()+1; i++){
	    dathistX->SetBinError(i, dathistX->GetMaximum()/40);
	  }
	  for (int i=1;i<dathistY->GetXaxis()->GetNbins()+1; i++){
            dathistY->SetBinError(i, dathistY->GetMaximum()/40);
          }

	  dathistX->SetMarkerColor(kBlack);
	  dathistY->SetMarkerColor(kBlack);
	  dathistX->SetMarkerStyle(2);
          dathistY->SetMarkerStyle(2);
	  //	  dathistX->SetMarkerSize(2);
	  // dathistY->SetMarkerSize(2);
	  dathistX->SetLineColor(kBlack);
	  dathistY->SetLineColor(kBlack);
	 
	  std::string mcname = "MC_" + samp[s];
	  std::cout << "getting " << mcname << std::endl;
          momname = mcname+"_x";
          thname = mcname+"_y";
          TH1D *mchistX = (TH1D*) file_->Get((momname).c_str())->Clone();
          TH1D *mchistY = (TH1D*) file_->Get((thname).c_str())->Clone();

          mchistX->SetLineColor(kRed);
          mchistY->SetLineColor(kRed);
	  mchistX->SetFillStyle(0);
          mchistY->SetFillStyle(0);
	  if(s==0){
	    leg->AddEntry(dathistX, "Data", "lp");
	    leg->AddEntry(mchistX, "MC", "l");
	  }

	  datname = "DATA_" + samp[s];
          momname = datname+"_x";
          thname = datname+"_y";
	  TH1D *datratioX = (TH1D*) file_->Get((momname).c_str())->Clone();
	  datratioX->Sumw2();
	  datratioX->Divide(mchistX);
	  TH1D *datratioY = (TH1D*) file_->Get((thname).c_str())->Clone();
	  datratioY->Sumw2();
	  datratioY->Divide(mchistY);
	  datratioX->SetMarkerColor(kBlack);
	  datratioY->SetMarkerColor(kBlack);
	  datratioX->SetLineColor(kBlack);
          datratioY->SetLineColor(kBlack);
	  datratioX->SetMarkerStyle(2);	  
	  datratioY->SetMarkerStyle(2);
	  
	  mcname = "MC_" + samp[s];
	  momname = mcname+"_x";
          thname = mcname+"_y";
	  TH1D *mcratioX = (TH1D*) file_->Get((momname).c_str())->Clone();
	  mcratioX->Divide(mchistX);
          TH1D *mcratioY = (TH1D*) file_->Get((thname).c_str())->Clone();
	  mcratioY->Divide(mchistY);
	  mcratioX->SetLineColor(kRed);
          mcratioY->SetLineColor(kRed);
	  mcratioX->SetFillStyle(0);
          mcratioY->SetFillStyle(0);
	}
	
	else {
	  std::string name = "MC_" + samp[s];
	  name = name + "_" + mode[m];
	  std::cout << "getting " << name << std::endl;
	  momname = name+"_x";
	  thname = name+"_y";
	  TH1D* projX = (TH1D*) file_->Get((momname).c_str())->Clone();
	  TH1D* projY = (TH1D*) file_->Get((thname).c_str())->Clone();
	  
	  projX->SetFillStyle(1001);
	  projY->SetFillStyle(1001);
	  projX->SetLineColor(0);
	  projY->SetLineColor(0);
	  projX->SetLineStyle(0);
	  projY->SetLineStyle(0);
	  projX->SetLineWidth(0);
          projY->SetLineWidth(0);

	  switch(m){
	   
	  case 1:
	    projX->SetFillColor(kGreen+4);
	    projY->SetFillColor(kGreen+4);
	    break;
	  case 2:
            projX->SetFillColor(kGreen+2);
            projY->SetFillColor(kGreen+2);
            break;
	  case 3:
            projX->SetFillColor(kRed+3);
            projY->SetFillColor(kRed+3);
            break;
	  case 4:
            projX->SetFillColor(kRed+1);
            projY->SetFillColor(kRed+1);
            break;
	  case 5:
            projX->SetFillColor(kRed-3);
            projY->SetFillColor(kRed-3);
            break;
	  case 6:
            projX->SetFillColor(kRed-6);
            projY->SetFillColor(kRed-6);
            break;
	  case 7:
            projX->SetFillColor(kRed-10);
            projY->SetFillColor(kRed-10);
            break;
	  case 8:
            projX->SetFillColor(kBlue+3);
            projY->SetFillColor(kBlue+3);
            break;
	  case 9:
            projX->SetFillColor(kBlue+1);
            projY->SetFillColor(kBlue+1);
            break;
	  case 10:
            projX->SetFillColor(kBlue-9);
            projY->SetFillColor(kBlue-9);
            break;
	  case 11:
            projX->SetFillColor(kBlue-10);
            projY->SetFillColor(kBlue-10);
            break;
	  case 12:
            projX->SetFillColor(kCyan);
            projY->SetFillColor(kCyan);
            break;
	  }

	  projX->SetLineColor(projX->GetFillColor());
	  projY->SetLineColor(projY->GetFillColor());

	  if(s==0){
	    leg->AddEntry(projX, title[m].c_str(), "f");
	  }
	  
	  MomStack->Add(projX);
	  ThStack->Add(projY);
	}
      }

      if(s==0){
	leg->Draw();
	c->Print((std::string("legend.pdf")).c_str());
      }

      gPad->Clear();

      c->cd();
      p1->Draw();
      p1->cd();
      dathistX->SetTitle("");
      dathistX->Draw("e");
      MomStack->Draw("same");
      mchistX->Draw("same");
      dathistX->Draw("e same");
      dathistX->SetTitle("");
      dathistX->GetXaxis()->SetTitle("p_{#mu} (MeV)");
      dathistX->GetYaxis()->SetTitle("Events/MeV");
      dathistX->GetYaxis()->SetTitleOffset(0.8);
      dathistX->GetXaxis()->SetRangeUser(0,5000);
      dathistX->GetYaxis()->SetRangeUser(0.1,1.25*dathistX->GetMaximum());
      dathistX->GetYaxis()->SetLabelSize(0.05);
      dathistX->GetYaxis()->SetTitleSize(0.07);
      c->Update();

      c->cd();
      p2->Draw();
      p2->cd();
      datratioX->SetTitle("");
      datratioX->GetXaxis()->SetTitle("p_{#mu} (MeV)");
      datratioX->GetYaxis()->SetTitle("Data/MC ");
      datratioX->GetYaxis()->SetTitleOffset(0.36);
      datratioX->GetXaxis()->SetRangeUser(0,5000);
      datratioX->GetYaxis()->SetRangeUser(0.01,1.99);
      datratioX->GetYaxis()->SetNdivisions(305);
      datratioX->GetYaxis()->SetLabelSize(0.1);
      datratioX->GetYaxis()->SetTitleSize(0.15);
      datratioX->GetXaxis()->SetLabelSize(0.12);
      datratioX->GetXaxis()->SetTitleSize(0.15);
      datratioX->Draw("e");
      TLine *line = new TLine(datratioX->GetXaxis()->GetXmin(), 1.0, 5000, 1.0);
      line->SetLineColor(kRed);
      line->Draw("same");
      c->Print((samp[s]+std::string("_p.pdf")).c_str());
      
      p1->cd();
      dathistY->SetTitle("");
      dathistY->Draw("e");
      ThStack->Draw("same");
      mchistY->Draw("e same");
      dathistY->Draw("e same");
      dathistY->SetTitle("");
      dathistY->GetXaxis()->SetTitle("cos #theta_{#mu}");
      dathistY->GetYaxis()->SetTitle("Events/1.0");
      dathistY->GetYaxis()->SetTitleOffset(0.8);
      dathistY->GetXaxis()->SetRangeUser(0.7,1.0);
      dathistY->GetYaxis()->SetRangeUser(0.1,1.25*dathistY->GetMaximum());
      dathistY->GetYaxis()->SetLabelSize(0.05);
      dathistY->GetYaxis()->SetTitleSize(0.07);
      c->Update();

      p2->cd();
      datratioY->SetTitle("");
      datratioY->GetXaxis()->SetTitle("cos #theta_{#mu}");
      datratioY->GetYaxis()->SetTitle("Data/MC ");
      datratioY->GetYaxis()->SetTitleOffset(0.36);
      datratioY->GetXaxis()->SetRangeUser(0.7,1.0);
      datratioY->GetYaxis()->SetRangeUser(0.01,1.99);
      datratioY->GetYaxis()->SetLabelSize(0.12);
      datratioY->GetYaxis()->SetNdivisions(305);
      datratioY->GetYaxis()->SetTitleSize(0.15);
      datratioY->GetXaxis()->SetLabelSize(0.12);
      datratioY->GetXaxis()->SetTitleSize(0.15);
      datratioY->Draw("e");
      TLine *line2 = new TLine(0.7, 1.0, 1.0, 1.0);
      line2->SetLineColor(kRed);
      line2->Draw("same");
      c->Print((samp[s]+std::string("_t.pdf")).c_str());
      
      std::cout << "done " << samp[s] << std::endl;
    }
}
