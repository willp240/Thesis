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


void PlotHistos(std::string file0, int run) {

  TCanvas *c = new TCanvas("canv", "canv", 1080, 1080);
  TFile *infile = TFile::Open(file0.c_str());

 
  c->SetTopMargin(0.07);
  c->SetBottomMargin(0.11);
  c->SetLeftMargin(0.10);
  c->SetRightMargin(0.18);
  c->Divide(2,1);
  c->Print((std::string("plots")+std::string(".pdf[")).c_str());

  gStyle->SetOptStat(0);
  gPad->SetLogz();

  for(int i_eve=0; i_eve<40; i_eve++){

    for(int i_cam=0; i_cam<4; i_cam++){


      TString hrawname;
      hrawname = Form("hRaw_run%d_cam%d_eve%d",run,i_cam,i_eve);
      TH2D *hraw = (TH2D*)(infile->Get(hrawname));
      if(!hraw)
	std::cout << "No " << hrawname << std::endl;
      c->cd(1);
      hraw->SetName(hrawname);
      hraw->SetTitle(hrawname);
      hraw->GetZaxis()->SetRangeUser(2000,5000);
      hraw->Draw("colz");
            
      TString hsubname;
      hsubname = Form("hBiasSub_run%d_cam%d_eve%d",run,i_cam,i_eve);
      TH2D *hsub = (TH2D*)(infile->Get(hsubname));
      if(!hsub)
	std::cout << "No " << hsubname << std::endl;
      c->cd(2);
      hsub->SetName(hsubname);
      hsub->SetTitle(hsubname);
      hsub->GetZaxis()->SetRangeUser(-50,250);
      hsub->Draw("colz");
      c->Print((std::string("plots")+std::string(".pdf")).c_str());
      std::cout << "Done for eve " << i_cam << " event " << i_eve << std::endl;
    }
  }

  
	  
  std::cout << "**********"<< std::endl;
   
  c->Print((std::string("plots")+std::string(".pdf]")).c_str());
}
