#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include "TFile.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TGaxis.h"
#include "TMarker.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TMath.h"

void CompareContours1Ddcp(char *cont1, char *cont2, std::string firstlegend,std::string secondlegend,std::string outname="",std::string cont3="none",std::string thirdlegend="none")
{
  std::cout << "in here" << std::endl;

  TGaxis::SetMaxDigits(3) ;
  gStyle->SetOptStat(0);
  // Get MaCh3 contours
  std::cout<<"Getting contour 1"<<endl;
  TFile *f_1 = new TFile(cont1,"open");

  TH1D* mdcp_input=(TH1D*)f_1->Get("h_dcp");
  if(!mdcp_input)mdcp_input=(TH1D*)f_1->Get("dcp_posterior");
  TH1D* mdcp;

  //TH1D* mdcp=(TH1D*)f_1->Get("h_dcp_chi2");
  if(mdcp_input){
    cout<<"Found contour 1"<<endl;
    mdcp=(TH1D*)mdcp_input->Clone("mdcp");
    double maxbin = mdcp_input->GetMaximumBin();
    double maxb = mdcp_input->GetBinContent(maxbin);
    for (int iBin = 1; iBin <= mdcp_input->GetNbinsX(); iBin++){
      double bincont = mdcp_input->GetBinContent(iBin);
      if (bincont == 0) bincont = 1;
      double likelihood = -2 * TMath::Log(bincont/maxb);
      mdcp->SetBinContent(iBin, likelihood);
    }
    mdcp=mdcp_input;
    mdcp->SetLineWidth(2);
    mdcp->SetLineColor(kBlue);
    mdcp->GetXaxis()->SetTitle("#delta_{CP}");
    mdcp->GetYaxis()->SetTitleOffset(1.0);
    mdcp->GetYaxis()->SetTitle("Posterior probability");
    mdcp->GetXaxis()->SetTitleOffset(1.0);
    mdcp->SetFillStyle(0);
  }
  else     cout<<"Didn't find contour 1"<<endl;

  //!!Make VALOR Contour 2
  TFile *f_2 = new TFile(cont2,"open");

  TH1D* vdcp_input=(TH1D*)f_2->Get("h_dcp");
  if(!vdcp_input)vdcp_input=(TH1D*)f_1->Get("dcp_posterior");
  TH1D* vdcp;

  //TH1D* mdcp=(TH1D*)f_1->Get("h_dcp_chi2");
  if(vdcp_input){
    cout<<"Found contour 2"<<endl;
    vdcp=(TH1D*)vdcp_input->Clone("vdcp");
    double maxbin = vdcp_input->GetMaximumBin();
    double maxb = vdcp_input->GetBinContent(maxbin);
    for (int iBin = 1; iBin <= vdcp_input->GetNbinsX(); iBin++){
      double bincont = vdcp_input->GetBinContent(iBin);
      if (bincont == 0) bincont = 1;
      double likelihood = -2 * TMath::Log(bincont/maxb);
      vdcp->SetBinContent(iBin, likelihood);
    }
    vdcp=vdcp_input;
    vdcp->GetXaxis()->SetTitle("#delta_{CP}");
    vdcp->GetYaxis()->SetTitle("#Delta#chi^{2}");
    vdcp->GetXaxis()->SetTitleOffset(1.0);
    vdcp->SetFillStyle(0);
  }
  else     cout<<"Didn't find contour 2"<<endl;

  TH1D* mdcp3;
  if(cont3.compare("none")!=0) {

    std::cout<<"Getting contour 3"<<endl;
    TFile *f_3 = new TFile(cont3.c_str(),"open");
  
    TH1D* mdcp3_input=(TH1D*)f_3->Get("h_dcp");
    if(!mdcp3_input)mdcp3_input=(TH1D*)f_3->Get("dcp_posterior");
  
    //TH1D* mdcp=(TH1D*)f_1->Get("h_dcp_chi2");
    if(mdcp3_input){
      cout<<"Found contour 3"<<endl;
      mdcp3=(TH1D*)mdcp3_input->Clone("mdcp3");
      double maxbin = mdcp3_input->GetMaximumBin();
      double maxb = mdcp3_input->GetBinContent(maxbin);
      for (int iBin = 1; iBin <= mdcp3_input->GetNbinsX(); iBin++){
        double bincont = mdcp3_input->GetBinContent(iBin);
        if (bincont == 0) bincont = 1;
        double likelihood = -2 * TMath::Log(bincont/maxb);
        mdcp3->SetBinContent(iBin, likelihood);
      }
      mdcp3=mdcp3_input;
      mdcp3->SetLineWidth(2);
      mdcp3->SetLineColor(kBlack);
      mdcp3->GetXaxis()->SetTitle("#delta_{CP}");
      mdcp3->GetYaxis()->SetTitleOffset(1.0);
      mdcp3->GetYaxis()->SetTitle("Posterior probability");
      mdcp3->GetXaxis()->SetTitleOffset(1.0);
      mdcp3->SetFillStyle(0);
    }
    else     cout<<"Didn't find contour 3"<<endl;

}

  if(vdcp) vdcp->SetLineWidth(2);
  if(vdcp) vdcp->SetLineColor(kRed);

  // Legend
  //TLegend *leg = new TLegend(0.5,0.6,0.89,0.89);//right hand side
  //TLegend *leg = new TLegend(0.5,0.5,0.95,0.89);//For Poly
  TLegend *leg = new TLegend(0.5,0.55,0.82,0.86);//Orig
  leg->SetNColumns(1);
  leg->SetFillColor(kWhite);
  leg->SetFillStyle(0);
  leg->SetLineColor(0);
  //leg->SetTextFont(132);
  if (mdcp) leg->AddEntry(mdcp,firstlegend.c_str(),"l");
  if (vdcp) leg->AddEntry(vdcp,secondlegend.c_str(),"l");
  if (mdcp3) leg->AddEntry(mdcp3,thirdlegend.c_str(),"l");

  // Draw
  TCanvas* c1=new TCanvas();
  c1->SetLeftMargin(0.13);
  // TH2D *dummy;
  // if(!do23){
  //   dummy= new TH2D("dummy","",50,0,0.08,50,-TMath::Pi(),TMath::Pi());
  //   dummy->GetXaxis()->SetTitle("#delta_{CP}");
  //   dummy->GetYaxis()->SetTitle("#delta_{CP}");
  //   dummy->Draw();
  //   m90_0->Draw("L same");
  // }
  // else{
  //   cout<<"doing 23 "<<"doing IH?"<<doIH<<endl;
  //   if(!doIH) dummy= new TH2D("dummy","",50,0.4,0.65,100,0.0023,0.003);
  //   else dummy= new TH2D("dummy","",50,0.4,0.65,100,-0.0028,-0.002);
  //   dummy->GetXaxis()->SetTitle("sin^{2}#theta_{23}");
  //   dummy->GetYaxis()->SetTitle("#Delta m^{2}_{32} (eV^{2})");
  //   dummy->Draw();
  //   if (m90_0) m90_0->Draw("L same");
  // }

  mdcp->Scale(1./mdcp->Integral(),"width");
  vdcp->Scale(1./vdcp->Integral(),"width");
  if (mdcp3) mdcp3->Scale(1./mdcp3->Integral(),"width");

  if(vdcp->GetMaximum()>mdcp->GetMaximum())mdcp->GetYaxis()->SetRangeUser(0,vdcp->GetMaximum()+0.1);
  else mdcp->GetYaxis()->SetRangeUser(0,mdcp->GetMaximum()+0.1);
  if (mdcp) mdcp->Draw("same");
  if (vdcp) vdcp->Draw("same");
  if (mdcp3) mdcp3->Draw("same");
  leg->Draw("same");
  c1->SaveAs((outname+".pdf").c_str());
  c1->SaveAs((outname+".png").c_str());
  c1->SaveAs((outname+".root").c_str());
}
