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


void make4Plots(std::string file0, int event, int run) {

  TCanvas *c = new TCanvas("CombinedCanvas", "CombinedCanvas", 600, 630);
  TFile *infile = TFile::Open(file0.c_str());
 
  //  c->SetTopMargin(0.07);
  //c->SetBottomMargin(0.11);
  //c->SetLeftMargin(0.10);
  //c->SetRightMargin(0.18);

  c->SetLeftMargin(0);
  c->SetRightMargin(0);
  c->SetBottomMargin(0);
  c->SetTopMargin(0);
  c->Divide(2,2,0,0);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(55);
  gPad->SetTicks();
  ////////////////////////////////////////////CCD0
  TString hname0;
  hname0 = Form("hBiasSub_run%d_cam0_eve%d",run,event);
  TH2D *hsub0 = (TH2D*)(infile->Get(hname0));
  if(!hsub0)
    std::cout << "No " << hname0 << std::endl;
  c->cd(1);
  gStyle->SetOptStat(0);
  gPad->SetTicks();

  hsub0->SetName(hname0);
  hsub0->SetTitle("");
  hsub0->GetZaxis()->SetRangeUser(-100,2000);
  hsub0->GetXaxis()->SetLabelSize(0);
  hsub0->GetYaxis()->SetLabelSize(0);
  hsub0->Draw("colz");
  gPad->Update();
  TPaletteAxis *palette = (TPaletteAxis*)hsub0->GetListOfFunctions()->FindObject("palette");
  if(!palette)
    std::cout << "No palette" << std::endl;
  palette->SetX1NDC(2.9);
  palette->SetX2NDC(3.9);
  gStyle->SetOptStat(0);
  hsub0->Draw("colz");
 
  TF1 *cam0source1up = new TF1("cam0source1up","20+sqrt(100-(x-135)^2)",125,145);
  cam0source1up->SetLineColor(kWhite);
  cam0source1up->Draw("same");
  TF1 *cam0source1down = new TF1("cam0source1down","20-sqrt(100-(x-135)^2)",125,145);
  cam0source1down->SetLineColor(kWhite);
  cam0source1down->Draw("same");
  TF1 *cam0source2up = new TF1("cam0source2up","42+sqrt(100-(x-355)^2)",345,365);
  cam0source2up->SetLineColor(kWhite);
  cam0source2up->Draw("same");
  TF1 *cam0source2down = new TF1("cam0source2down","42-sqrt(100-(x-355)^2)",345,365);
  cam0source2down->SetLineColor(kWhite);
  cam0source2down->Draw("same");
  TF1 *cam0source3up = new TF1("cam0source3up","150+sqrt(100-(x-350)^2)",340,360);
  cam0source3up->SetLineColor(kWhite);
  cam0source3up->Draw("same");
  TF1 *cam0source3down = new TF1("cam0source3down","150-sqrt(100-(x-350)^2)",340,360);
  cam0source3down->SetLineColor(kWhite);
  cam0source3down->Draw("same");

  ////////////////////////////////////////////CCD1
  TString hname1;
  hname1 = Form("hBiasSub_run%d_cam1_eve%d",run,event);
  TH2D *hsub1 = (TH2D*)(infile->Get(hname1));
  if(!hsub1)
    std::cout << "No " << hname1 << std::endl;
  c->cd(2);
  gStyle->SetOptStat(0);
  gPad->SetTicks();

  hsub1->SetName(hname1);
  hsub1->SetTitle("");
  hsub1->GetZaxis()->SetRangeUser(-100,2000);
  hsub1->GetXaxis()->SetLabelSize(0);
  hsub1->GetYaxis()->SetLabelSize(0);
  hsub1->Draw("colz");
  gPad->Update();
  TPaletteAxis *palette1 = (TPaletteAxis*)hsub1->GetListOfFunctions()->FindObject("palette");
  palette1->SetX1NDC(2.9);
  palette1->SetX2NDC(3.9);
  hsub1->Draw("colz");

  TF1 *cam1source1up = new TF1("cam1source1up","30+sqrt(100-(x-35)^2)",25,45);
  cam1source1up->SetLineColor(kWhite);
  cam1source1up->Draw("same");
  TF1 *cam1source1down = new TF1("cam1source1down","30-sqrt(100-(x-35)^2)",25,45);
  cam1source1down->SetLineColor(kWhite);
  cam1source1down->Draw("same");
  TF1 *cam1source2up = new TF1("cam1source2up","140+sqrt(100-(x-30)^2)",20,40);
  cam1source2up->SetLineColor(kWhite);
  cam1source2up->Draw("same");
  TF1 *cam1source2down = new TF1("cam1source2down","140-sqrt(100-(x-30)^2)",20,40);
  cam1source2down->SetLineColor(kWhite);
  cam1source2down->Draw("same");
  TF1 *cam1source3up = new TF1("cam1source3up","10+sqrt(100-(x-225)^2)",215,235);
  cam1source3up->SetLineColor(kWhite);
  cam1source3up->Draw("same");
  TF1 *cam1source3down = new TF1("cam1source3down","10-sqrt(100-(x-225)^2)",215,235);
  cam1source3down->SetLineColor(kWhite);
  cam1source3down->Draw("same");


  ////////////////////////////////////////////CCD2
  TString hname2;
  hname2 = Form("hBiasSub_run%d_cam2_eve%d",run,event);
  TH2D *hsub2 = (TH2D*)(infile->Get(hname2));
  if(!hsub2)
    std::cout << "No " << hname2 << std::endl;
  c->cd(3);
  gStyle->SetOptStat(0);
  gPad->SetTicks();

  hsub2->SetName(hname2);
  hsub2->SetTitle("");
  hsub2->GetZaxis()->SetRangeUser(-100,5000);
  hsub2->GetXaxis()->SetLabelSize(0);
  hsub2->GetYaxis()->SetLabelSize(0);
  hsub2->Draw("colz");
  gPad->Update();
  TPaletteAxis *palette2 = (TPaletteAxis*)hsub2->GetListOfFunctions()->FindObject("palette");
  palette2->SetX1NDC(2.9);
  palette2->SetX2NDC(3.9);
  hsub2->Draw("colz");

  TF1 *cam2source1up = new TF1("cam2source1up","340+sqrt(100-(x-130)^2)",120,140);
  cam2source1up->SetLineColor(kWhite);
  cam2source1up->Draw("same");
  TF1 *cam2source1down = new TF1("cam2source1down","340-sqrt(100-(x-130)^2)",120,140);
  cam2source1down->SetLineColor(kWhite);
  cam2source1down->Draw("same");
  TF1 *cam2source2up = new TF1("cam2source2up","353+sqrt(100-(x-345)^2)",335,355);
  cam2source2up->SetLineColor(kWhite);
  cam2source2up->Draw("same");
  TF1 *cam2source2down = new TF1("cam2source2down","353-sqrt(100-(x-345)^2)",335,355);
  cam2source2down->SetLineColor(kWhite);
  cam2source2down->Draw("same");
  TF1 *cam2source3up = new TF1("cam2source3up","210+sqrt(100-(x-345)^2)",335,355);
  cam2source3up->SetLineColor(kWhite);
  cam2source3up->Draw("same");
  TF1 *cam2source3down = new TF1("cam2source3down","210-sqrt(100-(x-345)^2)",335,355);
  cam2source3down->SetLineColor(kWhite);
  cam2source3down->Draw("same");


  ////////////////////////////////////////////CCD3
  TString hname3;
  hname3 = Form("hBiasSub_run%d_cam3_eve%d",run,event);
  TH2D *hsub3 = (TH2D*)(infile->Get(hname3));
  if(!hsub3)
    std::cout << "No " << hname3 << std::endl;
  c->cd(4);
  gStyle->SetOptStat(0);
  gPad->SetTicks();

  hsub3->SetName(hname3);
  hsub3->SetTitle("");
  hsub3->GetZaxis()->SetRangeUser(-100,10000);
  hsub3->GetXaxis()->SetLabelSize(0);
  hsub3->GetYaxis()->SetLabelSize(0);
  hsub3->Draw("colz");
  gPad->Update();
  TPaletteAxis *palette3 = (TPaletteAxis*)hsub3->GetListOfFunctions()->FindObject("palette");
  palette3->SetX1NDC(2.9);
  palette3->SetX2NDC(3.9);
  hsub3->Draw("colz");     
  
  TF1 *cam3source1up = new TF1("cam3source1up","180+sqrt(100-(x-20)^2)",10,30);
  cam3source1up->SetLineColor(kWhite);
  cam3source1up->Draw("same");
  TF1 *cam3source1down = new TF1("cam3source1down","180-sqrt(100-(x-20)^2)",10,30);
  cam3source1down->SetLineColor(kWhite);
  cam3source1down->Draw("same");
  TF1 *cam3source2up = new TF1("cam3source2up","350+sqrt(100-(x-20)^2)",10,30);
  cam3source2up->SetLineColor(kWhite);
  cam3source2up->Draw("same");
  TF1 *cam3source2down = new TF1("cam3source2down","350-sqrt(100-(x-20)^2)",10,30);
  cam3source2down->SetLineColor(kWhite);
  cam3source2down->Draw("same");
  TF1 *cam3source3up = new TF1("cam3source3up","325+sqrt(100-(x-230)^2)",220,240);
  cam3source3up->SetLineColor(kWhite);
  cam3source3up->Draw("same");
  TF1 *cam3source3down = new TF1("cam3source3down","325-sqrt(100-(x-230)^2)",220,240);
  cam3source3down->SetLineColor(kWhite);
  cam3source3down->Draw("same");


  std::cout << "**********"<< std::endl;

  TFile *outfile = TFile::Open("sparkyplots.root","RECREATE");
  c->Write();
  hsub0->Write();
  hsub1->Write();
  hsub2->Write();
  hsub3->Write();

}
