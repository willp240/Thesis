#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include "TFile.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TMarker.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TStyle.h"

void CompareContours(char *cont1D, char *cont2D,std::string firstlegend="1D",std::string secondlegend="2D",bool app=false,std::string outname="comparedcontours1d2d",bool dowrc=false,bool doIH=false,bool doconf=false)
{
  std::cout<<"Getting first contours"<<std::endl;
  // Get MaCh3 contours
  TFile *f_1D = new TFile(cont1D,"open");
  std::vector<TGraph*> m90;
  std::vector<TGraph*> m68;
  std::string prefix90;
  std::string prefix68;
  if(doconf){
    prefix90="g90_conf_";
    prefix68="g68_conf_";
  }
  else{
    prefix90="g90_";
    prefix68="g68_";
  }

  for(int iCont=0;true;iCont++){
    std::ostringstream name;
    name << prefix90<<iCont;
    if((TGraph*)f_1D->Get(name.str().c_str()))m90.push_back((TGraph*)f_1D->Get(name.str().c_str()));
    else break;
  }
  std::cout<<"Found "<<m90.size()<<" first file 90% contours"<<std::endl;
  for(int iCont=0;true;iCont++){
    std::ostringstream name;
    name << prefix68<<iCont;
    if((TGraph*)f_1D->Get(name.str().c_str()))m68.push_back((TGraph*)f_1D->Get(name.str().c_str()));
    else break;
  }
  std::cout<<"Found "<<m68.size()<<" first file 68% contours"<<std::endl;

  TGraph *mbf = (TGraph*)f_1D->Get("Bestfit");

  std::cout<<"Setting plot style for first file"<<std::endl;
  for(unsigned i90=0;i90<m90.size();i90++){
    m90.at(i90)->SetLineWidth(2);
    m90.at(i90)->SetLineColor(kBlue);
    m90.at(i90)->GetXaxis()->SetTitle("sin^{2}#theta_{23}");
    m90.at(i90)->GetYaxis()->SetTitleOffset(1.2);
    m90.at(i90)->GetYaxis()->SetTitle("#Delta m^{2}_{32} (eV^{2})");
  }
  for(unsigned i68=0;i68<m68.size();i68++){
    m68.at(i68)->SetLineWidth(2);
    m68.at(i68)->SetLineStyle(2);
    m68.at(i68)->SetLineColor(kBlue);
    m68.at(i68)->GetXaxis()->SetTitle("sin^{2}#theta_{23}");
    m68.at(i68)->GetYaxis()->SetTitleOffset(1.2);
    m68.at(i68)->GetYaxis()->SetTitle("#Delta m^{2}_{32} (eV^{2})");
  }

  std::cout<<"Getting 2D contours"<<std::endl;

  // // Get Valor contours
  TFile *f_2D = new TFile(cont2D,"open");

  std::vector<TGraph*> v90;
  std::vector<TGraph*> v68;

  for(int iCont=0;true;iCont++){
    std::ostringstream name;
    name << prefix90<<iCont;
    if((TGraph*)f_2D->Get(name.str().c_str()))v90.push_back((TGraph*)f_2D->Get(name.str().c_str()));
    else break;
  }
  std::cout<<"Found "<<v90.size()<<" second file 90% contours"<<std::endl;

  for(int iCont=0;true;iCont++){
    std::ostringstream name;
    name << prefix68<<iCont;
    if((TGraph*)f_2D->Get(name.str().c_str()))v68.push_back((TGraph*)f_2D->Get(name.str().c_str()));
    else break;
  }
  std::cout<<"Found "<<v68.size()<<" second file 68% contours"<<std::endl;

  TGraph *vbf = (TGraph*)f_2D->Get("Bestfit");

  std::cout<<"Setting 2D plotting style"<<std::endl;
  for(unsigned i90=0;i90<v90.size();i90++){
    v90.at(i90)->SetLineWidth(2);
    v90.at(i90)->SetLineColor(kRed);
    v90.at(i90)->GetXaxis()->SetTitle("sin^{2}#theta_{23}");
    v90.at(i90)->GetYaxis()->SetTitleOffset(1.2);
    v90.at(i90)->GetYaxis()->SetTitle("#Delta m^{2}_{32} (eV^{2})");
    v90.at(i90)->GetXaxis()->SetTitle("sin^{2}#theta_{23}");
    v90.at(i90)->GetYaxis()->SetTitle("#Delta m^{2}_{32} (eV^{2})");

  }
  std::cout<<"blah1"<<std::endl;//!!
  for(unsigned i68=0;i68<v68.size();i68++){
    v68.at(i68)->SetLineWidth(2);
    v68.at(i68)->SetLineStyle(2);
    v68.at(i68)->SetLineColor(kRed);
    v68.at(i68)->GetXaxis()->SetTitle("sin^{2}#theta_{23}");
    v68.at(i68)->GetYaxis()->SetTitleOffset(1.2);
    v68.at(i68)->GetYaxis()->SetTitle("#Delta m^{2}_{32} (eV^{2})");
  }  
  if (vbf) vbf->SetMarkerColor(kRed);
  if (vbf) vbf->SetMarkerStyle(22);

  std::cout<<"Setting up legend"<<std::endl;

  // Legend
  TLegend *leg; /*
  if(!app)leg= new TLegend(0.2,0.5,0.85,0.89);//poly
  //else leg= new TLegend(0.68,0.2,0.89,0.89);//poly
  else leg= new TLegend(0.45,0.15,0.99,0.89);//poly*/
 
  /*  if(!app)leg= new TLegend(0.3,0.55,0.85,0.89);//orig
  //else leg= new TLegend(0.68,0.2,0.89,0.89);//orig
  else leg= new TLegend(0.5,0.2,0.85,0.89);//orig*/

  if(!app)leg= new TLegend(0.3,0.54,0.85,0.88);//eb
  //else leg= new TLegend(0.68,0.2,0.89,0.89);//eb
  else leg= new TLegend(0.5,0.2,0.85,0.89);//eb

  leg->SetNColumns(1);
  leg->SetFillColor(kWhite);
  leg->SetLineColor(kBlack);
  //leg->SetTextFont(132);
  leg->SetBorderSize(0.);
  if (m90.at(0)) leg->AddEntry(m90.at(0),(firstlegend+" 90%").c_str(),"l");
  if (m68.at(0)) leg->AddEntry(m68.at(0),(firstlegend+" 68%").c_str(),"l");
  if (mbf) leg->AddEntry(mbf,(firstlegend+" Best Fit").c_str(),"p");
  if (v90.at(0)) leg->AddEntry(v90.at(0),(secondlegend+" 90%").c_str(),"l");
  if (v68.at(0)) leg->AddEntry(v68.at(0),(secondlegend+" 68%").c_str(),"l");
  if (vbf) leg->AddEntry(vbf,(secondlegend+" Best Fit").c_str(),"p");

  std::cout<<"Drawing"<<std::endl;
  // Draw
  TH2D *dummy;
  double yupper,ylower,xupper,xlower;
  TString xtitle,ytitle;
  if(app){
    xtitle="sin^{2}#theta_{13}";
    ytitle="#delta_{CP}";
    yupper=TMath::Pi();
    ylower=-TMath::Pi();
    if(dowrc){
      xupper=0.035;
      xlower=0.015;
    }
    else{
      xupper=0.075;
      xlower=0.01;
    }
  }
  else{
    xtitle="sin^{2}#theta_{23}";
    ytitle="#Delta M^{2}_{23}";
    xupper=0.65;
    xlower=0.4;
    if(!doIH){
      yupper=0.003;
      ylower=0.0023;
    }
    else{
      yupper=-0.0018;
      ylower=-0.0028;
    }
  }
    
  dummy= new TH2D("dummy","",50,xlower,xupper,50,ylower,yupper);
  dummy->GetXaxis()->SetTitle(xtitle);
  dummy->GetYaxis()->SetTitle(ytitle);
  gStyle->SetOptStat(0);

  //  if(!app) dummy->GetYaxis()->SetTitleOffset(1.4);
  //if(!app) TGaxis::SetMaxDigits(2);

  //if(dowrc) dummy->SetTitle("MaCh3 Asimov A Joint Fit, wRC");
  //else dummy->SetTitle("MaCh3 Asimov A Joint Fit, woRC");
  if(dowrc) dummy->SetTitle("");
  else dummy->SetTitle("");

  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleFillColor(0);

  dummy->GetXaxis()->SetTitleOffset(1.1);

  dummy->Draw();

  TCanvas* c1=new TCanvas("c1");
  c1->SetLeftMargin(0.13);
  dummy->Draw("");
  for(unsigned im90=0;im90<m90.size();im90++)m90.at(im90)->Draw("L same");
  for(unsigned im68=0;im68<m68.size();im68++)m68.at(im68)->Draw("L same");
  if (mbf) mbf->Draw("P same");
  for(unsigned iv90=0;iv90<v90.size();iv90++)v90.at(iv90)->Draw("L same");
  for(unsigned iv68=0;iv68<v68.size();iv68++)v68.at(iv68)->Draw("L same");
  if (vbf) vbf->Draw("P same");

  leg->Draw("same");
  c1->SaveAs((outname+".pdf").c_str());
  c1->SaveAs((outname+".png").c_str());
  c1->SaveAs((outname+".root").c_str());
}
