// Compare sigma variations to nominal
// .L compareSigVar.C+g // in root
// compareSigVar("filename.root" 
// ********************************

#include "TH1D.h"
#include "TH2Poly.h"
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


void compareSigVar(std::string filename) {

  //Get input file, make canvas and output file
  TCanvas *c = new TCanvas("canv", "canv", 1080, 1080);
  TFile *infile = TFile::Open(filename.c_str());
  
  c->SetTopMargin(0.07);
  c->SetBottomMargin(0.12);
  c->SetLeftMargin(0.18);
  c->SetRightMargin(0.23);
  
  //Get all entries in input file
  std::string dirname;

  // get directory names, ignore flux 
  dirname = std::string("FGD1_numuCC_0pi_EB_dial_C_nu");
  name = std::string("FGD1_numuCC_0pi_EB_dial_C_nu");
      
  //make -3,-1,0,1,3 polys
  TH2Poly *m3sig;
  TH2Poly *m1sig;
  TH2Poly *nom;
  TH2Poly *p1sig;
  TH2Poly *p3sig;

  name = dirname + "/" + dirname + "_-10";
  std::cout << "Getting " << name << std::endl;
  m3sig=(TH2Poly*)infile->Get(name.c_str())->Clone();
  
  name = dirname + "/" + dirname + "_-4";
  std::cout << "Getting " << name << std::endl;
  m1sig=(TH2Poly*)infile->Get(name.c_str())->Clone();

  name = dirname + "/" + dirname + "_2";
  std::cout << "Getting " << name << std::endl;
  nom=(TH2Poly*)infile->Get(name.c_str())->Clone();

  name = dirname + "/" + dirname + "_8";
  std::cout << "Getting " << name << std::endl;
  p1sig=(TH2Poly*)infile->Get(name.c_str())->Clone();

  name = dirname + "/" + dirname + "_15";
  std::cout << "Getting " << name << std::endl;
  p3sig=(TH2Poly*)infile->Get(name.c_str())->Clone();
  
  double maxz,minz;
  
  //Divide the plots for ratios and plot to a big old pdf
  m3sig->Divide(nom);
  m3sig->SetName((dirname + "_-3sigmaRatioToNom").c_str());
  m3sig->SetTitle((dirname + "_-3sigmaRatioToNom").c_str());
  c->Clear();
  maxz=m3sig->GetMaximum();
  minz=m3sig->GetMinimum();
  if (fabs(1-maxz)>fabs(1-minz))
    m3sig->GetZaxis()->SetRangeUser(1-fabs(1-maxz),1+fabs(1-maxz));
  else
    m3sig->GetZaxis()->SetRangeUser(1-fabs(1-minz),1+fabs(1-minz));
  m3sig->GetXaxis()->SetRangeUser(0,5000);
  m3sig->GetXaxis()->SetTitleOffset(1.1);
  m3sig->GetYaxis()->SetTitleOffset(1.1);
  m3sig->GetZaxis()->SetTitleOffset(1.5);
  m3sig->GetZaxis()->SetTitle("Ratio to Nominal");
  m3sig->Draw("colz");
  std::cout << "Drawn m3" << std::endl;
  // c->Print((outfilename).c_str());
  
  m1sig->Divide(nom);
  m1sig->SetName((dirname + "_-1sigmaRatioToNom").c_str());
  m1sig->SetTitle("");
  c->Clear();
  maxz=m1sig->GetMaximum();
  minz=m1sig->GetMinimum();
  if (fabs(1-maxz)>fabs(1-minz))
    m1sig->GetZaxis()->SetRangeUser(1-fabs(1-maxz),1+fabs(1-maxz));
  else
    m1sig->GetZaxis()->SetRangeUser(1-fabs(1-minz),1+fabs(1-minz));
  m1sig->GetXaxis()->SetRangeUser(0,5000);
  m1sig->GetXaxis()->SetTitleOffset(1.1);
  m1sig->GetYaxis()->SetTitleOffset(1.1);
  m1sig->GetZaxis()->SetTitleOffset(1.5);
  m1sig->GetZaxis()->SetTitle("Ratio to Nominal");
  m1sig->Draw("colz");
  std::cout << "Drawn m1" << std::endl;
  c->Print("EbNuCM1Ratio.pdf");
  
  p1sig->Divide(nom);
  p1sig->SetName((dirname + "_+1sigmaRatioToNom").c_str());
  p1sig->SetTitle("");
  c->Clear();
  maxz=p1sig->GetMaximum();
  minz=p1sig->GetMinimum();
  if (fabs(1-maxz)>fabs(1-minz))
    p1sig->GetZaxis()->SetRangeUser(1-fabs(1-maxz),1+fabs(1-maxz));
  else
    p1sig->GetZaxis()->SetRangeUser(1-fabs(1-minz),1+fabs(1-minz));
  p1sig->GetXaxis()->SetRangeUser(0,5000);
  p1sig->GetXaxis()->SetTitleOffset(1.1);
  p1sig->GetYaxis()->SetTitleOffset(1.1);
  p1sig->GetZaxis()->SetTitleOffset(1.5);
  p1sig->GetZaxis()->SetTitle("Ratio to Nominal");
  p1sig->Draw("colz");
  c->Print("EbNuCP1Ratio.pdf");
  std::cout << "Drawn p1" << std::endl;
  
  p3sig->Divide(nom);
  p3sig->SetName((dirname + "_+3sigmaRatioToNom").c_str());
  p3sig->SetTitle((dirname + "_+3sigmaRatioToNom").c_str());
  c->Clear();
  maxz=p3sig->GetMaximum();
  minz=p3sig->GetMinimum();
  if (fabs(1-maxz)>fabs(1-minz))
    p3sig->GetZaxis()->SetRangeUser(1-fabs(1-maxz),1+fabs(1-maxz));
  else
    p3sig->GetZaxis()->SetRangeUser(1-fabs(1-minz),1+fabs(1-minz));
  p3sig->GetXaxis()->SetRangeUser(0,5000);
  p3sig->GetXaxis()->SetTitleOffset(1.1);
  p3sig->GetYaxis()->SetTitleOffset(1.1);
  p3sig->GetZaxis()->SetTitleOffset(1.5);
  p3sig->GetZaxis()->SetTitle("Ratio to Nominal");
  p3sig->Draw("colz");
  std::cout << "Drawn p3"<< std::endl;
  //c->Print((outfilename).c_str());
  
  delete m3sig;
  delete m1sig;
  delete nom;
  delete p1sig;
  delete p3sig;
  
  //    if(count2 >42)
  //break;
}

