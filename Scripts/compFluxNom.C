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


void compFluxNom(std::string file1, std::string file2)
{

  TFile *f1 = TFile::Open(file1.c_str());
  TFile *f2 = TFile::Open(file2.c_str());
  TH2Poly* h1;
  TH2Poly* h2;
  TH2Poly* h3;
  TH2PolyBin* bin1;
  TH2PolyBin* bin2;
  
  h1 = (TH2Poly*)f1->Get("MC_FGD1_NuMuBkg_CC0pi_in_AntiNu_Mode");
  h2 = (TH2Poly*)f2->Get("MC_FGD1_NuMuBkg_CC0pi_in_AntiNu_Mode");
  h3 = (TH2Poly*)f2->Get("MC_FGD1_NuMuBkg_CC0pi_in_AntiNu_Mode");
  h3->Reset("");
  double maxx1, maxx2, maxy1, maxy2, minx1, minx2, miny1, miny2;

  for(int i=1; i<=h1->GetNumberOfBins(); i++){
  
    bin1 = (TH2PolyBin*)(((TH2Poly*)h1)->GetBins()->At(i-1)->Clone());
    bin2 = (TH2PolyBin*)(((TH2Poly*)h2)->GetBins()->At(i-1)->Clone());

    maxx1 = bin1->GetXMax();
    maxx2 = bin2->GetXMax();
    minx1 = bin1->GetXMin();
    minx2 = bin2->GetXMin();
    maxy1 = bin1->GetYMax();
    maxy2 = bin2->GetYMax();
    miny1 = bin1->GetYMin();
    miny2 = bin2->GetYMin();

    if(maxx1!=maxx2 || minx1!=minx2 || maxy1!=maxy2 || miny1!=miny2){
      std::cout << "Bin1: " << minx1 << "-" << maxx1 << ", " << miny1 << "-" << maxy1 << std::endl; 
      std::cout << "Bin2: " << minx2 << "-" << maxx2 << ", " << miny2 << "-" << maxy2 << std::endl;
      i=10000;
    }
    double binval1 = h1->GetBinContent(i);
    double binval2 = h2->GetBinContent(i);
    h3->SetBinContent(i,binval1/binval2);
    std::cout << "Bin " << i << " binval1: " << binval1 << " binval2: " << binval2 << " Ratio: " << binval1/binval2 << std::endl;
  }

  h3->GetZaxis()->SetRangeUser(0.7,1.3);
  h3->GetZaxis()->SetTitle("Ratio");
  h3->SetTitle("MC_FGD1_NuMuBkg_CC0pi_in_AntiNu_Mode Ratio av2:av6");
  h3->Draw("colz");
}
