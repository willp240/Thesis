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

void skerate(std::string file0) {

  TFile *file_ = TFile::Open(file0.c_str()); 
  
  TString spec[5] = {"numu", "numubar", "nue", "nuebar", "nue1pi"};

  TH1D *total = new TH1D("","", 1000,0,1000);
  double totvec[2000];

  for(int i=0; i<2000; i++)
    totvec[i]=0.0;

  for(int s=0; s<1; s++){

    TH1D *rates = new TH1D("","", 1000,0,1000); 

    for(int i=0; i<2000; i++){
      
      TString name = Form("sk_%s%d",spec[s].Data(),i);
      //      std::cout << name << std::endl;
      TH1D* hist = (TH1D*)file_->Get(name)->Clone();

      rates->Fill(hist->Integral());
      totvec[i]+=hist->Integral();
    }
    rates->Fit("gaus");
    std::cout << spec[s] << " Mean: " << rates->GetMean() << ", RMS: " << rates->GetRMS() << std::endl;
  }

  for(int i=0; i<2000; i++){
    total->Fill(totvec[i]);
  }

  total->Fit("gaus");
  std::cout < <"Total Mean: " << total->GetMean() << ", RMS: " << total->GetRMS() << std::endl;
  
}
