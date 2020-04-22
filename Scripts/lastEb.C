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


void lastEb(std::string file0) {

  TFile *zerofile = TFile::Open(file0.c_str());
  TFile *outfile = new TFile((file0+"_corr.root").c_str(), "RECREATE");
 
  TVectorD* FittedParameters0  = (TVectorD*)zerofile->Get("FittedParameters0");
  TMatrixT<double>* FittedCovariance0 = (TMatrixT<double>*)zerofile->Get("FittedCovariance0");

  int sz = FittedParameters0->GetNoElements();

  (*FittedParameters0)(sz-1)= -2.39;
  (*FittedCovariance0)(sz-1,sz-1)=23.81;
  (*FittedCovariance0)(sz-2,sz-2)=30.58;

  outfile->cd();
  FittedParameters0->Write("FittedParameters0");
  FittedCovariance0->Write("FittedCovariance0");
  

}
