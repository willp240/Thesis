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
#include "TColor.h"

void detbins(std::string file1, std::string file2) {

  TCanvas *c = new TCanvas("canv", "canv", 1580, 1380);
  TFile *file_1 = TFile::Open(file1.c_str());
  TFile *file_2 = TFile::Open(file2.c_str());
  
  c->SetTopMargin(0.05);
  c->SetBottomMargin(0.18);
  c->SetLeftMargin(0.18);
  c->SetRightMargin(0.22);

  file_1->cd();
  gDirectory->cd("BinVariations");

  Bin_5_NEvent_Spread->Print((std::string("detbin_allsysts5")+std::string(".pdf")).c_str());
  Bin_253_NEvent_Spread->Print((std::string("detbin_allsysts253")+std::string(".pdf")).c_str());
  Bin_430_NEvent_Spread->Print((std::string("detbin_allsysts430")+std::string(".pdf")).c_str());
  Bin_535_NEvent_Spread->Print((std::string("detbin_allsysts535")+std::string(".pdf")).c_str());
  
  Bin_258_NEvent_Spread->Print((std::string("detbin_allsysts258")+std::string(".pdf")).c_str());
  Bin_97_NEvent_Spread->Print((std::string("detbin_allsysts97")+std::string(".pdf")).c_str());
  Bin_103_NEvent_Spread->Print((std::string("detbin_allsysts103")+std::string(".pdf")).c_str());
  Bin_570_NEvent_Spread->Print((std::string("detbin_allsysts570")+std::string(".pdf")).c_str());

  file_2->cd();
  gDirectory->cd("BinVariations");

  Bin_258_NEvent_Spread->Print((std::string("detbin_nopisi258")+std::string(".pdf")).c_str());
  Bin_97_NEvent_Spread->Print((std::string("detbin_nopisi97")+std::string(".pdf")).c_str());
  Bin_103_NEvent_Spread->Print((std::string("detbin_nopisi103")+std::string(".pdf")).c_str());
  Bin_570_NEvent_Spread->Print((std::string("detbin_nopisi570")+std::string(".pdf")).c_str());
}
