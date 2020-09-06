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

void detbinsNoSI(std::string file1, std::string file2) {

  TCanvas *c = new TCanvas("canv", "canv", 1580, 1380);
  TFile *file_1 = TFile::Open(file1.c_str());
  TFile *file_2 = TFile::Open(file2.c_str());
  
  c->SetTopMargin(0.05);
  c->SetBottomMargin(0.18);
  c->SetLeftMargin(0.18);
  c->SetRightMargin(0.22);
 
  file_2->cd();
  gDirectory->cd("BinVariations");

  TLegend *leg9 = (TLegend*) Bin_258_NEvent_Spread->GetListOfPrimitives()->At(5)->Clone();
  std::cout <<"h1" << std::endl;
  TList *p9 = leg9->GetListOfPrimitives();
  TIter next9(p9);
  TObject *obj9;
  TLegendEntry *le9;
  std::cout <<"h1" << std::endl;
  i=0;
  leg9->SetX1NDC(0.58);
  leg9->SetX2NDC(0.9);
  leg9->SetY1NDC(0.58);
  leg9->SetY2NDC(0.9);
  while ((obj9 = next9())) {
    le9 = (TLegendEntry*)obj9;
    i++;
    if (i==1) le9->SetLabel("Bin Content from Throws");
    if (i==2) le9->SetLabel("Gauss without MC Stats");
    if (i==3) le9->SetLabel("Gauss with MC Stats");
    if (i==4) le9->SetLabel("Nominal Bin Content");
  }
  std::cout <<"h1" << std::endl;
  Bin_258_NEvent_Spread->Draw();
  Bin_258_NEvent_Spread->GetListOfPrimitives()->At(5)->Delete();
  leg9->SetTextSize(0.035);
  leg9->Draw();
  Bin_258_NEvent_Spread->Print((std::string("detbin_nopisi258")+std::string(".pdf")).c_str());
  
  TLegend *leg10 = (TLegend*) Bin_97_NEvent_Spread->GetListOfPrimitives()->At(5)->Clone();
  TList *p10 = leg10->GetListOfPrimitives();
  TIter next10(p10);
  TObject *obj10;
  TLegendEntry *le10;
  i=0;
  leg10->SetX1NDC(0.58);
  leg10->SetX2NDC(0.9);
  leg10->SetY1NDC(0.58);
  leg10->SetY2NDC(0.9);
  while ((obj10 = next10())) {
    le10 = (TLegendEntry*)obj10;
    i++;
    if (i==1) le10->SetLabel("Bin Content from Throws");
    if (i==2) le10->SetLabel("Gauss without MC Stats");
    if (i==3) le10->SetLabel("Gauss with MC Stats");
    if (i==4) le10->SetLabel("Nominal Bin Content");
  }
  Bin_97_NEvent_Spread->Draw();
  Bin_97_NEvent_Spread->GetListOfPrimitives()->At(5)->Delete();
  leg10->SetTextSize(0.035);
  leg10->Draw();
  Bin_97_NEvent_Spread->Print((std::string("detbin_nopisi97")+std::string(".pdf")).c_str());
 
  TLegend *leg11 = (TLegend*) Bin_103_NEvent_Spread->GetListOfPrimitives()->At(5)->Clone();
  TList *p11 = leg11->GetListOfPrimitives();
  TIter next11(p11);
  TObject *obj11;
  TLegendEntry *le11;
  i=0;
  leg11->SetX1NDC(0.58);
  leg11->SetX2NDC(0.9);
  leg11->SetY1NDC(0.58);
  leg11->SetY2NDC(0.9);
  while ((obj11 = next11())) {
    le11 = (TLegendEntry*)obj11;
    i++;
    if (i==1) le11->SetLabel("Bin Content from Throws");
    if (i==2) le11->SetLabel("Gauss without MC Stats");
    if (i==3) le11->SetLabel("Gauss with MC Stats");
    if (i==4) le11->SetLabel("Nominal Bin Content");
  }
  Bin_103_NEvent_Spread->Draw();
  Bin_103_NEvent_Spread->GetListOfPrimitives()->At(5)->Delete();
  leg11->SetTextSize(0.035);
  leg11->Draw();
  Bin_103_NEvent_Spread->Print((std::string("detbin_nopisi103")+std::string(".pdf")).c_str());

  TLegend *leg12 = (TLegend*) Bin_570_NEvent_Spread->GetListOfPrimitives()->At(5)->Clone();
  TList *p12 = leg12->GetListOfPrimitives();
  TIter next12(p12);
  TObject *obj12;
  TLegendEntry *le12;
  i=0;
  leg12->SetX1NDC(0.58);
  leg12->SetX2NDC(0.9);
  leg12->SetY1NDC(0.58);
  leg12->SetY2NDC(0.9);
  while ((obj12 = next12())) {
    le12 = (TLegendEntry*)obj12;
    i++;
    if (i==1) le12->SetLabel("Bin Content from Throws");
    if (i==2) le12->SetLabel("Gauss without MC Stats");
    if (i==3) le12->SetLabel("Gauss with MC Stats");
    if (i==4) le12->SetLabel("Nominal Bin Content");
  }
  Bin_570_NEvent_Spread->Draw();
  Bin_570_NEvent_Spread->GetListOfPrimitives()->At(5)->Delete();
  leg12->SetTextSize(0.035);
  leg12->Draw();
  Bin_570_NEvent_Spread->Print((std::string("detbin_nopisi570")+std::string(".pdf")).c_str());
}

