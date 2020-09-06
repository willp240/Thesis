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

  TLegend *leg = (TLegend*) Bin_5_NEvent_Spread->GetListOfPrimitives()->At(5)->Clone();
  TList *p = leg->GetListOfPrimitives();
  TIter next(p);
  TObject *obj;
  TLegendEntry *le;
  int i=0;
  leg->SetX1NDC(0.58);
  leg->SetX2NDC(0.9);
  leg->SetY1NDC(0.58);
  leg->SetY2NDC(0.9);
  while ((obj = next())) {
    le = (TLegendEntry*)obj;
    i++;
    if (i==1) le->SetLabel("Bin Content from Throws");
    if (i==2) le->SetLabel("Gauss without MC Stats");
    if (i==3) le->SetLabel("Gauss with MC Stats");
    if (i==4) le->SetLabel("Nominal Bin Content");
  }
  Bin_5_NEvent_Spread->Draw();
  Bin_5_NEvent_Spread->GetListOfPrimitives()->At(5)->Delete();
  leg->SetTextSize(0.035);
  leg->Draw();
  Bin_5_NEvent_Spread->Print((std::string("detbin_allsysts5")+std::string(".pdf")).c_str());

  TLegend *leg2 = (TLegend*) Bin_253_NEvent_Spread->GetListOfPrimitives()->At(5)->Clone();
  TList *p2 = leg2->GetListOfPrimitives();
  TIter next2(p2);
  TObject *obj2;
  TLegendEntry *le2;
  i=0;
  leg2->SetX1NDC(0.58);
  leg2->SetX2NDC(0.9);
  leg2->SetY1NDC(0.58);
  leg2->SetY2NDC(0.9);
  while ((obj2 = next2())) {
    le2 = (TLegendEntry*)obj2;
    i++;
    if (i==1) le2->SetLabel("Bin Content from Throws");
    if (i==2) le2->SetLabel("Gauss without MC Stats");
    if (i==3) le2->SetLabel("Gauss with MC Stats");
    if (i==4) le2->SetLabel("Nominal Bin Content");
  }
  Bin_253_NEvent_Spread->Draw();
  Bin_253_NEvent_Spread->GetListOfPrimitives()->At(5)->Delete();
  leg2->SetTextSize(0.035);
  leg2->Draw();
  Bin_253_NEvent_Spread->Print((std::string("detbin_allsysts253")+std::string(".pdf")).c_str());
  
  TLegend *leg3 = (TLegend*) Bin_430_NEvent_Spread->GetListOfPrimitives()->At(5)->Clone();
  TList *p3 = leg3->GetListOfPrimitives();
  TIter next3(p3);
  TObject *obj3;
  TLegendEntry *le3;
  i=0;
  leg3->SetX1NDC(0.58);
  leg3->SetX2NDC(0.9);
  leg3->SetY1NDC(0.58);
  leg3->SetY2NDC(0.9);
  while ((obj3 = next3())) {
    le3 = (TLegendEntry*)obj3;
    i++;
    if (i==1) le3->SetLabel("Bin Content from Throws");
    if (i==2) le3->SetLabel("Gauss without MC Stats");
    if (i==3) le3->SetLabel("Gauss with MC Stats");
    if (i==4) le3->SetLabel("Nominal Bin Content");
  }
  Bin_430_NEvent_Spread->Draw();
  Bin_430_NEvent_Spread->GetListOfPrimitives()->At(5)->Delete();
  leg3->SetTextSize(0.035);
  leg3->Draw();
  Bin_430_NEvent_Spread->Print((std::string("detbin_allsysts430")+std::string(".pdf")).c_str());
  
  TLegend *leg4 = (TLegend*) Bin_535_NEvent_Spread->GetListOfPrimitives()->At(5)->Clone();
  TList *p4 = leg4->GetListOfPrimitives();
  TIter next4(p4);
  TObject *obj4;
  TLegendEntry *le4;
  i=0;
  leg4->SetX1NDC(0.58);
  leg4->SetX2NDC(0.9);
  leg4->SetY1NDC(0.58);
  leg4->SetY2NDC(0.9);
  while ((obj4 = next4())) {
    le4 = (TLegendEntry*)obj4;
    i++;
    if (i==1) le4->SetLabel("Bin Content from Throws");
    if (i==2) le4->SetLabel("Gauss without MC Stats");
    if (i==3) le4->SetLabel("Gauss with MC Stats");
    if (i==4) le4->SetLabel("Nominal Bin Content");
  }
  Bin_535_NEvent_Spread->Draw();
  Bin_535_NEvent_Spread->GetListOfPrimitives()->At(5)->Delete();
  leg4->SetTextSize(0.035);
  leg4->Draw();
  Bin_535_NEvent_Spread->Print((std::string("detbin_allsysts535")+std::string(".pdf")).c_str());
  ////////////////////////////////////

  TLegend *leg5 = (TLegend*) Bin_258_NEvent_Spread->GetListOfPrimitives()->At(5)->Clone();
  TList *p5 = leg5->GetListOfPrimitives();
  TIter next5(p5);
  TObject *obj5;
  TLegendEntry *le5;
  i=0;
  leg5->SetX1NDC(0.58);
  leg5->SetX2NDC(0.9);
  leg5->SetY1NDC(0.58);
  leg5->SetY2NDC(0.9);
  while ((obj5 = next5())) {
    le5 = (TLegendEntry*)obj5;
    i++;
    if (i==1) le5->SetLabel("Bin Content from Throws");
    if (i==2) le5->SetLabel("Gauss without MC Stats");
    if (i==3) le5->SetLabel("Gauss with MC Stats");
    if (i==4) le5->SetLabel("Nominal Bin Content");
  }
  Bin_258_NEvent_Spread->Draw();
  Bin_258_NEvent_Spread->GetListOfPrimitives()->At(5)->Delete();
  leg5->SetTextSize(0.035);
  leg5->Draw();
  Bin_258_NEvent_Spread->Print((std::string("detbin_allsysts258")+std::string(".pdf")).c_str());

  TLegend *leg6 = (TLegend*) Bin_97_NEvent_Spread->GetListOfPrimitives()->At(5)->Clone();
  TList *p6 = leg6->GetListOfPrimitives();
  TIter next6(p6);
  TObject *obj6;
  TLegendEntry *le6;
  i=0;
  leg6->SetX1NDC(0.58);
  leg6->SetX2NDC(0.9);
  leg6->SetY1NDC(0.58);
  leg6->SetY2NDC(0.9);
  while ((obj6 = next6())) {
    le6 = (TLegendEntry*)obj6;
    i++;
    if (i==1) le6->SetLabel("Bin Content from Throws");
    if (i==2) le6->SetLabel("Gauss without MC Stats");
    if (i==3) le6->SetLabel("Gauss with MC Stats");
    if (i==4) le6->SetLabel("Nominal Bin Content");
  }
  Bin_97_NEvent_Spread->Draw();
  Bin_97_NEvent_Spread->GetListOfPrimitives()->At(5)->Delete();
  leg6->SetTextSize(0.035);
  leg6->Draw();
  Bin_97_NEvent_Spread->Print((std::string("detbin_allsysts97")+std::string(".pdf")).c_str());
  
  TLegend *leg7 = (TLegend*) Bin_103_NEvent_Spread->GetListOfPrimitives()->At(5)->Clone();
  TList *p7 = leg7->GetListOfPrimitives();
  TIter next7(p7);
  TObject *obj7;
  TLegendEntry *le7;
  i=0;
  leg7->SetX1NDC(0.58);
  leg7->SetX2NDC(0.9);
  leg7->SetY1NDC(0.58);
  leg7->SetY2NDC(0.9);
  while ((obj7 = next7())) {
    le7 = (TLegendEntry*)obj7;
    i++;
    if (i==1) le7->SetLabel("Bin Content from Throws");
    if (i==2) le7->SetLabel("Gauss without MC Stats");
    if (i==3) le7->SetLabel("Gauss with MC Stats");
    if (i==4) le7->SetLabel("Nominal Bin Content");
  }
  Bin_103_NEvent_Spread->Draw();
  Bin_103_NEvent_Spread->GetListOfPrimitives()->At(5)->Delete();
  leg7->SetTextSize(0.035);
  leg7->Draw();
  Bin_103_NEvent_Spread->Print((std::string("detbin_allsysts103")+std::string(".pdf")).c_str());
  
  TLegend *leg8 = (TLegend*) Bin_570_NEvent_Spread->GetListOfPrimitives()->At(5)->Clone();
  TList *p8 = leg8->GetListOfPrimitives();
  TIter next8(p8);
  TObject *obj8;
  TLegendEntry *le8;
  i=0;
  leg8->SetX1NDC(0.58);
  leg8->SetX2NDC(0.9);
  leg8->SetY1NDC(0.58);
  leg8->SetY2NDC(0.9);
  while ((obj8 = next8())) {
    le8 = (TLegendEntry*)obj8;
    i++;
    if (i==1) le8->SetLabel("Bin Content from Throws");
    if (i==2) le8->SetLabel("Gauss without MC Stats");
    if (i==3) le8->SetLabel("Gauss with MC Stats");
    if (i==4) le8->SetLabel("Nominal Bin Content");
  }
  Bin_570_NEvent_Spread->Draw();
  Bin_570_NEvent_Spread->GetListOfPrimitives()->At(5)->Delete();
  leg8->SetTextSize(0.035);
  leg8->Draw();
  Bin_570_NEvent_Spread->Print((std::string("detbin_allsysts570")+std::string(".pdf")).c_str());
  /////////////////////////
}
/* 
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
*/
