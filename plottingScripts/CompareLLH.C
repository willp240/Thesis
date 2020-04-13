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

TCanvas *c = new TCanvas("canv", "canv", 1080, 1080);
std::string globfile = "";
std::string globfile2 = "";
TFile *global = NULL;
TFile *global2 = NULL;
std::string OldName;

void LoopDirectory(TDirectory *dir) {

  TIter next(dir->GetListOfKeys());
  TKey *key;
  std::string name;
  std::string param;
  // Loop through all entries
  while ((key = (TKey*)next())) {
    // Get name of object
    name = std::string(key->GetName());
    if (name.find(OldName) != std::string::npos) continue;
    std::string ClassName = key->GetClassName();
    if (ClassName == "TDirectoryFile") {
      dir->cd(name.c_str());
      TDirectory *SubDir = gDirectory;
      LoopDirectory(SubDir);
    }

    if (name.find("sam") == std::string::npos) continue;
    if (name.find("ndd") != std::string::npos) continue;
    //    if (name.find("b_") != std::string::npos) continue;

    if((name.find("b_0") == std::string::npos)&&(name.find("b_18")== std::string::npos)&&((name.find("b_40")== std::string::npos))&&(name.find("2p2h_shape_O") == std::string::npos)&&(name.find("FEFCX") == std::string::npos)&&(name.find("FEFINEL") == std::string::npos) )continue;
    param=name;
    name = std::string(dir->GetName())+"/"+name;
    if (name.find(".root") != std::string::npos) continue;
    // The stack and data
    TH1D *Plot1 = (TH1D*)(global->Get(name.c_str()));
    TH1D *Plot2 = (TH1D*)(global2->Get(name.c_str()));
    if (Plot1 == NULL || Plot2 == NULL) continue;
    if (Plot1->Integral() == 0 && Plot2->Integral() == 0) {
      delete Plot1;
      delete Plot2;
      continue;
    }

    Plot1->SetLineColor(kRed);
    Plot1->SetLineWidth(2);
    Plot2->SetLineColor(kBlue);
    Plot2->SetLineWidth(2);
    Plot2->SetLineStyle(kDashed);

    TLegend *leg = new TLegend(0.35, 0.55, 0.85, 0.85);
    leg->AddEntry(Plot1, "Uniform Rectangular", "l");
    leg->AddEntry(Plot2, "Non-Uniform Rectangular", "l");

    TPad *p1 = new TPad("p1", "p1", 0., 0.3, 1.0, 1.0);
    TPad *p2 = new TPad("p2", "p2", 0., 0.0, 1.0, 0.3);
    p1->SetTopMargin(c->GetTopMargin());
    p1->SetBottomMargin(0);
    p1->SetLeftMargin(c->GetLeftMargin());
    p1->SetRightMargin(c->GetRightMargin());

    p2->SetTopMargin(0);
    p2->SetBottomMargin(0.25);
    p2->SetLeftMargin(c->GetLeftMargin());
    p2->SetRightMargin(c->GetRightMargin());


    TH1D *ratio = (TH1D*)Plot1->Clone();
    ratio->Divide(Plot2);
    //    if (ratio->GetMaximum() < 1.05 && ratio->GetMinimum() > 0.95) continue;
    Plot1->GetYaxis()->SetRangeUser(0-Plot1->GetMaximum()*0.02, Plot1->GetMaximum()*1.02);

    c->cd();
    p1->Draw();
    p1->cd();
    Plot1->GetYaxis()->SetTitle("-2LLH_{Sample}");
    Plot1->Draw("c");
    Plot2->Draw("c, same");
    leg->Draw("same");
    Plot1->GetYaxis()->SetLabelSize(0.);
    Plot1->GetYaxis()->SetTitleSize(30);
    Plot1->GetYaxis()->SetTitleFont(43);
    Plot1->GetYaxis()->SetTitleOffset(Plot1->GetYaxis()->GetTitleOffset()*1.7);
    Plot1->SetTitle("");
    c->Update();

    TGaxis *axis = new TGaxis(Plot1->GetXaxis()->GetBinLowEdge(Plot1->GetXaxis()->GetFirst()), gPad->GetUymin()-gPad->GetUymax()*0.02,
        Plot1->GetXaxis()->GetBinLowEdge(Plot1->GetXaxis()->GetFirst()), gPad->GetUymax(),
        gPad->GetUymin()-gPad->GetUymax()*0.02, gPad->GetUymax(),
        510, "");
    axis->SetLabelFont(43);
    axis->SetLabelSize(25);
    axis->Draw();

    c->cd();
    p2->Draw();
    p2->cd();

    ratio->SetLineWidth(2);
    ratio->SetTitle("");

    ratio->SetMarkerSize(2);
    ratio->SetMarkerStyle(20);

    ratio->GetYaxis()->SetTitleSize(30);
    ratio->GetYaxis()->SetTitleFont(43);
    ratio->GetYaxis()->SetTitleOffset(1.7);
    ratio->GetYaxis()->SetLabelFont(43);
    ratio->GetYaxis()->SetLabelSize(25);
    ratio->GetYaxis()->CenterTitle();
    ratio->GetYaxis()->SetNdivisions(5,3,0);
    ratio->GetXaxis()->SetTitle("Variation rel. nom.");

    ratio->GetXaxis()->SetTitleSize(30);
    ratio->GetXaxis()->SetTitleFont(43);
    ratio->GetXaxis()->SetTitleOffset(4.0);
    ratio->GetXaxis()->SetNdivisions(7,5,0);
    Plot1->GetXaxis()->SetNdivisions(7,5,0);
    ratio->GetXaxis()->SetLabelFont(43);
    ratio->GetXaxis()->SetLabelSize(25);
    ratio->GetYaxis()->SetTitle("Ratio");

    ratio->SetMinimum(-0.1);
    ratio->SetMaximum(2.1);
    ratio->Draw("c");

    c->Print((param+std::string("polyLLH.pdf")).c_str());

    c->Clear();
    delete Plot1;
    delete Plot2;
    delete ratio;
    delete leg;
  }

}

void CompareLLH(std::string file1, std::string file2) {

  TFile *rootfile = TFile::Open(file1.c_str());
  TFile *rootfile2 = TFile::Open(file2.c_str());

  std::cout << "2017 Runs 2-6" << file1 << std::endl;
  std::cout << "2020 Runs 2-9" << file2 << std::endl;

  global = rootfile;
  global2 = rootfile2;

  globfile = file1;
  globfile2 = file2;

  c->SetTopMargin(0.07);
  c->SetBottomMargin(0.08);
  c->SetLeftMargin(0.10);
  c->SetRightMargin(0.03);
  OldName = "empty";
  LoopDirectory(rootfile);
  
}

