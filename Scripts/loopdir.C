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
std::string globfile3 = "";
TFile *global = NULL;
TFile *global2 = NULL;
TFile *global3 = NULL;
std::string OldName;

void LoopDirectory(TDirectory *dir) {

  TIter next(dir->GetListOfKeys());
  TKey *key;
  std::string name;

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


    name = std::string(dir->GetName())+"/"+name;
    if (name.find(".root") != std::string::npos) continue;
    if (name.find("ndd") != std::string::npos) continue;
    // if (global2->GetListOfKeys()->Contains(name.c_str())) std::cout << "glob2 contains " << name <<std::endl;// || !global->GetListOfKeys()->Contains(name.c_str())) continue;
    std::cout << name << std::endl;

    // The stack and data
    TH1D *Plot1 = (TH1D*)(global->Get(name.c_str()));
    TH1D *Plot2 = (TH1D*)(global2->Get(name.c_str()));
    if(!Plot1 || !Plot2)continue;

    if (Plot1->Integral() == 0 && Plot2->Integral() == 0) {
      delete Plot1;
      delete Plot2;
      continue;
    }
    Plot1->SetLineColor(kRed);
    Plot1->SetLineWidth(2);
    Plot2->SetLineColor(kBlue);
    Plot2->SetLineWidth(2);
    Plot2->SetLineStyle(2);

    TLegend *leg = new TLegend(0.25, 0.25, 0.75, 0.75);
    leg->AddEntry(Plot1, "2017", "l");
    leg->AddEntry(Plot2, "2020", "l");

    c->Clear();
    Plot1->Draw();
    Plot2->Draw("same");
    leg->Draw("same");
    c->Print((globfile+std::string("_")+globfile2+std::string(".pdf")).c_str());

    delete Plot1;
    delete Plot2;
    delete leg;
  }

}

void loopdir(std::string file1, std::string file2) {

  TFile *rootfile = TFile::Open(file1.c_str());
  TFile *rootfile2 = TFile::Open(file2.c_str());

  std::cout << "2017" << file1 << std::endl;
  std::cout << "2020" << file2 << std::endl;

  global = rootfile;
  global2 = rootfile2;

  globfile = file1;
  globfile2 = file2;

  c->SetTopMargin(0.07);
  c->SetBottomMargin(0.08);
  c->SetLeftMargin(0.10);
  c->SetRightMargin(0.03);
  c->Print((globfile+std::string("_")+globfile2+std::string(".pdf[")).c_str());

  OldName = "empty";
  LoopDirectory(rootfile);

  c->Print((globfile+std::string("_")+globfile2+std::string(".pdf]")).c_str());
}
