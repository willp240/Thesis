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

    if (name.find("full") == std::string::npos) continue;
    if (name.find("ndd") != std::string::npos) continue;
    if (name.find("b_") != std::string::npos) continue;

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
    //std::cout << "plot1 xmin " << Plot1->GetXaxis()->GetXmin() << " plot1 xmax " << Plot1->GetXaxis()->GetXmax() << " plot2 xmin " << Plot2->GetXaxis()->GetXmin() << " plot2 xmax " << Plot2->GetXaxis()->GetXmax() << std::endl;

    //    Plot1->GetXaxis()->SetRangeUser(Plot2->GetXaxis()->GetXmin(),Plot2->GetXaxis()->GetXmax());

    //std::cout << "plot1 xmin " << Plot1->GetXaxis()->GetXmin() << " plot1 xmax " << Plot1->GetXaxis()->GetXmax() << " plot2 xmin " << Plot2->GetXaxis()->GetXmin() << " plot2 xmax " << Plot2->GetXaxis()->GetXmax() << std::endl;

    TLegend *leg = new TLegend(0.35, 0.55, 0.85, 0.85);
    leg->AddEntry(Plot1, "Runs 2-6 2017", "l");
    leg->AddEntry(Plot2, "Runs 2-9 2019", "l");

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

    //    double minx, maxx, p, r;
    //int nbins =150, pbin_plot1, pbin_plot2;

    /*    std::cout << "Plot1 xmin " << Plot1->GetXaxis()->GetXmin() << " Plot1 xmax " << Plot1->GetXaxis()->GetXmax() << " Plot2 xmin " << Plot2->GetXaxis()->GetXmin() << " Plot2 xmax " << Plot2->GetXaxis()->GetXmax() << std::endl;

    if(Plot1->GetXaxis()->GetXmin()>Plot2->GetXaxis()->GetXmin())
      minx = Plot1->GetXaxis()->GetXmin();
    else
      minx = Plot2->GetXaxis()->GetXmin();

    if(Plot1->GetXaxis()->GetXmax()<Plot2->GetXaxis()->GetXmax())
      maxx = Plot1->GetXaxis()->GetXmax();
    else
      maxx = Plot2->GetXaxis()->GetXmax();

    std::cout << "xmin " << minx << " xmax " << maxx << std::endl;

    //    TH1D *ratio = new TH1D("ratio","ratio",nbins,minx,maxx);
    */
    TH1D *ratio = (TH1D*)Plot1->Clone("ratio");
    /*
    for (int pbin=1;pbin<=nbins;pbin++){

      p=ratio->GetXaxis()->GetBinCenter(pbin);
      pbin_plot1 = Plot1->GetXaxis()->FindBin(p);
      pbin_plot2 = Plot2->GetXaxis()->FindBin(p);
      r = Plot2->GetBinContent(pbin_plot2)/Plot1->GetBinContent(pbin_plot1);
      ratio->SetBinContent(pbin,r);
      std::cout << "For bin " << pbin << ", p=" << p << std::endl;
      std::cout << "Plot1 bin " << pbin_plot1 << " Plot2 bin " << pbin_plot2 << std::endl;
      std::cout << "Plot1 content "<< Plot1->GetBinContent(pbin_plot1) << " Plot2 content " << Plot2->GetBinContent(pbin_plot2) <<std::endl;
      std::cout << "Ratio " << r << std::endl << std::endl;
    }
    */
    ratio->Divide(Plot2);
    //    if (ratio->GetMaximum() < 1.05 && ratio->GetMinimum() > 0.95) continue;
    Plot1->GetYaxis()->SetRangeUser(0-Plot1->GetMaximum()*0.02, Plot1->GetMaximum()*1.02);
    
    c->cd();
    p1->Draw();
    p1->cd();
    Plot1->GetYaxis()->SetTitle("#chi^{2}");
    Plot1->Draw("c");
    Plot2->Draw("c, same");
    leg->Draw("same");
    Plot1->GetYaxis()->SetLabelSize(0.);
    Plot1->GetYaxis()->SetTitleSize(30);
    Plot1->GetYaxis()->SetTitleFont(43);
    Plot1->GetYaxis()->SetTitleOffset(Plot1->GetYaxis()->GetTitleOffset()*1.7);
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
    ratio->GetYaxis()->SetTitle("Upd/Old");

    ratio->SetMinimum(-0.1);
    ratio->SetMaximum(2.1);
    ratio->Draw("c");

    c->Print((globfile+std::string("_")+globfile2+std::string(".pdf")).c_str());

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

  std::cout << "Runs2-6 " << file1 << std::endl;
  std::cout << "Runs2-9 " << file2 << std::endl;

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

