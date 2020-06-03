#include "TH2D.h"
#include "TH1D.h"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TStyle.h"
#include "TAxis.h"
#include "TGaxis.h"
#include "TGraph.h"
#include "TF1.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TMath.h"
#include "TVectorD.h"
#include <iostream>

void getInterval1D(TH1D *hist, double &p68, double &p90, double &p95, double &p99, double &p3sig)
{
  std::cout << "getting interval" << std::endl;
  TH1D *hCopy = (TH1D*)hist->Clone("hCopy");

  double integral = hCopy->Integral();
  double tsum = 0;
  
  //NB THIS IS NOT THE SET OF LEVELS USED IF YOU JUST CALL THIS ON THE COMMAND LINE!! LOOK LOWER FOR THAT
  double contlevel1=0.68;
  double contlevel2=0.9;
  double contlevel3=0.955;
  double contlevel4=0.99;
  double contlevel5=0.9973;

  double cont68lev = 0;
  double cont90lev = 0;
  double cont95lev = 0;
  double cont99lev = 0;
  double cont3siglev = 0;

  std::cout << integral << std::endl;

  while((tsum / integral) < contlevel5)
    {
      double tmax = hCopy->GetMaximum();
      tsum = tsum + tmax;
      int bin = hCopy->GetMaximumBin();
      if (tsum / integral < contlevel1)
	{
	  cont68lev = tmax;
	  hCopy->SetBinContent(bin, -1.0);
	}
      if ((tsum / integral < contlevel2) && (tsum / integral > contlevel1))
	{
	  cont90lev = tmax;
	  hCopy->SetBinContent(bin, -3.0);
	}
      if ((tsum / integral < contlevel3) && (tsum / integral > contlevel2))
	{
	  cont95lev = tmax;
	  hCopy->SetBinContent(bin, -5.0);
	}
      if ((tsum / integral < contlevel4) && (tsum / integral > contlevel3))
	{
	  cont99lev = tmax;
	  hCopy->SetBinContent(bin, -7.0);
	}
      if ((tsum / integral < contlevel5) && (tsum / integral > contlevel4))
	{
	  cont3siglev = tmax;
	  hCopy->SetBinContent(bin, -9.0);
	}
    }

  p3sig = cont3siglev;
  p99 = cont99lev;
  p95 = cont95lev;
  p90 = cont90lev;
  p68 = cont68lev;

  std::cout << "p68 = " << p68 << ", p90 = " << p90 << ", p95 = " << p95 << std::endl;
  std::cout << "p99 = " << p68 << ", p3sig = " << p3sig << std::endl;
}

void getInterval1D(TH1D *hist, TH1D &h68, TH1D &h90, TH1D &h95, TH1D &h99, TH1D &h3sig)
{
  std::cout << "getting interval" << std::endl;
  TH1D *hCopy = (TH1D*)hist->Clone("hCopy");
  TH1D *hCopy68 = (TH1D*)hist->Clone("hCopy68");
  TH1D *hCopy90 = (TH1D*)hist->Clone("hCopy90");
  TH1D *hCopy95 = (TH1D*)hist->Clone("hCopy95");
  TH1D *hCopy99 = (TH1D*)hist->Clone("hCopy99");
  TH1D *hCopy3sig = (TH1D*)hist->Clone("hCopy3sig");

  double integral = hCopy->Integral();
  double tsum = 0;
  double cont68lev = 0;
  double cont90lev = 0;
  double cont95lev = 0;
  double cont99lev = 0;
  double cont3siglev = 0;

  //!!NB THIS IS NOT THE SET OF LEVELS USED IF YOU JUST CALL THIS ON THE COMMAND LINE!! LOOK LOWER FOR THAT
  double contlevel1=0.68;
  double contlevel2=0.9;
  double contlevel3=0.955;
  double contlevel4=0.99;
  double contlevel5=0.9973;

  std::cout << integral << std::endl;

  while((tsum / integral) < contlevel5)
    {
      double tmax = hCopy->GetMaximum();
      tsum = tsum + tmax;
      int bin = hCopy->GetMaximumBin();
      if (tsum / integral < contlevel1)
	{
	  cont68lev = tmax;
	  hCopy->SetBinContent(bin, -1.0);
	  hCopy68->SetBinContent(bin, 0);
	  hCopy90->SetBinContent(bin, 0);
	  hCopy95->SetBinContent(bin, 0);
	  hCopy99->SetBinContent(bin, 0);
	  hCopy3sig->SetBinContent(bin, 0);
	}
      if ((tsum / integral < contlevel2) && (tsum / integral > contlevel1))
	{
	  cont90lev = tmax;
	  hCopy->SetBinContent(bin, -3.0);
	  hCopy90->SetBinContent(bin, 0);
	}
      if ((tsum / integral < contlevel3) && (tsum / integral > contlevel2))
	{
	  cont95lev = tmax;
	  hCopy->SetBinContent(bin, -5.0);
	  hCopy90->SetBinContent(bin, 0);
	}
    }

  h95 = (*hCopy95);
  h90 = (*hCopy90);
  h68 = (*hCopy68);
}

/// Set onehierarchyonly = 1 for NH, -1 for IH, 0 for all
/// Set axisX_sindcp to true if you want to fill with sin(dCP) instead of (dCP)
/// Set prior_sindcp to true if you want to apply flat prior in sin(dCP) (assume you used flat prior in dCP in the fit)
void contours_app1D(TString fs,int onehierarchyonly=0,int burnIn=0, bool doRCrw=false, bool doLogY = false, bool axisX_sindcp = false, bool prior_sindcp = false)
{
  //!!NB THIS IS THE SET OF LEVELS USED IF YOU JUST CALL THIS ON THE COMMAND LINE!!
  double contlevel1=0.68;
  double contlevel2=0.90;
  double contlevel3=0.954;
  double contlevel4=0.99;
  double contlevel5=0.9973;

  gStyle->SetPalette(51,0);
  gStyle->SetOptStat(0);
  TGaxis::SetMaxDigits(4);

  TFile *f = new TFile(fs.Data(), "READ");
  TTree *t = (TTree*)f->Get("osc_posteriors");

  TH1D *h_th13 = new TH1D("h_th13","",400,0,0.1);
  if(doRCrw) h_th13->SetBins(400,0.01,0.04);
  TH1D* h_dcp;
  if(axisX_sindcp) h_dcp = new TH1D("h_dcp","",/*100*/100,-1,1);
  else h_dcp = new TH1D("h_dcp","",/*100*/200,-1*TMath::Pi(),TMath::Pi());

  int max = t->GetEntries();

  double tth13, tth23, tdcp, tdm32, tbiaswgt = 1.0;
  double RCrw = 1.0;
  int stp;

  t->SetBranchAddress("dcp",&tdcp);
  t->SetBranchAddress("theta13",&tth13);
  t->SetBranchAddress("theta23",&tth23);
  t->SetBranchAddress("dm23", &tdm32);
  t->SetBranchAddress("biaswgt", &tbiaswgt);
  if(doRCrw=true) t->SetBranchAddress("RCreweight", &RCrw);
  t->SetBranchAddress("step", &stp);

  // Apply different prior to osc pars 
  // Note: this is assuming it is being run on the fit without reactor (so no prior to ‘undo’)
  
  // Flat in theta23 or theta13 (not sin^2 theta23 or sin^2 theta13) 
  TF1 *prior1 = new TF1("prior1","1.0/1.0/(2*TMath::Sqrt(x - x*x))");
  
  // Flat in sin^2 2theta13 (not sin^2 theta13) 
  TF1 *prior2 = new TF1("prior2", "4 - 8*x"); 
  
  // Flat in sin dcp (not dcp)
  TF1 *prior3 = new TF1("prior3","TMath::Abs(TMath::Cos(x))");
  
  // Flat in cos dcp (not dcp)
  TF1 *prior4 = new TF1("prior4","TMath::Abs(TMath::Sin(x))");
  
  for (int i=0; i<max; i++)
    {
      t->GetEntry(i);

      if(stp<burnIn) continue;

      if(onehierarchyonly==1){
	if (tdm32 < 0) continue; // NH only
      }
      if(onehierarchyonly==-1){
	if (tdm32 > 0) continue; // IH only
      }

      // Apply different prior to osc pars 
      // Note: this is assuming it is being run on the fit without reactor (so no prior to ‘undo’)
      // Comment out unless needed 

      // Flat in theta23 (not sin^2 theta23) 
      //tbiaswgt = prior1->Eval(tth23);

      // Flat in theta13 (not sin^2 theta13)
      //tbiaswgt = prior1->Eval(tth13);

      // Flat in sin^2 2theta13 (not sin^2 theta13) 
      //tbiaswgt = prior2->Eval(tth13);

      //Flat in sin dcp (not dcp)
      if(prior_sindcp)
	tbiaswgt = prior3->Eval(tdcp);  

      // Flat in cos dcp (not dcp)
      //tbiaswgt = prior4->Eval(tdcp);  
      double totWeight = tbiaswgt*RCrw;

      h_th13->Fill(tth13, totWeight);

      if(axisX_sindcp)
	h_dcp->Fill(TMath::Sin(tdcp), totWeight);
      else
	h_dcp->Fill(tdcp, totWeight);
    }

  if(doRCrw) h_th13->GetXaxis()->SetRangeUser(0.015,0.035);


  if(!doLogY) {
    //h_th13->Scale(1./h_th13->Integral());
    //h_dcp->Scale(1./h_dcp->Integral());
  }

  if(prior_sindcp) h_dcp->SetTitle("Prior: flat in sin#delta_{CP}");

  h_th13->GetXaxis()->SetTitle("sin^{2}#theta_{13}");
  //h_dcp->GetYaxis()->SetTitle("Postfit N_{steps}");
  h_th13-> GetXaxis() -> SetTitleOffset(1.2) ;
  h_th13-> GetXaxis() -> SetTitleFont(132) ;
  h_th13-> GetXaxis() -> SetTitleSize(0.06) ;
  h_th13-> GetXaxis() -> SetLabelFont(132) ;
  h_th13-> GetXaxis() -> SetLabelSize(0.05) ;
  h_th13->GetYaxis()->SetTitle("Posterior probability density");
  h_th13-> GetYaxis() -> SetTitleOffset(1.2) ;
  h_th13-> GetYaxis() -> SetTitleFont(132) ;
  h_th13-> GetYaxis() -> SetTitleSize(0.06) ;
  h_th13-> GetYaxis() -> SetLabelFont(132) ;
  h_th13-> GetYaxis() -> SetLabelSize(0.05) ;

  if(axisX_sindcp)
    h_dcp->GetXaxis()->SetTitle("sin(#delta_{CP})");
  else
    h_dcp->GetXaxis()->SetTitle("#delta_{CP} (rad.)");
  h_dcp-> GetXaxis() -> SetTitleOffset(1.2) ;
  h_dcp-> GetXaxis() -> SetTitleFont(132) ;
  h_dcp-> GetXaxis() -> SetTitleSize(0.06) ;
  h_dcp-> GetXaxis() -> SetLabelFont(132) ;
  h_dcp-> GetXaxis() -> SetLabelSize(0.05) ;
  h_dcp->GetYaxis()->SetTitle("Posterior probability density");
  h_dcp-> GetYaxis() -> SetTitleOffset(1.2) ;
  h_dcp-> GetYaxis() -> SetTitleFont(132) ;
  h_dcp-> GetYaxis() -> SetTitleSize(0.06) ;
  h_dcp-> GetYaxis() -> SetLabelFont(132) ;
  h_dcp-> GetYaxis() -> SetLabelSize(0.05) ;


  
  TH1D *h68th = (TH1D*)h_th13->Clone("h68th");
  TH1D *h90th = (TH1D*)h_th13->Clone("h90th");
  TH1D *h95th = (TH1D*)h_th13->Clone("h95th");
  TH1D *h99th = (TH1D*)h_th13->Clone("h99th");
  TH1D *h3sigth = (TH1D*)h_th13->Clone("h3sigth");

  TH1D *h68dcp = (TH1D*)h_dcp->Clone("h68dcp"); 
  TH1D *h90dcp = (TH1D*)h_dcp->Clone("h90dcp");
  TH1D *h95dcp = (TH1D*)h_dcp->Clone("h95dcp");
  TH1D *h99dcp = (TH1D*)h_dcp->Clone("h99dcp");
  TH1D *h3sigdcp = (TH1D*)h_dcp->Clone("h3sigdcp");

  TH1D *copyth = (TH1D*)h_th13->Clone("copyth");
  TH1D *copydcp = (TH1D*)h_dcp->Clone("copydcp");

  double integral, tsum=0;

  // Get intervals (th_13)
  std::cout << "getting th13 1D interval" << std::endl;

  integral = copyth->Integral();

  std::cout << integral << std::endl;

  while((tsum / integral) < contlevel5)
    {
      double tmax = copyth->GetMaximum();
      //tsum = tsum + tmax;
      int bin = copyth->GetMaximumBin();
      if (tsum / integral < contlevel1)
	{
	  copyth->SetBinContent(bin, -1.0);
	  h68th->SetBinContent(bin, 0);
	  h90th->SetBinContent(bin, 0);
	  h95th->SetBinContent(bin, 0);
	  h99th->SetBinContent(bin, 0);
	  h3sigth->SetBinContent(bin, 0);
	}
      if ((tsum / integral < contlevel2) && (tsum / integral > contlevel1))
	{
	  copyth->SetBinContent(bin, -3.0);
	  h90th->SetBinContent(bin, 0);
	  h95th->SetBinContent(bin, 0);
	  h99th->SetBinContent(bin, 0);
	  h3sigth->SetBinContent(bin, 0);
	}
      if ((tsum / integral < contlevel3) && (tsum / integral > contlevel2))
	{
	  copyth->SetBinContent(bin, -5.0);
	  h95th->SetBinContent(bin, 0);
	  h99th->SetBinContent(bin, 0);
	  h3sigth->SetBinContent(bin, 0);
	}
      if ((tsum / integral < contlevel4) && (tsum / integral > contlevel3))
	{
	  copyth->SetBinContent(bin, -7.0);
	  h99th->SetBinContent(bin, 0);
	  h3sigth->SetBinContent(bin, 0);
	}
      if ((tsum / integral < contlevel5) && (tsum / integral > contlevel4))
	{
	  copyth->SetBinContent(bin, -9.0);
	  h3sigth->SetBinContent(bin, 0);
	}
      tsum += tmax;
    }


  integral=0; 
  tsum=0;

  // Get intervals (dcp)
  std::cout << "getting dcp 1D interval" << std::endl;

  integral = copydcp->Integral();

  std::cout << "integral = " << integral << std::endl;

  while((tsum / integral) < contlevel5)
    {
      double tmax = copydcp->GetMaximum();
      //tsum = tsum + tmax;
      int bin = copydcp->GetMaximumBin();
      if (tsum / integral < contlevel1)
	{
	  copydcp->SetBinContent(bin, -1.0);
	  h68dcp->SetBinContent(bin, 0);
	  h90dcp->SetBinContent(bin, 0);
	  h95dcp->SetBinContent(bin, 0);
	  h99dcp->SetBinContent(bin, 0);
	  h3sigdcp->SetBinContent(bin, 0);
	}
      if ((tsum / integral < contlevel2) && (tsum / integral > contlevel1))
	{
	  copydcp->SetBinContent(bin, -3.0);
	  h90dcp->SetBinContent(bin, 0);
	  h95dcp->SetBinContent(bin, 0);
	  h99dcp->SetBinContent(bin, 0);
	  h3sigdcp->SetBinContent(bin, 0);
	}
      if ((tsum / integral < contlevel3) && (tsum / integral > contlevel2))
	{
	  copydcp->SetBinContent(bin, -5.0);
	  h95dcp->SetBinContent(bin, 0);
	  h99dcp->SetBinContent(bin, 0);
	  h3sigdcp->SetBinContent(bin, 0);
	}
      if ((tsum / integral < contlevel4) && (tsum / integral > contlevel3))
	{
	  copydcp->SetBinContent(bin, -7.0);
	  h99dcp->SetBinContent(bin, 0);
	  h3sigdcp->SetBinContent(bin, 0);
        }
      if ((tsum / integral < contlevel5) && (tsum / integral > contlevel4))
	{
	  copydcp->SetBinContent(bin, -9.0);
	  h3sigdcp->SetBinContent(bin, 0);
        }
      tsum += tmax;
    }


  TCanvas *c = new TCanvas("c","",1200,800);
  c->SetLeftMargin(0.15) ;
  c->SetBottomMargin(0.15) ;

  TCanvas *c1 = new TCanvas("c1","",1200,800);
  c1->SetLeftMargin(0.15) ;
  c1->SetBottomMargin(0.15) ;

  TCanvas *c_sig = new TCanvas("c_sig","",1200,800);
  c_sig->SetLeftMargin(0.15) ;
  c_sig->SetBottomMargin(0.15) ;

  TCanvas *c1_sig = new TCanvas("c1_sig","",1200,800);
  c1_sig->SetLeftMargin(0.15) ;
  c1_sig->SetBottomMargin(0.15) ;

  h_th13->SetLineColor(kBlack);
  h_th13->GetYaxis()->SetTitleOffset(1.3);
  h68th->SetLineColor(kBlack);
  h90th->SetLineColor(kBlack);
  h95th->SetLineColor(kBlack);
  h99th->SetLineColor(kBlack);
  h3sigth->SetLineColor(kBlack);

  //h68th->SetLineStyle(2);
  //h90th->SetLineStyle(3);

  // one plot will have 68, 90, 99...
  h68th->SetFillColor(13);
  h_th13->SetFillColor(12);
  h90th->SetFillColor(11);
  h99th->SetFillColor(10);
  // the other will have 68, 95, 99.73 (1,2,3 sigma)...
  h95th->SetFillColor(11);
  h3sigth->SetFillColor(10);

  h_dcp->SetLineColor(kBlack);
  h_dcp->GetYaxis()->SetTitleOffset(1.3);
  //if(!axisX_sindcp) h_dcp->GetYaxis()->SetRangeUser(0,0.05);
  //else h_dcp->GetYaxis()->SetRangeUser(0,2);

  h68dcp->SetLineColor(kBlack);
  h90dcp->SetLineColor(kBlack);
  h95dcp->SetLineColor(kBlack);
  h99dcp->SetLineColor(kBlack);
  h3sigdcp->SetLineColor(kBlack);

  //h68dcp->SetLineStyle(2);
  //h90dcp->SetLineStyle(3);

  // one plot will have 68, 90, 99...
  h_dcp->SetFillColor(13);
  h68dcp->SetFillColor(12);
  h90dcp->SetFillColor(11);
  h99dcp->SetFillColor(10);
  // the other will have 68, 95, 99.73 (1,2,3 sigma)...
  h95dcp->SetFillColor(11);
  h3sigdcp->SetFillColor(10);

  TString firstlegnum=TString::Format("%2.0f",contlevel1*100);
  TString secondlegnum=TString::Format("%2.0f",contlevel2*100);
  TString thirdlegnum=TString::Format("%2.2f",contlevel3*100);
  TString fourthlegnum=TString::Format("%2.0f",contlevel4*100);
  TString fifthlegnum=TString::Format("%2.2f",contlevel5*100);

  TLegend *leg = new TLegend(0.5,0.73,0.85,0.85);
  leg->SetFillColor(kWhite);
  leg->SetFillStyle(0);
  leg->SetLineColor(kWhite);
  leg->SetTextFont(132);
  leg->SetTextSize(0.05);

  leg->AddEntry(h90dcp,fourthlegnum+"% Credible Interval","f");
  leg->AddEntry(h68dcp,secondlegnum+"% Credible Interval","f");
  leg->AddEntry(h_dcp,firstlegnum+"% Credible Interval","f");

  TLegend *leg_sig = new TLegend(0.5,0.73,0.85,0.85);
  leg_sig->SetFillColor(kWhite);
  leg_sig->SetFillStyle(0);
  leg_sig->SetLineColor(kWhite);
  leg_sig->SetTextFont(132);
  leg_sig->SetTextSize(0.05);

  //KEVIN YOU ARE HERE!!!
  leg_sig->AddEntry(h90dcp,"3#sigma Credible Interval","f");
  leg_sig->AddEntry(h68dcp,"2#sigma redible Interval","f");
  leg_sig->AddEntry(h_dcp,"1#sigma Credible Interval","f");

  // Stop zero suppression of y axis
  h_th13->GetYaxis()->SetRange(h_th13->GetYaxis()->GetFirst(), h_th13->GetYaxis()->GetLast());
  //  h68th->GetYaxis()->SetRange(h68th->GetYaxis()->GetFirst(), h68th->GetYaxis()->GetLast());
  //  h90th->GetYaxis()->SetRange(h90th->GetYaxis()->GetFirst(), h90th->GetYaxis()->GetLast());
  h_dcp->GetYaxis()->SetRange(h_dcp->GetYaxis()->GetFirst(), h_dcp->GetYaxis()->GetLast());
  //  h68dcp->GetYaxis()->SetRange(h68dcp->GetYaxis()->GetFirst(), h68dcp->GetYaxis()->GetLast());
  //  h90dcp->GetYaxis()->SetRange(h90dcp->GetYaxis()->GetFirst(), h90dcp->GetYaxis()->GetLast());
  
  c->cd();
  //  gPad->SetLeftMargin(0.13);
  //  gPad->SetBottomMargin(0.13);
  gPad->SetTickx();
  gPad->SetTicky();
  h_th13->Draw();
  h68th->Draw("SAME");
  h90th->Draw("SAME");
  h99th->Draw("SAME");
  if(doLogY) {
    h_dcp->GetYaxis()->SetRangeUser(100.,h_dcp->GetMaximum()*5.);
    c->SetLogy(1);
  }
  leg->Draw("SAME");
  gPad->RedrawAxis();
  c->Update();
  gPad->Update();

  double h_dcp_tot = h_dcp->Integral();
  double h68dcp_tot = h68dcp->Integral();
  double h90dcp_tot = h90dcp->Integral();
  double h95dcp_tot = h95dcp->Integral();
  double h99dcp_tot = h99dcp->Integral();
  double h3sigdcp_tot = h3sigdcp->Integral();

  std::cout << "h_dcp->Integral() = " << h_dcp_tot << std::endl
	    << "h68dcp->Integral() = " << h68dcp_tot << "     " 
	    << (1.0-h68dcp_tot/h_dcp_tot)*100.0 << "%" << std::endl
	    << "h90dcp->Integral() = " << h90dcp_tot << "      "  
	    << (1.0-h90dcp_tot/h_dcp_tot)*100.0 << "%" << std::endl
	    << "h95dcp->Integral() = " << h95dcp_tot << "      "  
	    << (1.0-h95dcp_tot/h_dcp_tot)*100.0 << "%" << std::endl
	    << "h99dcp->Integral() = " << h99dcp_tot << "      "  
	    << (1.0-h99dcp_tot/h_dcp_tot)*100.0 << "%" << std::endl
	    << "h3sigdcp->Integral() = " << h3sigdcp_tot << "      "  
	    << (1.0-h3sigdcp_tot/h_dcp_tot)*100.0 << "%" << std::endl;

  c1->cd();
  gPad->SetTickx() ;
  gPad->SetTicky() ;
  //  gPad->SetLeftMargin(0.13);
  //  gPad->SetBottomMargin(0.13);
  h_dcp->GetXaxis()->LabelsOption("u");
  h_dcp->Draw();
  h68dcp->Draw("SAME");
  h90dcp->Draw("SAME");
  h99dcp->Draw("SAME");
  if(doLogY) {
    h_dcp->GetYaxis()->SetRangeUser(100.,h_dcp->GetMaximum()*5.);
    c1->SetLogy(1);
  }
   leg->Draw("SAME");
  gPad->RedrawAxis();
  c1->Update();
  gPad->Update();

  // draw 1, 2, 3, sigma contours on one plot:
  c_sig->cd();
  //  gPad->SetLeftMargin(0.13);
  //  gPad->SetBottomMargin(0.13);
  gPad->SetTickx() ;
  gPad->SetTicky() ;
  h_th13->Draw();
  h68th->Draw("SAME");
  h95th->Draw("SAME");
  h3sigth->Draw("SAME");
  if(doLogY) {
    h_th13->GetYaxis()->SetRangeUser(100.,h_dcp->GetMaximum()*5.);
    c_sig->SetLogy(1);
  }
  leg_sig->Draw("SAME");
  gPad->RedrawAxis();
  c_sig->Update();
  gPad->Update();

  c1_sig->cd();
  gPad->SetTickx() ;
  gPad->SetTicky() ;
  //  gPad->SetLeftMargin(0.13);
  //  gPad->SetBottomMargin(0.13);
  h_dcp->GetXaxis()->LabelsOption("u");
  h_dcp->Draw();
  h68dcp->Draw("SAME");
  h95dcp->Draw("SAME");
  h3sigdcp->Draw("SAME");
  if(doLogY) {
    h_dcp->GetYaxis()->SetRangeUser(100.,h_dcp->GetMaximum()*5.);
    c1_sig->SetLogy(1);
  }
  leg_sig->Draw("SAME");
  gPad->RedrawAxis();
  c1_sig->Update();
  gPad->Update();

  //TFile* outfile=new TFile("contours_app1D.root","recreate");
  //outfile->cd();
  //h_dcp->Write();
  //h_th13->Write();

  c->SaveAs("contours_1D_th13.pdf");
  c1->SaveAs("contours_1D_dcp.pdf");
  c->SaveAs("contours_1D_th13.root");
  c1->SaveAs("contours_1D_dcp.root");
  c_sig->SaveAs("contours_1D_th13_sigmas.pdf");
  c1_sig->SaveAs("contours_1D_dcp_sigmas.pdf");
  c_sig->SaveAs("contours_1D_th13_sigmas.root");
  c1_sig->SaveAs("contours_1D_dcp_sigmas.root");


  // Get limits
  int thbin_low=0, thbin_high=0, dcpbin_low=0, dcpbin_high=0;
  int thbin90_low=0, thbin90_high=0, dcpbin90_low=0, dcpbin90_high=0;
  int thbin95_low=0, thbin95_high=0, dcpbin95_low=0, dcpbin95_high=0;
  int thbin99_low=0, thbin99_high=0, dcpbin99_low=0, dcpbin99_high=0;
  int thbin3sig_low=0, thbin3sig_high=0, dcpbin3sig_low=0, dcpbin3sig_high=0;

  int nbinsth = h_th13->GetNbinsX();
  int nbinsdcp = h_dcp->GetNbinsX();
  int centre_th = h_th13->GetMaximumBin();
  int max_dcp = h_dcp->GetMaximumBin();
  int min_dcp = h_dcp->GetMinimumBin();

  //std::cout << centre_th << "  " << max_dcp << std::endl;

  for (int i=centre_th; i>0; i--)
    {
      //std::cout << i << "  " << h68th->GetBinContent(i) << std::endl;
      if (h68th->GetBinContent(i)!=0) 
	{
	  thbin_low = i;
	  break;
	}
    }
  for (int i=centre_th; i>0; i--)
    {
      //std::cout << i << "  " << h90th->GetBinContent(i) << std::endl;
      if (h90th->GetBinContent(i)!=0) 
	{
	  thbin90_low = i;
	  break;
	}
    }
  for (int i=centre_th; i>0; i--)
    {
      if (h95th->GetBinContent(i)!=0) 
	{
	  thbin95_low = i;
	  break;
	}
    }
  for (int i=centre_th; i>0; i--)
    {
      if (h99th->GetBinContent(i)!=0) 
	{
	  thbin99_low = i;
	  break;
	}
    }
  for (int i=centre_th; i>0; i--)
    {
      if (h3sigth->GetBinContent(i)!=0) 
	{
	  thbin3sig_low = i;
	  break;
	}
    }
  

  for (int i=centre_th; i<nbinsth+1; i++)
    {
      //std::cout << i << "  " << h68th->GetBinContent(i) << std::endl;
      if (h68th->GetBinContent(i)!=0) 
	{
	  thbin_high = i;
	  break;
	}
    }
  for (int i=centre_th; i<nbinsth+1; i++)
    {
      //std::cout << i << "  " << h90th->GetBinContent(i) << std::endl;
      if (h90th->GetBinContent(i)!=0) 
	{
	  thbin90_high = i;
	  break;
	}
    }
  for (int i=centre_th; i<nbinsth+1; i++)
    {
      if (h95th->GetBinContent(i)!=0) 
	{
	  thbin95_high = i;
	  break;
	}
    }
  for (int i=centre_th; i<nbinsth+1; i++)
    {
      if (h99th->GetBinContent(i)!=0) 
	{
	  thbin99_high = i;
	  break;
	}
    }
  for (int i=centre_th; i<nbinsth+1; i++)
    {
      if (h3sigth->GetBinContent(i)!=0) 
	{
	  thbin3sig_high = i;
	  break;
	}
    }



  // dcp: credible regions may be disconnected                             
  // Work down from maximum to minimum, wrapping dcp if necessary  

  // Down from maximum (towards -ve dcp)

  // h68dcp contains the beige histogram -- first non-zero bins out from the maximum are the 68% limits
  for (int i=max_dcp; i<nbinsdcp+1; i--)
    {
      //std::cout << i << "  " << h68dcp->GetBinContent(i) << std::endl;

      // wrap if necessary
      int bin = i;
      if (i<=0)
	{
	  bin = i+nbinsdcp+1;
	  //std::cout << bin << "  " << h68dcp->GetBinCenter(i) << "   " << h68dcp->GetBinCenter(i)+TMath::Pi() << "   " << h68dcp->GetBinCenter(bin) << std::endl;
	}
      if (h68dcp->GetBinContent(bin)!=0) 
	{
	  dcpbin_low = bin;
	  //std::cout << "dcpbin_low: " << h68dcp->GetXaxis()->GetBinLowEdge(bin) << " -- " << h68dcp->GetXaxis()->GetBinUpEdge(bin) << std::endl;
	  break;
	}
    }

  // Up from maximum (towards +ve dcp)
  // Don't need to wrap - this is covered in the dcpbin_low cal
  for (int i=max_dcp; i<nbinsdcp+1; i++)
    {
      //std::cout << i << "  " << h68dcp->GetBinContent(i) << std::endl;
      if (h68dcp->GetBinContent(i)!=0) 
	{
	  dcpbin_high = i;
	  //std::cout << "dcpbin_high: " << h68dcp->GetXaxis()->GetBinLowEdge(i) << " -- " << h68dcp->GetXaxis()->GetBinUpEdge(i) << std::endl;
	  break;
	}
    }

  // 90%: down from maximum
  // h90th contains the white historam -- first non-zero bins out from the maximum give the 90% limits
  for (int i=max_dcp; i<nbinsdcp+1; i--)
    {
      //std::cout << i << "  " << h90dcp->GetBinContent(i) << std::endl;
      // wrap if necessary                                                                      
      int bin= i;
      if (i<=0)
        {
          bin = i+nbinsdcp+1;
	  //std::cout << h90dcp->GetBinCenter(i) << "   "<< h90dcp->GetBinCenter(i)+TMath::Pi() << "   "<< h90dcp->GetBinCenter(bin) << std::endl;
        }

      if (h90dcp->GetBinContent(bin)!=0) 
	{
	  dcpbin90_low = bin;
	  //std::cout << "dcpbin90_low: " << h68dcp->GetXaxis()->GetBinLowEdge(bin) << " -- " << h68dcp->GetXaxis()->GetBinUpEdge(bin) << std::endl;
	  break;
	}
    }
  // 90%: up from maximum (no need to wrap)
  for (int i=max_dcp; i<nbinsdcp+1; i++)
    {
      //std::cout << i << "  " << h90dcp->GetBinContent(i) << std::endl;
      if (h90dcp->GetBinContent(i)!=0) 
	{
	  dcpbin90_high = i;
	  //std::cout << "dcpbin90_high: " << h68dcp->GetXaxis()->GetBinLowEdge(i) << " -- " << h68dcp->GetXaxis()->GetBinUpEdge(i) << std::endl;
	  break;
	}
    }

  for (int i=max_dcp; i<nbinsdcp+1; i--)
    {
      int bin= i;
      if (i<=0) bin = i+nbinsdcp+1;
      if (h95dcp->GetBinContent(bin)!=0) 
	{
	  dcpbin95_low = bin;
	  break;
	}
    }
  for (int i=max_dcp; i<nbinsdcp+1; i++)
    {
      if (h95dcp->GetBinContent(i)!=0) 
	{
	  dcpbin95_high = i;
	  break;
	}
    }

  for (int i=max_dcp; i<nbinsdcp+1; i--)
    {
      int bin= i;
      if (i<=0) bin = i+nbinsdcp+1;
      if (h99dcp->GetBinContent(bin)!=0) 
	{
	  dcpbin99_low = bin;
	  break;
	}
    }
  for (int i=max_dcp; i<nbinsdcp+1; i++)
    {
      if (h99dcp->GetBinContent(i)!=0) 
	{
	  dcpbin99_high = i;
	  break;
	}
    }

  for (int i=max_dcp; i<nbinsdcp+1; i--)
    {
      int bin= i;
      if (i<=0) bin = i+nbinsdcp+1;
      if (h3sigdcp->GetBinContent(bin)!=0) 
	{
	  dcpbin3sig_low = bin;
	  break;
	}
    }
  for (int i=max_dcp; i<nbinsdcp+1; i++)
    {
      if (h3sigdcp->GetBinContent(i)!=0) 
	{
	  dcpbin3sig_high = i;
	  break;
	}
    }


  double minusint=h_dcp->Integral(1,100);//!!??
  //double minusint2=h_dcp->Integral(89,100);//!!??
  double plusint=h_dcp->Integral(62,88);//!!??
  std::cout<<minusint<<" "<<plusint<<" "<<(minusint)/(plusint)<<std::endl;

  TVectorD* v68th = new TVectorD(2);
  TVectorD* v90th = new TVectorD(2);
  TVectorD* v95th = new TVectorD(2);
  TVectorD* v99th = new TVectorD(2);
  TVectorD* v3sigth = new TVectorD(2);
 
  TVectorD* v68dcp = new TVectorD(2);
  TVectorD* v90dcp = new TVectorD(2);
  TVectorD* v95dcp = new TVectorD(2);
  TVectorD* v99dcp = new TVectorD(2);
  TVectorD* v3sigdcp = new TVectorD(2);
 
  std::cout << "--- "<<contlevel1*100.<<"% Limits ---" <<std::endl;
  // thbin_low and thbin_high are on outside of 68% limits
  (*v68th)(0) = h_th13->GetXaxis()->GetBinUpEdge(thbin_low) ;
  (*v68th)(1) = h_th13->GetXaxis()->GetBinLowEdge(thbin_high) ;
  std::cout << "sin^2 theta_13: " << (*v68th)(0) 	    << " -- "	    << (*v68th)(1)	    << std::endl;
  
  // check if had to wrap dcp and correct if so
  // dcpbin_low and dcp_bin high are also on outside of 68% limits
  if (dcpbin_low > dcpbin_high) // had to wrap dcp
    {
      (*v68dcp)(0) = h_dcp->GetXaxis()->GetBinLowEdge(dcpbin_high) ;
      (*v68dcp)(1) = h_dcp->GetXaxis()->GetBinLowEdge(dcpbin_high) ;
      
      std::cout << "delta cp: -pi -- " << (*v68dcp)(0) << " & " << (*v68dcp)(1) << " -- pi" << std::endl;
    }
  else // as normal
    {
      (*v68dcp)(0) = h_dcp->GetXaxis()->GetBinUpEdge(dcpbin_low) ; 
      (*v68dcp)(1) = h_dcp->GetXaxis()->GetBinLowEdge(dcpbin_high) ;
      
      std::cout << "delta cp: " << (*v68dcp)(0) 		<< " -- "		<< (*v68dcp)(1) 		<< std::endl;
    }

  std::cout << "--- "<<contlevel2*100.<<"% Limits ---" <<std::endl;
  // thbin_low and thbin_high are on outside of 90% limits
  (*v90th)(0) = h_th13->GetXaxis()->GetBinUpEdge(thbin90_low) ;
  (*v90th)(1) = h_th13->GetXaxis()->GetBinLowEdge(thbin90_high) ;
               
  std::cout << "sin^2 theta_13: " << (*v90th)(0) 
	    << " -- "
	    << (*v90th)(1) 
	    << std::endl;
    
  // check if had to wrap dcp and correct if so 
  // dcpbin_low and dcp_bin high are also on outside of 90% limits
  if (dcpbin90_low > dcpbin90_high) // had to wrap dcp                                                                                                                    
    {
      (*v90dcp)(0) = h_dcp->GetXaxis()->GetBinLowEdge(dcpbin90_high) ;
      (*v90dcp)(1) = h_dcp->GetXaxis()->GetBinUpEdge(dcpbin90_low) ;
      std::cout << "delta cp: -pi -- " << (*v90dcp)(0) << " & " << (*v90dcp)(1) << " -- pi" << std::endl;
    }
  else // as normal                                                                                                                                  
    {
      (*v90dcp)(0) = h_dcp->GetXaxis()->GetBinUpEdge(dcpbin90_low) ; 
      (*v90dcp)(1) = h_dcp->GetXaxis()->GetBinLowEdge(dcpbin90_high) ;
      std::cout << "delta cp: " << (*v90dcp)(0) 		<< " -- "		<< (*v90dcp)(1) 		<< std::endl;
    }

  std::cout << "--- "<<contlevel3*100.<<"% Limits ---" <<std::endl;
  // thbin_low and thbin_high are on outside of 95% limits
  (*v95th)(0) = h_th13->GetXaxis()->GetBinUpEdge(thbin95_low) ;
  (*v95th)(1) = h_th13->GetXaxis()->GetBinLowEdge(thbin95_high) ;
               
  std::cout << "sin^2 theta_13: " << (*v95th)(0) 
	    << " -- "
	    << (*v95th)(1) 
	    << std::endl;
    
  // check if had to wrap dcp and correct if so 
  // dcpbin_low and dcp_bin high are also on outside of 95% limits
  if (dcpbin95_low > dcpbin95_high) // had to wrap dcp                                                                                                                    
    {
      (*v95dcp)(0) = h_dcp->GetXaxis()->GetBinLowEdge(dcpbin95_high) ;
      (*v95dcp)(1) = h_dcp->GetXaxis()->GetBinUpEdge(dcpbin95_low) ;
      std::cout << "delta cp: -pi -- " << (*v95dcp)(0) << " & " << (*v95dcp)(1) << " -- pi" << std::endl;
    }
  else // as normal                                                                                                                                  
    {
      (*v95dcp)(0) = h_dcp->GetXaxis()->GetBinUpEdge(dcpbin95_low) ; 
      (*v95dcp)(1) = h_dcp->GetXaxis()->GetBinLowEdge(dcpbin95_high) ;
      std::cout << "delta cp: " << (*v95dcp)(0) 		<< " -- "		<< (*v95dcp)(1) 		<< std::endl;
    }

  std::cout << "--- "<<contlevel4*100.<<"% Limits ---" <<std::endl;
  // thbin_low and thbin_high are on outside of 99% limits
  (*v99th)(0) = h_th13->GetXaxis()->GetBinUpEdge(thbin99_low) ;
  (*v99th)(1) = h_th13->GetXaxis()->GetBinLowEdge(thbin99_high) ;
               
  std::cout << "sin^2 theta_13: " << (*v99th)(0) 
	    << " -- "
	    << (*v99th)(1) 
	    << std::endl;
    
  // check if had to wrap dcp and correct if so 
  // dcpbin_low and dcp_bin high are also on outside of 99% limits
  if (dcpbin99_low > dcpbin99_high) // had to wrap dcp                                                                                                                    
    {
      (*v99dcp)(0) = h_dcp->GetXaxis()->GetBinLowEdge(dcpbin99_high) ;
      (*v99dcp)(1) = h_dcp->GetXaxis()->GetBinUpEdge(dcpbin99_low) ;
      std::cout << "delta cp: -pi -- " << (*v99dcp)(0) << " & " << (*v99dcp)(1) << " -- pi" << std::endl;
    }
  else // as normal                                                                                                                                  
    {
      (*v99dcp)(0) = h_dcp->GetXaxis()->GetBinUpEdge(dcpbin99_low) ; 
      (*v99dcp)(1) = h_dcp->GetXaxis()->GetBinLowEdge(dcpbin99_high) ;
      std::cout << "delta cp: " << (*v99dcp)(0) 		<< " -- "		<< (*v99dcp)(1) 		<< std::endl;
    }

  std::cout << "--- "<<contlevel5*100.<<"% Limits ---" <<std::endl;
  // thbin_low and thbin_high are on outside of 99.73% limits
  (*v3sigth)(0) = h_th13->GetXaxis()->GetBinUpEdge(thbin3sig_low) ;
  (*v3sigth)(1) = h_th13->GetXaxis()->GetBinLowEdge(thbin3sig_high) ;
               
  std::cout << "sin^2 theta_13: " << (*v3sigth)(0) 
	    << " -- "
	    << (*v3sigth)(1) 
	    << std::endl;
    
  // check if had to wrap dcp and correct if so 
  // dcpbin_low and dcp_bin high are also on outside of 99.73% limits
  if (dcpbin3sig_low > dcpbin3sig_high) // had to wrap dcp                                                                                                                    
    {
      (*v3sigdcp)(0) = h_dcp->GetXaxis()->GetBinLowEdge(dcpbin3sig_high) ;
      (*v3sigdcp)(1) = h_dcp->GetXaxis()->GetBinUpEdge(dcpbin3sig_low) ;
      std::cout << "delta cp: -pi -- " << (*v3sigdcp)(0) << " & " << (*v3sigdcp)(1) << " -- pi" << std::endl;
    }
  else // as normal                                                                                                                                  
    {
      (*v3sigdcp)(0) = h_dcp->GetXaxis()->GetBinUpEdge(dcpbin3sig_low) ; 
      (*v3sigdcp)(1) = h_dcp->GetXaxis()->GetBinLowEdge(dcpbin3sig_high) ;
      std::cout << "delta cp: " << (*v3sigdcp)(0) 		<< " -- "		<< (*v3sigdcp)(1) 		<< std::endl;
    }



  TFile *fout = new TFile("contours_app1D.root", "recreate");
  fout->cd();
  h_th13->Write("h_th13");
  h68th->Write("h68th");
  h90th->Write("h95th");
  h90th->Write("h99th");
  h90th->Write("h3sigth");
  v68th->Write("v68th");
  v90th->Write("v90th");
  v90th->Write("v95th");
  v90th->Write("v99th");
  v90th->Write("v3sigth");

  h_dcp->Write("h_dcp");
  h68dcp->Write("h68dcp");
  h90dcp->Write("h90dcp");
  h95dcp->Write("h95dcp");
  h99dcp->Write("h99dcp");
  h3sigdcp->Write("h3sigdcp");

  v68dcp->Write("v68dcp");
  v90dcp->Write("v90dcp");
  v90dcp->Write("v95dcp");
  v90dcp->Write("v99dcp");
  v90dcp->Write("v3sigdcp");

  leg->Write("legend") ;
  leg_sig->Write("legend_sig") ;


  // ------- chi2 for comparison to other groups ------- //

  // dcp only (this is all they do)
  TH1D *h_dcp_chi2 = (TH1D*)h_dcp->Clone("h_dcp_chi2");
  h_dcp_chi2->Reset();
  h_dcp_chi2->GetYaxis() -> SetTitle("#Delta #chi^{2}") ;
  std::cout << "-----------------" << std::endl;
  double maxbin = h_dcp->GetMaximumBin();
  double maxdcp = h_dcp->GetBinContent(maxbin);
  std::cout << "maxdcp = " << maxdcp << std::endl;
  for (int i=0; i<= h_dcp->GetNbinsX(); i++){
    double bc = h_dcp->GetBinContent(i);
    std::cout << "bc = " << bc << std::endl;
    if (bc==0){
      bc = 1.0;
    }
    double l = -2*TMath::Log(bc/maxdcp);
    std::cout << "-2LLH = " << l << std::endl;
    h_dcp_chi2->SetBinContent(i,l);
  }
  std::cout << "-----------------" << std::endl;

  h_dcp_chi2->Write("h_dcp_chi2");

}
