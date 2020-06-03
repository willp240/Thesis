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

void getInterval1D(TH1D *hist, double &p68, double &p90)
{
  std::cout << "getting interval" << std::endl;
  TH1D *hCopy = (TH1D*)hist->Clone("hCopy");

  double integral = hCopy->Integral();
  double tsum = 0;
  double cont68lev = 0;
  double cont90lev = 0;

  std::cout << integral << std::endl;

  while((tsum / integral) < 0.954)
    {
      double tmax = hCopy->GetMaximum();
      tsum = tsum + tmax;
      int bin = hCopy->GetMaximumBin();
      if (tsum / integral < 0.68)
	{
	  cont68lev = tmax;
	  hCopy->SetBinContent(bin, -1.0);
	}
      if ((tsum / integral < 0.954) && (tsum / integral > 0.68))
	{
	  cont90lev = tmax;
	  hCopy->SetBinContent(bin, -3.0);
	}
    }

  p90 = cont90lev;
  p68 = cont68lev;

  std::cout << "p68 = " << p68 << ", p90 = " << p90 << std::endl;
}

void getInterval1D(TH1D *hist, TH1D &h68, TH1D &h90)
{
  std::cout << "getting interval" << std::endl;
  TH1D *hCopy = (TH1D*)hist->Clone("hCopy");
  TH1D *hCopy68 = (TH1D*)hist->Clone("hCopy68");
  TH1D *hCopy90 = (TH1D*)hist->Clone("hCopy90");

  double integral = hCopy->Integral();
  double tsum = 0;
  double cont68lev = 0;
  double cont90lev = 0;

  std::cout << integral << std::endl;

  while((tsum / integral) < 0.954)
    {
      double tmax = hCopy->GetMaximum();
      tsum = tsum + tmax;
      int bin = hCopy->GetMaximumBin();
      if (tsum / integral < 0.68)
	{
	  cont68lev = tmax;
	  hCopy->SetBinContent(bin, -1.0);
	  hCopy68->SetBinContent(bin, 0);
	  hCopy90->SetBinContent(bin, 0);
	}
      if ((tsum / integral < 0.954) && (tsum / integral > 0.68))
	{
	  cont90lev = tmax;
	  hCopy->SetBinContent(bin, -3.0);
	  hCopy90->SetBinContent(bin, 0);
	}
    }

  h90 = (*hCopy90);
  h68 = (*hCopy68);
}

void contours_1D(TString fs,int onehierarchyonly=0,int burnIn=0, bool doRCrw=false)
{
  gStyle->SetPalette(51,0);
  gStyle->SetOptStat(0);

  TFile *f = new TFile(fs.Data(), "READ");
  TTree *t = (TTree*)f->Get("osc_posteriors");

  TH1D *h_th23 = new TH1D("h_th23","",300,0,1);
  TH1D *h_dm32 = new TH1D("h_dm32","",1500,-3e-3,3e-3);

  int max = t->GetEntries();

  double tth23, tdm32, tth13, tdcp, tbiaswgt=1.0;
  int stp;
  double RCrw=1.0;

  t->SetBranchAddress("dm23",&tdm32);
  t->SetBranchAddress("theta23",&tth23);
  t->SetBranchAddress("theta13",&tth13);
  t->SetBranchAddress("dcp",&tdcp);
  t->SetBranchAddress("step",&stp);
  if(doRCrw) t->SetBranchAddress("RCreweight", &RCrw);

  // Apply different prior to osc pars 
  // Note: this is assuming it is being run on the fit without reactor (so no prior to ‘undo’)
  
  // Flat in theta23 or theta13 (not sin^2 theta23 or sin^2 theta13) 
  TF1 *prior1 = new TF1("prior1","1.0/1.0/(2*TMath::Sqrt(x - x*x))");
  
  // Flat in sin^2 2theta13 (not sin^2 theta13) 
  TF1 *prior2 = new TF1("prior2", "4 - 8*x"); 
  
  // Flat in sin dcp (not dcp)
  TF1 *prior3 = new TF1("prior3","TMath::Abs(TMath::Cos(x))");
  
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
      //tbiaswgt = prior-1>Eval(tth23);

      // Flat in theta13 (not sin^2 theta13)
      //tbiaswgt = prior1->Eval(tth13);

      // Flat in sin^2 2theta13 (not sin^2 theta13) 
      //tbiaswgt = prior2->Eval(tth13);
      // Flat in sin dcp (not dcp)
      //      tbiaswgt = prior3->Eval(tdcp);  
 

      double totWeight = tbiaswgt * RCrw;
      h_th23->Fill(tth23, totWeight);
      h_dm32->Fill(tdm32, totWeight);
    }

  if(onehierarchyonly==-1) h_dm32->GetXaxis()->SetRangeUser(-0.003,-0.002);
  if(onehierarchyonly==1) h_dm32->GetXaxis()->SetRangeUser(0.002,0.003);
 

  //h_th23->Scale(1./h_th23->Integral(), "width") ;
  //h_dm32->Scale(1./h_dm32->Integral(), "width") ;

  h_th23->Scale(1./h_th23->Integral()) ;
  h_dm32->Scale(1./h_dm32->Integral()) ;
   
  h_th23->GetXaxis()->SetTitle("sin^{2}#theta_{23}");
  h_th23-> GetXaxis() -> SetTitleOffset(1.2) ;
  h_th23-> GetXaxis() -> SetTitleFont(132) ;
  h_th23-> GetXaxis() -> SetTitleSize(0.06) ;
  h_th23-> GetXaxis() -> SetLabelFont(132) ;
  h_th23-> GetXaxis() -> SetLabelSize(0.05) ;
  h_th23->GetYaxis()->SetTitle("Posterior probability density");
  h_th23-> GetYaxis() -> SetTitleOffset(1.2) ;
  h_th23-> GetYaxis() -> SetTitleFont(132) ;
  h_th23-> GetYaxis() -> SetTitleSize(0.06) ;
  h_th23-> GetYaxis() -> SetLabelFont(132) ;
  h_th23-> GetYaxis() -> SetLabelSize(0.05) ;

  h_dm32->GetXaxis()->SetTitle("#Delta m^{2}_{32} (eV^{2})");
  h_dm32-> GetXaxis() -> SetTitleOffset(1.2) ;
  h_dm32-> GetXaxis() -> SetTitleFont(132) ;
  h_dm32-> GetXaxis() -> SetTitleSize(0.06) ;
  h_dm32-> GetXaxis() -> SetLabelFont(132) ;
  h_dm32-> GetXaxis() -> SetLabelSize(0.05) ;
  h_dm32->GetYaxis()->SetTitle("Posterior probability density");
  h_dm32-> GetYaxis() -> SetTitleOffset(1.2) ;
  h_dm32-> GetYaxis() -> SetTitleFont(132) ;
  h_dm32-> GetYaxis() -> SetTitleSize(0.06) ;
  h_dm32-> GetYaxis() -> SetLabelFont(132) ;
  h_dm32-> GetYaxis() -> SetLabelSize(0.05) ;

  /*
  h_th23->GetYaxis()->SetTitle("Postfit N_{steps}");
  h_th23->GetXaxis()->SetTitle("sin^{2}#theta_{23}");
  h_dm32->GetYaxis()->SetTitle("Postfit N_{steps}");
  h_dm32->GetXaxis()->SetTitle("#Delta m_{32} (eV^{2})");
  */
  TH1D *h68th = (TH1D*)h_th23->Clone("h68th");
  TH1D *h90th = (TH1D*)h_th23->Clone("h90th");
  TH1D *h68dm = (TH1D*)h_dm32->Clone("h68dm"); 
  TH1D *h90dm = (TH1D*)h_dm32->Clone("h90dm");
  TH1D *copyth = (TH1D*)h_th23->Clone("copyth");
  TH1D *copydm = (TH1D*)h_dm32->Clone("copydm");

  double integral, tsum=0;

  // Get intervals (th_23)
  std::cout << "getting interval" << std::endl;

  integral = copyth->Integral();

  std::cout << integral << std::endl;

  while((tsum / integral) < 0.954)
    {
      double tmax = copyth->GetMaximum();
      tsum = tsum + tmax;
      int bin = copyth->GetMaximumBin();
      if (tsum / integral < 0.68)
	{
	  copyth->SetBinContent(bin, -1.0);
	  h68th->SetBinContent(bin, 0);
	  h90th->SetBinContent(bin, 0);
	}
      if ((tsum / integral < 0.954) && (tsum / integral > 0.68))
	{
	  copyth->SetBinContent(bin, -3.0);
	  h90th->SetBinContent(bin, 0);
	}
    }


  integral=0; 
  tsum=0;

  // Get intervals (dm_32)
  std::cout << "getting interval" << std::endl;

  integral = copydm->Integral();

  std::cout << integral << std::endl;

  while((tsum / integral) < 0.954)
    {
      double tmax = copydm->GetMaximum();
      tsum = tsum + tmax;
      int bin = copydm->GetMaximumBin();
      if (tsum / integral < 0.68)
	{
	  copydm->SetBinContent(bin, -1.0);
	  h68dm->SetBinContent(bin, -0.001);
	  h90dm->SetBinContent(bin, -0.001);
	}
      if ((tsum / integral < 0.954) && (tsum / integral > 0.68))
	{
	  copydm->SetBinContent(bin, -3.0);
	  h90dm->SetBinContent(bin, -0.001);
	}
    }


  TCanvas *c = new TCanvas("c","",1200,800);
  TCanvas *c1 = new TCanvas("c1","",1200,800);
  c->SetLeftMargin(0.15) ;
  c1->SetLeftMargin(0.15) ;
  c->SetBottomMargin(0.15) ;
  c1->SetBottomMargin(0.15) ;
  
  h_th23->SetLineColor(kBlack);
  //  h_th23->GetYaxis()->SetTitleOffset(1.7);
  h68th->SetLineColor(kBlack);
  h90th->SetLineColor(kBlack);
  //h68th->SetLineStyle(2);
  //h90th->SetLineStyle(3);
  h68th->SetFillColor(11);
  h_th23->SetFillColor(12);
  h90th->SetFillColor(10);

  h_dm32->SetLineColor(kBlack);
  //  h_dm32->GetYaxis()->SetTitleOffset(1.8);
  //  h_dm32->GetXaxis()->SetNdivisions(505);
  //h_dm32->GetXaxis()->SetRangeUser(0,5e-3);
  h68dm->SetLineColor(kBlack);
  h90dm->SetLineColor(kBlack);
  //h68dm->SetLineStyle(2);
  //h90dm->SetLineStyle(3);
  h68dm->SetFillColor(11);
  h_dm32->SetFillColor(12);
  h90dm->SetFillColor(10);

  TLegend *leg = new TLegend(0.4,0.6,0.6,0.85);
  if(onehierarchyonly==-1 || onehierarchyonly==1) h_dm32->SetMaximum(h_dm32->GetMaximum()*1.55);
  leg->SetFillColor(kWhite);
  leg->SetLineColor(kWhite);
  leg->SetTextFont(132);
  leg->SetTextSize(0.04);
  leg->AddEntry(h68dm,"90% Credible Interval","f");
  leg->AddEntry(h_dm32,"68% Credible Interval","f");

  TLegend *leg2 = new TLegend(0.2,0.6,0.4,0.85);
  leg2->SetFillColor(kWhite);
  leg2->SetLineColor(kWhite);
  leg2->SetTextFont(132);
  leg2->SetTextSize(0.04);
  leg2->AddEntry(h68dm,"90% Credible Interval","f");
  leg2->AddEntry(h_dm32,"68% Credible Interval","f");
  c->cd();
  //  gPad->SetLeftMargin(0.13);
  gPad->SetTickx() ;
  gPad->SetTicky() ;
  h_th23->Draw();
  h68th->Draw("SAME");
  h90th->Draw("SAME");
  leg2->Draw("SAME");
  gPad->RedrawAxis();
  c->Update();
  gPad->Update();

  c1->cd();
  gPad->SetTickx() ;
  gPad->SetTicky() ;
  //  gPad->SetLeftMargin(0.13);
  h_dm32->GetXaxis()->LabelsOption("u");
  h_dm32->Draw();
  h68dm->Draw("SAME");
  h90dm->Draw("SAME");
  leg->Draw("SAME");
  gPad->RedrawAxis();
  c1->Update();
  gPad->Update();

  c->SaveAs("contours_1D_th23.pdf");
  c1->SaveAs("contours_1D_dm23.pdf");
  c->SaveAs("contours_1D_th23.root");
  c1->SaveAs("contours_1D_dm23.root");
  //c->SaveAs("contours_1D_th23.png");
  //c1->SaveAs("contours_1D_dm23.png");

  // Get limits
  int thbin_low=0, thbin_high=0, dmbin_low_IH=0, dmbin_high_IH=0, dmbin_low_NH=0, dmbin_high_NH=0;
  int thbin90_low=0, thbin90_high=0, dmbin90_low_IH=0, dmbin90_high_IH=0, dmbin90_low_NH=0, dmbin90_high_NH=0;
  int nbinsth = h_th23->GetNbinsX();
  int nbinsdm = h_dm32->GetNbinsX();
  int centre_th = h_th23->GetMaximumBin();
  int centre_dm = h_dm32->FindBin(0);

  std::cout << centre_th << "  " << centre_dm << std::endl;

  // h68th contains the beige histogram -- first non-zero bins out from the centre give 68% limits
  for (int i=centre_th; i>0; i--)
    {
      //std::cout << i << "  " << h68th->GetBinContent(i) << std::endl;
      if (h68th->GetBinContent(i)!=0) 
	{
	  thbin_low = i;
	  break;
	}
    }
  // h90th contains the white histogram -- first non-zero bins out from the centre give 90% limits
  for (int i=centre_th; i>0; i--)
    {
      //std::cout << i << "  " << h90th->GetBinContent(i) << std::endl;
      if (h90th->GetBinContent(i)!=0) 
	{
	  thbin90_low = i;
	  break;
	}
    }
  // h90th contains the beige histogram -- first non-zero bins out from the centre give 68% limits
  for (int i=centre_th; i<nbinsth+1; i++)
    {
      //std::cout << i << "  " << h68th->GetBinContent(i) << std::endl;
      if (h68th->GetBinContent(i)!=0) 
	{
	  thbin_high = i;
	  break;
	}
    }
  // h90th contains the white histogram -- first non-zero bins out from the centre give 90% limits
  for (int i=centre_th; i<nbinsth+1; i++)
    {
      //std::cout << i << "  " << h68th->GetBinContent(i) << std::endl;
      if (h90th->GetBinContent(i)!=0) 
	{
	  thbin90_high = i;
	  break;
	}
    }

  // down from centre to find 68% limit
  // h68dm contains the beige histogram -- first zero bin after non-zero bins when counting down from centre gives upper IH 68% limit
  // Next non-zero bin if you keep counting gives lower IH 68% limit
  bool nonzero=false;
  bool zero=false;
 for (int i=centre_dm; i>0; i--)
    {
      //std::cout << i << "  " << h68dm->GetBinContent(i) << std::endl;
      if (h68dm->GetBinContent(i)!=0) nonzero=true;
      if (nonzero && !zero && h68dm->GetBinContent(i)==-0.001)
	{
	  dmbin_high_IH = i;
	  zero=true;
	  // std::cout << i << "  " << h68dm->GetXaxis()->GetBinLowEdge(i) << "  " << h68dm->GetXaxis()->GetBinUpEdge(i) << std::endl;
	}
      if (zero && h68dm->GetBinContent(i)!=-0.001)
	{
	  dmbin_low_IH = i;
	  //std::cout << i << "  " << h68dm->GetXaxis()->GetBinLowEdge(i) << "  " << h68dm->GetXaxis()->GetBinUpEdge(i) << std::endl;
	  break;
	}
    }

 // do the same thing going up from centre to find NH 68% limits
 nonzero=false;
 zero=false;
 for (int i=centre_dm; i<nbinsdm+1; i++)
   {
     if (h68dm->GetBinContent(i)!=0){
       nonzero=true;
     }
     if (nonzero && !zero && h68dm->GetBinContent(i)==-0.001)
       {
	 dmbin_low_NH = i;
	 zero=true;
       }
     if (zero && h68dm->GetBinContent(i)!=-0.001)
       {
	 dmbin_high_NH = i;
	 break;
       }
   }

 // Now do the same to find 90% limits
 // h90dm contains the white histogram -- first zero bin after non-zero bins gives the first 90% limit (counting out from centr, eg. low 90% limit in NH or high 90% limit in IH)
 // Next non-zero bin if you keep counting gives the other 90% limit
 // Need to check first zero bin isn
 nonzero=false;
 zero=false;
 for (int i=centre_dm; i>0; i--)
   {
     if (h90dm->GetBinContent(i)!=0){
       nonzero=true;
     }
     if (nonzero && !zero && h90dm->GetBinContent(i)==-0.001)
       {
         dmbin90_high_IH = i;
         zero=true;
       }
     if (zero && h90dm->GetBinContent(i)!=-0.001)
       {
         dmbin90_low_IH =i;
	 break;
       }
   }
 // Do the same thing going up from the centre to find the NH 90% limits
 nonzero=false;
 zero=false;
 for (int i=centre_dm; i<nbinsdm+1; i++)
   {
     if (h90dm->GetBinContent(i)!=0) nonzero=true;
     if (nonzero && !zero && h90dm->GetBinContent(i)==-0.001)
       {
         dmbin90_low_NH = i;
         zero=true;
       }
     if (zero && h90dm->GetBinContent(i)!=-0.001)
       {
         dmbin90_high_NH =i;
	 break;
       }
   }

 TVectorD* v68th = new TVectorD(2);
 TVectorD* v90th = new TVectorD(2);
 
 TVectorD* v90dm32_ih = new TVectorD(2);
 TVectorD* v68dm32_ih = new TVectorD(2);
 TVectorD* v90dm32_nh = new TVectorD(2);
 TVectorD* v68dm32_nh = new TVectorD(2);

 std::cout << "!!! be careful, there seem to be an issue with DM32 limits display. " << std::endl ;
  std::cout << "--- 68% Limits ---" <<std::endl;
  //thbin_low and thbin_high are on the outside of the 68% limits

  (*v68th)(0) = h_th23->GetXaxis()->GetBinUpEdge(thbin_low) ;
  (*v68th)(1) = h_th23->GetXaxis()->GetBinLowEdge(thbin_high) ;
  std::cout << "sin^2 theta_23: " << (*v68th)(0) 	    << " -- "	    << (*v68th)(1)	    << std::endl;
  
  // IH: dmbin_low is on the outside of the 68% limit and dmbin_high is on the inside of the 68% limit (want BinUpEdge for both)
  // NH: dmbin_low is on the inside of the 68% limit and dmbin_high is on the outside of the 68% limit (want BinLowEdge for both)

  (*v68dm32_ih)(0) = h_dm32->GetXaxis()->GetBinLowEdge(dmbin_low_IH) ;
  (*v68dm32_ih)(1) =  h_dm32->GetXaxis()->GetBinUpEdge(dmbin_high_IH) ;
  std::cout << "delta m^2_32 (IH): " 	    << (*v68dm32_ih)(0) 		<< " -- "		<< (*v68dm32_ih)(1) 		<< std::endl;


  (*v68dm32_nh)(0) = h_dm32->GetXaxis()->GetBinLowEdge(dmbin_low_NH) ;
  (*v68dm32_nh)(1) =  h_dm32->GetXaxis()->GetBinUpEdge(dmbin_high_NH) ;
  std::cout << "delta m^2_32 (NH): " 	    << (*v68dm32_nh)(0) 		<< " -- "		<< (*v68dm32_nh)(1) 		<< std::endl;

  std::cout << "--- 90% Limits ---" <<std::endl;
  // thbin90_low and thbin90_high are on the outside of the 90% limits

  (*v90th)(0) = h_th23->GetXaxis()->GetBinUpEdge(thbin90_low) ;
  (*v90th)(1) = h_th23->GetXaxis()->GetBinLowEdge(thbin90_high) ;
  std::cout << "sin^2 theta_23: " << (*v90th)(0) 	    << " -- "	    << (*v90th)(1) 	    << std::endl;
  // IH: dmbin90_low is on the outside of the 90% limit and dmbin90_high is on the inside of the 90% limit (want BinUpEdge for both) 
  // NH: dmbin90_low is on the inside of the 90% limit and dmbin90_high is on the outside of the 90% limit (want BinLowEdge for both)

  (*v68dm32_ih)(0) = h_dm32->GetXaxis()->GetBinLowEdge(dmbin90_low_IH) ;
  (*v68dm32_ih)(1) =  h_dm32->GetXaxis()->GetBinUpEdge(dmbin90_high_IH) ;
  std::cout << "delta m^2_32 (IH): " 	    << (*v68dm32_ih)(0) 		<< " -- "		<< (*v68dm32_ih)(1) 		<< std::endl;

  (*v90dm32_nh)(0) = h_dm32->GetXaxis()->GetBinLowEdge(dmbin90_low_NH) ;
  (*v90dm32_nh)(1) =  h_dm32->GetXaxis()->GetBinUpEdge(dmbin90_high_NH) ;
  std::cout << "delta m^2_32 (NH): " 	    << (*v90dm32_nh)(0) 		<< " -- "		<< (*v90dm32_nh)(1) 		<< std::endl;
  
  TFile *fout = new TFile("contours_1D.root", "recreate");
  fout->cd();
  h_th23->Write("h_th23");
  h68th->Write("h68th");
  h90th->Write("h90th");
  v68th->Write("v68th");
  v90th->Write("v90th");
  h_dm32->Write("h_dm32");
  h68dm->Write("h68dm");
  h90dm->Write("h90dm");
  ///!!! problem with those boundaries - save when filled properly
  /*  v68dm32_nh->Write("v68dm_nh");
  v90dm32_nh->Write("v90dm_nh");
  v68dm32_ih->Write("v68dm_ih");
  v90dm32_ih->Write("v90dm_ih");*/
  leg->Write("legend") ;
  
}
