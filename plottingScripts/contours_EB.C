#include "TH2D.h"
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
#include <iostream>

//******************************************************************************\\
// Plot Delta m_23^2 VS sin^2 theta_23 for Asimov fit
//
// To run :
// $ root
//  root[1] .L contours.C
//  root[2] llh_ratio("fit_output.root", "outFileName", true, 0.528, 2.4e-3)
//     where "fit_output.root" is the name of the input file, from reducedSyst
//           "outFileName"     is the name of the output file
//           true = fit is Re+Rmu, false = fit is Rmu only
//           0.528 = input sin^2 theta_23
//           2.4e-3 = input Delta m_23^2
//******************************************************************************\\

//TH2D *h68;
//TH2D *h90;

int colour_68 = 2; // 68% intervals are white
int colour_90 = 2; // 80% intervals are white also
int colour_95 = 2; // 95% intervals are white also
int colour_99 = 2; // 99% intervals are white also
int colour_3sig = 2; // 3-sig intervals are white also


double* getInterval2D(TH2D *hist, double &p68, double &p90, double &p95, double &p99 ,double &p3sig)//, TH2D &h68, TH2D &h90)
{
  std::cout << "getting interval" << std::endl;
  TH2D *hCopy = (TH2D*)hist->Clone("hCopy");
  TH2D *hCont68 = (TH2D*)hist->Clone("hCont68");
  TH2D *hCont90 = (TH2D*)hist->Clone("hCont90");
  TH2D *hCont95 = (TH2D*)hist->Clone("hCont95");
  TH2D *hCont99 = (TH2D*)hist->Clone("hCont99");
  TH2D *hCont3sig = (TH2D*)hist->Clone("hCont3sig");
  
  double integral = hCopy->Integral();
  double tsum = 0;
  double cont68lev = 0;// array('d', [0.0])
  double cont90lev = 0;//array('d', [0.0])
  double cont95lev = 0;//array('d', [0.0])
  double cont99lev = 0;//array('d', [0.0])
  double cont3siglev = 0;//array('d', [0.0])

  std::cout << integral << " " << tsum << std::endl;

  while ((tsum / float(integral)) < 0.9973)
    {
      double tmax = hCopy->GetMaximum();
      //std::cout << tmax << std::endl;
      tsum = tsum + tmax;
      int bin = hCopy->GetMaximumBin();
      if (tsum / float(integral) < 0.68)
	{
	  cont68lev = tmax;
	  hCopy->SetBinContent(bin, -1.0);
	}
      if ((tsum / float(integral) < 0.9) && (tsum / float(integral) > 0.68))
	{
	  cont90lev = tmax;
	  hCopy->SetBinContent(bin, -3.0);
	}
      if ((tsum / float(integral) < 0.95) && (tsum / float(integral) > 0.9))
	{
	  cont95lev = tmax;
	  hCopy->SetBinContent(bin, -5.0);
	}
      if ((tsum / float(integral) < 0.99) && (tsum / float(integral) > 0.95))
	{
	  cont99lev = tmax;
	  hCopy->SetBinContent(bin, -7.0);
        }	
      if ((tsum / float(integral) < 0.9973) && (tsum / float(integral) > 0.99))
	{
	  cont3siglev = tmax;
	  hCopy->SetBinContent(bin, -9.0);
        }	

    }
  
  //std::cout << " c "  << cont90lev << " " << cont68lev << std::endl;
  double quant[5];
  quant[0] = cont90lev;
  quant[1] = cont68lev;

  // adding these, may need to rethink ordering later
  quant[2] = cont95lev;
  quant[3] = cont99lev;
  quant[4] = cont3siglev;
  //hCont68->Smooth();
  //hCont90->Smooth();

  //  h68 = hCont68;
  //  h90 = hCont90;

  p90 = cont90lev;
  p68 = cont68lev;

  p95 = cont95lev;
  p99 = cont99lev;
  p3sig = cont3siglev;


  std::cout << "p68 = " << p68 << ", p90 = " << p90 << std::endl;
  std::cout << "p95 = " << p95 << ", p99 = " << p99 << std::endl;
  std::cout << "p3sig = " << p3sig << std::endl;

  return quant;//, hCont68, hCont90;
}


void contours_EB(TString fs, double truex=0, double truey=0, int max_steps=-999,int onehierarchyonly=0,int burnIn=0, bool doRCrw=false)
{
  gStyle->SetPalette(51,0);
  gStyle->SetOptStat(0);
  TGaxis::SetMaxDigits(4);

  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleFillColor(0);



  // get mach3
  TFile *f = new TFile(fs.Data(), "READ");
  TTree *t = (TTree*)f->Get("osc_posteriors");

  //this is actually th23dm23
  TH2D *th13_dcp = new TH2D("th13_dcp","", 30,0.4,0.65,200,-3.015e-3,2.985e-3);
  //TH2D *th13_dcp = new TH2D("th13_dcp","", 60,0.4,0.65,400,-3.015e-3,2.985e-3);
  if(onehierarchyonly==1) th13_dcp->GetYaxis()->SetRangeUser(2.25e-3,2.815e-3); // NH
  if(onehierarchyonly==-1) th13_dcp->GetYaxis()->SetRangeUser(-2.75e-3,-2.2e-3); // IH
  // numubar disapp
  //TH2D *th13_dcp = new TH2D("th13_dcp","", 51,-0.002,1.018,50,-5.015e-3,4.985e-3);

  // 51,-0.002,1.018,50,-0.041e-3,4.959e-3);
  //123, 0.196296, 0.8037037, 26, 1.71429e-3, 3.842857e-3); // Davide's binning
  //51,-0.002,1.018,100,-5.041e-3,4.959e-3);
  //25, 0, 1, 100, 0, 20e-3); 
  //101, -0.005, 1.005, 400, -0.005e-3, 19.995E-3); 
  // this is secretly th23 dm32
  // closest bin centre to true values (0.524,2.51e-3) is at (0.525,2.50e-3)

  TH2D *llh_th13_dcp = (TH2D*)th13_dcp->Clone();

  // Extra histograms to calculate posterior probabilities by octant/MH
  TH2D *lowoct_NH = (TH2D*)th13_dcp->Clone("lowoct_NH");
  TH2D *highoct_NH = (TH2D*)th13_dcp->Clone("highoct_NH");
  TH2D *lowoct_IH = (TH2D*)th13_dcp->Clone("lowoct_IH");
  TH2D *highoct_IH = (TH2D*)th13_dcp->Clone("highoct_IH");

  int max;
  if (max_steps == -999) max = t->GetEntries();
  else max = max_steps;

  double tdm32;
  double tth23;
  double tth13;
  double tdcp;
  double tbeta = 1.0;
  double tnu_dm32;
  double tnu_th23;
  double tnu_th13;
  int stp;
  double tbiaswgt = 1.0;
  double RCrw = 1.0;
  double xsec_20;

  t->SetBranchAddress("dm23", &tdm32);
  t->SetBranchAddress("theta23", &tth23);
  t->SetBranchAddress("theta13", &tth13);
  t->SetBranchAddress("dcp", &tdcp);
  t->SetBranchAddress("step", &stp);
  t->SetBranchAddress("xsec_20",&xsec_20);
  //t->SetBranchAddress("beta", &tbeta);
  //t->SetBranchAddress("nu_theta23", &tnu_th23);
  //t->SetBranchAddress("nu_dm23", &tnu_dm32);
  //t->SetBranchAddress("nu_theta13", &tnu_th13);
  if(doRCrw)t->SetBranchAddress("RCreweight", &RCrw);

  // Apply different prior to osc pars 
  // Note: this is assuming it is being run on the fit without reactor (so no prior to ‘undo’)
  // Comment out unless needed
  
  // Flat in theta23 or theta13 (not sin^2 theta23 or sin^2 theta13) 
  TF1 *prior1 = new TF1("prior1","1.0/1.0/(2*TMath::Sqrt(x - x*x))");
  
  // Flat in sin^2 2theta13 (not sin^2 theta13) 
  TF1 *prior2 = new TF1("prior2", "4 - 8*x"); 
  
  // Flat in sin dcp (not dcp)
  TF1 *prior3 = new TF1("prior3","TMath::Abs(TMath::Cos(x))");
  
  // Flat in cos dcp (not dcp)
  TF1 *prior4 = new TF1("prior4","TMath::Abs(TMath::Sin(x))");
  
  
  for (int i = 0; i < t->GetEntries(); i++)
    {
      t->GetEntry(i);

      if(stp<burnIn) continue;
      if(stp>max) continue;

      if(onehierarchyonly==1){
	if (tdm32 < 0) continue; // select NH only
      }
      if(onehierarchyonly==-1){
	if (tdm32 > 0) continue;
      }
      if(xsec_20<-3.2||xsec_20>-1.6) continue;

      // Valor marginalisation: 
      // nu_theta23: 0.3 - 0.7
      /* if(tnu_th23 > 0.7 || tnu_th23 < 0.3) continue;
      // nu_theta13: 0 - 0.11
      if (tnu_th13 > 0.11) continue;
      // nu_dm32: 2-3 e-3 (NH), -3--2 e-3 (IH)
      if (!((tnu_dm32>2e-3 && tnu_dm32<3e-3) || (tnu_dm32>-3e-3 && tnu_dm32<-2e-3))) continue;
      // theta13: 0 - 0.11
      if (tth13 > 0.11) continue;

      // Valor fit
      // theta23: 0.246875 - 0.753125
      if(tth23 > 0.753125 || tth23 < 0.246875) continue;
      // dm32: 1.765-3 - 3.235e-3 (NH)
      if (!((tdm32>1.765e-3 && tdm32<3.235e-3) || (tdm32>-3e-3 && tdm32<-2e-3))) continue;*/


      // beta (comment out if not doing nuebar appearance!)
      //if (tbeta==0) continue; // accept beta=1 steps only
      //if (tbeta==1) continue; // accept beta=0 steps only

      // Apply different prior to osc pars 
      // Note: this is assuming it is being run on the fit without reactor (so no prior to ‘undo’)
      // Comment out unless needed 

      // Flat in theta23 (not sin^2 theta23) 
      //tbiaswgt = prior1->Eval(tth23);

      // Flat in theta13 (not sin^2 theta13)
      //tbiaswgt = prior1->Eval(tth13);

      // Flat in sin^2 2theta13 (not sin^2 theta13) 
      //tbiaswgt = prior2->Eval(tth13);

      // Flat in sin dcp (not dcp)
      //tbiaswgt = prior3->Eval(tdcp);  

      // Flat in cos dcp (not dcp)
      //tbiaswgt = prior4->Eval(tdcp);  

      //if (tbiaswgt<100){
	// Fill main histogram (this one to be used to make contours and plot)
	th13_dcp->Fill(tth23, tdm32, RCrw);

      // Fill histograms to tell you the posterior probability by octant/MH
	if (tdm32>0) // if NH
	  {
	    if (tth23>0.5) {highoct_NH->Fill(tth23,tdm32,RCrw);} // high octant, NH
	    else { lowoct_NH->Fill(tth23,tdm32,RCrw);} // low octant, NH
	  }
	else // if IH
	  {
	    if (tth23>0.5) {highoct_IH->Fill(tth23,tdm32,RCrw);} // high octant, IH
	    else { lowoct_IH->Fill(tth23,tdm32,RCrw);} // low octant, IH
	  }
      //}
    }

  //th13_dcp->Rebin2D(1,5);
  //th13_dcp->Smooth(1);

  TCanvas *c = new TCanvas("c","",400/**2*/,400);

  
  //c->Divide(2);
  //c->cd(1);
  gPad->SetLogz();
  gPad->SetLeftMargin(0.13);
  if(doRCrw && onehierarchyonly==1)th13_dcp->SetTitle("With Reactor Constraint, NH;sin^{2}#theta_{23};#Delta m^{2}_{32} (eV^{2})");
  if(doRCrw && onehierarchyonly==-1)th13_dcp->SetTitle("With Reactor Constraint, IH;sin^{2}#theta_{23};#Delta m^{2}_{32} (eV^{2})");
  if(doRCrw && onehierarchyonly==0)th13_dcp->SetTitle("With Reactor Constraint, Both Hierarchies;sin^{2}#theta_{23};#Delta m^{2}_{32} (eV^{2})");
  if(!doRCrw && onehierarchyonly==1)th13_dcp->SetTitle("Without Reactor Constraint, NH;sin^{2}#theta_{23};#Delta m^{2}_{32} (eV^{2})");
  if(!doRCrw && onehierarchyonly==-1)th13_dcp->SetTitle("Without Reactor Constraint, IH;sin^{2}#theta_{23};#Delta m^{2}_{32} (eV^{2})");
  if(!doRCrw && onehierarchyonly==0)th13_dcp->SetTitle("Without Reactor Constraint, Both Hierarchies;sin^{2}#theta_{23};#Delta m^{2}_{32} (eV^{2})");
  th13_dcp->GetYaxis()->SetTitleOffset(1);
  th13_dcp->GetYaxis()->SetTitleSize(0.04);
  //th13_dcp->GetXaxis()->SetTitleOffset(1.03);
  th13_dcp->GetXaxis()->SetTitleSize(0.04);
  //th13_dcp->Draw("COLZ");

  // make contours
  double p68, p90, p95, p99, p3sig;
  double *p = getInterval2D(th13_dcp, p68, p90, p95, p99, p3sig);
  TH2D* a68th23 = (TH2D*)th13_dcp->Clone("cl_68");
  TH2D* b90th23 = (TH2D*)th13_dcp->Clone("cl_90");
  TH2D* b95th23 = (TH2D*)th13_dcp->Clone("cl_95");
  TH2D* b99th23 = (TH2D*)th13_dcp->Clone("cl_99");
  TH2D* b3sigth23 = (TH2D*)th13_dcp->Clone("cl_3sig");

  double tpp68[1];
  double tpp90[1];
  double tpp95[1];
  double tpp99[1];
  double tpp3sig[1];

  tpp68[0] = p68;
  tpp90[0] = p90;
  tpp95[0] = p95;
  tpp99[0] = p99;
  tpp3sig[0] = p3sig;

  th13_dcp->SetContour(255);

  a68th23->Smooth(1);
  b90th23->Smooth(1);
  b95th23->Smooth(1);
  b99th23->Smooth(1);
  b3sigth23->Smooth(1);

  a68th23->SetContour(1, tpp68);
  b90th23->SetContour(1, tpp90);
  b95th23->SetContour(1, tpp95);
  b99th23->SetContour(1, tpp99);
  b3sigth23->SetContour(1, tpp3sig);

  a68th23->SetLineColor(colour_68); //colour_68);
  b90th23->SetLineColor(colour_90); //colour_90);
  b95th23->SetLineColor(colour_95); //colour_90);
  b99th23->SetLineColor(colour_99); //colour_68);
  b3sigth23->SetLineColor(colour_3sig); //colour_90);

  a68th23->SetLineWidth(2);
  a68th23->SetLineStyle(2);

  b90th23->SetLineWidth(2);
  //b90th23->SetLineStyle(2);

  b95th23->SetLineWidth(2);
  //b95th23->SetLineStyle(2);

  b99th23->SetLineWidth(2);
  b99th23->SetLineStyle(7);

  b3sigth23->SetLineWidth(2);
  b3sigth23->SetLineStyle(7);

  // extract contours
  TGraph *gr_68, *gr_90, *gr_95, *gr_99,*gr_3sig;

  TH2D *hconts = (TH2D*)th13_dcp->Clone("hconts");
  //hconts->Smooth(3);

  // Save contours and best fit to file
  TFile *contours_out = new TFile("contours_th23dm23.root","RECREATE");

  double contlev[2];
  contlev[0] = p90;
  contlev[1] = p68;

/*
  contlev[2] = p95;
  contlev[3] = p99;
  contlev[4] = p3sig;
*/

  hconts->SetContour(2,contlev);
  hconts->Draw("CONT LIST");//!!
  c->Update();
  c->Clear();//!!
  TH2D* dummyhist=(TH2D*)hconts->Clone("dummyhist");//
  for(unsigned xbin=0;xbin<=dummyhist->GetNbinsX();xbin++){
    for(unsigned ybin=0;ybin<=dummyhist->GetNbinsY();ybin++){
      dummyhist->SetBinContent(xbin,ybin,0);
    }
  }
  dummyhist->Draw();
  gPad->Update();
  TObjArray *contours = (TObjArray*) gROOT->GetListOfSpecials()->FindObject("contours");
  TList *list68 = (TList*)contours->At(1);
  TList *list90 = (TList*)contours->At(0);

/*
  TList *list95 = (TList*)contours->At(2);
  TList *list99 = (TList*)contours->At(3);
  TList *list3sig = (TList*)contours->At(4);
*/

  TIter next(list90);
  int it=0;
  TGraph *g;
  while ((g = (TGraph*)next()))
    {
      char name[10];
      sprintf(name, "g90_%d",it);
      g->SetLineWidth(2);
      g->SetLineColor(kBlue);
      g->Write(name);
      it++;
    }
  TIter next2(list68);
  it = 0;
  while ((g = (TGraph*)next2()))
    {
      char name[10];
      sprintf(name, "g68_%d",it);
      g->SetLineWidth(2);
      g->SetLineColor(kBlue);
      g->SetLineStyle(2);
      g->Write(name);
      it++;
    }

/*  
  TIter next3(list95);
  it = 0;
  while ((g = (TGraph*)next3()))
    {
      char name[10];
      sprintf(name, "g95_%d",it);
      g->SetLineWidth(2);
      g->SetLineColor(kBlue);
      //g->SetLineStyle(2);
      g->Write(name);
      it++;
    }

  
  TIter next4(list99);
  it = 0;
  while ((g = (TGraph*)next4()))
    {
      char name[10];
      sprintf(name, "g99_%d",it);
      g->SetLineWidth(2);
      g->SetLineColor(kBlue);
      g->SetLineStyle(7);
      g->Write(name);
      it++;
    }

  
  TIter next5(list3sig);
  it = 0;
  while ((g = (TGraph*)next5()))
    {
      char name[10];
      sprintf(name, "g3sig_%d",it);
      g->SetLineWidth(2);
      g->SetLineColor(kBlue);
      g->SetLineStyle(7);
      g->Write(name);
      it++;
    }
*/

  // back to plotting things
  //th13_dcp->Draw("colz");
  a68th23->Draw("same CONT3");
  b90th23->Draw("SAME CONT3");
  b99th23->Draw("SAME CONT3");

  TH2D *cl_a68 = (TH2D*)a68th23->Clone();
  TH2D *cl_b90 = (TH2D*)b90th23->Clone();
  TH2D *cl_b99 = (TH2D*)b99th23->Clone();
  cl_a68->SetLineColor(kBlack);
  cl_b90->SetLineColor(kBlack);
  cl_b99->SetLineColor(kBlack);

  // Make Legend
  TLegend *leg1;
  if(onehierarchyonly!=0) leg1 = new TLegend(0.55,0.68,0.94,0.87);
  else leg1 = new TLegend(0.55,0.4,0.94,0.6);
  leg1->SetLineColor(0);
  leg1->SetTextFont(132);
  leg1->AddEntry(cl_a68,"68% credible interval","l");
  leg1->AddEntry(cl_b90,"90% credible interval","l");
  leg1->AddEntry(cl_b99,"99% credible interval","l");

  // find best fit point
  int Mbx, Mby, Mbz, MaCh3bin;
  MaCh3bin = th13_dcp->GetMaximumBin();
  th13_dcp->GetBinXYZ(MaCh3bin,Mbx,Mby,Mbz);
  double Mx, My;
  Mx = th13_dcp->GetXaxis()->GetBinCenter(Mbx);
  My = th13_dcp->GetYaxis()->GetBinCenter(Mby);
  TGraph *bestfitM = new TGraph(1);
  bestfitM->SetPoint(0,Mx,My);
  bestfitM->SetMarkerStyle(22);
  bestfitM->SetMarkerSize(1);
  bestfitM->SetMarkerColor(kBlack);
  bestfitM->Draw("SAME.P");

  leg1->AddEntry(bestfitM,"MaCh3 best fit","p");
  leg1->Draw("SAME");

  TCanvas *c1 = new TCanvas("c1","",400,400);
  gPad->Update();

  dummyhist->Draw();
  gPad->Update(); 

  a68th23->Draw("same CONT3");
  b95th23->Draw("SAME CONT3");
  b3sigth23->Draw("SAME CONT3");
  bestfitM->Draw("SAME.P");

  TH2D *cl_b95 = (TH2D*)b95th23->Clone();
  TH2D *cl_b3sig = (TH2D*)b3sigth23->Clone();
  cl_b95->SetLineColor(kBlack);
  cl_b3sig->SetLineColor(kBlack);

  // Make Legend
  TLegend *leg2;
  if(onehierarchyonly!=0) leg2 = new TLegend(0.55,0.68,0.94,0.87);
  else leg2 = new TLegend(0.55,0.4,0.94,0.6);
  leg2->SetLineColor(0);
  leg2->SetTextFont(132);
  leg2->AddEntry(cl_a68,"1#sigma credible interval","l");
  leg2->AddEntry(cl_b95,"2#sigma credible interval","l");
  leg2->AddEntry(cl_b3sig,"3#sigma credible interval","l");
  leg2->AddEntry(bestfitM,"MaCh3 best fit","p");
  leg2->Draw("SAME");


  std::cout << " -------------------- 2D Best-fit point ------------------- " << std::endl;
  std::cout << "MaCh3 best fit point: sin^2 theta_23 = " << Mx << ", Delta m^2_32 = " << My << endl;
  std::cout << "Posterior value: " << th13_dcp->GetBinContent(MaCh3bin) << std::endl << std::endl;

  // make graph to plot true values used to generate the fake data
  TGraph *trueg = new TGraph(1);
  if (!(truex == 0 && truey == 0))
    {
      trueg->SetPoint(0,truex,truey);
      std::cout << "True parameter values: sin^2 theta_23 = " << truex << ", Delta m^2_32 = " << truey << std::endl;
      std::cout << "Posterior value: " << th13_dcp->GetBinContent(th13_dcp->FindBin(truex,truey)) << std::endl;
      std::cout << "Posterior (Best fit) / Posterior (true): " << th13_dcp->GetBinContent(MaCh3bin)/th13_dcp->GetBinContent(th13_dcp->FindBin(truex,truey)) << std::endl;
      trueg->SetMarkerStyle(34);
      trueg->SetMarkerSize(1);
      trueg->SetMarkerColor(kBlue);
    }


  // ------- CONFIDENCE INTERVALS ------- //
  std::cout<<"Getting confidence intervals:"<<std::endl;
  double maxb = th13_dcp->GetMaximum();
  std::cout<<"Max bin is: "<<maxb<<std::endl;
  for (int i = 1; i <= th13_dcp->GetNbinsX(); i++) 
    {                                                                
      for (int j = 1; j <= th13_dcp->GetNbinsY(); j++)                
        { 
	  double bc = th13_dcp->GetBinContent(i,j);
	  if (bc == 0) bc = 1;
	  double l = -2 * TMath::Log(bc/maxb);
	  llh_th13_dcp->SetBinContent(i,j, l);
	}
    }



  //c->cd(2);
  TCanvas *cc1 = new TCanvas("cc1","",400,400);
  gPad->Update();

  // make contours                                                                               
  
  double minchi2 = llh_th13_dcp->GetMinimum();
  std::cout << minchi2 << std::endl;
  double param68[1] = {minchi2+2.30};
  double param90[1] = {(minchi2+4.61)};
  TH2D *cllh_th13_dcp = (TH2D*)llh_th13_dcp->Clone("cloned");
  TH2D *cllh_th13_dcp90 = (TH2D*)llh_th13_dcp->Clone("cloned90");
  cllh_th13_dcp->Smooth(1);
  cllh_th13_dcp90->Smooth(1);
  cllh_th13_dcp->SetContour(1, param68);
  cllh_th13_dcp90->SetContour(1, param90);


  cllh_th13_dcp->SetTitle("Confidence Intervals;sin^{2}#theta_{23};#Delta m^{2}_{32}");
  gPad->SetLeftMargin(0.13);
  cllh_th13_dcp->GetYaxis()->SetTitleOffset(1.3);
  cllh_th13_dcp->GetYaxis()->SetTitleSize(0.06);
  //cllh_th13_dcp->GetXaxis()->SetTitleOffset(1.03);
  cllh_th13_dcp->GetXaxis()->SetTitleSize(0.06);

  //valor->SetLineWidth(1);

  // Save confidence intervals too

  double contlev2[2] = {(minchi2+2.30), (minchi2+4.61)};;
  llh_th13_dcp->SetContour(2,contlev2);
  llh_th13_dcp->SetLineColor(kBlue);
  llh_th13_dcp->SetLineWidth(2);
  llh_th13_dcp->Draw("CONT LIST");//!!
  cc1->Update();
  gPad->Update();
  TObjArray *contours_conf = (TObjArray*) gROOT->GetListOfSpecials()->FindObject("contours");
  TList *list68_conf = (TList*)contours_conf->At(0);
  TList *list90_conf = (TList*)contours_conf->At(1);
  TIter next_conf(list90_conf);
  int it_conf=0;
  while ((g = (TGraph*)next_conf()))
    {
      char name[10];
      sprintf(name, "g90_conf_%d",it_conf);
      g->SetLineColor(kBlue);
      g->SetLineWidth(2);
      g->Write(name);
      it_conf++;
    }
  TIter next2_conf(list68_conf);
  it_conf = 0;
  while ((g = (TGraph*)next2_conf()))
    {
      char name[10];
      sprintf(name, "g68_conf_%d",it_conf);
      g->SetLineColor(kBlue);
      g->SetLineWidth(2);
      g->SetLineStyle(2);
      g->Write(name);
      it_conf++;
    }

  // Now back to making plots

  cllh_th13_dcp->SetLineColor(kBlue);
  cllh_th13_dcp->SetLineWidth(2);
  cllh_th13_dcp->SetLineStyle(2);
  //cllh_th13_dcp->Draw("CONT3.LIST");//!!
  cllh_th13_dcp90->SetLineColor(kBlue);
  cllh_th13_dcp90->SetLineWidth(2);
  //cllh_th13_dcp90->Draw("CONT3.SAME.LIST");//!!
  bestfitM->SetMarkerColor(kBlue);
  bestfitM->Draw("SAME.P");
  if (!(truex == 0 && truey == 0)) trueg->Draw("SAME.P");


  // Make Legend
  // Fake graphs for VaLOR and p-theta in legend
/*
  TGraph *fakev = new TGraph(1);
  TGraph *fakev2 = new TGraph(1);
  fakev->SetLineWidth(2);
  fakev->SetLineColor(kRed);
  fakev2->SetMarkerStyle(23);
  fakev2->SetMarkerColor(kRed);
*/
  TGraph *fakep = new TGraph(1);
  TGraph *fakep2 = new TGraph(1);
  fakep->SetLineWidth(2);
  fakep->SetLineColor(kGreen+2);
  fakep2->SetMarkerStyle(22);
  fakep2->SetMarkerColor(kGreen+2);

  TLegend *leg = new TLegend(0.6,0.8,0.93,0.93);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(132);
  if (truex == 0 && truey == 0)
    {
      // do nothing 
    }
  else
    {
      leg->AddEntry(trueg,"Input values","p");
    }
  leg->AddEntry(bestfitM,"MaCh3 best fit","p");
  leg->AddEntry(cllh_th13_dcp,"MaCh3 confidence interval","l");
  //leg->AddEntry(fakev2,"VaLOR best fit","p");
  //leg->AddEntry(fakev,"VaLOR confidence interval","l");
  leg->AddEntry(fakep2,"P-Theta best fit","p");
  leg->AddEntry(fakep,"P-Theta confidence interval","l");
  leg->Draw("SAME");

/*
 *
  c->cd();
  TGraph *bestfitM2 = (TGraph*)bestfitM->Clone();
  bestfitM2->SetMarkerColor(kBlack); //kWhite);
  bestfitM2->Draw("SAME.P");
  if (truex == 0 && truey == 0)
    {
      // do nothing
    }
  else
    trueg->Draw("SAME.P");

  leg1->AddEntry(bestfitM2,"MaCh3 best fit","p");
  if (truex ==0 && truey == 0)
    {
      // do nothing
    }
  else
    leg1->AddEntry(trueg,"Input values","p");

  leg1->Draw("SAME");

*/

  contours_out->cd();
  bestfitM->Write("Bestfit");

  c->SaveAs("contours_th23dm23.pdf");
  c->SaveAs("contours_th23dm23forroot.root");
  c1->SaveAs("contours_th23dm23_sigmas.pdf");
  c1->SaveAs("contours_th23dm23_sigmas.root");
  cc1->SaveAs("contours_th23dm23_conf.pdf");
  cc1->SaveAs("contours_th23dm23_conf.root");


  // Print out posterior probability by octant/hierarchy to the screen

  // First normalise histograms
  double total = th13_dcp->Integral();
  lowoct_NH->Scale(1.0/total);
  highoct_NH->Scale(1.0/total);
  lowoct_IH->Scale(1.0/total);
  highoct_IH->Scale(1.0/total);

  std::cout.precision(3);

  std::cout << " ------------ Posterior probability by octant/mass hierarchy --------------- " << std::endl
	    << "Low octant, IH: " << lowoct_IH->Integral() << std::endl
	    << "High octant, IH: " << highoct_IH->Integral() << std::endl
	    << "Low octant, NH: " << lowoct_NH->Integral() << std::endl
	    << "High octant, NH: " << highoct_NH->Integral() << std::endl
	    << std::endl
	    << "Total low octant: " << lowoct_IH->Integral() + lowoct_NH->Integral() << std::endl
	    << "Total high octant: " << highoct_IH->Integral() + highoct_NH->Integral() << std::endl
	    << std::endl
	    << "Total NH: " << lowoct_NH->Integral() + highoct_NH->Integral() << std::endl
	    << "Total IH: " << lowoct_IH->Integral() + highoct_IH->Integral() << std::endl
	    << std::endl;

  th13_dcp->Write("hist2D");
  lowoct_NH->Write("lowoct_NH");
  highoct_NH->Write("highoct_NH");
  lowoct_IH->Write("lowoct_IH");
  highoct_IH->Write("highoct_IH");

}


/*
void llh_ratio()
{
  contours_moreIntervals("antinu_withsysts_FDSA_31Jan.root", true);
  //llh_ratio("rdf_0_reactor_reduced.root", true);
}


// Version for use with a histogram
void llh_ratio(TH2D *hh, double truex, double truey, int max_steps=-999)
{
  gStyle->SetPalette(51,0);
  gStyle->SetOptStat(0);

  TH2D *th13_dcp = (TH2D*)hh->Clone(); // this is secretly th23 dm32
  // closest bin centre to true values (0.524,2.51e-3) is at (0.525,2.50e-3)


  double min = th13_dcp->GetMinimum();
  if (min<0)
    {
      std::cout << "min = " << min << std::endl;
      for (int i=0; i<th13_dcp->GetXaxis()->GetNbins(); i++)
	{
	  for (int j=0; j<th13_dcp->GetYaxis()->GetNbins(); j++)
	    {
	      double val = th13_dcp->GetBinContent(i+1,j+1);
	      th13_dcp->SetBinContent(i+1,j+1,val-min);
	    }
	}
    }


  TH2D *llh_th13_dcp = (TH2D*)th13_dcp->Clone();

    TCanvas *c = new TCanvas("c","",400,400);
  
  //c->Divide(2);
  //c->cd(1);
  //gPad->SetLogz();
  gPad->SetLeftMargin(0.13);
  th13_dcp->SetTitle(";sin^{2}#theta_{23};#Delta m^{2}_{32}");
  th13_dcp->GetYaxis()->SetTitleOffset(1.4);
  th13_dcp->GetYaxis()->SetTitleSize(0.06);
  //th13_dcp->GetXaxis()->SetTitleOffset(1.03);
  th13_dcp->GetXaxis()->SetTitleSize(0.06);
  th13_dcp->Draw("COLZ");//!!

  // make contours
  double p68, p90;
  double *p = getInterval2D(th13_dcp, p68, p90);
  TH2D* a68th23 = (TH2D*)th13_dcp->Clone();
  TH2D* b90th23 = (TH2D*)th13_dcp->Clone();

  double tpp68[1];
  double tpp90[1];

  tpp68[0] = p68;
  tpp90[0] = p90;

  th13_dcp->SetContour(255);

  a68th23->SetContour(1, tpp68);
  b90th23->SetContour(1, tpp90);

  a68th23->SetLineColor(colour_68);
  b90th23->SetLineColor(colour_90);
  a68th23->SetLineWidth(1);
  a68th23->SetLineStyle(2);
  b90th23->SetLineWidth(1);

  a68th23->Draw("SAME.CONT3.LIST");//!!
  b90th23->Draw("SAME.CONT3.LIST");//!!

  TH2D *cl_a68 = (TH2D*)a68th23->Clone();
  TH2D *cl_b90 = (TH2D*)b90th23->Clone();
  cl_a68->SetLineColor(kBlack);
  cl_b90->SetLineColor(kBlack);

  // Make Legend
  TLegend *leg1 = new TLegend(0.6,0.8,0.93,0.93);
  leg1->SetFillColor(kWhite);
  leg1->SetTextFont(132);
  leg1->AddEntry(cl_a68,"68% credible interval","l");
  leg1->AddEntry(cl_b90,"90% credible interval","l");
  //leg1->Draw("SAME");

  double maxb = th13_dcp->GetMaximum();
  for (int i = 1; i <= th13_dcp->GetNbinsX(); i++)                                                                                                                                                                                                                            
    {                                                                                                                                                                                                                                                                          
      for (int j = 1; j <= th13_dcp->GetNbinsY(); j++)                                                                                                                                                                                                                        
        { 
	  double bc = th13_dcp->GetBinContent(i,j);
	  if (bc == 0)
	    bc = 1;
	  double l = -2 * TMath::Log(bc/maxb);
	  llh_th13_dcp->SetBinContent(i,j, l);
	}
    }


  //c->cd(2);
  
  // make contours                                                                               
  
//  double minchi2 = llh_th13_dcp->GetMinimum();
//  std::cout << minchi2 << std::endl;
//  double param[2] = {(minchi2+4.61)};//2.30), (minchi2+4.61)};
//  TH2D *cllh_th13_dcp = (TH2D*)llh_th13_dcp->Clone("cloned");
//  cllh_th13_dcp->SetContour(1, param); // change to 2
//
//  cllh_th13_dcp->SetTitle("Confidence Intervals;sin^{2}#theta_{23};#Delta m^{2}_{32}");
//  gPad->SetLeftMargin(0.13);
//  cllh_th13_dcp->GetYaxis()->SetTitleOffset(1.3);
//  cllh_th13_dcp->GetYaxis()->SetTitleSize(0.06);
//  //cllh_th13_dcp->GetXaxis()->SetTitleOffset(1.03);
//  cllh_th13_dcp->GetXaxis()->SetTitleSize(0.06);

  //valor->SetLineWidth(1);

  // find best fit point
  int Mbx, Mby, Mbz, MaCh3bin;
  MaCh3bin = th13_dcp->GetMaximumBin();
  th13_dcp->GetBinXYZ(MaCh3bin,Mbx,Mby,Mbz);
  double Mx, My;
  Mx = th13_dcp->GetXaxis()->GetBinCenter(Mbx);
  My = th13_dcp->GetYaxis()->GetBinCenter(Mby);
  TGraph *bestfitM = new TGraph(1);
  bestfitM->SetPoint(0,Mx,My);
  bestfitM->SetMarkerStyle(22);
  bestfitM->SetMarkerColor(kBlack);
  std::cout << "MaCh3 best fit point: sin^2 theta_23 = " << Mx << ", Delta m^2_32 = " << My << endl;

  // make graph to plot true values used to generate the fake data
  TGraph *trueg = new TGraph(1);
  trueg->SetPoint(0,truex,truey);
  std::cout << "True parameter values: sin^2 theta_23 = " << truex << ", Delta m^2_32 = " << truey << std::endl;
  trueg->SetMarkerStyle(34);
  trueg->SetMarkerColor(kBlue);

//  cllh_th13_dcp->SetLineColor(1);
//  cllh_th13_dcp->SetLineWidth(2);
//  cllh_th13_dcp->Draw("CONT2");
//  bestfitM->Draw("SAME.P");
//  trueg->Draw("SAME.P");
//
//  // Make Legend
//  TLegend *leg = new TLegend(0.6,0.8,0.93,0.93);
//  leg->SetFillColor(kWhite);
//  leg->SetTextFont(132);
//  leg->AddEntry(trueg,"Input values","p");
//  leg->AddEntry(bestfitM,"MaCh3 best fit","p");
//  leg->AddEntry(cllh_th13_dcp,"MaCh3 90% confidence interval","l");
//  leg->Draw("SAME");

  //c->cd(1);
  TGraph *bestfitM2 = (TGraph*)bestfitM->Clone();
  bestfitM2->SetMarkerColor(kBlack);
  bestfitM2->Draw("SAME.P");
  trueg->Draw("SAME.P");

  leg1->AddEntry(bestfitM2,"MaCh3 best fit");
  leg1->AddEntry(trueg,"Input values","p");

  leg1->Draw("SAME");

}


// -------------------------------------------------
TH1D* dchi2_1d(TH1D *hist) 
{
  TH1D *dchi2 = (TH1D*)hist->Clone("dchi2");
  dchi2->Reset();

  double maxb = hist->GetMaximum();
  for (int i=1; i<hist->GetNbinsX()+1; i++)
    {
      double bc = hist->GetBinContent(i);
      if (bc == 0)
	bc = 1;
      double l = -2 * TMath::Log(bc/maxb);
      dchi2->SetBinContent(i,l);
    }

  return dchi2;
}
*/
