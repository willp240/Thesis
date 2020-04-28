// Compile with something like
//
// g++ `root-config --cflags` -g -o DrawComp DrawComp.cpp -I`root-config --incdir` `root-config --glibs --libs`
//
// Run on ND280 only fits to get xsec and beam posteriors prettified. Essentially a compiled makeTransferMatrixAll
//
// ./DrawComp will tell you how
//
// Should probably prettify this into something managable like classes, or at least some structs

#include <iostream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <algorithm>

#include "TObjArray.h"
#include "TChain.h"
#include "TFile.h"
#include "TBranch.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TLegend.h"
#include "TString.h"
#include "TStyle.h"
#include "TH1.h"
#include "TF1.h"
#include "TH2.h"
#include "TGraphErrors.h"
#include "TVectorD.h"
#include "TColor.h"

// The input cross-section covariance matrix
std::string XsecCov = "";
std::string FluxCov = "";

// Check that xsecCov is set
void CheckXsecCov();

// Helper to put a name on the cross-section plots
std::string GetXsecName(int xsecParamNo);
// Helper to get the limits for each cross-section parameter
void GetXsecLimits(int param, double &central, double &prior, double &down_error, double &up_error);
// Helper to plot the prefit limits
TH1D* MakePrefit(int ndraw);

// Get the highest posterior density numbers from a 1D posterior
void GetHPD(TH1D * const post, double &central, double &error, double &error_pos, double &error_neg);
void GetArithmetic(TH1D * const hpost, double &mean, double &error);
void GetGaussian(TH1D *& hpost, TF1 *& gauss, double &central, double &error);

// The main comparison functions
void DrawComp(std::string inputfile, bool drawcorr = false);
void DrawComp(std::string inputFile1, std::string title1, std::string inputFile2, std::string title2);
void DrawComp(std::string inputFile1, std::string title1, std::string inputFile2, std::string title2, std::string inputFile3, std::string title3);


int main(int argc, char *argv[]) {

  if (argc != 2 && argc != 3 && argc !=5 && argc != 7) {
    std::cerr << "./DrawComp root_file_to_analyse.root drawcorr (root_file_to_compare_to1.root title root_file_to_compare_to2.root title)" << std::endl;
    exit(-1);
  }

  // Have implemented some good old hard-coding for comparing many fits; argc == 2 is 1 fit, argc == 3 is 2 fits, etc
  // Could additionally do fancy stuff; do we want beam+fsi+xsec params, do we want do draw a covariance matrix, etc
  if (argc == 2 || argc == 3) {

    std::cout << "Producing single fit output" << std::endl;
    std::string filename = argv[1];
    //bool drawCorr = false;
    if (argc == 3) {
      //drawCorr = true;
    }
    //    DrawComp(filename, drawCorr);

    // If we want to compare two fits (e.g. binning changes or introducing new params/priors)
  } else if (argc == 3+2*1) {

    std::cout << "Producing two fit comparison" << std::endl;
    std::string filename  = argv[1];
    std::string title1    = argv[2];
    std::string filename2 = argv[3];
    std::string title2    = argv[4];
    //DrawComp(filename, title1, filename2, title2);

    // If we want to compare three fits (e.g. binning changes or introducing new params/priors)
  } else if (argc == 4+3*1) {

    std::cout << "Producing three fit comparison" << std::endl;
    std::string filename = argv[1];
    std::string title1    = argv[2];
    std::string filename2 = argv[3];
    std::string title2    = argv[4];
    std::string filename3 = argv[5];
    std::string title3    = argv[6];
    DrawComp(filename, title1, filename2, title2, filename3, title3);

  }

  return 0;
}


void DrawComp(std::string inputFile1, std::string title1, std::string inputFile2, std::string title2, std::string inputFile3, std::string title3) {

  std::cout << "File for study:       " << inputFile1 << std::endl;
  std::cout << "File for study2:      " << inputFile2 << std::endl;
  std::cout << "File for study3:      " << inputFile3 << std::endl;

  // Open the chain for file 1
  TChain* chain = new TChain("posteriors","");
  chain->Add(inputFile1.c_str());
  // Open the chain for file 2
  TChain* chain2 = new TChain("posteriors","");
  chain2->Add(inputFile2.c_str());
  // Open the chain for file 3
  TChain* chain3 = new TChain("posteriors","");
  chain3->Add(inputFile3.c_str());

  //  int nEntries = chain->GetEntries();
  //int nEntries2 = chain2->GetEntries();
  ///int nEntries3 = chain3->GetEntries();

  // Use first 1/4 as burn-in
  int cut = 50000;
  int cut2 = 50000;
  int cut3 = 50000;
  std::stringstream ss;
  ss << "step > " << cut;
  std::stringstream ss2;
  ss2 << "step > " << cut2;
  std::stringstream ss3;
  ss3 << "step > " << cut3;

  // What entries are we interested in?
  // Here specifically 20k steps in, and events that have a logL < 5000
  std::string stepcut = ss.str();
  std::string stepcut2 = ss2.str();
  std::string stepcut3 = ss3.str();

  // Loop over interesting params; only need to do this for one of the files granted we have same params
  //
  // Get the list of branches
  TObjArray* brlis = (TObjArray*)chain->GetListOfBranches();
  // Get the number of branches
  int nbr = brlis->GetEntries();
  std::cout << "# of branches: " << nbr << std::endl;
  // Make an array of TStrings
  TString bnames[nbr];
  // Array of doubles containing the nominal values (first step for xsec)
  std::vector<double> nom;
  nom.resize(nbr);

  // Have a counter for how many of each systematic we have
  int ndraw = 0;
  int nflux = 0;
  int nxsec = 0;


  // Have bool to remember what parameters we're plotting
  bool plotFlux = false;
  //  bool plotXsec = false;

  // Loop over the number of branches
  // Find the name and how many of each systematic we have
  chain->SetBranchStatus("*", false);
  chain->SetBranchStatus("step", true);
  chain2->SetBranchStatus("*", false);
  chain2->SetBranchStatus("step", true);
  chain3->SetBranchStatus("*", false);
  chain3->SetBranchStatus("step", true);
  for(int i=0; i < nbr; i++) {

    // Get the TBranch and its name
    TBranch* br = (TBranch*)brlis->At(i);
    TString bname = br->GetName();

    // Get first entry of this branch (i.e. systematic); this is the first reconfigure so should be the fake-data which it was generated for (for xsec this is true...!)
    if(bname.BeginsWith("LogL") || bname.BeginsWith("ndd_")) continue;

    // If we're on beam systematics
    if(bname.BeginsWith("b_")) {

      plotFlux = true;
      chain->SetBranchStatus(bname, true);
      chain->SetBranchAddress(bname, &(nom[i]));
      chain2->SetBranchStatus(bname, true);
      chain2->SetBranchAddress(bname, &(nom[i]));
      chain3->SetBranchStatus(bname, true);
      chain3->SetBranchAddress(bname, &(nom[i]));
      bnames[ndraw]=bname;
      nflux++;
      ndraw++;

    } else if(bname.BeginsWith("xsec_")) {

      //      plotXsec = true;
      chain->SetBranchStatus(bname, true);
      chain->SetBranchAddress(bname, &(nom[i]));
      chain2->SetBranchStatus(bname, true);
      chain2->SetBranchAddress(bname, &(nom[i]));
      chain3->SetBranchStatus(bname, true);
      chain3->SetBranchAddress(bname, &(nom[i]));
      bnames[ndraw]=bname;
      ndraw++;
      nxsec++;

    }
  }
  TFile *TempFile = new TFile(inputFile1.c_str(), "open");
  TTree* Settings = (TTree*)(TempFile->Get("Settings"));
  // Get the xsec covariance matrix
  std::string *XsecInput = 0;
  Settings->SetBranchAddress("XsecCov", &XsecInput);
  // And the flux covariance matrix
  std::string *FluxInput = 0;
  Settings->SetBranchAddress("FluxCov", &FluxInput);
  Settings->GetEntry(0);
  delete Settings;
  TempFile->Close();
  delete TempFile;

  XsecCov = "../"+*XsecInput;
  FluxCov = "../"+*FluxInput;
  chain->GetEntry(0);
  nom.resize(nbr);

  std::cout << "# useful entries (flux, xsec, fsi): " << ndraw << std::endl;
  std::cout << "************************************************" << std::endl;

  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);

  // Open a TCanvas to write the posterior onto
  TCanvas* c0 = new TCanvas("c0", "c0", 0, 0, 1024, 1024);
  c0->SetGrid();
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  c0->SetBottomMargin(0.1);
  c0->SetTopMargin(0.02);
  c0->SetRightMargin(0.03);
  c0->SetLeftMargin(0.12);


  // Write to a PDF file
  // Strip out the initial ../ if can find
  // Then strip out any additional /
  // Strip out the last .root
  inputFile1 = inputFile1.substr(0, inputFile1.find(".root"));
  std::cout << inputFile1 << std::endl;
  while (inputFile2.find("/") != std::string::npos) {
    inputFile2 = inputFile2.substr(inputFile2.find("/")+1, inputFile2.find(".root")-1);
  }
  std::cout << inputFile2 << std::endl;
  while (inputFile3.find("/") != std::string::npos) {
    inputFile3 = inputFile3.substr(inputFile3.find("/")+1, inputFile3.find(".root")-1);
  }
  std::cout << inputFile3 << std::endl;
  TString canvasname = inputFile1+"_"+inputFile2+"_"+inputFile3+".pdf[";
  // c0->Print(canvasname);
  // Once the pdf file is open no longer need to bracket
  canvasname.ReplaceAll("[","");

  // We fit with this gaussian
  TF1 *gauss = new TF1("gauss","[0]/sqrt(2.0*3.14159)/[2]*TMath::Exp(-0.5*pow(x-[1],2)/[2]/[2])",   -5, 5);
  TF1 *gauss2 = new TF1("gauss2","[0]/sqrt(2.0*3.14159)/[2]*TMath::Exp(-0.5*pow(x-[1],2)/[2]/[2])", -5, 5);
  TF1 *gauss3 = new TF1("gauss3","[0]/sqrt(2.0*3.14159)/[2]*TMath::Exp(-0.5*pow(x-[1],2)/[2]/[2])", -5, 5);

  TVectorD* mean_vec = new TVectorD(ndraw);
  TVectorD* err_vec = new TVectorD(ndraw); 
  TVectorD* gaus_mean_vec = new TVectorD(ndraw);
  TVectorD* gaus_err_vec = new TVectorD(ndraw); 

  TVectorD* mean_vec2 = new TVectorD(ndraw);
  TVectorD* err_vec2 = new TVectorD(ndraw); 
  TVectorD* gaus_mean_vec2 = new TVectorD(ndraw);
  TVectorD* gaus_err_vec2 = new TVectorD(ndraw); 

  TVectorD* mean_vec3 = new TVectorD(ndraw);
  TVectorD* err_vec3 = new TVectorD(ndraw); 
  TVectorD* gaus_mean_vec3 = new TVectorD(ndraw);
  TVectorD* gaus_err_vec3 = new TVectorD(ndraw); 

  TVectorD* xsec_nom  = new TVectorD(nxsec);

  // Output file to write to
  TString rootfilename = inputFile1+"_"+inputFile2+"_"+inputFile3+".root";
  // TFile* file = new TFile(rootfilename, "RECREATE");
  //file->cd();
  //TDirectory *params = file->mkdir("hpost");

  int nbins = 70;
  bool isXsec = false;

  // ndraw is number of draws we want to do
  for(int i=118; i<122; ++i) {

    //    file->cd();

    double asimovLine = 0.0;

    std::string tempString;

    // Make a line with the asimov on for xsec params
    if (bnames[i].BeginsWith("xsec_")) {
      asimovLine = nom.at(i);
    } else {
      asimovLine = 1.0;
    }


    // Set up slightly different binning
    // If beam systematic
    if (bnames[i].BeginsWith("b_")) {
      tempString = bnames[i];
      // Also if we find xsec systematic make a string WITH THE ACTUAL PARAMETER NAME FOR 2016a WOOP
    } else if (bnames[i].BeginsWith("xsec_")) {
      isXsec = true;

      int paramNo = 0;
      if (plotFlux) {
        paramNo = i - nflux;
      } 
      tempString = GetXsecName(paramNo);

      double central;
      double prior;
      double down;
      double up;
      GetXsecLimits(paramNo, central, prior, down, up);
      (*xsec_nom)(paramNo) = central;
    }
    //file->cd();
    std::cout << bnames[i] << std::endl;

    double maxi = chain->GetMaximum(bnames[i]);
    double mini = chain->GetMinimum(bnames[i]);
    double maxi2 = chain2->GetMaximum(bnames[i]);
    double mini2 = chain2->GetMinimum(bnames[i]);
    double maxi3 = chain3->GetMaximum(bnames[i]);
    double mini3 = chain3->GetMinimum(bnames[i]);

    // This holds the posterior density
    TH1D *hpost = new TH1D(bnames[i], bnames[i], nbins, mini, maxi);
    hpost->SetMinimum(0);
    hpost->GetYaxis()->SetTitle("Steps (area norm.)");
    hpost->GetYaxis()->SetTitleOffset(1.2);
    hpost->GetYaxis()->SetNoExponent(false);
    TH1D *hpost2 = new TH1D(bnames[i]+"_2", bnames[i]+"_2", nbins, mini2, maxi2);
    hpost2->SetMinimum(0);
    TH1D *hpost3 = new TH1D(bnames[i]+"_3", bnames[i]+"_3", nbins, mini3, maxi3);
    hpost3->SetMinimum(0);

    hpost->SetTitle(bnames[i]);
    hpost2->SetTitle(bnames[i]+"_2");
    hpost3->SetTitle(bnames[i]+"_3");

    // Project bnames[i] onto hpost, applying stepcut
    chain->Project(bnames[i], bnames[i], stepcut.c_str());
    chain2->Project(bnames[i]+"_2", bnames[i], stepcut2.c_str());
    chain3->Project(bnames[i]+"_3", bnames[i], stepcut3.c_str());

    // Area normalise the distributions
    hpost->Scale(1./hpost->Integral(), "width");
    hpost2->Scale(1./hpost2->Integral(), "width");
    hpost3->Scale(1./hpost3->Integral(), "width");

    hpost->Smooth();
    hpost2->Smooth();
    hpost3->Smooth();

    hpost->SetTitle(tempString.c_str());
    hpost->GetXaxis()->SetTitle(tempString.c_str());
    hpost2->SetTitle(tempString.c_str());
    hpost2->GetXaxis()->SetTitle(tempString.c_str());
    hpost3->SetTitle(tempString.c_str());
    hpost3->GetXaxis()->SetTitle(tempString.c_str());

    // Get the characteristics of the hpost
    double mean = hpost->GetMean();
    double rms = hpost->GetRMS();
    double peakval = hpost->GetBinCenter(hpost->GetMaximumBin());
    double gauss_mean = gauss->GetParameter(1);
    double gauss_rms = gauss->GetParameter(2);

    double mean2 = hpost2->GetMean();
    double rms2 = hpost2->GetRMS();
    double peakval2 = hpost2->GetBinCenter(hpost2->GetMaximumBin());
    double gauss_mean2 = gauss2->GetParameter(1);
    double gauss_rms2 = gauss2->GetParameter(2);

    double mean3 = hpost3->GetMean();
    double rms3 = hpost3->GetRMS();
    double peakval3 = hpost3->GetBinCenter(hpost3->GetMaximumBin());
    double gauss_mean3 = gauss3->GetParameter(1);
    double gauss_rms3 = gauss3->GetParameter(2);


    double sigma_p, sigma_m, sigma_hpd;
    GetHPD(hpost, peakval, sigma_hpd, sigma_p, sigma_m);
    double sigma_p2, sigma_m2, sigma_hpd2;
    GetHPD(hpost2, peakval2, sigma_hpd2, sigma_p2, sigma_m2);
    double sigma_p3, sigma_m3, sigma_hpd3;
    GetHPD(hpost3, peakval3, sigma_hpd3, sigma_p3, sigma_m3);

    // Set the range for the Gaussian fit
    gauss->SetRange(mean - 1.5*rms , mean + 1.5*rms);
    gauss2->SetRange(mean2 - 1.5*rms2 , mean2 + 1.5*rms2);
    gauss3->SetRange(mean3 - 1.5*rms3 , mean3 + 1.5*rms3);

    // Set the starting parameters close to RMS and peaks of the histograms
    gauss->SetParameters(hpost->GetMaximum()*rms*sqrt(2*3.14), peakval, rms);
    gauss2->SetParameters(hpost2->GetMaximum()*rms2*sqrt(2*3.14), peakval2, rms2);
    gauss3->SetParameters(hpost3->GetMaximum()*rms3*sqrt(2*3.14), peakval3, rms3);

    // Set some nice colours
    hpost->SetLineColor(kBlue);
    hpost->SetLineWidth(2);
    gauss->SetLineColor(kBlue);
    gauss->SetLineStyle(kDashed);

    hpost2->SetLineColor(kBlack);
    hpost2->SetLineWidth(2);
    gauss2->SetLineColor(kBlack);
    gauss2->SetLineStyle(kDashed);

    hpost3->SetLineColor(kGreen+2);
    hpost3->SetLineWidth(2);
    gauss3->SetLineColor(kGreen+2);
    gauss3->SetLineStyle(kDashed);

    // Perform the fit
    hpost->Fit("gauss","Rq");
    hpost2->Fit("gauss2","Rq");
    hpost3->Fit("gauss3","Rq");

    // Now make the TLine for the asimov
    TLine *asimov = new TLine(asimovLine, hpost->GetMinimum(), asimovLine, 1.5*hpost->GetMaximum());
    asimov->SetLineColor(kRed);
    asimov->SetLineWidth(2);
    asimov->SetLineStyle(kDashed);

    // Make a nice little TLegend
    TLegend *leg = new TLegend(0.13, 0.7, 0.8, 0.97);
    leg->SetTextSize(0.03);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->SetLineColor(0);
    leg->SetLineStyle(0);
    TString nombinLeg = Form("#splitline{%s}{HPD = %.2f, #sigma = %.2f}", title1.c_str(), peakval, sigma_hpd);
    TString rebinLeg = Form("#splitline{%s}{HPD = %.2f, #sigma = %.2f}", title2.c_str(), peakval2, sigma_hpd2);
    TString rebinLeg2 = Form("#splitline{%s}{HPD = %.2f, #sigma = %.2f}", title3.c_str(), peakval3, sigma_hpd3);
    TString asimovLeg = Form("Prior, #mu = %.2f", asimovLine);
    leg->AddEntry(asimov, asimovLeg, "l");
    leg->AddEntry(hpost, nombinLeg, "l");
    leg->AddEntry(hpost2, rebinLeg, "l");
    leg->AddEntry(hpost3, rebinLeg2, "l");

    TLine *hpd = new TLine(peakval, hpost->GetMinimum(), peakval, hpost->GetMaximum());
    hpd->SetLineColor(hpost->GetLineColor());
    hpd->SetLineWidth(2);
    hpd->SetLineStyle(kSolid);

    TLine *hpd2 = new TLine(peakval2, hpost2->GetMinimum(), peakval2, hpost2->GetMaximum());
    hpd2->SetLineColor(hpost2->GetLineColor());
    hpd2->SetLineWidth(2);
    hpd2->SetLineStyle(kSolid);

    TLine *hpd3 = new TLine(peakval3, hpost3->GetMinimum(), peakval3, hpost3->GetMaximum());
    hpd3->SetLineColor(hpost3->GetLineColor());
    hpd3->SetLineWidth(2);
    hpd3->SetLineStyle(kSolid);
    
    if (isXsec && ((*xsec_nom)(i-nflux)) != 0) {
      mean = mean / ((*xsec_nom)(i-nflux));
      rms = rms / ((*xsec_nom)(i-nflux));
      gauss_mean = gauss_mean / ((*xsec_nom)(i-nflux));
      gauss_rms = gauss_rms / ((*xsec_nom)(i-nflux));

      mean2 = mean2 / ((*xsec_nom)(i-nflux));
      rms2 = rms2 / ((*xsec_nom)(i-nflux));
      gauss_mean2 = gauss_mean2 / ((*xsec_nom)(i-nflux));
      gauss_rms2 = gauss_rms2 / ((*xsec_nom)(i-nflux));

      mean3 = mean3 / ((*xsec_nom)(i-nflux));
      rms3 = rms3 / ((*xsec_nom)(i-nflux));
      gauss_mean3 = gauss_mean3 / ((*xsec_nom)(i-nflux));
      gauss_rms3 = gauss_rms3 / ((*xsec_nom)(i-nflux));
    } else if (isXsec && ((*xsec_nom)(i-nflux)) == 0) {
      mean = mean + 1.0;
      gauss_mean = gauss_mean + 1.0;
      mean2 = mean2 + 1.0;
      gauss_mean2 = gauss_mean2 + 1.0;
      mean3 = mean3 + 1.0;
      gauss_mean3 = gauss_mean3 + 1.0;
    }

    // Save the mean and RMS for histo and Gaussian fit
    (*mean_vec)(i)      = mean;
    (*err_vec)(i)       = rms;
    (*gaus_mean_vec)(i) = gauss_mean;
    (*gaus_err_vec)(i)  = gauss_rms;

    (*mean_vec2)(i)     = mean2;
    (*err_vec2)(i)      = rms2;
    (*gaus_mean_vec2)(i)= gauss_mean2;
    (*gaus_err_vec2)(i) = gauss_rms2;

    (*mean_vec3)(i)     = mean3;
    (*err_vec3)(i)      = rms3;
    (*gaus_mean_vec3)(i)= gauss_mean3;
    (*gaus_err_vec3)(i) = gauss_rms3;


    // Do the difference between 1 and 2
    double maximum = 0;
    maximum = std::max(hpost->GetMaximum(), hpost3->GetMaximum());
    maximum = std::max(maximum, hpost3->GetMaximum());

    hpost->SetMaximum(1.5*maximum);
    hpost2->SetMaximum(hpost->GetMaximum());
    hpost3->SetMaximum(hpost->GetMaximum());

    //    file->cd();
    //    hpost->GetYaxis()->SetTitle("titleeee");
    hpost->GetYaxis()->SetTitleOffset(1.5);
    hpost->Draw();
    hpost2->Draw("same");
    hpost3->Draw("same");
    asimov->Draw("same");
    hpd->Draw("same");
    hpd2->Draw("same");
    hpd3->Draw("same");
    leg->Draw("same");

    // cd into params directory in root file
    c0->SetName(hpost->GetName());
    c0->SetTitle(hpost->GetTitle());
    c0->Print((tempString+".pdf").c_str());

    delete hpost;
    delete asimov;
    delete leg;
  } // End the for ndraw loop

}


// *****************************
// Make the prefit plots
TH1D* MakePrefit(int ndraw) {
  // *****************************

  // Get the flux covariance matrix
  TFile *FluxFile = new TFile(FluxCov.c_str(), "open");
  FluxFile->cd();
  TMatrixDSym *FluxMatrix = (TMatrixDSym*)(FluxFile->Get("total_flux_cov"));
  int nFlux = FluxMatrix->GetNrows();
  double *FluxPrefitArray = new double[nFlux];
  double *FluxErrorArray = new double[nFlux];
  for (int i = 0; i < nFlux; ++i) {
    // The nominals are equal to the priors are equal to 1 for flux params
    FluxPrefitArray[i] = 1.0;
    FluxErrorArray[i] = sqrt((*FluxMatrix)(i,i));
  }
  FluxFile->Close();
  delete FluxFile;
  delete FluxMatrix;

  // Do the same for the cross-section file
  TFile *XsecFile = new TFile(XsecCov.c_str(), "open");
  XsecFile->cd();
  // Get the matrix
  TMatrixDSym *XsecMatrix = (TMatrixDSym*)(XsecFile->Get("xsec_cov"));
  // And the central priors
  TVectorD *XsecPrior = (TVectorD*)(XsecFile->Get("xsec_param_prior"));
  // And the nominal values
  TVectorD *XsecNominal = (TVectorD*)(XsecFile->Get("xsec_param_nom"));
  TMatrixT<double> *XsecID = (TMatrixD*)(XsecFile->Get("xsec_param_id"));

  // Now make a TH1D of it
  int nXsec = XsecPrior->GetNrows();

  // Now drive in the actually fill
  double* XsecPrefitArray = new double[nXsec];
  // And the errors from the covariance matrix
  double* XsecErrorArray = new double[nXsec];

  for (int i = 0; i < nXsec; ++i) {
    // Normalise the prior relative the nominal: this is just by choice for BANFF comparisons
    if (((*XsecNominal)(i)) != 0) {
      XsecPrefitArray[i] = ((*XsecPrior)(i)) / ((*XsecNominal)(i));
      XsecErrorArray[i] = sqrt((*XsecMatrix)(i,i))/((*XsecNominal)(i));
    } else {
      XsecPrefitArray[i] = ((*XsecPrior)(i)) + 1.0;
      XsecErrorArray[i] = sqrt((*XsecMatrix)(i,i));
    }
  }

  if (ndraw != nFlux+nXsec) {
    std::cerr << "Number of requested drawn parameters are not equal to flux and xsec input" << std::endl;
    std::cerr << "ndraw = " << ndraw << std::endl;
    std::cerr << "nFlux = " << nFlux << "   " << "nXsec = " << nXsec << std::endl;
    std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
    throw;
  }

  TH1D *PreFitPlot = new TH1D("Prefit", "Prefit", ndraw, 0, ndraw);
  for (int i = 0; i < PreFitPlot->GetNbinsX() + 1; ++i) {
    PreFitPlot->SetBinContent(i+1, 0);
    PreFitPlot->SetBinError(i+1, 0);
  }
  // Set the flux values
  for (int i = 0; i < nFlux; ++i) {
    PreFitPlot->SetBinContent(i+1, FluxPrefitArray[i]);
    PreFitPlot->SetBinError(i+1, FluxErrorArray[i]);

    if (i % 5 == 0) PreFitPlot->GetXaxis()->SetBinLabel(i+1, Form("Flux %i", i));
  }
  for (int i = nFlux; i < nXsec+nFlux; ++i) {
    PreFitPlot->SetBinContent(i+1, XsecPrefitArray[i-nFlux]);
    PreFitPlot->SetBinError(i+1, XsecErrorArray[i-nFlux]);

    PreFitPlot->GetXaxis()->SetBinLabel(i+1, GetXsecName(i-nFlux).c_str());
    // Also scale xsec entries relative the error
  }

  PreFitPlot->GetYaxis()->SetTitle("Variation rel. nominal");

  PreFitPlot->SetDirectory(0);

  PreFitPlot->SetFillStyle(1001);
  PreFitPlot->SetFillColor(kRed-3);
  PreFitPlot->SetMarkerStyle(21);
  PreFitPlot->SetMarkerSize(2.7);
  PreFitPlot->SetMarkerColor(kWhite);
  PreFitPlot->SetLineColor(PreFitPlot->GetFillColor());

  PreFitPlot->GetXaxis()->LabelsOption("v");

  delete XsecMatrix;
  delete XsecPrior;
  delete XsecNominal;
  delete XsecID;
  delete[] XsecPrefitArray;
  delete[] XsecErrorArray;

  XsecFile->Close();
  delete XsecFile;

  return PreFitPlot;
}


// **************************
// Function to get limits for a parameter from the input
// Hard-code to +/-2 sigma
void GetXsecLimits(const int param, double &central, double &prior, double &down_error, double &up_error) {
  // **************************

  central = 0.0;
  double error = 0.0;
  double sigmas = 1.0;

  // Open the input covariance file to get the pre-fit error
  TFile *covFile = new TFile(XsecCov.c_str(), "OPEN");
  covFile->cd();
  TMatrixDSym *covMatrix = (TMatrixDSym*)(covFile->Get("xsec_cov"));

  TVectorD* prior_a = (TVectorD*)(covFile->Get("xsec_param_prior"));
  TVectorD* nominal = (TVectorD*)(covFile->Get("xsec_param_nom"));
  TVectorD* lower   = (TVectorD*)(covFile->Get("xsec_param_lb"));
  TVectorD* upper   = (TVectorD*)(covFile->Get("xsec_param_ub"));
  TMatrixT<double> *id = (TMatrixD*)(covFile->Get("xsec_param_id"));

  if ( ((*nominal)(param)) != 0 ) {
    // Set the central to the nominal
    central = ((*nominal)(param));
    prior = ((*prior_a)(param));
    error = sigmas * sqrt(((*covMatrix)(param,param))) / central;
  } else {
    central = ((*nominal)(param));
    prior = ((*prior_a)(param));
    error = sigmas * sqrt(((*covMatrix)(param,param)));
  }

  // We might be passed the valid range of parameter
  // Do a check to see if this is true
  if (central - error <  ((*lower)(param))) {
    down_error = (*lower)(param);
  } else {
    down_error = central - error;
  }

  if (central + error > ((*upper)(param))) {
    up_error = (*upper)(param);
  } else {
    up_error = central + error;
  }

  // Normalisation parameter are unlikely of going above 2 sigma
  if ((*id)(param,0) == -1) {
    up_error = central + 2*sqrt((*covMatrix)(param, param));
  }

  delete covMatrix;
  delete prior_a;
  delete nominal;
  delete id;

  covFile->Close();
  delete covFile;
}



// **************************
// Just a converter from xsec_i to a real parameter name
// Use this when the conversion wasn't written to file :(
std::string GetXsecName(int xsecParamNo) {
  // **************************

  CheckXsecCov();
  // Open the input covariance file to get the pre-fit error
  TFile *covFile = new TFile(XsecCov.c_str(), "OPEN");
  covFile->cd();

  // Get the array of cross-section parameter names
  TObjArray *xsec_param_names = (TObjArray*)(covFile->Get("xsec_param_names"));

  std::string returnString;
  if (xsecParamNo >= xsec_param_names->GetEntries()) {
    std::cerr << "You gave me a cross-section parameter number which was greater than what the xsec covariance matrix has dimensions" << std::endl;
    std::cerr << "Please revise!" << std::endl;
    std::cerr << "xsecParamNo = " << xsecParamNo << std::endl;
    std::cerr << "GetEntries  = " << xsec_param_names->GetEntries() << std::endl;
    returnString = "INVALID";
  } else {
    returnString = std::string(((TObjString*)xsec_param_names->At(xsecParamNo))->GetString());
  }

  covFile->Close();

  return returnString;

} // End function


// **************************
// Function to check the cross-section covariance matrix
void CheckXsecCov() {
  // **************************

  if (XsecCov.empty()) {
    std::cerr << "Have not got a name for the cross-section covariance matrix" << std::endl;
    std::cerr << "Please specify as a global variable at the top of " << __FILE__ << std::endl;
    throw;
  }

}

// Get the highest posterior density from a TH1D
void GetHPD(TH1D * const hpost, double &central, double &error, double &error_pos, double &error_neg) {

  // Get the bin which has the largest posterior density
  int MaxBin = hpost->GetMaximumBin();
  // And it's value
  double peakval = hpost->GetBinCenter(MaxBin);

  // The total integral of the posterior
  double integral = hpost->Integral();

  // Keep count of how much area we're covering
  double sum = 0.0;

  // Counter for current bin
  int CurrBin = MaxBin;
  while (sum/integral < 0.6827/2.0 && CurrBin < hpost->GetNbinsX()+1) {
    sum += hpost->GetBinContent(CurrBin);
    CurrBin++;
  }
  double sigma_p = fabs(hpost->GetBinCenter(MaxBin)-hpost->GetBinCenter(CurrBin));
  // Reset the sum
  sum = 0.0;

  // Reset the bin counter
  CurrBin = MaxBin;
  // Counter for current bin
  while (sum/integral < 0.6827/2.0 && CurrBin >= 0) {
    sum += hpost->GetBinContent(CurrBin);
    CurrBin--;
  }
  double sigma_m = fabs(hpost->GetBinCenter(CurrBin)-hpost->GetBinCenter(MaxBin));

  // Now do the double sided HPD
  sum = 0.0;
  int LowBin = MaxBin-1;
  int HighBin = MaxBin+1;
  double LowCon = 0.0;
  double HighCon = 0.0;
  while (sum/integral < 0.6827 && LowBin >= 0 && HighBin < hpost->GetNbinsX()+1) {

    // Get the slice
    LowCon = hpost->GetBinContent(LowBin);
    HighCon = hpost->GetBinContent(HighBin);

    // If we're on the last slice and the lower contour is larger than the upper
    if ((sum+LowCon+HighCon)/integral > 0.6827 && LowCon > HighCon) {
      sum += LowCon;
      break;
      // If we're on the last slice and the upper contour is larger than the lower
    } else if ((sum+LowCon+HighCon)/integral > 0.6827 && HighCon >= LowCon) {
      sum += HighCon;
      break;
    } else {
      sum += LowCon + HighCon;
    }

    LowBin--;
    HighBin++;
  }

  double sigma_hpd = 0.0;
  if (LowCon > HighCon) {
    sigma_hpd = fabs(hpost->GetBinCenter(LowBin)-hpost->GetBinCenter(MaxBin));
  } else {
    sigma_hpd = fabs(hpost->GetBinCenter(HighBin)-hpost->GetBinCenter(MaxBin));
  }

  central = peakval;
  error = sigma_hpd;
  error_pos = sigma_p;
  error_neg = sigma_m;

}

// **************************
// Get the mean and RMS of a 1D posterior
void GetArithmetic(TH1D * const hpost, double &mean, double &error) {
  // **************************
  mean = hpost->GetMean();
  error = hpost->GetRMS();
}

// **************************
// Get Gaussian characteristics
void GetGaussian(TH1D *& hpost, TF1 *& gauss, double &central, double &error) {
  // **************************

  double mean = hpost->GetMean();
  double err = hpost->GetRMS();
  double peakval = hpost->GetBinCenter(hpost->GetMaximumBin());

  // Set the range for the Gaussian fit
  gauss->SetRange(mean - 1.5*err , mean + 1.5*err);
  // Set the starting parameters close to RMS and peaks of the histograms
  gauss->SetParameters(hpost->GetMaximum()*err*sqrt(2*3.14), peakval, err);

  // Perform the fit
  hpost->Fit(gauss->GetName(),"Rq");
  hpost->SetStats(0);

  central = gauss->GetParameter(1);
  error = gauss->GetParameter(2);

}
