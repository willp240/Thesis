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
    bool drawCorr = false;
    if (argc == 3) {
      drawCorr = true;
    }
    DrawComp(filename, drawCorr);

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
    //    DrawComp(filename, title1, filename2, title2, filename3, title3);

  }

  return 0;
}




// All right, let's re-write this in a faster more efficient way!
void DrawComp(std::string inputFile, bool drawCorr) {

  // drawCorr decideds if we want to draw correlations or not
  std::cout << "File for study:       " << inputFile << std::endl;
  std::cout << "Draw correlations?    " << drawCorr << std::endl;

  // Open the chain
  TChain* chain = new TChain("posteriors","");
  chain->Add(inputFile.c_str());
  //int nEntries = chain->GetEntries();
  // Use first 1/5 as burn-in
  int cut = 200000;// nEntries/4;
  std::stringstream ss;
  ss << "step > " << cut;

  // What entries are we interested in?
  // Here specifically 20k steps in, and events that have a logL < 5000
  std::string stepcut = ss.str();

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
  // Probably take these as user input instead!
  bool plotFlux = false;

  // Loop over the number of branches
  // Find the name and how many of each systematic we have
  chain->SetBranchStatus("*", false);
  chain->SetBranchStatus("step", true);
  for (int i = 0; i < nbr; i++) {

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
      bnames[ndraw]=bname;
      ndraw++;
      nflux++;
    } else if(bname.BeginsWith("xsec_")) {
      chain->SetBranchStatus(bname, true);
      chain->SetBranchAddress(bname, &(nom[i]));
      bnames[ndraw]=bname;
      ndraw++;
      nxsec++;
    }
  }

  TFile *TempFile = new TFile(inputFile.c_str(), "open");
  TempFile->ls();
  TTree* Settings = (TTree*)(TempFile->Get("Settings"));
  Settings->ls();
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

  // Get first entry in chain
  chain->GetEntry(0);
  nom.resize(nbr);

  std::cout << "# useful entries (flux, xsec, fsi): " << ndraw << std::endl;
  std::cout << "************************************************" << std::endl;

  gStyle->SetOptFit(111);

  // Open a TCanvas to write the posterior onto
  TCanvas* c0 = new TCanvas("c0", "c0", 0, 0, 1024, 1024);
  c0->SetGrid();
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  c0->SetTickx();
  c0->SetTicky();
  c0->SetBottomMargin(0.1);
  c0->SetTopMargin(0.04);
  c0->SetRightMargin(0.03);
  c0->SetLeftMargin(0.10);

  // Make sure we can read files located anywhere and strip the .root ending
  inputFile = inputFile.substr(0, inputFile.find(".root"));

  TString canvasname = inputFile;
  // Append if we're drawing correlations
  // Open bracket means we want to make a pdf file
  if (drawCorr) {
    canvasname += "_drawCorr.pdf[";
  } else {
    canvasname += "_drawPar.pdf[";
  }

  // Once the pdf file is open no longer need to bracket
  canvasname.ReplaceAll("[","");

  // We fit with this Gaussian
  TF1 *gauss    = new TF1("gauss","[0]/sqrt(2.0*3.14159)/[2]*TMath::Exp(-0.5*pow(x-[1],2)/[2]/[2])",-5,5);
  gauss->SetLineWidth(2);
  gauss->SetLineColor(kOrange-5);

  // Some TVectors with the means, errors and Gaussian parameters of the PDFs
  TVectorD* mean_vec  = new TVectorD(ndraw);
  TVectorD* err_vec   = new TVectorD(ndraw); 
  TVectorD* gaus_mean_vec = new TVectorD(ndraw);
  TVectorD* gaus_err_vec  = new TVectorD(ndraw); 
  TVectorD* HPD_mean_vec  = new TVectorD(ndraw);
  TVectorD* HPD_err_p_vec   = new TVectorD(ndraw); 
  TVectorD* HPD_err_vec   = new TVectorD(ndraw); 
  TVectorD* HPD_err_m_vec   = new TVectorD(ndraw); 

  TVectorD* xsec_nom  = new TVectorD(nxsec);

  // Only want to draw covariance of xsec parameters!
  TMatrixT<double>* covariance = new TMatrixT<double>(ndraw,ndraw);
  TMatrixT<double>* correlation = new TMatrixT<double>(ndraw,ndraw);
  for (int i = 0; i < ndraw; ++i) {
    for (int j = 0; j < ndraw; ++j) {
      (*covariance)(i,j) = 0.;
      (*correlation)(i,j) = 0.;
    }
  }

  bool isXsec = false;

  // Remember which parameters we're varying
  bool *vary = new bool[ndraw]();
  for (int i=0; i < ndraw; ++i) {
    vary[i] = false;
  }

  // ndraw is number of draws we want to do
  for(int i = 118; i < 122; ++i) {
    isXsec = false;

      int nbins = 70;

      double asimovLine = 1.0;

      std::string tempString = std::string(bnames[i]);
      if (bnames[i].BeginsWith("xsec_")) {
        isXsec = true;

        int paramNo = 0;
        if (plotFlux && i >= nflux) {
          paramNo = i - nflux;
        } 
        tempString = GetXsecName(paramNo);

        double central, prior, down, up;
        GetXsecLimits(paramNo, central, prior, down, up);
        (*xsec_nom)(paramNo) = central;
        asimovLine = prior;
      }

      double maxi = chain->GetMaximum(bnames[i]);
      double mini = chain->GetMinimum(bnames[i]);
      // This holds the posterior density
      TH1D *hpost = new TH1D(bnames[i], bnames[i], nbins, mini, maxi);
      hpost->SetMinimum(0);
      hpost->GetYaxis()->SetTitle("Steps");
      hpost->GetYaxis()->SetNoExponent(false);

      // Project bnames[i] onto hpost, applying stepcut
      chain->Project(bnames[i], bnames[i], stepcut.c_str());

      hpost->Smooth();

      // Get the characteristics of the hpost
      double mean, rms;
      GetArithmetic(hpost, mean, rms);
      double peakval, sigma_p, sigma_m, sigma_hpd;
      GetHPD(hpost, peakval, sigma_hpd, sigma_p, sigma_m);
      double gauss_mean, gauss_rms;
      GetGaussian(hpost, gauss, gauss_mean, gauss_rms);

      std::cout << mean << " +/- " << rms << " (" << peakval << "+/-" << sigma_hpd << " + " << sigma_p << " - " << sigma_m << ")" << " (" << gauss_mean << "+/-" << gauss_rms << ")" << std::endl;

      TLine *hpd = new TLine(peakval, hpost->GetMinimum(), peakval, hpost->GetMaximum());
      hpd->SetLineColor(kBlack);
      hpd->SetLineWidth(2);
      hpd->SetLineStyle(kSolid);

      TLegend *leg = new TLegend(0.12, 0.6, 0.5, 0.95);
      leg->SetTextSize(0.04);
      leg->AddEntry(hpost, Form("#splitline{PDF}{#mu = %.2f, #sigma = %.2f}", hpost->GetMean(), hpost->GetRMS()), "l");
      leg->AddEntry(gauss, Form("#splitline{Gauss}{#mu = %.2f, #sigma = %.2f}", gauss->GetParameter(1), gauss->GetParameter(2)), "l");
      leg->AddEntry(hpd, Form("#splitline{HPD}{#mu = %.2f, #sigma = %.2f (+%.2f-%.2f)}", peakval, sigma_hpd,sigma_p, sigma_m), "l");

      if (isXsec && ((*xsec_nom)(i-nflux)) != 0) {
        mean = mean / ((*xsec_nom)(i-nflux));
        rms = rms / ((*xsec_nom)(i-nflux));
        gauss_mean = gauss_mean / ((*xsec_nom)(i-nflux));
        gauss_rms = gauss_rms / ((*xsec_nom)(i-nflux));
        peakval = peakval / ((*xsec_nom)(i-nflux));
        sigma_hpd = sigma_hpd / ((*xsec_nom)(i-nflux));
        sigma_p = sigma_p / ((*xsec_nom)(i-nflux));
        sigma_m = sigma_m / ((*xsec_nom)(i-nflux));
      } else if (isXsec && ((*xsec_nom)(i-nflux)) == 0) {
        mean = mean + 1.0;
        gauss_mean = gauss_mean + 1.0;
        peakval = peakval + 1.0;
      }

      (*mean_vec)(i)      = mean;
      (*err_vec)(i)       = rms;
      (*gaus_mean_vec)(i) = gauss_mean;
      (*gaus_err_vec)(i)  = gauss_rms;
      (*HPD_mean_vec)(i) = peakval;
      (*HPD_err_p_vec)(i)  = sigma_p;
      (*HPD_err_vec)(i)   = sigma_hpd;
      (*HPD_err_m_vec)(i)  = sigma_m;
      (*covariance)(i,i)  = rms*rms;
      (*correlation)(i,i)  = 1.0;

      hpost->SetLineWidth(2);
      hpost->SetMaximum(hpost->GetMaximum()*1.5);
      hpost->SetTitle(tempString.c_str());
      hpost->GetXaxis()->SetTitle((std::string(hpost->GetTitle())+" rel. nom").c_str());

      // Now make the TLine for the asimov
      TLine *asimov = new TLine(asimovLine, hpost->GetMinimum(), asimovLine, hpost->GetMaximum());
      asimov->SetLineColor(kRed-3);
      asimov->SetLineWidth(2);
      asimov->SetLineStyle(kDashed);
      hpost->GetYaxis()->SetTitleOffset(1.3);
      hpost->Draw();
      hpd->Draw("same");
      asimov->Draw("same");

      // Make the legend
      leg->AddEntry(asimov, Form("#splitline{Input}{x = %.2f}", asimovLine), "l");
      leg->SetLineColor(0);
      leg->SetLineStyle(0);
      leg->SetFillColor(0);
      leg->SetFillStyle(0);
      leg->Draw("same");

      if (hpost->GetMaximum() == hpost->Integral()*1.5) {
        vary[i] = false;
        delete hpost;
        delete asimov;
        delete hpd;
        delete leg;
        continue;
      }
      vary[i] = true;

      // Write to file
      c0->SetName(hpost->GetName());
      c0->SetTitle(hpost->GetTitle());
      c0->Print((tempString+".pdf").c_str());

      delete hpost;
      delete asimov;
      delete hpd;
      delete leg;

  } // End for ndraw

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
