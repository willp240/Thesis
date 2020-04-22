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
#include "TH2Poly.h"
#include "TGraphErrors.h"
#include "TVectorD.h"
#include "TColor.h"

// Get the highest posterior density numbers from a 1D posterior
void GetHPD(TH1D * const post, double &central, double &error, double &error_pos, double &error_neg);
void GetArithmetic(TH1D * const hpost, double &mean, double &error);
int GetTH2DBinNoOverflow(TH2D* histo, int xbin, int ybin );
// The main comparison functions
void extractHPD(std::string inputfile);

int main(int argc, char *argv[]) {

  if (argc != 2) {
    std::cerr << "./DrawComp root_file_to_analyse.root" << std::endl;
    exit(-1);
  }

  // Have implemented some good old hard-coding for comparing many fits; argc == 2 is 1 fit, argc == 3 is 2 fits, etc
  // Could additionally do fancy stuff; do we want beam+fsi+xsec params, do we want do draw a covariance matrix, etc
  if (argc == 2) {

    std::cout << "Producing single fit output" << std::endl;
    std::string filename = argv[1];
    extractHPD(filename);
  }

  return 0;
}


// All right, let's re-write this in a faster more efficient way!
void extractHPD(std::string inputFile) {

  // drawCorr decideds if we want to draw correlations or not
  std::cout << "File for study:       " << inputFile << std::endl;
 
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
  int nddet = 0;
  bool isXsec = false;
 
  // Loop over the number of branches
  // Find the name and how many of each systematic we have
  chain->SetBranchStatus("*", false);
  chain->SetBranchStatus("step", true);

  for (int i = 0; i < nbr; i++) {
    
    // Get the TBranch and its name
    TBranch* br = (TBranch*)brlis->At(i);
    TString bname = br->GetName();

    // Get first entry of this branch (i.e. systematic); this is the first reconfigure so should be the fake-data which it was generated for (for xsec this is true...!)
    if(bname.BeginsWith("LogL")) continue;

    // If we're on beam systematics
    if(bname.BeginsWith("b_")) {
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
    } else if(bname.BeginsWith("ndd_")) {
      chain->SetBranchStatus(bname, true);
      chain->SetBranchAddress(bname, &(nom[i]));
      bnames[ndraw]=bname;
      ndraw++;
      nddet++;
    }
  }
  // Get first entry in chain
  chain->GetEntry(0);
  nom.resize(nbr);

  std::cout << "# useful entries (flux, xsec, fsi, nddet): " << ndraw << std::endl;
  std::cout << "************************************************" << std::endl;

 
  // Make sure we can read files located anywhere and strip the .root ending
  inputFile = inputFile.substr(0, inputFile.find(".root"));

  // Some TVectors with the means, errors and Gaussian parameters of the PDFs
  TVectorD* HPD_mean_vec  = new TVectorD(ndraw);
  TVectorD* HPD_err_vec   = new TVectorD(ndraw);
  TVectorD* xsec_nom  = new TVectorD(nxsec);

  // Only want to draw covariance of xsec parameters!
  TMatrixT<double>* covariance = new TMatrixT<double>(ndraw,ndraw);
  TMatrixT<double>* correlation = new TMatrixT<double>(ndraw,ndraw);
  for (int i = 0; i < ndraw-2; ++i) {
    for (int j = 0; j < ndraw-2; ++j) {
      (*covariance)(i,j) = 0.;
      (*correlation)(i,j) = 0.;
    }
  }

  TString rootfilename = inputFile;
  rootfilename += "_HPDForEd.root";

  // The output file
  TFile* file = new TFile(rootfilename, "RECREATE");
  file->cd();
 
   // Remember which parameters we're varying
  bool *vary = new bool[ndraw]();
  for (int i=0; i < ndraw; ++i) {
    vary[i] = false;
  }

  int postEbindex = -1;

  // ndraw is number of draws we want to do
  for(int i = 0; i < ndraw; ++i) {
    
    postEbindex++;
    isXsec = false;
       
    file->cd();
    int nbins = 70;

    std::string tempString = std::string(bnames[i]);
    if (bnames[i].BeginsWith("xsec_")) {
      isXsec = true;
    }    

    file->cd();
      
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
      
    double peakval, sigma_p, sigma_m, sigma_hpd;
    GetHPD(hpost, peakval, sigma_hpd, sigma_p, sigma_m);
    double mean, rms;
    GetArithmetic(hpost, mean, rms);
      
    std::cout << i << " " << peakval << "+/-" << sigma_hpd << " + " << sigma_p << " - " << sigma_m << std::endl;

    if (isXsec && ((*xsec_nom)(i-nflux)) != 0) {
      mean = mean / ((*xsec_nom)(i-nflux));
      rms = rms / ((*xsec_nom)(i-nflux));
      peakval = peakval / ((*xsec_nom)(i-nflux));
      sigma_hpd = sigma_hpd / ((*xsec_nom)(i-nflux));
    } else if (isXsec && ((*xsec_nom)(i-nflux)) == 0) {
      mean = mean + 1.0;
      peakval = peakval + 1.0;
    }

    (*HPD_mean_vec)(postEbindex) = peakval;
    (*HPD_err_vec)(postEbindex)   = sigma_hpd;
    (*covariance)(postEbindex,postEbindex)  = rms*rms;
    (*correlation)(postEbindex,postEbindex)  = 1.0;
    
    if (hpost->GetMaximum() == hpost->Integral()*1.5) {
      vary[i] = false;
      delete hpost;
      continue;
    }
    vary[i] = true;
    
    // cd into params directory in root file
    file->cd();
      
    delete hpost;   
  } // End for ndraw

  // cd into the output file
  file->cd();

  

  // Write all the nice vectors
  file->cd();
  HPD_mean_vec->Write("HPD_mean_vec");
  HPD_err_vec->Write("HPD_err_vec");
    
  file->Close();
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

int GetTH2DBinNoOverflow(TH2D* histo, int xbin, int ybin ) {

  int nbinsx = histo->GetXaxis()->GetNbins();
  return (ybin-1)*nbinsx + (xbin-1);
}
