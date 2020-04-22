// Convert ND280 binning
// Compiled because that means we can link to psyche and avoid some hard-coding!
// g++ `root-config --cflags` -g -o ConvertND280Bin ConvertND280Bin.cpp -I`root-config --incdir` -I../psycheinc `root-config --glibs --libs` -L../lib
// ********************************
// ______WARNING______
// THE ORDER OF SAMPLES IN THE FOR LOOP MATTERS AND ____IS HARD CODED____
// FOR 2017 ANALYSES THE SAMPLE ORDERING WAS FGD1 SAMPLES, FGD2 SAMPLES, GOING __AGAINST__ THE PSYCHE NUMBERING
// FOR 2018 ANALYSES THE SAMPLE ORDERING FOLLOWED PSYCHE NUMBERING
// CHECK GIT LOG TO SEE WHAT WAS DONE
// ______END WARNING______

#include <iostream>
#include <iomanip>
#include <string>

#include "TFile.h"
#include "TAxis.h"
#include "TObjArray.h"
#include "TH1I.h"
#include "TVectorD.h"
#include "TMatrixDSym.h"

#include "CoreDataClasses.hxx"
#include "CoreUtils.hxx"
#include "BaseDataClasses.hxx"
#include "SampleId.hxx"

void ConvertND280Binning(const std::string &);


int main(int argc, char **argv) {

  if (argc != 2) {
    std::cout << "./ConvertND280Bin inputFile.root" << std::endl;
    exit(-1);
  }

  std::string infile = argv[1];
  ConvertND280Binning(infile);

  return 0;
};


void ConvertND280Binning(const std::string & inputFile) {

  TFile *f = new TFile(inputFile.c_str(), "READ");
  if (f->IsZombie()) {
    std::cerr << "Can't find file " << inputFile << std::endl;
    std::cerr << "I'll just have to hang up!" << std::endl;
    exit(-1);
  }

  // Check the input file: unfortunately need to change this code for every new input covariance
  // Any of the input files would be fine to check, I've just chosen the GOF for no reason
  TH1D *temp = (TH1D*)(f->Get("Nominal"));
  if (temp->GetNbinsX() != 556) {
    std::cerr << "This executable is hard-coded for 2016/2017 analyses" << std::endl;
    std::cerr << "Am expecting 556 bins in total, I found " << temp->GetNbinsX() << std::endl;
    //   throw;
  }

  const int psycheSamples = SampleId::kNSamples;

  // Let's hard-code the whole bloody binning here; it's the only way of getting the binning correct
  // Unfortunately the input file has no information on what momentum/angle combination the bin corresponds too: there are only bin numbers
  //
  // This should probably be fixed!
  //
  // From Simon Bienstock's talk on 25 Oct 2016. Slide 4 was decided to be the best binning for the detector covariance
  /*
  ////////// Old binning
  // FHC CC0pi muon momentum
  const int n_cc0pi_pmu = 6;
  const double cc0pi_pmu[n_cc0pi_pmu+1] = {0, 1000, 1250, 2000, 3000, 5000, 30000};
  // FHC CC0pi muon cos theta
  const int n_cc0pi_thmu = 7;
  const double cc0pi_thmu[n_cc0pi_thmu+1] = {-1.0, 0.6, 0.7, 0.8, 0.85, 0.94, 0.96, 1.0};
  TAxis *axis_cc0pi_pmu = new TAxis(n_cc0pi_pmu, cc0pi_pmu);
  TAxis *axis_cc0pi_thmu = new TAxis(n_cc0pi_thmu, cc0pi_thmu);

  // FHC CC1pi muon momentum
  const int n_cc1pi_pmu = 5;
  const double cc1pi_pmu[n_cc1pi_pmu+1] = {0, 300, 1250, 1500, 5000, 30000};
  const int n_cc1pi_thmu = 8;
  double cc1pi_thmu[n_cc1pi_thmu+1] = {-1.0, 0.7, 0.85, 0.9, 0.92, 0.96, 0.98, 0.99, 1.0};
  TAxis *axis_cc1pi_pmu = new TAxis(n_cc1pi_pmu, cc1pi_pmu);
  TAxis *axis_cc1pi_thmu = new TAxis(n_cc1pi_thmu, cc1pi_thmu);

  // FHC CCOth muon momentum
  const int n_ccOth_pmu = 5;
  const double ccOth_pmu[n_ccOth_pmu+1] = {0, 1500, 2000, 3000, 5000, 30000};
  const int n_ccOth_thmu = 8;
  double ccOth_thmu[n_ccOth_thmu+1] = {-1.0, 0.8, 0.85, 0.9, 0.92, 0.96, 0.98, 0.99, 1.0};
  TAxis *axis_ccOth_pmu = new TAxis(n_ccOth_pmu, ccOth_pmu);
  TAxis *axis_ccOth_thmu = new TAxis(n_ccOth_thmu, ccOth_thmu);

  // RHC antinu CC1track
  const int n_cc1Trk_nubar_pmu = 5;
  const double cc1Trk_nubar_pmu[n_cc1Trk_nubar_pmu+1] = {0, 400, 900, 1100, 2000, 10000};
  const int n_cc1Trk_nubar_thmu = 8;
  double cc1Trk_nubar_thmu[n_cc1Trk_nubar_thmu+1] = {-1.0, 0.6, 0.7, 0.88, 0.95, 0.97, 0.98, 0.99, 1.00};
  TAxis *axis_cc1Trk_nubar_pmu = new TAxis(n_cc1Trk_nubar_pmu, cc1Trk_nubar_pmu);
  TAxis *axis_cc1Trk_nubar_thmu = new TAxis(n_cc1Trk_nubar_thmu, cc1Trk_nubar_thmu);

  // RHC antinu CCNtrack
  const int n_ccNTrk_nubar_pmu = 6;
  const double ccNTrk_nubar_pmu[n_ccNTrk_nubar_pmu+1] = {0, 700, 1200, 1500, 2000, 3000, 10000};
  const int n_ccNTrk_nubar_thmu = 6;
  double ccNTrk_nubar_thmu[n_ccNTrk_nubar_thmu+1] = {-1.0, 0.85, 0.88, 0.93, 0.98, 0.99, 1.00};
  TAxis *axis_ccNTrk_nubar_pmu = new TAxis(n_ccNTrk_nubar_pmu, ccNTrk_nubar_pmu);
  TAxis *axis_ccNTrk_nubar_thmu = new TAxis(n_ccNTrk_nubar_thmu, ccNTrk_nubar_thmu);

  // RHC nu in anu CC1track
  const int n_cc1Trk_nu_pmu = 5;
  const double cc1Trk_nu_pmu[n_cc1Trk_nu_pmu+1] = {0, 400, 800, 1100, 2000, 10000};
  const int n_cc1Trk_nu_thmu = 8;
  double cc1Trk_nu_thmu[n_cc1Trk_nu_thmu+1] = {-1.0, 0.7, 0.85, 0.90, 0.93, 0.96, 0.98, 0.99, 1.00};
  TAxis *axis_cc1Trk_nu_pmu = new TAxis(n_cc1Trk_nu_pmu, cc1Trk_nu_pmu);
  TAxis *axis_cc1Trk_nu_thmu = new TAxis(n_cc1Trk_nu_thmu, cc1Trk_nu_thmu);

  // RHC nu in an CCNtrack
  const int n_ccNTrk_nu_pmu = 5;
  const double ccNTrk_nu_pmu[n_ccNTrk_nu_pmu+1] = {0, 1000, 1500, 2000, 3000, 10000};
  const int n_ccNTrk_nu_thmu = 8;
  double ccNTrk_nu_thmu[n_ccNTrk_nu_thmu+1] = {-1.0, 0.8, 0.9, 0.93, 0.95, 0.96, 0.97, 0.99, 1.00};
  TAxis *axis_ccNTrk_nu_pmu = new TAxis(n_ccNTrk_nu_pmu, ccNTrk_nu_pmu);
  TAxis *axis_ccNTrk_nu_thmu = new TAxis(n_ccNTrk_nu_thmu, ccNTrk_nu_thmu);
*/


  //////////////////////////////////// New 574 Binning
  // FHC CC0pi muon mom
  const int n_cc0pi_pmu = 8;
  const double cc0pi_pmu[n_cc0pi_pmu+1] = {0., 300., 1000., 1250., 1500., 2000., 3000., 5000., 30000.};
  // FHC CC0pi muon cos theta
  const int n_cc0pi_thmu = 8;
  const double cc0pi_thmu[n_cc0pi_thmu+1] = {-1.0, 0.6, 0.8, 0.85, 0.9, 0.92, 0.98, 0.99, 1.0};
  TAxis *axis_cc0pi_pmu = new TAxis(n_cc0pi_pmu, cc0pi_pmu);
  TAxis *axis_cc0pi_thmu = new TAxis(n_cc0pi_thmu, cc0pi_thmu);

  // FHC CC1pi muon momentum
  const int n_cc1pi_pmu = 9;
  const double cc1pi_pmu[n_cc1pi_pmu+1] = {0., 300., 400., 700., 800., 1000., 1500., 2000., 5000., 30000.};
  const int n_cc1pi_thmu = 9;
  double cc1pi_thmu[n_cc1pi_thmu+1] = {-1.0, 0.6, 0.8, 0.9, 0.92, 0.94, 0.96, 0.98, 0.99, 1.0};
  TAxis *axis_cc1pi_pmu = new TAxis(n_cc1pi_pmu, cc1pi_pmu);
  TAxis *axis_cc1pi_thmu = new TAxis(n_cc1pi_thmu, cc1pi_thmu);

  // FHC CCOth muon momentum
  const int n_ccOth_pmu = 10;
  const double ccOth_pmu[n_ccOth_pmu+1] = {0., 300., 400., 700., 800., 900., 1250., 2000., 3000., 5000., 30000.};
  const int n_ccOth_thmu = 9;
  double ccOth_thmu[n_ccOth_thmu+1] = {-1.0, 0.6, 0.8, 0.85, 0.9, 0.92, 0.96, 0.98, 0.99, 1.0};
  TAxis *axis_ccOth_pmu = new TAxis(n_ccOth_pmu, ccOth_pmu);
  TAxis *axis_ccOth_thmu = new TAxis(n_ccOth_thmu, ccOth_thmu);

  // RHC CC0pi muon momentum
  const int n_cc0pi_nubar_pmu = 4;
  const double cc0pi_nubar_pmu[n_cc0pi_nubar_pmu+1] = {0., 300., 2000., 4000., 30000.};
  // RHC CC0pi muon cos theta
  const int n_cc0pi_nubar_thmu = 5;
  const double cc0pi_nubar_thmu[n_cc0pi_nubar_thmu+1] = {-1., 0.6, 0.8, 0.9, 0.96, 1.};
  TAxis *axis_cc0pi_nubar_pmu = new TAxis(n_cc0pi_nubar_pmu, cc0pi_nubar_pmu);
  TAxis *axis_cc0pi_nubar_thmu = new TAxis(n_cc0pi_nubar_thmu, cc0pi_nubar_thmu);

  // RHC CC1pi muon momentum
  const int n_cc1pi_nubar_pmu = 2;
  const double cc1pi_nubar_pmu[n_cc1pi_nubar_pmu+1] = {0., 500., 30000.};
  const int n_cc1pi_nubar_thmu = 2;
  double cc1pi_nubar_thmu[n_cc1pi_nubar_thmu+1] = {-1, 0.7, 1.};
  TAxis *axis_cc1pi_nubar_pmu = new TAxis(n_cc1pi_nubar_pmu, cc1pi_nubar_pmu);
  TAxis *axis_cc1pi_nubar_thmu = new TAxis(n_cc1pi_nubar_thmu, cc1pi_nubar_thmu);

  // RHC CCOth muon momentum
  const int n_ccOth_nubar_pmu = 3;
  const double ccOth_nubar_pmu[n_ccOth_nubar_pmu+1] = {0., 600., 800., 30000.};
  const int n_ccOth_nubar_thmu = 4;
  double ccOth_nubar_thmu[n_ccOth_nubar_thmu+1] = {-1., 0.7, 0.95, 0.97, 1.};
  TAxis *axis_ccOth_nubar_pmu = new TAxis(n_ccOth_nubar_pmu, ccOth_nubar_pmu);
  TAxis *axis_ccOth_nubar_thmu = new TAxis(n_ccOth_nubar_thmu, ccOth_nubar_thmu);

  // RHCnu CC0pi muon momentum
  const int n_cc0pi_nubarnu_pmu = 3;
  const double cc0pi_nubarnu_pmu[n_cc0pi_nubarnu_pmu+1] = {0., 300., 1500., 30000.};
  // RHCnu CC0pi muon cos theta
  const int n_cc0pi_nubarnu_thmu = 2;
  const double cc0pi_nubarnu_thmu[n_cc0pi_nubarnu_thmu+1] = {-1., 0.7, 1.};
  TAxis *axis_cc0pi_nubarnu_pmu = new TAxis(n_cc0pi_nubarnu_pmu, cc0pi_nubarnu_pmu);
  TAxis *axis_cc0pi_nubarnu_thmu = new TAxis(n_cc0pi_nubarnu_thmu, cc0pi_nubarnu_thmu);

  // RHCnu CC1pi muon momentum
  const int n_cc1pi_nubarnu_pmu = 3;
  const double cc1pi_nubarnu_pmu[n_cc1pi_nubarnu_pmu+1] = {0., 600., 800., 30000.};
  const int n_cc1pi_nubarnu_thmu = 2;
  double cc1pi_nubarnu_thmu[n_cc1pi_nubarnu_thmu+1] = {-1, 0.7, 1.};
  TAxis *axis_cc1pi_nubarnu_pmu = new TAxis(n_cc1pi_nubarnu_pmu, cc1pi_nubarnu_pmu);
  TAxis *axis_cc1pi_nubarnu_thmu = new TAxis(n_cc1pi_nubarnu_thmu, cc1pi_nubarnu_thmu);

  // RHCnu CCOth muon momentum
  const int n_ccOth_nubarnu_pmu = 2;
  const double ccOth_nubarnu_pmu[n_ccOth_nubarnu_pmu+1] = {0., 600., 30000.};
  const int n_ccOth_nubarnu_thmu = 2;
  double ccOth_nubarnu_thmu[n_ccOth_nubarnu_thmu+1] = {-1., 0.7, 1.};
  TAxis *axis_ccOth_nubarnu_pmu = new TAxis(n_ccOth_nubarnu_pmu, ccOth_nubarnu_pmu);
  TAxis *axis_ccOth_nubarnu_thmu = new TAxis(n_ccOth_nubarnu_thmu, ccOth_nubarnu_thmu);

  
  // RHC antinu CC1track
  const int n_cc1Trk_nubar_pmu = 5;
  const double cc1Trk_nubar_pmu[n_cc1Trk_nubar_pmu+1] = {0, 400, 900, 1100, 2000, 10000};
  const int n_cc1Trk_nubar_thmu = 8;
  double cc1Trk_nubar_thmu[n_cc1Trk_nubar_thmu+1] = {-1.0, 0.6, 0.7, 0.88, 0.95, 0.97, 0.98, 0.99, 1.00};
  TAxis *axis_cc1Trk_nubar_pmu = new TAxis(n_cc1Trk_nubar_pmu, cc1Trk_nubar_pmu);
  TAxis *axis_cc1Trk_nubar_thmu = new TAxis(n_cc1Trk_nubar_thmu, cc1Trk_nubar_thmu);

  // RHC antinu CCNtrack
  const int n_ccNTrk_nubar_pmu = 6;
  const double ccNTrk_nubar_pmu[n_ccNTrk_nubar_pmu+1] = {0, 700, 1200, 1500, 2000, 3000, 10000};
  const int n_ccNTrk_nubar_thmu = 6;
  double ccNTrk_nubar_thmu[n_ccNTrk_nubar_thmu+1] = {-1.0, 0.85, 0.88, 0.93, 0.98, 0.99, 1.00};
  TAxis *axis_ccNTrk_nubar_pmu = new TAxis(n_ccNTrk_nubar_pmu, ccNTrk_nubar_pmu);
  TAxis *axis_ccNTrk_nubar_thmu = new TAxis(n_ccNTrk_nubar_thmu, ccNTrk_nubar_thmu);

  // RHC nu in anu CC1track
  const int n_cc1Trk_nu_pmu = 5;
  const double cc1Trk_nu_pmu[n_cc1Trk_nu_pmu+1] = {0, 400, 800, 1100, 2000, 10000};
  const int n_cc1Trk_nu_thmu = 8;
  double cc1Trk_nu_thmu[n_cc1Trk_nu_thmu+1] = {-1.0, 0.7, 0.85, 0.90, 0.93, 0.96, 0.98, 0.99, 1.00};
  TAxis *axis_cc1Trk_nu_pmu = new TAxis(n_cc1Trk_nu_pmu, cc1Trk_nu_pmu);
  TAxis *axis_cc1Trk_nu_thmu = new TAxis(n_cc1Trk_nu_thmu, cc1Trk_nu_thmu);

  // RHC nu in an CCNtrack
  const int n_ccNTrk_nu_pmu = 5;
  const double ccNTrk_nu_pmu[n_ccNTrk_nu_pmu+1] = {0, 1000, 1500, 2000, 3000, 10000};
  const int n_ccNTrk_nu_thmu = 8;
  double ccNTrk_nu_thmu[n_ccNTrk_nu_thmu+1] = {-1.0, 0.8, 0.9, 0.93, 0.95, 0.96, 0.97, 0.99, 1.00};
  TAxis *axis_ccNTrk_nu_pmu = new TAxis(n_ccNTrk_nu_pmu, ccNTrk_nu_pmu);
  TAxis *axis_ccNTrk_nu_thmu = new TAxis(n_ccNTrk_nu_thmu, ccNTrk_nu_thmu);
  

  // All the types that MaCh3 need
  // This holds e.g. lepton momentum (x axis)
  TObjArray* k1_axes = new TObjArray();
  // This holds e.g. lepton cos theta (y axis)
  TObjArray* k2_axes = new TObjArray();
  // This holds the offset in going from a bin on x and y axis to a global bin
  TH1I* offset = new TH1I("offset", "offset", psycheSamples, 0, psycheSamples);

  // This holds the actual covariance, just read straight out of the file
  TMatrixDSym *nddet_cov = (TMatrixDSym*)(f->Get("Covariance_Matrix")->Clone());
  // This holds the central values, just read straight out of the file
  TVectorD* det_weight = (TVectorD*)(f->Get("Mean_Value")->Clone());

  // This is the dummy for when we have a psyche sample but it isn't included in the switch below (so isn't included as a sample in ND280 fits!)
  TAxis *dummy = NULL;

  // Have a counter that stores the offset for each sample
  int offsetCounter = 0;

  // Loop over the psyche samples
  // According to Simon the ordering of the bins is FGD1 CC0pi, CC1pi... FGD2 CC0pi, CC1pi...
  //
  // So am hacking this together somewhat: loop over the psyche sample, picking out the FGD1 selections only
  // Then do the same loop again but pick out the FGD2 selections instead
  // This calls two for loops which are identical _EXCEPT_ for the cases in the switch

  for (int i = 0; i < psycheSamples; ++i) {

    // Make the x axis from the momentum of lepton
    TAxis* xaxis;
    // Make the y axis from the cos theta of lepton
    TAxis* yaxis;
    
    switch (i) {

      // FGD1 and FGD2 selections share same binning
      // CC 0pi
      case SampleId::kFGD1NuMuCC0Pi:
      case SampleId::kFGD2NuMuCC0Pi:
        xaxis = axis_cc0pi_pmu;
        yaxis = axis_cc0pi_thmu;
        break;

      // CC 1pi
      case SampleId::kFGD1NuMuCC1Pi:
      case SampleId::kFGD2NuMuCC1Pi:
        xaxis = axis_cc1pi_pmu;
        yaxis = axis_cc1pi_thmu;
        break;

      // CC Oth
      case SampleId::kFGD1NuMuCCOther:
      case SampleId::kFGD2NuMuCCOther:
        xaxis = axis_ccOth_pmu;
        yaxis = axis_ccOth_thmu;
        break;
	
      // FGD1 and FGD2 selections share same binning
      // CC 0pi
      case SampleId::kFGD1AntiNuMuCC0Pi:
      case SampleId::kFGD2AntiNuMuCC0Pi:
        xaxis = axis_cc0pi_nubar_pmu;
        yaxis = axis_cc0pi_nubar_thmu;
        break;

      // CC 1pi
      case SampleId::kFGD1AntiNuMuCC1Pi:
      case SampleId::kFGD2AntiNuMuCC1Pi:
        xaxis = axis_cc1pi_nubar_pmu;
        yaxis = axis_cc1pi_nubar_thmu;
        break;

      // CC Oth
      case SampleId::kFGD1AntiNuMuCCOther:
      case SampleId::kFGD2AntiNuMuCCOther:
        xaxis = axis_ccOth_nubar_pmu;
        yaxis = axis_ccOth_nubar_thmu;
        break;

      // FGD1 and FGD2 selections share same binning
      // Neutrino in anti-neutrino CC 0pi
      case SampleId::kFGD1NuMuBkgInAntiNuModeCC0Pi:
      case SampleId::kFGD2NuMuBkgInAntiNuModeCC0Pi:
        xaxis = axis_cc0pi_nubarnu_pmu;
        yaxis = axis_cc0pi_nubarnu_thmu;
        break;

      // Neutrino in anti-neutrino CC 1pi
      case SampleId::kFGD1NuMuBkgInAntiNuModeCC1Pi:
      case SampleId::kFGD2NuMuBkgInAntiNuModeCC1Pi:
        xaxis = axis_cc1pi_nubarnu_pmu;
        yaxis = axis_cc1pi_nubarnu_thmu;
        break;

      // Neutrino in anti-neutrino CC Oth
      case SampleId::kFGD1NuMuBkgInAntiNuModeCCOther:
      case SampleId::kFGD2NuMuBkgInAntiNuModeCCOther:
        xaxis = axis_ccOth_nubarnu_pmu;
        yaxis = axis_ccOth_nubarnu_thmu;
        break;

      // Anti-neutrino CC 1 track
      case SampleId::kFGD1AntiNuMuCC1Track:
      case SampleId::kFGD2AntiNuMuCC1Track:
        xaxis = axis_cc1Trk_nubar_pmu;
        yaxis = axis_cc1Trk_nubar_thmu;
        break;

      // Anti-neutrino CC N track
      case SampleId::kFGD1AntiNuMuCCNTracks:
      case SampleId::kFGD2AntiNuMuCCNTracks:
        xaxis = axis_ccNTrk_nubar_pmu;
        yaxis = axis_ccNTrk_nubar_thmu;
        break;

      // Neutrino in anti-neutrino CC 1 track
      case SampleId::kFGD1NuMuBkgInAntiNuModeCC1Track:
      case SampleId::kFGD2NuMuBkgInAntiNuModeCC1Track:
        xaxis = axis_cc1Trk_nu_pmu;
        yaxis = axis_cc1Trk_nu_thmu;
        break;

      // Neutrino in anti-neutrino CC N track
      case SampleId::kFGD1NuMuBkgInAntiNuModeCCNTracks:
      case SampleId::kFGD2NuMuBkgInAntiNuModeCCNTracks:
        xaxis = axis_ccNTrk_nu_pmu;
        yaxis = axis_ccNTrk_nu_thmu;
        break;

      default:
        xaxis = dummy;
        yaxis = dummy;
        break;
    }


    // Need the number of total bins for the selection
    // This is needed because there's an "offset" histogram in MaCh3
    // This histogram stores where each selection's bins start in the global bin numbering
    // The global bin numbering runs from 0~580, where the first n bins are one selection, next bins are another selection, and so on
    std::string name = SampleId::ConvertSample(SampleId::SampleEnum(i));

    if (yaxis != NULL && xaxis != NULL) {
      int size = xaxis->GetNbins()*yaxis->GetNbins();
      std::cout << name << ":" << std::endl;
      std::cout << std::setw(10) << "offset: " << offsetCounter << std::endl;
      std::cout << std::setw(10) << "xaxis: " << xaxis->GetNbins() << std::endl;
      std::cout << std::setw(10) << "yaxis: " << yaxis->GetNbins() << std::endl;
      std::cout << std::setw(10) << "nbins: " << size << std::endl;
      std::cout << std::endl;

      offset->SetBinContent(i+1, offsetCounter);
      offsetCounter += size;

      // Chuck axes into the TObjArray
      k1_axes->AddAtAndExpand(xaxis, i);
      k2_axes->AddAtAndExpand(yaxis, i);
    }
    offset->GetXaxis()->SetBinLabel(i+1, name.c_str());

  }

  std::string outputName = inputFile;
  while (outputName.find(".") != std::string::npos) {
    outputName = outputName.substr(0, outputName.find("."));
  }
  outputName += "_convertedToMaCh3.root";

  TFile *output = new TFile(outputName.c_str(), "RECREATE");

  k1_axes->Write("k1_axes", 1);
  k2_axes->Write("k2_axes", 1);
  offset->Write("offset");
  nddet_cov->Write("nddet_cov");
  det_weight->Write("det_weights");

  output->Close();
  f->Close();

};
