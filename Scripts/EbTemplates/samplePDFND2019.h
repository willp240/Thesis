#ifndef _samplePDFND2019_h_
#define _samplePDFND2019_h_

// Do we use TF1 or TSpline3* for spline evaluations
#define USE_TSpline3 1
#define USE_TSpline3_red 2
#define USE_TF1 3
#define USE_TF1_red 4

// Can use:
//  TSpline3 (third order spline in ROOT)
//  TSpline3_red (reduced class third order spline in ROOT)
//  TF1 (fifth order poly in ROOT)
//  TF1_red (reduced class fifth order poly)

#define USE_SPLINE USE_TSpline3

// Set the __SPLINE_TYPE__ accordingly
#if USE_SPLINE == USE_TSpline3
#define __SPLINE_TYPE__ TSpline3
#elif USE_SPLINE == USE_TSpline3_red
#define __SPLINE_TYPE__ TSpline3_red
#elif USE_SPLINE == USE_TF1
#define __SPLINE_TYPE__ TF1
#elif USE_SPLINE == USE_TF1_red
#define __SPLINE_TYPE__ TF1_red
#endif

// Do we want to debug the TF1 treatment? This involves writing TSpline3* and TF1* to file and comparing how well they fit TGraph*
//#define DEBUG_TF1


#ifdef MULTITHREAD
#include "omp.h"
#endif

// C++ includes
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <sstream>

// ROOT include
#include <TH1D.h>
#include <TH2D.h>
#include <TH2Poly.h>
#include <TSpline.h>
#include <TMath.h>
#include <TFile.h>
#include <TTree.h>
#include <TRandom3.h>
#include <TStopwatch.h>
#include <TString.h>
#include <TMatrixDSym.h>
#include <TLegend.h>
#include <TLine.h>

// Psyche includes
#include "psycheinc/AnalysisManager.hxx"
#include "psycheinc/AnalysisUtils.hxx"
#include "psycheinc/ToyMakerExample.hxx"
#include "psycheinc/ND280AnalysisUtils.hxx"

// MaCh3 samplePDF includes
#include "samplePDFBase.h"
#include "Structs.h"

// MaCh3 covariance includes
#include "covariance/covarianceFlux.h"
#include "covariance/covarianceXsec.h"
#include "covariance/covarianceXsec2015.h"
#include "covariance/covarianceNDFSI.h"
#include "covariance/covariancePsyche.h"
#include "covariance/covarianceNDDet_2019Poly.h"

// MaCh3 effective rpa for 2017
#include "eRPA/eRPA.h"
// MINERvA fake-data studies for 2017
#include "eRPA/MINERvAq0q3.h"

// Include the manager for an alternate constructor
#include "manager/manager.h"

//NIWGReweight includes
#include "NIWGReWeight/NIWGFSLepEbMomentumVariation_continuous.h"
#include "NIWGReWeight/NIWGFSLepEbMomentumVariation_utils.h"
#include "NIWGReWeight/NIWGFSLepEbMomentumVariation.h"

// Make an enum of the test statistic that we're using
enum TestStatistic {
  kPoisson,
  kBarlowBeeston,
  kIceCube
};

// Do same as BANFF here: replaces the ToyMakerExample in psycheSteering
// It's identical to BANFF
class NominalToyMaker : public ToyMaker {
  public:
    NominalToyMaker(void): ToyMaker(), _sigma(0.0), _weight(1.0) {}
    // Need to override destructor because ToyMaker destructor doesn't work...
    ~NominalToyMaker() {};

    void FillToyExperiment(ToyExperiment& toy) {
      std::cout << "Number of systematics in psyche = " << _nSystematics << std::endl;
      for (UInt_t isyst = 0; isyst<NMAXSYSTEMATICS; isyst++ ) {
        SystematicBase* syst = _systematics[isyst];
        if (!syst) continue;
        std::cout << "Making toy variation for systematic " << syst->Name() << " which has " << syst->GetNParameters() << " params" << std::endl;

        // Create the proper structure for the ToyExperiment adding each of the systematics
        toy.AddToyVariation(syst->GetIndex(), syst->GetNParameters());

        // Loop over parameters for this systematic
        for (UInt_t ipar = 0; ipar < syst->GetNParameters(); ipar++ ) {
          toy.SetToyVariation(syst->GetIndex(), ipar, _sigma, _weight);
        }
      }

    }
  private:
    Float_t _sigma;
    Float_t _weight;
};


class samplePDFND2019 : public samplePDFBase {
  public:
      samplePDFND2019(manager *FitManager);
      ~samplePDFND2019();

      // The different covariance matrices to be associated with the samplePDF
      void setFluxCov(covarianceFlux * const flx) { flux = flx; }
      void setSimpleDetCov(covarianceNDDet_2019Poly * const indet) { detsimple = indet; }
      void setXsecCov(covarianceXsec2015 * const xs, bool norms = false);

      covarianceFlux * const GetFluxCov() const { return flux; }
      covarianceNDDet_2019Poly * const GetSimpleDetCov() const { return detsimple; }
      covarianceXsec2015 * const GetXsecCov() const { return xs2015; }

      // Function to set the Asimov fake-data from within samplePDF
      // Put this back in
      void setAsimovFakeData(bool CustomReWeight = false);
      void setAsimovFakeData_FromFile(std::string &FileName);
      void setAsimovFakeDataThrow();

      double getLikelihood();
      double getTestStatLLH(double data, double mc, double w2);
      double getSampleLikelihood(int isample);
      double getEventRate(int index); 

      // The reweighting functions
      void reweight(double *oscpar);
      void reweight(double *oscpar, double *oscpar2) { reweight(oscpar); };

      // DEPRECATED
      std::vector< std::vector<double> > generate2D(int index);
      std::vector<double> generate(int index);
      std::vector<double>* getDataSample() {return NULL;}
      
      void addData(std::vector<double> &dat) {return;}
      void addData(std::vector< vector <double> > &dat) {return;}
      void addData(TH1D* binneddata){return;}
      void addData(TH2Poly* binneddata){return;}
      
      void addData(std::vector<double> &dat, int index);
      void addData(std::vector< std::vector <double> > &dat, int index);
      void addData(TH1D* binneddata, int index);
      void addData(TH2Poly* binneddata, int index);

      // ADD method!
      void printPosteriors() {return;}
      void addSelection(SampleId::SampleEnum selection, ND280KinematicTypes var1=kLeptonMomentum, ND280KinematicTypes var2=kLeptonCosTheta);

      // Randomize the starting position
      void RandomStart();

      // Provide a setter for the test-statistic
      void SetTestStatistic(TestStatistic test_stat) {
        fTestStatistic = test_stat;
        if (fTestStatistic == kPoisson) std::cout << "Using Poisson likelihood in ND280" << std::endl;
        else if (fTestStatistic == kBarlowBeeston) std::cout << "Using Barlow-Beeston likelihood in ND280" << std::endl;
        else if (fTestStatistic == kIceCube) std::cout << "Using IceCube likelihood in ND280" << std::endl;
        else {
          std::cerr << "UNKNOWN LIKELHOOD SPECIFIED TO ND280!" << std::endl;
          std::cerr << "You gave test-statistic " << fTestStatistic << std::endl;
          std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
          throw;
        }
      }

      // Getters for the event histograms
      TH1* getPDF(int Selection)  { return (TH1*)samplepdfs->At(Selection); }
      TH1* getData(int Selection) { return (TH1*)datapdfs->At(Selection); }
      TH1* getPDFMode(int Selection, int Mode) { return (TH1*)((TObjArray*)samplemodepdfs->At(Selection))->At(Mode); }

      void GetKinVars(ND280KinematicTypes &TypeX, ND280KinematicTypes &TypeY);

      void fillDataFromSamples();
      void fillReweigtingBins();
      void fillFluxandXsec();

      virtual void enableModeHistograms();

      void printMCRates();
      void printDataRates();
      void printRates();

      covariancePsyche* getDetectorCovariance() { return det; }

      void savePDF(std::string saveLoc);

      // DEPRECATED
      // {
      std::vector< std::vector <double> > generate2D(TH2Poly* pdf) {std::vector< std::vector <double> > null; return null;}
      std::vector<double> generate() {std::vector<double> null; return null;}
      // }

  protected:
      // Using TF1 or TSpline3 for spline evaluation?
      bool tf1;

      // Helper function to find normalisation bins for normalisation parameters that aren't simply a mode scaling
      // e.g. different normalisation parameters for carbon and oxygen, neutrino and anti-neutrino, etc
      inline void FindNormBins();

      // Set each events normalisation bin
      void SetEventNormBin(const int EventNumber);

      // Make the TChain which connects our psyche events and cross-section model
      TChain* MakeXsecChain();

      // Setup the binning for a given selection
      void SetupBinning(SampleId::SampleEnum Selection, std::vector<double> &BinningX, std::vector<double> &BinningY);

      void SetupPoly(std::string, TH2Poly*&);

      // Reserve spline memory 
      // Virtual for GPU
      void ReserveMemory(int nEventMC);
      void ShrinkMemory(int nEventMC);

      // Quick helpers to load psyche and samples from the fit manager
      inline void LoadPsyche();
      inline void LoadSamples();
      inline void PrintPOT();

      // virtual for GPU
      void SetSplines(TGraph** &xsecgraph, const int i);
#if USE_SPLINE < USE_TF1
      void SetSplines_Reduced(TGraph** &xsecgraph, const int i);
#endif
      //void SetFunct(TGraph** &xsecgraph, const int);

      // Helper function to check if the covariances have been set
      inline void CheckCovariances();

      // Redirect std::cout to silence psyche
      void QuietPlease();
      void NowTalk();

      // DEPRECATED BUT NEEDED BECAUSE OF INHERITANCE, POOP
      // PLEASE FIX IF YOU HAVE TIME
      // {
      double getCovLikelihood() {return 0.0;}
      double getLikelihood_kernel(std::vector<double> &data) {return 0.0;}
      void fill1DHist() {return;}
      void fill2DHist() {return;}
      // }

      // Get a kinematic variable from a psyche AnaEventB object
      inline double getKinematicVariable(ND280KinematicTypes variable, AnaEventB* event);
      // Get highest momentum particle
      inline TLorentzVector GetHMParticle(AnaEventB const * event, int PDG);

      // Perform the main reweight loop
      void ReWeight_NoPsyche();
#ifdef MULTITHREAD
      void ReWeight_NoPsyche_MP();
#endif
      // Prepare weights
      // virtual for GPU
      virtual void PrepareWeights();

      // Helper function to reset histograms
      inline void ResetHistograms();

      // Calculate the cross-section weight for a given event for NIWG 2017
      double CalcXsecWeight(const int EventNumber);
      // Calculate the spline weight for a given event for NIWG 2017
      // virtual for GPU
      virtual double CalcXsecWeight_Spline(const int EventNumber);
      // Calculate the tf1 weight for a given event for NIWG 2017
      //virtual double CalcXsecWeight_Funct(const int EventNumber);
      // Calculate the spline weight for a given event for NIWG 2017
      double CalcXsecWeight_Norm(const int EventNumber);

#ifdef DEBUG
      // Dump spline information
      inline void DumpSplines(double xsecw, double beamw, double detw, int EventNumber);
      // Print information about the xsec and ndobj structs
      inline void PrintStructs(double xsecw, double beamw, double detw, int i);
#endif

      bool isRHC(int runperiod);
      bool HaveIRandomStart; // Have I random started?

      // The covariance classes
      covarianceFlux * flux;
      covarianceXsec2015 * xs2015;
      covariancePsyche * det;
      covarianceNDDet_2019Poly * detsimple;

      // This is the number of cross-section splines we're loading
      // It's set in the constructor/init, looking only at TString prod
      int nXsecSplines;
      // This one is read from setCovMatrix function which looks at the input covariance and counts the number of cross-section splines
      int nXsecSplinesFromCov;

      // Needed for openMP
      // Contains how many samples we've got
      __int__ nSamples;
      // what is the maximum number of bins we have
      __int__ maxBins;

      // Struct containing the cross-section info
      xsec2018<__SPLINE_TYPE__*>* xsecInfo;

      // Struct containing the ND280 information
      ND280MC* ndobj;

      //WP: The EB varier
      niwg::var::NIWGFSLepEbMomentumVariation_continuous* varier_eb;

      // String with what production we want
      TString *prod;

      // Psyche classes
      AnalysisManager _man;
      ToyExperiment *toy;
      Experiment *exp;
      NominalToyMaker ToyMaker;

      // Dimensions of the ith psyche selection (2D or 1D)
      int ndims[SampleId::kNSamples];
      // The kinematic type we're plotting
      ND280KinematicTypes kinvars[SampleId::kNSamples][2];

      // The pdfs we're storing
      // MC pdfs
      TObjArray* samplepdfs;
      // data pdfs
      TObjArray* datapdfs;
      // mode of MC pdfs
      TObjArray* samplemodepdfs;
      TH2PolyBin ***polybins;
      TH2Poly***** ebtemplates;

      TH1D* NumuErec;
      TH1D* NumubarErec;
      TH1D* NueErec;
      TH1D* NuebarErec;
      
      // do we want mode MC pdf to be save
      bool modepdf;
      TObjArray* modeobjarray[SampleId::kNSamples];

      // A nice random
      TRandom3* rand;

      // number of cross-section normalisation params?
      int nxsec_norm_modes;
      int *xsec_norm_modes;
      int *xsec_norm_startbin;
      std::vector<std::string> xsec_norm_names;
      // Information about the normliastion parameters
      std::vector<XsecNorms2> xsec_norms;

      // Number of normalisation parameters
      // Used in the !DoItSmart version
      // The old version hard-codes these instead
      int nNormParams;
      std::vector<int> normParsModes;
      std::vector< std::vector<int> > normParsTarget;
      std::vector<int> normParsIndex;
      std::vector< std::vector<int> > normParsPDG;

      // Number of spline parameters
      int nSplineParams;
      std::vector<int> splineParsIndex;
      std::vector<std::string> splineParsNames;

      // Number of spline parameters that aren't repeated
      int nSplineParamsUniq;
      std::vector<int> splineParsUniqIndex;
      std::vector<std::string> splineParsUniqNames;

      // Number of spline parameters that are repeated
      int nSplineParamsShare;
      std::vector<int> splineParsShareIndex;
      std::vector<int> splineParsShareToUniq;
      std::vector<std::string> splineParsShareNames;

      // Number of function parameters
      int nFuncParams;
      std::vector<int> funcParsIndex;
      std::vector<std::string> funcParsNames;

      // Bit-field comparison
      std::vector<int> xsecBitField;

      // Takes a spline number and maps it to the correct global cross-section parameter
      std::vector<int> splineToGlobal;
      // Takes a norm. number and maps it to the correct global cross-section parameter
      std::vector<int> normToGlobal;

      // Use covariance matrix for detector
      bool simple;

      // noRFG?
      bool norfg;

      // NIWG xsec model
      int xsecmodel;

      // The position of the neutrino/anti-neutrino C/O scaling parameters
      int mec_norm_nu_pos;
      int mec_norm_anu_pos;
      int mec_norm_CtoO_pos;
      int coh_c_pos;
      int coh_o_pos;
      int nue_norm_nu_pos;
      int nue_norm_anu_pos;

      //Position of the Q2 parameters
      int Q2_norm_0_pos;
      int Q2_norm_1_pos;
      int Q2_norm_2_pos;
      int Q2_norm_3_pos;
      int Q2_norm_4_pos;

      //Position of the Eb Dial parameters
      int Eb_C_nu_pos;
      int Eb_C_nubar_pos;
      int Eb_O_nu_pos;
      int Eb_O_nubar_pos;

      // Position of the 2018 CC parameters
      int CC_nu_pos;
      int CC_nubar_pos;

      //Position of CC Misc DIS and MultPi normalisations
      int CC_misc_norm_pos;
      int CC_dis_multipi_norm_nu_pos;
      int CC_dis_multipi_norm_nubar_pos;

      // Array of FastSplineInfo structs: keeps information on each xsec spline for fast evaluation
      // Method identical to TSpline3::Eval(double) but faster because less operations
#if USE_SPLINE < USE_TF1
      FastSplineInfo *SplineInfoArray;
      template <class T> 
        double FastSplineEval(T* spline, const int SplineNumber);
      void FindSplineSegment();
#endif

#ifdef DEBUG_TF1
      TFile *dumpfile;
#endif

      inline double ApplyCoulombShift(AnaEventB* Summary, int target, int reaction, ND280KinematicTypes kKintype, double InitialValue = 0.0);

      inline double CalculateQ2(double PLep, double PUpd, double EnuTrue, double InitialQ2 = 0.0);
      inline double CalculateEnu(double PLep, double cosTheta, double EB, bool neutrino);
      // Pointer to fit manager
      manager *FitManager;

  private:
      int getValue();
      int parseLine(char* line);
      std::streambuf *buf; // Keep the cout buffer
      std::streambuf *errbuf; // Keep the cerr buffer

      TestStatistic fTestStatistic;

      // A vector of TH2Ds to hold the weight^2 for e.g. Barlow Beeston
      std::vector<TH2Poly*> W2Hist;

#if DEBUG > 0
      unsigned int nEventsInTree;
      // Vector of POT weights needed for pot+xsec weights
      std::vector<double> POTWeights;

      // 2D histograms which match BANFF's numbers
      std::vector<TH2Poly*> MCOnly;
      std::vector<TH2Poly*> POTOnly;
      std::vector<TH2Poly*> FluxOnly;
      std::vector<TH2Poly*> XsecOnly;
      std::vector<TH2Poly*> DetOnly;
      std::vector<TH2Poly*> NDCovOnly;
      std::vector<TH2Poly*> BeamCovOnly;
      std::vector<TH2Poly*> AllOnly;

      // TH1D which matches NuMu group's binning
      std::vector<TH1D*> NuMuPmu;

      // Events per normalisation parameter
      std::vector<unsigned int> EventsPerNorm;
      // Events per normalisation mode
      std::vector<unsigned int> EventsPerMode;
#endif

};

#if USE_SPLINE < USE_TF1
// ***************************************************************************
// Fast spline evaluation
// Inspired by TSpline3::Eval, which we can speed up considerably
// Main reason is that we know that for one parameter (e.g. MAQE) we will have the same number of points, x min, xmax, etc for all MAQE splines, so we can signficantly reduce number of operations
// The curious can find very similar GPU code in splines/gpuSplineUtils.cu and CPU code in spline/SplineMonolith.cpp::Eval
// I've included it here for more transparency: this kind of eval should be possible for SK splines too
template <class T>
double samplePDFND2019::FastSplineEval(T* spline, const int SplineNumber) {
  // ***************************************************************************

  // Check if the spline is NULL
  if (spline == NULL) return 1.0;

  // The segment has already been found in FindSplineSegment()
  int segment = SplineInfoArray[SplineNumber].CurrSegment;

  // These are what we can extract from the TSpline3
  double x = -999.99;
  double y = -999.99;
  double b = -999.99;
  double c = -999.99;
  double d = -999.99;

  // Now write the coefficients
  spline->GetCoeff(segment, x, y, b, c, d);

  // Get the variation for this reconfigure for the ith parameter
  int GlobalIndex = splineParsIndex[SplineNumber];
  double xvar = xs2015->calcReWeightFrac(GlobalIndex);

  // The Delta(x)
  double dx = xvar - x;

  // The spline weight to return
  double weight = y+dx*(b+dx*(c+d*dx));

  // Check that eval on the TSpline3 is the same as our weight
#ifdef DEBUG
  // Difference between eval and weight
  double diff = fabs(spline->Eval(xs2015->calcReWeightFrac(GlobalIndex)) - weight);

  if (diff > 1.E-7) {

    std::cerr << "TSpline3->Eval() != custom eval: Something is wrong with FastSplineEval!" << std::endl;
    std::cerr << "Difference in spline evaluation > 1.E-5!" << std::endl;

    std::cerr << "Eval      = " << spline->Eval(xs2015->calcReWeightFrac(GlobalIndex)) << std::endl;
    std::cerr << "Cust      = " << weight << std::endl;
    std::cerr << "diff      = " << diff << std::endl;
    std::cerr << "param     = " << splineParsNames[SplineNumber] << std::endl;
    std::cerr << "variation = " << xs2015->calcReWeightFrac(GlobalIndex) << std::endl;
    std::cerr << "paramVal  = " << xs2015->getParProp(GlobalIndex) << std::endl;

    // Check we've found the right segment
    int klow = spline->FindX(xvar);
    if (klow >= spline->GetNp()-1 && spline->GetNp() > 1) klow = spline->GetNp()-2;

    std::cerr << "segment   = " << segment << std::endl;
    std::cerr << "spl segm  = " << klow << std::endl;
    std::cerr << "nPoints   = " << SplineInfoArray[SplineNumber].nPts << std::endl;
    std::cerr << "nPoints sp= " << spline->GetNp() << std::endl;


    std::cerr << "Printing x information:" << std::endl;
    for (int i = 0; i < SplineInfoArray[SplineNumber].nPts; ++i) {
      std::cerr << "   " << i << " = " << SplineInfoArray[SplineNumber].xPts[i] << std::endl;
    }
    std::cerr << "From spline: " << std::endl;
    for (int i = 0; i < spline->GetNp(); ++i) {
      double xtmp, ytmp;
      spline->GetKnot(i, xtmp, ytmp);
      std::cerr << "   " << i << " = " << xtmp << std::endl;
    }

    if (klow != segment) {
      std::cerr << "Have found wrong segment in FindSplineSegment!" << std::endl;
      std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
      throw;
    }

    std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
    throw;
  }
#endif
  return weight;
}
#endif


#endif
