#include <iostream>
#include <vector>

#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <THStack.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TApplication.h>

#include "samplePDF/samplePDFND2018.h"
#include "samplePDF/samplePDFND2018GPU.h"
#include "covariance/covarianceFlux.h"
#include "covariance/covarianceXsec2015.h"
#include "covariance/covarianceNDDet_2014Simple.h"

#include "mcmc/mcmc.h"
#include "manager/manager.h"

// *******************
// Class to hold the variations, throws and so on
class SampleSummary {
// *******************
  public:
    SampleSummary(int nSamples, const std::string &Outputfile, int nChainSteps);
    ~SampleSummary();

    void AddData(std::vector<TH2D*> &DataHist);
    void AddNominal(std::vector<TH2D*> &NominalHist);
    void AddThrow(std::vector<TH2D*> &MCHist, int DrawNumber = 0, double LLHPenalty = 0.0);

    void Write();

  private:
    TRandom3* rnd;
    bool first_pass;

    // Finalise the distributions from the thrown samples
    void inline MakePosteriorPredictive();

    void inline PrepareTree();

    // Helper functions to calculate likelihoods for TH1D and TH2Ds
    void inline CalcLLH(TH2D * const & Data, TH2D * const & MC);
    void inline CalcLLH(TH1D * const & Data, TH1D * const & MC);

    double inline GetLLH(TH2D * const & Data, TH2D * const & MC);
    double inline GetLLH(TH1D * const & Data, TH1D * const & MC);

    // Helper functions to change titles etc of finished plots, calculate pvalues etc
    void inline MakeCutLLH();
    void inline MakeCutLLH1D(TH1D *Histogram, TH1D *Reference = NULL);
    void inline MakeCutLLH2D(TH2D *Histogram);
    void inline MakeChi2Hists();

    // Check the length of psyche samples agrees
    bool inline CheckPsycheSamples(int Length);

    // Helper to make ratio histograms
    template<class HistType> HistType* RatioHists(HistType* NumHist, HistType* DenomHist);
    // Helper to project TH2D onto axis
    TH1D* ProjectHist(TH2D* Histogram, bool ProjectX);
    // Helper to Normalise histograms
    template<class HistType> HistType* NormaliseHist(HistType* Histogram);

    // Vector of vectors which holds the loaded MC histograms
    std::vector<std::vector<TH2D*> > MCVector;
    // Vector to hold the penalty term
    std::vector<double> LLHPenaltyVector;
    // Vectors to hold exact LLH
    std::vector<double> LLH_DrawFluc_V;
    std::vector<double> LLH_DrawData_V;
    std::vector<double> LLH_DrawFlucDraw_V;

    // A map to keep track of what psyche sample indices we want. Save some read time
    std::vector<int> PsycheSampleMap;
    std::vector<std::string> SampleNames;

    // The posterior predictive for the whole selection: this gets built after adding in the toys
    TH3D **PosteriorHist;
    // The data histogram for the selection
    TH2D **DataHist;
    // The nominal histogram for the selection
    TH2D **NominalHist;

    // The histogram containing the lnL for each throw
    TH1D *lnLHist;
    // The lnLhist for the draw vs MC
    TH1D *lnLHist_draw;
    // The lnLhist for the draw vs MC fluctuated
    TH1D *lnLHist_drawfluc;
    // The lnLhist for the draw vs MC fluctuated
    TH1D *lnLHist_drawflucdraw;
    // The lnLhist for the draw vs data
    TH1D *lnLHist_drawdata;
    // The 2D lnLhist, showing (draw vs data) and (draw vs fluct), anything above y=x axis is the p-value
    TH2D *lnLDrawHist;
    // The 2D lnLHist, showing (draw vs data) and (draw vs draw fluct), anything above y=x axis is the p-value
    TH2D *lnLFlucHist;

    // The histogram containing the lnL for each throw for each sample
    TH1D **lnLHist_Sample;

    // The LLH distribution in pmu cosmu for using the mean in each bin
    TH2D **lnLHist_Mean;
    // The LLH distribution in pmu cosmu for using the mode in each bin
    TH2D **lnLHist_Mode;
    // The posterior predictive distribution in pmu cosmu using the mean
    TH2D **MeanHist;
    // The posterior predictive distribution in pmu cosmu using the mode
    TH2D **ModeHist;

    // Holds the event rate for the distribution
    TH1D **SumHist;
    // Holds the total event rate
    TH1D *EventHist;
    // Holds the bin-by-bin LLH for the mean posterior predictive vs the data
    TH1D **lnLHist_Mean1D;
    // Holds the bin-by-bin LLH for the mode posterior predictive vs the data
    TH1D **lnLHist_Mode1D;

    // Holds the history of which entries have been drawn in the MCMC file
    TH1D *RandomHist;

    // Number of samples
    int nSamples;

    // Number of throws by user
    int nChainSteps;

    // Number of throws
    int nThrows;

    // Total LLH for the posterior predictive distribution
    double llh_total;

    // Output filename
    std::string OutputName;
    TFile *Outputfile;

    // TTree which we save useful information to
    TTree *OutputTree;
    // Data vs Draw
    double *llh_data_draw;
    // Data vs Fluctuated Draw
    double *llh_data_drawfluc;
    // Data vs Fluctuated Predictive
    double *llh_data_predfluc;

    // Draw vs Predictive
    double *llh_draw_pred;

    // Fluctuated Draw vs Predictive
    double *llh_drawfluc_pred;
    // Fluctuated Draw vs Fluctuated Predictive
    double *llh_drawfluc_predfluc;
    // Fluctuated Draw vs Draw
    double *llh_drawfluc_draw;

    // Fluctuated Predictive vs Predictive
    double *llh_predfluc_pred;
    // Fluctuated Predicitve vs Draw
    double *llh_predfluc_draw;

    // LLH penalty for each throw
    double llh_penalty;

    // Event rate for each throw
    double event_rate;

    // Data vs Draw
    double total_llh_data_draw;
    // Data vs Fluctuated Draw
    double total_llh_data_drawfluc;
    // Data vs Fluctuated Predictive
    double total_llh_data_predfluc;

    // Draw vs Predictive
    double total_llh_draw_pred;

    // Fluctuated Draw vs Predictive
    double total_llh_drawfluc_pred;
    // Fluctuated Draw vs Fluctuated Predictive
    double total_llh_drawfluc_predfluc;
    // Fluctuated Draw vs Draw
    double total_llh_drawfluc_draw;

    // Fluctuated Predictive vs Predictive
    double total_llh_predfluc_pred;
    // Fluctuated Predicitve vs Draw
    double total_llh_predfluc_draw;
};

// *******************
// The constructor
SampleSummary::SampleSummary(const int PsycheSamples, const std::string &Filename, int nSteps) {
// *******************

  std::cout << "Making sample summary class..." << std::endl;

  nChainSteps = nSteps;
  OutputName = Filename;
  nSamples = PsycheSamples;
  nThrows = 0;
  first_pass = true;
  Outputfile = NULL;
  OutputTree = NULL;

  rnd = new TRandom3();

  DataHist = new TH2D*[nSamples];
  NominalHist = new TH2D*[nSamples];
  PosteriorHist = new TH3D*[nSamples];

  lnLHist_Mean = new TH2D*[nSamples];
  lnLHist_Mode = new TH2D*[nSamples];
  MeanHist = new TH2D*[nSamples];
  ModeHist = new TH2D*[nSamples];
  SumHist = new TH1D*[nSamples];
  lnLHist_Mean1D = new TH1D*[nSamples];
  lnLHist_Mode1D = new TH1D*[nSamples];
  lnLHist_Sample = new TH1D*[nSamples];

  lnLHist = new TH1D("lnLHist_predpredfluc", "lnLHist_predpredfluc", 100, 0, 100);
  lnLHist->GetXaxis()->SetTitle("-2LLH (Pred Fluc, Pred)");
  lnLHist->GetYaxis()->SetTitle("Counts");

  lnLHist_draw = new TH1D("lnLHist_preddraw", "lnLHist_preddraw", 100, 0, 100);
  lnLHist_draw->GetXaxis()->SetTitle("-2LLH (Draw, Pred)");
  lnLHist_draw->GetYaxis()->SetTitle("Counts");

  lnLHist_drawdata = new TH1D("lnLHist_drawdata", "lnLHist_drawdata", 100, 0, 100);
  lnLHist_drawdata->GetXaxis()->SetTitle("-2LLH (Data, Draw)");
  lnLHist_drawdata->GetYaxis()->SetTitle("Counts");

  lnLHist_drawfluc = new TH1D("lnLHist_drawpredfluc", "lnLHist_drawpredfluc", 100, 0, 100);
  lnLHist_drawfluc->GetXaxis()->SetTitle("-2LLH (Pred Fluc, Draw)");
  lnLHist_drawfluc->GetYaxis()->SetTitle("Counts");

  lnLHist_drawflucdraw = new TH1D("lnLHist_drawflucdraw", "lnLHist_drawflucdraw", 100, 0, 100);
  lnLHist_drawflucdraw->GetXaxis()->SetTitle("-2LLH (Draw Fluc, Draw)");
  lnLHist_drawflucdraw->GetYaxis()->SetTitle("Counts");

  lnLDrawHist = new TH2D("lnLDrawHist", "lnLDrawHist", 100, 0, 100, 100, 0, 100);
  lnLDrawHist->GetXaxis()->SetTitle("-2LLH_{Draw, Pred Fluc}");
  lnLDrawHist->GetYaxis()->SetTitle("-2LLH_{Draw, Data}");

  lnLFlucHist = new TH2D("lnLFlucHist", "lnLFlucHist", 100, 0, 100, 100, 0, 100);
  lnLFlucHist->GetXaxis()->SetTitle("-2LLH_{Draw, Draw Fluc}");
  lnLFlucHist->GetYaxis()->SetTitle("-2LLH_{Draw, Data}");

  EventHist = new TH1D("EventHist", "EventHist", 100, 0, 100);
  EventHist->GetXaxis()->SetTitle("Total event rate");
  EventHist->GetYaxis()->SetTitle("Counts");
  EventHist->SetLineWidth(2);

  // Holds the hist of random number draws
  RandomHist = new TH1D("RandomHist", "RandomHist", 100, 0, nChainSteps);
  RandomHist->GetXaxis()->SetTitle("Step");
  double binwidth = nChainSteps/RandomHist->GetNbinsX();
  std::stringstream ss;
  ss << "Draws/" << binwidth;
  RandomHist->GetYaxis()->SetTitle(ss.str().c_str());
  RandomHist->SetLineWidth(2);

  for (int i = 0; i < nSamples; ++i) {
    DataHist[i] = NULL;
    NominalHist[i] = NULL;
    PosteriorHist[i] = NULL;
    SumHist[i] = NULL;
    MeanHist[i] = NULL;
    lnLHist_Mean[i] = NULL;
    lnLHist_Mode[i] = NULL;
    lnLHist_Mean1D[i] = NULL;
    lnLHist_Mode1D[i] = NULL;
    lnLHist_Sample[i] = NULL;
  }
}

// *******************
// Empty destructor
// Maybe should delete all the new?
SampleSummary::~SampleSummary() {
// *******************
}

// *******************
// Check size of psyche sample against size of vectors
bool SampleSummary::CheckPsycheSamples(int Length) {
// *******************
  bool ok = (SampleId::kNSamples == Length);
  if (!ok) {
    std::cerr << "Size of SampleVector input != number of psyche samples" << std::endl;
    std::cout << "Size of SampleVector: " << Length << std::endl;
    std::cout << "Size of Psyche samples: " << SampleId::kNSamples << std::endl;
    std::cerr << "Something has gone wrong with making the Samples" << std::endl;
    std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
    throw;
  }
  return ok;
}

// *******************
// Add a data histogram to the list (will have N_psyche_samples of these)
// Since the data doesn't change with varying the MC
void SampleSummary::AddData(std::vector<TH2D*> &Data) {
// *******************
  int Length = Data.size();
  // Check length of psyche samples are OK
  if (!CheckPsycheSamples(Length)) throw;
  for (int i = 0; i < Length; ++i) {
    if (Data[i] == NULL) {
      DataHist[i] = NULL;
    } else {
      PsycheSampleMap.push_back(i);
      DataHist[i] = (TH2D*)(Data[i]->Clone());
    }
  }
}

// *******************
// Add the nominal histograms to the list (will have N_psyche_samples of these)
void SampleSummary::AddNominal(std::vector<TH2D*> &Nominal) {
// *******************

  // Initialise the z-axis bins for the TH3D* PosteriorPredictive
  const int nZbins = 5000;
  double zBins[nZbins+1];
  // Go from 0 to 2000 in 5000 bins in z
  for (int i = 0; i < nZbins+1; ++i) {
    zBins[i] = 0.0;
  }

  int Length = Nominal.size();
  if (!CheckPsycheSamples(Length)) throw;
  // Loop over the lenght of nominal and set the initial distributions up
  for (int i = 0; i < Length; ++i) {

  // If NULL it indicates the selection was turned off, so initialise all the hists to NULL
    if (Nominal[i] == NULL) {
      PosteriorHist[i] = NULL;
      NominalHist[i] = NULL;
      lnLHist_Mean[i] = NULL;
      lnLHist_Mode[i] = NULL;
      MeanHist[i] = NULL;
      ModeHist[i] = NULL;
      SumHist[i] = NULL;
      lnLHist_Sample[i] = NULL;
      // If not NULL it indicates the selection was turned on, so initialise the privates
    } else {
      NominalHist[i] = (TH2D*)(Nominal[i]->Clone());
      std::string name = std::string(NominalHist[i]->GetName());
      name = name.substr(0, name.find("_nom"));

      double Maximum = NominalHist[i]->GetMaximum()*1.3;
      double Minimum = NominalHist[i]->GetMinimum()*0.7;
      for (int z = 0; z < nZbins+1; ++z) {
        zBins[z] = double(Minimum+(Maximum-Minimum))*(double(z)/double(nZbins));
      }

      PosteriorHist[i] = new TH3D((name+"_PP").c_str(), (name+"_PP").c_str(), 
          NominalHist[i]->GetXaxis()->GetNbins(), NominalHist[i]->GetXaxis()->GetXbins()->GetArray(),
          NominalHist[i]->GetYaxis()->GetNbins(), NominalHist[i]->GetYaxis()->GetXbins()->GetArray(),
          nZbins, zBins);
      PosteriorHist[i]->GetXaxis()->SetTitle(NominalHist[i]->GetXaxis()->GetTitle());
      PosteriorHist[i]->GetYaxis()->SetTitle(NominalHist[i]->GetYaxis()->GetTitle());
      PosteriorHist[i]->GetZaxis()->SetTitle("Events in bin");

      lnLHist_Mean[i] = (TH2D*)(NominalHist[i]->Clone());
      lnLHist_Mean[i]->SetNameTitle((name+"_MeanlnL").c_str(), (name+"_MeanlnL").c_str());
      lnLHist_Mean[i]->Reset();
      lnLHist_Mean[i]->GetZaxis()->SetTitle("-2lnL_{sample} #times sign(MC-data)");

      lnLHist_Mode[i] = (TH2D*)(NominalHist[i]->Clone());
      lnLHist_Mode[i]->SetNameTitle((name+"_ModelnL").c_str(), (name+"_ModelnL").c_str());
      lnLHist_Mode[i]->Reset();
      lnLHist_Mode[i]->GetZaxis()->SetTitle("-2lnL_{sample} #times sign(MC-data)");

      MeanHist[i] = (TH2D*)(NominalHist[i]->Clone());
      MeanHist[i]->SetNameTitle((name+"_mean").c_str(), (name+"_mean").c_str());
      MeanHist[i]->Reset();
      MeanHist[i]->GetZaxis()->SetTitle("Mean");

      ModeHist[i] = (TH2D*)(NominalHist[i]->Clone());
      ModeHist[i]->SetNameTitle((name+"_mode").c_str(), (name+"_mode").c_str());
      ModeHist[i]->Reset();
      ModeHist[i]->GetZaxis()->SetTitle("Mode");

      SumHist[i] = new TH1D((name+"_sum").c_str(),(name+"_sum").c_str(), 100, NominalHist[i]->Integral()*0.95, NominalHist[i]->Integral()*1.05);
      SumHist[i]->GetXaxis()->SetTitle("N_{events}");
      SumHist[i]->GetYaxis()->SetTitle("Counts");
      double Integral = DataHist[i]->Integral();
      std::stringstream ss;
      ss << Integral;
      SumHist[i]->SetTitle((std::string(SumHist[i]->GetTitle())+"_"+ss.str()).c_str());

      // Declare the lnL histograms
      // Get the likelihood for data and nominal hist
      double llh = GetLLH(DataHist[i], NominalHist[i]);

      lnLHist_Mean1D[i] = new TH1D((name+"_MeanlnL1D").c_str(),(name+"_MeanlnL1D").c_str(), 50, 0, 20);
      lnLHist_Mean1D[i]->GetXaxis()->SetTitle("-2LLH");
      lnLHist_Mean1D[i]->GetYaxis()->SetTitle("Counts");

      lnLHist_Mode1D[i] = new TH1D((name+"_ModelnL1D").c_str(),(name+"_ModelnL1D").c_str(), 50, 0, 20);
      lnLHist_Mode1D[i]->GetXaxis()->SetTitle("-2LLH");
      lnLHist_Mode1D[i]->GetYaxis()->SetTitle("Counts");

      lnLHist_Sample[i] = new TH1D((name+"_lnL").c_str(),(name+"_lnL").c_str(), 100, llh*0.9, llh*1.1);
      lnLHist_Sample[i]->GetXaxis()->SetTitle("-2LLH");
      lnLHist_Sample[i]->GetYaxis()->SetTitle("Counts");
    }
  }

}

// *******************
// Add a throw from the MCMC to the posterior predictive
// The input here is nSamples long
void SampleSummary::AddThrow(std::vector<TH2D*> &SampleVector, int DrawNumber, double LLHPenalty) {
  // *******************

  nThrows++;
  RandomHist->Fill(DrawNumber);

  int size = SampleVector.size();
  if (!CheckPsycheSamples(size)) throw;

  // Push back the throw
  MCVector.push_back(SampleVector);
  LLHPenaltyVector.push_back(LLHPenalty);

  // Loop over the sameples
  for (int z = 0; z < SampleId::kNSamples; ++z) {

    if (SampleVector[z] == NULL) continue;

    // Rebin the posterior hist
    if (first_pass) {
      // Initialise the z-axis bins for the TH3D* PosteriorPredictive
      const int nZbins = 5000;
      double zBins[nZbins+1];
      double Maximum = SampleVector[z]->GetMaximum()*1.2;
      double Minimum = SampleVector[z]->GetMinimum()*0.8;
      for (int k = 0; k < nZbins+1; ++k) {
        zBins[k] = double(Minimum+(Maximum-Minimum))*(double(k)/double(nZbins));
      }

      PosteriorHist[z]->SetBins(NominalHist[z]->GetXaxis()->GetNbins(), NominalHist[z]->GetXaxis()->GetXbins()->GetArray(), NominalHist[z]->GetYaxis()->GetNbins(), NominalHist[z]->GetYaxis()->GetXbins()->GetArray(), nZbins, zBins);
      double Minimum2 = SampleVector[z]->Integral()*0.95;
      double Maximum2 = SampleVector[z]->Integral()*1.05;
      if (Minimum2 > DataHist[z]->Integral()*0.95) {
        Minimum2 = DataHist[z]->Integral()*0.95;
      }
      if (Maximum2 < DataHist[z]->Integral()*1.05) {
        Maximum2 = DataHist[z]->Integral()*1.05;
      }
      SumHist[z]->SetBins(SumHist[z]->GetNbinsX(), Minimum2, Maximum2);
    }

    // Fill the sum histogram with the integral of the sampled distribution
    SumHist[z]->Fill(SampleVector[z]->Integral());
    // Loop over the distribution and fill the posterior predictive
    for (int i = 0; i < SampleVector[z]->GetXaxis()->GetNbins(); ++i) {
      double xCenter = SampleVector[z]->GetXaxis()->GetBinCenter(i+1);
      for (int j = 0; j < SampleVector[z]->GetYaxis()->GetNbins(); ++j) {
        double yCenter = SampleVector[z]->GetYaxis()->GetBinCenter(j+1);
        double Content = SampleVector[z]->GetBinContent(i+1, j+1);
        PosteriorHist[z]->Fill(xCenter, yCenter, Content);
      } // end y-axis loop
    } // end x-axis loop
  } // end nd280 samples loop

  first_pass = false;
} // end AddThrow

// **********************
void SampleSummary::PrepareTree() {
// **********************

  // Number of good samples
  int nPsycheSamples = PsycheSampleMap.size();

  // The array of doubles we write to the TTree
  // Data vs Draw
  llh_data_draw = new double[nPsycheSamples];
  // Data vs Fluctuated Draw
  llh_data_drawfluc = new double[nPsycheSamples];
  // Data vs Fluctuated Predictive
  llh_data_predfluc = new double[nPsycheSamples];

  // Draw vs Predictive
  llh_draw_pred = new double[nPsycheSamples];

  // Fluctuated Draw vs Predictive
  llh_drawfluc_pred = new double[nPsycheSamples];
  // Fluctuated Draw vs Fluctuated Predictive
  llh_drawfluc_predfluc = new double[nPsycheSamples];
  // Fluctuated Draw vs Draw
  llh_drawfluc_draw = new double[nPsycheSamples];
  // Fluctuated Predictive vs Predictive
  llh_predfluc_pred = new double[nPsycheSamples];
  // Fluctuated Predicitve vs Draw
  llh_predfluc_draw = new double[nPsycheSamples];

  // The output tree we're going to write to
  OutputTree = new TTree("LLH_draws", "LLH_draws");

  // Loop over the samples and set the addresses of the variables to write to file
  for (int i = 0; i < nPsycheSamples; ++i) {
    // Get the psyche sample number
    int PsycheSample = PsycheSampleMap[i];
    // Get the name
    std::string SampleName = SampleId::ConvertSample(SampleId::SampleEnum(PsycheSample));
    // Strip out spaces
    while (SampleName.find(" ") != std::string::npos) {
      SampleName.replace(SampleName.find(" "), 1, std::string("_"));
    }
    // Also strip out - signs because it messes up TBranches
    while (SampleName.find("-") != std::string::npos) {
      SampleName.replace(SampleName.find("-"), 1, std::string("_"));
    }

    SampleNames.push_back(SampleName);

    OutputTree->Branch((SampleName+"_data_draw").c_str(),     &llh_data_draw[i]);
    OutputTree->Branch((SampleName+"_data_drawfluc").c_str(), &llh_data_drawfluc[i]);
    OutputTree->Branch((SampleName+"_data_predfluc").c_str(), &llh_data_predfluc[i]);
    OutputTree->Branch((SampleName+"_draw_pred").c_str(),     &llh_draw_pred[i]);
    OutputTree->Branch((SampleName+"_drawfluc_pred").c_str(), &llh_drawfluc_pred[i]);
    OutputTree->Branch((SampleName+"_drawfluc_predfluc").c_str(), &llh_drawfluc_predfluc[i]);
    OutputTree->Branch((SampleName+"_drawfluc_draw").c_str(), &llh_drawfluc_draw[i]);
    OutputTree->Branch((SampleName+"_predfluc_pred").c_str(), &llh_predfluc_pred[i]);
    OutputTree->Branch((SampleName+"_predfluc_draw").c_str(), &llh_predfluc_draw[i]);
  }

  OutputTree->Branch("LLH_Penalty",         &llh_penalty);
  OutputTree->Branch("Total_LLH_Data_Draw", &total_llh_data_draw);
  OutputTree->Branch("Total_LLH_Data_DrawFluc", &total_llh_data_drawfluc);
  OutputTree->Branch("Total_LLH_Data_PredFluc", &total_llh_data_predfluc);
  OutputTree->Branch("Total_LLH_Draw_Pred",     &total_llh_draw_pred);
  OutputTree->Branch("Total_LLH_DrawFluc_Pred", &total_llh_drawfluc_pred);
  OutputTree->Branch("Total_LLH_DrawFluc_PredFluc", &total_llh_drawfluc_predfluc);
  OutputTree->Branch("Total_LLH_DrawFluc_Draw", &total_llh_drawfluc_draw);
  OutputTree->Branch("Total_LLH_PredFluc_Pred", &total_llh_predfluc_pred);
  OutputTree->Branch("Total_LLH_PredFluc_Draw", total_llh_predfluc_draw);
  OutputTree->Branch("Event_Rate", &event_rate);
}

// *******************
// Write the contents to the file
void SampleSummary::Write() {
  // *******************

  // Make the output file (MakePosterioPredictive call writes to this)
  std::string TempString = OutputName;
  TempString.replace(TempString.find(".root"), 5, std::string("_procs.root"));
  Outputfile = new TFile(TempString.c_str(), "RECREATE");
  // Prepare the output tree
  PrepareTree();

  std::cout << "Summarising " << nThrows << " throws..." << std::endl;
  // After all the throws are added finalise the sample
  MakePosteriorPredictive();
  std::cout << "Made posterior predictive, now writing..." << std::endl;

  Outputfile->cd();
  TDirectory *Dir[SampleId::kNSamples];

  OutputTree->Write();

  // Make the various distributions
  lnLHist->Write();
  lnLHist_draw->Write();
  lnLHist_drawfluc->Write();
  lnLHist_drawflucdraw->Write();
  lnLHist_drawdata->Write();
  lnLDrawHist->Write();
  lnLFlucHist->Write();
  RandomHist->Write();
  EventHist->Write();

  // Loop over each sample and write to file
  int nPsycheSamples = PsycheSampleMap.size();
  for (int j = 0; j < nPsycheSamples; ++j) {
    int i = PsycheSampleMap[j];
    // Skip the null histograms
    if (DataHist[i] == NULL || DataHist[i]->Integral() == 0) continue;

    // Make a new direcotry
    Dir[i] = Outputfile->mkdir(ConvertSample(SampleId::SampleEnum(i)).c_str());
    Dir[i]->cd();

    // Make the data/MC ratio histogram
    TH2D *RatioHistMean = RatioHists(DataHist[i], MeanHist[i]);
    RatioHistMean->GetZaxis()->SetTitle("Data/Mean");
    TH2D *RatioHistMode = RatioHists(DataHist[i], ModeHist[i]);
    RatioHistMode->GetZaxis()->SetTitle("Data/Mode");
    TH2D *RatioHistNom = RatioHists(DataHist[i], NominalHist[i]);
    RatioHistNom->GetZaxis()->SetTitle("Data/Nom");

    // And the normalised data histogram
    TH2D *DataNormHist = NormaliseHist(DataHist[i]);
    // Last true refers to if project along x or y
    TH1D *DataProjectX = ProjectHist(DataHist[i], true);
    TH1D *DataProjectY = ProjectHist(DataHist[i], false);

    TH2D *MeanNormHist = NormaliseHist(MeanHist[i]);
    TH2D *ModeNormHist = NormaliseHist(ModeHist[i]);
    TH1D *MeanProjectX = ProjectHist(MeanHist[i], true);
    TH1D *MeanProjectY = ProjectHist(MeanHist[i], false);
    TH1D *ModeProjectX = ProjectHist(ModeHist[i], true);
    TH1D *ModeProjectY = ProjectHist(ModeHist[i], false);

    TH2D *NomNormHist = NormaliseHist(NominalHist[i]);
    TH1D *NomProjectX = ProjectHist(NominalHist[i], true);
    TH1D *NomProjectY = ProjectHist(NominalHist[i], false);

    // Same for the TH2Ds
    CalcLLH(DataHist[i], NominalHist[i]);
    CalcLLH(DataHist[i], MeanHist[i]);
    CalcLLH(DataHist[i], ModeHist[i]);

    // Calculate the log likelihood for the 1D dists
    // Sets the title of the second TH1D to the -2LLH
    CalcLLH(DataProjectX, NomProjectX);
    CalcLLH(DataProjectX, MeanProjectX);
    CalcLLH(DataProjectX, ModeProjectX);
    CalcLLH(DataProjectY, NomProjectY);
    CalcLLH(DataProjectY, MeanProjectY);
    CalcLLH(DataProjectY, ModeProjectY);

    // For the event rate histogram add a TLine to the data rate
    TLine *TempLine = new TLine(DataHist[i]->Integral(), SumHist[i]->GetMinimum(), DataHist[i]->Integral(), SumHist[i]->GetMaximum());
    TempLine->SetLineColor(kRed);
    TempLine->SetLineWidth(2);
    // Also fit a Gaussian because why not?
    TF1 *Fitter = new TF1("Fit", "gaus", SumHist[i]->GetBinLowEdge(1), SumHist[i]->GetBinLowEdge(SumHist[i]->GetNbinsX()+1));
    SumHist[i]->Fit(Fitter, "RQ");
    Fitter->SetLineColor(kRed-5);
    // Calculate a p-value
    double Above = 0.0;
    for (int z = 0; z < SumHist[i]->GetNbinsX(); ++z) {
      double xvalue = SumHist[i]->GetBinCenter(z+1);
      if (xvalue >= DataHist[i]->Integral()) {
        Above += SumHist[i]->GetBinContent(z+1);
      }
    }
    double pvalue = Above/SumHist[i]->Integral();
    TLegend *Legend = new TLegend(0.4, 0.67, 0.98, 0.93);
    Legend->SetFillColor(0);
    Legend->SetFillStyle(0);
    Legend->SetLineWidth(0);
    Legend->SetLineColor(0);
    Legend->AddEntry(TempLine, Form("Data, %.0f, p-value=%.2f", DataHist[i]->Integral(), pvalue), "l");
    Legend->AddEntry(SumHist[i], Form("MC, #mu=%.1f#pm%.1f", SumHist[i]->GetMean(), SumHist[i]->GetRMS()), "l");
    Legend->AddEntry(Fitter, Form("Gauss, #mu=%.1f#pm%.1f", Fitter->GetParameter(1), Fitter->GetParameter(2)), "l");
    std::string TempTitle = std::string(SumHist[i]->GetName());
    TempTitle += "_canv";
    TCanvas *TempCanvas = new TCanvas(TempTitle.c_str(), TempTitle.c_str(), 1024, 1024);
    TempCanvas->SetGridx();
    TempCanvas->SetGridy();
    TempCanvas->SetRightMargin(0.03);
    TempCanvas->SetBottomMargin(0.08);
    TempCanvas->SetLeftMargin(0.10);
    TempCanvas->SetTopMargin(0.06);
    TempCanvas->cd();
    SumHist[i]->Draw();
    TempLine->Draw("same");
    Fitter->Draw("same");
    Legend->Draw("same");
    TempCanvas->Write();
    delete TempLine;
    delete TempCanvas;
    delete Fitter;
    delete Legend;
    SumHist[i]->Write();

    // Also write the 2D histograms for the p-value
    OutputTree->Draw((SampleNames[j]+"_data_draw:"+SampleNames[j]+"_predfluc_pred>>htemp").c_str());
    TH2D *TempHistogram = (TH2D*)((gDirectory->Get("htemp"))->Clone());
    TempHistogram->GetXaxis()->SetTitle("-2LLH(Pred Fluc, Pred)");
    TempHistogram->GetYaxis()->SetTitle("-2LLH(Data, Draw)");
    TempHistogram->SetNameTitle((SampleNames[j]+"_predfluc_pred").c_str(), (SampleNames[j]+"_predfluc_pred").c_str());
    MakeCutLLH2D(TempHistogram);
    TempHistogram->Write();
    delete TempHistogram;

    OutputTree->Draw((SampleNames[j]+"_data_draw:"+SampleNames[j]+"_drawfluc_draw>>htemp").c_str());
    TempHistogram = (TH2D*)((gDirectory->Get("htemp"))->Clone());
    TempHistogram->GetXaxis()->SetTitle("-2LLH(Draw Fluc, Draw)");
    TempHistogram->GetYaxis()->SetTitle("-2LLH(Data, Draw)");
    TempHistogram->SetNameTitle((SampleNames[j]+"_drawfluc_draw").c_str(), (SampleNames[j]+"_drawfluc_draw").c_str());
    MakeCutLLH2D(TempHistogram);
    TempHistogram->Write();
    delete TempHistogram;

    // Write the Histgorams to each folder
    DataHist[i]->Write();
    NominalHist[i]->Write();
    MeanHist[i]->Write();
    ModeHist[i]->Write();
    RatioHistMean->Write();
    RatioHistMode->Write();
    RatioHistNom->Write();

    DataNormHist->Write();
    NomNormHist->Write();
    MeanNormHist->Write();
    ModeNormHist->Write();

    DataProjectX->Write();
    NomProjectX->Write();
    MeanProjectX->Write();
    ModeProjectX->Write();

    DataProjectY->Write();
    NomProjectY->Write();
    MeanProjectY->Write();
    ModeProjectY->Write();

    PosteriorHist[i]->Write();
    lnLHist_Mean[i]->Write();
    lnLHist_Mode[i]->Write();
    lnLHist_Mean1D[i]->Write();
    lnLHist_Mode1D[i]->Write();
    lnLHist_Sample[i]->Write();

    // Delete temporary objects
    delete RatioHistMean;
    delete RatioHistMode;
    delete RatioHistNom;

    delete DataNormHist;
    delete MeanNormHist;
    delete ModeNormHist;
    delete NomNormHist;

    delete DataProjectX;
    delete MeanProjectX;
    delete ModeProjectX;
    delete NomProjectX;

    delete DataProjectY;
    delete MeanProjectY;
    delete ModeProjectY;
    delete NomProjectY;

    std::cout << std::endl;
  }

  std::cout << "Wrote to " << Outputfile->GetName() << std::endl;
  Outputfile->Close();
}

// *******************
// Make the posterior predictive distributions: fit Poissons etc
void SampleSummary::MakePosteriorPredictive() {
  // *******************

  // First make the projection on the z axis of the TH3D* for every pmu cosmu bin
  llh_total = 0.0;
  // Loop over the psyche samples
  for (int i = 0; i < SampleId::kNSamples; ++i) {

    // Skip disabled samples
    if (DataHist[i] == NULL || DataHist[i]->Integral() == 0) continue;

    // Count the -2LLH for each histogram
    double negLogL_Mean = 0.0;
    double negLogL_Mode = 0.0;

    // Loop over each pmu cosmu bin
    for (int j = 0; j < PosteriorHist[i]->GetXaxis()->GetNbins(); ++j) {
      for (int k = 0; k < PosteriorHist[i]->GetYaxis()->GetNbins(); ++k) {

        // Just make a little fancy name
        std::string name = PosteriorHist[i]->GetName();
        std::stringstream ss;
        ss << "_px_" << j << "_" << k;

        // Make the posterior predictive projection on z
        // The z axis of PosteriorPredictive is the bin content
        // Essentailly zooming in on one bin and looking at the mean and mode of that bin
        TH1D *Projection = PosteriorHist[i]->ProjectionZ((name+ss.str()).c_str(), j+1, j+1, k+1, k+1);
        Projection->SetAxisRange(Projection->GetMean()-16*Projection->GetRMS(), Projection->GetMean()+16*Projection->GetRMS());

        // Make the title
        std::stringstream ss2;
        ss2 << "p_{#mu} (" << PosteriorHist[i]->GetXaxis()->GetBinLowEdge(j+1) << "-" << PosteriorHist[i]->GetXaxis()->GetBinLowEdge(j+2) << ")";
        ss2 << " cos#theta_{#mu} (" << PosteriorHist[i]->GetYaxis()->GetBinLowEdge(k+1) << "-" << PosteriorHist[i]->GetYaxis()->GetBinLowEdge(k+2) << ")";
        Projection->SetNameTitle(ss2.str().c_str(), ss2.str().c_str());

        // Data content for the j,kth bin
        double nData = DataHist[i]->GetBinContent(j+1, k+1);

        // Get the mean for this projection for all the samples
        // This is the mean prediction for this given j,k bin
        double nMean = Projection->GetMean();
        double nMode = Projection->GetBinCenter(Projection->GetMaximumBin());

        double TempLLH_Mean = 0.0;
        double TempLLH_Mode = 0.0;
        if (nData > 0 && nMean > 0) {
          TempLLH_Mean = nMean - nData + nData * log(nData/nMean);
        } else if (nData == 0 && nMean > 0) {
          TempLLH_Mean = nMean;
        }

        if (nData > 0 && nMode > 0) {
          TempLLH_Mode = nMode - nData + nData * log(nData/nMode);
        } else if (nData == 0 && nMode > 0) {
          TempLLH_Mode = nMode;
        }

        // Increment -2LLH
        negLogL_Mean += 2*TempLLH_Mean;
        negLogL_Mode += 2*TempLLH_Mode;
        // Set the content to the mean in the bin
        MeanHist[i]->SetBinContent(j+1, k+1, nMean);
        // Set the content to the mode in the bin
        ModeHist[i]->SetBinContent(j+1, k+1, nMode);

        // Set the mean and average LLH for this given bin
        // Can use these hists to see where the largest -2LLH hists come from
        lnLHist_Mean[i]->SetBinContent(j+1, k+1, fabs(2.0*TempLLH_Mean));
        lnLHist_Mode[i]->SetBinContent(j+1, k+1, fabs(2.0*TempLLH_Mode));

        lnLHist_Mean1D[i]->Fill(fabs(2.0*TempLLH_Mean));
        lnLHist_Mode1D[i]->Fill(fabs(2.0*TempLLH_Mode));

        delete Projection;
      } // End loop over y bins
    } // End loop over x bins
    llh_total += negLogL_Mean;
  } // End loop over samples

  // Now we have our posterior predictive histogram and it's LLH
  std::cout << "Posterior predictive LLH mean (sample only) = " << llh_total << std::endl;
  std::stringstream ss;
  ss << llh_total;
  lnLHist->SetTitle((std::string(lnLHist->GetTitle())+"_"+ss.str()).c_str());
  lnLHist_draw->SetTitle((std::string(lnLHist_draw->GetTitle())+"_"+ss.str()).c_str());

  // Now make the fluctuated hists of the MeanHist and ModeHist
  MakeChi2Hists();

  // Get the 1D LLH dists
  MakeCutLLH();

} // End MakePosteriorPredictive() function

// *******************
// Make the fluctuated histograms (2D and 1D) for the chi2s
// Essentially taking the MCMC draws and calculating their LLH to the Posterior predicitve distribution
// And additionally taking the data histogram and calculating the LLH to the predictive distribution
// Additionally we calculate the chi2 of the draws (fluctuated) of  the MC with the posterior predictive and plot it vs the chi2 from the draws of MCMC and the data
void SampleSummary::MakeChi2Hists() {
  // *******************

  std::cout << "Making the chi2 histograms..." << std::endl;
  // Have this to signify if we're doing the first pass
  first_pass = true;

  double AveragePenalty = 0;

  // Loop over the draws
  // Should look into multi-threading this. Would require temporary THxx structures from arrays
  for (int i = 0; i < nThrows; ++i) {

    if (i % (nThrows/10) == 0) {
      std::cout << "   On throw " << i << "/" << nThrows << " (" << int(double(i)*100.0/double(nThrows)) << "%)"<< std::endl;
    }

    // Set the total LLH to zero to initialise
    total_llh_data_draw = 0.0;
    total_llh_data_drawfluc = 0.0;
    total_llh_data_predfluc = 0.0;
    total_llh_draw_pred = 0.0;
    total_llh_drawfluc_pred = 0.0;
    total_llh_drawfluc_predfluc = 0.0;
    total_llh_drawfluc_draw = 0.0;
    total_llh_predfluc_pred = 0.0;
    total_llh_predfluc_draw = 0.0;
    event_rate = 0.0;

    // Save the double that gets written to file
    llh_penalty = LLHPenaltyVector[i];
    AveragePenalty += llh_penalty;

    // Loop over the samples
    int StdIt = 0;
    for (std::vector<int>::iterator it = PsycheSampleMap.begin(); it != PsycheSampleMap.end(); ++it, StdIt++) {

      int j = *it;
      // Get the ith draw for the jth psyche sample
      TH2D *DrawHist = MCVector[i][j];
      // Skip empty samples
      if (DrawHist == NULL) continue;

      // Data vs Draw
      llh_data_draw[StdIt] = llh_penalty;
      // Data vs Fluctuated Draw
      llh_data_drawfluc[StdIt] = llh_penalty;
      // Data vs Fluctuated Predictive
      llh_data_predfluc[StdIt] = 0.0;

      // Draw vs Predictive
      llh_draw_pred[StdIt] = llh_penalty;

      // Fluctuated Draw vs Predictive
      llh_drawfluc_pred[StdIt] = llh_penalty;
      // Fluctuated Draw vs Fluctuated Predictive
      llh_drawfluc_predfluc[StdIt] = llh_penalty;
      // Fluctuated Draw vs Draw
      llh_drawfluc_draw[StdIt] = llh_penalty;
      // Fluctuated Predictive vs Predictive
      llh_predfluc_pred[StdIt] = 0.0;
      // Fluctuated Predicitve vs Draw
      llh_predfluc_draw[StdIt] = llh_penalty;

      TH2D *FluctHist = (TH2D*)(MeanHist[j]->Clone());
      TH2D *FluctDrawHist = (TH2D*)(MeanHist[j]->Clone());

      // Make the Poisson fluctuated hist
      FluctHist->Reset();
      // Also Poisson fluctuate the drawn MCMC hist
      FluctDrawHist->Reset();

      for (int x = 0; x < FluctHist->GetXaxis()->GetNbins(); ++x) {
        for (int y = 0; y < FluctHist->GetYaxis()->GetNbins(); ++y) {
          // Get the posterior predictive bin content
          double MeanContent = MeanHist[j]->GetBinContent(x+1, y+1);
          // Get a Poisson fluctuation of the content
          double Random = rnd->PoissonD(MeanContent);
          // Set the fluctuated histogram content to the Poisson variation of the posterior predictive histogram
          FluctHist->SetBinContent(x+1, y+1, Random);

          // Do the same but fluctuate the drawn MCMC step
          double DrawContent = DrawHist->GetBinContent(x+1, y+1);
          double RandomDraw = rnd->PoissonD(DrawContent);
          FluctDrawHist->SetBinContent(x+1, y+1, RandomDraw);
        }
      }

      // Likelihood between the drawn histogram and the data
      double DataDrawLLH = GetLLH(DataHist[j], DrawHist);
      llh_data_draw[StdIt] = DataDrawLLH;
      total_llh_data_draw += DataDrawLLH;

      // Likelihood between the fluctuated drawn histogram and the data
      double DataDrawFlucLLH = GetLLH(DataHist[j], FluctDrawHist);
      llh_data_drawfluc[StdIt] += DataDrawFlucLLH;
      total_llh_data_drawfluc += DataDrawFlucLLH;

      // Likelihood between the drawn histogram and the data
      double DataPredFlucLLH = GetLLH(DataHist[j], FluctHist);
      llh_data_predfluc[StdIt] += DataPredFlucLLH;
      total_llh_data_predfluc += DataPredFlucLLH;

      // Likelihood between the drawn hist and the Posterior Predictive
      double DrawPredLLH = GetLLH(DrawHist, MeanHist[j]);
      llh_draw_pred[StdIt] += DrawPredLLH;
      total_llh_draw_pred += DrawPredLLH;

      // Likelihood between fluctuated drawn and predictive
      double DrawFlucPredLLH = GetLLH(FluctDrawHist, MeanHist[j]);
      llh_drawfluc_pred[StdIt]  += DrawFlucPredLLH;
      total_llh_drawfluc_pred  += DrawFlucPredLLH;

      // Likelihood between drawn histogram and fluctuated drawn histogram
      double DrawFlucPredFlucLLH = GetLLH(FluctDrawHist, FluctHist);
      llh_drawfluc_predfluc[StdIt] += DrawFlucPredFlucLLH;
      total_llh_drawfluc_predfluc += DrawFlucPredFlucLLH;

      // Likelihood between drawn histogram and fluctuated drawn histogram
      double DrawFlucDrawLLH = GetLLH(FluctDrawHist, DrawHist);
      llh_drawfluc_draw[StdIt] += DrawFlucDrawLLH;
      total_llh_drawfluc_draw += DrawFlucDrawLLH;

      // Likelihood between the fluctuated drawn histogram and the posterior predictive
      double PredFlucPredLLH = GetLLH(FluctHist, MeanHist[j]);
      llh_predfluc_pred[StdIt] += PredFlucPredLLH;
      total_llh_predfluc_pred += PredFlucPredLLH;

      // Likelihood between drawn histogram and fluctuated posterior predictive distribution
      double PredFlucDrawLLH = GetLLH(FluctHist, DrawHist);
      llh_predfluc_draw[StdIt] += PredFlucDrawLLH;
      total_llh_predfluc_draw += PredFlucDrawLLH;

      event_rate += DrawHist->Integral();

      // Fill the random draw and posterior predictive likelihood
      if (first_pass) {
        lnLHist_Sample[j]->SetBins(lnLHist_Sample[j]->GetXaxis()->GetNbins(), DataDrawLLH*0.2, DataDrawLLH*2.5);
      }
      lnLHist_Sample[j]->Fill(DataDrawLLH);

      // Delete the fluctuated histograms
      delete FluctHist;
      delete FluctDrawHist;
    } // End loop over samples (still looping throws)


    // Add LLH penalties from the systematics to the LLH that use the drawn histogram
    total_llh_data_draw     += llh_penalty;
    total_llh_data_drawfluc += llh_penalty;
    total_llh_draw_pred     += llh_penalty;
    total_llh_drawfluc_pred += llh_penalty;
    total_llh_drawfluc_predfluc += llh_penalty;
    total_llh_drawfluc_draw += llh_penalty;
    total_llh_predfluc_draw += llh_penalty;


    if (first_pass) {
      first_pass = false;
      // Set the 1D hists
      lnLHist->SetBins(lnLHist->GetXaxis()->GetNbins(), 0.8*total_llh_predfluc_pred, 1.2*total_llh_predfluc_pred);
      lnLHist_draw->SetBins(lnLHist_draw->GetXaxis()->GetNbins(), 0.8*total_llh_draw_pred, 1.4*total_llh_draw_pred);
      lnLHist_drawdata->SetBins(lnLHist_drawdata->GetXaxis()->GetNbins(), 0.9*total_llh_data_draw, 1.1*total_llh_data_draw);
      lnLHist_drawfluc->SetBins(lnLHist_drawfluc->GetXaxis()->GetNbins(), 0.8*total_llh_predfluc_draw, 1.2*total_llh_predfluc_draw);
      lnLHist_drawflucdraw->SetBins(lnLHist_drawflucdraw->GetXaxis()->GetNbins(), 0.8*total_llh_drawfluc_draw, 1.2*total_llh_drawfluc_draw);

      // Set the 1D hists
      lnLDrawHist->SetBins(lnLDrawHist->GetXaxis()->GetNbins(), total_llh_predfluc_draw*0.8, total_llh_predfluc_draw*1.2, lnLDrawHist->GetYaxis()->GetNbins(), total_llh_data_draw*0.9, total_llh_data_draw*1.1);
      lnLFlucHist->SetBins(lnLFlucHist->GetXaxis()->GetNbins(), total_llh_drawfluc_draw*0.8, total_llh_drawfluc_draw*1.2, lnLFlucHist->GetYaxis()->GetNbins(), total_llh_data_draw*0.9, total_llh_data_draw*1.1);
      EventHist->SetBins(EventHist->GetXaxis()->GetNbins(), event_rate*0.98, event_rate*1.02);
    }

    lnLHist->Fill(total_llh_predfluc_pred);
    lnLHist_draw->Fill(total_llh_draw_pred);
    lnLHist_drawdata->Fill(total_llh_data_draw);
    lnLHist_drawfluc->Fill(total_llh_predfluc_draw);
    lnLHist_drawflucdraw->Fill(total_llh_drawfluc_draw);

    lnLDrawHist->Fill(total_llh_predfluc_draw, total_llh_data_draw);
    lnLFlucHist->Fill(total_llh_drawfluc_draw, total_llh_data_draw);

    EventHist->Fill(event_rate);

    // Also save to arrays to make sure we have the utmost super accuracy
    LLH_DrawFluc_V.push_back(total_llh_predfluc_draw);
    LLH_DrawData_V.push_back(total_llh_data_draw);
    LLH_DrawFlucDraw_V.push_back(total_llh_drawfluc_draw);

    // Write to the output tree
    OutputTree->Fill();
  } // End loop over throws

  AveragePenalty = AveragePenalty/double(nThrows);
  std::cout << "Average LLH penalty over toys is " << AveragePenalty << std::endl;

  // Calculate exact p-value instead of binned
  unsigned int Steps = LLH_DrawFluc_V.size();
  unsigned int Accept = 0;
  for (size_t i = 0; i < Steps; ++i) {
    if (LLH_DrawData_V[i] > LLH_DrawFlucDraw_V[i]) {
      Accept++;
    }
  }
  double pvalue = double(Accept)/double(Steps);
  std::stringstream ss;
  ss << pvalue << std::endl;
  std::cout << "Calculated exact p-value is " << pvalue << std::endl;
  lnLFlucHist->SetTitle((std::string(lnLFlucHist->GetTitle())+"_"+ss.str()).c_str());

  // Make the TCanvas for the total event rate
  TCanvas *TempCanv = new TCanvas((std::string(EventHist->GetName())+"_canv").c_str(), (std::string(EventHist->GetName())+"_canv").c_str(), 1024, 1024);
  double DataRate = 0.0;
  for (int i = 0; i < SampleId::kNSamples; ++i) {
    if (DataHist[i] == NULL) continue;
    DataRate += DataHist[i]->Integral();
  }
  TLine *Line = new TLine(DataRate, EventHist->GetMinimum(), DataRate, EventHist->GetMaximum());
  Line->SetLineColor(kRed);
  Line->SetLineWidth(2);
  TF1 *Fitter = new TF1("Fit", "gaus", EventHist->GetBinLowEdge(1), EventHist->GetBinLowEdge(EventHist->GetNbinsX()+1));
  EventHist->Fit(Fitter, "RQ");
  Fitter->SetLineColor(kRed-5);
  // Calculate a p-value
  double Above = 0.0;
  for (int z = 0; z < EventHist->GetNbinsX(); ++z) {
    double xvalue = EventHist->GetBinCenter(z+1);
    if (xvalue >= EventHist->Integral()) {
      Above += EventHist->GetBinContent(z+1);
    }
  }
  pvalue = Above/EventHist->Integral();
  TLegend *Legend = new TLegend(0.6, 0.6, 0.9, 0.9);
  Legend->SetFillColor(0);
  Legend->SetFillStyle(0);
  Legend->SetLineWidth(0);
  Legend->SetLineColor(0);
  Legend->AddEntry(Line, Form("Data rate, %.0f, p-value=%.2f", DataRate, pvalue), "l");
  Legend->AddEntry(EventHist, Form("Event rate, #mu=%.1f#pm%.1f", EventHist->GetMean(), EventHist->GetRMS()), "l");
  Legend->AddEntry(Fitter, Form("Gauss, #mu=%.1f#pm%.1f", Fitter->GetParameter(1), Fitter->GetParameter(2)), "l");
  TempCanv->SetGridx();
  TempCanv->SetGridy();
  TempCanv->SetRightMargin(0.03);
  TempCanv->SetBottomMargin(0.08);
  TempCanv->SetLeftMargin(0.10);
  TempCanv->SetTopMargin(0.06);
  TempCanv->cd();
  EventHist->Draw();
  Line->Draw("same");
  Fitter->Draw("same");
  Legend->Draw("same");
  TempCanv->Write();
  delete TempCanv;
  delete Line;
  delete Legend;
  delete Fitter;
}

// *******************
// Make the cut LLH histogram
void SampleSummary::MakeCutLLH() {
  // *******************
  MakeCutLLH1D(lnLHist);
  MakeCutLLH1D(lnLHist_draw);
  MakeCutLLH1D(lnLHist_drawfluc);
  MakeCutLLH1D(lnLHist_drawdata);
  MakeCutLLH1D(lnLHist_drawflucdraw);
  MakeCutLLH2D(lnLDrawHist);
  MakeCutLLH2D(lnLFlucHist);
}

// ****************
// Make the 1D cut distribution and give the 1D p-value
void SampleSummary::MakeCutLLH1D(TH1D *Histogram, TH1D* Reference) {
  // ****************
  double TotalIntegral = Histogram->Integral();
  double Above = 0.0;
  // Get the LLH reference from total llh or some reference histogram
  double llh_reference = 0.0;
  if (Reference != NULL) {
    llh_reference = GetLLH(Reference, Histogram);
  } else {
    llh_reference = llh_total;
  }
  for (int i = 0; i < Histogram->GetXaxis()->GetNbins(); ++i) {
    double xvalue = Histogram->GetBinCenter(i+1);
    if (xvalue >= llh_reference) {
      Above += Histogram->GetBinContent(i+1);
    }
  }
  double pvalue = Above/TotalIntegral;
  std::stringstream ss;
  ss << int(Above) << "/" << int(TotalIntegral) << "=" << pvalue;
  Histogram->SetTitle((std::string(Histogram->GetTitle())+"_"+ss.str()).c_str());

  // Write a TCanvas and make a line and a filled histogram
  TLine *TempLine = new TLine(llh_reference , Histogram->GetMinimum(), llh_reference, Histogram->GetMaximum());
  TempLine->SetLineColor(kBlack);
  TempLine->SetLineWidth(2);

  // Make the fill histogram
  TH1D *TempHistogram = (TH1D*)(Histogram->Clone());
  TempHistogram->SetFillStyle(1001);
  TempHistogram->SetFillColor(kRed);
  for (int i = 0; i < TempHistogram->GetNbinsX(); ++i) {
    if (TempHistogram->GetBinCenter(i+1) < llh_reference) {
      TempHistogram->SetBinContent(i+1, 0.0);
    }
  }

  std::string Title = Histogram->GetName();
  Title += "_canv";
  TCanvas *TempCanvas = new TCanvas(Title.c_str(), Title.c_str(), 1024, 1024);
  TempCanvas->SetGridx();
  TempCanvas->SetGridy();
  TempCanvas->cd();
  Histogram->Draw();
  TempHistogram->Draw("same");
  TempLine->Draw("same");

  Outputfile->cd();
  TempCanvas->Write();

  delete TempLine;
  delete TempHistogram;
  delete TempCanvas;
}


// ****************
// Make the 2D cut distribution and give the 2D p-value
void SampleSummary::MakeCutLLH2D(TH2D *Histogram) {
  // ****************

  double TotalIntegral = Histogram->Integral();
  // Count how many fills are above y=x axis
  // This is the 2D p-value
  double Above = 0.0;
  for (int i = 0; i < Histogram->GetXaxis()->GetNbins(); ++i) {
    double xvalue = Histogram->GetXaxis()->GetBinCenter(i+1);
    for (int j = 0; j < Histogram->GetYaxis()->GetNbins(); ++j) {
      double yvalue = Histogram->GetYaxis()->GetBinCenter(j+1);
      // We're only interested in being _ABOVE_ the y=x axis
      if (xvalue >= yvalue) {
        Above += Histogram->GetBinContent(i+1, j+1);
      }
    }
  }
  double pvalue = Above/TotalIntegral;
  std::stringstream ss;
  ss << int(Above) << "/" << int(TotalIntegral) << "=" << pvalue;
  Histogram->SetTitle((std::string(Histogram->GetTitle())+"_"+ss.str()).c_str());

  // Now add the TLine going diagonally
  double minimum = Histogram->GetXaxis()->GetBinLowEdge(1);
  if (Histogram->GetYaxis()->GetBinLowEdge(1) > minimum) {
    minimum = Histogram->GetYaxis()->GetBinLowEdge(1);
  }
  double maximum = Histogram->GetXaxis()->GetBinLowEdge(Histogram->GetXaxis()->GetNbins());
  if (Histogram->GetYaxis()->GetBinLowEdge(Histogram->GetYaxis()->GetNbins()) < maximum) {
    maximum = Histogram->GetYaxis()->GetBinLowEdge(Histogram->GetYaxis()->GetNbins());
  }
  TLine *TempLine = new TLine(minimum, minimum, maximum, maximum);
  TempLine->SetLineColor(kRed);
  TempLine->SetLineWidth(2);

  std::string Title = Histogram->GetName();
  Title += "_canv";
  TCanvas *TempCanvas = new TCanvas(Title.c_str(), Title.c_str(), 1024, 1024);
  TempCanvas->SetGridx();
  TempCanvas->SetGridy();
  TempCanvas->cd();
  gStyle->SetPalette(81);
  Histogram->SetMinimum(-0.01);
  Histogram->Draw("colz");
  TempLine->Draw("same");

  TempCanvas->Write();
  delete TempLine;
  delete TempCanvas;
}

// ****************
// Make a ratio histogram
template<class HistType>
HistType* SampleSummary::RatioHists(HistType *NumHist, HistType *DenomHist) {
  // ****************

  HistType *NumCopy = (HistType*)(NumHist->Clone());
  std::string title = std::string(DenomHist->GetName()) + "_ratio";
  NumCopy->SetNameTitle(title.c_str(), title.c_str());

  NumCopy->Divide(DenomHist);
  NumCopy->SetMinimum(0.4);
  NumCopy->SetMaximum(1.6);

  return NumCopy;
}

// ****************
// Normalise a histogram
template<class HistType>
HistType* SampleSummary::NormaliseHist(HistType *Histogram) {
  // ****************

  HistType* HistCopy = (HistType*)(Histogram->Clone());
  HistCopy->Scale(1./HistCopy->Integral("width"), "width");
  std::string title = std::string(HistCopy->GetName())+"_norm";
  HistCopy->SetNameTitle(title.c_str(), title.c_str());

  return HistCopy;
}


// ****************
// Calculate the LLH for TH1D, set the LLH to title of MCHist
void SampleSummary::CalcLLH(TH1D * const & DataHist, TH1D * const & MCHist) {
  // ****************
  double llh = GetLLH(DataHist, MCHist);
  std::stringstream ss;
  ss << "__n2LLH=" << llh;
  MCHist->SetTitle((std::string(MCHist->GetTitle())+ss.str()).c_str());
  std::cout << std::setw(45) << std::left << MCHist->GetName() << std::setw(10) << DataHist->Integral() << std::setw(10) << MCHist->Integral() << std::setw(10) << llh << std::endl;
}

// ****************
// Calculate the LLH for TH1D, set the LLH to title of MCHist
void SampleSummary::CalcLLH(TH2D * const & DataHist, TH2D * const & MCHist) {
  // ****************
  double llh = GetLLH(DataHist, MCHist);
  std::stringstream ss;
  ss << "__n2LLH_" << llh;
  MCHist->SetTitle((std::string(MCHist->GetTitle())+ss.str()).c_str());
  std::cout << std::setw(45) << std::left << MCHist->GetName() << std::setw(10) << DataHist->Integral() << std::setw(10) << MCHist->Integral() << std::setw(10) << llh << std::endl;
}

// ****************
double SampleSummary::GetLLH(TH2D * const & DataHist, TH2D * const & MCHist) {
  // ****************
  double llh = 0.0;

  for (int i = 0; i < DataHist->GetNbinsX(); ++i) {
    for (int j = 0; j < DataHist->GetNbinsY(); ++j) {
      double nData = DataHist->GetBinContent(i+1, j+1);
      double nMC = MCHist->GetBinContent(i+1, j+1);
      if (nMC == 0) nMC = 1E-8;
      if (nData > 0 && nMC > 0) {
        llh += 2*(nMC - nData + nData * log(nData/nMC));
      } else if (nData == 0 && nMC > 0) {
        llh += 2*nMC;
      }
    }
  }
  return llh;
}

// ****************
double SampleSummary::GetLLH(TH1D * const & DataHist, TH1D * const & MCHist) {
  // ****************
  double llh = 0.0;
  for (int i = 0; i < DataHist->GetXaxis()->GetNbins(); ++i) {
    double nData = DataHist->GetBinContent(i+1);
    double nMC = MCHist->GetBinContent(i+1);
    if (nMC == 0) nMC = 1E-8;
    if (nData > 0 && nMC > 0) {
      llh += 2*(nMC - nData + nData * log(nData/nMC));
    } else if (nData == 0 && nMC > 0) {
      llh += 2*nMC;
    }
  }
  return llh;
}

// ****************
// Make a projection
TH1D* SampleSummary::ProjectHist(TH2D* Histogram, bool ProjectX) {
  // ****************

  TH1D* Projection = NULL;
  std::string name;
  if (ProjectX) {
    name = std::string(Histogram->GetName()) + "_x";
    Projection = Histogram->ProjectionX(name.c_str(), 1, Histogram->GetYaxis()->GetNbins(), "e");
  } else {
    name = std::string(Histogram->GetName()) + "_y";
    Projection = Histogram->ProjectionY(name.c_str(), 1, Histogram->GetXaxis()->GetNbins(), "e");
  }

  return Projection;
}

// *******************
// Finally, the main
int main (int argc, char **argv) {
  // *******************

  if (argc != 2) {
    std::cout << "ND280_PosteriorPredictive_2016 runs MCMC at ND280" << std::endl;
    std::cout << "Syntax is $: ND280_PosteriorPredictive_2016 config.cfg" << std::endl;
    std::cout << "Where config.cfg is a valid config file, compatible with the manager class (manager/manager.cpp/h); an example is provided in example_config.cfg" << std::endl;
    exit(-1);
  }

  manager *fitMan = new manager(argv[1]);

  // there's a check inside the manager class that does this; left here for demonstrative purposes
  if (fitMan->getGoodConfig() == false) {
    std::cerr << "Didn't find a good config in input configuration" << std::endl;
    throw;
  }

  // Number of toys
  int Ntoys = fitMan->getNtoy();
  std::vector<double> Eb = fitMan->getOscParameters();//WP doing this to get range of parameter peaks to s\elect data from
 
 // Use the full likelihood for the posterior predictive pvalue
  bool FullLLH = fitMan->getFullLLH();

  // Check if we've already made all the data and MC
  std::string OutputName = "/vols/t2k/users/wparker2/OA2019/Eb/3OctEbData100Cut/3OctEbData100Cut.root";//fitMan->getOutputFilename();
  std::string PosteriorName = fitMan->getPosteriorFiles();
  if (PosteriorName == "EMPTY") {
    PosteriorName = OutputName;
    std::stringstream eblow, ebhigh;
    eblow << Eb[0];
    ebhigh << Eb[1];
    std::string suff = "_PostPredStore" + eblow.str() +"to" + ebhigh.str() + ".root";
    PosteriorName.replace(PosteriorName.find(".root"), 5, suff);
  }
  TFile* Outfile = new TFile(PosteriorName.c_str(), "READ");

  std::vector<TH2D*> NominalVector;
  std::vector<TH2D*> DataVector;
  std::vector< std::vector<TH2D*> > MasterSampleVector;
  // Know that these guys are going to be exist for every psyche sample, sometimes just NULL though
  NominalVector.resize(SampleId::kNSamples);
  DataVector.resize(SampleId::kNSamples);
  int nEntries = -1;
  std::vector<double> DrawVector;
  std::vector<double> PenaltyVector;

  // If the output file doesn't exist
  if (Outfile->IsZombie()) {
    delete Outfile;
    Outfile = NULL;

    MCMCProcessor Processor(OutputName, false);

    // Get inputted systematic parameters covariance matrices
    TString xsecCovMatrixFile = Processor.GetXSecCov();
    TString fluxCovMatrixFile = Processor.GetFluxCov();
    TString ndDetCovMatrixFile = Processor.GetND280Cov();

    // Check the ND runs in the run MCMC vs the given psyche
    std::string MCMCNDruns = Processor.GetNDruns();
    std::string ManNDruns = std::string(fitMan->getNDRuns());
    if (MCMCNDruns != ManNDruns) {
      std::cout << "*** WARNING: " << std::endl;
      std::cout << "    Running Posterior predictive with different ND280 runs than the MCMC in " << OutputName << " was sampled with" << std::endl;
      std::cout << "    MCMC runs: " << MCMCNDruns << std::endl;
      std::cout << "    ManNDruns: " << ManNDruns << std::endl;
      std::cout << "    This is _ENTIRELY FINE_ if you wanted to do this!" << std::endl;
      sleep(3);
    }
    // Check the ND selections
    std::vector<std::string> MCMCSel = Processor.GetNDsel();
    std::vector<std::string> ManSel = fitMan->getNDsel();
    size_t sizet = MCMCSel.size();
    std::cout << sizet << std::endl;
    if (sizet != ManSel.size()) {
      std::cout << "*** WARNING: " << std::endl;
      std::cout << "    Running Posterior predictive with different ND280 selections in the MCMC in " << OutputName << " was sampled with" << std::endl;
      std::cout << "    MCMC sels size: " << MCMCSel.size() << std::endl;
      std::cout << "    Man sels size:  " << ManSel.size() << std::endl;
      std::cout << "    This is _ENTIRELY FINE_ if you wanted to do this!" << std::endl;
      sleep(3);
    }
    for (int i = 0; i < sizet; ++i) {
      std::cout << MCMCSel[i] << " " << ManSel[i] << std::endl;
      if (MCMCSel[i] != ManSel[i]) {
        std::cout << "*** WARNING: " << std::endl;
        std::cout << "    Running Posterior predictive with different ND280 selections in the MCMC in " << OutputName << " was sampled with" << std::endl;
        std::cout << "    MCMC sels: " << MCMCSel[i] << std::endl;
        std::cout << "    Man sels:  " << ManSel[i] << std::endl;
        std::cout << "    This is _ENTIRELY FINE_ if you wanted to do this!" << std::endl;
        sleep(3);
      }
    }


#ifdef CUDA
    samplePDFND2018GPU* sample = new samplePDFND2018GPU(fitMan);
#else
    samplePDFND2018* sample = new samplePDFND2018(fitMan);
#endif

    // Setup the covariance matrices
    covarianceXsec2015 *xsec = new covarianceXsec2015("xsec_cov", xsecCovMatrixFile);
    covarianceFlux* flux = new covarianceFlux("total_flux_cov", fluxCovMatrixFile);
    covarianceNDDet_2014Simple* det = new covarianceNDDet_2014Simple("nddet_cov", ndDetCovMatrixFile);

    // Get the flat parameters
    std::vector<int> NDDetFlatParams = fitMan->getNdDetFlat();
    std::vector<int> FluxFlatParams  = fitMan->getFluxFlat();
    std::vector<int> XsecFlatParams  = fitMan->getXsecFlat();
    // Flat flux parameters loop
    if (FluxFlatParams.size() == 1 && FluxFlatParams.at(0) == -1) {
      for (int j = 0; j < flux->getSize(); j++) {
        flux->setEvalLikelihood(j, false);
      }
    } else {
      for (int j = 0; j < FluxFlatParams.size(); j++) {
        flux->setEvalLikelihood(FluxFlatParams.at(j), false);
      }
    }
    // Flat ND parameters loop
    if (NDDetFlatParams.size() == 1 && NDDetFlatParams.at(0) == -1) {
      for (int j = 0; j < det->getSize(); j++) {
        det->setEvalLikelihood(j, false);
      }
    } else {
      for (int j = 0; j < NDDetFlatParams.size(); j++) {
        det->setEvalLikelihood(NDDetFlatParams.at(j), false);
      }
    }
    // Flat xsec parameters loop
    if (XsecFlatParams.size() == 1 && XsecFlatParams.at(0) == -1) {
      for (int j = 0; j < xsec->getSize(); j++) {
        xsec->setEvalLikelihood(j, false);
      }
    } else {
      for (int j = 0; j < XsecFlatParams.size(); j++) {
        xsec->setEvalLikelihood(XsecFlatParams.at(j), false);
      }
    }

    xsec->setParameters();
    flux->setParameters();
    det->setParameters();

    sample->setXsecCov(xsec);
    sample->setFluxCov(flux);
    sample->setSimpleDetCov(det);

    sample->fillDataFromSamples();
    sample->fillReweigtingBins();
#ifdef CUDA
    sample->fillGPUSplines();
#endif

    // Do a reweight and save the nominal distributions
    double *fake = NULL;
    sample->reweight(fake);

    // If we're running fake-data fits, we should load up this as data
    bool fakeDataFit = fitMan->getFakeDataFitFlag();
    if (fakeDataFit) {
      // Maybe make a method that takes some parameters?
      sample->setAsimovFakeData(fitMan->getCustomReWeight());
    }

    for (int j = 0; j < SampleId::kNSamples; ++j) {
      if (sample->getPDF(j) == NULL) {
        NominalVector[j] = NULL;
        DataVector[j] = NULL;
        continue;
      }
      std::string name = std::string(sample->getPDF(j)->GetName());
      TH2D *TempDist = (TH2D*)(sample->getPDF(j)->Clone());
      TH2D *TempDataDist = (TH2D*)(sample->getData(j)->Clone());
      TempDist->SetNameTitle((name+"_nom").c_str(), (name+"_nom").c_str());
      TempDataDist->SetNameTitle((name+"_data").c_str(), (name+"_data").c_str());
      NominalVector[j] = TempDist;
      DataVector[j] = TempDataDist;
    }
    // Now have the nominal vectors in NominalVector

    // Open up the TChain
    TChain* mcmc = new TChain("posteriors");
    mcmc->Add(OutputName.c_str());

    TObjArray* brlis = (TObjArray*)mcmc->GetListOfBranches();
    int nbr = brlis->GetEntries();

    // Initialise the variables
    std::vector<TString> bnames(nbr,"");
    int nFlux = 0;
    int nXSec = 0;
    int nND280 = 0;
    int nDraw = 0;
    std::vector<double> fluxpars = flux->getNominalArray();
    std::vector<double> xsecpars = xsec->getNominalArray();
    std::vector<double> nddpars = det->getNominalArray();
    double LLH = 0.0;
    int step = 0;

    // Loop over the branches
    mcmc->SetBranchStatus("*", false);
    for (int i = 0; i < nbr; i++) {
      TBranch* br = (TBranch*)brlis->At(i);
      TString bname = br->GetName();

      // Read in the flux variables
      if (bname.BeginsWith("b_")) {
        mcmc->SetBranchStatus(bname, true);
        mcmc->SetBranchAddress(bname, &(fluxpars[nFlux]));
        nFlux++;
      }

      // Read in the cross-section variables
      if (bname.BeginsWith("xsec_")) {
        mcmc->SetBranchStatus(bname, true);
        mcmc->SetBranchAddress(bname, &(xsecpars[nXSec]));
        nXSec++;
      }

      // Read in the ND280 variables
      if (bname.BeginsWith("ndd_")) {
        mcmc->SetBranchStatus(bname, true);
        mcmc->SetBranchAddress(bname, &(nddpars[nND280]));
        nND280++;
      }

      // This LogL is the current LogL for each accepted step
      if (bname == "LogL") {
        mcmc->SetBranchStatus(bname, true);
        mcmc->SetBranchAddress(bname, &LLH);
      }

      // Read in the step
      if (bname == "step") {
        mcmc->SetBranchStatus(bname, true);
        mcmc->SetBranchAddress(bname, &step);
      }
    }

    nEntries = mcmc->GetEntries();
    int BurnIn = nEntries/16;
    
    // Make a vector of the TH2Ds instead of writing for every iteration
    // Holds the drawn entry numbers
    DrawVector.reserve(Ntoys);
    PenaltyVector.reserve(Ntoys);
    TRandom3 *rnd = new TRandom3();

    std::cout << "Generating " << Ntoys << std::endl;
    for (int i = 0; i < Ntoys; i++) {

      // Get a random entry after burn in
      int entry = (int)(nEntries*rnd->Rndm());
      mcmc->GetEntry(entry);

      // If we have combined chains by hadd need to check the step in the chain
      // Note, entry is not necessarily same as step due to merged ROOT files, so can't choose entry in the range BurnIn - nEntries :(
      if (step < BurnIn || xsecpars[32]<Eb[0] || xsecpars[33]>Eb[1]) {
        i--;
        continue;
      }

      // Set the flux xsec and ND280 parameters to what the MCMC stepped
      flux->setParameters(fluxpars);
      xsec->setParameters(xsecpars);
      det->setParameters(nddpars);

      // Reweight the sample with the new systematic variations
      sample->reweight(fake);

      // Calculate the likelihood
      double Sample_LLH = sample->getLikelihood();
      double Flux_LLH = flux->getLikelihood();
      double Xsec_LLH = xsec->getLikelihood();
      double ND280_LLH = det->getLikelihood();

      // Output some info for the user
      if (i % (Ntoys/10) == 0) {
        std::cout << "On toy " << i << "/" << Ntoys << " (" << int(double(i)/double(Ntoys)*100.0) << "%)" << std::endl;
        std::cout << "   Getting random entry " << entry << std::endl;
	}

      // Calculate the total combined likelihood which we can compared to the saved likelihood in the ROOT file
      double PenaltyLLH = Flux_LLH + Xsec_LLH + ND280_LLH;

      DrawVector[i] = entry;
      PenaltyVector[i] = 2.0*PenaltyLLH;

      std::stringstream ss;
      ss << "_" << i;

      std::vector<TH2D*> SampleVector;
      SampleVector.resize(SampleId::kNSamples);
      // Loop over the psyche samples and write them to file
      for (int j = 0; j < SampleId::kNSamples; ++j) {
        // Skip unloaded samples
        if (sample->getPDF(j) == NULL) {
          SampleVector[j] = NULL;
          continue;
        }
        std::string name = sample->getPDF(j)->GetName() + ss.str();
        TH2D *TempDist = (TH2D*)(sample->getPDF(j)->Clone());
        TempDist->SetNameTitle(name.c_str(), name.c_str());
        SampleVector[j] = TempDist;
      }
      MasterSampleVector.push_back(SampleVector);
    }

    // The output file
    TFile* Outfile = new TFile(PosteriorName.c_str(), "RECREATE");
    // Write the penalty for the jth toy to file
    double Penalty;
    int Draw;
    TTree *tree = new TTree("ToySummary", "ToySummary");
    tree->Branch("Penalty", &Penalty, "Penalty/D");
    tree->Branch("Draw", &Draw, "Draw/I");
    for (int j = 0; j < Ntoys; ++j) {
      Penalty = PenaltyVector[j];
      Draw = DrawVector[j];
      tree->Fill();
    }

    // Write all the templates to file
    for (int i = 0; i < SampleId::kNSamples; ++i) {
      if (DataVector[i] == NULL || DataVector[i]->Integral() == 0) continue;
      DataVector[i]->Write();
      NominalVector[i]->Write();
      for (int j = 0; j < Ntoys; ++j) {
        MasterSampleVector[j][i]->Write();
      }
    }

    Outfile->cd();
    tree->Write();
    std::cout << "Wrote " << Ntoys << " output toys to " << Outfile->GetName() << std::endl;
    Outfile->Close();
    Outfile = NULL;

  } else {

    std::cout << "Found saved posterior file " << PosteriorName << ", loading toys..." << std::endl;

    // This MasterSampleVector has indices [sample][toy] whereas the one above has [toy][sample]
    // We later sort the MasterSampleVector to take this into account so we can call AddThrow properly
    MasterSampleVector.resize(SampleId::kNSamples);
    Ntoys = -1;
    for (int i = 0; i < SampleId::kNSamples; ++i) {
      std::string Title = SampleId::ConvertSample(SampleId::SampleEnum(i));
      // Replace the spaces with underscores
      std::replace(Title.begin(), Title.end(), ' ', '_');
      TH2D *Data = (TH2D*)Outfile->Get((Title+"_data").c_str());
      TH2D *Nom = (TH2D*)Outfile->Get((Title+"_nom").c_str());
      // Pushing back NULL pointers are fine here, it's what we want
      if (Data != NULL) Data = (TH2D*)Data->Clone();
      if (Nom != NULL) Nom = (TH2D*)Nom->Clone();
      DataVector[i] = Data;
      NominalVector[i] = Nom;
      std::vector<TH2D*> SampleVector;
      // Can reserve the SampleVector if we know how many toys we have from having counted a previous sample
      if (Ntoys != -1) SampleVector.reserve(Ntoys);
      TH2D *Toy = (TH2D*)Outfile->Get((Title+"_0").c_str());
      int nToys = 0;
      if (Data != NULL) {
        std::cout << "Loading data, nominal and toys for " << Title << std::endl;
        while (Toy != NULL) {
          nToys++;
          SampleVector.push_back(Toy);
          std::stringstream ss;
          ss << "_" << nToys;
          Toy = (TH2D*)Outfile->Get((Title+ss.str()).c_str());
          if (Toy != NULL) Toy = (TH2D*)Toy->Clone();
          if (nToys % 1000 == 0) {
            std::cout << "   Loaded toy " << nToys << std::endl;
          }
        }
        Ntoys = nToys;
        MasterSampleVector[i] = SampleVector;
      }
    }

    std::cout << "Found total " << Ntoys << " pregenerated toys in " << Outfile->GetName() << std::endl;
    if (Ntoys != fitMan->getNtoy()) {
      std::cout << "*** WARNING found different number of toys in saved file than asked to run!" << std::endl;
      std::cout << "    Will read _ALL_ toys in the file" << std::endl;
      std::cout << "    Ntoys in file: " << Ntoys << std::endl;
      std::cout << "    Ntoys specified: " << fitMan->getNtoy() << std::endl;
      sleep(3);
    }

    // Then loop over and fill the MasterSampleVector (holding the toys) with NULLs where DataVector has it
    for (int i = 0; i < SampleId::kNSamples; ++i) {
      if (DataVector[i] != NULL) continue;
      std::vector<TH2D*> SampleVector;
      SampleVector.resize(Ntoys);
      for (int j = 0; j < Ntoys; ++j) {
        SampleVector[j] = NULL;
      }
      MasterSampleVector[i] = SampleVector;
    }

    // Now our MasterSampleVector has indices [sample][toy] but we want [toy][sample] for AddThrow method
    std::vector<std::vector<TH2D*> > MasterSampleCopy;
    MasterSampleCopy.resize(Ntoys);
    for (int i = 0; i < Ntoys; ++i) {
      std::vector<TH2D*> SampleVector;
      SampleVector.resize(SampleId::kNSamples);
      for (int j = 0; j < SampleId::kNSamples; ++j) {
        SampleVector[j] = MasterSampleVector[j][i];
      }
      MasterSampleCopy[i] = SampleVector;
    }
    MasterSampleVector = MasterSampleCopy;

    // Finally get the TTree branch with the penalty vectors for each of the toy throws
    TTree *PenaltyTree = (TTree*)Outfile->Get("ToySummary");
    double Penalty;
    int Draw;
    nEntries = 0;
    PenaltyVector.resize(Ntoys);
    DrawVector.resize(Ntoys);
    // If the TTree exists we can read the penalty and draw number from there
    if (PenaltyTree != NULL) {
      PenaltyTree->SetBranchAddress("Penalty", &Penalty);
      PenaltyTree->SetBranchAddress("Draw", &Draw);
      const int TreeEntries = PenaltyTree->GetEntries();
      if (TreeEntries != Ntoys) {
        std::cerr << "Number of entries in penalty tree != number of generated toys found" << std::endl;
        std::cerr << "TreeEntries: " << TreeEntries << std::endl;
        std::cerr << "nToys: " << Ntoys << std::endl;
        throw;
      }
      // Loop over the file and load entries into the penalty vector
      for (int i = 0; i < Ntoys; ++i) {
        PenaltyTree->GetEntry(i);
        if (FullLLH) {
          PenaltyVector[i] = Penalty;
        } else {
          PenaltyVector[i] = 0.0;
        }
        DrawVector[i] = Draw;
        // Find the highest draw number to make the draw plots
        if (Draw > nEntries) nEntries = Draw;
      }
    } else {
      for (int i = 0; i < Ntoys; ++i) {
        PenaltyVector[i] = 0.;
        DrawVector[i] = 1;
      }
    }
  }

  std::cout << "Done making templates and reweighting... Now making the SampleSummary instance" << std::endl;

  // Set the output name depending on if we've done the full LLH or not
  if (FullLLH) {
    PosteriorName.replace(PosteriorName.find(".root"), 5, std::string("_FullLLH.root"));
  } else {
    PosteriorName.replace(PosteriorName.find(".root"), 5, std::string("_SampLLH.root"));
  }
  SampleSummary Samples(SampleId::kNSamples, PosteriorName, nEntries);

  // Add in the Data and nominal to the SampleSummary class
  Samples.AddData(DataVector);
  Samples.AddNominal(NominalVector);

  // Now we have the MasterSampleVector filled with our toys
  // Add the throws
  for (int i = 0; i < Ntoys; ++i) {
    Samples.AddThrow(MasterSampleVector[i], DrawVector[i], PenaltyVector[i]);
  }
  std::cout << "Finished adding throws to SampleSummary, now writing..." << std::endl;

  // Calculate the posterior predictive p-values and write the samples to file
  Samples.Write();

  return 0;
}

