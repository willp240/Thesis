// g++ -o Mergebins $(root-config --cflags) MergeBins.cpp $(root-config --libs)

#include <vector>
#include <iostream>
#include <iomanip>

#include <TAxis.h>
#include <TFile.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <map>

TFile *Infile = NULL;

//////////////////////////////////////
TAxis *FGD1_CC0pi_p = NULL;
TAxis *FGD1_CC0pi_th = NULL;
TAxis *FGD1_CC1pi_p = NULL;
TAxis *FGD1_CC1pi_th = NULL;
TAxis *FGD1_CCOth_p = NULL;
TAxis *FGD1_CCOth_th = NULL;

TAxis *FGD2_CC0pi_p = NULL;
TAxis *FGD2_CC0pi_th = NULL;
TAxis *FGD2_CC1pi_p = NULL;
TAxis *FGD2_CC1pi_th = NULL;
TAxis *FGD2_CCOth_p = NULL;
TAxis *FGD2_CCOth_th = NULL;

TAxis *FGD1_nubar_CC0pi_p = NULL;
TAxis *FGD1_nubar_CC0pi_th = NULL;
TAxis *FGD1_nubar_CC1pi_p = NULL;
TAxis *FGD1_nubar_CC1pi_th = NULL;
TAxis *FGD1_nubar_CCOth_p = NULL;
TAxis *FGD1_nubar_CCOth_th = NULL;

TAxis *FGD2_nubar_CC0pi_p = NULL;
TAxis *FGD2_nubar_CC0pi_th = NULL;
TAxis *FGD2_nubar_CC1pi_p = NULL;
TAxis *FGD2_nubar_CC1pi_th = NULL;
TAxis *FGD2_nubar_CCOth_p = NULL;
TAxis *FGD2_nubar_CCOth_th = NULL;

TAxis *FGD1_nubar_nu_CC0pi_p = NULL;
TAxis *FGD1_nubar_nu_CC0pi_th = NULL;
TAxis *FGD1_nubar_nu_CC1pi_p = NULL;
TAxis *FGD1_nubar_nu_CC1pi_th = NULL;
TAxis *FGD1_nubar_nu_CCOth_p = NULL;
TAxis *FGD1_nubar_nu_CCOth_th = NULL;

TAxis *FGD2_nubar_nu_CC0pi_p = NULL;
TAxis *FGD2_nubar_nu_CC0pi_th = NULL;
TAxis *FGD2_nubar_nu_CC1pi_p = NULL;
TAxis *FGD2_nubar_nu_CC1pi_th = NULL;
TAxis *FGD2_nubar_nu_CCOth_p = NULL;
TAxis *FGD2_nubar_nu_CCOth_th = NULL;

int nBinsFGD1CC0pi;
int nBinsFGD1CC1pi;
int nBinsFGD1CCOth;

int nBinsFGD2CC0pi;
int nBinsFGD2CC1pi;
int nBinsFGD2CCOth;

int nBinsFGD1CC0pi_nubar;
int nBinsFGD1CC1pi_nubar;
int nBinsFGD1CCOth_nubar;

int nBinsFGD2CC0pi_nubar;
int nBinsFGD2CC1pi_nubar;
int nBinsFGD2CCOth_nubar;

int nBinsFGD1CC0pi_nubarnu;
int nBinsFGD1CC1pi_nubarnu;
int nBinsFGD1CCOth_nubarnu;

int nBinsFGD2CC0pi_nubarnu;
int nBinsFGD2CC1pi_nubarnu;
int nBinsFGD2CCOth_nubarnu;

int nBins;

std::map<int, std::vector<bool> >mergemomboolmap;
std::map<int, std::vector<bool> >mergethetaboolmap;

std::vector<int> BinVector;
std::vector<TAxis*> MomVector;
std::vector<TAxis*> CosVector;

TH2D* Covariance = NULL;
TH1D *Nominal = NULL;

void FindBins();
bool PmuMerge(int i, int sample, double pmu, double cosmu);
bool CosMerge(int i, int sample, double pmu, double cosmu);
void FindSamplePmuCosmu(int bin, int &selection, double &pmu, double &cosmu);
double threshold;

//////////////////////////////////////

void MergeBins(std::string File) {

  Infile = new TFile(File.c_str(), "READ");

  Covariance = (TH2D*)(Infile->Get("Total_Covariance_Matrix")->Clone());
  Nominal = (TH1D*)(Infile->Get("Nominal")->Clone());

  FindBins();

  // Get number of bins
  nBins = Nominal->GetXaxis()->GetNbins();
  
  TCanvas *Canv = new TCanvas("Canv", "Canv", 1024, 1024);
  TH1D* TheOutputCov = new TH1D("outputcov", "outputcov", nBins, 0, nBins);

  int nMergesPmu = 0;
  int nMergesCos = 0;
  int nMerges = 0;

  for (int i = 0; i < nBins; ++i) {
    // Check the sigma of each bin
    double pmu, cosmu;
    int sample;
    // Get sample, pmu, cosmu from i
    FindSamplePmuCosmu(i, sample, pmu, cosmu);
    TheOutputCov->SetBinContent(i+1, (Covariance->GetBinContent(i, i)));
    //std::cout << i << ": " << sample << "  " << std::setw(5) << pmu << "  " << std::setw(5) << cosmu << std::endl;
    // Now check if we can merge all the bins in one sub-bin
    // Try pmu first
    bool okpmu = CosMerge(i, sample, pmu, cosmu);
    if (okpmu) {
      std::cout << "Ok to merge sample " << sample << " pmu: " << pmu << " cosmu: " << cosmu << " in cosmu UPWARDS" << std::endl;
      nMergesPmu++;
      mergemomboolmap[sample].push_back(false);
      mergethetaboolmap[sample].push_back(true);
    }

    bool okcosmu = false;
    // Try cosmu first
    if (!okpmu) {
      okcosmu = PmuMerge(i, sample, pmu, cosmu);
      if (okcosmu) {
	std::cout << "Ok to merge sample " << sample << " pmu: " << pmu << " cosmu: " << cosmu << " in mom UPWARDS" << std::endl;
        nMergesCos++;
	mergemomboolmap[sample].push_back(true);
	mergethetaboolmap[sample].push_back(false);
      }
    }
    if (!okpmu && !okcosmu) {
      std::cout << " Do not merge sample " << sample << " pmu: " << pmu << " cosmu: " << cosmu << " in momentum UPWARDS" << std::endl;
      mergemomboolmap[sample].push_back(false);
      mergethetaboolmap[sample].push_back(false);
    }
  }

  std::cout << "***************************************************************" << std::endl;


  int counter=0, prevsample = -999;

  for (int i = 0; i < nBins; ++i) {

    double pmu, cosmu;
    int sample;
    // Get sample, pmu, cosmu from i
    FindSamplePmuCosmu(i, sample, pmu, cosmu);
    if(sample!=prevsample)
      counter = 0;

    if((sample<3 || (sample<9&&sample>5) || (sample<15&&sample>11))  && mergemomboolmap[sample][counter] && mergemomboolmap[sample+3][counter]){
      std::cout << "Ok to merge sample " << sample << " pmu: " << pmu << " cosmu: " << cosmu << " in mom UPWARDS" << std::endl; 
      nMerges++;
    }
    else if((sample<3 || (sample<9&&sample>5) || (sample<15&&sample>11))  && mergethetaboolmap[sample][counter] && mergethetaboolmap[sample+3][counter]){
      std::cout << "Ok to merge sample " << sample << " pmu: " << pmu << " cosmu: " << cosmu << " in cos UPWARDS" << std::endl;
      nMerges++;
    }
    else if(sample<3 || (sample<9&&sample>5) || (sample<15&&sample>11)){
      std::cout << "Do not merge sample " << sample << " pmu: " << pmu << " cosmu: " << cosmu << " in momentum UPWARDS" << std::endl;
    }

    counter++;
    prevsample = sample;
  }


  Canv->cd();
  TheOutputCov->Draw();
  Canv->Print("outputcov.pdf");
  std::cout << "***************************************************************" << std::endl;
  std::cout << "Found " << nMergesPmu << "/" << nBins << " possible merges in pmu" << std::endl;
  std::cout << "Found " << nMergesCos << "/" << nBins << " possible merges in cos" << std::endl;
  std::cout << "Found " << nMergesPmu+nMergesCos << "/" << nBins << " possible merges in cos" << std::endl;
  std::cout << "Bringing total bin number to " << nBins-nMergesPmu-nMergesCos << std::endl;
  std::cout << "With threshold " << threshold << std::endl;
  std::cout << "***************************************************************" << std::endl;
  std::cout << "Combining FGD1 and 2: " << std::endl << "Found " << 2*nMerges << "/" << nBins << std::endl; 
  std::cout << "Bringing total bin number to " << nBins-(2*nMerges) << std::endl;
}

void FindBins() {

  //////////////////////////////////////
  FGD1_CC0pi_p = (TAxis*)(Infile->Get( "AxisForCov/FGD1_numuCC_0pi_p_axis"));
  FGD1_CC0pi_th = (TAxis*)(Infile->Get("AxisForCov/FGD1_numuCC_0pi_th_axis"));
  FGD1_CC1pi_p = (TAxis*)(Infile->Get( "AxisForCov/FGD1_numuCC_1pi_p_axis"));
  FGD1_CC1pi_th = (TAxis*)(Infile->Get("AxisForCov/FGD1_numuCC_1pi_th_axis"));
  FGD1_CCOth_p = (TAxis*)(Infile->Get( "AxisForCov/FGD1_numuCC_other_p_axis"));
  FGD1_CCOth_th = (TAxis*)(Infile->Get("AxisForCov/FGD1_numuCC_other_th_axis"));

  FGD2_CC0pi_p = (TAxis*)(Infile->Get( "AxisForCov/FGD2_numuCC_0pi_p_axis"));
  FGD2_CC0pi_th = (TAxis*)(Infile->Get("AxisForCov/FGD2_numuCC_0pi_th_axis"));
  FGD2_CC1pi_p = (TAxis*)(Infile->Get( "AxisForCov/FGD2_numuCC_1pi_p_axis"));
  FGD2_CC1pi_th = (TAxis*)(Infile->Get("AxisForCov/FGD2_numuCC_1pi_th_axis"));
  FGD2_CCOth_p = (TAxis*)(Infile->Get( "AxisForCov/FGD2_numuCC_other_p_axis"));
  FGD2_CCOth_th = (TAxis*)(Infile->Get("AxisForCov/FGD2_numuCC_other_th_axis"));
  FGD1_nubar_CC0pi_p = (TAxis*)(Infile->Get( "AxisForCov/FGD1_anti-numuCC_0pi_p_axis"));
  FGD1_nubar_CC0pi_th = (TAxis*)(Infile->Get("AxisForCov/FGD1_anti-numuCC_0pi_th_axis"));
  FGD1_nubar_CC1pi_p = (TAxis*)(Infile->Get( "AxisForCov/FGD1_anti-numuCC_1pi_p_axis"));
  FGD1_nubar_CC1pi_th = (TAxis*)(Infile->Get("AxisForCov/FGD1_anti-numuCC_1pi_th_axis"));
  FGD1_nubar_CCOth_p = (TAxis*)(Infile->Get( "AxisForCov/FGD1_anti-numuCC_other_p_axis"));
  FGD1_nubar_CCOth_th = (TAxis*)(Infile->Get("AxisForCov/FGD1_anti-numuCC_other_th_axis"));

  FGD2_nubar_CC0pi_p = (TAxis*)(Infile->Get( "AxisForCov/FGD2_anti-numuCC_0pi_p_axis"));
  FGD2_nubar_CC0pi_th = (TAxis*)(Infile->Get("AxisForCov/FGD2_anti-numuCC_0pi_th_axis"));
  FGD2_nubar_CC1pi_p = (TAxis*)(Infile->Get( "AxisForCov/FGD2_anti-numuCC_1pi_p_axis"));
  FGD2_nubar_CC1pi_th = (TAxis*)(Infile->Get("AxisForCov/FGD2_anti-numuCC_1pi_th_axis"));
  FGD2_nubar_CCOth_p = (TAxis*)(Infile->Get( "AxisForCov/FGD2_anti-numuCC_other_p_axis"));
  FGD2_nubar_CCOth_th = (TAxis*)(Infile->Get("AxisForCov/FGD2_anti-numuCC_other_th_axis"));

  FGD1_nubar_nu_CC0pi_p = (TAxis*)(Infile->Get( "AxisForCov/FGD1_NuMuBkg_CC0pi_in_AntiNu_Mode_p_axis"));
  FGD1_nubar_nu_CC0pi_th = (TAxis*)(Infile->Get("AxisForCov/FGD1_NuMuBkg_CC0pi_in_AntiNu_Mode_th_axis"));
  FGD1_nubar_nu_CC1pi_p = (TAxis*)(Infile->Get( "AxisForCov/FGD1_NuMuBkg_CC1pi_in_AntiNu_Mode_p_axis"));
  FGD1_nubar_nu_CC1pi_th = (TAxis*)(Infile->Get("AxisForCov/FGD1_NuMuBkg_CC1pi_in_AntiNu_Mode_th_axis"));
  FGD1_nubar_nu_CCOth_p = (TAxis*)(Infile->Get( "AxisForCov/FGD1_NuMuBkg_CCother_in_AntiNu_Mode_p_axis"));
  FGD1_nubar_nu_CCOth_th = (TAxis*)(Infile->Get("AxisForCov/FGD1_NuMuBkg_CCother_in_AntiNu_Mode_th_axis"));

  FGD2_nubar_nu_CC0pi_p = (TAxis*)(Infile->Get( "AxisForCov/FGD2_NuMuBkg_CC0pi_in_AntiNu_Mode_p_axis"));
  FGD2_nubar_nu_CC0pi_th = (TAxis*)(Infile->Get("AxisForCov/FGD2_NuMuBkg_CC0pi_in_AntiNu_Mode_th_axis"));
  FGD2_nubar_nu_CC1pi_p = (TAxis*)(Infile->Get( "AxisForCov/FGD2_NuMuBkg_CC1pi_in_AntiNu_Mode_p_axis"));
  FGD2_nubar_nu_CC1pi_th = (TAxis*)(Infile->Get("AxisForCov/FGD2_NuMuBkg_CC1pi_in_AntiNu_Mode_th_axis"));
  FGD2_nubar_nu_CCOth_p = (TAxis*)(Infile->Get( "AxisForCov/FGD2_NuMuBkg_CCother_in_AntiNu_Mode_p_axis"));
  FGD2_nubar_nu_CCOth_th = (TAxis*)(Infile->Get("AxisForCov/FGD2_NuMuBkg_CCother_in_AntiNu_Mode_th_axis"));

  MomVector.push_back(FGD1_CC0pi_p);
  MomVector.push_back(FGD1_CC1pi_p);
  MomVector.push_back(FGD1_CCOth_p);
  MomVector.push_back(FGD2_CC0pi_p);
  MomVector.push_back(FGD2_CC1pi_p);
  MomVector.push_back(FGD2_CCOth_p);
  MomVector.push_back(FGD1_nubar_CC0pi_p);
  MomVector.push_back(FGD1_nubar_CC1pi_p);
  MomVector.push_back(FGD1_nubar_CCOth_p);
  MomVector.push_back(FGD2_nubar_CC0pi_p);
  MomVector.push_back(FGD2_nubar_CC1pi_p);
  MomVector.push_back(FGD2_nubar_CCOth_p);
  MomVector.push_back(FGD1_nubar_nu_CC0pi_p);
  MomVector.push_back(FGD1_nubar_nu_CC1pi_p);
  MomVector.push_back(FGD1_nubar_nu_CCOth_p);
  MomVector.push_back(FGD2_nubar_nu_CC0pi_p);
  MomVector.push_back(FGD2_nubar_nu_CC1pi_p);
  MomVector.push_back(FGD2_nubar_nu_CCOth_p);

  CosVector.push_back(FGD1_CC0pi_th);
  CosVector.push_back(FGD1_CC1pi_th);
  CosVector.push_back(FGD1_CCOth_th);
  CosVector.push_back(FGD2_CC0pi_th);
  CosVector.push_back(FGD2_CC1pi_th);
  CosVector.push_back(FGD2_CCOth_th);
  CosVector.push_back(FGD1_nubar_CC0pi_th);
  CosVector.push_back(FGD1_nubar_CC1pi_th);
  CosVector.push_back(FGD1_nubar_CCOth_th);
  CosVector.push_back(FGD2_nubar_CC0pi_th);
  CosVector.push_back(FGD2_nubar_CC1pi_th);
  CosVector.push_back(FGD2_nubar_CCOth_th);
  CosVector.push_back(FGD1_nubar_nu_CC0pi_th);
  CosVector.push_back(FGD1_nubar_nu_CC1pi_th);
  CosVector.push_back(FGD1_nubar_nu_CCOth_th);
  CosVector.push_back(FGD2_nubar_nu_CC0pi_th);
  CosVector.push_back(FGD2_nubar_nu_CC1pi_th);
  CosVector.push_back(FGD2_nubar_nu_CCOth_th);

  nBinsFGD1CC0pi = FGD1_CC0pi_p->GetNbins() * FGD1_CC0pi_th->GetNbins();
  nBinsFGD1CC1pi = FGD1_CC1pi_p->GetNbins() * FGD1_CC1pi_th->GetNbins();
  nBinsFGD1CCOth = FGD1_CCOth_p->GetNbins() * FGD1_CCOth_th->GetNbins();

  nBinsFGD2CC0pi = FGD2_CC0pi_p->GetNbins() * FGD2_CC0pi_th->GetNbins();
  nBinsFGD2CC1pi = FGD2_CC1pi_p->GetNbins() * FGD2_CC1pi_th->GetNbins();
  nBinsFGD2CCOth = FGD2_CCOth_p->GetNbins() * FGD2_CCOth_th->GetNbins();

  nBinsFGD1CC0pi_nubar = FGD1_nubar_CC0pi_p->GetNbins() * FGD1_nubar_CC0pi_th->GetNbins();
  nBinsFGD1CC1pi_nubar = FGD1_nubar_CC1pi_p->GetNbins() * FGD1_nubar_CC1pi_th->GetNbins();
  nBinsFGD1CCOth_nubar = FGD1_nubar_CCOth_p->GetNbins() * FGD1_nubar_CCOth_th->GetNbins();

  nBinsFGD2CC0pi_nubar = FGD2_nubar_CC0pi_p->GetNbins() * FGD2_nubar_CC0pi_th->GetNbins();
  nBinsFGD2CC1pi_nubar = FGD2_nubar_CC1pi_p->GetNbins() * FGD2_nubar_CC1pi_th->GetNbins();
  nBinsFGD2CCOth_nubar = FGD2_nubar_CCOth_p->GetNbins() * FGD2_nubar_CCOth_th->GetNbins();

  nBinsFGD1CC0pi_nubarnu = FGD1_nubar_nu_CC0pi_p->GetNbins() * FGD1_nubar_nu_CC0pi_th->GetNbins();
  nBinsFGD1CC1pi_nubarnu = FGD1_nubar_nu_CC1pi_p->GetNbins() * FGD1_nubar_nu_CC1pi_th->GetNbins();
  nBinsFGD1CCOth_nubarnu = FGD1_nubar_nu_CCOth_p->GetNbins() * FGD1_nubar_nu_CCOth_th->GetNbins();

  nBinsFGD2CC0pi_nubarnu = FGD2_nubar_nu_CC0pi_p->GetNbins() * FGD2_nubar_nu_CC0pi_th->GetNbins();
  nBinsFGD2CC1pi_nubarnu = FGD2_nubar_nu_CC1pi_p->GetNbins() * FGD2_nubar_nu_CC1pi_th->GetNbins();
  nBinsFGD2CCOth_nubarnu = FGD2_nubar_nu_CCOth_p->GetNbins() * FGD2_nubar_nu_CCOth_th->GetNbins();

  BinVector.push_back(0);
  BinVector.push_back(nBinsFGD1CC0pi);
  BinVector.push_back(nBinsFGD1CC1pi);
  BinVector.push_back(nBinsFGD1CCOth);
  BinVector.push_back(nBinsFGD2CC0pi);
  BinVector.push_back(nBinsFGD2CC1pi);
  BinVector.push_back(nBinsFGD2CCOth);

  BinVector.push_back(nBinsFGD1CC0pi_nubar);
  BinVector.push_back(nBinsFGD1CC1pi_nubar);
  BinVector.push_back(nBinsFGD1CCOth_nubar);
  BinVector.push_back(nBinsFGD2CC0pi_nubar);
  BinVector.push_back(nBinsFGD2CC1pi_nubar);
  BinVector.push_back(nBinsFGD2CCOth_nubar);

  BinVector.push_back(nBinsFGD1CC0pi_nubarnu);
  BinVector.push_back(nBinsFGD1CC1pi_nubarnu);
  BinVector.push_back(nBinsFGD1CCOth_nubarnu);
  BinVector.push_back(nBinsFGD2CC0pi_nubarnu);
  BinVector.push_back(nBinsFGD2CC1pi_nubarnu);
  //BinVector.push_back(nBinsFGD2CCOth_nubarnu);

  // Sum up to make additive
  int oldcontent = 0;
  for (std::vector<int>::iterator it = BinVector.begin(); it!=BinVector.end(); ++it) {
    (*it) += oldcontent;
    oldcontent = *it;
  }
}

void FindSamplePmuCosmu(int bin, int &selection, double &pmu, double &cosmu) {
  // Find the sample
  int sample = 0;
  while (bin >= BinVector[sample+1]) {
    sample++;
  }
  selection = sample;

  int nbinsp = MomVector[sample]->GetNbins();
  int nbinsth = CosVector[sample]->GetNbins();
  int bincount = BinVector[sample];

  // Get the local index
  int LocalIndex = bin - bincount;
  int MomBin = LocalIndex % nbinsp + 1;
  int CosBin = LocalIndex / nbinsp + 1;
  pmu = MomVector[sample]->GetBinLowEdge(MomBin);
  cosmu = CosVector[sample]->GetBinLowEdge(CosBin);
}

// Check if we can merge along pmu
bool CosMerge(int GlobalBin, int sample, double pmu, double cosmu) {

  int nbinp  = MomVector[sample]->GetNbins();
  int CosBin = CosVector[sample]->FindBin(cosmu);

  // The lowest momentum bin in this slice of pmu
  int LowIndex = BinVector[sample] + nbinp*(CosBin-1);
  // The highest momentum bin in this slice of theta
  int TopIndex = LowIndex+nbinp-1;

  //std::cout << "Looping from " << LowIndex << "-" << TopIndex << std::endl;
  for (int j = LowIndex; j < TopIndex; ++j) {
    int GlobalIndex = j;
    // Get the variance in the current binning
    //int sel;
    //double testpmu, testcosmu;
    //FindSamplePmuCosmu(GlobalIndex, sel, testpmu, testcosmu);
    //std::cout << "    GlobalIndex " << GlobalIndex << " sel: " << sel << " pmu: " << testpmu << " cosmu: " << testcosmu << std::endl;
    //FindSamplePmuCosmu(GlobalIndex+nbinp, sel, testpmu, testcosmu);
    //std::cout << "    GlobalIndex " << GlobalIndex+nbinp << " sel: " << sel << " pmu: " << testpmu << " cosmu: " << testcosmu << std::endl;

    // Get the variance in the next bins
    double oldsigma = (Covariance->GetBinContent(j, j));
    double sigma = (Covariance->GetBinContent(j+nbinp, j+nbinp));

    if (fabs(sigma/oldsigma-1) > threshold  && Nominal->GetBinContent(j) > 1 && sqrt(sigma)*Nominal->GetBinContent(j) > 1 && sqrt(sigma) > threshold) {
    //if (fabs(sigma/oldsigma-1) > threshold && sqrt(sigma) > 0.05 && Nominal->GetBinContent(j) > 1 && sqrt(sigma)*Nominal->GetBinContent(j) > 1) {
    //if (fabs(sigma/oldsigma-1) > threshold && sqrt(sigma)*Nominal->GetBinContent(j) > 0.5) {
      return false;
    }
  }
  return true;
}

// Check if we can merge along pmu
bool PmuMerge(int GlobalBin, int sample, double pmu, double cosmu) {

  int PmuBin = MomVector[sample]->FindBin(pmu);
  int nbinp = MomVector[sample]->GetNbins();
  int nbinth = CosVector[sample]->GetNbins();

  int LocalIndex = GlobalBin;
  // The lowest cosmu bin in this slice of pmu
  int LowLocalIndex = BinVector[sample]+(PmuBin%nbinp)-1;
  // The highest momentum bin in this slice of theta
  int TopLocalIndex = LowLocalIndex + (nbinth-1)*nbinp;

  //std::cout << "Looping from " << LowLocalIndex << "-" << TopLocalIndex << std::endl;
  for (int j = LowLocalIndex; j < TopLocalIndex; j = j + nbinp) {
    int GlobalIndex = j;
    //int sel;
    //double testpmu, testcosmu;
    //FindSamplePmuCosmu(GlobalIndex, sel, testpmu, testcosmu);
    //std::cout << "     GlobalIndex " << GlobalIndex << " sel: " << sel << " pmu: " << testpmu << " cosmu: " << testcosmu << std::endl;
    //FindSamplePmuCosmu(GlobalIndex+1, sel, testpmu, testcosmu);
    //std::cout << "     GlobalIndex " << GlobalIndex+1 << " sel: " << sel << " pmu: " << testpmu << " cosmu: " << testcosmu << std::endl;
    double oldsigma = (Covariance->GetBinContent(j,j));
    double sigma = (Covariance->GetBinContent(j+1, j+1));
    if (fabs(sigma/oldsigma-1) > threshold  && Nominal->GetBinContent(j) > 1 && sqrt(sigma)*Nominal->GetBinContent(j) > 1 && sqrt(sigma) > threshold) {
    //if (fabs(sigma/oldsigma-1) > threshold && sqrt(sigma) > 0.05 && Nominal->GetBinContent(j) > 1 && sqrt(sigma)*Nominal->GetBinContent(j) > 1) {
    //if (fabs(sigma/oldsigma-1) > threshold && sqrt(sigma)*Nominal->GetBinContent(j) > 0.5) {
      return false;
    }
  }
  return true;
}

int main(int argc, char** argv) {
  if (argc != 3) {
    std::cerr << "Need two arguments, input file and tolerance" << std::endl;
    std::cerr << "I got " << argc << std::endl;
    for (int i = 0; i < argc; ++i) {
      std::cout << argv[i] << " ";
    }
    std::cout << std::endl;
    throw;
  }

  threshold = std::atof(argv[2]);
  std::string File = argv[1];
  MergeBins(File);
}

