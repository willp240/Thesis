#include <iostream>
#include <fstream>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TH2D.h"
#include "TMatrixTSym.h"
#include "TVectorT.h"
#include "TLine.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TStyle.h"
#include "TColor.h"
#include "TFitResultPtr.h"
#include "TDirectory.h"
#include "TROOT.h"
#include "TKey.h"
#include "THn.h"

#include "BinningDefinition.hxx"
#include "ND280AnalysisUtils.hxx"
#include "Parameters.hxx"

bool ApplyNIWG = false;
bool Apply1p1h = false;
const int nThrows = 2000;

void MakeCovariance(char* MCFilesListStr, char* NIWGFilesListStr,
    char* Cov1P1HStr,     char* OutFileStr){

  ND::params().SetReadParamOverrideFilePointPassed();
  BANFF::BinningDefinition& bd = BANFF::BinningDefinition::Get();

  std::vector<std::string> name_vector = bd.GetListOfBins();
  std::vector<std::string> name_vector_det = bd.GetListOfBins_Det();
  /*
  for (std::vector<std::string>::iterator it = name_vector.begin(); it != name_vector.end(); ++it) {
    std::cout << *it << std::endl;
  }
  for (std::vector<std::string>::iterator it = name_vector_det.begin(); it != name_vector_det.end(); ++it) {
    std::cout << *it << std::endl;
  }
  */

  // Count the number of active samples
  int nCountedSamples = 0;
  for (int i = 0; i < SampleId::kNSamples; ++i) {
    if (bd.IsActiveSample(SampleId::SampleEnum(i))) nCountedSamples++;
  }

  std::cout << "Found " << nCountedSamples << " enabled samples in BinningDefinition" << std::endl;

  TFile* fake = NULL;
  TH1D* fake_cov = NULL;

  int nSamples = 14;

  const int nSamplesMultiTrack = 14;
  const int nSamplesMultiPi = 18;

  if (nSamplesMultiTrack != nCountedSamples && nSamplesMultiPi != nCountedSamples) {
    std::cerr << "Hard-coded number of samples does not match found samples in BinningDefinition" << std::endl;
    std::cerr << "nSamples (hard-coded): " << nSamples << std::endl;
    std::cerr << "nFoundSamples: " << nCountedSamples << std::endl;
    throw;
  }
  nSamples = nCountedSamples;

  TH2D** Nominal_Plots = new TH2D*[nSamples];
  TH2D** Fake_Plots = new TH2D*[nSamples];

  if(Apply1p1h){
    std::cout << "Fake1P1H File " << Cov1P1HStr << std::endl;
    fake = new TFile(Cov1P1HStr);

    THnD* Nominal_Plots_In[nSamples];
    THnD* Fake_Plots_In[nSamples];

    Fake_Plots_In[0] =  (THnD*) fake->Get("FGD1_numuCC_0pi_fakedata");
    Fake_Plots_In[1] =  (THnD*) fake->Get("FGD1_numuCC_1pi_fakedata");
    Fake_Plots_In[2] =  (THnD*) fake->Get("FGD1_numuCC_other_fakedata");
    Fake_Plots_In[3] =  (THnD*) fake->Get("FGD2_numuCC_0pi_fakedata");
    Fake_Plots_In[4] =  (THnD*) fake->Get("FGD2_numuCC_1pi_fakedata");
    Fake_Plots_In[5] =  (THnD*) fake->Get("FGD2_numuCC_other_fakedata");

    if (nSamples == 14) {
      Fake_Plots_In[6] =  (THnD*) fake->Get("FGD1_anti-numuCC_QE_fakedata");
      Fake_Plots_In[7] =  (THnD*) fake->Get("FGD1_anti-numuCC_nQE_fakedata");
      Fake_Plots_In[8] =  (THnD*) fake->Get("FGD2_anti-numuCC_1_track_fakedata");
      Fake_Plots_In[9] =  (THnD*) fake->Get("FGD2_anti-numuCC_N_tracks_fakedata");
      Fake_Plots_In[10] = (THnD*) fake->Get("FGD1_NuMuBkg_CCQE_in_AntiNu_Mode__fakedata");
      Fake_Plots_In[11] = (THnD*) fake->Get("FGD1_NuMuBkg_CCnQE_in_AntiNu_Mode_fakedata");
      Fake_Plots_In[12] = (THnD*) fake->Get("FGD2_NuMuBkg_CCQE_in_AntiNu_Mode__fakedata");
      Fake_Plots_In[13] = (THnD*) fake->Get("FGD2_NuMuBkg_CCnQE_in_AntiNu_Mode_fakedata");
    } else {
      Fake_Plots_In[6] =  (THnD*) fake->Get("FGD1_anti-numuCC_0pi_fakedata");
      Fake_Plots_In[7] =  (THnD*) fake->Get("FGD1_anti-numuCC_1pi_fakedata");
      Fake_Plots_In[8] =  (THnD*) fake->Get("FGD1_anti-numuCC_other_fakedata");

      Fake_Plots_In[9] =  (THnD*) fake->Get("FGD2_anti-numuCC_0pi_fakedata");
      Fake_Plots_In[10] = (THnD*) fake->Get("FGD2_anti-numuCC_1pi_fakedata");
      Fake_Plots_In[11] = (THnD*) fake->Get("FGD2_anti-numuCC_other_fakedata");

      Fake_Plots_In[12] = (THnD*) fake->Get("FGD1_NuMuBkg_CC0pi_in_AntiNu_Mode_fakedata");
      Fake_Plots_In[13] = (THnD*) fake->Get("FGD1_NuMuBkg_CC1pi_in_AntiNu_Mode_fakedata");
      Fake_Plots_In[14] = (THnD*) fake->Get("FGD1_NuMuBkg_CCother_in_AntiNu_Mode_fakedata");

      Fake_Plots_In[15] = (THnD*) fake->Get("FGD2_NuMuBkg_CC0pi_in_AntiNu_Mode_fakedata");
      Fake_Plots_In[16] = (THnD*) fake->Get("FGD2_NuMuBkg_CC1pi_in_AntiNu_Mode_fakedata");
      Fake_Plots_In[17] = (THnD*) fake->Get("FGD2_NuMuBkg_CCother_in_AntiNu_Mode_fakedata");
    }


    Nominal_Plots_In[0] =  (THnD*) fake->Get("FGD1_numuCC_0pi_prefit");
    Nominal_Plots_In[1] =  (THnD*) fake->Get("FGD1_numuCC_1pi_prefit");
    Nominal_Plots_In[2] =  (THnD*) fake->Get("FGD1_numuCC_other_prefit");
    Nominal_Plots_In[3] =  (THnD*) fake->Get("FGD2_numuCC_0pi_prefit");
    Nominal_Plots_In[4] =  (THnD*) fake->Get("FGD2_numuCC_1pi_prefit");
    Nominal_Plots_In[5] =  (THnD*) fake->Get("FGD2_numuCC_other_prefit");

    if (nSamples == 14) {
      Nominal_Plots_In[6] =  (THnD*) fake->Get("FGD1_anti-numuCC_QE_prefit");
      Nominal_Plots_In[7] =  (THnD*) fake->Get("FGD1_anti-numuCC_nQE_prefit");
      Nominal_Plots_In[8] =  (THnD*) fake->Get("FGD2_anti-numuCC_1_track_prefit");
      Nominal_Plots_In[9] =  (THnD*) fake->Get("FGD2_anti-numuCC_N_tracks_prefit");
      Nominal_Plots_In[10] = (THnD*) fake->Get("FGD1_NuMuBkg_CCQE_in_AntiNu_Mode__prefit");
      Nominal_Plots_In[11] = (THnD*) fake->Get("FGD1_NuMuBkg_CCnQE_in_AntiNu_Mode_prefit");
      Nominal_Plots_In[12] = (THnD*) fake->Get("FGD2_NuMuBkg_CCQE_in_AntiNu_Mode__prefit");
      Nominal_Plots_In[13] = (THnD*) fake->Get("FGD2_NuMuBkg_CCnQE_in_AntiNu_Mode_prefit");
    } else {
      Nominal_Plots_In[6] =  (THnD*) fake->Get("FGD1_anti-numuCC_0pi_prefit");
      Nominal_Plots_In[7] =  (THnD*) fake->Get("FGD1_anti-numuCC_1pi_prefit");
      Nominal_Plots_In[8] =  (THnD*) fake->Get("FGD1_anti-numuCC_other_prefit");

      Nominal_Plots_In[9] =  (THnD*) fake->Get("FGD2_anti-numuCC_0pi_prefit");
      Nominal_Plots_In[10] = (THnD*) fake->Get("FGD2_anti-numuCC_1pi_prefit");
      Nominal_Plots_In[11] = (THnD*) fake->Get("FGD2_anti-numuCC_other_prefit");

      Nominal_Plots_In[12] = (THnD*) fake->Get("FGD1_NuMuBkg_CC0pi_in_AntiNu_Mode_prefit");
      Nominal_Plots_In[13] = (THnD*) fake->Get("FGD1_NuMuBkg_CC1pi_in_AntiNu_Mode_prefit");
      Nominal_Plots_In[14] = (THnD*) fake->Get("FGD1_NuMuBkg_CCother_in_AntiNu_Mode_prefit");

      Nominal_Plots_In[15] = (THnD*) fake->Get("FGD2_NuMuBkg_CC0pi_in_AntiNu_Mode_prefit");
      Nominal_Plots_In[16] = (THnD*) fake->Get("FGD2_NuMuBkg_CC1pi_in_AntiNu_Mode_prefit");
      Nominal_Plots_In[17] = (THnD*) fake->Get("FGD2_NuMuBkg_CCother_in_AntiNu_Mode_prefit");
    }

    std::cout << "Got fake data from " << Cov1P1HStr << std::endl;

    for(int i = 0; i < nSamples; ++i){
      Fake_Plots[i] = (TH2D*) Fake_Plots_In[i]->Projection(1,0)->Clone(Form("%s_%i_2D",Fake_Plots_In[i]->GetName(), i));
      Nominal_Plots[i] = (TH2D*) Nominal_Plots_In[i]->Projection(1,0)->Clone(Form("%s_%i_2D",Nominal_Plots_In[i]->GetName(), i));
    }

    // Write the fake-data and nominal distributions to file
    std::string titlename = OutFileStr;
    titlename += "_fakedata.root";
    TFile *tempfile = new TFile(titlename.c_str(), "recreate");
    for (int i = 0; i < nSamples; ++i) {
      Fake_Plots[i]->Write();
      Nominal_Plots[i]->Write();
    }
    tempfile->Close();
  }

  if(Apply1p1h){
    fake_cov = new TH1D("fake_error","fake_error",bd.GetNbins(), 0, bd.GetNbins());
    for(int i = 0; i < nSamples; ++i){
      SampleId::SampleEnum sample = SampleId::kUnassigned;
      switch(i) {
        case 0:
          sample = SampleId::kFGD1NuMuCC0Pi;
          break;
        case 1:
          sample = SampleId::kFGD1NuMuCC1Pi;
          break;
        case 2:
          sample = SampleId::kFGD1NuMuCCOther;
          break;
        case 3:
          sample = SampleId::kFGD2NuMuCC0Pi;
          break;
        case 4:
          sample = SampleId::kFGD2NuMuCC1Pi;
          break;
        case 5:
          sample = SampleId::kFGD2NuMuCCOther;
          break;
        case 6:
          if (nSamples == 14) {
            sample = SampleId::kFGD1AntiNuMuCC1Track;
          } else {
            sample = SampleId::kFGD1AntiNuMuCC0Pi;
          }
          break;
        case 7:
          if (nSamples == 14) {
            sample = SampleId::kFGD1AntiNuMuCCNTracks;
          } else {
            sample = SampleId::kFGD1AntiNuMuCC1Pi;
          }
          break;
        case 8:
          if (nSamples == 14) {
            sample = SampleId::kFGD2AntiNuMuCC1Track;
          } else {
            sample = SampleId::kFGD1AntiNuMuCCOther;
          }
          break;
        case 9:
          if (nSamples == 14) {
            sample = SampleId::kFGD2AntiNuMuCCNTracks;
          } else {
            sample = SampleId::kFGD2AntiNuMuCC0Pi;
          }
          break;
        case 10:
          if (nSamples == 14) {
            sample = SampleId::kFGD1NuMuBkgInAntiNuModeCC1Track;
          } else {
            sample = SampleId::kFGD2AntiNuMuCC1Pi;
          }
          break;
        case 11:
          if (nSamples == 14) {
            sample = SampleId::kFGD1NuMuBkgInAntiNuModeCCNTracks;
          } else {
            sample = SampleId::kFGD2AntiNuMuCCOther;
          }
          break;
        case 12:
          if (nSamples == 14) {
            sample = SampleId::kFGD2NuMuBkgInAntiNuModeCC1Track;
          } else {
            sample = SampleId::kFGD1NuMuBkgInAntiNuModeCC0Pi;
          }
          break;
        case 13:
          if (nSamples == 14) {
            sample = SampleId::kFGD2NuMuBkgInAntiNuModeCCNTracks;
          } else {
            sample = SampleId::kFGD1NuMuBkgInAntiNuModeCC1Pi;
          }
          break;
        case 14:
          sample = SampleId::kFGD1NuMuBkgInAntiNuModeCCOther;
          break;
        case 15:
          sample = SampleId::kFGD2NuMuBkgInAntiNuModeCC0Pi;
          break;
        case 16:
          sample = SampleId::kFGD2NuMuBkgInAntiNuModeCC1Pi;
          break;
        case 17:
          sample = SampleId::kFGD2NuMuBkgInAntiNuModeCCOther;
          break;
        default:
          sample = SampleId::kUnassigned;
          break;
      }

      for(int j = 1; j <= Nominal_Plots[i]->GetNbinsX(); ++j){
        for(int k = 1; k <= Nominal_Plots[i]->GetNbinsY(); ++k){
          int global_bin = bd.GetGlobalBinMatrixMomCos(sample, Nominal_Plots[i]->GetXaxis()->GetBinCenter(j), Nominal_Plots[i]->GetYaxis()->GetBinCenter(k));
          if (global_bin == -1) {
            std::cerr << "Couldn't find global bin for " << SampleId::ConvertSample(sample) << std::endl;
            std::cerr << "With " << Nominal_Plots[i]->GetXaxis()->GetBinCenter(j) << std::endl;
            std::cerr << "With " << Nominal_Plots[i]->GetYaxis()->GetBinCenter(k) << std::endl;
            throw;
          }
          fake_cov->Fill(global_bin + 1, (Fake_Plots[i]->GetBinContent(j,k) - Nominal_Plots[i]->GetBinContent(j,k))/Nominal_Plots[i]->GetBinContent(j,k));
        }
      }
    }
  } // Have now made the correlation matrix

  ifstream MCFilesList;
  MCFilesList.open(MCFilesListStr);

  std::map<std::string, TChain*> MCChains;

  if(MCFilesList.is_open()){
    std::string line;
    while(getline(MCFilesList, line)){
      if(MCChains.size() == 0){
        std::cout << "Reading list of systematics trees from file" << std::endl;
        TFile temp(line.c_str(), "READ");

        TIter next(temp.GetListOfKeys());
        TKey *key;
        while ((key = (TKey*)next())) {
          TClass *cl = gROOT->GetClass(key->GetClassName());
          if (!cl->InheritsFrom("TTree")) continue;
          TTree *h = (TTree*)key->ReadObj();
          std::string name = h->GetName();
          // Skip the NRooTrackerVtx chain
          if(name.find("RooTrackerVtx")!=std::string::npos) continue;
          std::cout << "Adding syst " << name << " to TChain map" << std::endl;
          MCChains.insert(std::pair<std::string, TChain*>(name, new TChain(name.c_str())));
        } 
      }
      std::cout << "Adding the file " << line << " to the TChain" << std::endl;
      for(std::map<std::string, TChain*>::iterator it = MCChains.begin(); it != MCChains.end(); ++it){
        it->second->Add(line.c_str());
      }
    }
    MCFilesList.close();
  }else{
    std::cerr << "List of MC files not open!!" << std::endl;
    throw;
  }

  ifstream NIWGFilesList;
  NIWGFilesList.open(NIWGFilesListStr);
  TChain* NIWGTune = NULL;

  if(ApplyNIWG){
    if(NIWGFilesList.is_open()){
      NIWGTune = new TChain("NIWGWeightTree");
      std::string line;
      while(getline(NIWGFilesList, line)){
        std::cout << "NIWG tune File " << line << std::endl;
        NIWGTune->Add(line.c_str());
      }
      NIWGFilesList.close();
    }else{
      std::cerr << "List of MC files not open!!" << std::endl;
      throw;
    }
  }

  std::cout << "Added all the Psyche and NIWG files to the file" << std::endl;

  int    Run            ;
  int    SubRun         ;
  int    EventNumber    ;
  int    SelectionNom   ;
  double TrueEnuNom     ;
  int    TrueNuPDGNom   ;
  int    TrueVertexIDNom;
  double LeptonMomNom   ;
  double LeptonCosNom   ;
  double WeightNom      ;
  double FluxWeightNom  ;
  double NIWGWeightNom  ;
  int    nToys          ;
  int    Toy            [nThrows];
  int    TrueVertexIDToy[nThrows];
  int    SelectionToy   [nThrows];
  double TrueEnuToy     [nThrows];
  int    TrueNuPDGToy   [nThrows];
  double LeptonMomToy   [nThrows];
  double LeptonCosToy   [nThrows];
  double WeightToy      [nThrows];
  double FluxWeightToy  [nThrows];
  double NIWGWeightToy  [nThrows];


  TChain* MCChain;

  for(std::map<std::string, TChain*>::iterator it = MCChains.begin(); it != MCChains.end(); ++it){

    MCChain = it->second;
    MCChain->SetBranchAddress("Run"            , &Run            );
    MCChain->SetBranchAddress("SubRun"         , &SubRun         );
    MCChain->SetBranchAddress("EventNumber"    , &EventNumber    );
    MCChain->SetBranchAddress("SelectionNom"   , &SelectionNom   );
    MCChain->SetBranchAddress("TrueEnuNom"     , &TrueEnuNom     );
    MCChain->SetBranchAddress("TrueNuPDGNom"   , &TrueNuPDGNom   );
    MCChain->SetBranchAddress("TrueVertexIDNom", &TrueVertexIDNom);
    MCChain->SetBranchAddress("LeptonMomNom"   , &LeptonMomNom   );
    MCChain->SetBranchAddress("LeptonCosNom"   , &LeptonCosNom   );
    MCChain->SetBranchAddress("WeightNom"      , &WeightNom      );
    MCChain->SetBranchAddress("FluxWeightNom"  , &FluxWeightNom  );
    MCChain->SetBranchAddress("nToys"          , &nToys          );
    MCChain->SetBranchAddress("Toy"            ,  Toy            );
    MCChain->SetBranchAddress("TrueVertexIDToy",  TrueVertexIDToy);
    MCChain->SetBranchAddress("SelectionToy"   ,  SelectionToy   );
    MCChain->SetBranchAddress("TrueEnuToy"     ,  TrueEnuToy     );
    MCChain->SetBranchAddress("TrueNuPDGToy"   ,  TrueNuPDGToy   );
    MCChain->SetBranchAddress("LeptonMomToy"   ,  LeptonMomToy   );
    MCChain->SetBranchAddress("LeptonCosToy"   ,  LeptonCosToy   );
    MCChain->SetBranchAddress("WeightToy"      ,  WeightToy      );
    MCChain->SetBranchAddress("FluxWeightToy"  ,  FluxWeightToy  );
    if(ApplyNIWG){
      NIWGTune->SetBranchAddress("T2KRW_NIWGWeightNom",    &NIWGWeightNom);
      NIWGTune->SetBranchAddress("T2KRW_NIWGWeightToy",     NIWGWeightToy);
      MCChain->AddFriend(NIWGTune);

      /*
      // Create some test projections of Enu from the full MCChain
      TH2D* check1 = new TH2D("EnuCheck", "EnuCheck;Enu;Enu", 1000, 0, 5000, 1000, 0, 5000);
      MCChain->Project("EnuCheck", "TrueEnuToy:T2KRW_EnuToy");

      TH2D* check2 = new TH2D("EnuCheck2", "EnuCheck2;Enu;Enu", 1000, 0, 5000, 1000, 0, 5000);
      MCChain->Project("EnuCheck2", "TrueEnuNom:T2KRW_EnuNom");

      TCanvas* c = new TCanvas();
      check1->Draw();
      c->SaveAs("check1.png");

      TCanvas* b = new TCanvas();
      check2->Draw();
      b->SaveAs("check2.png");
      */
    }else{
      std::cout << "WARNING!! NOT APPLYING THE NIWG WEIGHTING" << std::endl;
      for(int i = 0;i < nThrows; i++) NIWGWeightToy[i] = 1;
      NIWGWeightNom = 1.;
    }

    std::cout << "Creating the file " << OutFileStr << ".root" << std::endl;
    TFile* OutputFile = new TFile(Form("%s_%s.root",OutFileStr,(it->first).c_str()), "RECREATE");

    int nBinsRedu = bd.GetNbins_Det();
    int nBinsFull = bd.GetNbins();

    std::cout << "We have " << nBinsRedu <<  " bins in the covariance matrix" << std::endl;

    TMatrixTSym<double> matrix(nBinsRedu);
    TVectorT   <double> vector(nBinsRedu);
    TH2D* corr               = new TH2D("Total_Correlation_Matrix", "Total Correlation Matrix",                      nBinsRedu, 0, nBinsRedu, nBinsRedu, 0, nBinsRedu);
    TH2D* corr_nofake        = new TH2D("Total_Correlation_Matrix_No1p1h", "Total Correlation Matrix without 1p1h error",                      nBinsRedu, 0, nBinsRedu, nBinsRedu, 0, nBinsRedu);
    TH2D* cov                = new TH2D("Total_Covariance_Matrix",  "Total Covariance Matrix",                       nBinsRedu, 0, nBinsRedu, nBinsRedu, 0, nBinsRedu);
    TH2D* cov_var            = new TH2D("Covariance_Matrix_NoMCStats",    "Covariance Matrix without MC stats error",               nBinsRedu, 0, nBinsRedu, nBinsRedu, 0, nBinsRedu);
    TH2D* corr_var           = new TH2D("Correlation_Matrix_NoMCStats",   "Correlation Matrix without MC stats error",              nBinsRedu, 0, nBinsRedu, nBinsRedu, 0, nBinsRedu);
    TH2D* cov_zero           = new TH2D("Correlation_Matrix_1p1h",  "Correlation Matrix for 1p1h error",             nBinsRedu, 0, nBinsRedu, nBinsRedu, 0, nBinsRedu);
    TH2D* cov_nofake         = new TH2D("Covariance_Matrix_NoFake", "Covariance Matrix without 1p1h error",          nBinsRedu, 0, nBinsRedu, nBinsRedu, 0, nBinsRedu);
    TH2D* fake_error_cov     = new TH2D("1p1h_fake_Error_cov",      "1p1h systematic error",                         nBinsRedu, 0, nBinsRedu, nBinsRedu, 0, nBinsRedu);
    TH1D* fake_error_redu    = new TH1D("1p1h_fake_Error_redu",     "1p1h systematic error",                         nBinsRedu, 0, nBinsRedu);
    TH1D* nominal            = new TH1D("Nominal",                  "Nominal Events",                                nBinsRedu, 0, nBinsRedu);
    TH1D* mc_sys_error       = new TH1D("MC_Sys_Error",             "MC systematic error",                           nBinsRedu, 0, nBinsRedu);
    TH1D* mc_stat_error_redu = new TH1D("MC_Stat_Error_Redu",       "MC statistical error",                          nBinsRedu, 0, nBinsRedu); 
    TH1D* var_mean           = new TH1D("Varied_Mean",              "Varied Mean",                                   nBinsRedu, 0, nBinsRedu);
    TH1D* zero_point         = new TH1D("Zero_Point",               "Zero Point",                                    nBinsRedu, 0, nBinsRedu);
    TH1D* var_mean_offset    = new TH1D("Varied_Mean_over_Nominal", "Varied Mean normalised by nominal",             nBinsRedu, 0, nBinsRedu);
    TH1D* number_events_redu = new TH1D("number_events",            "Number of Events",                              nBinsRedu, 0, nBinsRedu);
    TH1D* GOF                = new TH1D("GOF",                      "Chi-square between gaussian and bin variation", nBinsRedu, 0, nBinsRedu);
    TH1D* GOFnDOF            = new TH1D("GOFnDOF",                      "Chi-square/bins between gaussian and bin variation", nBinsRedu, 0, nBinsRedu);
    GOF->GetYaxis()->SetTitle("#chi^{2}");
    GOFnDOF->GetYaxis()->SetTitle("#chi^{2}_{red}");
    TH1D* DiffMean           = new TH1D("DiffMean",                      "Difference in Mean between gaussian and binned relative nominal", nBinsRedu, 0, nBinsRedu);
    DiffMean->GetYaxis()->SetTitle("#mu_{Gauss} - #mu_{Binned}");
    TH1D* DiffError          = new TH1D("DiffError",                      "Difference in Error between gaussian and binned relative nominal", nBinsRedu, 0, nBinsRedu);
    DiffError->GetYaxis()->SetTitle("#sigma_{Gauss} - #sigma_{Binned}");
    TH1D* varied[nThrows];
    for (int i = 0; i < nThrows; ++i){
      varied[i] = new TH1D(Form("Varied_Events_%d", i), Form("Varied Events %d", i), nBinsRedu, 0, nBinsRedu);
    }

    TH1D* mc_stat_error_full = new TH1D("MC_Stat_Error_full",   "MC statistical error",   nBinsFull, 0, nBinsFull);
    TH1D* mc_stat_nom        = new TH1D("MC_Stat_Nom",          "MC statistical nominal", nBinsFull, 0, nBinsFull);
    TH1D* number_events_full = new TH1D("number_events_full",   "Number of Events",       nBinsFull, 0, nBinsFull);
    TH1D* fake_error_full    = new TH1D("1p1h_fake_Error_full", "1p1h systematic error",  nBinsFull, 0, nBinsFull);

    if(fake_cov){
      for(int i = 1; i <= fake_cov->GetNbinsX(); ++i){
        fake_error_full->SetBinContent(i,fake_cov->GetBinContent(i));
      }
    }

    TH1D* Selection = new TH1D("Selection", "Selection;SampleId;N_{events} MC", SampleId::kNSamples, 0, SampleId::kNSamples);

    for(int iBin = 1; iBin < nBinsRedu+1; ++iBin){
      for(int jBin = 1; jBin < nBinsRedu+1; ++jBin){
        cov       ->SetBinContent(iBin, jBin, 0);
        cov_var   ->SetBinContent(iBin, jBin, 0);
        cov_zero  ->SetBinContent(iBin, jBin, 0);
        cov_nofake->SetBinContent(iBin, jBin, 0);
      }
      nominal           ->SetBinContent(iBin, 0);
      mc_stat_error_redu->SetBinContent(iBin, 0);
      mc_sys_error      ->SetBinContent(iBin, 0);
      var_mean          ->SetBinContent(iBin, 0);
      var_mean_offset   ->SetBinContent(iBin, 0);
      zero_point        ->SetBinContent(iBin, 0);

      for(int iToy = 0; iToy < nThrows; ++iToy){
        varied[iToy]->SetBinContent(iBin, 0);
      }
    }

    int nEntries = MCChain->GetEntries();

    double POT_Weight = 1;
    int countwidth=double(nEntries)/double(20)+1;

    for(int i = 0; i < nEntries; ++i){
      MCChain->GetEntry(i);
      //std::cout << MCChain->GetFile()->GetName() << std::endl;
      if(nThrows != nToys){
        std::cerr << "nThrows should be equal to nToys, change this in the beginning MakeND280Cov and recompile" << std::endl;
        throw;
      }
      if (i % countwidth == 0) {
        std::cout << "On event " << i << "/" << nEntries << " (" << int(100. * ((double) i) / double(nEntries)) << "%)" << std::endl;
      }
      bool isSand = 0;
      Int_t runp = anaUtils::GetRunPeriod(Run, SubRun);

      if(anaUtils::GetSandMode(Run) >= 0)
        isSand = true;
      switch (runp){
        case 0://run1 water
          POT_Weight = 0.018057359;
          break;
        case 1://run2 water
          POT_Weight = 4.33765e+19/1.20341e+21;
          if(isSand) POT_Weight = 4.33765e+19/4.00035e+19;
          break;
        case 2://run2 air
          POT_Weight = 3.59337e+19/9.23937e+20;
          if(isSand) POT_Weight = 3.59337e+19/3.7132e+19;
          break;       
        case 3://run3b air -- edited from 4.44873
          POT_Weight = 2.1705e+19/4.47864e+20;
          if(isSand) POT_Weight = 2.1705e+19/2.35053e+19;
          break;       
        case 4://run3c air -- edited from 2.633027
          POT_Weight = 1.36398e+20/2.63227e+21;
          if(isSand) POT_Weight = 1.36398e+20/1.31337e+20;
          break;
        case 5://run4 water
          POT_Weight = 1.64277e+20/2.26216e+21;	
          if(isSand) POT_Weight = 1.64277e+20/1.59801e+20;
          break;
        case 6://run4 air
          POT_Weight = 1.78271e+20/3.4996e+21;
          if(isSand) POT_Weight = 1.78271e+20/1.74125e+20;
          break;
        case 8://run5 antinu-water
          POT_Weight = 4.3468e+19/2.29627e+21;
          if(isSand) POT_Weight = 4.3468e+19/9.07403e+19;
          break;
        case 9://run6b antinu-air -- edited from 3.41655 sand
          POT_Weight = 1.27301e+20/1.4174e+21;
          if(isSand) POT_Weight = 1.27301e+20/2.59187e+20;
          break;
        case 10://run6c antinu-air
          POT_Weight = 5.07819e+19/5.27562e+20;
          if(isSand) POT_Weight = 5.07819e+19/1.04626e+20;
          break;
        case 11://run6d antinu-air
          POT_Weight = 7.75302e+19/6.883e+20;
          if(isSand) POT_Weight = 7.75302e+19/1.58059e+20;
          break;
        case 12://run6e antinu-air -- edited from 1.75435 sand
          POT_Weight = 8.51429e+19/8.59439e+20;
          if(isSand) POT_Weight = 8.51429e+19/1.72691e+20;
          break;
        case 15://run7b antinu-water
          POT_Weight = 2.43683e+20/3.37059e+21;
          if(isSand) POT_Weight = 2.43683e+20/5.03961e+20;
          break;
        case 17://run8 nu-water
          POT_Weight = 1.58053e+20/2.64115e+21;
          if(isSand) POT_Weight = 1.58053e+20/1.61263e+20;
          break;
        case 18://run8 nu-air -- edited from 4.04241 sand
          POT_Weight = 4.14909e+20/3.63054e+21;
          if(isSand) POT_Weight = 4.14909e+20/4.01875e+20;
          break;
        default:
          std::cout << "NOT FOUND?!?!" << std::endl;
          throw;
      }
      if(isSand){
        for(int i = 0;i < nThrows; i++) NIWGWeightToy[i] = 1;
        NIWGWeightNom = 1.;
      }

      if(POT_Weight > 2){
        std::cout << "BAD POT WEIGHT" << std::endl;
        throw;
      }

      if(bd.IsActiveSample(SelectionNom)){
        Selection->Fill(SelectionNom);
        int BinRedu = bd.GetGlobalBinMatrixMomCos_Det((SampleId::SampleEnum)SelectionNom, LeptonMomNom, LeptonCosNom) + 1;
        if(BinRedu > nBinsRedu){
          std::cout << "BinRedu > nBinsRedu" << std::endl;
          throw;
        }
        if (BinRedu == -1) {
          std::cout << "Couldn't find a bin in GetGlobalBinmatrixMomCos_Det for " << SelectionNom << " (" << SampleId::ConvertSample(SampleId::SampleEnum(SelectionNom)) << ")" << std::endl;
          throw;
        }
        int BinFull = bd.GetGlobalBinMatrixMomCos((SampleId::SampleEnum)SelectionNom,
            LeptonMomNom, LeptonCosNom) + 1;

        //std::cout << SampleId::ConvertSample(SampleId::SampleEnum(SelectionNom)) << " detbin: " << BinRedu << " fitbin: " << BinFull << std::endl;
        if(BinFull > nBinsFull){
          std::cout << "BinFull > nBinsFull" << std::endl;
          throw;
        }
        if (BinFull == -1) {
          std::cout << "Couldn't find a bin in GetGlobalBinmatrixMomCos for " << SelectionNom << " (" << SampleId::ConvertSample(SampleId::SampleEnum(SelectionNom)) << ")" << std::endl;
          throw;
        }
        if(TMath::IsNaN(WeightNom    ) || !TMath::Finite(WeightNom    ) || WeightNom      == -999 || WeightNom      > 1000){ std::cerr << "WeightNom isn't finite: "     << WeightNom     << std::endl; throw;}
        if(TMath::IsNaN(NIWGWeightNom) || !TMath::Finite(NIWGWeightNom) || NIWGWeightNom  == -999 || NIWGWeightNom  > 1000){ std::cerr << "NIWGWeightNom isn't finite: " << NIWGWeightNom << std::endl; throw;}
        if(TMath::IsNaN(FluxWeightNom) || !TMath::Finite(FluxWeightNom) || FluxWeightNom  == -999 || FluxWeightNom  > 1000){ std::cerr << "FluxWeightNom isn't finite: " << FluxWeightNom << std::endl; throw;}
        if(TMath::IsNaN(POT_Weight   ) || !TMath::Finite(POT_Weight   ) || POT_Weight     == -999 || POT_Weight     > 1000){ std::cerr << "POT_Weight isn't finite: "    << POT_Weight    << std::endl; throw;}

        number_events_redu->Fill(BinRedu);
        number_events_full->Fill(BinFull);

        nominal           ->SetBinContent(BinRedu, nominal           ->GetBinContent(BinRedu) + WeightNom * NIWGWeightNom * FluxWeightNom * POT_Weight);
        mc_stat_nom       ->SetBinContent(BinFull, mc_stat_nom       ->GetBinContent(BinFull) + WeightNom * FluxWeightNom * NIWGWeightNom * POT_Weight);
        mc_stat_error_full->SetBinContent(BinFull, mc_stat_error_full->GetBinContent(BinFull) + TMath::Power(WeightNom * FluxWeightNom * NIWGWeightNom * POT_Weight, 2.));
      }

      for(int iToy = 0; iToy < nThrows; ++iToy){
        if(bd.IsActiveSample(SelectionToy[iToy])){
          int BinRedu = bd.GetGlobalBinMatrixMomCos_Det((SampleId::SampleEnum)SelectionToy[iToy], LeptonMomToy[iToy], LeptonCosToy[iToy])+1;

          if(BinRedu > nBinsRedu){
            std::cout << "BinRedu > nBinsRedu" << std::endl;
            throw;
          }
          if (BinRedu == -1) {
            std::cout << "Couldn't find a bin in GetGlobalBinmatrixMomCos_Det for " << SelectionNom << " (" << SampleId::ConvertSample(SampleId::SampleEnum(SelectionNom)) << ")" << std::endl;
            throw;
          }
          if(!TMath::IsNaN(WeightToy[iToy]) && TMath::Finite(WeightToy[iToy])){

            if(WeightToy[iToy] < 100){
              varied[iToy]->SetBinContent(BinRedu, varied[iToy]->GetBinContent(BinRedu) + POT_Weight * FluxWeightToy[iToy] * NIWGWeightToy[iToy] * WeightToy[iToy]);
              if(TMath::IsNaN(WeightToy[iToy]    ) || !TMath::Finite(WeightToy[iToy]    ) || WeightToy[iToy]      == -999 || WeightNom      > 1000){ std::cerr << "WeightToy["     << iToy << "] isn't finite: " << WeightToy[iToy]     << std::endl; throw;}
              if(TMath::IsNaN(NIWGWeightToy[iToy]) || !TMath::Finite(NIWGWeightToy[iToy]) || NIWGWeightToy[iToy]  == -999 || NIWGWeightNom  > 1000){ std::cerr << "NIWGWeightToy[" << iToy << "] isn't finite: " << NIWGWeightToy[iToy] << std::endl; throw;}
              if(TMath::IsNaN(FluxWeightToy[iToy]) || !TMath::Finite(FluxWeightToy[iToy]) || FluxWeightToy[iToy]  == -999 || FluxWeightNom  > 1000){ std::cerr << "FluxWeightToy[" << iToy << "] isn't finite: " << FluxWeightToy[iToy] << std::endl; throw;}
              if(TMath::IsNaN(POT_Weight         ) || !TMath::Finite(POT_Weight         ) || POT_Weight           == -999 || POT_Weight     > 1000){ std::cerr << "POT_Weight isn't finite: "                    << POT_Weight          << std::endl; throw;}
            }else{
              varied[iToy]->SetBinContent(BinRedu, varied[iToy]->GetBinContent(BinRedu) + FluxWeightToy[iToy] * NIWGWeightToy[iToy] * POT_Weight);
            }
          }
        }
      }
    } // End the total for loop over events
    // Now summaries

    for(int iToy = 0; iToy < nThrows; ++iToy){
      for(int iBin = 1; iBin < nBinsRedu+1; ++iBin){
        var_mean->SetBinContent(iBin, var_mean->GetBinContent(iBin) + varied[iToy]->GetBinContent(iBin)/double(nThrows));
      }
    }

    for(int iToy = 0; iToy < nThrows; ++iToy){
      for(int iBin = 1; iBin < nBinsRedu+1; ++iBin){
        for(int jBin = 1; jBin < nBinsRedu+1; ++jBin){
          cov_var->SetBinContent(iBin, jBin, cov_var->GetBinContent(iBin,jBin) + ((varied[iToy]->GetBinContent(iBin) - var_mean->GetBinContent(iBin))*(varied[iToy]->GetBinContent(jBin) - var_mean->GetBinContent(jBin))));
        }
      }
    }

    for(int iBin = 1; iBin < nBinsFull+1; ++iBin){
      if(mc_stat_nom->GetBinContent(iBin)>0){
        mc_stat_error_full->SetBinContent(iBin, mc_stat_error_full->GetBinContent(iBin) / (mc_stat_nom->GetBinContent(iBin) * mc_stat_nom->GetBinContent(iBin)));
      }
    }

    int lastk = 1;
    for(int j = 0; j < nBinsFull; ++j){
      double staterr= mc_stat_error_full->GetBinContent(j+1);
      double fakeerr= fake_error_full->GetBinContent(j+1);
      SampleId::SampleEnum sample1 = bd.GetSampleFromIndex(j);
      int localbin = j - bd.GetGlobalBinMatrix(sample1);
      double p    = bd.GetBinning_Mom(sample1)->GetBinUpEdge(int(localbin%bd.GetNbins_Mom(sample1))+1);
      double c    = bd.GetBinning_Cos(sample1)->GetBinUpEdge(int(localbin/bd.GetNbins_Mom(sample1))+1);
      double plow = bd.GetBinning_Mom(sample1)->GetBinLowEdge(int(localbin%bd.GetNbins_Mom(sample1))+1);
      double clow = bd.GetBinning_Cos(sample1)->GetBinLowEdge(int(localbin/bd.GetNbins_Mom(sample1))+1);

      std::cout << "p = " << p << " plow = " << plow << " c = " << c << " clow = " << clow << " sample = " << SampleId::ConvertSample(sample1) << std::endl;
      for(int k = bd.GetGlobalBinMatrix_Det(sample1); k < bd.GetGlobalBinMatrix_Det(sample1) + bd.GetNbins_Det(sample1); ++k){

        localbin = k - bd.GetGlobalBinMatrix_Det(sample1);

        // If the fit bin does not overlap at all with the systematics obsnorm bin then skip it.
        if (bd.GetBinning_Mom_Det(sample1)->GetBinUpEdge(int(localbin%bd.GetNbins_Mom_Det(sample1))+1) <= plow || 
            bd.GetBinning_Cos_Det(sample1)->GetBinUpEdge(int(localbin/bd.GetNbins_Mom_Det(sample1))+1) <= clow || 
            bd.GetBinning_Mom_Det(sample1)->GetBinLowEdge(int(localbin%bd.GetNbins_Mom_Det(sample1))+1) >= p ||
            bd.GetBinning_Cos_Det(sample1)->GetBinLowEdge(int(localbin/bd.GetNbins_Mom_Det(sample1))+1) >= c)
          continue;

        std::cout << " plow+det = " << bd.GetBinning_Mom_Det(sample1)->GetBinLowEdge(int(localbin%bd.GetNbins_Mom_Det(sample1))+1) << std::endl;
        std::cout << " p+det = " << bd.GetBinning_Mom_Det(sample1)->GetBinUpEdge(int(localbin%bd.GetNbins_Mom_Det(sample1))+1) << std::endl;
        std::cout << " clow+det = " << bd.GetBinning_Cos_Det(sample1)->GetBinLowEdge(int(localbin/bd.GetNbins_Mom_Det(sample1))+1) << std::endl;
        std::cout << " c+det = " << bd.GetBinning_Cos_Det(sample1)->GetBinUpEdge(int(localbin/bd.GetNbins_Mom_Det(sample1))+1) << std::endl;
        std::cout << " local bin = " << localbin << std::endl;
        std::cout << " Det Bin = " << k << std::endl;

        if(staterr > mc_stat_error_redu->GetBinContent(k + 1)){
          mc_stat_error_redu->SetBinContent(k + 1, staterr);
        } 
        if(Apply1p1h){
          if(fabs(fakeerr) > fabs(fake_error_redu->GetBinContent(k + 1))){
            fake_error_redu->SetBinContent(k + 1, fakeerr);
          }
        }
        break;
      }
    }


    if(Apply1p1h){
      fake_error_redu->Write("Fake_Error_Lin");

      TH1D* fake_lin = new TH1D("Lin_1p1h","Lin_1p1h",nBinsRedu,0,nBinsRedu); 

      for(int j = 1; j <= nBinsRedu; ++j){
        for(int k = j+1; k <= nBinsRedu; ++k){
          fake_error_cov->SetBinContent(j,k, fake_error_redu->GetBinContent(j) * fake_error_redu->GetBinContent(k));
        }
        fake_error_cov->SetBinContent(j,j, fake_error_redu->GetBinContent(j) * fake_error_redu->GetBinContent(j));
        fake_lin->SetBinContent(j,fake_error_cov->GetBinContent(j,j));
      }

      fake_lin->Write();
      fake_error_full->Write();
      fake_cov->Write();

      for(int j = 1; j <= nBinsRedu; ++j){
        for(int k = 1; k < j; ++k){
          fake_error_cov->SetBinContent(j,k, fake_error_cov->GetBinContent(k,j));
        }
      }
    }

    for(int iBin = 1; iBin <= nBinsRedu; ++iBin){
      if(nominal->GetBinContent(iBin)>0){
        vector(iBin-1) = var_mean->GetBinContent(iBin) / nominal->GetBinContent(iBin);
        var_mean_offset->SetBinContent(iBin, var_mean->GetBinContent(iBin)/nominal->GetBinContent(iBin));
      }else{
        vector(iBin-1) = 1.0;
        var_mean_offset->SetBinContent(iBin,1);
      }

      for(int jBin = 1; jBin <= nBinsRedu; ++jBin){
        if(fabs(cov_var->GetBinContent(iBin,jBin)) > 0 && var_mean->GetBinContent(iBin) > 0 && var_mean->GetBinContent(jBin) > 0){
          cov_var->SetBinContent(iBin, jBin, cov_var->GetBinContent(iBin,jBin)/(double(nThrows)*var_mean->GetBinContent(iBin)*var_mean->GetBinContent(jBin)));
        }

        cov->SetBinContent(iBin, jBin, cov_var->GetBinContent(iBin,jBin));

        if(jBin == iBin){
          mc_sys_error->SetBinContent(iBin, cov_var->GetBinContent(iBin,jBin));
          cov->SetBinContent(iBin, jBin, cov->GetBinContent(iBin,jBin) + mc_stat_error_redu->GetBinContent(iBin));//WP: Taking out so MC stats not added to cov
        }

        if(Apply1p1h) matrix(iBin-1,jBin-1) = cov->GetBinContent(iBin,jBin) + fake_error_cov->GetBinContent(iBin,jBin);
        else          matrix(iBin-1,jBin-1) = cov->GetBinContent(iBin,jBin);
        cov_nofake->SetBinContent(iBin, jBin, cov->GetBinContent(iBin,jBin));
      }
      mc_stat_error_redu->SetBinContent(iBin, sqrt(mc_stat_error_redu->GetBinContent(iBin)));
    }

    for(int iBin = 1; iBin <= nBinsRedu; ++iBin){
      for(int jBin = 1; jBin <= nBinsRedu; ++jBin){
        if(cov->GetBinContent(iBin,iBin) > 0 && cov->GetBinContent(jBin,jBin) > 0)
          corr_nofake->SetBinContent(iBin, jBin, cov->GetBinContent(iBin,jBin) / sqrt(cov->GetBinContent(iBin,iBin) * cov->GetBinContent(jBin,jBin)));
        if(cov_var->GetBinContent(iBin,iBin) > 0 && cov_var->GetBinContent(jBin,jBin) > 0)
          corr_var->SetBinContent(iBin, jBin, cov_var->GetBinContent(iBin,jBin) / sqrt(cov_var->GetBinContent(iBin,iBin) * cov_var->GetBinContent(jBin,jBin)));
      }
    }

    if(Apply1p1h) cov->Add(&*fake_error_cov);

    for(int iBin = 1; iBin <= nBinsRedu; ++iBin){
      for(int jBin = 1; jBin <= nBinsRedu; ++jBin){
        if(cov->GetBinContent(iBin,iBin) > 0 && cov->GetBinContent(jBin,jBin) > 0){
          corr->SetBinContent(iBin, jBin, cov->GetBinContent(iBin,jBin) / sqrt(cov->GetBinContent(iBin,iBin) * cov->GetBinContent(jBin,jBin)));
          cov_zero->SetBinContent(iBin, jBin, fake_error_cov->GetBinContent(iBin,jBin) / sqrt(fake_error_cov->GetBinContent(iBin,iBin) * fake_error_cov->GetBinContent(jBin,jBin)));
        }
      }
    }

    OutputFile->cd();

    vector.Write("Mean_Value");
    matrix.Write("Covariance_Matrix");
    std::string CurSampleStr = "";
    for(int iBin = 1; iBin < nBinsRedu+1; ++iBin){
      std::string SampleStr = ConvertSample(bd.GetSampleFromIndex_Det(iBin));
      if(SampleStr != CurSampleStr){
        cov        ->GetXaxis()->SetBinLabel(iBin,SampleStr.c_str());
        corr       ->GetXaxis()->SetBinLabel(iBin,SampleStr.c_str());
        corr_nofake->GetXaxis()->SetBinLabel(iBin,SampleStr.c_str());
        cov_zero   ->GetXaxis()->SetBinLabel(iBin,SampleStr.c_str());
        cov_nofake ->GetXaxis()->SetBinLabel(iBin,SampleStr.c_str());
        cov        ->GetYaxis()->SetBinLabel(iBin,SampleStr.c_str());
        corr       ->GetYaxis()->SetBinLabel(iBin,SampleStr.c_str());
        corr_nofake->GetYaxis()->SetBinLabel(iBin,SampleStr.c_str());
        cov_zero   ->GetYaxis()->SetBinLabel(iBin,SampleStr.c_str());
        cov_nofake ->GetYaxis()->SetBinLabel(iBin,SampleStr.c_str());
        fake_error_redu    ->GetXaxis()->SetBinLabel(iBin,SampleStr.c_str());
        nominal            ->GetXaxis()->SetBinLabel(iBin,SampleStr.c_str());
        mc_sys_error       ->GetXaxis()->SetBinLabel(iBin,SampleStr.c_str());
        mc_stat_error_redu ->GetXaxis()->SetBinLabel(iBin,SampleStr.c_str());

        var_mean           ->GetXaxis()->SetBinLabel(iBin,SampleStr.c_str());
        var_mean_offset    ->GetXaxis()->SetBinLabel(iBin,SampleStr.c_str());
        number_events_redu ->GetXaxis()->SetBinLabel(iBin,SampleStr.c_str());
        GOF                ->GetXaxis()->SetBinLabel(iBin,SampleStr.c_str());
        GOFnDOF            ->GetXaxis()->SetBinLabel(iBin,SampleStr.c_str());
        DiffMean           ->GetXaxis()->SetBinLabel(iBin,SampleStr.c_str());
        DiffError          ->GetXaxis()->SetBinLabel(iBin,SampleStr.c_str());
        for (int j = 0; j < nThrows; ++j) {
          varied[j]->GetYaxis()->SetTitle("N_{events} MC");
          varied[j]->GetXaxis()->SetBinLabel(iBin, SampleStr.c_str());
        }
        CurSampleStr = SampleStr;
      }
    }

    for (int iBin = 1; iBin < nBinsFull+1; ++iBin) {
      std::string SampleStr = ConvertSample(bd.GetSampleFromIndex(iBin));
      if(SampleStr != CurSampleStr){
        mc_stat_error_full    ->GetXaxis()->SetBinLabel(iBin, SampleStr.c_str());
        mc_stat_nom           ->GetXaxis()->SetBinLabel(iBin, SampleStr.c_str());
        number_events_full    ->GetXaxis()->SetBinLabel(iBin, SampleStr.c_str());
        fake_error_full       ->GetXaxis()->SetBinLabel(iBin, SampleStr.c_str());
        CurSampleStr = SampleStr;
      }
    }

    for(int iBin = 1; iBin < Selection->GetXaxis()->GetNbins()+1; ++iBin){
      std::string SampleStr = SampleId::ConvertSample((SampleId::SampleEnum)(iBin-1));
      Selection->GetXaxis()->SetBinLabel(iBin, SampleStr.c_str());
    }

    Selection          ->Write();
    corr               ->Write();
    corr_nofake        ->Write();
    cov                ->Write();
    cov_var            ->Write();
    cov_nofake         ->Write();
    nominal            ->Write();
    mc_sys_error       ->Write();
    mc_stat_error_redu ->Write(); 
    var_mean            ->Write();
    var_mean_offset    ->Write();
    number_events_redu ->Write();
    // Take root of MC stats plot to convert it from variance into an error
    for(int i = 1; i <= mc_stat_error_full->GetNbinsX(); ++i){
      mc_stat_error_full->SetBinContent(i, sqrt(mc_stat_error_full->GetBinContent(i)));
    }
    mc_stat_error_full ->Write();
    mc_stat_nom        ->Write();
    number_events_full ->Write();
    Selection          ->Write();
    fake_error_cov     ->Write();
    cov_zero           ->Write();
    corr_var->Write();

    const int NCont = 255;
    int MyPalette[NCont];
    const Int_t NRGBs = 5;
    Double_t stops[NRGBs] = { 0.00, 0.25, 0.50, 0.75, 1.00 };
    Double_t red[NRGBs]   = { 0.00, 0.25, 1.00, 1.00, 0.50 };
    Double_t green[NRGBs] = { 0.00, 0.25, 1.00, 0.25, 0.00 };
    Double_t blue[NRGBs]  = { 0.50, 1.00, 1.00, 0.25, 0.00 };
    int start = TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    for (int i = 0; i < NCont; ++i) MyPalette[i] = start+i;
    gStyle->SetPalette(NCont, MyPalette);  // use the nice blue->white->red palette

    TCanvas c1;
    cov->Draw("COLZ");
    c1.SetTitle("Covariance");
    c1.SetName ("canCovariance");
    c1.Write();
    cov_zero->Draw("COLZ");
    c1.SetTitle("CovarianceZero");
    c1.SetName ("canCovarianceZero");
    c1.Write();
    cov_nofake->Draw("COLZ");
    c1.SetTitle("CovarianceNoFake");
    c1.SetName ("canCovarianceNoFake");
    c1.Write();


    corr->SetMaximum( 1);
    corr->SetMinimum(-1);
    corr->Draw("COLZ");
    c1.SetTitle("Correlation");
    c1.SetName ("canCorrelation");
    c1.Write();

    TDirectory* VarU = OutputFile->mkdir("VariedUniverse");
    VarU->cd();
    for (int i = 0; i < nThrows; ++i) {
      varied[i]->Write();
    }

    TCanvas c2;
    OutputFile->cd();
    TDirectory* BinVar = OutputFile->mkdir("BinVariations");
    BinVar->cd();

    for(int iBin = 1; iBin < nBinsRedu+1; ++iBin){

      c2.SetName(Form("Bin_%d_NEvent_Spread",iBin));
      c2.SetTitle(Form("Bin_%d_NEvent_Spread",iBin));
      int nbins = nominal->GetBinContent(iBin)*80;
      if (nbins > 1000) nbins = 1000;
      double lower = 0.0;
      double upper = nominal->GetBinContent(iBin)*3.0;
      std::string title = Form("Bin %d", iBin);
      // Get the sample
      SampleId::SampleEnum sample = SampleId::kUnassigned;
      int mombin = -1;
      int cosbin = -1;
      bd.GetLocalBinMatrixMomCos_Det(iBin, sample, mombin, cosbin);
      if (sample != SampleId::kUnassigned) {
        TAxis* xaxis = bd.GetBinning_Mom_Det(sample);
        TAxis* yaxis = bd.GetBinning_Cos_Det(sample);
        double momlow = xaxis->GetBinLowEdge(mombin+1);
        double momhi = xaxis->GetBinLowEdge(mombin+2);
        double coslow = yaxis->GetBinLowEdge(cosbin+1);
        double coshi = yaxis->GetBinLowEdge(cosbin+2);
        std::stringstream ss;
        ss << momlow << "-" << momhi;
        std::stringstream ss2;
        ss2 << coslow << "-" << coshi;
        title += " (" + SampleId::ConvertSample(sample) + ", {p_{#mu}, cos#theta_{#mu}} = {" + ss.str() + ", " + ss2.str() + "})";
      }

      TH1F h1(Form("Bin_%d_NEvent_Spread_Hist", iBin), title.c_str(), nbins, lower, upper);

      TF1 f1("f1","gaus",nominal->GetBinContent(iBin)*0.5, nominal->GetBinContent(iBin)*2.0);
      TF1 f2("f2","gaus",nominal->GetBinContent(iBin)*0.5, nominal->GetBinContent(iBin)*2.0);

      f1.SetNameTitle(Form("Bin_%d_NEvent_Varying_Gaussian",iBin), "Varying Gaussian");
      f2.SetNameTitle(Form("Bin_%d_NEvent_Total_Gaussian",iBin), "Total Gaussian");

      for(int iToy = 0; iToy < nThrows; ++iToy){
        h1.Fill(varied[iToy]->GetBinContent(iBin));
      }

      h1.GetXaxis()->SetTitle("N_{Events}^{MC} selected");
      h1.GetYaxis()->SetTitle("N_{toys}");
      h1.SetStats(kFALSE);
      h1.Draw();

      f1.SetParameters(h1.Integral()/(sqrt(2.0*TMath::Pi()*cov_var->GetBinContent(iBin,iBin))*nominal->GetBinContent(iBin)),
          var_mean->GetBinContent(iBin),
          sqrt(cov_var->GetBinContent(iBin,iBin))*nominal->GetBinContent(iBin));

      f2.SetParameters(h1.Integral()/(sqrt(2.0*TMath::Pi()*cov->GetBinContent(iBin,iBin))*nominal->GetBinContent(iBin)),
          var_mean->GetBinContent(iBin),
          sqrt(cov->GetBinContent(iBin,iBin))*nominal->GetBinContent(iBin));

      h1.GetYaxis()->SetRangeUser(0,h1.GetMaximum()*1.1);
      h1.GetXaxis()->SetRange(h1.FindFirstBinAbove(0)-5,h1.FindLastBinAbove(0)+5);

      //for (int i = 0; i < 3; ++i) {
      //std::cout << i << ": " << f1.GetParameter(i) << std::endl;
      //}
      double CalcMean = f1.GetParameter(1);
      double CalcError = f1.GetParameter(2);

      //f1.FixParameter(0,f1.GetParameter(0));
      //f1.FixParameter(1,f1.GetParameter(1));
      //f1.FixParameter(2,f1.GetParameter(2));
      f2.FixParameter(1,f2.GetParameter(1));
      f2.FixParameter(2,f2.GetParameter(2));

      h1.Fit(&f1,"BNQ");
      f2.FixParameter(0, f1.GetParameter(0));
      //for (int i = 0; i < 3; ++i) {
      //std::cout << i << ": " << f1.GetParameter(i) << std::endl;
      //}
      double FitMean = f1.GetParameter(1);
      double FitError = f1.GetParameter(2);

      double chi2 = f1.GetChisquare();
      int lowbin = h1.FindBin(nominal->GetBinContent(iBin)*0.5);
      int hibin = h1.FindBin(nominal->GetBinContent(iBin)*1.5);

      GOF->SetBinContent(iBin,chi2);
      GOFnDOF->SetBinContent(iBin, chi2/h1.GetNbinsX());
      DiffMean->SetBinContent(iBin,(FitMean-CalcMean)/nominal->GetBinContent(iBin));
      DiffError->SetBinContent(iBin,(FitError-CalcError)/nominal->GetBinContent(iBin));

      TLine l1(nominal->GetBinContent(iBin),0,nominal->GetBinContent(iBin),h1.GetMaximum());
      //TLine l2(zero_point->GetBinContent(iBin),0,zero_point->GetBinContent(iBin),h1.GetMaximum());

      f1.SetLineWidth(2);
      f1.SetLineColor(kRed);
      f2.SetLineWidth(2);
      f2.SetLineColor(kGreen);

      f1.Draw("SAME");
      f2.Draw("SAME");
      h1.SetLineWidth(2);

      l1.SetLineWidth(2);
      l1.SetLineStyle(2);
      l1.SetLineColor(kBlack);
      //l2.SetLineWidth(2);
      //l2.SetLineStyle(2);
      //l2.SetLineColor(kOrange);

      l1.Draw("SAME");
      //l2.Draw("SAME");

      TLegend* leg = new TLegend(0.65, 0.65, 0.98, 0.94);
      double mean = h1.GetMean();
      double rms = h1.GetRMS();
      leg->AddEntry(&h1, Form("Bin content from throws, #mu=%.1f#pm%.1f", mean, rms), "l");
      leg->AddEntry(&f1, Form("Gaussian from thrown cov, #chi^{2}_{nbins}=%.1f, #mu=%.1f#pm%.1f", chi2/h1.GetNbinsX(), FitMean, FitError),"l");
      leg->AddEntry(&f2, Form("Gaussian from total cov, #mu=%.1f#pm%.1f", f2.GetParameter(1), f2.GetParameter(2)), "l");
      leg->AddEntry(&l1, Form("Nominal bin content, %.1f events", nominal->GetBinContent(iBin)), "l");
      //leg->AddEntry(&l2, "Zero weight bin content","l");
      leg->SetFillColor(kWhite);
      leg->Draw();
      c2.Update();
      c2.Write();
    }

    OutputFile->cd();
    GOF->Write();
    GOFnDOF->Write();
    DiffMean->Write();
    DiffError->Write();
    TDirectory* AxisCovDir = OutputFile->mkdir("AxisForCov");
    OutputFile->cd();
    TDirectory* AxisFitDir = OutputFile->mkdir("AxisForFit");

    for(int i = SampleId::kUnassigned; i < SampleId::kNSamples; i++){
      SampleId::SampleEnum smpl = static_cast<SampleId::SampleEnum>(i);

      if(bd.IsActiveSample(smpl)){
        std::string name = SampleId::ConvertSample(smpl);
        std::replace (name.begin(), name.end(), ' ', '_');
        AxisCovDir->cd();
        bd.GetBinning_Mom_Det(smpl)->Write((name+"_p_axis").c_str());
        bd.GetBinning_Cos_Det(smpl)->Write((name+"_th_axis").c_str());
        AxisFitDir->cd();
        bd.GetBinning_Mom(smpl)->Write((name+"_p_axis").c_str());
        bd.GetBinning_Cos(smpl)->Write((name+"_th_axis").c_str());
      }

    }
    OutputFile->Close();
  }
}


void MakeCovariance4Arg(char* MCFilesListStr, char* NIWGFilesListStr,
    char* Cov1P1HStr,     char* OutFileStr){
  MakeCovariance(MCFilesListStr, NIWGFilesListStr, Cov1P1HStr, OutFileStr);
}

void MakeCovariance3Arg(char* MCFilesListStr, char* NIWGFilesListStr, char* OutFileStr){
  MakeCovariance(MCFilesListStr, NIWGFilesListStr, "", OutFileStr);
}

void MakeCovariance2Arg(char* MCFilesListStr, char* OutFileStr){
  MakeCovariance(MCFilesListStr, "", "", OutFileStr);
}

int main(int argn, char** arg){

  switch (argn){
    case 5:{
             ApplyNIWG = true;
             Apply1p1h = true;
             MakeCovariance4Arg(arg[1], arg[2], arg[3], arg[4]);
             break;
           }
    case 4:{
             ApplyNIWG = true;
             Apply1p1h = false;
             MakeCovariance3Arg(arg[1], arg[2], arg[3]);
             break;
           }
    case 3:{
             ApplyNIWG = false;
             Apply1p1h = false;
             MakeCovariance2Arg(arg[1], arg[2]);
             break;
           }
    default:{
              std::cout << "\n\n\nUsage:\nMakeND280DetectorCov.exe ListOfFilesWithNDThrows ListOfFileWithNIWGWeight 1p1hCov.root outfile.root" << std::endl;
              std::cout << "or:\nMakeND280DetectorCov.exe ListOfFilesWithNDThrows ListOfFileWithNIWGWeight outfile.root" << std::endl;
              std::cout << "or:\nMakeND280DetectorCov.exe ListOfFilesWithNDThrows outfile.root" << std::endl;
              throw;
            }
  }

  std::cout << "Finished" << std::endl;

  return 0;

}


