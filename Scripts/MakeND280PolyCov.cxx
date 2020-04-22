#include <iostream>
#include <fstream>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TH2Poly.h"
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

#include "PolyBinning.hxx"
#include "ND280AnalysisUtils.hxx"
#include "Parameters.hxx"

bool ApplyNIWG = false;
bool Apply1p1h = false;
const int nThrows = 2000;

void MakeCovariance(char* MCFilesListStr, char* NIWGFilesListStr,
    char* Cov1P1HStr,     char* OutFileStr){

  ND::params().SetReadParamOverrideFilePointPassed();
  BANFF::PolyBinning& pb = BANFF::PolyBinning::Get();
  pb.SetupPoly();
  //WP will need to get binning details

  // Count the number of active samples
  int nCountedSamples = 0;
  for (int i = 0; i < SampleId::kNSamples; ++i) {
    if (pb.IsActiveSample(SampleId::SampleEnum(i))) nCountedSamples++;
  }//WP can hack something to count samples we're using

  std::cout << "Found " << nCountedSamples << " enabled samples in PolyBinning" << std::endl;

  TFile* fake = NULL;
  TH1D* fake_cov = NULL;

  int nSamples = 18;
  //WP hard code 18
  
  const int nSamplesMultiTrack = 14;
  const int nSamplesMultiPi = 18;

  if (nSamplesMultiTrack != nCountedSamples && nSamplesMultiPi != nCountedSamples) {
    std::cerr << "Hard-coded number of samples does not match found samples in PolyBinning" << std::endl;
    std::cerr << "nSamples (hard-coded): " << nSamples << std::endl;
    std::cerr << "nFoundSamples: " << nCountedSamples << std::endl;
    throw;
  }
  nSamples = nCountedSamples;

  //WP probably want these to be polys
  TH2Poly** Nominal_Plots = new TH2Poly*[nSamples];
  TH2Poly** Fake_Plots = new TH2Poly*[nSamples];

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

    }else{
      std::cout << "WARNING!! NOT APPLYING THE NIWG WEIGHTING" << std::endl;
      for(int i = 0;i < nThrows; i++) NIWGWeightToy[i] = 1;
      NIWGWeightNom = 1.;
    }

    std::cout << "Creating the file " << OutFileStr << ".root" << std::endl;
    TFile* OutputFile = new TFile(Form("%s_%s.root",OutFileStr,(it->first).c_str()), "RECREATE");

    //WP just get total number of bins for each
    int nBinsRedu = pb.GetNbins_Det();
    int nBinsFull = pb.GetNbins();

    std::cout << "We have " << nBinsRedu <<  " bins in the covariance matrix" << std::endl;

    TMatrixTSym<double> matrix(nBinsRedu);
    TVectorT   <double> vector(nBinsRedu);
    //WP make polys, same names, just a clone of inputted polys prolly?
    //WP Actually no most of these stay as th2ds :)
    TH2D* corr               = new TH2D("Total_Correlation_Matrix", "Total Correlation Matrix",                      nBinsRedu, 0, nBinsRedu, nBinsRedu, 0, nBinsRedu);
    TH2D* corr_nofake        = new TH2D("Total_Correlation_Matrix_No1p1h", "Total Correlation Matrix without 1p1h error",                      nBinsRedu, 0, nBinsRedu, nBinsRedu, 0, nBinsRedu);
    TH2D* cov                = new TH2D("Total_Covariance_Matrix",  "Total Covariance Matrix",                       nBinsRedu, 0, nBinsRedu, nBinsRedu, 0, nBinsRedu);
    TH2D* cov_var            = new TH2D("Covariance_Matrix_NoMCStats",    "Covariance Matrix without MC stats error",               nBinsRedu, 0, nBinsRedu, nBinsRedu, 0, nBinsRedu);
    TH2D* corr_var           = new TH2D("Correlation_Matrix_NoMCStats",   "Correlation Matrix without MC stats error",              nBinsRedu, 0, nBinsRedu, nBinsRedu, 0, nBinsRedu);
    TH2D* cov_zero           = new TH2D("Correlation_Matrix_1p1h",  "Correlation Matrix for 1p1h error",             nBinsRedu, 0, nBinsRedu, nBinsRedu, 0, nBinsRedu);
    TH2D* cov_nofake         = new TH2D("Covariance_Matrix_NoFake", "Covariance Matrix without 1p1h error",          nBinsRedu, 0, nBinsRedu, nBinsRedu, 0, nBinsRedu);
    TH2D* fake_error_cov     = new TH2D("1p1h_fake_Error_cov",      "1p1h systematic error",                         nBinsRedu, 0, nBinsRedu, nBinsRedu, 0, nBinsRedu);
    //WP 1Ds are ok
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

      //WP check sample again
      if(pb.IsActiveSample(SelectionNom)){
        Selection->Fill(SelectionNom);
	//WP Have function to get global
        int BinRedu = pb.GetGlobalBinMatrixMomCos_Det((SampleId::SampleEnum)SelectionNom, LeptonMomNom, LeptonCosNom);
        if(BinRedu > nBinsRedu){
          std::cout << "BinRedu > nBinsRedu" << std::endl;
          throw;
        }

	//WP again get global
        int BinFull = pb.GetGlobalBinMatrixMomCos((SampleId::SampleEnum)SelectionNom,
            LeptonMomNom, LeptonCosNom);

        //std::cout << SampleId::ConvertSample(SampleId::SampleEnum(SelectionNom)) << " detbin: " << BinRedu << " fitbin: " << BinFull << std::endl;
        if(BinFull > nBinsFull){
          std::cout << "BinFull > nBinsFull" << std::endl;
          throw;
        }

        if(TMath::IsNaN(WeightNom    ) || !TMath::Finite(WeightNom    ) || WeightNom      == -999 || WeightNom      > 1000){ std::cerr << "WeightNom isn't finite: "     << WeightNom     << std::endl; throw;}
        if(TMath::IsNaN(NIWGWeightNom) || !TMath::Finite(NIWGWeightNom) || NIWGWeightNom  == -999 || NIWGWeightNom  > 1000){ std::cerr << "NIWGWeightNom isn't finite: " << NIWGWeightNom << std::endl; throw;}
        if(TMath::IsNaN(FluxWeightNom) || !TMath::Finite(FluxWeightNom) || FluxWeightNom  == -999 || FluxWeightNom  > 1000){ std::cerr << "FluxWeightNom isn't finite: " << FluxWeightNom << std::endl; throw;}
        if(TMath::IsNaN(POT_Weight   ) || !TMath::Finite(POT_Weight   ) || POT_Weight     == -999 || POT_Weight     > 1000){ std::cerr << "POT_Weight isn't finite: "    << POT_Weight    << std::endl; throw;}

	//WP this ok as 1ds 
        number_events_redu->Fill(BinRedu);
        number_events_full->Fill(BinFull);
	
	//WP I think these still ok
        nominal           ->SetBinContent(BinRedu, nominal           ->GetBinContent(BinRedu) + WeightNom * NIWGWeightNom * FluxWeightNom * POT_Weight);
        mc_stat_nom       ->SetBinContent(BinFull, mc_stat_nom       ->GetBinContent(BinFull) + WeightNom * FluxWeightNom * NIWGWeightNom * POT_Weight);
        mc_stat_error_full->SetBinContent(BinFull, mc_stat_error_full->GetBinContent(BinFull) + TMath::Power(WeightNom * FluxWeightNom * NIWGWeightNom * POT_Weight, 2.));
      }

      for(int iToy = 0; iToy < nThrows; ++iToy){
        if(pb.IsActiveSample(SelectionToy[iToy])){
          //WP Get global again
	  int BinRedu = pb.GetGlobalBinMatrixMomCos_Det((SampleId::SampleEnum)SelectionToy[iToy], LeptonMomToy[iToy], LeptonCosToy[iToy]) ;

          if(BinRedu > nBinsRedu){
            std::cout << "BinRedu > nBinsRedu" << std::endl;
            //throw;
          }

          if(!TMath::IsNaN(WeightToy[iToy]) && TMath::Finite(WeightToy[iToy])){

            if(WeightToy[iToy] < 100){
	      //WP varied 1d so ok
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

    //WP 1d ok
    for(int iToy = 0; iToy < nThrows; ++iToy){
      for(int iBin = 1; iBin < nBinsRedu+1; ++iBin){
        var_mean->SetBinContent(iBin, var_mean->GetBinContent(iBin) + varied[iToy]->GetBinContent(iBin)/double(nThrows));
      }
    }

    //WP will be poly but setbincontent ok
    for(int iToy = 0; iToy < nThrows; ++iToy){
      for(int iBin = 1; iBin < nBinsRedu+1; ++iBin){
        for(int jBin = 1; jBin < nBinsRedu+1; ++jBin){
          cov_var->SetBinContent(iBin, jBin, cov_var->GetBinContent(iBin,jBin) + ((varied[iToy]->GetBinContent(iBin) - var_mean->GetBinContent(iBin))*(varied[iToy]->GetBinContent(jBin) - var_mean->GetBinContent(jBin))));
        }
      }
    }

    //WP fine, we wont use but can keep in here
    for(int iBin = 1; iBin < nBinsFull+1; ++iBin){
      if(mc_stat_nom->GetBinContent(iBin)>0){
        mc_stat_error_full->SetBinContent(iBin, mc_stat_error_full->GetBinContent(iBin) / (mc_stat_nom->GetBinContent(iBin) * mc_stat_nom->GetBinContent(iBin)));
      }
    }
    
    int lastk = 1;
    for(int j = 0; j < nBinsFull; ++j){
      //WP ok
      double staterr= mc_stat_error_full->GetBinContent(j+1);
      double fakeerr= fake_error_full->GetBinContent(j+1);
      //WP need a sample to go from global bin num to sample
      SampleId::SampleEnum sample1 = pb.GetSampleFromIndex(j);
      //WP subtract sample offset (need to make function or do on fly)
      int localbin = j - pb.GetGlobalBinMatrixOffset(sample1);
      //WP get p and c edges of Fit bin, I think can do
      TH2PolyBin* pbin = (TH2PolyBin*)pb.getPoly(sample1)->GetBins()->At(localbin)->Clone();
      double p    = pbin->GetXMax();
      double c    = pbin->GetYMax();
      double plow = pbin->GetXMin();
      double clow = pbin->GetYMin();

      std::cout << "p = " << p << " plow = " << plow << " c = " << c << " clow = " << clow << " sample = " << SampleId::ConvertSample(sample1) << std::endl;
      //WP loop over bins in this sample
      for(int k = pb.GetGlobalBinMatrixOffset_Det(sample1); k < pb.GetGlobalBinMatrixOffset_Det(sample1) + pb.GetNbins_Det(sample1); ++k){
	
        localbin = k - pb.GetGlobalBinMatrixOffset_Det(sample1);

	TH2PolyBin* pbin_det = (TH2PolyBin*)pb.getPoly_det(sample1)->GetBins()->At(localbin)->Clone();
	double p_det    = pbin_det->GetXMax();
	double c_det    = pbin_det->GetYMax();
	double plow_det = pbin_det->GetXMin();
	double clow_det = pbin_det->GetYMin();

        // If the fit bin does not overlap at all with the systematics obsnorm bin then skip it.
	//WP just get bin numbers again
        if (p_det <= plow || c_det <= clow || plow_det >= p || clow_det >= c)
          continue;

        std::cout << " plow_det = " << plow_det << " p_det = " << p_det << " clow_det = " << clow_det << " c_det = " << c_det << " local bin = " << localbin << " Det Bin = " << k << std::endl;

	//WP not doing but ok anyway
        if(staterr > mc_stat_error_redu->GetBinContent(k + 1)){
          mc_stat_error_redu->SetBinContent(k + 1, staterr);
        } 
        break;
      }
    }
    
    //WP I think ok
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
          cov->SetBinContent(iBin, jBin, cov->GetBinContent(iBin,jBin));// + mc_stat_error_redu->GetBinContent(iBin));//WP: Taking out so MC stats not added to cov. Yea take that out
        }
	
	matrix(iBin-1,jBin-1) = cov->GetBinContent(iBin,jBin);
        cov_nofake->SetBinContent(iBin, jBin, cov->GetBinContent(iBin,jBin));//WP ok i think
      }
      mc_stat_error_redu->SetBinContent(iBin, sqrt(mc_stat_error_redu->GetBinContent(iBin)));
    }

    //WP I think ok as well!
    for(int iBin = 1; iBin <= nBinsRedu; ++iBin){
      for(int jBin = 1; jBin <= nBinsRedu; ++jBin){
        if(cov->GetBinContent(iBin,iBin) > 0 && cov->GetBinContent(jBin,jBin) > 0)
          corr_nofake->SetBinContent(iBin, jBin, cov->GetBinContent(iBin,jBin) / sqrt(cov->GetBinContent(iBin,iBin) * cov->GetBinContent(jBin,jBin)));
        if(cov_var->GetBinContent(iBin,iBin) > 0 && cov_var->GetBinContent(jBin,jBin) > 0)
          corr_var->SetBinContent(iBin, jBin, cov_var->GetBinContent(iBin,jBin) / sqrt(cov_var->GetBinContent(iBin,iBin) * cov_var->GetBinContent(jBin,jBin)));
      }
    }

    //WP same as above if no MC so ok
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
      //WP these ok, just get sample from global det bin
      std::string SampleStr = ConvertSample(pb.GetSampleFromIndex_Det(iBin));
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

    //WP these ok, just need get sample from global bin
    for (int iBin = 1; iBin < nBinsFull+1; ++iBin) {
      std::string SampleStr = ConvertSample(pb.GetSampleFromIndex(iBin));
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

    for(int iBin = 0; iBin < nBinsRedu; ++iBin){

      c2.SetName(Form("Bin_%d_NEvent_Spread",iBin));
      c2.SetTitle(Form("Bin_%d_NEvent_Spread",iBin));

      int nbins = nominal->GetBinContent(iBin+1)*80;
      if (nbins > 1000) nbins = 1000;
      double lower = 0.0;
      double upper = nominal->GetBinContent(iBin+1)*3.0;
      std::string title = Form("Bin %d", iBin);
      // Get the sample
      SampleId::SampleEnum sample = SampleId::kUnassigned;
      int bin = -1;
      //WP get local bin
      pb.GetLocalBinMatrix_Det(iBin, sample, bin);
      if (sample != SampleId::kUnassigned) {
	//WP get axis wont work, but we can just get polybin edges I think
	TH2PolyBin* pbin_det = (TH2PolyBin*)pb.getPoly_det(sample)->GetBins()->At(bin)->Clone();
        double momlow = pbin_det->GetXMin();
        double momhi  = pbin_det->GetXMax();
        double coslow = pbin_det->GetYMin();
        double coshi  = pbin_det->GetYMax();
        std::stringstream ss;
        ss << momlow << "-" << momhi;
        std::stringstream ss2;
        ss2 << coslow << "-" << coshi;
        title += " (" + SampleId::ConvertSample(sample) + ", {p_{#mu}, cos#theta_{#mu}} = {" + ss.str() + ", " + ss2.str() + "})";
      }

      TH1F h1(Form("Bin_%d_NEvent_Spread_Hist", iBin), title.c_str(), nbins, lower, upper);

      TF1 f1("f1","gaus",nominal->GetBinContent(iBin+1)*0.5, nominal->GetBinContent(iBin+1)*2.0);
      TF1 f2("f2","gaus",nominal->GetBinContent(iBin+1)*0.5, nominal->GetBinContent(iBin+1)*2.0);

      f1.SetNameTitle(Form("Bin_%d_NEvent_Varying_Gaussian",iBin), "Varying Gaussian");
      f2.SetNameTitle(Form("Bin_%d_NEvent_Total_Gaussian",iBin), "Total Gaussian");

      for(int iToy = 0; iToy < nThrows; ++iToy){
        h1.Fill(varied[iToy]->GetBinContent(iBin+1));
      }
      //WP can still do this ok
      h1.GetXaxis()->SetTitle("N_{Events}^{MC} selected");
      h1.GetYaxis()->SetTitle("N_{toys}");
      h1.SetStats(kFALSE);
      h1.Draw();

      f1.SetParameters(h1.Integral()/(sqrt(2.0*TMath::Pi()*cov_var->GetBinContent(iBin+1,iBin+1))*nominal->GetBinContent(iBin+1)),
          var_mean->GetBinContent(iBin+1),
          sqrt(cov_var->GetBinContent(iBin+1,iBin+1))*nominal->GetBinContent(iBin+1));

      f2.SetParameters(h1.Integral()/(sqrt(2.0*TMath::Pi()*cov->GetBinContent(iBin+1,iBin+1))*nominal->GetBinContent(iBin+1)),
          var_mean->GetBinContent(iBin+1),
          sqrt(cov->GetBinContent(iBin+1,iBin+1))*nominal->GetBinContent(iBin+1));
      //WP think ok
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
      int lowbin = h1.FindBin(nominal->GetBinContent(iBin+1)*0.5);
      int hibin = h1.FindBin(nominal->GetBinContent(iBin+1)*1.5);

      GOF->SetBinContent(iBin+1,chi2);
      GOFnDOF->SetBinContent(iBin+1, chi2/h1.GetNbinsX());
      DiffMean->SetBinContent(iBin+1,(FitMean-CalcMean)/nominal->GetBinContent(iBin+1));
      DiffError->SetBinContent(iBin+1,(FitError-CalcError)/nominal->GetBinContent(iBin+1));

      TLine l1(nominal->GetBinContent(iBin+1),0,nominal->GetBinContent(iBin+1),h1.GetMaximum());
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
      leg->AddEntry(&l1, Form("Nominal bin content, %.1f events", nominal->GetBinContent(iBin+1)), "l");
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
    TDirectory* PolyCovDir = OutputFile->mkdir("PolyForCov");
    OutputFile->cd();
    TDirectory* PolyFitDir = OutputFile->mkdir("PolyForFit");

    //WP no axis so maybe just don't do this? will have to handle in mach3. Let's save the poly templates 
    for(int i = SampleId::kUnassigned; i < SampleId::kNSamples; i++){
      SampleId::SampleEnum smpl = static_cast<SampleId::SampleEnum>(i);

      if(pb.IsActiveSample(smpl)){
        std::string name = SampleId::ConvertSample(smpl);
        std::replace (name.begin(), name.end(), ' ', '_');
        PolyCovDir->cd();
        pb.getPoly_det(smpl)->Write((name+"_poly").c_str());
        PolyFitDir->cd();
        pb.getPoly(smpl)->Write((name+"_poly").c_str());
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


