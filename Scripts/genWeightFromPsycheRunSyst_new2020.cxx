#include <stdlib.h>
#include <cstdlib>

#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TString.h"
#include "TGraph2DErrors.h"

#include "T2KReWeight.h"
#include "T2KSyst.h"

#include "T2KGenieReWeight.h" 
#include "T2KNeutReWeight.h"
#include "T2KJNuBeamReWeight.h"
#include "T2KNIWGReWeight.h"
#include "T2KNIWGUtils.h"

#include "SampleId.hxx"

using namespace t2krew;

std::string FileInStr;
std::string FileOutStr;
bool IndividualThrow;
bool Sand;
int nEvents;

double calcRPA(double Q2, double A, double B, double D, double E, double U) {
  // Callum's fitted nominals
  // A = 0.59 +/- 20%
  // B = 1.05 +/- 20%
  // D = 1.13 +/- 15%
  // E = 0.88 +/- 40%
  // U = 1.2

  // Kept for convenience
  double eRPA = 1.;
 
  // Q2 transition; less than U -> polynominal
  if (Q2 < U) {
    // xprime as prescribed by Callum
    double xprime = Q2/U;
    double C      = D + U*E*(D-1)/3.;
    eRPA          = A*(1-xprime)*(1-xprime)*(1-xprime) + 3*B*(1-xprime)*(1-xprime)*xprime + 3*C*(1-xprime)*xprime*xprime + D*xprime*xprime*xprime;
  } else {
    eRPA = 1 + (D-1)*exp(-E*(Q2-U));
  }
  return eRPA;
}

double calcRPA(ND::NRooTrackerVtx *vtx=NULL) {
  if(vtx != NULL){
    if(fabs(atoi(((vtx->EvtCode)->String()).Data() )) != 1) return 1.;
    if(vtx->StdHepN < 3) return 1;
    double Q2=T2KNIWGUtils::Q2(vtx);
    return calcRPA(Q2, 0.59, 1.05, 1.13, 0.88, 1.2);
  }else{
    return 1;
  }
}

double calcQ2Tuning(ND::NRooTrackerVtx *vtx=NULL) {

  if(vtx != NULL){
    double Q2=T2KNIWGUtils::Q2(vtx);
    
    if(Q2>=0.0 && Q2<0.05) {
      return 0.505;
    }
    else if (Q2>=0.05 && Q2<0.1) {
      return 0.305;
    }
    else if (Q2>=0.1 && Q2<0.15) {
      return 0.22;
    }
    else if (Q2>=0.15 && Q2<0.2) {
      return 0.11;
    }
    else if (Q2>=0.2 && Q2<0.25) {
      return 0.07;
    }
    else {
      return 1.0;
    }
  }
  else{
    return 1;
  }
}

void ReweightTree(){

  std::vector< t2krew::ET2KSyst > dials;
  t2krew::T2KReWeight rw; 
  rw.AdoptWghtEngine("neut_rw", new t2krew::T2KNeutReWeight());
  rw.AdoptWghtEngine("niwg_rw", new t2krew::T2KNIWGReWeight());

  // Uncertainties
  // CCQE:
  std::cout << "Setting the NIWG dials" << std::endl;
  //  rw.Systematics().Include(t2krew::kNIWG2014a_pF_C12);
  //rw.Systematics().Include(t2krew::kNIWG2014a_pF_O16);//WP gone
  //rw.Systematics().Include(t2krew::kNIWG2014a_Eb_C12);
  //rw.Systematics().Include(t2krew::kNIWG2014a_Eb_O16);
  rw.Systematics().Include(t2krew::kNIWGMEC_Norm_C12);
  rw.Systematics().Include(t2krew::kNIWGMEC_Norm_O16);  
  rw.Systematics().Include(t2krew::kNIWGMEC_q3Cut);  
  rw.Systematics().Include(t2krew::kNXSec_MaCCQE);
  //rw.Systematics().Include(t2krew::kNXSec_VecFFCCQE); // used to set MAQE to act according to for RFG MC (2) or SF MC (402)

  // CC and NC single pion resonance:
  rw.Systematics().Include(t2krew::kNXSec_CA5RES);
  rw.Systematics().Include(t2krew::kNXSec_MaRES);//WP now mares 
  rw.Systematics().Include(t2krew::kNXSec_BgSclRES);
  rw.Systematics().Include(t2krew::kNXSec_BgSclLMCPiBarRES);//WP new
  // nue/numu uncertainties
  rw.Systematics().Include(t2krew::kNXSec_SCCVecQE);
  rw.Systematics().Include(t2krew::kNXSec_SCCAxlQE);
  rw.Systematics().Include(t2krew::kNIWG2012a_ccnueE0);

  // All other CC and NC
  //rw.Systematics().Include(t2krew::kNIWG2012a_dismpishp);//WP Gone, add new dis
  rw.Systematics().Include(t2krew::kNIWG_DIS_BY_corr);
  rw.Systematics().Include(t2krew::kNIWG_MultiPi_BY_corr);
  rw.Systematics().Include(t2krew::kNIWG_MultiPi_Xsec_AGKY);
  rw.Systematics().Include(t2krew::kNIWG2012a_cccohE0);
  rw.Systematics().Include(t2krew::kNIWG2012a_nccohE0);
  rw.Systematics().Include(t2krew::kNIWG2012a_ncotherE0); 

  //  rw.Systematics().Include(t2krew::kNIWG2014a_SF_RFG);//WP gone ?
  rw.Systematics().Include(t2krew::kNIWG_rpaCCQE_norm);
  rw.Systematics().Include(t2krew::kNIWG_rpaCCQE_shape);

  // Add FSI dials
  rw.Systematics().Include(t2krew::kNCasc_FrAbs_pi);
  rw.Systematics().Include(t2krew::kNCasc_FrCExLow_pi);
  rw.Systematics().Include(t2krew::kNCasc_FrInelLow_pi);
  rw.Systematics().Include(t2krew::kNCasc_FrPiProd_pi);
  rw.Systematics().Include(t2krew::kNCasc_FrCExHigh_pi);
  rw.Systematics().Include(t2krew::kNCasc_FrInelHigh_pi);

  //-- PDD Weights, New Eb dial
  rw.Systematics().Include(t2krew::kNIWGMEC_PDDWeight_C12);
  rw.Systematics().Include(t2krew::kNIWGMEC_PDDWeight_O16);
  //rw.Systematics().Include(t2krew::kNNucl_CCQEBindingEnergy_C12);
  //rw.Systematics().Include(t2krew::kNNucl_CCQEBindingEnergy_O16);

  //WP Add 2p2h E dep dials
  rw.Systematics().Include(t2krew::kNIWG_2p2hEdep_lowEnu);
  rw.Systematics().Include(t2krew::kNIWG_2p2hEdep_highEnu);
  rw.Systematics().Include(t2krew::kNIWG_2p2hEdep_lowEnubar);
  rw.Systematics().Include(t2krew::kNIWG_2p2hEdep_highEnubar);

  //WP add Q2 ?

  //  dials.push_back(t2krew::kNIWG2014a_pF_C12);//WP gone
  //dials.push_back(t2krew::kNIWG2014a_pF_O16);//WP gone
  //dials.push_back(t2krew::kNIWG2014a_Eb_C12);
  //dials.push_back(t2krew::kNIWG2014a_Eb_O16);
  dials.push_back(t2krew::kNIWGMEC_Norm_C12);
  dials.push_back(t2krew::kNIWGMEC_Norm_O16);  
  dials.push_back(t2krew::kNIWGMEC_q3Cut);  
  dials.push_back(t2krew::kNXSec_MaCCQE);
  //dials.push_back(t2krew::kNXSec_VecFFCCQE); // used to set MAQE to act according to for RFG MC (2) or SF MC (402)
  dials.push_back(t2krew::kNXSec_CA5RES);
  dials.push_back(t2krew::kNXSec_MaRES);
  dials.push_back(t2krew::kNXSec_BgSclRES);
  dials.push_bacl(t2krew::kNXSec_BgSclLMCPiBarRES);//WP new
  dials.push_back(t2krew::kNXSec_SCCVecQE);
  dials.push_back(t2krew::kNXSec_SCCAxlQE);
  dials.push_back(t2krew::kNIWG2012a_ccnueE0); 
  //  dials.push_back(t2krew::kNIWG2012a_dismpishp);//WP replace
  dials.push_back(t2krew::kNIWG_DIS_BY_corr);//WP check name
  dials.push_back(t2krew::kNIWG_MultiPi_BY_corr);//WP check name
  dials.push_back(t2krew::kNIWG_MultiPi_Xsec_AGKY);//WP check name
  dials.push_back(t2krew::kNIWG2012a_cccohE0);
  dials.push_back(t2krew::kNIWG2012a_nccohE0);
  dials.push_back(t2krew::kNIWG2012a_ncotherE0); 
  dials.push_back(t2krew::kNIWG2014a_SF_RFG);//WP sf->rfg stay for now
  dials.push_back(t2krew::kNIWG_rpaCCQE_norm);
  dials.push_back(t2krew::kNIWG_rpaCCQE_shape);
  dials.push_back(t2krew::kNCasc_FrAbs_pi);
  dials.push_back(t2krew::kNCasc_FrCExLow_pi);
  dials.push_back(t2krew::kNCasc_FrInelLow_pi);
  dials.push_back(t2krew::kNCasc_FrPiProd_pi);
  dials.push_back(t2krew::kNCasc_FrCExHigh_pi);
  dials.push_back(t2krew::kNCasc_FrInelHigh_pi);
  dials.push_back(t2krew::kNIWGMEC_PDDWeight_C12);
  dials.push_back(t2krew::kNIWGMEC_PDDWeight_O16);
  //dials.push_back(t2krew::kNNucl_CCQEBindingEnergy_C12);
  //dials.push_back(t2krew::kNNucl_CCQEBindingEnergy_O16);
  //WP Add 2p2h E dep dials
  dials.push_back(t2krew::kNIWG_2p2hEdep_lowEnu);
  dials.push_back(t2krew::kNIWG_2p2hEdep_highEnu);
  dials.push_back(t2krew::kNIWG_2p2hEdep_lowEnubar);
  dials.push_back(t2krew::kNIWG_2p2hEdep_highEnubar);

  
  // CCQE:
  //  rw.Systematics().SetAbsTwk(t2krew::kNIWG2014a_pF_C12);
  //rw.Systematics().SetAbsTwk(t2krew::kNIWG2014a_pF_O16);
  //rw.Systematics().SetAbsTwk(t2krew::kNIWG2014a_Eb_C12);
  //rw.Systematics().SetAbsTwk(t2krew::kNIWG2014a_Eb_O16);
  rw.Systematics().SetAbsTwk(t2krew::kNIWGMEC_Norm_C12);
  rw.Systematics().SetAbsTwk(t2krew::kNIWGMEC_Norm_O16);  
  rw.Systematics().SetAbsTwk(t2krew::kNIWGMEC_q3Cut);  
  rw.Systematics().SetAbsTwk(t2krew::kNXSec_MaCCQE);

  // CC and NC single pion resonance:
  rw.Systematics().SetAbsTwk(t2krew::kNXSec_CA5RES);
  rw.Systematics().SetAbsTwk(t2krew::kNXSec_MaRES);//WP gone
  rw.Systematics().SetAbsTwk(t2krew::kXSecTwkDial_BgSclRES);//wp add new
  rw.Systematics().SetAbsTwk(t2krew::kXSecTwkDial_BgSclLMCPiBarRES);//wp add new

  // nue/numu uncertainties
  rw.Systematics().SetAbsTwk(t2krew::kNXSec_SCCVecQE);
  rw.Systematics().SetAbsTwk(t2krew::kNXSec_SCCAxlQE);
  rw.Systematics().SetAbsTwk(t2krew::kNIWG2012a_ccnueE0); 

  // All other CC and NC
  rw.Systematics().SetAbsTwk(t2krew::kNIWG_DIS_BY_corr);
  rw.Systematics().SetAbsTwk(t2krew::kNIWG_MultiPi_BY_corr);
  rw.Systematics().SetAbsTwk(t2krew::kNIWG_MultiPi_Xsec_AGKY);
  rw.Systematics().SetAbsTwk(t2krew::kNIWG2012a_cccohE0);
  rw.Systematics().SetAbsTwk(t2krew::kNIWG2012a_nccohE0);
  rw.Systematics().SetAbsTwk(t2krew::kNIWG2012a_ncotherE0);

  rw.Systematics().SetAbsTwk(t2krew::kNIWG_rpaCCQE_norm);
  rw.Systematics().SetAbsTwk(t2krew::kNIWG_rpaCCQE_shape);

  // FSI dials
  rw.Systematics().SetAbsTwk(t2krew::kNCasc_FrAbs_pi);
  rw.Systematics().SetAbsTwk(t2krew::kNCasc_FrCExLow_pi);
  rw.Systematics().SetAbsTwk(t2krew::kNCasc_FrInelLow_pi);
  rw.Systematics().SetAbsTwk(t2krew::kNCasc_FrPiProd_pi);
  rw.Systematics().SetAbsTwk(t2krew::kNCasc_FrCExHigh_pi);
  rw.Systematics().SetAbsTwk(t2krew::kNCasc_FrInelHigh_pi);

  rw.Systematics().SetAbsTwk(t2krew::kNIWGMEC_Nu_High);//WP Check Names
  rw.Systematics().SetAbsTwk(t2krew::kNIWGMEC_Nubar_High);//WP Check Names
  rw.Systematics().SetAbsTwk(t2krew::kNIWGMEC_Nu_Low);//WP Check Names
  rw.Systematics().SetAbsTwk(t2krew::kNIWGMEC_Nubar_Low);//WP Check Names

  //-- PDD Weights, New Eb dial
  rw.Systematics().SetAbsTwk(t2krew::kNIWGMEC_PDDWeight_C12);
  rw.Systematics().SetAbsTwk(t2krew::kNIWGMEC_PDDWeight_O16);
  //rw.Systematics().SetAbsTwk(t2krew::kNNucl_CCQEBindingEnergy_C12);
  //rw.Systematics().SetAbsTwk(t2krew::kNNucl_CCQEBindingEnergy_O16);
  


  // Initialise dials to nominal
  for(unsigned int sys_iter1=0; sys_iter1<dials.size(); sys_iter1++){
    rw.Systematics().SetTwkDial(dials[sys_iter1], 0.);
  }
  // Apply nominal tuning
  //  rw.Systematics().SetTwkDial(t2krew::kNIWG2014a_SF_RFG,0); // SF->RFG tuning// WP gone
  //rw.Systematics().SetTwkDial(t2krew::kNXSec_VecFFCCQE, 2); // Note that when we set MAQE, we need to set NXSec_VecFFCCQE to 2 for SF->RFG MC //WP now don;t need this as SF?
  //  rw.Systematics().SetTwkDial(t2krew::kNIWG_rpaCCQE_norm, 1);//WP should these stay as they are? 
  //rw.Systematics().SetTwkDial(t2krew::kNIWG_rpaCCQE_shape,0);

  //WP Iso1/2, CA5, MaRES all different from neut nominal. 
  rw.Systematics().SetTwkDial(t2krew::kNXSec_MaCCQE, 0.85124);//WP tweak MAQE: 1.03/1.21=0.85124
  rw.Systematics().SetTwkDial(t2krew::kNXSec_CA5RES, 0.9505);//WP tweak CA5RES: 0.96/1.01=0.9505
  rw.Systematics().SetTwkDial(t2krew::kNXSec_BgSclRES,0.7385);//WP tweak I1/2: 0.96/1.3=0.7385
  rw.Systematics().SetTwkDial(t2krew::kNXSec_BgSclLMCPiBarRES,0.7385);//WP tweak I1/2 Low M Pi: 0.96/1.3-0.7385 


  std::cout << "Setting the NIWG dials -done-" << std::endl;
  rw.Reconfigure();
  std::cout << "Reconfigured NIWG dials -done-" << std::endl;


  TFile *FileIn = new TFile(FileInStr.c_str(), "READ");
  TTree* RTV = NULL;
  TClonesArray *nRooVtxs = new TClonesArray("ND::NRooTrackerVtx",100);
  int NRooVtx = 0;
  if(!Sand){
    RTV = (TTree*)FileIn->Get("NRooTrackerVtx");
    RTV->SetBranchAddress("Vtx",  &nRooVtxs);
    RTV->SetBranchAddress("NVtx", &NRooVtx);
  }
  TTree* all_syst = (TTree*)FileIn->Get("all");
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
  int    TrueVertexIDToy[2000];
  int    SelectionToy   [2000];
  double TrueEnuToy     [2000];
  int    TrueNuPDGToy   [2000];
  double LeptonMomToy   [2000];
  double LeptonCosToy   [2000];
  double WeightToy      [2000];
  double FluxWeightToy  [2000];
    
  all_syst->SetBranchAddress("Run"            , &Run            );
  all_syst->SetBranchAddress("SubRun"         , &SubRun         );
  all_syst->SetBranchAddress("EventNumber"    , &EventNumber    );
  all_syst->SetBranchAddress("SelectionNom"   , &SelectionNom   );
  all_syst->SetBranchAddress("TrueEnuNom"     , &TrueEnuNom     );
  all_syst->SetBranchAddress("TrueNuPDGNom"   , &TrueNuPDGNom   );
  all_syst->SetBranchAddress("TrueVertexIDNom", &TrueVertexIDNom);
  all_syst->SetBranchAddress("LeptonMomNom"   , &LeptonMomNom   );
  all_syst->SetBranchAddress("LeptonCosNom"   , &LeptonCosNom   );
  all_syst->SetBranchAddress("WeightNom"      , &WeightNom      );
  all_syst->SetBranchAddress("FluxWeightNom"  , &FluxWeightNom  );
  all_syst->SetBranchAddress("TrueVertexIDToy",  TrueVertexIDToy);
  all_syst->SetBranchAddress("SelectionToy"   ,  SelectionToy   );
  all_syst->SetBranchAddress("TrueEnuToy"     ,  TrueEnuToy     );
  all_syst->SetBranchAddress("TrueNuPDGToy"   ,  TrueNuPDGToy   );
  all_syst->SetBranchAddress("LeptonMomToy"   ,  LeptonMomToy   );
  all_syst->SetBranchAddress("LeptonCosToy"   ,  LeptonCosToy   );
  all_syst->SetBranchAddress("WeightToy"      ,  WeightToy      );
  all_syst->SetBranchAddress("FluxWeightToy"  ,  FluxWeightToy  );

  int prevRTV, rtvi=0, prevTruthID;
  double prevNIWGWeight;
  TFile* FileOut = new TFile(FileOutStr.c_str(), "RECREATE");
  TTree* NIWGWeightTree = new TTree("NIWGWeightTree", "NIWGWeightTree");
  double NIWGWeightNom = 1;
  double EnuNom        = 1;
  // double LeptonMom1p1hFDS = 1;//WP comment out fds
  double NIWGWeightToy[2000];
  double EnuToy[2000];
  int nToy=2000;
  int Toy[2000];

  // TFile* FileFDS_1p1h_carbon = new TFile("/neut/data20/jmgwalker/banff-work/T2KReWeight/app/numu_C12_MA1030_MDLQE0002_theta.root", "READ");//WP comment everything out about the 1p1h fake dataset
  //TFile* FileFDS_1p1h_oxygen = new TFile("/neut/data20/jmgwalker/banff-work/T2KReWeight/app/numu_O16_MA1030_MDLQE0002_theta.root", "READ");
  //  TGraph2DErrors* FDS_1p1h_carbon = (TGraph2DErrors*)FileFDS_1p1h_carbon->Get("nieves_minus_neut_graph");
  //TGraph2DErrors* FDS_1p1h_oxygen = (TGraph2DErrors*)FileFDS_1p1h_oxygen->Get("nieves_minus_neut_graph");
  
  NIWGWeightTree->Branch("T2KRW_NIWGWeightNom"   , &NIWGWeightNom   , "T2KRW_NIWGWeightNom/D");
  NIWGWeightTree->Branch("T2KRW_EnuNom"          , &EnuNom          , "T2KRW_EnuNom/D");
  //  NIWGWeightTree->Branch("T2KRW_LeptonMom1p1hFDS", &LeptonMom1p1hFDS, "T2KRW_LeptonMom1p1hFDS/D");
  NIWGWeightTree->Branch("T2KRW_nToys"           , &nToy            , "T2KRW_nToys/I");
  NIWGWeightTree->Branch("T2KRW_Toy"             ,  Toy             , "T2KRW_Toy[T2KRW_nToys]/I");
  NIWGWeightTree->Branch("T2KRW_NIWGWeightToy"   ,  NIWGWeightToy   , "T2KRW_NIWGWeightToy[T2KRW_nToys]/D");
  NIWGWeightTree->Branch("T2KRW_EnuToy"          ,  EnuToy          , "T2KRW_EnuToy[T2KRW_nToys]/D");
  
  for(int iEntry = 0; iEntry < (int)all_syst->GetEntries(); ++iEntry){
    NIWGWeightNom = 1;
    EnuNom = -999;
    if(iEntry%100==0) std::cout << "Progress " << (double)iEntry*100./(int)all_syst->GetEntries() << "%" << std::endl;
    for(int iToy = 0; iToy < 2000; ++iToy){
      Toy[iToy]           = -999;
      NIWGWeightToy[iToy] = 1;
      EnuToy       [iToy] = -999;
    }
    
    all_syst->GetEntry(iEntry);
    
    if(Sand){
      NIWGWeightTree->Fill();
      continue;
    }
    

    if(SelectionNom > SampleId::kUnassigned){
      bool foundvtx = false;
      ND::NRooTrackerVtx * vtx = NULL;
      if(TrueVertexIDNom != -999){
        prevTruthID = TrueVertexIDNom;
        prevRTV = rtvi;
        if(!foundvtx){
          while(vtx == NULL){
            //Pull out the correct RooTrackerVtx tree entry.
            RTV->GetEntry(rtvi);
            for(int i = 0; i < NRooVtx; ++i){
              vtx = (ND::NRooTrackerVtx*)(*nRooVtxs)[i];
              if(vtx->TruthVertexID == TrueVertexIDNom &&
                 fabs(vtx->StdHepP4[0][3]*1000 - TrueEnuNom) < 0.1){
                foundvtx = true;
                break;
              }
              vtx = NULL;
            }
            if(vtx==NULL) rtvi++;
            if(rtvi == RTV->GetEntries()){
              std::cout << "Looping to find correct vertex Nominal. VtxID: " << TrueVertexIDNom << std::endl;
              rtvi = 0;
              //break;
            }
            if(rtvi == prevRTV - 1){
              std::cout << "Event failed to find ANY vertex Nominal. VtxID: " << TrueVertexIDNom << std::endl;
              break;
            }
          }
        }else{
          //may not work, depending on how RTV->GetEntry(rvti) works
          for(int iPart = 0; iPart > NRooVtx; ++iPart){
            vtx = (ND::NRooTrackerVtx*)(*nRooVtxs)[iPart];
            if(vtx->TruthVertexID == TrueVertexIDNom &&
               fabs(vtx->StdHepP4[0][3]*1000 - TrueEnuNom) < 0.1){
              break;
            } 
          }
        }
        if(vtx != NULL){
          NIWGWeightNom = rw.CalcWeight(vtx);
          if(NIWGWeightNom>100 || NIWGWeightNom<0){
             std::cout << "NIWGWeightNom = " << NIWGWeightNom << " setting it to 1." << std::endl;
            NIWGWeightNom=1;
          }
          EnuNom        = vtx->StdHepP4[0][3]*1000;
          double MomDiff = 0;
          if(vtx->StdHepN > 2){
            int MuonIndex = T2KNIWGUtils::IdxFsl(*vtx);
            double mom = TMath::Sqrt(vtx->StdHepP4[MuonIndex][0]*vtx->StdHepP4[MuonIndex][0]+
                                     vtx->StdHepP4[MuonIndex][1]*vtx->StdHepP4[MuonIndex][1]+
                                     vtx->StdHepP4[MuonIndex][2]*vtx->StdHepP4[MuonIndex][2]);
            double angle = (TMath::ACos(vtx->StdHepP4[MuonIndex][2] / mom) * 180.0 / 3.14159265);
	    /* if      (vtx->StdHepPdg[1] == 1000060120 && TMath::Abs(vtx->StdHepPdg[MuonIndex]) == 13){
              MomDiff = FDS_1p1h_carbon->Interpolate(EnuNom,angle);
            }else if(vtx->StdHepPdg[1] == 1000080160 && TMath::Abs(vtx->StdHepPdg[MuonIndex]) == 13){
              MomDiff = FDS_1p1h_oxygen->Interpolate(EnuNom,angle);
	      }*/ //WP comment out 1p1h FDS
          }
	  //          LeptonMom1p1hFDS = LeptonMomNom + MomDiff;//WP comment out 1p1h1 fds
          //Add Kendall's pion tuning stuff
	  /*          double piEnergy = -1.0;
          if(fabs(atoi(((vtx->EvtCode)->String()).Data() )) == 16){
            for(int iPart = 0; iPart < vtx->StdHepN; ++iPart){
              if(fabs(vtx->StdHepPdg[iPart]) != 211) continue;
              piEnergy = vtx->StdHepP4[iPart][3];
            }
            if     (piEnergy < 0   ) { NIWGWeightNom = 1     * NIWGWeightNom; }
            else if(piEnergy < 0.25) { NIWGWeightNom = 0.135 * NIWGWeightNom; }
            else if(piEnergy < 0.5 ) { NIWGWeightNom = 0.4   * NIWGWeightNom; }
            else if(piEnergy < 0.75) { NIWGWeightNom = 0.294 * NIWGWeightNom; }
            else if(piEnergy < 1.0 ) { NIWGWeightNom = 1.206 * NIWGWeightNom; }
	    }*/
          NIWGWeightNom = NIWGWeightNom * calcQ2Tuning(vtx);
        }else{
          std::cout << "Didnt find vertex for this event nom" << std::endl;
        }
        prevNIWGWeight = NIWGWeightNom;
      }else{
        std::cout << "No vertex id saved for this event nom" << std::endl;
        NIWGWeightNom   =1;
        //LeptonMom1p1hFDS=LeptonMomNom; // WP comment out 1p1h fds
      }
    
      if(!IndividualThrow){
        NIWGWeightTree->Fill();
        continue;
      }
    }
    for(int iToy=0; iToy < 2000; ++iToy){
      Toy[iToy] = iToy;
      EnuToy[iToy] = -999;
      NIWGWeightToy[iToy] = 1;

      if(SelectionToy[iToy] > SampleId::kUnassigned){
        bool foundvtx = false;
        ND::NRooTrackerVtx * vtx = NULL;
        if(TrueVertexIDToy[iToy] != -999){
          prevTruthID = TrueVertexIDToy[iToy];
          prevRTV = rtvi;
          if(!foundvtx){
            while(vtx == NULL){
              //Pull out the correct RooTrackerVtx tree entry.
              RTV->GetEntry(rtvi);
              for(int i = 0; i < NRooVtx; ++i){
                vtx = (ND::NRooTrackerVtx*)(*nRooVtxs)[i];
                if(vtx->TruthVertexID == TrueVertexIDToy[iToy] &&
                   fabs(vtx->StdHepP4[0][3]*1000 - TrueEnuToy[iToy]) < 0.1){
                  foundvtx = true;
                  break;
                }
                vtx = NULL;
              }
              if(vtx==NULL) rtvi++;
              if(rtvi == RTV->GetEntries()){
                std::cout << "Looping to find correct vertex Toy" << iToy << ". VtxID: " << TrueVertexIDToy[iToy] << std::endl;
                rtvi = 0;
                //break;
              }
              if(rtvi == prevRTV - 1){
                std::cout << "Event failed to find ANY vertex Toy" << iToy << ". VtxID: " << TrueVertexIDToy[iToy] << std::endl;
                break;
              }
            }
          }else{
            //may not work, depending on how RTV->GetEntry(rvti) works
            for(int i = 0; i < NRooVtx; ++i){
              vtx = (ND::NRooTrackerVtx*)(*nRooVtxs)[i];
              if(vtx->TruthVertexID == TrueVertexIDToy[iToy] &&
                 fabs(vtx->StdHepP4[0][3]*1000 - TrueEnuToy[iToy]) < 0.1){
                break;
              } 
            }
          }
          if(vtx != NULL){
            NIWGWeightToy[iToy] = rw.CalcWeight(vtx);
            if(NIWGWeightToy[iToy]>100 || NIWGWeightToy[iToy] < 0){
              std::cout << "NIWGWeightToy[" << iToy << "] = " << NIWGWeightToy[iToy] << " setting it to 1." << std::endl;
              NIWGWeightToy[iToy]=1;
            }
            EnuToy       [iToy] = vtx->StdHepP4[0][3]*1000;
	    /*          
            //Add Kendall's pion tuning stuff
            double piEnergy = -1.0;
            if(fabs(atoi(((vtx->EvtCode)->String()).Data() )) == 16){
              for(int iPart = 0; iPart < vtx->StdHepN; ++iPart){
                if(fabs(vtx->StdHepPdg[iPart]) != 211) continue;
                piEnergy = vtx->StdHepP4[iPart][3];
              }
            }
            if     (piEnergy < 0   ) { NIWGWeightToy[iToy] = 1     * NIWGWeightToy[iToy]; }
            else if(piEnergy < 0.25) { NIWGWeightToy[iToy] = 0.135 * NIWGWeightToy[iToy]; }
            else if(piEnergy < 0.5 ) { NIWGWeightToy[iToy] = 0.4   * NIWGWeightToy[iToy]; }
            else if(piEnergy < 0.75) { NIWGWeightToy[iToy] = 0.294 * NIWGWeightToy[iToy]; }
            else if(piEnergy < 1.0 ) { NIWGWeightToy[iToy] = 1.206 * NIWGWeightToy[iToy]; }
	    */
            NIWGWeightToy[iToy] = NIWGWeightToy[iToy] * calcQ2Tuning(vtx);//WP replace with Q2
          }else{
            std::cout << "Didnt find vertex for this toy event" << std::endl;
          }
          prevNIWGWeight = NIWGWeightToy[iToy];
        }else{
          std::cout << "No vertex id saved for this event toy" << iToy << std::endl;
          NIWGWeightToy[iToy]   =1;
        }
      }
    }
    NIWGWeightTree->Fill();
  }
  FileOut->cd();
  NIWGWeightTree->Write();
  FileOut->Close();
  
}


// Print the cmd line syntax
void Usage(){
  std::cout << "Cmd line syntax should be:" << std::endl;
  std::cout << "genWeightFromPsycheRunSyst_new.exe -i RunSyst_New_file.root "
    "-o outputfile.root [-n nevents] [-s if running on sand]" << std::endl;
}

int ParseArgs(int argc, char **argv){

  IndividualThrow = true;
  Sand = false;
  for (;;) {
    int c = getopt(argc, argv, "sn:i:o:");
    if (c < 0)
      break;
    switch (c){
    case 'n':{
      std::istringstream tmp(optarg);
      tmp >> nEvents;
      std::cout << "Looping over " << nEvents << " events." << std::endl;
      break;
    }
    case 'i':{
      FileInStr = optarg;
      break;
    }  
    case 'o':{
      FileOutStr = optarg;
      break;
    }
    case 's':{
      Sand=true;
      break;
    }
    default:{
      std::cerr << "Option " << c << " is not recognised." << std::endl;
      Usage();
      throw;
      break;
    }
    }
  }
  return 0;
}


int main(int argc, char *argv[]){
  ParseArgs(argc,argv);
  ReweightTree();
  return 1;
}

