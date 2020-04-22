#include "ND280AnalysisUtils.hxx"
#include <TVector3.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <typeinfo>
#include "BaseDataClasses.hxx"
#include "DetectorDefinition.hxx"
#include <ND280BeamBunching.hxx>
#include "PionInteractionSystematic.hxx"
#include "Parameters.hxx"

ND280BeamBunching bunching;


//**************************************************
Int_t anaUtils::GetRunPeriod(Int_t run, Int_t subrun){
  //**************************************************

  if      (             run<6000)   return 0;//run1  water                
  else if (run>=6462  && run<7664)  return 1;//run2  water                
  else if (run>=7665  && run<7755)  return 2;//run2  air                  
  else if (run>=8000  && run<8550)  return 3;//run3b air => mc associated to be the one of run2                  
  else if (run>=8550  && run<8800)  return 4;//run3c air                  
  else if (run>=8983  && run<9414)  return 5;//run4  water                
  else if (run>=9426  && run<9800)  return 6;//run4  air                  
  else if ((run>=10252 && run<10334)              //                    + 10252_10 to 10333
      || (run>=10519 && run<10522)) return 7;//run5  water         + 10518_28 to 10521_13             
  else if (run>=10334 && run<10518) return 8;//run5  water (ANTINU) up to 10334_11 to 10518_9
  else if (run == 10518){
    if(subrun < 10) return 8; //run5  water (ANTINU)
    else return 7; //run5  water
  }
  else if (run>=10828 && run<10954) return 13; //run6  air
  else if (run == 10954){
    if(subrun > -1 && subrun < 10) return 13; //run6  air
    else return 9; // run6b  air (ANTINU) 
  }
  else if (run>=10955 && run<11240) return 9;//run6b  air (ANTINU)
  else if (run>=11305 && run<11443) return 13; //run6  air
  else if (run>=11449 && run<11492) return 10;//run6c  air (ANTINU)
  else if (run>=11492 && run<11564) return 11;//run6d  air (ANTINU)
  else if (run>=11619 && run<11687) return 12;//run6e  air (ANTINU)
  else if (run == 11687){
    if(subrun < 7) return 12;//run6e  air (ANTINU)
    else return 14; //run6  air
  }
  else if(run>=12080 && run<12556) return 15; //run7b water (ANTINU) 
  else if(run>=12563 && run<12590) return 16; //run7c water

  else if(run>=12683 && run<13159) return 17; //run8 water  
  else if(run>=13160 && run<13740) return 18; //run8 air

  //  NEUT
  else if (run>=90110000 && run<=90119999) return 0;  //run1 water          
  else if (run>=90210000 && run<=90219999) return 1;  //run2 water          
  else if (run>=90200000 && run<=90209999) return 2;  //run2 air            
  else if (run>=90300000 && run<=90300015) return 3;  //run3b air separated based on the pot of data (run3b/(run3b+run3c)=0.13542             
  else if (run>=90300016 && run<=90300110) return 4;  //run3c air separated based on the pot of data            
  else if (run>=90410000 && run<=90419999) return 5;  //run4 water          
  else if (run>=90400000 && run<=90409999) return 6;  //run4 air            
  else if (run>=80510000 && run<=80519999) return 8;  //run5 antinu-water
  else if (run>=80600000 && run<=80600159) return 9;  //run6b antinu-air
  else if (run>=80600160 && run<=80600219) return 10; //run6c antinu-air - have to split Run 6 due to different flux tunings for the different parts
  else if (run>=80600220 && run<=80600299) return 11; //run6d antinu-air
  else if (run>=80600300 && run<=80600399) return 12; //run6e antinu-air
  else if (run>=80730000 && run<=80730119) return 15; //run7b antinu-water 
  else if (run>=90830000 && run<=90830095) return 17; //run8 water
  else if (run>=90820000 && run<=90820129) return 18; //run8 air 
  else if (run>=90820130 && run<=90829999) return 18; //run8 air - set this as default, to catch any files that are not currently at TRIUMF 

  else if (run>=80307000 && run<=80307999) return 8;  //sand muons antinu
  else if (run>=90307000 && run<=90307999) return 4;  //sand muons neut
  else if (run>=92307000 && run<=92307999) return 4;  //sand muons old neut

  // GENIE
  else if (run>=91110000 && run<=91119999) return 0;  //run1 water   
  else if (run>=91210000 && run<=91219999) return 1;  //run2 water          
  else if (run>=91200000 && run<=91209999) return 2;  //run2 air            
  else if (run>=91300000 && run<=91300015) return 3;  //run3 air separated based on the pot of data (run3b/(run3b+run3c)=0.13542             
  else if (run>=91300016 && run<=91300110) return 4;  //run3 air separated based on the pot of data (run3b/(run3b+run3c)=0.13542            
  else if (run>=91410000 && run<=91419999) return 5;  //run4 water          
  else if (run>=91400000 && run<=91409999) return 6;  //run4 air            
  else if (run>=81510000 && run<=81519999) return 8;  //run5 antinu-water
  else if (run>=81600000 && run<=81600186) return 9;  //run6b antinu-air
  else if (run>=81600187 && run<=81600260) return 10; //run6c antinu-air - have to split Run 6 due to different flux tunings for the different parts
  else if (run>=81600261 && run<=81600369) return 11; //run6d antinu-air
  else if (run>=81600370 && run<=81600490) return 12; //run6e antinu-air
  else if (run>=81710000 && run<=81710124) return 15; //run7b antinu-water 
  else if (run>=91810000 && run<=91810095) return 17; //run8 water
  else if (run>=91800000 && run<=91800129) return 18; //run8 air
  else if (run>=91800130 && run<=91800999) return 18; //run8 air - set this as default, to catch any files that are not currently at TRIUMF

  else return -1;

}


//**************************************************
int anaUtils::GetSandMode(int run){
  //**************************************************

  if      (run>=80307000 && run<=80307999) return 0;	// antineutrino flux
  else if (run>=90307000 && run<=90307999) return 1;	// neutrino flux
  else if (run>=92307000 && run<=92307999) return 2;	// neutrino flux, old NEUT

  else return -1;
}


//**************************************************
bool anaUtils::IsRHC(int run){
  //**************************************************

  // First check for sand
  if     (GetSandMode(run)  ==  0) return true;
  // Now check for standard MC and Data
  else if(GetRunPeriod(run) ==  8 ||
      GetRunPeriod(run) ==  9 ||
      GetRunPeriod(run) == 10 ||
      GetRunPeriod(run) == 11 ||
      GetRunPeriod(run) == 12 ||
      GetRunPeriod(run) == 15) return true;

  return false;

}


//*****************************************************************************
int anaUtils::GetFGDMichelElectrons(const AnaEventB& event, const SubDetId::SubDetEnum det, AnaFgdTimeBinB** arr, bool prod5Cut){
  //*****************************************************************************  

  int nSel = 0;

  Int_t ibunch = event.Bunch;

  Float_t selTime = bunching.GetBunchCentralTime(event, ibunch); 
  if(selTime < 0) return 0;

  Float_t numWidthsReqInSorting = 4; 
  Float_t binBunchWidth = 25.; 
  Float_t binTimeOffsetToBunch = -13.;
  Float_t minRecBinQ = 0.; 
  Int_t minHitsFGD1 = 6;
  Int_t minHitsFGD2 = 5;

  if(prod5Cut){
    minRecBinQ = 200.0;
    minHitsFGD1 = 0;
    minHitsFGD2 = 0;
  }

  int index=-1;

  selTime += numWidthsReqInSorting*binBunchWidth;

  for( Int_t i = 0; i < event.nFgdTimeBins; i++ ){

    AnaFgdTimeBinB *FgdTimeBin = event.FgdTimeBins[i];

    if( det == SubDetId::kFGD1 && FgdTimeBin->NHits[0] == 0 ) continue; 
    if( det == SubDetId::kFGD2 && FgdTimeBin->NHits[1] == 0 ) continue;
    if( det == SubDetId::kFGD  && FgdTimeBin->NHits[0] == 0 && FgdTimeBin->NHits[1] == 0 ) continue;

    Float_t binTime = FgdTimeBin->MinTime - binTimeOffsetToBunch;

    bool inBunch = false;

    for( UInt_t k = 0 ; k < bunching.GetNBunches(event); k++ ){
      Float_t time =  bunching.GetBunchCentralTime(event, k); 
      if( TMath::Abs( binTime - time ) < numWidthsReqInSorting*binBunchWidth )  inBunch = true;
    }

    if( inBunch )  continue;

    //find bin with min raw charge and greater than the required number of hits
    if( binTime >= selTime && FgdTimeBin->RawChargeSum[0] > minRecBinQ && FgdTimeBin->NHits[0] > minHitsFGD1 &&  (det == SubDetId::kFGD || det == SubDetId::kFGD1) ){
      index = i;
    }
    if( binTime >= selTime && FgdTimeBin->RawChargeSum[1] > minRecBinQ && FgdTimeBin->NHits[1] > minHitsFGD2 && (det == SubDetId::kFGD || det == SubDetId::kFGD2) ){
      index = i;
    }

  }
  if(index > -1){
    arr[nSel] = event.FgdTimeBins[index];
    nSel++;
  }

  return nSel; 
}


// //*****************************************************************************
// bool anaUtils::IsFGD1Analysis(const AnaEventB& event){
// //*****************************************************************************

//   if(!event.Summary)              return false;
//   if(!event.Summary->EventSample) return false;

//   if((event.Summary->EventSample >= SampleId::kFGD1NuMuCC
//       && event.Summary->EventSample<=SampleId::kFGD1NuMuCCOther) ||
//      (event.Summary->EventSample >= SampleId::kFGD1AntiNuMuCC
//       && event.Summary->EventSample<=SampleId::kFGD1AntiNuMuCCnQE)||
//      (event.Summary->EventSample >= SampleId::kFGD1NuMuBkgInAntiNuModeCC
//       && event.Summary->EventSample<=SampleId::kFGD1NuMuBkgInAntiNuModeCCnQE))
//     return true;
//   return false;
// }
// //*****************************************************************************
// bool anaUtils::IsFGD2Analysis(const AnaEventB& event){
// //*****************************************************************************
//   if(!event.Summary)              return false;
//   if(!event.Summary->EventSample) return false;

//   if((event.Summary->EventSample   >= SampleId::kFGD2NuMuCC
//       && event.Summary->EventSample<= SampleId::kFGD2NuMuCCOther) || 
//      (event.Summary->EventSample   >= SampleId::kFGD2AntiNuMuCC 
//       && event.Summary->EventSample<= SampleId::kFGD2AntiNuMuCCNTracks) ||
//      (event.Summary->EventSample    >= SampleId::kFGD2NuMuBkgInAntiNuModeCC
//       && event.Summary->EventSample <= SampleId::kFGD2NuMuBkgInAntiNuModeCCnQE))
//     return true;
//   return false;
// }

//*********************************************************
double anaUtils::GetNucleonsPerNucleus(massComponentEnum massComponent, int target_nucleons) {
  //*********************************************************
  // target_nucleons: 1 for protons, 2 for neutrons, 3 for both

  double nNucleonsPerNucleus = 0., nProtonsPerNucleus=0;

  // Areal densities of each element in FGDs (=>relative abundance), taken from TN 91 and TN 198 (weighted averaged for FGD2)
  // those commented are the values used in TN 117 and TN 210, where 2nd and 3rd column are rounded fractions, i.e. divided by the total areal density
  //     element[4]   = {Nprotons,
  //                     Nnucleons (initialized later),
  //                     'typical' areal densities of XY modules in g/cm^2,
  //                     weighted averaged areal densities of water modules in mg/cm^2}
  double carbon[4]    = {6,  0, 1.849,  421.0704607 }; // {6,  0, 0.8613, 0.1505};
  double oxygen[4]    = {8,  0, 0.0794, (12397./6.) }; // {8,  0, 0.0370, 0.7383};
  double hydrogen[4]  = {1,  0, 0.1579, (1759./6.)  }; // {1,  0, 0.0736, 0.1048};
  double titanium[4]  = {22, 0, 0.0355, 0.          }; // {22, 0, 0.0165, 0.    };
  double silicon[4]   = {14, 0, 0.0218, 11          }; // {14, 0, 0.0101, 0.0039};
  double nitrogen[4]  = {7,  0, 0.0031, 0.          }; // {7,  0, 0.0014, 0.    };
  double magnesium[4] = {12, 0, 0.,     7.          }; // {12, 0, 0.,     0.0025};
  // totArealDensityXY = 2.1467 g/cm^2 summing, while in TN 91 is reported 2.147 and in NIM FGDs 2.1463
  double totArealDensityXY = carbon[2] + oxygen[2] + hydrogen[2] + titanium[2] + silicon[2] + nitrogen[2] + magnesium[2];
  double totArealDensityWater = carbon[3] + oxygen[3] + hydrogen[3] + titanium[3] + silicon[3] + nitrogen[3] + magnesium[3];

  // Nnucleons with isotopic abundance
  // values taken from www.ciaaw.org/isotopic-abundances.htm: Elements published in Pure and Applied Chemistry in 2011
  // those commented are the values used in TN 117 and TN 210
  carbon[1]    = 12 * 0.9893   + 13 * 0.0107;                                             // 12*0.989   + 13*0.011;
  oxygen[1]    = 16 * 0.99757  + 17 * 0.00038 + 18 * 0.00205;                             // 16*0.99762 + 17*0.00038 + 18*0.002;
  hydrogen[1]  = 1  * 0.999885 + 2  * 0.000115;                                           // 1*0.99985  + 2*0.00015;
  titanium[1]  = 46 * 0.0825   + 47 * 0.0744  + 48 * 0.7372  + 49 * 0.0541 + 50 * 0.0518; // 46*0.08    + 47*0.075   + 48*0.738  + 49*0.055 + 50*0.054;
  silicon[1]   = 28 * 0.92223  + 29 * 0.04685 + 30 * 0.03092;                             // 28*0.9222  + 29*0.0468  + 30*0.0309;
  nitrogen[1]  = 14 * 0.996346 + 15 * 0.00364;                                            // 14*0.99634 + 15*0.00366;
  magnesium[1] = 24 * 0.7899   + 25 * 0.1000  + 26 * 0.1101;                              // 24*0.7899  + 25*0.1     + 26*0.11;

  if(massComponent == kFGD1 || massComponent == kFGD2xymodules || massComponent == kFGD2xylike) {

    // average nucleons/protons per nucleus in FGD1
    nProtonsPerNucleus  = carbon[0]    * carbon[2] +
      oxygen[0]    * oxygen[2] +
      hydrogen[0]  * hydrogen[2] +
      titanium[0]  * titanium[2] +
      silicon[0]   * silicon[2] +
      nitrogen[0]  * nitrogen[2] +
      magnesium[0] * magnesium[2] ;
    nProtonsPerNucleus  = nProtonsPerNucleus / totArealDensityXY;

    nNucleonsPerNucleus = carbon[1]    * carbon[2] +
      oxygen[1]    * oxygen[2] +
      hydrogen[1]  * hydrogen[2] +
      titanium[1]  * titanium[2] +
      silicon[1]   * silicon[2] +
      nitrogen[1]  * nitrogen[2] +
      magnesium[1] * magnesium[2] ;
    nNucleonsPerNucleus = nNucleonsPerNucleus / totArealDensityXY;
  }
  else if(massComponent == kFGD2watermodules || massComponent == kFGD2waterlike) {

    // average nucleons/protons per nucleus in FGD2
    nProtonsPerNucleus  = carbon[0]    * carbon[3] +
      oxygen[0]    * oxygen[3] +
      hydrogen[0]  * hydrogen[3] +
      titanium[0]  * titanium[3] +
      silicon[0]   * silicon[3] +
      nitrogen[0]  * nitrogen[3] +
      magnesium[0] * magnesium[3] ;
    nProtonsPerNucleus = nProtonsPerNucleus / totArealDensityWater;

    nNucleonsPerNucleus = carbon[1]    * carbon[3] +
      oxygen[1]    * oxygen[3] +
      hydrogen[1]  * hydrogen[3] +
      titanium[1]  * titanium[3] +
      silicon[1]   * silicon[3] +
      nitrogen[1]  * nitrogen[3] +
      magnesium[1] * magnesium[3] ;
    nNucleonsPerNucleus = nNucleonsPerNucleus / totArealDensityWater;
  }
  else if (massComponent == kP0Dwater) {

    nProtonsPerNucleus = (8*1/3. +  // Oxygen
        1*2/3.);  // Hydrogen

    nNucleonsPerNucleus = ((16*0.99762+17*0.00038+18*0.002)*1/3. + // Oxygen
        (1*0.99985+2*0.00015)*2/3.); // Hydrogen
  }
  else {
    std::cout << "Unknown target component selected." << std::endl;
    return 0;
  }

  if     (target_nucleons == 1) return nProtonsPerNucleus;
  else if(target_nucleons == 2) return (nNucleonsPerNucleus - nProtonsPerNucleus);
  else if(target_nucleons == 3) return nNucleonsPerNucleus;
  else {
    std::cout << "Unknown target nucleons selected." << std::endl;
    return 0;
  }
}

//*********************************************************
double anaUtils::GetAreaFV(massComponentEnum massComponent) {
  //*********************************************************
  TVector3 fgd1min = anaUtils::ArrayToTVector3(DetDef::fgd1min);
  TVector3 fgd2min = anaUtils::ArrayToTVector3(DetDef::fgd2min);
  TVector3 fgd1max = anaUtils::ArrayToTVector3(DetDef::fgd1max);
  TVector3 fgd2max = anaUtils::ArrayToTVector3(DetDef::fgd2max);
  TVector3 FVdefminFGD1 = anaUtils::ArrayToTVector3(FVDef::FVdefminFGD1);
  TVector3 FVdefminFGD2 = anaUtils::ArrayToTVector3(FVDef::FVdefminFGD2);
  TVector3 FVdefmaxFGD1 = anaUtils::ArrayToTVector3(FVDef::FVdefmaxFGD1);
  TVector3 FVdefmaxFGD2 = anaUtils::ArrayToTVector3(FVDef::FVdefmaxFGD2);

  TVector3 FVfgd1 = (fgd1max - FVdefmaxFGD1)-(fgd1min + FVdefminFGD1);
  TVector3 FVfgd2 = (fgd2max - FVdefmaxFGD2)-(fgd2min + FVdefminFGD2);

  double areaFV = 0;
  if (massComponent == kFGD1) {
    areaFV = FVfgd1.x() * FVfgd1.y() / 100.;  // convert from mm^2 to cm^2
  } else if (massComponent == kFGD2xymodules || massComponent == kFGD2watermodules || massComponent == kFGD2waterlike || massComponent == kFGD2xylike) {
    areaFV = FVfgd2.x() * FVfgd2.y() / 100.;  // convert from mm^2 to cm^2
  }
  else {
    std::cout << "Unknown target component selected." << std::endl;
    return 0;
  }

  std::cout << "XY-plane area of the FV: " << areaFV << " cm^2" << std::endl;
  return areaFV;
}

//*********************************************************
double anaUtils::GetNTargets(massComponentEnum massComponent, int target_nucleons) {
  //*********************************************************
  // target_nucleons: 1 for protons, 2 for neutrons, 3 for both

  // get number of XY modules in the FV
  double nXYModulesFGD1FV = 15.; // whole volume
  double nXYModulesFGD2FV = 7.; // whole volume
  // avoid dealing with air between modules by looking only at FVdef and taking the integer of the division
  double discardedModulesFGD1 = (FVDef::FVdefminFGD1[2] + FVDef::FVdefmaxFGD1[2]) / DetDef::fgdXYModuleWidth;
  double discardedModulesFGD2 = (FVDef::FVdefminFGD2[2] + FVDef::FVdefmaxFGD2[2]) / DetDef::fgdXYModuleWidth;
  // multiply (and then divide) by 2 to allow approximation to half module
  int nFGD1 = (int)(discardedModulesFGD1 * 2) ;
  int nFGD2 = (int)(discardedModulesFGD2 * 2) ;
  nXYModulesFGD1FV -= (double)(nFGD1/2.);
  nXYModulesFGD2FV -= (double)(nFGD2/2.);

  // get number of water modules in the FV
  double xyplusair = DetDef::fgdXYModuleWidth + DetDef::fgdXYAirWidth;
  if (FVDef::FVdefminFGD2[2] > xyplusair || FVDef::FVdefmaxFGD2[2] > xyplusair)
    std::cout << "ERROR: FGD2 FV doesn't include the 6 water modules --> GetNTargets will be wrong for FGD2!" << std::endl;
  double nWaterModules = 6.;

  double massFV = 0.;
  if(massComponent == kFGD1 || massComponent == kFGD2xymodules) {

    double nXYModulesFV = 0;
    if (massComponent == kFGD1) {
      nXYModulesFV = nXYModulesFGD1FV;
      std::cout << "FGD1 (" << nXYModulesFV << " modules):" << std::endl;
    } else if (massComponent == kFGD2xymodules) {
      nXYModulesFV = nXYModulesFGD2FV;
      std::cout << "FGD2 XY (" << nXYModulesFV << " modules):" << std::endl;
    }

    // areal density from tab 1 in TN 91, even though the sum seems wrong (vertically 2.1463, horizontally 2.1467)
    // but air is missing (which is 0.000258 from TN 122)
    massFV = 2.147 * GetAreaFV(massComponent) * nXYModulesFV;
    //    double density = 2.147 / (((double)DetDef::fgdXYModuleWidth + (double)DetDef::fgdXYAirWidth) / 10.); // convert mm to cm
    //    std::cout << "   density (included air): " << density << " g/cm^3" << std::endl;
  }

  else if(massComponent == kFGD2watermodules) {
    std::cout << "FGD2 water modules (" << nWaterModules << " modules):" << std::endl;
    // areal density from tab 9 in TN 198 (weighted average of 'total'), as built
    // in tab 5 of TN 122 (bottom three rows) it is instead 2.791, but it is the MC one
    // (in both air is missing but from TN 198 it is only 0.000266011)
    //    massFV = 2.79867 * GetAreaFV(massComponent) * nWaterModules;
    //    double density = 2.79867 / (((double)DetDef::fgdWaterModuleWidth + (double)DetDef::fgdWaterAirWidth) / 10.); // convert mm to cm
    //    std::cout << "   density (included air): " << density << " g/cm^3" << std::endl;
    double arealDensity = 13.838 + 2.931;
    massFV = arealDensity * GetAreaFV(massComponent);
  }

  else if(massComponent == kFGD2waterlike) {
    std::cout << "FGD2 water-like (" << nWaterModules << " modules):" << std::endl;
    // areal density from tab 78 in TN 212, already for 6 modules, which is the MC one
    // it corresponds to the as built one, from tab 9 and 10 in TN 198 (weighted average of 'water' + 'virtual water')
    massFV = 13.838 * GetAreaFV(massComponent);
  }

  else if(massComponent == kFGD2xylike) { // this is XY + virtual XY
    std::cout << "FGD2 XY-like (" << nXYModulesFGD2FV << " modules):" << std::endl;
    // areal density from tab 78 in TN 212
    // it corresponds to the as built one, from tab 9 and 10 in TN 198 (weighted average of 'water' + 'virtual water')
    double arealDensity = (2.147 * nXYModulesFGD2FV) + 2.931;
    massFV = arealDensity * GetAreaFV(massComponent);
  }

  else if(massComponent == kFGD2) {
    double Nwater = GetNTargets(kFGD2watermodules,target_nucleons);
    double Nscint = GetNTargets(kFGD2xymodules,target_nucleons);
    return (Nwater + Nscint);
  }

  else if(massComponent == kFGDs) {
    double Nfgd1 = GetNTargets(kFGD1,target_nucleons);
    double Nfgd2 = GetNTargets(kFGD2,target_nucleons);
    return (Nfgd1 + Nfgd2);
  }

  else if (massComponent == kP0Dwater) {
    // Water mass from TN82v2 and p0dasbuilt spreadsheet
    massFV = 1902.*1.e3; // g
    std::cout << "P0D:" << std::endl;
  }

  else {
    std::cout << "Unknown target component selected." << std::endl;
    return 0;
  }

  std::cout << "   mass FV: " << massFV / 1000. << " kg" << std::endl; // convert g to kg

  // number of nucleons (assuming 1g/mol/nucleon)
  double nNucleons = massFV * units::Avogadro;
  double nTargets;
  if (target_nucleons == 3) nTargets = nNucleons;
  else {
    double nNucleonsPerNucleus = GetNucleonsPerNucleus(massComponent, 3);
    nTargets = nNucleons * GetNucleonsPerNucleus(massComponent, target_nucleons) / nNucleonsPerNucleus;
  }
  std::cout << "Number of targets: " << nTargets << std::endl;
  return nTargets;
}

//********************************************************************
void anaUtils::IncrementPOTBySpill(const AnaSpillB& spill, Header& header) {
  //********************************************************************

  const AnaBeamB& beam = *spill.Beam;
  const AnaEventInfoB& info = *spill.EventInfo;

  if (beam.POTSincePreviousSavedSpill > 1e+16) {
    std::cout << "WARNING: POT is suspiciously large for run " << info.Run << ", subrun " << info.SubRun << ", event " << info.Event << ": " << beam.POTSincePreviousSavedSpill << ". POT will not be counted for this event." << std::endl;
    return;
  }

  if (beam.POTSincePreviousSavedSpill < 0) {
    std::cout << "WARNING: POT is negative for run " << info.Run << ", subrun " << info.SubRun << ", event " << info.Event << ": " << beam.POTSincePreviousSavedSpill << ". POT will not be counted for this event." << std::endl;
    return;
  }

  // For real data
  if (!spill.GetIsMC()) {

    header._POTInfo[kNoCut]        += beam.POTSincePreviousSavedSpill;
    header._SpillInfo[kSpillNoCut] += beam.SpillsSincePreviousSavedSpill;

    if (beam.GoodSpill == 0) {
      header._POTInfo[kBadBeam]        += beam.POTSincePreviousSavedSpill;
      header._SpillInfo[kSpillBadBeam] += beam.SpillsSincePreviousSavedSpill;
      return;
    }

    if (!spill.DataQuality->GoodDaq) {
      header._POTInfo[kBadND280]        += beam.POTSincePreviousSavedSpill;
      header._SpillInfo[kSpillBadND280] += beam.SpillsSincePreviousSavedSpill;
      return;
    }


    if (beam.GoodSpill == 100) {
      header._POTInfo[k0KA] += beam.POTSincePreviousSavedSpill;
    } else if (beam.GoodSpill == 1) {
      header._POTInfo[k250KA] += beam.POTSincePreviousSavedSpill;
    } else if (beam.GoodSpill == 2) {
      header._POTInfo[k200KA] += beam.POTSincePreviousSavedSpill;
    } else if (beam.GoodSpill == -1) {
      header._POTInfo[km250KA] += beam.POTSincePreviousSavedSpill;
    } else {
      header._POTInfo[kOtherKA] += beam.POTSincePreviousSavedSpill;
    }

    header._POTInfo[kGoodBeamGoodND280]        += beam.POTSincePreviousSavedSpill;
    header._SpillInfo[kSpillGoodBeamGoodND280] += beam.SpillsSincePreviousSavedSpill;

    if (!spill.DataQuality->GoodTimeDaq) {
      header._POTInfo[kBadND280Time]        += beam.POTSincePreviousSavedSpill;
      header._SpillInfo[kSpillBadND280Time] += beam.SpillsSincePreviousSavedSpill;
      return;
    }


    if (beam.GoodSpill == 1) {
      header._POTInfo[k250KA_ND280Time] += beam.POTSincePreviousSavedSpill;
    } else if (beam.GoodSpill == -1) {
      header._POTInfo[km250KA_ND280Time] += beam.POTSincePreviousSavedSpill;
    }

    header._POTInfo[kGoodBeamGoodND280Time]        += beam.POTSincePreviousSavedSpill;
    header._SpillInfo[kSpillGoodBeamGoodND280Time] += beam.SpillsSincePreviousSavedSpill;


  }
  else{
    // For MC there is no information about magnet current, Spill and DQ are always good
    header._SpillInfo[kSpillNoCut]                  += beam.SpillsSincePreviousSavedSpill;
    header._SpillInfo[kSpillGoodBeamGoodND280]      += beam.SpillsSincePreviousSavedSpill;
    header._SpillInfo[kSpillGoodBeamGoodND280Time]  += beam.SpillsSincePreviousSavedSpill;
    
    header._POTInfo[kNoCut]                 += beam.POTSincePreviousSavedSpill;
    header._POTInfo[kGoodBeamGoodND280]     += beam.POTSincePreviousSavedSpill;
    header._POTInfo[kGoodBeamGoodND280Time] += beam.POTSincePreviousSavedSpill;
  
  }

}

//********************************************************************
Float_t anaUtils::ComputeInversePT(const AnaDetCrossingB &cross, bool entrance) {
  //********************************************************************

  Float_t p = entrance ? GetEntranceMomentum(cross): GetExitMomentum(cross);

  if (p <= 0) return -999; 

  Float_t ut = entrance ? sqrt(1-pow(cross.EntranceMomentum[0]/p,2)) :
    sqrt(1-pow(cross.ExitMomentum[0]/p,2));

  Float_t pt = ut * p;

  if (pt!= 0)
    return 1/pt;
  else
    return -999;
}

//********************************************************************
Float_t anaUtils::ComputeInversePTFlip(const AnaTrackB& track){
  //********************************************************************

  Float_t u_t = sqrt(1-pow(track.DirectionEnd[0],2));
  Float_t pt = u_t*(track.MomentumFlip);
  if (pt!=0)
    return 1/pt;
  else
    return -999;
}

//********************************************************************
Float_t anaUtils::ComputeMomentumFromInversePTFlip(const AnaParticleB &part, Float_t PTinv) {
  //********************************************************************

  if (PTinv==0) return -999;
  Float_t pt = 1/PTinv;  
  Float_t u_t = sqrt(1-pow(part.DirectionEnd[0],2));
  if (u_t==0) return -999;

  Float_t p = pt/u_t;
  return p;
}

