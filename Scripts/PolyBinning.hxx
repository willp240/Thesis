#include <map>
#include "TFile.h"
#include "TAxis.h"
#include "TH2Poly.h"
#include "SampleId.hxx"
#include "Parameters.hxx"


namespace BANFF{
  
  class PolyBinning{
 
  private:
    PolyBinning();
    PolyBinning(const PolyBinning& binning);
    static PolyBinning* fPolyBinning;

    std::map<SampleId::SampleEnum, bool> ActiveSample;
    
    std::map<SampleId::SampleEnum, TH2Poly*> poly;
    std::map<SampleId::SampleEnum, TH2Poly*> poly_det;

    std::map<int, int> FullToReduMap;

    bool Do4PiFHC    ;
    bool DoNue       ;
    bool DoOnlyNue   ;
    bool DoMultiPiRHC;

    bool Do1DCheckMom;
    bool Do1DCheckCos;
    bool UseFlatNuMuBinning;
    
  public:
    ~PolyBinning(){
      for (std::map<SampleId::SampleEnum, TH2Poly*>::iterator it=poly.begin(); it!=poly.end(); ++it){
        if(it->second) delete it->second;
        poly.erase(it);
      }

      for (std::map<SampleId::SampleEnum, TH2Poly*>::iterator it=poly_det.begin(); it!=poly_det.end(); ++it){
        if(it->second) delete it->second;
        poly_det.erase(it);
      }

    };
    int GetNbins    (const SampleId::SampleEnum sample){
      if(ActiveSample[sample]) return poly        [sample]->GetNumberOfBins();
      else return 0;
    };
    int GetNbins_Det(const SampleId::SampleEnum sample){
      if(ActiveSample[sample]) return poly_det    [sample]->GetNumberOfBins();
      else return 0;
    };

    void SetupPoly();
    
    int GetGlobalBinMatrixOffset_Det(const SampleId::SampleEnum sample);
    int GetGlobalBinMatrixOffset    (const SampleId::SampleEnum sample);
    void GetLocalBinMatrix          (const int GlobalBin, SampleId::SampleEnum& Sample, int& Bin);
    void GetLocalBinMatrix_Det      (const int GlobalBin, SampleId::SampleEnum& Sample, int& Bin);

    int GetNbins_Det();
    int GetNbins    ();
    
    int GetGlobalBinMatrixMomCos_Det(const SampleId::SampleEnum sample, const double Mom, const double Cos);
    int GetGlobalBinMatrixMomCos    (const SampleId::SampleEnum sample, const double Mom, const double Cos);

    SampleId::SampleEnum GetSampleFromIndex    (const int index);
    SampleId::SampleEnum GetSampleFromIndex_Det(const int index);

    TH2Poly* getPoly(const SampleId::SampleEnum sample){
      return poly[sample];
    };

    TH2Poly* getPoly_det(const SampleId::SampleEnum sample){
      return poly_det[sample];
    };

    bool IsActiveSample(const SampleId::SampleEnum sample){
      return ActiveSample[sample];
    };

    bool IsActiveSample(const int sample){
      if(sample >= SampleId::kUnassigned && sample < SampleId::kNSamples){
        SampleId::SampleEnum s = static_cast<SampleId::SampleEnum>(sample);
        return ActiveSample[s];
      }else return false;
    };

    std::vector<std::string> GetListOfBins();
    std::vector<std::string> GetListOfBins_Det();
    std::map<int,int> GetFullToReduMap();

    
    static BANFF::PolyBinning& Get() {
      if (!fPolyBinning) fPolyBinning = new BANFF::PolyBinning;
      return *(fPolyBinning);
    };
    
  };


}
