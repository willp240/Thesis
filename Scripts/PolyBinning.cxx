#include <PolyBinning.hxx>

// The static member pointer to the singleton.
BANFF::PolyBinning* BANFF::PolyBinning::fPolyBinning = NULL;


BANFF::PolyBinning::PolyBinning(){
  Do4PiFHC     = (bool)ND::params().GetParameterI("BANFF.Do4PiFHC");
  DoNue        = (bool)ND::params().GetParameterI("BANFF.DoNueSelections");
  DoOnlyNue    = (bool)ND::params().GetParameterI("BANFF.DoOnlyNueSelections");
  DoMultiPiRHC = (bool)ND::params().GetParameterI("BANFF.DoMultiPiRHC");

  Do1DCheckMom       = (bool)ND::params().GetParameterI("BANFF.Do1DCheckMom");
  Do1DCheckCos       = (bool)ND::params().GetParameterI("BANFF.Do1DCheckCos");
  UseFlatNuMuBinning = (bool)ND::params().GetParameterI("BANFF.UseFlatNuMuBinning");

  if(Do1DCheckMom && Do1DCheckCos){
    std::cerr << "You cannot check 1D Cov Mat for cos and mom in the same time." << std::endl;
    std::cerr << "Change the parameter file (Do1DCheck)" << std::endl;
    throw;
  }

  FullToReduMap.clear();
}

void BANFF::PolyBinning::SetupPoly(){ 

  std::string polyFile =  "inputs/convertedtoPolyJacobsBinning.root";  
  TFile* fp = TFile::Open(polyFile.c_str());
  if(fp==NULL){
    std::cerr << "Input binning file not found. " << std::endl;
    throw;
  } 

  for(int i = SampleId::kUnassigned; i < SampleId::kNSamples; i++){
    SampleId::SampleEnum sample = static_cast<SampleId::SampleEnum>(i);
    std::string samplename = ConvertSample(sample);
    ActiveSample[sample]=false;
    std::string temp = samplename;
    while (temp.find(" ") != std::string::npos) {
      temp.replace(temp.find(" "), 1, std::string("_"));
    }

    if(fp->GetListOfKeys()->Contains(temp.c_str())){
      poly[sample]=(TH2Poly*)(fp->Get(temp.c_str()))->Clone();
      poly[sample]->Reset("");
      ActiveSample[sample]=true;
      std::cout << temp << " Active in poly file " << std::endl;
    }  
  }

  std::string polyFile_det = "inputs/convertedtoPolyJacobsDetBinning.root";
  TFile* fp_det = TFile::Open(polyFile_det.c_str());
  if(fp_det==NULL){
    std::cerr << "Input binning file for det not found. " << std::endl;
    throw;
  }
  
  for(int i = SampleId::kUnassigned; i < SampleId::kNSamples; i++){
    SampleId::SampleEnum sample = static_cast<SampleId::SampleEnum>(i);
    std::string samplename = ConvertSample(sample);
    std::string temp = samplename;
    while (temp.find(" ") != std::string::npos) {
      temp.replace(temp.find(" "), 1, std::string("_"));
    }

    if(ActiveSample[sample]==true){
      if(fp_det->GetListOfKeys()->Contains(temp.c_str())){
	poly_det[sample]=(TH2Poly*)fp_det->Get(temp.c_str())->Clone();
      }
      else{
	std::cout << "Sample " << samplename << " active in input poly file, but not input poly for det file. " << std::endl;
	throw;
      }	  
      poly_det[sample]->Reset("");
    }
  }
}


int BANFF::PolyBinning::GetGlobalBinMatrixMomCos_Det(SampleId::SampleEnum sample, double Mom, double Cos){
  /// This gets the global bin number in the detector matrix
  int offset = 0;
  if(Mom < 0     || Cos < -1) return -1;
  if(Mom > 30000 || Cos >  1) return -1;

  for (std::map<SampleId::SampleEnum, TH2Poly*>::const_iterator it=poly_det.begin(); it!=poly_det.end(); ++it){
        
    if(sample == (*it).first){
      TH2Poly* samplePoly = poly_det[(*it).first];
      int bin = samplePoly->FindBin(Mom,Cos);
      int nbins = samplePoly->GetNumberOfBins();
      if(bin==0 || bin>nbins){ return -1; }
      return offset + bin;
    }else{
      offset = offset + GetNbins_Det((*it).first);
    }
  }
  return -1; 
}


 int BANFF::PolyBinning::GetGlobalBinMatrixMomCos(SampleId::SampleEnum sample, double Mom, double Cos){
   /// This gets the global bin number in the detector matrix
   int offset = 0;
   if(Mom < 0     || Cos < -1) return -1;
   if(Mom > 30000 || Cos >  1) return -1;

   for (std::map<SampleId::SampleEnum, TH2Poly*>::const_iterator it=poly_det.begin(); it!=poly.end(); ++it){

     if(sample == (*it).first){
       TH2Poly* samplePoly = poly[(*it).first];
       int bin = samplePoly->FindBin(Mom,Cos);
       int nbins = samplePoly->GetNumberOfBins();
       if(bin==0 || bin>nbins){ return -1; }
       return offset + bin;
     }else{
       offset = offset + GetNbins((*it).first);
     }
   }
   return -1;
 }

void BANFF::PolyBinning::GetLocalBinMatrix(const int GlobalBin, SampleId::SampleEnum& Sample, int& Bin){
  
  Sample = SampleId::kUnassigned;
  Bin = -1;
    
  if(GlobalBin > GetNbins()){
    return;
  }
  Sample = GetSampleFromIndex(GlobalBin);
  
  // Gets the first index in the matrix with this sample
  int Ind = GetGlobalBinMatrixOffset(Sample);
  int LocalIndex = GlobalBin - Ind;
  Bin = LocalIndex; 

  return;
}

 void BANFF::PolyBinning::GetLocalBinMatrix_Det(const int GlobalBin, SampleId::SampleEnum& Sample, int& Bin){

   Sample = SampleId::kUnassigned;
   Bin = -1;

   if(GlobalBin > GetNbins_Det()){
     return;
   }
   Sample = GetSampleFromIndex_Det(GlobalBin);

   // Gets the first index in the matrix with this sample
   int Ind = GetGlobalBinMatrixOffset_Det(Sample);
   int LocalIndex = GlobalBin - Ind;
   Bin = LocalIndex;

   return;
 }


int BANFF::PolyBinning::GetGlobalBinMatrixOffset_Det(const SampleId::SampleEnum sample){
  /// This gets the global bin number in the detector matrix
  int row = 0;

  for (std::map<SampleId::SampleEnum, TH2Poly*>::iterator it=poly_det.begin(); it!=poly_det.end(); ++it){
    if(sample == it->first)
      return row;
    row = row + GetNbins_Det(it->first);
  }
  return -1; 
}


int BANFF::PolyBinning::GetGlobalBinMatrixOffset(const SampleId::SampleEnum sample){
  /// This gets the global bin number in the matrix
  int row = 0;

  for (std::map<SampleId::SampleEnum, TH2Poly*>::iterator it=poly.begin(); it!=poly.end(); ++it){
    if(sample == it->first)
      return row;
    row = row + GetNbins(it->first);
  }
  return -1; 

}


int BANFF::PolyBinning::GetNbins_Det(){
  /// This gets the total number of bins in the detector matrix
  int row = 0;

  for (std::map<SampleId::SampleEnum, TH2Poly*>::iterator it=poly_det.begin(); it!=poly_det.end(); ++it){
    row = row + GetNbins_Det(it->first);
  }
  return row; 

}


int BANFF::PolyBinning::GetNbins(){
  /// This gets the total number of bins in the matrix
  int row = 0;

  for (std::map<SampleId::SampleEnum, TH2Poly*>::iterator it=poly.begin(); it!=poly.end(); ++it){
    row = row + GetNbins(it->first);
  }
  return row; 

}

SampleId::SampleEnum BANFF::PolyBinning::GetSampleFromIndex(const int index){

  int row = 0;

  for (std::map<SampleId::SampleEnum, TH2Poly*>::iterator it=poly.begin(); it!=poly.end(); ++it){
    row = row + GetNbins(it->first);
    if(index < row)
      return it->first;
  }
  return SampleId::kUnassigned; 
}

SampleId::SampleEnum BANFF::PolyBinning::GetSampleFromIndex_Det(const int index){

  int row = 0;

  for (std::map<SampleId::SampleEnum, TH2Poly*>::iterator it=poly_det.begin(); it!=poly_det.end(); ++it){
    row = row + GetNbins_Det(it->first);
    if(index < row)
      return it->first;
  }
  return SampleId::kUnassigned; 
  
}

std::vector<std::string> BANFF::PolyBinning::GetListOfBins(){

  std::vector<std::string> BinList;
  
  for(int i = SampleId::kUnassigned; i < SampleId::kNSamples; i++){
    SampleId::SampleEnum sample = static_cast<SampleId::SampleEnum>(i);
    if(ActiveSample[sample]){
      std::string samplename = ConvertSample(sample);
      for(int ibin = 0; ibin < GetNbins(sample); ibin++)
	BinList.push_back(std::string(Form("%s_%d", samplename.c_str(), ibin)));
    }
  }
  
  return BinList;
  
}

std::vector<std::string> BANFF::PolyBinning::GetListOfBins_Det(){

  std::vector<std::string> BinList;
  
  for(int i = SampleId::kUnassigned; i < SampleId::kNSamples; i++){
    SampleId::SampleEnum sample = static_cast<SampleId::SampleEnum>(i);
    if(ActiveSample[sample]){
      std::string samplename = ConvertSample(sample);
      for(int ibin = 0; ibin < GetNbins_Det(sample); ibin++)
	BinList.push_back(std::string(Form("%s_%d", samplename.c_str(), ibin)));
    }
  }
  
  return BinList;
  
}

// A map to be able to go from full -> reduced binning
std::map<int,int> BANFF::PolyBinning::GetFullToReduMap(){

  if(FullToReduMap.size() > 0)
    return FullToReduMap;
  
  //int nBinsRedu = GetNbins_Det();
  int nBinsFull = GetNbins();

  //Start with some check that the actual binning is ok.
  for(int j = 0; j < nBinsFull; ++j){
    int bin;
    SampleId::SampleEnum sample;
    // First get the local bin (in the block of the sample)
    GetLocalBinMatrix(j, sample, bin);
    TH2PolyBin* pbin = (TH2PolyBin*)poly[sample]->GetBins()->At(j+1)->Clone();
    double ph = pbin->GetXMax();
    double pl = pbin->GetXMin();
    double ch = pbin->GetYMax();
    double cl = pbin->GetYMin();

    // Check that these match up they should be the same
    if(j != GetGlobalBinMatrixMomCos(sample, 0.5*(ph+pl), 0.5*(ch+cl))){
      std::cout << "Closure test failed!! " << std::endl;
      std::cout << "Global Bin: " << j << std::endl;
      std::cout << "Should be identical to GetGlobalBinMatrixMomCos(" << sample << ", " << 0.5*(ph+pl) << ", "<< 0.5*(ch+cl) << ") " <<  GetGlobalBinMatrixMomCos(sample, 0.5*(ph+pl), 0.5*(ch+cl)) << std::endl;
      // If they dont exit
      throw;
    }
    if(j != GetGlobalBinMatrixOffset(sample)+bin){
      std::cout << "Closure test failed!! " << std::endl;
      std::cout << "Global Bin: " << j << std::endl;
      std::cout << "Should be identical to GetGlobalBinMatrix(" << sample << ") + bin: " <<  GetGlobalBinMatrixOffset(sample) << " + " << bin << std::endl;
      // If they dont exit
      throw;
    }

    double epsilon = 0.00001; // Make sure this is smaller that the cos binning!!
    int bin_det  = GetGlobalBinMatrixMomCos_Det(sample, ph-epsilon, ch-epsilon);
    int bin2_det = GetGlobalBinMatrixMomCos_Det(sample, pl+epsilon, cl+epsilon);
    int localbin_det;
    SampleId::SampleEnum sample_det;
    GetLocalBinMatrix_Det(bin_det, sample_det, localbin_det);
    
    if(bin_det != bin2_det){
      std::cerr << "The fine binning is more coarse than the detector binning!!!" << std::endl;
      std::cerr << "Maybe epilson is too large for the Cos binning?" << std::endl;
      
      throw;
    }

    TH2PolyBin* pbin_det = (TH2PolyBin*)poly_det[sample]->GetBins()->At(j+1)->Clone();
    double ph_det = pbin_det->GetXMax();
    double pl_det = pbin_det->GetXMin();
    double ch_det = pbin_det->GetYMax();
    double cl_det = pbin_det->GetYMin();
    
    if(bin_det != GetGlobalBinMatrixMomCos_Det(sample_det, 0.5*(ph_det+pl_det), 0.5*(ch_det+cl_det))){
      std::cout << "Closure test failed!! " << std::endl;
      std::cout << "Global Bin: " << bin_det << std::endl;
      std::cout << "Should be identical to GetGlobalBinMatrixMomCos_Det(" << sample_det << ", " << 0.5*(ph_det+pl_det) << ", "<< 0.5*(ch_det+cl_det) << ") " <<  GetGlobalBinMatrixMomCos_Det(sample_det, 0.5*(ph_det+pl_det), 0.5*(ch_det+cl_det)) << std::endl;
      throw;
    }
    if(bin_det != GetGlobalBinMatrixOffset_Det(sample_det)+localbin_det){
      std::cout << "Closure test failed!! " << std::endl;
      std::cout << "Global Bin: " << bin_det << std::endl;
      std::cout << "Should be identical to GetGlobalBinMatrix(" << sample_det << ", " << 0.5*(ph_det+pl_det) << ", "<< 0.5*(ch_det+cl_det) << ") " <<  GetGlobalBinMatrixOffset_Det(sample_det) << " + localbin_det" << std::endl;
      throw;
    }
    FullToReduMap[j] = bin_det;
  }
  return FullToReduMap;
}




