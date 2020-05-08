// Compile with something like
//
// g++ `root-config --cflags` -g -o DrawComp DrawComp.cpp -I`root-config --incdir` `root-config --glibs --libs`
//
// Run on ND280 only fits to get xsec and beam posteriors prettified. Essentially a compiled makeTransferMatrixAll
//
// ./DrawComp will tell you how
//
// Should probably prettify this into something managable like classes, or at least some structs

#include <iostream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <algorithm>

#include "TObjArray.h"
#include "TChain.h"
#include "TFile.h"
#include "TBranch.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TLegend.h"
#include "TString.h"
#include "TStyle.h"
#include "TH1.h"
#include "TF1.h"
#include "TH2.h"
#include "TH2Poly.h"
#include "TGraphErrors.h"
#include "TVectorD.h"
#include "TColor.h"


// All right, let's re-write this in a faster more efficient way!
void convertDetCovBinNumbers() {

  TFile *histofile = TFile::Open("../inputs/574BinDetectorTH2DBinning.root");
  TFile *polyfile = TFile::Open("../inputs/convertedtoPoly574BinDetector.root");
  TFile *outfile = new TFile("detcovVector.root", "RECREATE");

  std::string samples[18] = {"FGD1_numuCC_0pi", "FGD1_numuCC_1pi", "FGD1_numuCC_other", "FGD2_numuCC_0pi", "FGD2_numuCC_1pi", "FGD2_numuCC_other", "FGD1_anti-numuCC_0pi", "FGD1_anti-numuCC_1pi", "FGD1_anti-numuCC_other", "FGD2_anti-numuCC_0pi", "FGD2_anti-numuCC_1pi", "FGD2_anti-numuCC_other", "FGD1_NuMuBkg_CC0pi_in_AntiNu_Mode", "FGD1_NuMuBkg_CC1pi_in_AntiNu_Mode", "FGD1_NuMuBkg_CCother_in_AntiNu_Mode", "FGD2_NuMuBkg_CC0pi_in_AntiNu_Mode;1", "FGD2_NuMuBkg_CC1pi_in_AntiNu_Mode", "FGD2_NuMuBkg_CCother_in_AntiNu_Mode"};

  TVectorT   <double> polyIndex(574);
  TVectorT   <double> th2dIndex(574);

  int bincount =0;
  int globalOffset=0;
  for (int sam=0; sam<18; sam++){

    std::cout << "sample " << samples[sam] << std::endl;

    TH2D* samplehist = (TH2D*)histofile->Get(("MC_"+samples[sam]).c_str());
    TH2Poly* samplepoly = (TH2Poly*)polyfile->Get((samples[sam]).c_str());

    int numxbins = samplehist->GetXaxis()->GetNbins();
    int numybins = samplehist->GetYaxis()->GetNbins();

    for (int binx=1; binx<=numxbins; binx++){
      for(int biny=1; biny<=numybins; biny++){

	double binmaxx = samplehist->GetXaxis()->GetBinUpEdge(binx);
	double binminx = samplehist->GetXaxis()->GetBinLowEdge(binx);
	double binmaxy = samplehist->GetYaxis()->GetBinUpEdge(biny);
	double binminy = samplehist->GetYaxis()->GetBinLowEdge(biny);
	std::cout << "Bin xmax: " << binmaxx << ", Bin xmin: " << binminx << " Bin ymax: " << binmaxy << ", Bin ymin: " << binminy << std::endl;
	
	double x = (binmaxx-binminx)/2 + binminx;
	double y = (binmaxy-binminy)/2 + binminy;

	int polybin = samplepoly->FindBin(x,y) - 1 + globalOffset;
	int histbin = (biny-1)*numxbins + (binx-1) + globalOffset;
	
	std::cout << x << ", " << y << ". poly bin: " << polybin << ", histbin: " << histbin << ", globalOffset: " << globalOffset << std::endl;

	polyIndex(histbin) = polybin;
	th2dIndex(polybin) = histbin;

      }
    }
    globalOffset+=numxbins*numybins;
  }
  
  // Write all the nice vectors
  outfile->cd();
  polyIndex.Write("TH2PolyIndex");
  th2dIndex.Write("TH2DIndex");
    
  outfile->Close();
}
 
