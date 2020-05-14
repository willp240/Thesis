#include "TH1D.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TFile.h"
#include "TKey.h"
#include "THStack.h"
#include "TROOT.h"
#include "TGaxis.h"
#include <iostream>
#include <algorithm>

/* This code throws the kinematic variables for each bin. This was originally done for looking at Eb peaks in fluctuated asimov data. So the idea is you've already ran ND280_MCMC_2019 with FAKEDATAFLUCTUATEDFIT in the config. So you already have some fluctuated asimov p-theta distributions (i.e like the normal asimovs, but with an integer number of events in each bin.

This code then runs over that file, and loops over every bin in every sample. It then throws p and theta inside that bin range N_bincontent number of times. That gives a set of p theta pairs (i.e events) so that we can use different binnings on the same fluctuated asimov dataset.

So then once we've got the p theta pairs, we bin them up in the poly binning we want, which is read in from a file.

Then you can do a fake data fit if you like, I know I do!
*/

void throwBins(std::string file0){

  //open up poly template file
  //TFile *polyfile = TFile::Open("../inputs/polybinning_minwidth100MeV_template.root");
  TFile *polyfile = TFile::Open("../inputs/polybinning_BANFFTH2Ds_template.root");
  //open th2d file with events
  TFile *flucfile = TFile::Open(file0.c_str());
  TFile *outfileNew = new TFile("th2dFlucAsimov_binnedTH2D2.root", "RECREATE");
  TFile *outfileOrig = new TFile("th2dFlucAsimov_binnedPoly2.root","RECREATE");

  bool debug =false;

  int nsamples=18;
  //define sample names
  std::string sampleName[18] = {"FGD1_numuCC_0pi", "FGD1_numuCC_1pi", "FGD1_numuCC_other", "FGD1_anti-numuCC_0pi", "FGD1_anti-numuCC_1pi", "FGD1_anti-numuCC_other", "FGD1_NuMuBkg_CC0pi_in_AntiNu_Mode", "FGD1_NuMuBkg_CC1pi_in_AntiNu_Mode", "FGD1_NuMuBkg_CCother_in_AntiNu_Mode", "FGD2_numuCC_0pi", "FGD2_numuCC_1pi", "FGD2_numuCC_other", "FGD2_anti-numuCC_0pi", "FGD2_anti-numuCC_1pi", "FGD2_anti-numuCC_other", "FGD2_NuMuBkg_CC0pi_in_AntiNu_Mode", "FGD2_NuMuBkg_CC1pi_in_AntiNu_Mode", "FGD2_NuMuBkg_CCother_in_AntiNu_Mode"};

  TRandom3 *rnd = new TRandom3();

  //loop over samples
  for(int s=0; s<nsamples; s++){

    std::cout << "Doing sample MC_" << sampleName[s] << std::endl;

    // clone polys from template file
    TH2Poly* fluchist = (TH2Poly*) flucfile->Get(("DATA_"+sampleName[s]).c_str());
    TH2Poly* rebinned = (TH2Poly*) polyfile->Get((sampleName[s]).c_str());
    rebinned->Reset("");

    // loop over TH2D bins 
    for(int bin=1; bin<=fluchist->GetNumberOfBins(); bin++){ 

      if(debug)
	std::cout << "Bin " << bin << " of " << fluchist->GetNumberOfBins() << std::endl;

      TH2PolyBin* pbin = (TH2PolyBin*)((TH2Poly*)(fluchist)->GetBins()->At(bin-1))->Clone();
      double maxmom = pbin->GetXMax();
      double minmom = pbin->GetXMin();
      double maxang = pbin->GetYMax();
      double minang = pbin->GetYMin();

      if(debug)
	std::cout << "Bin " << bin << " " << minmom << "-" << maxmom << ", " << minang << "-" << maxang << std::endl;

      // loop over events in bin
      for(int eve=0; eve<fluchist->GetBinContent(bin); eve++){
	
	// throw p and theta in the bin range
	double p = (double)((maxmom-minmom)*rnd->Rndm()) + minmom;
	double t = (double)((maxang-minang)*rnd->Rndm()) + minang;
	
	int newbin = rebinned->FindBin(p,t);

	if(newbin<=0 || newbin>rebinned->GetNumberOfBins()){
	  std::cout << "WARNING: Bin " << newbin << ". Eve " << eve << " " << p << ", " << t << std::endl;
	  return 0;
	}

	if(debug)
	  std::cout << "Eve " << eve << " " << p << ", " << t << " " << newbin << std::endl;

	// fill the poly
	rebinned->SetBinContent(newbin,rebinned->GetBinContent(newbin)+1);

      }
    }

    outfileNew->cd();
    rebinned->Write(("MC_"+sampleName[s]).c_str());
    outfileOrig->cd();
    fluchist->Write(("MC_"+sampleName[s]).c_str());
    if(debug)
      s=18;
  }
  // Write the polys to file (then fit th2ds and polys!)
 
  for(int s=0; s<nsamples; s++){
      
    std::cout << "Compare sample MC_" << sampleName[s] << std::endl;
    TH2Poly* fluchist2 = (TH2Poly*) flucfile->Get(("DATA_"+sampleName[s]).c_str());
    TH2Poly* rebinned2 = (TH2Poly*) outfileNew->Get(("MC_"+sampleName[s]).c_str());

    int integral1 =0, integral2 =0;
    for (int i1=1; i1<=fluchist2->GetNumberOfBins(); i1++){
      integral1+= fluchist2->GetBinContent(i1);
    }
    for(int i1=1; i1<=rebinned2->GetNumberOfBins(); i1++){
      integral2+= rebinned2->GetBinContent(i1);
    }
    std::cout << "Orig integral " << integral1 << " new integral " << integral2 << std::endl;

    for(int b=1; b<=fluchist2->GetNumberOfBins(); b++){
      
      TH2PolyBin* pbin2 = (TH2PolyBin*)((TH2Poly*)(fluchist2)->GetBins()->At(b-1))->Clone();
      double maxmom2 = pbin2->GetXMax();
      double minmom2 = pbin2->GetXMin();
      double maxang2 = pbin2->GetYMax();
      double minang2 = pbin2->GetYMin();
      double midmom = ((maxmom2-minmom2)/2) + minmom2;
      double midang = ((maxang2-minang2)/2) + minang2;

      int newbinnum = rebinned2->FindBin(midmom,midang);

      if(debug)
	std::cout << "Orig bin: " << b << " " << minmom2 << "-" << maxmom2 << " " << minang2 << "-" << maxang2 << " " << std::endl << midmom << " " << midang << " " << fluchist2->GetBinContent(b) << std::endl;

      //      if(debug)
      //std::cout << "Bin " << b << ": " << fluchist2->GetBinContent(b) << ", " << newbinnum << ": " << rebinned2->GetBinContent(newbinnum) << std::endl << std::endl;

      TH2PolyBin* pbin3 = (TH2PolyBin*)((TH2Poly*)(rebinned2)->GetBins()->At(newbinnum-1))->Clone();
      double maxmom3 = pbin3->GetXMax();
      double minmom3 = pbin3->GetXMin();
      double maxang3 = pbin3->GetYMax();
      double minang3 = pbin3->GetYMin();
      double midmom3 = ((maxmom3-minmom3)/2) + minmom3;
      double midang3 = ((maxang3-minang3)/2) + minang3;

      if(debug)
	std::cout << "New bin:" << newbinnum << " " << minmom3 << "-" << maxmom3 << " " << minang3 << "-" << maxang3 << " " << std::endl << midmom3 << " " << midang3 << " " << rebinned2->GetBinContent(newbinnum) << std::endl << std::endl;

      if(debug)
	s=18;

    }

  }

}
