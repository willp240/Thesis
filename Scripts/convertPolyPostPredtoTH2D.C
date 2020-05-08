#include <TH2Poly.h>
#include <TH2D.h>
#include <vector>
#include <TMath.h>
#include <TFile.h>
#include <TTree.h>
#include <THStack.h>
#include <iostream> 

/*Script takes in a file of th2polys, prolly the output of ND280_MCMC_2019 when savenom is true in config, and turns it into a th2d for comparison with banff if the th2poly had uniform rectangular binning (this isn't checked so stringently so be careful). 
We loop over samples, then bins in the poly, filling a vector of x and y bin edges for each sample as we go. Then we make the th2d with those vectors, and finally loop over poly bins and fill the corresponding th2d bin.
 */

void convertPolyPostPredtoTH2D (std::string PolyFileName){
   
  std::string TH2DFileName = PolyFileName;
  TH2DFileName.replace(TH2DFileName.find(".root"), 5, std::string("_ConvertedToTH2D.root"));
  TFile* filePoly = TFile::Open(PolyFileName.c_str());
  TFile * fileTH2D = new TFile(TH2DFileName.c_str(),"RECREATE");
  
  const double nsamples = 18;
  std::string sampleName[nsamples] = {"FGD1_numuCC_0pi", "FGD1_numuCC_1pi", "FGD1_numuCC_other", "FGD1_anti-numuCC_0pi", "FGD1_anti-numuCC_1pi", "FGD1_anti-numuCC_other", "FGD1_NuMuBkg_CC0pi_in_AntiNu_Mode", "FGD1_NuMuBkg_CC1pi_in_AntiNu_Mode", "FGD1_NuMuBkg_CCother_in_AntiNu_Mode", "FGD2_numuCC_0pi", "FGD2_numuCC_1pi", "FGD2_numuCC_other", "FGD2_anti-numuCC_0pi", "FGD2_anti-numuCC_1pi", "FGD2_anti-numuCC_other", "FGD2_NuMuBkg_CC0pi_in_AntiNu_Mode", "FGD2_NuMuBkg_CC1pi_in_AntiNu_Mode", "FGD2_NuMuBkg_CCother_in_AntiNu_Mode"};

  std::string sampleNamespace[nsamples] = {"FGD1 numuCC 0pi", "FGD1 numuCC 1pi", "FGD1 numuCC other", "FGD1 anti-numuCC 0pi", "FGD1 anti-numuCC 1pi", "FGD1 anti-numuCC other", "FGD1 NuMuBkg CC0pi in AntiNu Mode", "FGD1 NuMuBkg CC1pi in AntiNu Mode", "FGD1 NuMuBkg CCother in AntiNu Mode", "FGD2 numuCC 0pi", "FGD2 numuCC 1pi", "FGD2 numuCC other", "FGD2 anti-numuCC 0pi", "FGD2 anti-numuCC 1pi", "FGD2 anti-numuCC other", "FGD2 NuMuBkg CC0pi in AntiNu Mode", "FGD2 NuMuBkg CC1pi in AntiNu Mode", "FGD2 NuMuBkg CCother in AntiNu Mode"};

  //loop over all samples
  for(int s = 0; s<nsamples; s++){//nsamples; s++){

    std::cout << "Doing for " << sampleNamespace[s] << std::endl;
    //get the poly
    std::string getname = sampleNamespace[s] + "/" + sampleName[s] + "_nom_mean";
    std::string outname = sampleName[s] + "_nom_mean";

    TH2Poly *poly=(TH2Poly*)(filePoly->Get(getname.c_str()));
    if(!poly)
      std::cout << "No hist " << getname << std::endl;
    poly->SetName(("p"+getname).c_str());
    poly->SetTitle(("p"+getname).c_str());

    double xlow, xup, ylow, yup;
    std::vector<double> xedges;
    std::vector<double> yedges;

    //loop over bins in the poly
    for(int i=0; i<poly->GetNumberOfBins(); i++)
      {
	//get bin and its edges
	TH2PolyBin* bin = poly->GetBins()->At(i)->Clone();
	xlow=bin->GetXMin();
	xup=bin->GetXMax();
	ylow=bin->GetYMin();
	yup=bin->GetYMax();

	//loop over vectors of edges to see if we already saved the edge value in the vector
	int xlowIn=0, xupIn=0;
	if(xedges.size()!=0){
	  for(int j=0; j<xedges.size(); j++)
	    {
	      if(abs(xedges.at(j)-xlow)<0.0001){
		xlowIn=1;
	      }
	      if(abs(xedges.at(j)-xup)<0.0001){
		xupIn=1;
	      }
	    }
	}
	//If it's not already in the vector, add it!
	if(!xlowIn){
	  xedges.push_back(xlow);
	}
	if(!xupIn){
	  xedges.push_back(xup);
	}
	
	//exactly the same for y
	int ylowIn=0, yupIn=0;
	if(yedges.size()!=0){
	  for(int j=0; j<yedges.size(); j++)
	    {
	      if(fabs(yedges.at(j)-ylow)<0.0001){
		ylowIn=1;
	      }
	      if(fabs(yedges.at(j)-yup)<0.0001){
		yupIn=1;
	      }
	    }
	}
        if(!ylowIn){
          yedges.push_back(ylow);
	}
	if(!yupIn){
          yedges.push_back(yup);
	}
      }//end loop over bins

    //Let's check if xbins*ybins = total bins
    if((xedges.size()-1)*(yedges.size()-1) != poly->GetNumberOfBins())
      {
	std::cout << "Num X bins ("<< xedges.size()-1 << ") * Num Y bins (" << yedges.size()-1 << ") does not equal number of bins in poly (" << poly->GetNumberOfBins() << "). You're input poly probably is not uniform rectangularly binned" << std::endl;
      }
    else
      {
	std::cout << "Num X bins ("<< xedges.size()-1 << ") * Num Y bins (" << yedges.size()-1 << ") equals number of bins in poly ("<< poly->GetNumberOfBins() << "). So input poly probably is uniform rectangularly binned" << std::endl;
      }
    //put vector in array of doubles for easier making of th2d (probably a better way of doing this)

    const int numxedges = xedges.size();
    const int numyedges = yedges.size();

    double* xbins = &xedges[0];
    double* ybins = &yedges[0];
    
    std::cout << "Xbins = ";
    for(int i=0; i<xedges.size(); i++){
      std::cout << xedges.at(i) << ", ";
    } 
    std::cout << std::endl;
    std::cout << "Ybins = ";
    for(int i=0; i<yedges.size(); i++){
      std::cout << yedges.at(i) << ", ";
    }
    std::cout << std::endl;

    //make the th2d
    TH2D* hist = new TH2D(outname.c_str(), outname.c_str(), xedges.size()-1, &xbins[0], yedges.size()-1, &ybins[0]);

    //Now loop over poly bins, find the corresponding th2d and setbincontent!
    for(int i=0; i<poly->GetNumberOfBins(); i++){
      TH2PolyBin* polybin = poly->GetBins()->At(i)->Clone();
      xlow=polybin->GetXMin();
      xup=polybin->GetXMax();
      ylow=polybin->GetYMin();
      yup=polybin->GetYMax();
      int xbin, ybin;
      xbin=hist->GetXaxis()->FindBin(xlow+(xup-xlow)/2);
      ybin=hist->GetYaxis()->FindBin(ylow+(yup-ylow)/2);

      std::cout << "Poly bin " << i << ", xlow: " << xlow << ", xup: " << xup << ", ylow: " << ylow << ", yup: " << yup << ". Finding bin for (" << (xlow+(xup-xlow)/2) << "," << (ylow+(yup-ylow)/2) << ")" << ". Found Bin (" << xbin << "," << ybin << ") with content " << polybin->GetContent() << ". But Poly content: " << poly->GetBinContent(i) << std::endl;
      hist->SetBinContent(xbin,ybin,polybin->GetContent());
    }


    /*For now ignoring overflow. BInned/unbinned event rates are calculated and 

    //Finally let's handle the overflow. Bins -1, -2, -3 are above Y max. Bins -7, -8, -9 are below Ymin. Bins -1, -4 -7 are below Xmin. Bins -3, -6, -9 are above Xmax. For the bins in both we'll just put half in each for want of a better solution.
    double m[9];
    double xunderflow=0, xoverflow=0, yunderflow=0, yoverflow=0;
    for (int i=-9;i<0;i++)
      {
	TH2PolyBin* polybin = poly->GetBins()->At(i)->Clone();
	m[-1*i]=polybin->GetContent();
      }
    xunderflow+=m[1*/


    std::cout << "Made the TH2D for " << sampleName[s] << "!" << std::endl;
    //save the th2ds
    fileTH2D->cd();
    hist->Write();
    //hist->Write();
  }
}//end convertPolyToTH2D

