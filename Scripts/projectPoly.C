#include <TH2Poly.h>
#include <TH2D.h>
#include <vector>
#include <TMath.h>
#include <TFile.h>
#include <TTree.h>
#include <THStack.h>
#include <iostream> 

/*
Script takes in TH2Polys, projects them onto x and y. Need a user defined binning for each axis for each sample. This is read in from a file of th2ds. I.e take the banff th2d binning and project onto these. If a poly bin overlaps with a th2d bin, we assume the events are uniformly distributed within a bin. i.e. if a poly bin is 50% in one th2d bin, and 50% in another, we put 50 % of events from the poly bin in each th2d bin.
 */

void projectPoly (std::string PolyFileName){
   
  std::string OutFileName = PolyFileName;
  OutFileName.replace(OutFileName.find(".root"), 5, std::string("_Projected.root"));
  TFile* filePoly = TFile::Open(PolyFileName.c_str());
  TFile * fileOut = new TFile(OutFileName.c_str(),"RECREATE");
  
  const double nsamples = 18;
  std::string sampleName[nsamples] = {"FGD1_numuCC_0pi", "FGD1_numuCC_1pi", "FGD1_numuCC_other", "FGD1_anti-numuCC_0pi", "FGD1_anti-numuCC_1pi", "FGD1_anti-numuCC_other", "FGD1_NuMuBkg_CC0pi_in_AntiNu_Mode", "FGD1_NuMuBkg_CC1pi_in_AntiNu_Mode", "FGD1_NuMuBkg_CCother_in_AntiNu_Mode", "FGD2_numuCC_0pi", "FGD2_numuCC_1pi", "FGD2_numuCC_other", "FGD2_anti-numuCC_0pi", "FGD2_anti-numuCC_1pi", "FGD2_anti-numuCC_other", "FGD2_NuMuBkg_CC0pi_in_AntiNu_Mode", "FGD2_NuMuBkg_CC1pi_in_AntiNu_Mode", "FGD2_NuMuBkg_CCother_in_AntiNu_Mode"};

  double numMomBins[nsamples]= {29,18,18,29,18,18,17,6,8,17,6,8,10,4,6,10,4,6};
  double momBins[nsamples][30]=  {
    {0.   , 200. , 300. , 400. , 450. , 500. , 550. , 600. ,
     650. , 700. , 750. , 800. , 850. , 900. , 950. , 1000.,
     1050., 1100., 1200., 1300., 1400., 1500., 1600., 1700.,
     1800., 2000., 2500., 3000., 5000., 30000.}, //FGD1 CC0Pi
    {0.   , 300. , 350. , 400. , 500. , 600. , 650. , 700. ,
     750. , 800. , 900. , 1000., 1100., 1200., 1500., 2000.,
     3000., 5000., 30000.},//FGD1 CC1Pi
    {0.   , 300. , 400. , 500. , 600. , 650. , 700. , 750. ,
     800. , 900. , 1000., 1100., 1250., 1500., 1750., 2000.,
     3000., 5000., 30000.},//FGD1 CCOther
    {0.   , 200. , 300. , 400. , 450. , 500. , 550. , 600. ,
     650. , 700. , 750. , 800. , 850. , 900. , 950. , 1000.,
     1050., 1100., 1200., 1300., 1400., 1500., 1600., 1700.,
     1800., 2000., 2500., 3000., 5000., 30000.}, //FGD2 CC0Pi
    {0.   , 300. , 350. , 400. , 500. , 600. , 650. , 700. ,
     750. , 800. , 900. , 1000., 1100., 1200., 1500., 2000.,
     3000., 5000., 30000.},//FGD2 CC1Pi
    {0.   , 300. , 400. , 500. , 600. , 650. , 700. , 750. ,
     800. , 900. , 1000., 1100., 1250., 1500., 1750., 2000.,
     3000., 5000., 30000.},//FGD2 CCOther
    {0.   , 300. , 400. , 500. , 550. , 600. , 650. , 700. ,
     750. , 800. , 900. , 1000., 1100., 1200., 1500., 2000.,
     4000., 30000.}, //FGD1 nubar CC0Pi
    {0.   , 500. , 700. , 900. , 1300., 2500., 30000.}, //FGD1 nubar CC1Pi
    {0. , 600. , 800. , 1000., 1250., 1500., 2000., 4000., 30000.}, //FGD1 nubar CCOther
    {0.   , 300. , 400. , 500. , 550. , 600. , 650. , 700. ,
     750. , 800. , 900. , 1000., 1100., 1200., 1500., 2000.,
     4000., 30000.}, //FGD2 nubar CC0Pi
    {0.   , 500. , 700. , 900. , 1300., 2500., 30000.}, //FGD2 nubar CC1Pi
    {0. , 600. , 800. , 1000., 1250., 1500., 2000., 4000., 30000.}, //FGD2 nubar CCOther
    {0. , 300. , 500. , 700. , 800. , 900. , 1250., 1500., 2000., 4000., 30000.}, //FGD1 bkg CC0Pi
    {0. , 600. , 800. , 1500., 30000.}, //FGD1 bkg CC1Pi
    {0. , 600. , 1000., 1250., 2000., 4000., 30000.}, //FGD1 bkg CCOther
    {0. , 300. , 500. , 700. , 800. , 900. , 1250., 1500., 2000., 4000., 30000.}, //FGD2 bkg CC0Pi
    {0. , 600. , 800. , 1500., 30000.}, //FGD2 bkg CC1Pi
    {0. , 600. , 1000., 1250., 2000., 4000., 30000.}, //FGD2 bkg CCOther
  };

  double numThetaBins[nsamples]= {29,16,19,29,16,19,18,8,10,18,8,10,12,10,9,12,10,9};
  double thetaBins[nsamples][30]= {
    {-1.0,0.5,0.6,0.7,0.76,0.78,0.8,0.83,0.85,0.88,0.89,0.9,0.91,0.92,0.925,0.93,0.935,0.94,0.945,0.95,0.955,0.96,0.965,0.97,0.975,0.98,0.985,0.99,0.995,1.0}, //FGD1 CC0Pi
    {-1.0,0.6,0.7,0.8,0.85,0.88,0.9,0.92,0.93,0.94,0.95,0.96,0.97,0.98,0.99,0.995,1.0},//FGD1 CC1Pi
    {-1.0,0.6,0.7,0.76,0.8,0.85,0.88,0.89,0.9,0.91,0.92,0.93,0.94,0.95,0.96,0.97,0.98,0.99,0.995,1.0},//FGD1 CCOther
    {-1.0,0.5,0.6,0.7,0.76,0.78,0.8,0.83,0.85,0.88,0.89,0.9,0.91,0.92,0.925,0.93,0.935,0.94,0.945,0.95,0.955,0.96,0.965,0.97,0.975,0.98,0.985,0.99,0.995,1.0}, //FGD2 CC0Pi
    {-1.0,0.6,0.7,0.8,0.85,0.88,0.9,0.92,0.93,0.94,0.95,0.96,0.97,0.98,0.99,0.995,1.0},//FGD2 CC1Pi
    {-1.0,0.6,0.7,0.76,0.8,0.85,0.88,0.89,0.9,0.91,0.92,0.93,0.94,0.95,0.96,0.97,0.98,0.99,0.995,1.0},//FGD2 CCOther
    {-1.0,0.6,0.7,0.8,0.85,0.9,0.92,0.93,0.94,0.95,0.96,0.965,0.97,0.975,0.98,0.985,0.99,0.995,1.0}, //FGD2 nubar CC0Pi
    {-1.0,0.7,0.8,0.9,0.94,0.96,0.98,0.99,1.0}, //FGD2 nubar CC1Pi
    {-1.0,0.7,0.8,0.85,0.9,0.93,0.95,0.97,0.98,0.99,1.0}, //FGD2 nubar CCOther
    {-1.0,0.6,0.7,0.8,0.85,0.9,0.92,0.93,0.94,0.95,0.96,0.965,0.97,0.975,0.98,0.985,0.99,0.995,1.0}, //FGD2 nubar CC0Pi
    {-1.0,0.7,0.8,0.9,0.94,0.96,0.98,0.99,1.0}, //FGD2 nubar CC1Pi
    {-1.0,0.7,0.8,0.85,0.9,0.93,0.95,0.97,0.98,0.99,1.0}, //FGD2 nubar CCOther
    {-1.0,0.7,0.8,0.85,0.88,0.9,0.92,0.94,0.96,0.97,0.98,0.99,1.0}, //FGD1 bkg CC0Pi
    {-1.0,0.7,0.8,0.86,0.9,0.94,0.96,0.97,0.98,0.99,1.0}, //FGD1 bkg CC1Pi
    {-1.0,0.7,0.8,0.86,0.9,0.93,0.95,0.97,0.99,1.0}, //FGD1 bkg CCOther
    {-1.0,0.7,0.8,0.85,0.88,0.9,0.92,0.94,0.96,0.97,0.98,0.99,1.0}, //FGD2 bkg CC0Pi
    {-1.0,0.7,0.8,0.86,0.9,0.94,0.96,0.97,0.98,0.99,1.0}, //FGD2 bkg CC1Pi
    {-1.0,0.7,0.8,0.86,0.9,0.93,0.95,0.97,0.99,1.0}, //FGD2 bkg CCOther
  };

  //loop over all samples
  for(int s = 0; s<1; s++){//nsamples; s++){

    std::cout << "Doing for " << sampleName[s] << std::endl;
    //get the poly
    std::string getname = sampleName[s];//"MC_" + sampleName[s];

    TH2Poly *poly=(TH2Poly*)(filePoly->Get(getname.c_str()));
    if(!poly)
      std::cout << "No hist " << getname << std::endl;
    poly->SetName(("p"+getname).c_str());
    poly->SetTitle(("p"+getname).c_str());

    TH1D* hProjX = new TH1D((getname+"_x").c_str(),(getname+"_x").c_str(),numMomBins[s],&momBins[s][0]);
    TH1D* hProjY = new TH1D((getname+"_y").c_str(),(getname+"_y").c_str(),numThetaBins[s],&thetaBins[s][0]);

    double xlow, xup, ylow, yup, binlow, binup, frac=0;

    //loop over bins in the poly
    for(int i=0; i<poly->GetNumberOfBins(); i++)
      {
	//get bin and its edges
	TH2PolyBin* bin = poly->GetBins()->At(i)->Clone();
	xlow=bin->GetXMin();
	xup=bin->GetXMax();
	ylow=bin->GetYMin();
	yup=bin->GetYMax();

	std::cout << "Bin " << xlow << "-" << xup << ", " << ylow << "-" << yup << std::endl;  
	int acase =0;
	for(int dx=0; dx<numMomBins[s]; dx++)
	  {	 
	    if(momBins[s][dx+1]<=xlow || momBins[s][dx]>=xup)
	      {
		frac=0;
		acase=0;
	      }
	    else if(momBins[s][dx]<=xlow && momBins[s][dx+1]>=xup)
	      {
		frac=1;
		acase=1;
	      }
	    else if(momBins[s][dx]<=xlow && momBins[s][dx+1]<=xup)
	      {
		frac=(momBins[s][dx+1]-xlow)/(xup-xlow);
		acase=2;
	      }
	    else if(momBins[s][dx]>=xlow && momBins[s][dx+1]>=xup)
	      {
                frac=(xup-momBins[s][dx])/(xup-xlow);
		acase=3;
              }
	    else if(momBins[s][dx]>=xlow && momBins[s][dx+1]<=xup)
              {
                frac=(momBins[s][dx+1]-momBins[s][dx])/(xup-xlow);
		acase=4;
              }
	    else
	      {
		frac=0;
		acase=5;
	      }
	    std::cout << "        1D Bin:  " << momBins[s][dx] << "-" << momBins[s][dx+1] << " so Frac=" << frac << ", case = " << acase << std::endl; 
	    hProjX->SetBinContent(dx+1,hProjX->GetBinContent(dx+1)+frac*bin->GetContent());
	  }

        for(int dy=0; dy<numThetaBins[s]; dy++)
          {
            if(thetaBins[s][dy+1]<=ylow || thetaBins[s][dy]>=yup)
              {
                frac=0;
                acase=0;
              }
            else if(thetaBins[s][dy]<=ylow && thetaBins[s][dy+1]>=yup)
              {
                frac=1;
                acase=1;
              }
            else if(thetaBins[s][dy]<=ylow && thetaBins[s][dy+1]<=yup)
              {
                frac=(thetaBins[s][dy+1]-ylow)/(yup-ylow);
                acase=2;
              }
            else if(thetaBins[s][dy]>=ylow && thetaBins[s][dy+1]>=yup)
              {
                frac=(yup-thetaBins[s][dy])/(yup-ylow);
                acase=3;
              }
            else if(thetaBins[s][dy]>=ylow && thetaBins[s][dy+1]<=yup)
              {
                frac=(thetaBins[s][dy+1]-thetaBins[s][dy])/(yup-ylow);
                acase=4;
              }
            else
              {
                frac=0;
                acase=5;
              }
	    std::cout << "        1D Bin:  " << thetaBins[s][dy] << "-" << thetaBins[s][dy+1] << " so Frac=" << frac << ", case = " << acase << std::endl;
            hProjY->SetBinContent(dy+1,hProjY->GetBinContent(dy+1)+frac*bin->GetContent());
          }

      }//end loop over bins

    poly->Write();
    hProjX->Write();
    hProjY->Write();
    //save the th1ds
    fileOut->cd();
  }
}//end projectPoly

