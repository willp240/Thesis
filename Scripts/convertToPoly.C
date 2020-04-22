#include <TH2Poly.h>
#include <TH2D.h>
#include <vector>
#include <TMath.h>
#include <TFile.h>
#include <TTree.h>
#include <THStack.h>
#include <iostream> 

void convertToPoly (){
   
  //  TFile* fileTH2 = TFile::Open("/scratch4/wparker2/files/PreFarmegeddon/144403391079/12DecRHCAsimovRuns2to6/12DecRHCAsimovRuns2to6_0_ND280_nom.root");
  TFile* fileTH2 = TFile::Open("JacFitBinDetectorTH2DBinning.root");
  TFile * filePoly = new TFile("convertedtoPolyJacFitBinDetector.root","RECREATE");
  const double nsamples = 18;
  TString sampleName[nsamples] = {"FGD1_numuCC_0pi", "FGD1_numuCC_1pi", "FGD1_numuCC_other", "FGD1_anti-numuCC_0pi", "FGD1_anti-numuCC_1pi", "FGD1_anti-numuCC_other", "FGD1_NuMuBkg_CC0pi_in_AntiNu_Mode", "FGD1_NuMuBkg_CC1pi_in_AntiNu_Mode", "FGD1_NuMuBkg_CCother_in_AntiNu_Mode", "FGD2_numuCC_0pi", "FGD2_numuCC_1pi", "FGD2_numuCC_other", "FGD2_anti-numuCC_0pi", "FGD2_anti-numuCC_1pi", "FGD2_anti-numuCC_other", "FGD2_NuMuBkg_CC0pi_in_AntiNu_Mode", "FGD2_NuMuBkg_CC1pi_in_AntiNu_Mode", "FGD2_NuMuBkg_CCother_in_AntiNu_Mode"};

  for(int s = 0; s<nsamples; s++){

    TString histname = sampleName[s];
    TString getname = "MC_" + sampleName[s];
    TH2D *hist=(TH2D*)(fileTH2->Get(getname));
    if(!hist)
      std::cout << "No hist " << sampleName[s] << std::endl;
    
    hist->SetName("h_"+histname);
    hist->SetTitle("h_"+histname);
    TH2Poly *poly = new TH2Poly();
    poly->SetName(histname);
    poly->SetTitle(histname);
    double xmax, xmin, ymax, ymin;
    
    for(int i=1; i<=hist->GetXaxis()->GetNbins(); i++)
      {
	xmax = hist->GetXaxis()->GetBinUpEdge(i);
	xmin = hist->GetXaxis()->GetBinLowEdge(i);
	for(int j=1; j<=hist->GetYaxis()->GetNbins(); j++)
	  {
	    ymax = hist->GetYaxis()->GetBinUpEdge(j);
	    ymin = hist->GetYaxis()->GetBinLowEdge(j);
	    double binofx[] = {xmin, xmax, xmax, xmin};
	    double binofy[] = {ymin, ymin, ymax, ymax};
	  
	    poly->AddBin(4,binofx,binofy);
	  }
      }
    /*    int globBin =0;
    for(int i=1; i<=hist->GetXaxis()->GetNbins(); i++)
      {
        for(int j=1; j<=hist->GetYaxis()->GetNbins(); j++)
          {
	    globBin++;
	    poly->SetBinContent(globBin,hist->GetBinContent(i,j));
          }
	  }*/
    

    filePoly->cd();
    poly->Write();
    //hist->Write();
  }


}//end converToPoly

