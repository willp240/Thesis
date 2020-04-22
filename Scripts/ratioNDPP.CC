#include <TH2Poly.h>
#include <TH2D.h>
#include <vector>
#include <TMath.h>
#include <TFile.h>
#include <TTree.h>
#include <THStack.h>
#include <iostream> 

void ratioNDPP(std::string file1, std::string file2){
   
  TFile* f1 = TFile::Open(file1.c_str());
  TFile* f2 = TFile::Open(file2.c_str());
  std::string outname = file1+file2+"_ratio";
  TFile * fileOut = new TFile((outname+".root").c_str(),"RECREATE");
  
  const double nsamples = 18;
  std::string sampleName[nsamples] = {"FGD1 numuCC 0pi", "FGD1 numuCC 1pi", "FGD1 numuCC other", "FGD1 anti-numuCC 0pi", "FGD1 anti-numuCC 1pi", "FGD1 anti-numuCC other", "FGD1 NuMuBkg CC0pi in AntiNu Mode", "FGD1 NuMuBkg CC1pi in AntiNu Mode", "FGD1 NuMuBkg CCother in AntiNu Mode", "FGD2 numuCC 0pi", "FGD2 numuCC 1pi", "FGD2 numuCC other", "FGD2 anti-numuCC 0pi", "FGD2 anti-numuCC 1pi", "FGD2 anti-numuCC other", "FGD2 NuMuBkg CC0pi in AntiNu Mode", "FGD2 NuMuBkg CC1pi in AntiNu Mode", "FGD2 NuMuBkg CCother in AntiNu Mode"};

  std::string polyName[nsamples] = {"FGD1_numuCC_0pi", "FGD1_numuCC_1pi", "FGD1_numuCC_other", "FGD1_anti-numuCC_0pi", "FGD1_anti-numuCC_1pi", "FGD1_anti-numuCC_other", "FGD1_NuMuBkg_CC0pi_in_AntiNu_Mode", "FGD1_NuMuBkg_CC1pi_in_AntiNu_Mode", "FGD1_NuMuBkg_CCother_in_AntiNu_Mode", "FGD2_numuCC_0pi", "FGD2_numuCC_1pi", "FGD2_numuCC_other", "FGD2_anti-numuCC_0pi", "FGD2_anti-numuCC_1pi", "FGD2_anti-numuCC_other", "FGD2_NuMuBkg_CC0pi_in_AntiNu_Mode", "FGD2_NuMuBkg_CC1pi_in_AntiNu_Mode", "FGD2_NuMuBkg_CCother_in_AntiNu_Mode"};

  TCanvas *c = new TCanvas("canv", "canv", 1080, 1080);
  c->SetTopMargin(0.07);
  c->SetBottomMargin(0.11);
  c->SetLeftMargin(0.10);
  c->SetRightMargin(0.18);

  c->Print((outname+".pdf[").c_str());

  //loop over all samples
  for(int s = 0; s<nsamples; s++){//nsamples; s++){

    std::cout << "Doing for " << sampleName[s] << "/" << polyName[s] << "_nom_mean" << std::endl;
    //get the poly
    std::string getname = sampleName[s] + "/" + polyName[s] + "_nom_mean";
    std::string p1name = polyName[s] + "_p1";
    std::string p2name = polyName[s] + "_p2";
    TH2Poly *poly1=(TH2Poly*)(f1->Get(getname.c_str()))->Clone();
    if(!poly1)
      std::cout << "No hist " << getname << " in f1" << std::endl;
    poly1->SetName((p1name).c_str());
    poly1->SetTitle((p1name).c_str());
    TH2Poly *poly2=(TH2Poly*)(f2->Get(getname.c_str()))->Clone();
    if(!poly2)
      std::cout << "No hist " << getname << " in f1" << std::endl;
    poly2->SetName((p2name).c_str());
    poly2->SetTitle((p2name).c_str());

    TH2Poly *ratio= (TH2Poly*)(f1->Get(getname.c_str()))->Clone();
    std::string rationame = polyName[s] + "_mean_posterior_prediction";
    ratio->Reset("");
    
    double xlow1, xup1, ylow1, yup1, xlow2, xup2, ylow2, yup2;
    
    for(int i=1; i<=poly1->GetNumberOfBins(); i++){
      TH2PolyBin* polybin1 = poly1->GetBins()->At(i-1)->Clone();
      xlow1=polybin1->GetXMin();
      xup1=polybin1->GetXMax();
      ylow1=polybin1->GetYMin();
      yup1=polybin1->GetYMax();
      int binR, bin2;
      bin2=poly2->FindBin(xlow1+(xup1-xlow1)/2, ylow1+(yup1-ylow1)/2);
      binR=ratio->FindBin(xlow1+(xup1-xlow1)/2, ylow1+(yup1-ylow1)/2);
      ratio->SetBinContent(binR,poly1->GetBinContent(i)/poly2->GetBinContent(bin2));
    }
    ratio->SetName((rationame).c_str());
    ratio->SetTitle((rationame).c_str());
    ratio->GetZaxis()->SetRangeUser(0.85,1.15);
    ratio->GetXaxis()->SetRangeUser(0,5000);
    ratio->GetZaxis()->SetTitle("Ratio");
    poly1->GetXaxis()->SetRangeUser(0,5000);
    poly2->GetXaxis()->SetRangeUser(0,5000);    


    fileOut->cd();
    poly1->Draw("colz");
    c->Print((outname+".pdf").c_str());

    poly2->Draw("colz");
    c->Print((outname+".pdf").c_str());

    ratio->Draw("colz");
    c->Print((outname+".pdf").c_str());

  }//end loop over samples
  c->Print((outname+".pdf]").c_str());

}//end ratioNDPP

