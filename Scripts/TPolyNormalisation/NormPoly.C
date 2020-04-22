
void NormPoly(){
  TFile *outfile = new TFile("TPolyNorm_5.root", "RECREATE");
  TFile *file1 = TFile::Open("5JunTPolyAsimov_0_ND280_nom.root");

  TCanvas *c = new TCanvas("canv", "canv", 1080, 1080);
  c->SetTopMargin(0.07);
  c->SetBottomMargin(0.11);
  c->SetLeftMargin(0.11);
  c->SetRightMargin(0.18);

  c->Print((std::string("TPolyNorm_5")+std::string(".pdf[")).c_str());

  TString sample[18] = {"MC_FGD1_numuCC_0pi",
			"MC_FGD1_numuCC_1pi",
			"MC_FGD1_numuCC_other",
			"MC_FGD2_numuCC_0pi",
			"MC_FGD2_numuCC_1pi",
			"MC_FGD2_numuCC_other",
			"MC_FGD1_anti-numuCC_0pi",
			"MC_FGD1_anti-numuCC_1pi",
			"MC_FGD1_anti-numuCC_other",
			"MC_FGD2_anti-numuCC_0pi",
			"MC_FGD2_anti-numuCC_1pi",
			"MC_FGD2_anti-numuCC_other",
			"MC_FGD1_NuMuBkg_CC0pi_in_AntiNu_Mode",
			"MC_FGD1_NuMuBkg_CC1pi_in_AntiNu_Mode",
			"MC_FGD1_NuMuBkg_CCother_in_AntiNu_Mode",
			"MC_FGD2_NuMuBkg_CC0pi_in_AntiNu_Mode",
			"MC_FGD2_NuMuBkg_CC1pi_in_AntiNu_Mode",
			"MC_FGD2_NuMuBkg_CCother_in_AntiNu_Mode"};
  
  outfile->cd();

  for(int sam =0; sam < 18; sam++){
    TH2Poly *poly = (TH2Poly*)(file1->Get(sample[sam]))->Clone();
    TH2Poly *polyNorm = (TH2Poly*)(file1->Get(sample[sam]))->Clone();

    for(int i = 1; i< poly->GetNumberOfBins(); i++)
      {
	TH2PolyBin *polybin;
	polybin = (TH2PolyBin*)poly->GetBins()->At(i)->Clone();
	double norm = polybin->GetContent()/polybin->GetArea();
	polyNorm->SetBinContent(i,norm);
    }
    
    polyNorm->GetZaxis()->SetTitle("Events MeV^{-1} cos #theta^{-1}");
    polyNorm->GetXaxis()->SetRangeUser(0,5000);
    polyNorm->Draw("colz");
    c->Print((std::string("TPolyNorm_5.pdf")).c_str());
    polyNorm->Write();
  }
  c->Print((std::string("TPolyNorm_5")+std::string(".pdf]")).c_str());
}
