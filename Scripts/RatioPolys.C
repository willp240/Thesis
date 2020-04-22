#include <TH2D.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <string.h>
#include <TH2Poly.h>

void RatioPolys (std::string file1name) {
 
  TCanvas c1("c1");
  c1.SetRightMargin(0.18);
  TString name = file1name "_RatioPlots";

  TFile *file1 = TFile::Open(file1name.c_str());
   c1.Print(name + ".pdf[");
  TFile* outfile = new TFile(name + ".root", "RECREATE");

  int numSamples = 36;
  TString samples[36] = {"DATA_FGD1_numuCC_0pi",
			 "DATA_FGD1_numuCC_1pi",
			 "DATA_FGD1_numuCC_other",
			 "DATA_FGD2_numuCC_0pi",
			 "DATA_FGD2_numuCC_1pi",
			 "DATA_FGD2_numuCC_other",
			 "DATA_FGD1_anti-numuCC_QE",
			 "DATA_FGD1_anti-numuCC_nQE",
			 "DATA_FGD2_anti-numuCC_1_track",
			 "DATA_FGD2_anti-numuCC_N_tracks",
			 "MC_FGD1_numuCC_0pi",
			 "MC_FGD1_numuCC_1pi",
			 "MC_FGD1_numuCC_other",
			 "MC_FGD2_numuCC_0pi",
			 "MC_FGD2_numuCC_1pi",
			 "MC_FGD2_numuCC_other",
			 "MC_FGD1_anti-numuCC_QE",
			 "MC_FGD1_anti-numuCC_nQE",
			 "MC_FGD2_anti-numuCC_1_track",
			 "MC_FGD2_anti-numuCC_N_tracks"};

  for (int k = 0; k < numSamples; k++)
    {
      
      std::cout << "Getting " << samples[k] << std::endl;
      
      TH2D *h1 = (TH2D*)file1->Get(samples[k]);
      if(!h1)
	std::cout << "No h1!" << std::endl;
      TString h1name = samples[k]+"_1";
      h1->SetName(h1name);
      h1->SetTitle(h1name);

      TH2D *h2 = (TH2D*)file2->Get(samples[k]);
      if(!h2)
	std::cout << "No h2!" << std::endl;
      TString h2name = samples[k]+"_2";
      h2->SetName(h2name);
      h2->SetTitle(h2name);

      TH2D *hRatio = (TH2D*)file1->Get(samples[k]);
      TString hRationame = samples[k]+"_Ratio";
      hRatio->SetName(hRationame);
      hRatio->SetTitle(hRationame);

      hRatio->Divide(h2);
      TString zname = "Ratio";
      TString title = samples[k] + " Ratio"; 
      hRatio->GetZaxis()->SetTitle(zname);
      hRatio->GetXaxis()->SetRangeUser(0,5000);
      hRatio->GetZaxis()->SetRangeUser(0.7,1.3);
      gPad->Update();
      
      if((k>5&&k<10)||k>15){
	hRatio->GetXaxis()->SetRangeUser(0,3000);
      }

      h1->Draw("colz");
      gPad->Update();
      //c1.SetRightMargin(1.3);
      palette = (TPaletteAxis*) h1->GetListOfFunctions()->FindObject("palette");
      palette->SetX1NDC(0.83);
      palette->SetX2NDC(0.88);
      h1->GetXaxis()->SetRangeUser(0,5000);
      h1->GetZaxis()->SetLabelSize(0.03);
      c1.Print(name + ".pdf");
      h1->Write(h1name);

      h2->Draw("colz");
      gPad->Update();
      //c1.SetRightMargin(1.3);
      palette = (TPaletteAxis*) h2->GetListOfFunctions()->FindObject("palette");
      palette->SetX1NDC(0.83);
      palette->SetX2NDC(0.88);
      h2->GetXaxis()->SetRangeUser(0,5000);
      h2->GetZaxis()->SetLabelSize(0.03);
      c1.Print(name + ".pdf");
      h2->Write(h2name);

      hRatio->Draw("colz");
      gPad->Update();
      //c1.SetRightMargin(1.3);
      palette = (TPaletteAxis*) hRatio->GetListOfFunctions()->FindObject("palette");
      palette->SetX1NDC(0.83);
      palette->SetX2NDC(0.88);
      hRatio->GetXaxis()->SetLabelSize(0.03);
      hRatio->GetZaxis()->SetLabelSize(0.03);
      c1.Print(name + ".pdf");
      hRatio->Write(hRationame);
     
      delete h1;
      delete h2;
      delete hRatio;

    }

  c1.Print(name + ".pdf]");
 
}
