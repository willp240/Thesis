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
#include "TH2Poly.h"
#include <vector>
#include "../samplePDF/Structs.h"

// *******************
// Template to make vector out of an array of any length
//template< typename T, size_t N >
//std::vector<T> MakeVector( const T (&data)[N] ) {
  // *******************
  //return std::vector<T>(data, data+N);
//}

TH1D* PolyProjectionX(TObject* poly, std::string TempName, std::vector<double> xbins);
TH1D* PolyProjectionY(TObject* poly, std::string TempName, std::vector<double> ybins);

void SetBinning(int Selection, std::vector<double> & BinningX,std::vector<double> & BinningY);

int main(int argc, char **argv) {

    if (argc != 2) {
      std::cout << "thstack inputFile" << std::endl;
      exit(-1);
    }

    TFile *outfile = new TFile("thstacks.root", "RECREATE");

    TCanvas *c = new TCanvas("canv", "canv", 1580, 1080);
    TFile *file_ = TFile::Open(argv[1]);
    
    c->SetTopMargin(0.05);
    c->SetBottomMargin(0.14);
    c->SetLeftMargin(0.1);
    c->SetRightMargin(0.16);
    
    std::string samp[18] = {"FGD1_numuCC_0pi", "FGD1_numuCC_1pi", "FGD1_numuCC_other", "FGD1_anti-numuCC_0pi", "FGD1_anti-numuCC_1pi", "FGD1_anti-numuCC_other", "FGD1_NuMuBkg_CC0pi_in_AntiNu_Mode", "FGD1_NuMuBkg_CC1pi_in_AntiNu_Mode", "FGD1_NuMuBkg_CCother_in_AntiNu_Mode", "FGD2_numuCC_0pi", "FGD2_numuCC_1pi", "FGD2_numuCC_other", "FGD2_anti-numuCC_0pi","FGD2_anti-numuCC_1pi", "FGD2_anti-numuCC_other", "FGD2_NuMuBkg_CC0pi_in_AntiNu_Mode", "FGD2_NuMuBkg_CC1pi_in_AntiNu_Mode", "FGD2_NuMuBkg_CCother_in_AntiNu_Mode"};  
    
    std::string mode[13] = {"", "CCQE", "CC1pi", "CCcoh", "CCMpi", "CCDIS", "NC1pi0", "NC1pipm", "NCcoh", "NCoth", "2p2h", "NC1gam", "CCMisc"};

    for(int s=0; s<18; s++){
      
      // Stacks for the projected samples
      THStack *MomStack = new THStack((samp[s]+"_x").c_str(), (samp[s]+"_x").c_str());
      THStack *ThStack = new THStack((samp[s]+"_y").c_str(), (samp[s]+"_y").c_str());
      
      for(int m=0; m<13; m++){
	
	if(m==0){
	  
	  std::string datname = "DATA_" + samp[s];
	  std::cout << "getting " <<datname << std::endl;
	  
	  std::vector<double> datxbins;
	  std::vector<double> datybins;
	  
	  SetBinning(s, datxbins, datybins);
	  
	  TH2Poly *dathist = (TH2Poly*) file_->Get(datname.c_str())->Clone();
	  TH1D* dataprojX = PolyProjectionX(dathist, datname.c_str(), datxbins);
	  TH1D* dataprojY = PolyProjectionY(dathist, datname.c_str(), datybins);

	  dataprojX->Scale(1., "width");
	  dataprojY->Scale(1., "width");

	  MomStack->Add(dataprojX);
	  ThStack->Add(dataprojY);

	  outfile->cd();
	  dataprojX->Write();
	  dataprojY->Write();
	  
	}
	
	std::string name = "MC_" + samp[s];
	if(m>0)
	  name = name + "_" + mode[m];
	std::cout << "getting " << name << std::endl;
	
	std::vector<double> xbins;
	std::vector<double> ybins;
	
	SetBinning(s, xbins, ybins);
	
	TH2Poly *hist = (TH2Poly*) file_->Get(name.c_str())->Clone();
	TH1D* projX = PolyProjectionX(hist, name.c_str(), xbins);
	TH1D* projY = PolyProjectionY(hist, name.c_str(), ybins);
	
	projX->Scale(1., "width");
	projY->Scale(1., "width");
	
	projX->SetFillColor(m+1);
	projX->SetLineStyle(0);
	projY->SetFillColor(m+1);
	projY->SetLineStyle(0);
	projX->SetFillStyle(3001);
	projY->SetFillStyle(3001);
	
	MomStack->Add(projX);
	ThStack->Add(projY);

	outfile->cd();
	projX->Write();
	projY->Write();
	
      }

      outfile->cd();

      MomStack->Draw();
      MomStack->GetXaxis()->SetTitle("p_{#mu} (MeV)");
      MomStack->GetYaxis()->SetTitle("Events/(bin width)");
      c->Print((samp[s]+std::string("_p.pdf")).c_str());
      
      ThStack->Draw();
      ThStack->GetXaxis()->SetTitle("cos #theta_{#mu}");
      ThStack->GetYaxis()->SetTitle("Events/(bin width)");
      c->Print((samp[s]+std::string("_t.pdf")).c_str());

      MomStack->Write();
      ThStack->Write();
      
      std::cout << "done " << samp[s] << std::endl;
    }
}

//Basically copy SetupBinning from samplePDF. Used for the 1d axis in projections (have to define your axis if you want to project a poly). We have functions to do it (assuming uniform distribution of events in a bin when they overlap), but need to hard code the axis
void SetBinning(int Selection, std::vector<double> & BinningX,std::vector<double> & BinningY) {

  switch (Selection) {

  case (0):
  case (9):
    {
      const double BinArray_x[] = {0., 200., 300., 400., 450., 500., 550., 600., 650., 700., 750., 800., 850., 900., 950., 1000., 1050., 1100., 1200., 1300., 1400., 1500., 1600., 1700., 1800., 2000., 2500., 3000., 5000., 30000.};
      const double BinArray_y[] = {-1., 0.5, 0.6, 0.7, 0.76, 0.78, 0.8, 0.83, 0.85, 0.88, 0.89, 0.9, 0.91, 0.92, 0.925, 0.93, 0.935, 0.94, 0.945, 0.95, 0.955, 0.96, 0.965, 0.97, 0.975, 0.98, 0.985, 0.99, 0.995, 1.};
      BinningX = MakeVector(BinArray_x);
      BinningY = MakeVector(BinArray_y);
      break;
    }

  case (1):
  case (10):
    {
      const double BinArray_x[] = {0., 300., 350., 400., 500., 600., 650., 700., 750., 800., 900., 1000., 1100., 1200., 1500., 2000., 3000., 5000., 30000.};
      const double BinArray_y[] = {-1., 0.6, 0.7, 0.8, 0.85, 0.88, 0.9, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 0.995, 1.};
      BinningX = MakeVector(BinArray_x);
      BinningY = MakeVector(BinArray_y);
      break;
    }

  case (2):
  case (11):
    {
      const double BinArray_x[] = {0., 300., 400., 500., 600., 650., 700., 750., 800., 900., 1000., 1100., 1250., 1500., 1750., 2000., 3000., 5000., 30000.};
      const double BinArray_y[] = {-1., 0.6, 0.7, 0.76, 0.8, 0.85, 0.88, 0.89, 0.9, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 0.995, 1.};
      BinningX = MakeVector(BinArray_x);
      BinningY = MakeVector(BinArray_y);
      break;
    }

  case (3):
  case (12):
    {
      const double BinArray_x[] = {0., 300., 400., 500., 550., 600., 650., 700., 750., 800., 900., 1000., 1100., 1200., 1500., 2000., 4000., 30000.};
      const double BinArray_y[] = {-1., 0.6, 0.7, 0.8, 0.85, 0.9, 0.92, 0.93, 0.94, 0.95, 0.96, 0.965, 0.97, 0.975, 0.98, 0.985, 0.99, 0.995, 1.};
      BinningX = MakeVector(BinArray_x);
      BinningY = MakeVector(BinArray_y);
      break;
    }

  case (4):
  case (13):
    {
      const double BinArray_x[] = {0., 500., 700., 900., 1300., 2500., 30000.};
      const double BinArray_y[] = {-1, 0.7, 0.8, 0.9, 0.94, 0.96, 0.98, 0.99, 1};
      BinningX = MakeVector(BinArray_x);
      BinningY = MakeVector(BinArray_y);
      break;
    }

  case (5):
  case (14):
    {
      const double BinArray_x[] = {0., 600., 800., 1000., 1250., 1500., 2000., 4000., 30000.};
      const double BinArray_y[] = {-1., 0.7, 0.8, 0.85, 0.9, 0.93, 0.95, 0.97, 0.98, 0.99, 1.};
      BinningX = MakeVector(BinArray_x);
      BinningY = MakeVector(BinArray_y);
      break;
    }

  case (6):
  case (15):
    {
      const double BinArray_x[] = {0., 300., 500., 700., 800., 900., 1250., 1500., 2000., 4000., 30000.};
      const double BinArray_y[] = {-1., 0.7, 0.8, 0.85, 0.88, 0.9, 0.92, 0.94, 0.96, 0.97, 0.98, 0.99, 1.};
      BinningX = MakeVector(BinArray_x);
      BinningY = MakeVector(BinArray_y);
      break;
    }

  case (7):
  case (16):
    {
      const double BinArray_x[] = {0., 600., 800., 1500., 30000.};
      const double BinArray_y[] = {-1, 0.7, 0.8, 0.86, 0.9, 0.94, 0.96, 0.97, 0.98, 0.99, 1};
      BinningX = MakeVector(BinArray_x);
      BinningY = MakeVector(BinArray_y);
      break;
    }

  case (8):
  case (17):
    {
      const double BinArray_x[] = {0., 600., 1000., 1250., 2000., 4000., 30000.};
      const double BinArray_y[] = {-1., 0.7, 0.8, 0.86, 0.9, 0.93, 0.95, 0.97, 0.99, 1.};
      BinningX = MakeVector(BinArray_x);
      BinningY = MakeVector(BinArray_y);
      break;
    }

  default:
    std::cerr << "You gave me a selection I don't have a binning for, whaaa?" << std::endl;
    std::cerr << "Selection: " << Selection << " = " << std::endl;
    throw;
    break;
  }
}

// **************************************************
// Helper function for projecting TH2Poly onto the X axis
TH1D* PolyProjectionX(TObject* poly, std::string TempName, std::vector<double> xbins) {
  // **************************************************

  TH1D* hProjX = new TH1D((TempName+"_x").c_str(),(TempName+"_x").c_str(),xbins.size()-1,&xbins[0]);
  double xlow, xup, frac=0;

  //loop over bins in the poly
  for(int i=0; i<((TH2Poly*)poly)->GetNumberOfBins(); i++)
    {
      //get bin and its edges
      TH2PolyBin* bin = (TH2PolyBin*)((TH2Poly*)poly)->GetBins()->At(i)->Clone();
      xlow=bin->GetXMin();
      xup=bin->GetXMax();

      //Loop over projected bins, find fraction of poly bin in each
      for(int dx=0; dx<int(xbins.size()); dx++)
        {
          if(xbins[dx+1]<=xlow || xbins[dx]>=xup)
            {
              frac=0;
            }
          else if(xbins[dx]<=xlow && xbins[dx+1]>=xup)
            {
              frac=1;
            }
          else if(xbins[dx]<=xlow && xbins[dx+1]<=xup)
            {
              frac=(xbins[dx+1]-xlow)/(xup-xlow);
            }
          else if(xbins[dx]>=xlow && xbins[dx+1]>=xup)
            {
              frac=(xup-xbins[dx])/(xup-xlow);
            }
          else if(xbins[dx]>=xlow && xbins[dx+1]<=xup)
            {
              frac=(xbins[dx+1]-xbins[dx])/(xup-xlow);
            }
          else
            {
              frac=0;
            }
          hProjX->SetBinContent(dx+1,hProjX->GetBinContent(dx+1)+frac*bin->GetContent());
	  hProjX->SetBinError(dx+1,hProjX->GetBinError(dx+1)+(frac*poly->GetBinError(i+1))*(frac*poly->GetBinError(i+1)));
        }
    }
  
  for(int i=1; i<=hProjX->GetXaxis()->GetNbins(); i++){
    hProjX->SetBinError(i, sqrt(hProjX->GetBinError(i)));
  }

  return hProjX;
} // end project poly X function

// **************************************************
// Helper function for projecting TH2Poly onto the Y axis
TH1D* PolyProjectionY(TObject* poly, std::string TempName, std::vector<double> ybins) {
  // **************************************************

  TH1D* hProjY = new TH1D((TempName+"_y").c_str(),(TempName+"_y").c_str(),ybins.size()-1,&ybins[0]);
  double ylow, yup, frac=0;

  //loop over bins in the poly
  for(int i=0; i<((TH2Poly*)poly)->GetNumberOfBins(); i++)
    {
      //get bin and its edges
      TH2PolyBin* bin = (TH2PolyBin*)((TH2Poly*)poly)->GetBins()->At(i)->Clone();
      ylow=bin->GetYMin();
      yup=bin->GetYMax();

      //Loop over projected bins, find fraction of poly bin in each
      for(int dy=0; dy<int(ybins.size()); dy++)
        {
          if(ybins[dy+1]<=ylow || ybins[dy]>=yup)
            {
              frac=0;
            }
          else if(ybins[dy]<=ylow && ybins[dy+1]>=yup)
            {
              frac=1;
            }
          else if(ybins[dy]<=ylow && ybins[dy+1]<=yup)
            {
              frac=(ybins[dy+1]-ylow)/(yup-ylow);
            }
          else if(ybins[dy]>=ylow && ybins[dy+1]>=yup)
            {
              frac=(yup-ybins[dy])/(yup-ylow);
            }
          else if(ybins[dy]>=ylow && ybins[dy+1]<=yup)
            {
              frac=(ybins[dy+1]-ybins[dy])/(yup-ylow);
            }
          else
            {
              frac=0;
            }
          hProjY->SetBinContent(dy+1,hProjY->GetBinContent(dy+1)+frac*bin->GetContent());
	  hProjY->SetBinError(dx+1,hProjY->GetBinError(dx+1)+(frac*poly->GetBinError(i+1))*(frac*poly->GetBinError(i+1)))
        }
    }

  for(int i=1; i<=hProjY->GetXaxis()->GetNbins(); i++){
    hProjY->SetBinError(i, sqrt(hProjY->GetBinError(i)));
  }

  return hProjY;
} // end project poly Y function

