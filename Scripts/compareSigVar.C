// Compare sigma variations to nominal
// .L compareSigVar.C+g // in root
// compareSigVar("filename.root" 
// ********************************

#include "TH1D.h"
#include "TH2Poly.h"
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


void compareSigVar(std::string filename) {

  std::string OldName;
  TCanvas *c = new TCanvas("canv", "canv", 1080, 1080);
  TFile *infile = TFile::Open(filename.c_str());
  std::string outfilename = filename.substr(0, filename.find(".root"));
  outfilename = outfilename + "_RatioPlots.pdf";
  
  c->SetTopMargin(0.07);
  c->SetBottomMargin(0.11);
  c->SetLeftMargin(0.18);
  c->SetRightMargin(0.23);
  
  c->Print((outfilename+"[").c_str());

  OldName = "empty";
 
  std::string dirname;
  TIter next(infile->GetListOfKeys());
  TKey *key;
  
  int counter=0;
  TDirectory *dir;

  int count2=0;
  // Loop through all entries
  while ((key = (TKey*)next())) {
    
    // get names, ignore flux 
    std::string classname = std::string(key->GetClassName());   
    dirname = std::string(key->GetName());  
    if (dirname.find("_b_") != std::string::npos) continue;

    // If directory go into it
    if (classname == "TDirectoryFile") {
      count2++;
      std::cout << "Entering " << dirname << std::endl;
      infile->cd(dirname.c_str());
      dir=gDirectory;
      counter=0;
      
      TH2Poly *m3sig;
      TH2Poly *m1sig;
      TH2Poly *nom;
      TH2Poly *p1sig;
      TH2Poly *p3sig;
      TIter nextsub(dir->GetListOfKeys());
      TKey *subkey;
    
      //loop over items in directory, hard code which th2poly we want
      while ((subkey = (TKey*)nextsub())) {
	std::string name = std::string(subkey->GetName());
	classname = std::string(subkey->GetClassName());
	if(name.substr( name.length() - 5 ) == "_2p2h") continue;
	
	if (classname == "TH2Poly") {
	  counter++;
	  
	  //Order saved in output in src/ND280_SigmaVar_2019.cpp is 2,4,6,8,10 are -3,-1,0,+1,+3 sigma
	  if(counter==1){
	    std::cout << "Getting " << name << std::endl;
	    name = dirname + "/" + name;
	    m3sig=(TH2Poly*)infile->Get(name.c_str())->Clone();
	  }
	  if(counter==2){
	    std::cout << "Getting " << name << std::endl;
	    name = dirname + "/" + name;
	    m1sig=(TH2Poly*)infile->Get(name.c_str())->Clone();
	  }
	  if(counter==3){
	    std::cout << "Getting " << name << std::endl;
	    name = dirname + "/" + name;
	    nom=(TH2Poly*)infile->Get(name.c_str())->Clone();
	  }
	  if(counter==4){
	    std::cout << "Getting " << name << std::endl;
	    name = dirname + "/" + name;
	    p1sig=(TH2Poly*)infile->Get(name.c_str())->Clone();
	  }
	  if(counter==5){
	    std::cout << "Getting " << name << std::endl;
	    name = dirname + "/" + name;
	    p3sig=(TH2Poly*)infile->Get(name.c_str())->Clone();
	  }
	}
      }	
      
      //Divide the plots for ratios and plot to a big old pdf
      m3sig->Divide(nom);
      m3sig->SetName((dirname + "_-3sigmaRatioToNom").c_str());
      m3sig->SetTitle((dirname + "_-3sigmaRatioToNom").c_str());
      c->Clear();
      m3sig->GetXaxis()->SetRangeUser(0,5000);
      m3sig->GetXaxis()->SetTitleOffset(1.1);
      m3sig->GetYaxis()->SetTitleOffset(1.1);
      m3sig->GetZaxis()->SetTitleOffset(1.5);
      m3sig->GetZaxis()->SetRangeUser(0.7,1.3);
      m3sig->GetZaxis()->SetTitle("Ratio to Nominal");
      m3sig->Draw("colz");
      c->Print((outfilename).c_str());
	
      m1sig->Divide(nom);
      m1sig->SetName((dirname + "_-1sigmaRatioToNom").c_str());
      m1sig->SetTitle((dirname + "_-1sigmaRatioToNom").c_str());
      c->Clear();
      m1sig->GetXaxis()->SetRangeUser(0,5000);
      m1sig->GetXaxis()->SetTitleOffset(1.1);
      m1sig->GetYaxis()->SetTitleOffset(1.1);
      m1sig->GetZaxis()->SetTitleOffset(1.5);
      m1sig->GetZaxis()->SetRangeUser(0.7,1.3);
      m1sig->GetZaxis()->SetTitle("Ratio to Nominal");
      m1sig->Draw("colz");
      c->Print((outfilename).c_str());

      p1sig->Divide(nom);
      p1sig->SetName((dirname + "_+1sigmaRatioToNom").c_str());
      p1sig->SetTitle((dirname + "_+1sigmaRatioToNom").c_str());
      c->Clear();
      p1sig->GetXaxis()->SetRangeUser(0,5000);
      p1sig->GetXaxis()->SetTitleOffset(1.1);
      p1sig->GetYaxis()->SetTitleOffset(1.1);
      p1sig->GetZaxis()->SetTitleOffset(1.5);
      p1sig->GetZaxis()->SetRangeUser(0.7,1.3);
      p1sig->GetZaxis()->SetTitle("Ratio to Nominal");
      p1sig->Draw("colz");
      c->Print((outfilename).c_str());

      p3sig->Divide(nom);
      p3sig->SetName((dirname + "_+3sigmaRatioToNom").c_str());
      p3sig->SetTitle((dirname + "_+3sigmaRatioToNom").c_str());
      c->Clear();
      p3sig->GetXaxis()->SetRangeUser(0,5000);
      p3sig->GetXaxis()->SetTitleOffset(1.1);
      p3sig->GetYaxis()->SetTitleOffset(1.1);
      p3sig->GetZaxis()->SetTitleOffset(1.5);
      p3sig->GetZaxis()->SetRangeUser(0.7,1.3);
      p3sig->GetZaxis()->SetTitle("Ratio to Nominal");
      p3sig->Draw("colz");
      c->Print((outfilename).c_str());
   
      gDirectory->cd("..");

      delete m3sig;
      delete m1sig;
      delete nom;
      delete p1sig;
      delete p3sig;
      
    }
    else{
      continue;
    }
    //    if(count2 >42)
    //break;
  }
  c->Print((outfilename+"]").c_str());
}

