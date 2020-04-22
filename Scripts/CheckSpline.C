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


void CheckSpline(std::string fileloc) {

  //  std::string runs[11] = {"2a","2w","3","4a","4w","5","6","7","8a","8w","9"};//2020
  std::string runs[13] = {"2a","2w","3b","3c","4a","5","6b","6c","6d","6e","7","8a","8w"};

  //count events w maqe !=1.0
  int runcount=0, totalcount=0;
  //count events w maqe !=1.0 + reac not -1
  int countnonn1=0, totalcountnonn1=0;
  //count good events w reac -1
  int okcountn1=0, totalokcountn1=0;
  //count events w maqe !=1.0 + target not H
  int countnonH=0, totalcountnonH=0;
  //count good events + target H 
  int okcountH=0, totalokcountH=0;

  for (int r=0; r<13; r++){
    std::string filename = fileloc+"/MCRun"+runs[r]+"_Splines.root";
    TFile *f0 = TFile::Open(filename.c_str());
    
    runcount=0;
    countnonn1=0;
    okcountn1=0;
    countnonH=0;
    okcountH=0;
    double weight=0;
    int code=0, target=0;
    TTree* ssum = (TTree*)f0->Get("sample_sum");
    TClonesArray* arr = 0;
    ssum->SetBranchAddress("MAQEGraph",&arr);
    //ssum->SetBranchAddress("ReactionCode",&code);
    //   ssum->SetBranchAddress("TgtMat",&target);
    
    for(int i=0; i<ssum->GetEntries(); i++){
      ssum->GetEntry(i);
      weight = ((TGraph*)arr->At(0))->Eval(0);

      if(weight<0.9999||weight>1.0001){
	runcount++;
	if(code!=-1 ){
	  countnonn1++;
	}
	if(target!=1)
	  countnonH++;
      }
      else {
	if (code==-1){
	  okcountn1++;
	}
	if(target==1 && code==-1){
	  okcountH++;
	}
      }
    }
    
    totalcountnonH+=countnonH;
    totalokcountH+=okcountH;
    totalokcountn1+=okcountn1;
    totalcount+=runcount;
    totalcountnonn1+=countnonn1;
    totalcountnonH+=countnonH;

    std::cout << "Run " << runs[r] << ": " << std::endl;
    std::cout << runcount << "/" << ssum->GetEntries() << " MAQE !=1.0 at 0. "<< std::endl;
    std::cout << countnonn1  << " for mode not -1 " << std::endl;
    std::cout << okcountn1 << " reac -1 with ok maqe spline. " << std::endl; 
    std::cout << countnonH << " for target not H" << std::endl;
    std::cout << okcountH << " ok splines on H" << std::endl;

    f0->Close();
    delete f0;
  }

  std::cout << "Total: " << std::endl << totalcount << " MAQE !=1.0 at 0. " << std::endl;
  std::cout << totalcountnonn1  << " for mode not -1 " << std::endl;
  std::cout << totalokcountn1 << " reac -1 with ok maqe spline. " << std::endl;
  std::cout << totalokcountH << "target H " << std::endl;
  std::cout << totalcountnonH  << " maqe bad on not H. " << std::endl;
  
  std::cout << "done" << std::endl;
}
