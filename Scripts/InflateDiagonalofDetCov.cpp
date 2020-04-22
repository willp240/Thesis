// Convert ND280 binning
// Compiled because that means we can link to psyche and avoid some hard-coding!
// g++ `root-config --cflags` -g -o ConvertND280Bin ConvertND280Bin.cpp -I`root-config --incdir` -I../psycheinc `root-config --glibs --libs` -L../lib
// ********************************
// ______WARNING______
// THE ORDER OF SAMPLES IN THE FOR LOOP MATTERS AND ____IS HARD CODED____
// FOR 2017 ANALYSES THE SAMPLE ORDERING WAS FGD1 SAMPLES, FGD2 SAMPLES, GOING __AGAINST__ THE PSYCHE NUMBERING
// FOR 2018 ANALYSES THE SAMPLE ORDERING FOLLOWED PSYCHE NUMBERING
// CHECK GIT LOG TO SEE WHAT WAS DONE
// ______END WARNING______

#include <iostream>
#include <iomanip>
#include <string>

#include "TFile.h"
#include "TAxis.h"
#include "TObjArray.h"
#include "TH1I.h"
#include "TVectorD.h"
#include "TMatrixDSym.h"

#include "CoreDataClasses.hxx"
#include "CoreUtils.hxx"
#include "BaseDataClasses.hxx"

void InflateDiagonalofDetCov(const std::string &);


int main(int argc, char **argv) {

  if (argc != 2) {
    std::cout << "./InflateDiagonalofDetCov inputFile.root" << std::endl;
    exit(-1);
  }

  std::string infile = argv[1];
  InflateDiagonalofDetCov(infile);

  return 0;
};


void InflateDiagonalofDetCov(const std::string & inputFile) {

  TFile *f = new TFile(inputFile.c_str(), "READ");
  if (f->IsZombie()) {
    std::cerr << "Can't find file " << inputFile << std::endl;
    std::cerr << "I'll just have to hang up!" << std::endl;
    exit(-1);
  }

  TObjArray* k1_axes_orig = (TObjArray*)(f->Get("k1_axes")->Clone());
  TObjArray* k2_axes_orig = (TObjArray*)(f->Get("k2_axes")->Clone());
  TH1I* offset_orig = (TH1I*)(f->Get("offset")->Clone());
  TMatrixDSym* nddet_cov_infl = (TMatrixDSym*)(f->Get("nddet_cov")->Clone());
  TVectorD* det_weights_orig = (TVectorD*)(f->Get("det_weights")->Clone());

  for(int i=0; i<nddet_cov_infl->GetNrows(); i++){
    (*nddet_cov_infl)[i][i]+= 1.0e-7;
  }

  std::string outputName = inputFile;
  while (outputName.find(".") != std::string::npos) {
    outputName = outputName.substr(0, outputName.find("."));
  }
  outputName += "_inflated.root";

  TFile *output = new TFile(outputName.c_str(), "RECREATE");

  k1_axes_orig->Write("k1_axes", 1);
  k2_axes_orig->Write("k2_axes", 1);
  offset_orig->Write("offset");
  nddet_cov_infl->Write("nddet_cov");
  det_weights_orig->Write("det_weights");

  output->Close();
  f->Close();

};
