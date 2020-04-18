#include <iostream>
#include <vector>
#include <string>
#include <sstream>

#include "TFile.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TStyle.h"

#include "TVectorD.h"
#include "TVectorT.h"
#include "TVector.h"

#include "TMatrixD.h"
#include "TMatrixT.h"
#include "TMatrix.h"

#include "TObjArray.h"
#include "TObjString.h"
#include "TString.h"

#include "TLegend.h"

// Declare these ahead
std::vector<std::string> MakeBANFFNames();
std::vector<std::string> MakeMaCh3Names();
std::vector<double> GetXSecNominal();

void MaCh3_LLH(std::string &);

int main(int argc, char ** argv) {

  if (argc != 2) {
    std::cerr << "Need 1 arguments to run..." << std::endl;
    std::cerr << "MaCh3_file" << std::endl;
    throw;
  }

  std::string MaCh3 = std::string(argv[1]);

  MaCh3_LLH(MaCh3);

  return 0;
}


// *******************
// The main function
void MaCh3_LLH(std::string &MaCh3File) {
  // *******************
  std::cout << "MaCh3 file: " << MaCh3File << std::endl;

  TFile *MaCh3 = new TFile(MaCh3File.c_str());

  TCanvas *Canvas = new TCanvas("canv", "canv", 1024, 1024);
  Canvas->SetBottomMargin(0.08);
  Canvas->SetTopMargin(0.06);
  Canvas->SetLeftMargin(0.15);
  Canvas->SetRightMargin(0.03);
  Canvas->SetGridx();
  Canvas->SetGridy();
  std::string pdfname = MaCh3File;
  //Canvas->Print((pdfname+".pdf[").c_str());

  // Get all the BANFF names, woop
  std::vector<std::string> MaCh3Names = MakeMaCh3Names();

  // BANFF and MaCh3 have different values for nominal
  // In BANFF we always use the spline value for the parameter value, so MAQE = 0 means MAQE nominal
  // In MaCh3 we instead have everything relative nominal, so MAQE = 1 means nominal
  // To get the nominal values we read from file instead of hard-coding our life
  // These are also re-ordered to put FSI before the xsec, like BANFF do
  //std::vector<double> XSecNominals = GetXSecNominal();

  // Now need to map these onto MaCh3 names somehow...
  unsigned int Params = MaCh3Names.size();
  //  int nXsec = 0;
  for (unsigned int i = 0; i < Params; ++i) {

    TH1D *FluxGraph = (TH1D*)(MaCh3->Get((std::string("Flux_LLH/")+MaCh3Names[i]+"_flux").c_str())->Clone());
    TH1D *XsecGraph = (TH1D*)(MaCh3->Get((std::string("XSec_LLH/")+MaCh3Names[i]+"_xs").c_str())->Clone());
    TH1D *ND280Graph = (TH1D*)(MaCh3->Get((std::string("ND280_LLH/")+MaCh3Names[i]+"_nd").c_str())->Clone());
    TH1D *SampleGraph = (TH1D*)(MaCh3->Get((std::string("Sample_LLH/")+MaCh3Names[i]+"_sam").c_str())->Clone());
    TH1D *TotalGraph = (TH1D*)(MaCh3->Get((std::string("Total_LLH/")+MaCh3Names[i]+"_full").c_str())->Clone());

    TotalGraph->SetTitle("");
    TotalGraph->SetName("");
    ND280Graph->SetTitle("");
    XsecGraph->SetTitle("");
    FluxGraph->SetTitle("");
    SampleGraph->SetTitle("");
    if (TotalGraph == NULL) std::cout << MaCh3Names[i] << std::endl;
    TotalGraph->Draw("c");
    if (FluxGraph->Integral() > 1E-6) FluxGraph->Draw("same,c");
    if (XsecGraph->Integral() > 1E-6) XsecGraph->Draw("same,c");
    if (ND280Graph->Integral() > 1E-6) ND280Graph->Draw("same,c");
    if (SampleGraph->Integral() > 1E-6) SampleGraph->Draw("same,c");
    gStyle->SetOptStat(0);

    TotalGraph->SetMinimum(0);
    TotalGraph->SetMaximum(TotalGraph->GetMaximum()*1.1);
    TotalGraph->SetFillStyle(0);
    TotalGraph->SetFillColor(0);

    //    TotalGraph->SetTitle(MaCh3Names[i].c_str());
    TotalGraph->SetLineStyle(kSolid);
    TotalGraph->SetLineWidth(3);
    TotalGraph->SetLineColor(kRed+2);
    TotalGraph->GetXaxis()->SetTitle("Variation rel. nom.");
    TotalGraph->GetYaxis()->SetTitle("#chi^{2}");
    TotalGraph->GetYaxis()->SetTitleOffset(2.0);

    FluxGraph->SetLineColor(kGreen+3);
    FluxGraph->SetFillStyle(0);
    FluxGraph->SetFillColor(0);
    FluxGraph->SetLineWidth(3);
    XsecGraph->SetLineColor(kGreen+3);
    XsecGraph->SetFillStyle(0);
    XsecGraph->SetFillColor(0);
    XsecGraph->SetLineWidth(3);
    ND280Graph->SetLineColor(kGreen+3);
    ND280Graph->SetFillStyle(0);
    ND280Graph->SetFillColor(0);
    ND280Graph->SetLineWidth(3);

    FluxGraph->SetLineStyle(kDashed);
    XsecGraph->SetLineStyle(kDashed);
    ND280Graph->SetLineStyle(kDashed);

    SampleGraph->SetLineColor(kBlue+2);
    SampleGraph->SetFillStyle(0);
    SampleGraph->SetFillColor(0);
    SampleGraph->SetLineWidth(3);
    SampleGraph->SetLineStyle(kDashed);

    // Add a legend
    TLegend *leg = new TLegend(0.45, 0.5, 0.8, 0.75);
    leg->AddEntry(TotalGraph, "Total");
    if (FluxGraph->Integral() > 1E-6) leg->AddEntry(FluxGraph, "Prior, Flux");
    if (XsecGraph->Integral() > 1E-6) leg->AddEntry(XsecGraph, "Prior, Int.");
    if (ND280Graph->Integral() > 1E-6) leg->AddEntry(ND280Graph, "Prior, ND280");
    if (SampleGraph->Integral() > 1E-6) leg->AddEntry(SampleGraph, "Samples");
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->SetLineWidth(0);
    leg->SetLineColor(0);
    leg->SetBorderSize(0);
    leg->Draw("same");

    Canvas->Print(("llh/"+MaCh3Names[i]+"_llh.pdf").c_str());

    delete leg;
  }
  //  Canvas->Print((pdfname+".pdf]").c_str());
}


// *******************
// Make vector for the BANFF TGraph names
std::vector<std::string> MakeBANFFNames() {
  // *******************
  std::vector<std::string> NameVector;
  // ND280 numu flux, 11 of them
  for (int i = 0; i <= 10; ++i) {
    std::stringstream ss;
    ss << "NDNuModeNumu" << i;
    NameVector.push_back(ss.str());
  }

  // ND280 numub flux, 5 of them
  for (int i = 0; i <= 4; ++i) {
    std::stringstream ss;
    ss << "NDNuModeNumub" << i;
    NameVector.push_back(ss.str());
  }

  // ND nue flux, 7 of them
  for (int i = 0; i <= 6; ++i) {
    std::stringstream ss;
    ss << "NDNuModeNue" << i;
    NameVector.push_back(ss.str());
  }

  // ND nueb flux, 2 of them
  for (int i = 0; i <= 1; ++i) {
    std::stringstream ss;
    ss << "NDNuModeNueb" << i;
    NameVector.push_back(ss.str());
  }

  // ND anu numu, 5 of them
  for (int i = 0; i <= 4; ++i) {
    std::stringstream ss;
    ss << "NDANuModeNumu" << i;
    NameVector.push_back(ss.str());
  }

  // ND anu numub, 11 of them
  for (int i = 0; i <= 10; ++i) {
    std::stringstream ss;
    ss << "NDANuModeNumub" << i;
    NameVector.push_back(ss.str());
  }

  // ND anu nue, 2 of them
  for (int i = 0; i <= 1; ++i) {
    std::stringstream ss;
    ss << "NDANuModeNue" << i;
    NameVector.push_back(ss.str());
  }

  // ND anu nueb, 7 of them
  for (int i = 0; i <= 6; ++i) {
    std::stringstream ss;
    ss << "NDANuModeNueb" << i;
    NameVector.push_back(ss.str());
  }

  // Now get the size, and replace the "ND" part with "SK" to make the SK names
  unsigned int TempSize = NameVector.size();
  for (unsigned int i = 0; i < TempSize; ++i) {
    std::string TempName = NameVector[i];
    TempName.replace(TempName.find("ND"), 2, std::string("SK"));
    NameVector.push_back(TempName);
  }

  // Now make the xsec, fully hardcoded hurray (BANFF numbering do not necessarily follow a convention?)
  NameVector.push_back("FEFQE");
  NameVector.push_back("FEFQEH");
  NameVector.push_back("FEFINEL");
  NameVector.push_back("FEFABS");
  NameVector.push_back("FEFCX");

  NameVector.push_back("MAQE");

  NameVector.push_back("2p2h_norm_nu");
  NameVector.push_back("2p2h_norm_nubar");
  NameVector.push_back("2p2h_normCtoO");
  NameVector.push_back("2p2h_shape_C");
  NameVector.push_back("2p2h_shape_O");

  NameVector.push_back("2p2h_Edep_lowEnu");
  NameVector.push_back("2p2h_Edep_highEnu");
  NameVector.push_back("2p2h_Edep_lowEnubar");
  NameVector.push_back("2p2h_Edep_highEnubar");

  NameVector.push_back("Q2_norm_0");
  NameVector.push_back("Q2_norm_1");
  NameVector.push_back("Q2_norm_2");
  NameVector.push_back("Q2_norm_3");
  NameVector.push_back("Q2_norm_4");

  NameVector.push_back("EB_dial_C_nu");
  NameVector.push_back("EB_dial_O_nu");
  NameVector.push_back("EB_dial_C_nubar");
  NameVector.push_back("EB_dial_O_nubar");

  NameVector.push_back("CA5");
  NameVector.push_back("MARES");
  NameVector.push_back("ISO_BKG_LowPPi");
  NameVector.push_back("ISO_BKG");

  NameVector.push_back("CC_norm_nu");
  NameVector.push_back("CC_norm_nubar");

  NameVector.push_back("nue_numu");
  NameVector.push_back("nuebar_numubar");

  NameVector.push_back("CC_BY_DIS");
  NameVector.push_back("CC_BY_MPi");
  NameVector.push_back("CC_AGKY_Mult");
  NameVector.push_back("CC_Misc");

  NameVector.push_back("CC_DIS_MultPi_Norm_Nu");
  NameVector.push_back("CC_DIS_MultPi_Norm_Nubar");

  NameVector.push_back("CC_Coh_C");
  NameVector.push_back("CC_Coh_O");

  NameVector.push_back("NC_Coh");
  NameVector.push_back("NC_1gamma");
  NameVector.push_back("NC_other_near");
  NameVector.push_back("NC_other_far");

  // Now make the OBSNORM, fairly easy
  for (int i = 0; i <= 574; ++i) {
    std::stringstream ss;
    ss << "OBSNORM_" << i;
    NameVector.push_back(ss.str());
  }

  return NameVector;

}

// *******************
// Make vector for the MaCh3 TH1D names
std::vector<std::string> MakeMaCh3Names() {
  // *******************

  std::vector<std::string> NameVector;
  // The beam parameters
  for (int i = 0; i <= 99; ++i) {
    std::stringstream ss;
    ss << "b_" << i;
    NameVector.push_back(ss.str());
  }

  // The cross-section parameters
  // First FSI to match BANFF
  NameVector.push_back("FEFQE");
  NameVector.push_back("FEFQEH");
  NameVector.push_back("FEFINEL");
  NameVector.push_back("FEFABS");
  NameVector.push_back("FEFCX");

  NameVector.push_back("MAQE");

  NameVector.push_back("2p2h_norm_nu");
  NameVector.push_back("2p2h_norm_nubar");
  NameVector.push_back("2p2h_normCtoO");
  NameVector.push_back("2p2h_shape_C");
  NameVector.push_back("2p2h_shape_O");

  NameVector.push_back("2p2h_Edep_lowEnu");
  NameVector.push_back("2p2h_Edep_highEnu");
  NameVector.push_back("2p2h_Edep_lowEnubar");
  NameVector.push_back("2p2h_Edep_highEnubar");

  NameVector.push_back("Q2_norm_0");
  NameVector.push_back("Q2_norm_1");
  NameVector.push_back("Q2_norm_2");
  NameVector.push_back("Q2_norm_3");
  NameVector.push_back("Q2_norm_4");
  NameVector.push_back("Q2_norm_5");
  NameVector.push_back("Q2_norm_6");
  NameVector.push_back("Q2_norm_7");

  NameVector.push_back("EB_dial_C_nu");
  NameVector.push_back("EB_dial_O_nu");
  NameVector.push_back("EB_dial_C_nubar");
  NameVector.push_back("EB_dial_O_nubar");

  NameVector.push_back("CA5");
  NameVector.push_back("MARES");
  NameVector.push_back("ISO_BKG_LowPPi");
  NameVector.push_back("ISO_BKG");

  NameVector.push_back("CC_norm_nu");
  NameVector.push_back("CC_norm_nubar");

  NameVector.push_back("nue_numu");
  NameVector.push_back("nuebar_numubar");

  NameVector.push_back("CC_BY_DIS");
  NameVector.push_back("CC_BY_MPi");
  NameVector.push_back("CC_AGKY_Mult");
  NameVector.push_back("CC_Misc");

  NameVector.push_back("CC_DIS_MultPi_Norm_Nu");
  NameVector.push_back("CC_DIS_MultPi_Norm_Nubar");

  NameVector.push_back("CC_Coh_C");
  NameVector.push_back("CC_Coh_O");

  NameVector.push_back("NC_Coh");
  NameVector.push_back("NC_1gamma");
  NameVector.push_back("NC_other_near");
  NameVector.push_back("NC_other_far");

  // The ND280 parameters
  for (int i = 0; i <= 573; ++i) {
    std::stringstream ss;
    ss << "ndd_" << i;
    NameVector.push_back(ss.str());
  }

  return NameVector;
}

// *******************
// Get the xsec nominals
std::vector<double> GetXSecNominal() {
  // *******************

  std::string XsecCovName = "../inputs/xsec_covariance_2020a_v7.root";
  TFile *XsecFile = new TFile(XsecCovName.c_str());
  // Get the nominal vector so we can shift all the MaCh3 values to BANFF values
  TVectorT<double> *xsec_param_nom = (TVectorD*)(XsecFile->Get("xsec_param_nom")->Clone());
  // Also need to find which parameters are spline parameters
  TMatrixT<double> *xsec_param_id = (TMatrixD*)(XsecFile->Get("xsec_param_id")->Clone());
  // Doesn't hurt to get the names whilst we're at it
  TObjArray *xsec_param_names = (TObjArray*)(XsecFile->Get("xsec_param_names")->Clone());

  int size = xsec_param_nom->GetNrows();
  std::vector<double> ParamNominals(size, 0.0);
  for (int i = 0; i < size; ++i) {
    int id = int((*xsec_param_id)(i, 0));
    if (id < 0) {
      ParamNominals[i] = 0.0;
    } else {
      double NominalValue = (*xsec_param_nom)(i);
      ParamNominals[i] = NominalValue;
    }
  }

  // Then set the FSI parameters to be in the front, yuck
  std::vector<double> ReOrderNominals(size, 0.0);
  std::vector<std::string> ReOrderNames(size, "");

  for (int i = 0; i < size; ++i) {
    std::string ParamName = std::string(((TObjString*)(xsec_param_names->At(i)))->GetString());
    if (ParamName == std::string("FEFQE")) {
      ReOrderNominals[0] = ParamNominals[i];
      ReOrderNames[0] = ParamName;
    } else if (ParamName == std::string("FEFQEH")) {
      ReOrderNominals[1] = ParamNominals[i];
      ReOrderNames[1] = ParamName;
    } else if (ParamName == std::string("FEFINEL")) {
      ReOrderNominals[2] = ParamNominals[i];
      ReOrderNames[2] = ParamName;
    } else if (ParamName == std::string("FEFABS")) {
      ReOrderNominals[3] = ParamNominals[i];
      ReOrderNames[3] = ParamName;
    } else if (ParamName == std::string("FEFCX")) {
      ReOrderNominals[4] = ParamNominals[i];
      ReOrderNames[4] = ParamName;
    //} else if (ParamName == std::string("FSI_CEX_HI")) {
      //ReOrderNominals[5] = ParamNominals[i];
      //ReOrderNames[5] = ParamName;
    } else {
      ReOrderNominals[i+6] = ParamNominals[i];
      ReOrderNames[i+6] = ParamName;
    }
  }

  return ReOrderNominals;
}
