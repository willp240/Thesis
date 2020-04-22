#include "TArrayD.h"
void SKSpectra_Divide2Panel4Spectra(TString file_FD, TString file_FD2, TString file_FD3, TString file_Nom, TString outfile, TString spectrum, int counter)
{
  gStyle->SetOptStat(0);
  //gStyle->SetOptTitle(0);
  //gStyle->SetTextFont(132);
  TGaxis::SetMaxDigits(3) ;

  TString nominal = "2198 Bins";
  TString fd1 = "1347 Bins";
  TString fd2 = "179 Bins";
  TString fd3 = "Peak 3: EbC 33-37MeV";
  /*  TString nominal = "2198 Bins";
  TString fd1 = "1825 Bins";
  TString fd2 = "574 Bins";
  TString fd3 = "574 Bins";*/

  TFile *fout = new TFile(outfile,"recreate");

  TString name = "hist_1d_"+spectrum;
  std::cout << "name " << name << std::endl;
  // Get fake data spectrum
  TFile *f_FD = new TFile(file_FD,"open");
  TString FDName = "FD" + name;
  std::cout << "name " << name << std::endl;
  TH1D *h_FD = (TH1D*)f_FD->Get(name);
  h_FD->SetName(FDName);
  h_FD->Sumw2();

  TFile *f_FD2 = new TFile(file_FD2,"open");
  TString FD2Name = "FD2" + name;

  TH1D *h_FD2 = (TH1D*)f_FD2->Get(name);
  h_FD2->SetName(FD2Name);
  h_FD2->Sumw2();

  TFile *f_FD3 = new TFile(file_FD3,"open");
  TString FD3Name = "FD3" + name;

  TH1D *h_FD3 = (TH1D*)f_FD3->Get(name);
  h_FD3->SetName(FD3Name);
  h_FD3->Sumw2();

  // Get Nom spectrum
  TFile *f_Nom = new TFile(file_Nom,"open");
  TString Nomname = "Nom" + name;

  TH1D *h_Nom = (TH1D*)f_Nom->Get(name);
  h_Nom->SetName(Nomname);
  h_Nom->Sumw2();
  
  if (h_Nom->GetXaxis()->GetNbins() != h_FD->GetXaxis()->GetNbins() != h_FD2->GetXaxis()->GetNbins() != h_FD3->GetXaxis()->GetNbins())
    std::cerr << "Error! " << h_Nom->GetXaxis()->GetNbins() << " bins in Nom histogram and " << h_FD->GetXaxis()->GetNbins() << " bins in FD histogram, " << h_FD2->GetXaxis()->GetNbins() << " bins in FD2 histogram and " << h_FD3->GetXaxis()->GetNbins() << " bins in FD3 histogram" << std::endl;

  // Set draw options to be applied to all plots
  h_FD->SetLineColor(kRed);
  h_FD->SetLineWidth(2);
  h_FD->SetMarkerStyle(20);
  h_FD2->SetLineColor(kBlue);
  h_FD2->SetLineWidth(2);
  h_FD2->SetMarkerStyle(20);
  h_FD3->SetLineColor(kGreen);
  h_FD3->SetLineWidth(2);
  h_FD3->SetMarkerStyle(20);
  h_Nom->SetLineColor(kBlack);
  h_Nom->SetLineWidth(2);

  h_FD->GetXaxis()->SetTitleSize(0.07);
  h_FD->GetXaxis()->SetLabelSize(0.06);
  h_FD2->GetYaxis()->SetTitleSize(0.07);
  h_FD2->GetYaxis()->SetLabelSize(0.06);
  h_FD3->GetYaxis()->SetTitleSize(0.07);
  h_FD3->GetYaxis()->SetLabelSize(0.06);

  // True histogram for dividing
  TH1D *Nom_copy = (TH1D*)h_Nom->Clone("Nom_copy");

  // Histograms for dividing by True
  TH1D *h_FDn = (TH1D*)h_FD->Clone("FDn");
  TH1D *h_FD2n = (TH1D*)h_FD2->Clone("FD2n");
  TH1D *h_FD3n = (TH1D*)h_FD3->Clone("FD3n");
  TH1D *h_Nomn = (TH1D*)h_Nom->Clone("Nomn");
  
  // Histograms for normalising by bin width
  TH1D *h_FDn_un = (TH1D*)h_FD->Clone("FDn_un");
  TH1D *h_FD2n_un = (TH1D*)h_FD2->Clone("FD2n_un");
  TH1D *h_FD3n_un = (TH1D*)h_FD3->Clone("FD3n_un");
  TH1D *h_Nomn_un = (TH1D*)h_Nom->Clone("Nomn_un");
  
  // Histograms for normalising difference to nominal by uncertainty
  TH1D *h_FDn_diff = (TH1D*)h_FD->Clone("FDn_diff");
  TH1D *h_FD2n_diff = (TH1D*)h_FD2->Clone("FD2n_diff");
  TH1D *h_FD3n_diff = (TH1D*)h_FD3->Clone("FD3n_diff");
  double FDerror =0, FDdiff=0, FDval=0, FD2error =0, FD2diff=0, FD2val=0, FD3error =0, FD3diff=0, FD3val=0, Nomval=0;

  for (int i=0; i<h_FD->GetXaxis()->GetNbins(); i++)
    {
      // Divide hist
      bool dummy;
      dummy = h_FDn->Divide(h_FD,Nom_copy);
      dummy = h_FD2n->Divide(h_FD2,Nom_copy);
      dummy = h_FD3n->Divide(h_FD3,Nom_copy);
      dummy = h_Nomn->Divide(h_Nom,Nom_copy);
     
      FDval = h_FDn_diff->GetBinContent(i);
      FD2val = h_FD2n_diff->GetBinContent(i);
      FD3val = h_FD3n_diff->GetBinContent(i);
      FDerror = h_FDn_diff->GetBinError(i);
      FD2error = h_FD2n_diff->GetBinError(i);
      FD3error = h_FD3n_diff->GetBinError(i);
      Nomval = h_Nom->GetBinContent(i);
      FDdiff = fabs(FDval-Nomval);
      FD2diff = fabs(FD2val-Nomval);
      FD3diff = fabs(FD3val-Nomval);

      if (FDerror==0)
	FDerror=0.001;
      if (FD2error==0)
        FD2error=0.001;
      if (FD3error==0)
        FD3error=0.001;
      
      h_FDn_diff->SetBinContent(i,FDdiff/FDerror);
      h_FD2n_diff->SetBinContent(i,FD2diff/FD2error);
      h_FD3n_diff->SetBinContent(i,FD3diff/FD3error);
    }

  fout->cd();

  // Draw un-divided, bin-normalised plot
  TCanvas *c2 = new TCanvas("c2", "c2", 1000, 1000);
  h_FDn_un->GetYaxis()->SetTitle("Events/bin");
  h_FDn_un->GetYaxis()->SetTitleOffset(1.4);
  h_FDn_un->GetXaxis()->SetTitle("Neutrino Energy (GeV)");
  h_FDn_un->SetTitle("SK PP Comparison : " + spectrum);
  h_FDn_un->Draw("p e0");
  h_FD2n_un->Draw("p e0 same");
  //h_FD3n_un->Draw("p e0 same");
  h_Nomn_un->Draw("c same");
  c2->Write("_undiv_norm");
  gPad->SetTickx() ;
  gPad->SetTicky() ;

  // Draw divided plot
  h_FDn->GetXaxis()->SetTitle("Neutrino Energy (GeV)");
  h_FDn->GetYaxis()->SetRangeUser(0,3.9);
  h_FD2n->GetYaxis()->SetRangeUser(0,3.9);
  h_FD3n->GetYaxis()->SetRangeUser(0,3.9);
  h_FDn->GetYaxis()->SetTitle("Events (ratio to nominal)");
  h_FDn->SetTitle("");
  h_FDn->Draw("c");
  h_FD2n->Draw("c same");
  //h_FD3n->Draw("c same");
  h_Nomn->Draw("c same");
  c2->Write("divided");

  // Draw difference plot
  h_FDn_diff->SetLineColor(kRed);
  h_FD2n_diff->SetLineColor(kBlue);
  h_FD3n_diff->SetLineColor(kGreen);
  h_FDn_diff->GetYaxis()->SetRangeUser(0,3.9);
  h_FD2n_diff->GetYaxis()->SetRangeUser(0,3.9);
  h_FD3n_diff->GetYaxis()->SetRangeUser(0,3.9);
  h_FDn_diff->GetYaxis()->SetTitle("#Delta/#delta");
  h_FDn_diff->SetTitle("");
  h_FDn_diff->Draw("c");
  h_FD2n_diff->Draw("c same");
  //h_FD3n_diff->Draw("c same");
  h_Nomn->Draw("c same");
  c2->Write("difference");

  // Make one canvas with both plots
  double uppermin=0.35;
  double middlemin=0.237;
  TPad *lower = new TPad("lower","pad",0,0,1,middlemin);
  TPad *middle = new TPad("middle", "pad",0,middlemin,1,uppermin);
  TPad *upper = new TPad("upper","pad",0,uppermin,1,1);
  upper->SetBottomMargin(0.01);
  middle->SetTopMargin(0.0);
  middle->SetBottomMargin(0.02);
  lower->SetTopMargin(0.01);
  lower->SetBottomMargin(0.58);
  upper->Draw();
  middle->Draw();
  lower->Draw();
  c2->cd();

  // set axis ranges
  // set axis ranges
  double maxy[5] = {30,10,15,3,2};
  double maxx[5] = {2,2,1.3,1.3,1.3};

  // divided
  h_FDn->GetYaxis()->SetRangeUser(0.951,1.049);
  h_FDn->GetXaxis()->SetRangeUser(0,2);
  h_FD2n->GetYaxis()->SetRangeUser(0.951,1.049);
  h_FD2n->GetXaxis()->SetRangeUser(0,2);
  h_FD3n->GetYaxis()->SetRangeUser(0.951,1.049);
  h_FD3n->GetXaxis()->SetRangeUser(0,2);
  // not divided 
  h_FD->GetXaxis()->SetRangeUser(0,maxx[counter]);
  h_FD->GetYaxis()->SetRangeUser(0.1,maxy[counter]);
  h_FD2->GetXaxis()->SetRangeUser(0,maxx[counter]);
  h_FD2->GetYaxis()->SetRangeUser(0.1,maxy[counter]);
  h_FD3->GetXaxis()->SetRangeUser(0,maxx[counter]);
  h_FD3->GetYaxis()->SetRangeUser(0.1,maxy[counter]);
  // Difference
  h_FDn_diff->GetYaxis()->SetRangeUser(0,1.9);
  h_FDn_diff->GetXaxis()->SetRangeUser(0,2);
  h_FD2n_diff->GetYaxis()->SetRangeUser(0,1.9);
  h_FD2n_diff->GetXaxis()->SetRangeUser(0,2);  
  h_FD3n_diff->GetYaxis()->SetRangeUser(0,1.9);
  h_FD3n_diff->GetXaxis()->SetRangeUser(0,2);

  //difference
  lower->cd();
  h_FDn_diff->SetLineWidth(2);
  h_FDn_diff->SetLineColor(kRed);
  h_FD2n_diff->SetLineColor(kBlue);
  h_FD3n_diff->SetLineColor(kGreen);
  h_FDn_diff->SetTitle("");
  h_FDn_diff->GetXaxis()->SetTitle("Neutrino Energy (GeV)");
  h_FDn_diff->GetYaxis()->SetTitle("#Delta/#delta");
  h_FDn_diff->GetYaxis()->CenterTitle();
  h_FDn_diff->GetYaxis()->SetTitleSize(0.16);
  h_FDn_diff->GetYaxis()->SetTitleOffset(0.25);
  h_FDn_diff->GetYaxis()->SetNdivisions(305,kTRUE);
  h_FDn_diff->GetYaxis()->SetLabelSize(0.1);
  h_FDn_diff->GetXaxis()->SetTitleSize(0.17);
  h_FDn_diff->GetXaxis()->SetTitleOffset(1.2);
  h_FDn_diff->GetXaxis()->SetLabelSize(0.13);
  h_FDn_diff->Draw("L hist");
  h_FD2n_diff->Draw("L hist same");
  //h_FD3n_diff->Draw("L hist same");
  h_Nomn->Draw("L hist same");

  // divided
  middle->cd();
  h_FDn->SetTitle("");

  h_FDn->GetYaxis()->SetTitle("Ratio");
  h_FDn->GetYaxis()->CenterTitle();
  h_FDn->GetYaxis()->SetTitleSize(0.31);
  h_FDn->GetYaxis()->SetTitleOffset(0.13);
  h_FDn->GetYaxis()->SetNdivisions(305,kTRUE);
  h_FDn->GetYaxis()->SetLabelSize(0.19);
  h_FDn->GetXaxis()->SetLabelSize(1);
  h_FDn->Draw("L hist");
  h_FD2n->Draw("L hist same");
  //h_FD3n->Draw("L hist same");
  h_Nomn->Draw("L hist same");

  // un-divided, not- bin normalised
  upper->cd();
  upper->SetTitle("SK PP Comparison : " + spectrum);
  h_FD->GetYaxis()->SetTitle("Events");
  h_FD->GetXaxis()->SetLabelOffset(1);
  h_FD->GetYaxis()->SetTitleOffset(0.8);
  h_FD->GetYaxis()->SetTitleSize(0.05);
  h_FD->SetTitle("SK PP Comparison : " + spectrum);
  h_FD2->SetTitle("SK PP Comparison : " + spectrum);
  h_FD3->SetTitle("SK PP Comparison : " + spectrum);
  h_Nom->SetTitle("SK PP Comparison : " + spectrum);
  //  h_FD->GetXaxis()->SetLabelOffset(1);
  h_FD->Draw("p e0");
  h_FD2->Draw("p e0 same");
  //h_FD3->Draw("p e0 same");
  h_Nom->Draw("hist same");

  TPaveText *preliminary = new TPaveText(0.3,0.9,0.7,0.99, "NDC");
  preliminary -> SetBorderSize(0) ;
  preliminary -> SetFillColor(kWhite) ;
  preliminary -> SetTextFont(62) ;
  preliminary -> SetTextColor(kBlack);
  preliminary -> SetTextSizePixels(52) ;
  preliminary -> SetTextSize(0.04);

  TLegend *l = new TLegend(0.5,0.6,0.85,0.87);
  l->AddEntry(h_Nom,nominal,"lp");
  l->AddEntry(h_FD, fd1,"lp");
  l->AddEntry(h_FD2, fd2,"lp");
  //l->AddEntry(h_FD3, fd3,"lp");
  l->SetFillColor(0);
  l->Draw();
  l->Write("legend");

  c2->Write("combined_canvas");

  h_FD->Write("FD");
  h_FD2->Write("FD2");
  h_FD3->Write("FD3");
  h_Nom->Write("Nom");
  h_FDn->Write("FD_divNom");
  h_FD2n->Write("FD2_divNom");
  h_FD3n->Write("FD3_divNom");
  h_FDn_un->Write("FD_norm");
  h_FD2n_un->Write("FD2_norm");
  h_FD3n_un->Write("FD3_norm");  
  h_Nomn->Write("Nom_divNom");
  h_Nomn_un->Write("Nom_norm");
  h_FDn_diff->Write("FD_diff");
  h_FD2n_diff->Write("FD2_diff");
  h_FD3n_diff->Write("FD3_diff");
  std::cout << "Closing " <<  outfile  <<  std::endl;
  fout->Close();
}



// Function to do this for all spectra

void DivideUnosc_all(TString file_osc, TString file_osc2, TString file_osc3, TString file_unosc, TString outfile)
{
  SKSpectra_Divide2Panel4Spectra(file_osc, file_osc2, file_osc3, file_unosc, outfile+"_numu.root", "numu", 0);
  SKSpectra_Divide2Panel4Spectra(file_osc, file_osc2, file_osc3, file_unosc, outfile+"_numubar.root", "numubar",1);
  SKSpectra_Divide2Panel4Spectra(file_osc, file_osc2, file_osc3, file_unosc, outfile+"_nue.root", "nue",2);
  SKSpectra_Divide2Panel4Spectra(file_osc, file_osc2, file_osc3, file_unosc, outfile+"_nuebar.root", "nuebar",3);
  SKSpectra_Divide2Panel4Spectra(file_osc, file_osc2, file_osc3,file_unosc, outfile+"_nue1pi.root", "nue1pi",4);
}
