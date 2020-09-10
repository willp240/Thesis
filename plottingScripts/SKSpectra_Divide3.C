#include "TArrayD.h"
void SKSpectra_Divide3(char *file_prior, char *file_post1, char *file_post2, char *file_post3, char *outfile)
{

  TGaxis::SetMaxDigits(3) ;
  TString spec[5] = {"numu", "numubar", "nue", "nuebar", "nue1pi"};
  TString label[4] = {"Prior", "Uniform Fit Bins, 574 Det. Bins", "Non-Uniform Fit Bins, 574 Det. Bins", "Non-Uniform Fit and Det. Bins"};
  TString title[5] = {"1R_{#mu}", "RHC 1R_{#mu}", "1R_{e}", "RHC R_{e}", "1R_{e} 1de"};
  double MAXY[5] = {25,8,12,2,2};
  double MAXX[5] = {2,2,1.3,1.3,1.3};
  
  double ratiomax[5] = {2.0, 2.0, 1.25, 1.25, 1.25};
  double ratiomin[5] = {0.25, 0.25, 0.1, 0.15, 0.45};

  for(int s=0; s<5; s++){

    TString name = "hist_1d_" + spec[s];

    // Get posterior spectrum
    TFile *f_post1 = new TFile(file_post1, "open");
    TFile *f_post2 = new TFile(file_post2, "open");
    TFile *f_post3 = new TFile(file_post3, "open");
    TString postname1 = "posterior1_" + name;
    TString postname2 = "posterior2_" + name;
    TString postname3 = "posterior3_" + name;

    TH1D *h_post1 = (TH1D*)f_post1->Get(name);
    h_post1->SetName(postname1);
    h_post1->Sumw2();
    TH1D *h_post2 = (TH1D*)f_post2->Get(name);
    h_post2->SetName(postname2);
    h_post2->Sumw2();
    TH1D *h_post3 = (TH1D*)f_post3->Get(name);
    h_post3->SetName(postname3);
    h_post3->Sumw2();

    // Get prior spectrum
    TFile *f_prior = new TFile(file_prior,"open");
    TString priorname = "prior_" + name;

    TH1D *h_prior = (TH1D*)f_prior->Get(name);
    h_prior->SetName(priorname);
    h_prior->Sumw2();
  
    if (h_post1->GetXaxis()->GetNbins() != h_prior->GetXaxis()->GetNbins())
      std::cerr << "Error! " << h_post->GetXaxis()->GetNbins() << " bins in posterior1 histogram and " << h_prior->GetXaxis()->GetNbins() << " bins in prior histogram" << std::endl;
    if (h_post2->GetXaxis()->GetNbins() != h_prior->GetXaxis()->GetNbins())
      std::cerr << "Error! " << h_post->GetXaxis()->GetNbins() << " bins in posterior2 histogram and " << h_prior->GetXaxis()->GetNbins() << " bins in prior histogram" << std::endl;
    if (h_post3->GetXaxis()->GetNbins() != h_prior->GetXaxis()->GetNbins())
      std::cerr << "Error! " << h_post3->GetXaxis()->GetNbins() << " bins in posterior3 histogram and " << h_prior->GetXaxis()->GetNbins() << " bins in prior histogram" << std::endl;

    // Set draw options to be applied to all plots
    h_post1->SetLineColor(kBlue);
    h_post1->SetLineWidth(2);
    h_post1->SetMarkerStyle(20);
    h_post2->SetLineColor(kBlack);
    h_post2->SetLineWidth(2);
    h_post2->SetMarkerStyle(20);
    h_post3->SetLineColor(kGreen);
    h_post3->SetLineWidth(2);
    h_post3->SetMarkerStyle(20);
    h_prior->SetLineColor(kRed);
    h_prior->SetLineWidth(2);
    h_prior->SetMarkerStyle(20); 

    h_post1->GetXaxis()->SetTitleSize(0.07);
    h_post1->GetXaxis()->SetLabelSize(0.06);
    h_post1->GetYaxis()->SetTitleSize(0.07);
    h_post1->GetYaxis()->SetLabelSize(0.06);

    // copy histogram for dividing
    TH1D *prior_copy = (TH1D*)h_prior->Clone("prior_copy");
    
    // Histograms for dividing 
    TH1D *h_post1n = (TH1D*)h_post1->Clone("post1n");
    TH1D *h_post2n = (TH1D*)h_post2->Clone("post2n");
    TH1D *h_post3n = (TH1D*)h_post3->Clone("post3n");
    TH1D *h_priorn = (TH1D*)h_prior->Clone("priorn");
  
    //    for (int i=0; i<h_post1->GetXaxis()->GetNbins(); i++)
    // {
	// Divide hist
    //bool dummy;
    h_post1n->Divide(prior_copy);
    h_post2n->Divide(prior_copy);
    h_post3n->Divide(prior_copy);
    h_priorn->Divide(prior_copy);
    //      }

    TCanvas *c2 = new TCanvas("c2", "c2", 1000, 1000);
    
    // Draw divided plot
    h_post1n->GetXaxis()->SetTitle("E_{rec} (GeV)");
    h_post1n->GetYaxis()->SetRangeUser(0,3.9);
    h_post1n->GetYaxis()->SetTitle("Ratio");
    h_post1n->SetTitle("");
    h_post2n->GetXaxis()->SetTitle("E_{rec} (GeV)");
    h_post2n->GetYaxis()->SetRangeUser(0,3.9);
    h_post2n->GetYaxis()->SetTitle("Ratio");
    h_post2n->SetTitle("");
    h_post3n->GetXaxis()->SetTitle("E_{rec} (GeV)");
    h_post3n->GetYaxis()->SetRangeUser(0,3.9);
    h_post3n->GetYaxis()->SetTitle("Ratio");
    h_post3n->SetTitle("");
    //h_post1n->Draw("L hist");
    //h_post2n->Draw("L hist same");
    //h_post3n->Draw("c same");
    //h_priorn->Draw("c same");

    // Make one canvas with both plots
    double uppermin=0.35;
    double middlemin=0.237;
    TPad *lower = new TPad("lower","pad",0,0,1,middlemin);
    TPad *upper = new TPad("upper","pad",0,middlemin,1,1);
    upper->SetBottomMargin(0.01);
    lower->SetTopMargin(0.01);
    lower->SetBottomMargin(0.58);
    upper->Draw();
    lower->Draw();
    c2->cd();
    
    // set axis ranges
    double maxy = MAXY[s];
    // divided
    h_post1n->GetYaxis()->SetRangeUser(0.01,1.99);
    h_post1n->GetXaxis()->SetRangeUser(0,2);
    h_post2n->GetYaxis()->SetRangeUser(0.01,1.99);
    h_post2n->GetXaxis()->SetRangeUser(0,2);
    h_post3n->GetYaxis()->SetRangeUser(0.01,1.99);
    h_post3n->GetXaxis()->SetRangeUser(0,2);
    // not divided 
    h_post1->GetXaxis()->SetRangeUser(0,2);
    h_post1->GetYaxis()->SetRangeUser(0.1,maxy);
    h_post2->GetXaxis()->SetRangeUser(0,2);
    h_post2->GetYaxis()->SetRangeUser(0.1,maxy);   
    h_post3->GetXaxis()->SetRangeUser(0,2);
    h_post3->GetYaxis()->SetRangeUser(0.1,maxy);

    // divided
    lower->cd();

    for (int index=0; index<h_post1n->GetXaxis()->GetNbins()+2; index++){
      if(h_post1n->GetBinContent(index)>1.5 || h_post1n->GetBinContent(index)<0.5)
        h_post1n->SetBinContent(index,1.0);
      if(h_post2n->GetBinContent(index)>1.5 || h_post2n->GetBinContent(index)<0.5)
        h_post2n->SetBinContent(index,1.0);
      if(h_post3n->GetBinContent(index)>1.5 || h_post3n->GetBinContent(index)<0.5)
        h_post3n->SetBinContent(index,1.0);
    }
    
    h_post1n->SetTitle("");
    h_post1n->GetXaxis()->SetTitle("E_{rec} (GeV)");
    h_post1n->GetYaxis()->SetTitle("Ratio");
    h_post1n->GetYaxis()->CenterTitle();
    h_post1n->GetYaxis()->SetTitleSize(0.16);
    h_post1n->GetYaxis()->SetTitleOffset(0.25);
    h_post1n->GetYaxis()->SetNdivisions(305,kTRUE);
    h_post1n->GetYaxis()->SetLabelSize(0.1);
    h_post1n->GetXaxis()->SetTitleSize(0.17);
    h_post1n->GetXaxis()->SetTitleOffset(1.2);
    h_post1n->GetXaxis()->SetLabelSize(0.13);
    h_post1n->SetLineColor(kBlue);
    h_post2n->SetTitle("");
    h_post2n->GetXaxis()->SetTitle("E_{rec} (GeV)");
    h_post2n->GetYaxis()->SetTitle("Ratio");
    h_post2n->GetYaxis()->CenterTitle();
    h_post2n->GetYaxis()->SetTitleSize(0.16);
    h_post2n->GetYaxis()->SetTitleOffset(0.25);
    h_post2n->GetYaxis()->SetNdivisions(305,kTRUE);
    h_post2n->GetYaxis()->SetLabelSize(0.1);
    h_post2n->GetXaxis()->SetTitleSize(0.17);
    h_post2n->GetXaxis()->SetTitleOffset(1.2);
    h_post2n->GetXaxis()->SetLabelSize(0.13);
    h_post2n->SetLineColor(kBlack);
    h_post3n->SetTitle("");
    h_post3n->GetXaxis()->SetTitle("E_{rec} (GeV)");
    h_post3n->GetYaxis()->SetTitle("Ratio");
    h_post3n->GetYaxis()->CenterTitle();
    h_post3n->GetYaxis()->SetTitleSize(0.16);
    h_post3n->GetYaxis()->SetTitleOffset(0.25);
    h_post3n->GetYaxis()->SetNdivisions(305,kTRUE);
    h_post3n->GetYaxis()->SetLabelSize(0.1);
    h_post3n->GetXaxis()->SetTitleSize(0.17);
    h_post3n->GetXaxis()->SetTitleOffset(1.2);
    h_post3n->GetXaxis()->SetLabelSize(0.13);
    h_post3n->SetLineColor(kGreen);
    //    h_post1n->GetXaxis()->SetRangeUser(ratiomin[s],ratiomax[s]);
    //    h_post2n->GetXaxis()->SetRangeUser(ratiomin[s],ratiomax[s]);
    //    h_post3n->GetXaxis()->SetRangeUser(ratiomin[s],ratiomax[s]);
    h_post1n->Draw("L hist");
    h_post2n->Draw("L hist same");
    h_post3n->Draw("L hist same");    

    TLine *line = new TLine(0, 1.0, MAXX[s], 1.0);
    line->SetLineWidth(2);
    line->SetLineColor(kRed+2);
    line->Draw();
    
    // un-divided, not- bin normalised
    upper->cd();
    h_post1->SetTitle("");
    h_post1->GetYaxis()->SetTitle("Events/bin");
    h_post1->GetYaxis()->SetTitleOffset(0.87);
    h_post1->GetYaxis()->SetTitleSize(0.05);
    h_post1->GetYaxis()->SetLabelSize(0.04);
    h_post1->SetFillStyle(3224);
    h_post1->SetFillColor(kBlue);
    h_post2->SetTitle("");
    h_post2->GetYaxis()->SetTitle("Events/bin");
    h_post2->GetYaxis()->SetTitleOffset(0.87);
    h_post2->GetYaxis()->SetTitleSize(0.05);
    h_post2->GetYaxis()->SetLabelSize(0.04);
    h_post2->SetFillStyle(3224);
    h_post2->SetFillColor(kBlack);
    h_post3->SetTitle("");
    h_post3->GetYaxis()->SetTitle("Events/bin");
    h_post3->GetYaxis()->SetTitleOffset(0.87);
    h_post3->GetYaxis()->SetTitleSize(0.05);
    h_post3->GetYaxis()->SetLabelSize(0.04);
    h_post3->SetFillStyle(3224);
    h_post3->SetFillColor(kGreen);
    h_prior->SetFillColor(kRed);
    h_prior->SetFillStyle(3003);
    h_prior->SetMarkerStyle(7);
    h_prior->SetMarkerColor(kRed);
    h_prior->SetLineColor(kRed);
    h_post1->SetMarkerStyle(1);
    h_post1->SetMarkerColor(kBlue);
    h_post1->SetLineColor(kBlue);
    h_post2->SetMarkerStyle(1);
    h_post2->SetMarkerColor(kBlack);
    h_post2->SetLineColor(kBlack);
    h_post3->SetMarkerStyle(1);
    h_post3->SetMarkerColor(kGreen);
    h_post3->SetLineColor(kGreen);
    h_post1->Draw("e1");
    h_post2->Draw("e1 same");
    h_post3->Draw("e1 same");
    h_prior->Draw("p e5 same");
    
    TPaveText *preliminary = new TPaveText(0.3,0.9,0.7,0.99, "NDC");
    preliminary -> SetBorderSize(0) ;
    preliminary -> SetFillColor(kWhite) ;
    preliminary -> SetTextFont(62) ;
    preliminary -> SetTextColor(kBlack);
    preliminary -> SetTextSizePixels(52) ;
    preliminary -> SetTextSize(0.04);

    c2->Print(outfile+spec[s]+".pdf");
    c2->Clear();
    TLegend *l = new TLegend(0.,0.,1.0,1.0);
    l->AddEntry(h_prior,label[0],"lpf");
    l->AddEntry(h_post1,label[1],"lp");
    l->AddEntry(h_post2,label[2],"lp");
    l->AddEntry(h_post3,label[3],"lp");
    l->SetFillColor(0);
    l->Draw();
    c2->Print("polyskspecleg.pdf");
  }
}
