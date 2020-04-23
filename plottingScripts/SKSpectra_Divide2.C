#include "TArrayD.h"
void SKSpectra_Divide2(char *file_prior, char *file_post, char *outfile)
{

  TGaxis::SetMaxDigits(3) ;
  TString spec[5] = {"numu", "numubar", "nue", "nuebar", "nue1pi"}; 
  TString label[2] = {"Pre ND Fit", "Post ND Fit"};
  TString title[5] = {"1R_{#mu}", "RHC 1R_{#mu}", "1R_{e}", "RHC R_{e}", "1R_{e} 1de"};
  double MAXY[5] = {25,8,12,2,2};
  double MAXX[5] = {2,2,1.3,1.3,1.3};
  
  for(int s=0; s<5; s++){

    TString name = "hist_1d_" + spec[s];

    // Get posterior spectrum
    TFile *f_post = new TFile(file_post, "open");
    TString postname = "posterior_" + name;

    TH1D *h_post = (TH1D*)f_post->Get(name);
    h_post->SetName(postname);
    h_post->Sumw2();

    // Get prior spectrum
    TFile *f_prior = new TFile(file_prior,"open");
    TString priorname = "prior_" + name;

    TH1D *h_prior = (TH1D*)f_prior->Get(name);
    h_prior->SetName(priorname);
    h_prior->Sumw2();
  
    if (h_post->GetXaxis()->GetNbins() != h_prior->GetXaxis()->GetNbins())
      std::cerr << "Error! " << h_post->GetXaxis()->GetNbins() << " bins in posterior histogram and " << h_prior->GetXaxis()->GetNbins() << " bins in prior histogram" << std::endl;

    // Set draw options to be applied to all plots
    h_post->SetLineColor(kBlue);
    h_post->SetLineWidth(2);
    h_post->SetMarkerStyle(20);
    h_prior->SetLineColor(kBlack);
    h_prior->SetLineWidth(2);
    h_prior->SetMarkerStyle(20); 

    h_post->GetXaxis()->SetTitleSize(0.07);
    h_post->GetXaxis()->SetLabelSize(0.06);
    h_post->GetYaxis()->SetTitleSize(0.07);
    h_post->GetYaxis()->SetLabelSize(0.06);

    // copy histogram for dividing
    TH1D *prior_copy = (TH1D*)h_prior->Clone("prior_copy");
    
    // Histograms for dividing 
    TH1D *h_postn = (TH1D*)h_post->Clone("postn");
    TH1D *h_priorn = (TH1D*)h_prior->Clone("priorn");
  
    for (int i=0; i<h_post->GetXaxis()->GetNbins(); i++)
      {
	// Divide hist
	bool dummy;
	dummy = h_postn->Divide(h_post,prior_copy);
	dummy = h_priorn->Divide(h_prior,prior_copy);
      }
    
    TCanvas *c2 = new TCanvas("c2", "c2", 1000, 1000);
    
    // Draw divided plot
    h_postn->GetXaxis()->SetTitle("E_{rec} (GeV)");
    h_postn->GetYaxis()->SetRangeUser(0,3.9);
    h_postn->GetYaxis()->SetTitle("Ratio");
    h_postn->SetTitle("");
    h_postn->Draw("c");
    h_priorn->Draw("c same");
    
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
    h_postn->GetYaxis()->SetRangeUser(0.01,1.99);
    h_postn->GetXaxis()->SetRangeUser(0,2);
    // not divided 
    h_post->GetXaxis()->SetRangeUser(0,2);
    h_post->GetYaxis()->SetRangeUser(0.1,maxy);
   
    // divided
    lower->cd();
    h_postn->SetTitle("");
    
    h_postn->SetTitle("");
    h_postn->GetXaxis()->SetTitle("E_{rec} (GeV)");
    h_postn->GetYaxis()->SetTitle("Ratio");
    h_postn->GetYaxis()->CenterTitle();
    h_postn->GetYaxis()->SetTitleSize(0.16);
    h_postn->GetYaxis()->SetTitleOffset(0.25);
    h_postn->GetYaxis()->SetNdivisions(305,kTRUE);
    h_postn->GetYaxis()->SetLabelSize(0.1);
    h_postn->GetXaxis()->SetTitleSize(0.17);
    h_postn->GetXaxis()->SetTitleOffset(1.2);
    h_postn->GetXaxis()->SetLabelSize(0.13);
    h_postn->SetLineColor(kBlue+2);
    h_postn->Draw("L hist");
   
    TLine *line = new TLine(0, 1.0, MAXX[s], 1.0);
    line->SetLineWidth(2);
    line->SetLineColor(kRed+2);
    line->Draw();
    
    // un-divided, not- bin normalised
    upper->cd();
    //    h_post->SetTitle(title[s]);
    h_post->GetYaxis()->SetTitle("Events/bin");
    h_post->GetYaxis()->SetTitleOffset(0.87);
    h_post->GetYaxis()->SetTitleSize(0.05);
    h_post->GetYaxis()->SetLabelSize(0.04);
    h_post->SetFillStyle(3224);
    h_post->SetFillColor(kBlue+2);
    h_prior->SetFillColor(kRed+2);
    h_prior->SetFillStyle(3003);
    h_prior->SetMarkerStyle(7);
    h_prior->SetMarkerColor(kRed+2);
    h_prior->SetLineColor(kRed+2);
    h_post->SetMarkerStyle(1);
    h_post->SetMarkerColor(kBlue+2);
    h_post->SetLineColor(kBlue+2);
    h_post->Draw("e1");
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
    l->AddEntry(h_post,label[1],"lp");
    l->SetFillColor(0);
    l->Draw();
    c2->Print("skspecleg.pdf");
  }
}
