
void plotLLHContAcc(std::string fileA, std::string fileD ){

  TCanvas *c = new TCanvas("canv", "canv", 900, 600);
  c->SetTopMargin(0.07);
  c->SetBottomMargin(0.11);
  c->SetLeftMargin(0.18);
  c->SetRightMargin(0.1);

  TFile *fA= new TFile(fileA.c_str(),"READ");
  TFile *fD= new TFile(fileD.c_str(),"READ");

  TH1D* hAx = new TH1D("hAx","hAx",100,0,100);
  TH1D* hDx = new TH1D("hDx","hDx",100,0,100);
  TH1D* hAb = new TH1D("hAb","hAb",100,0,100);
  TH1D* hDb = new TH1D("hDb","hDb",100,0,100);
  TH1D* hAnd = new TH1D("hAnd","hAnd",1000,0,1000);
  TH1D* hDnd = new TH1D("hDnd","hDnd",1000,0,1000);
  TH1D* hAs = new TH1D("hAs","hAs",5000,0,5000);
  TH1D* hDs = new TH1D("hDs","hDs",5000,0,5000);  

  fA->cd();
  TTree *posA = (TTree*)fA->Get("posteriors");
  fD->cd();
  TTree *posD = (TTree*)fD->Get("posteriors");

  double xsLLHA=0, xsLLHD=0, xs_0A=0, xs_0D=0, ndLLHA=0, ndLLHD=0, bLLHA=0, bLLHD=0, sLLHA=0, sLLHD=0;
  double prevxs_0A=-1, prevxs_0D=-1, prevxs_LLHA=-1, prevxs_LLHD=-1, prevnd_LLHA=-1, prevnd_LLHD=-1, prevb_LLHA=-1, prevb_LLHD=-1, prevs_LLHA=-1, prevs_LLHD=-1;
  int stepA, stepD;

  posA->SetBranchAddress("LogL_systematic_xsec_cov",&xsLLHA);
  posA->SetBranchAddress("xsec_0",&xs_0A);
  posD->SetBranchAddress("LogL_systematic_xsec_cov",&xsLLHD);
  posD->SetBranchAddress("xsec_0",&xs_0D);
  posA->SetBranchAddress("LogL_systematic_total_flux_cov",&bLLHA);
  posD->SetBranchAddress("LogL_systematic_total_flux_cov",&bLLHD);
  posA->SetBranchAddress("LogL_systematic_nddet_cov",&ndLLHA);
  posD->SetBranchAddress("LogL_systematic_nddet_cov",&ndLLHD);
  posA->SetBranchAddress("LogL_sample_0",&sLLHA);
  posD->SetBranchAddress("LogL_sample_0",&sLLHD);
  posA->SetBranchAddress("step", &stepA);
  posD->SetBranchAddress("step", &stepD);

  for(int i=0; i<posA->GetEntries(); i++){

    posA->GetEntry(i);
    if(stepA<200000) continue;

    //if xs_0A hasn't changed step was rejected
    if(xs_0A==prevxs_0A){
      hAx->Fill(prevxs_LLHA);
      hAb->Fill(prevb_LLHA);
      hAnd->Fill(prevnd_LLHA);
      hAs->Fill(prevs_LLHA);
    }
    else{
      hAx->Fill(xsLLHA);
      prevxs_LLHA=xsLLHA;
      hAb->Fill(bLLHA);
      prevb_LLHA=bLLHA;
      hAnd->Fill(ndLLHA);
      prevnd_LLHA=ndLLHA;
      hAs->Fill(sLLHA);
      prevs_LLHA=sLLHA;
      prevxs_0A=xs_0A;
    }
  }

  std::cout << "Done Asimov, now data " << std::endl;

  for(int i=0; i<posD->GetEntries(); i++){

    posD->GetEntry(i);
    if(stepD<200000) continue;

    if(xs_0D==prevxs_0D){
      hDx->Fill(prevxs_LLHD);
      hDb->Fill(prevb_LLHD);
      hDnd->Fill(prevnd_LLHD);
      hDs->Fill(prevs_LLHD);
    }
    else{
      hDx->Fill(xsLLHD);
      prevxs_LLHD=xsLLHD;
      hDb->Fill(bLLHD);
      prevb_LLHD=bLLHD;
      hDnd->Fill(ndLLHD);
      prevnd_LLHD=ndLLHD;
      hDs->Fill(sLLHD);
      prevs_LLHD=sLLHD;
      prevxs_0D=xs_0D;
    }
  }

  TFile *outfile = new TFile("fullout.root", "RECREATE");

  hAx->SetLineColor(kRed+2);
  hDx->SetLineColor(kBlue+2);
  hAb->SetLineColor(kRed+2);
  hDb->SetLineColor(kBlue+2);
  hAnd->SetLineColor(kRed+2);
  hDnd->SetLineColor(kBlue+2);
  hAs->SetLineColor(kRed+2);
  hDs->SetLineColor(kBlue+2);

  hAx->Scale(1/(hAx->Integral()));
  hDx->Scale(1/(hDx->Integral()));
  hAnd->Scale(1/(hAnd->Integral()));
  hDnd->Scale(1/(hDnd->Integral()));
  hAb->Scale(1/(hAb->Integral()));
  hDb->Scale(1/(hDb->Integral()));
  hAs->Scale(1/(hAs->Integral()));
  hDs->Scale(1/(hDs->Integral()));

  hAx->GetXaxis()->SetTitle("LogL_systematic_xsec_cov");
  hAx->GetYaxis()->SetTitle("Steps");
  hAs->GetXaxis()->SetTitle("LogL_sample_0");
  hAs->GetYaxis()->SetTitle("Steps");
  hAnd->GetXaxis()->SetTitle("LogL_systematic_nddet_cov");
  hAnd->GetYaxis()->SetTitle("Steps");
  hAb->GetXaxis()->SetTitle("LogL_systematic_flux_cov");
  hAb->GetYaxis()->SetTitle("Steps");

  hAx->SetTitle("LogL_systematic_xsec_cov");
  hAx->Draw();
  hDx->Draw("same");
  hAx->Fit("gaus");
  hDx->Fit("gaus");
  TLegend* legxs = new TLegend(0.5,0.6,0.9,0.9);
  legxs->AddEntry(hAx,Form("Asimov %.2f +/- %.2f", hAx->GetFunction("gaus")->GetParameter(1), hAx->GetFunction("gaus")->GetParameter(2)),"l");
  legxs->AddEntry(hDx,Form("Data %.2f +/- %.2f", hDx->GetFunction("gaus")->GetParameter(1), hDx->GetFunction("gaus")->GetParameter(2)),"l");
  legxs->Draw();
  c->Write("xscanv");
  hAx->Write();
  hDx->Write();

  hAnd->SetTitle("LogL_systematic_nddet_cov");
  hAnd->Draw();
  hDnd->Draw("same");
  hAnd->Fit("gaus");
  hDnd->Fit("gaus");
  hAnd->GetXaxis()->SetRangeUser(0,750);
  TLegend* legnd = new TLegend(0.5,0.6,0.9,0.9);
  legnd->AddEntry(hAnd,Form("Asimov %.2f +/- %.2f", hAnd->GetFunction("gaus")->GetParameter(1), hAnd->GetFunction("gaus")->GetParameter(2)),"l");
  legnd->AddEntry(hDnd,Form("Data %.2f +/- %.2f", hDnd->GetFunction("gaus")->GetParameter(1), hDnd->GetFunction("gaus")->GetParameter(2)),"l");
  legnd->Draw();
  c->Write("ndcanv");
  hAnd->Write();
  hDnd->Write();

  hAb->Draw();
  hDb->Draw("same");
  hAb->Fit("gaus");
  hDb->Fit("gaus");
  TLegend* legb = new TLegend(0.5,0.6,0.9,0.9);
  legb->AddEntry(hAb,Form("Asimov %.2f +/- %.2f", hAb->GetFunction("gaus")->GetParameter(1), hAb->GetFunction("gaus")->GetParameter(2)),"l");
  legb->AddEntry(hDb,Form("Data %.2f +/- %.2f", hDb->GetFunction("gaus")->GetParameter(1), hDb->GetFunction("gaus")->GetParameter(2)),"l");
  legb->Draw();
  c->Write("bcanv");
  hAb->Write();
  hDb->Write();

  hAs->SetTitle("LogL_sample_0");
  hAs->Draw();
  hDs->Draw("same");
  hAs->Fit("gaus");
  hDs->Fit("gaus");
  hAs->GetXaxis()->SetRangeUser(0,3000);
  TLegend* legs = new TLegend(0.3,0.6,0.7,0.9);
  legs->AddEntry(hAs,Form("Asimov %.2f +/- %.2f", hAs->GetFunction("gaus")->GetParameter(1), hAs->GetFunction("gaus")->GetParameter(2)),"l");
  legs->AddEntry(hDs,Form("Data %.2f +/- %.2f", hDs->GetFunction("gaus")->GetParameter(1), hDs->GetFunction("gaus")->GetParameter(2)),"l");
  hAs->SetTitle("LogL_sample_0");
  legs->Draw();
  c->Write("scanv");
  hAs->Write();
  hDs->Write();

}
