TFile *File = NULL;
TCanvas *c = new TCanvas("canv", "canv");
std::string filename = "";

void LoopDirectory(TDirectory *dir) {
  TIter next(dir->GetListOfKeys());
  TKey *key;
  std::string name;
  int count = 0;
  
  //  TH2Poly *p1sigma;
  //TH2Poly *m1sigma;
  //TH2Poly *nominal;

  // Loop through all entries
  while ((key = (TKey*)next())) {

    std::string name = std::string(key->GetName());
    std::string classname = std::string(key->GetClassName());
    if (classname == "TDirectoryFile") {
      count = 0;
      dir->cd(name.c_str());
      TDirectory *SubDir = gDirectory;
      LoopDirectory(SubDir);
    } 
    std::string dirname = std::string(dir->GetName());

    if (dirname.find(".root") != std::string::npos) continue;

    name = dirname + "/" + name;
    if (classname == "TH2Poly") {
      count++;
      if (count==1){
	TH2Poly *m1sigma = File->Get(name.c_str())->Clone();
	m1sigma->SetName((name+"_m1").c_str());
        m1sigma->SetTitle((name+" Minus 1Sigma Ratio to Nominal").c_str());
      }
      if (count==2)
	{
	  TH2Poly *nominal = File->Get(name.c_str())->Clone();
	  nominal->SetName((name+"_Nominal").c_str());
	  nominal->SetTitle((name+" Nominal").c_str());
	}
      if (count==3){
        TH2Poly *p1sigma = File->Get(name.c_str())->Clone();
	p1sigma->SetName((name+"_p1").c_str());
	p1sigma->SetTitle((name+" Plus 1Sigma Ratio to Nominal").c_str());
      
	c->cd();
	nominal->Draw("colz");
	nominal->GetXaxis()->SetRangeUser(0,5000);
	c->Print((filename+".pdf").c_str());
	m1sigma->Divide(nominal);
	m1sigma->GetZaxis()->SetRangeUser(0.9,1.1);
	m1sigma->GetXaxis()->SetRangeUser(0,5000);
	m1sigma->Draw("colz");
	c->Print((filename+".pdf").c_str());
	p1sigma->Divide(nominal);
	p1sigma->GetZaxis()->SetRangeUser(0.9,1.1);
	p1sigma->GetXaxis()->SetRangeUser(0,5000);
	p1sigma->Draw("colz");
	c->Print((filename+".pdf").c_str());
	 }
    }

  }
}

void sigmaVarRatiosPoly(std::string FileName) {
  File = new TFile(FileName.c_str(), "OPEN");
  filename = FileName;

  c->SetTopMargin(0.1);
  c->SetBottomMargin(0.13);
  c->SetLeftMargin(0.13);
  c->SetRightMargin(0.19);

  c->Print((filename+".pdf[").c_str());
  c->cd();

  TPaveText *pt = new TPaveText(0.05, 0.05, 0.9, 0.9);
  pt->AddText(filename.c_str());
  pt->Draw();
  c->Print((filename+".pdf").c_str());
  LoopDirectory(File);

  c->Print((filename+".pdf]").c_str());
}

