TFile *File = NULL;
TCanvas *c = new TCanvas("canv", "canv", 1080, 1080);
std::string filename = "";

void LoopDirectory(TDirectory *dir) {
  TIter next(dir->GetListOfKeys());
  TKey *key;
  std::string name;

  // Loop through all entries
  while ((key = (TKey*)next())) {

    std::string name = std::string(key->GetName());
    std::string classname = std::string(key->GetClassName());
    if (classname == "TDirectoryFile") {
      dir->cd(name.c_str());
      TDirectory *SubDir = gDirectory;
      LoopDirectory(SubDir);
    } 
    std::string dirname = std::string(dir->GetName());
    if (dirname.find(".root") != std::string::npos)
      ;
    else
      name = dirname + "/" + name;
    // Select only data plots
    if (classname == "TH2Poly") {
      TH2Poly *plot = File->Get(name.c_str())->Clone();
      plot->GetXaxis()->SetRangeUser(0,5000);
      c->cd();
      c->SetTopMargin(0.06);
      c->SetRightMargin(0.18);
      c->SetLeftMargin(0.13);
      c->SetBottomMargin(0.06)
      plot->Draw("colz");
      c->Print((filename+".pdf").c_str());
      delete plot;
    }

  }
}

void plotPolyMode(std::string FileName) {
  File = new TFile(FileName.c_str(), "OPEN");
  filename = FileName;

  c->SetTopMargin(0.07);
  c->SetBottomMargin(0.08);
  c->SetLeftMargin(0.10);
  c->SetRightMargin(0.03);

  c->Print((filename+".pdf[").c_str());
  c->cd();

  TPaveText *pt = new TPaveText(0.05, 0.05, 0.9, 0.9);
  pt->AddText(filename.c_str());
  pt->Draw();
  c->Print((filename+".pdf").c_str());

  LoopDirectory(File);

  c->Print((filename+".pdf]").c_str());
}

