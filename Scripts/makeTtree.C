

void makeTtree(){

  TFile *file1 = TFile::Open("testIsohalf.root");
 
  TTree *ssum = (TTree*)(file1->Get("sample_sum"))->Clone();
  int sample;
  double x,y;
  TObjArray* gMAQE;
  TObjArray* gDIS;
  TGraph* g1;
  ssum->SetBranchAddress("BgSclLMCPiResGraph",&gDIS);
  ssum->SetBranchAddress("MAQEGraph",&gMAQE);
  ssum->SetBranchAddress("SelectedSample",&sample);
  for(int i=0; i<ssum->GetEntries(); i++)
    {
      x=0, y=0;
      ssum->GetEntry(i);
      g1 = (TGraph*)(gDIS->At(sample)->Clone());
      g1->GetPoint(1,x,y);
      if(y>0.0)
	std::cout << "Entry " << i << ", Sample: " << sample << std::endl; 
    }
  ssum->GetEntry(4);
  g1 = (TGraph*)(gDIS->At(sample)->Clone());
  g1->SetTitle("Run3b_000 Event 4 Sample 5");
  g1->GetXaxis()->SetTitle("kNIWG_HadMult_AGKY");
  g1->GetYaxis()->SetTitle("Weight");
  g1->Draw();
}
