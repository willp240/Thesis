// Chain reweighting script, should work for both reduced/non-reduced chain,
// and for both flat and gaussian priors, for any parameter you wish. This
// script will not overwrite your precious chains, it will make a new file with
// a "_reweighted.root" postfix.
// USAGE:
//
// updatePrior("chain.root") 
//   reweight sinth13 prior from flat to PDG 2018
//
// updatePrior("chain.root", "parameter", central_new, error_new) 
//   reweight "parameter" prior from flat to a gaussian centered around
//   central_new with error_new uncertainty
//
// updatePrior("chain.root", "parameter", central_new, error_new, central_old, error_old)
//   reweight "parameter" prior from gaussian centered around central_old with
//   error_old uncertainty to a gaussian centered around central_new with
//   error_new uncertainty

// TODO:
//  * More cerr if e.g. branch not loaded (wrong branch name given)
//  * Reweight from gauss prior to flat (if you need it)
//  * Sort out some memory issues, I'm sure it is possible have only one TFile
//    and TTree at the time... Probably because of CloneTree()?

void updatePrior(std::string sin, std::string param="theta13", double new_central = 0.0212, double new_error = 0.0008,  double old_central = -999, double old_error = -999 ){

  // First get the in/out name
  std::size_t pos;
  pos = sin.find(".root");
  // This will be your file output name
  std::string sout = sin.substr(0,pos) + "_reweighted.root";

  // Weights
  Double_t reweight;
  // Parameter value
  Double_t parameter;

  // Load the file in
  TFile *fin= new TFile(sin.c_str(),"READ");

  // Load the ttree
  TTree *treeold; 
  TTree *tree;
  if (fin->Get("osc_posteriors")){
    std::cout << "Loading a non-reduced chain (osc_posteriors)" << std::endl;
    treeold = (TTree*)fin->Get("osc_posteriors");
  }
  else if (fin->Get("posteriors")) {
    std::cout << "Loading a reduced chain (posteriors)" << std::endl;
    treeold = (TTree*)fin->Get("posteriors");
  }
  else {
    std::cerr << "Make sure you're loading a chain from MaCh3..." << std::endl;
    std::exit(1);
  }

  // Make a new file and clone the old TTree
  TFile *fout= new TFile(sout.c_str(),"recreate");
  tree = treeold ->CloneTree();
  fin->Close();
  delete fin;
  // Can't delete treeold, it crashes :/ Is there's somethig better than CloneTree()?

  // Set the new branch containing our weights
  TBranch *brout = tree->Branch("reweight", &reweight);

  // Load the parameter that we want to reweight the prior for
  tree->SetBranchAddress(param.c_str(), &parameter);
  double old_chi;
  double old_prior;
  double new_chi;
  double new_prior;

  // Now let's iterate through the chain and reweight
  Long64_t nentries = tree->GetEntries();
  for (Long64_t i = 0; i < nentries; ++i) {
    tree->GetEntry(i);

    // Calculate contribution from the new prior
    new_chi = (parameter - new_central)/new_error;
    new_prior = std::exp(-0.5 * new_chi * new_chi );

    // Calculate contribution from the old prior (1 if flat)
    if(old_central == -999)
      old_prior = 1.0;
    else{
      old_chi = (parameter - old_central)/old_error;
      old_prior = std::exp(-0.5 * old_chi * old_chi);
    }

    // Filling weight branch
    reweight = new_prior/old_prior;
    brout->Fill();
  }

  // Save. The end.
  tree->Write();
  fout->Close();
}
