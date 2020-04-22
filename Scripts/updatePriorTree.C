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

void updatePriorTree(std::string file1, std::string file2, std::string file3 ){

  TFile *f1= new TFile(file1.c_str(),"READ");
  TFile *f2= new TFile(file2.c_str(),"READ");
  TFile *f3= new TFile(file3.c_str(),"NEW");
  
  TTree *pos = (TTree*)f1->Get("posteriors");
  TTree *set = (TTree*)f2->Get("Settings");

  f3->cd();
  pos->CloneTree()->Write();
  set->CloneTree()->Write();
  
}
