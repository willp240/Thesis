

TString checkTree(TString filename){

  TFile *f = new TFile(filename, "READ");
  
  TString result = "";

  if(f->GetListOfKeys()->Contains("rap_header"))
    result = "header exists";

  return result;

}
