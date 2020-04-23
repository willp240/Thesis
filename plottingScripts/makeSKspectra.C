void makeSKspectra_sample(char *infile = "sk_nominal.root", TString name = "hist", TString outname = "output.root", int maxy = 1, int nthrows = 1, bool do_var = false, TString var ="")
{

  TFile* file = new TFile(infile, "open");
  std::cout << "infile is " << infile << std::endl;
  int nbins = maxy*10000;
  //std::cout << "Setting nbinsy to be " << nbins << std::endl;
  //std::cout << "Max bins is " << maxy << std::endl;
  TString tempname = name+"0";
  //std::cout << "Got temp name as " << tempname << std::endl;
  //std::cout << "Got name as " << name << std::endl;
  TH1D* temp = (TH1D*)file->Get(tempname);
  TH1D* nominal = (TH1D*)file->Get(name);

  int nbinsx = temp->GetNbinsX();
  double binedges[100];

  std::cout << "BLAH" << std::endl;
  printf("%i \n",nbinsx);

  for(int i=0; i<=nbinsx; i++){binedges[i]=temp->GetXaxis()->GetBinLowEdge(i+1);}

  std::cout << "Made Th2D spectra" << std::endl;
  TH2D* spectra = new TH2D("test","test",nbinsx,binedges,nbins+1,0,maxy); //400,0,/*20*/24/*100*/); //500,0,10
  std::cout << "Made Th2D spectra" << std::endl;

  double average=0;

  // I don't remember what these histograms were for or why they were useful,
  // and I never use them anyway, so commented out
  //TH1D* hAvg = new TH1D("hAvg","hAvg",200,10,30);

  int nsteps=nthrows -1;
  std::cout << "Number of throws is " << nthrows << std::endl;

  TRandom3* rnd = new TRandom3(0);
  for(int i=0; i<nsteps; i++)
  {
	TString name1=name;
	name1+=i;

	temp = (TH1D*)file->Get(name1);
	average+=temp->Integral();
	//hAvg->Fill(temp->Integral());

	for(int j=1; j<=temp->GetNbinsX(); j++)
	  spectra->Fill(temp->GetXaxis()->GetBinCenter(j),temp->GetBinContent(j));
  }
  std::cout << "Filled spectra " << std::endl;

  average/=(double)nsteps;
  //printf("%f %f %f\n",average,hAvg->GetMean(),hAvg->GetRMS());
  printf("%f \n",average);


  TCanvas* c = new TCanvas("c","c",0,0,700,700);
  spectra->Draw("colz");
  c->Update();
  //hAvg->Draw();
  //c->Update();

  TFile* outfile = new TFile(outname+var+".root","UPDATE");


  double x[100], y[100], ex[100], ey[100];
  c->Clear();
  TH1D* proj=0;
  for(int i=1; i<=nbinsx; i++)
  {
	TString namepy=name+"_py";
	namepy+=i;
	proj=spectra->ProjectionY(namepy,i,i);

	if(proj->GetEntries() > 0){

	  int first_bin = 99999;
	  int last_bin = 0;

	  for(int bin_i = 1 ; bin_i < proj->GetNbinsX() - 1 ; bin_i++){
		double val = proj->GetBinContent(bin_i);
		if(val > 0 && bin_i < first_bin){
		  first_bin = bin_i;
		}

		if(val > 0){
		  last_bin = bin_i;
		}
	  }

	  //Get the range the filled bins cover
	  double diff = proj->GetBinCenter(last_bin) - proj->GetBinCenter(first_bin); 
	  double opt_bin_width = diff/25;
	  double current_bin_width = proj->GetBinWidth(last_bin);
	  int rebin = 1;
	  if(diff > 0){rebin= (opt_bin_width/current_bin_width);}
	  //std::cout << "rebin is " << rebin << std::endl;
	  //std::cout << "diff is " << diff << std::endl;
	  if(rebin==0){rebin+=1;} 
	  while((proj->GetNbinsX() % rebin) != 0){
		rebin--;
	  }

	  if(rebin > 0){
		proj->Rebin(rebin);
	  }

	  if(diff / proj->GetBinWidth(last_bin) < 50){
   	/*	std::cout << "This hasn't worked!! " << std::endl;
		std::cout << "Current width is " << current_bin_width << std::endl;
		std::cout << "Omptimal width is " << opt_bin_width << std::endl;
		std::cout << "Bin width now is " << proj->GetBinWidth(last_bin) << std::endl;
		std::cout << "Rebin is " << rebin << std::endl;
		std::cout << "First bin is " << first_bin << std::endl;
		std::cout << "Last bin is " << last_bin << std::endl;
	*/	
	  }

	  //Get the bin centre array
	  x[i-1]=spectra->GetXaxis()->GetBinCenter(i);
	  //Set the error on as being half the bin width
	  ex[i-1]=spectra->GetXaxis()->GetBinWidth(i)/2.0;
	  //Set the centre of the bin to be 
	  //y[i-1]=proj->GetBinCenter(proj->GetMaximumBin());
	  y[i-1]=proj->GetMean();
	  //Set the error on y to be RMS of the bin
	  ey[i-1]=proj->GetRMS();
	  //TFitResultPtr r = proj->Fit("gaus","QS","goff",y[i-1]-ey[i-1],y[i-1]+ey[i-1]);
	  //if(y[i-1]>0.1 && r->Parameter(1)>0){y[i-1]=r->Parameter(1);}
	  //      c->Update();
	}
	proj->Write(namepy);


  }

  TGraphErrors* errorbars = new TGraphErrors(nbinsx,x,y,ex,ey);
  errorbars->Draw("A E2");

  TString dummyname;
  if(do_var){
	if (name=="sk_nue_"+var){dummyname = "nue";}
	if (name=="sk_nuebar_"+var){dummyname = "nuebar";}
	if (name=="sk_nue1pi_"+var){dummyname = "nue1pi";}
	if (name=="sk_numu_"+var){dummyname = "numu";}
	if (name=="sk_numubar_"+var){dummyname = "numubar";}
  }
  else{
	if (name=="sk_nue"){dummyname = "nue";}
	if (name=="sk_nuebar"){dummyname = "nuebar";}
	if (name=="sk_nue1pi"){dummyname = "nue1pi";}
	if (name=="sk_numu"){dummyname = "numu";}
	if (name=="sk_numubar"){dummyname = "numubar";}
  }

  // make TH1D with y-errors only
  TString histname = "hist_1d_"+dummyname;
  TH1D *hist1d = (TH1D*)temp->Clone(histname);
  hist1d->Reset();
  for (int i=0; i<hist1d->GetXaxis()->GetNbins(); i++)
  {
	hist1d->SetBinContent(i+1,y[i]);
	hist1d->SetBinError(i+1,ey[i]);
  }

  errorbars->Write(name);
  histname="spectra2d_"+dummyname;
  spectra->Write(histname);
  histname="hist_1d_"+dummyname;
  hist1d->Write(histname);
  histname=name+"_nominal";
  nominal->Write(histname);

  outfile->Close();
}

// -------------------------------------------------------------------------------- //



// -------------------------------------------------------------------------------- //

void makeSKspectra(char *infile="sk_nominal.root", TString outname = "output.root", int nthrows = 1, int do_var = -1)
{
  bool do_var_bool = false;
  TString var;
  switch(do_var){
  case 0:
	do_var_bool = true;
	var = "plep";
	break;	
  case 1:
	do_var_bool = true;
	var = "costheta_lep";
	break;
  case 2:
	do_var_bool = true;
	var = "trueEnu";
	break;
  case 3:
	do_var_bool = true;
	var = "trueQ2";
	break;
  case 4:
	do_var_bool = true;
	var = "ErecQE";
	break;
  case 5:
	do_var_bool = true;
	var = "Q2QE";
	break;
  default:
	break;
  }

  if(do_var_bool){
	std::cout << std::endl << " ---------------- numu ----------------" <<std::endl << std::endl;
	makeSKspectra_sample(infile, "sk_numu_"+var, outname, 100, nthrows, do_var, var);
	std::cout << std::endl << " ---------------- nue ----------------" <<std::endl << std::endl;
	makeSKspectra_sample(infile, "sk_nue_"+var, outname, 25, nthrows, do_var, var);
	std::cout << std::endl << " ---------------- nue1pi ----------------" <<std::endl << std::endl;
	makeSKspectra_sample(infile, "sk_nue1pi_"+var, outname, 5, nthrows, do_var, var);
	std::cout << std::endl << " ---------------- numubar ----------------" <<std::endl << std::endl;
	makeSKspectra_sample(infile, "sk_numubar_"+var, outname, 100, nthrows, do_var, var);
	std::cout << std::endl << " ---------------- nuebar ----------------" <<std::endl << std::endl;
	makeSKspectra_sample(infile, "sk_nuebar_"+var, outname, 10, nthrows, do_var, var);
  }
  else{
	std::cout << std::endl << " ---------------- numu ----------------" <<std::endl << std::endl;
	makeSKspectra_sample(infile, "sk_numu", outname, 25, nthrows);
	std::cout << std::endl << " ---------------- nue ----------------" <<std::endl << std::endl;
	makeSKspectra_sample(infile, "sk_nue", outname, 12, nthrows);
	std::cout << std::endl << " ---------------- nue1pi ----------------" <<std::endl << std::endl;
	makeSKspectra_sample(infile, "sk_nue1pi", outname, 2, nthrows);
	std::cout << std::endl << " ---------------- numubar ----------------" <<std::endl << std::endl;
	makeSKspectra_sample(infile, "sk_numubar", outname, 25, nthrows);
	std::cout << std::endl << " ---------------- nuebar ----------------" <<std::endl << std::endl;
	makeSKspectra_sample(infile, "sk_nuebar", outname, 5, nthrows);
  }
}

// -------------------------------------------------------------------------------- //
void makeSKspectra_beta_sample(char *infile="sk_nominal_beta.root", TString name = "hist", int nbinsy = 1, int maxy = 1)
{

  TFile* file = new TFile(infile, "open");

  TString tmpname = name+"0";
  TString tmpname_beta0 = name+"_beta0_0";
  TString tmpname_beta1 = name+"_beta1_1";
  TH1D* temp = (TH1D*)file->Get(tmpname);
  TH1D* temp_beta0 = (TH1D*)file->Get(tmpname_beta0);
  TH1D* temp_beta1 = (TH1D*)file->Get(tmpname_beta1);
  int nbinsx = temp->GetNbinsX(); // Assume binning is the same for all three
  const int n = nbinsx;
  double binedges[n+1];

  printf("%i \n",nbinsx);

  for(int i=0; i<=nbinsx; i++)
  {
	binedges[i]=temp->GetXaxis()->GetBinLowEdge(i+1);
  }
  binedges[nbinsx+1] = temp->GetXaxis()->GetBinUpEdge(nbinsx);

  TH2D* spectra = new TH2D("spectra","spectra",nbinsx,binedges,nbinsy,0,maxy);
  TH2D* spectra_beta0 = new TH2D("spectra_beta0","spectra_beta0",nbinsx,binedges,nbinsy,0,maxy);
  TH2D* spectra_beta1 = new TH2D("spectra_beta1","spectra_beta1",nbinsx,binedges,nbinsy,0,maxy);

  double average=0, average_beta0=0, average_beta1=0;

  // I don't remember what these histograms were for or why they were useful,
  // and I never use them anyway, so commented out
  //TH1D* hAvg = new TH1D("hAvg","hAvg",200,0,10);
  //TH1D* hAvg_beta0 = new TH1D("hAvg_beta0","hAvg_beta0",200,0,10);
  //TH1D* hAvg_beta1 = new TH1D("hAvg_beta1","hAvg_beta1",200,0,10);
  TRandom3* rnd = new TRandom3(0);

  for(int i=0; i<1000; i++)
  {
	TString name1=name;
	TString name_beta0=name+"_beta0_";
	TString name_beta1=name+"_beta1_";
	name1+=i;
	name_beta0+=i;
	name_beta1+=i;

	temp = (TH1D*)file->Get(name1);
	average+=temp->Integral();
	//hAvg->Fill(temp->Integral());

	temp_beta0 = (TH1D*)file->Get(name_beta0);
	average_beta0+=temp_beta0->Integral();
	//hAvg_beta0->Fill(temp_beta0->Integral());

	temp_beta1 = (TH1D*)file->Get(name_beta1);
	average_beta1+=temp_beta1->Integral();
	//hAvg_beta1->Fill(temp_beta1->Integral());

	for(int j=1; j<=temp->GetNbinsX(); j++)
	{
	  spectra->Fill(temp->GetXaxis()->GetBinCenter(j),temp->GetBinContent(j));
	  spectra_beta0->Fill(temp_beta0->GetXaxis()->GetBinCenter(j),temp_beta0->GetBinContent(j));
	  spectra_beta1->Fill(temp_beta1->GetXaxis()->GetBinCenter(j),temp_beta1->GetBinContent(j));
	}
  }

  average/=2000.0;
  average_beta0/=2000.0;
  average_beta1/=2000.0;
  //   printf("combined: %f %f %f\n",average,hAvg->GetMean(),hAvg->GetRMS());
  //   printf("beta=0: %f %f %f\n",average_beta0,hAvg_beta0->GetMean(),hAvg_beta0->GetRMS());
  //   printf("beta=1: %f %f %f\n",average_beta1,hAvg_beta1->GetMean(),hAvg_beta1->GetRMS());
  printf("combined average: %f\n",average);
  printf("beta=0 average: %f\n",average_beta0);
  printf("beta=1 average: %f\n",average_beta1);


  TCanvas* c = new TCanvas("c","c",0,0,700,700);
  spectra->Draw("colz");
  c->Update();
  //   hAvg->Draw();
  //c->Update();
  spectra_beta0->Draw("colz");
  c->Update();
  //   hAvg_beta0->Draw();
  //c->Update();
  spectra_beta1->Draw("colz");
  c->Update();
  //   hAvg_beta1->Draw();
  //c->Update();

  TFile* outfile = new TFile("spectra_beta.root","UPDATE");

  double x[100], y[100], ex[100], ey[100], x_beta0[100], y_beta0[100], ex_beta0[100], ey_beta0[100], x_beta1[100], y_beta1[100], ex_beta1[100], ey_beta1[100];
  c->Clear();

  TH1D* proj=0;
  for(int i=1; i<=nbinsx; i++)
  {
	// Combined
	std::cout << "-------------------- combined (" << i << ")-------------------" << std::endl;
	TString namepy_beta=name+"_py";
	namepy_beta+=i;
	proj=spectra->ProjectionY(namepy_beta,i,i);

	x[i-1]=spectra->GetXaxis()->GetBinCenter(i);
	ex[i-1]=spectra->GetXaxis()->GetBinWidth(i)/2.0;
	y[i-1]=proj->GetMean();
	ey[i-1]=proj->GetRMS();
	TFitResultPtr r = proj->Fit("gaus","QS","goff",y[i-1]-ey[i-1],y[i-1]+ey[i-1]);
	if(y[i-1]>0.02 && r->Parameter(1)>0)
	  y[i-1]=r->Parameter(1);
	proj->Write(namepy_beta);

	// beta=0
	std::cout << "-------------------- beta = 0 (" << i << ")-------------------" << std::endl;
	TString namepy_beta0=name+"_py_beta0";
	namepy_beta0+=i;
	proj=spectra_beta0->ProjectionY(namepy_beta0,i,i);

	x_beta0[i-1]=spectra_beta0->GetXaxis()->GetBinCenter(i);
	ex_beta0[i-1]=spectra_beta0->GetXaxis()->GetBinWidth(i)/2.0;
	y_beta0[i-1]=proj->GetMean();
	ey_beta0[i-1]=proj->GetRMS();
	TFitResultPtr r_beta0 = proj->Fit("gaus","QS","goff",y_beta0[i-1]-ey_beta0[i-1],y_beta0[i-1]+ey_beta0[i-1]);
	if(y_beta0[i-1]>0.02 && r_beta0->Parameter(1)>0)
	  y_beta0[i-1]=r_beta0->Parameter(1);
	proj->Write(namepy_beta0);

	// beta=1
	std::cout << "-------------------- beta = 1 (" << i << ")-------------------" << std::endl;
	TString namepy_beta1=name+"_py_beta1";
	namepy_beta1+=i;
	proj=spectra_beta1->ProjectionY(namepy_beta1,i,i);


	x_beta1[i-1]=spectra_beta1->GetXaxis()->GetBinCenter(i);
	ex_beta1[i-1]=spectra_beta1->GetXaxis()->GetBinWidth(i)/2.0;
	y_beta1[i-1]=proj->GetMean();
	ey_beta1[i-1]=proj->GetRMS();
	TFitResultPtr r_beta1 = proj->Fit("gaus","QS","goff",y_beta1[i-1]-ey_beta1[i-1],y_beta1[i-1]+ey_beta1[i-1]);
	if(y_beta1[i-1]>0.02 && r_beta1->Parameter(1)>0)
	  y_beta1[i-1]=r_beta1->Parameter(1);
	proj->Write(namepy_beta1);


  }

  TGraphErrors* errorbars = new TGraphErrors(nbinsx,x,y,ex,ey);
  errorbars->Draw("A E2");

  TGraphErrors* errorbars_beta0 = new TGraphErrors(nbinsx,x_beta0,y_beta0,ex_beta0,ey_beta0);
  errorbars_beta0->Draw("A E2");

  TGraphErrors* errorbars_beta1 = new TGraphErrors(nbinsx,x_beta1,y_beta1,ex_beta1,ey_beta1);
  errorbars_beta1->Draw("A E2");

  // make TH1D with y-errors only
  TString histname = name+"_hist1d";
  TH1D *hist1d = (TH1D*)temp->Clone(histname);
  hist1d->Reset();
  for (int i=0; i<hist1d->GetXaxis()->GetNbins(); i++)
  {
	hist1d->SetBinContent(i+1,y[i]);
	hist1d->SetBinError(i+1,ey[i]);
  }

  histname = name+"_hist1d_beta0";
  TH1D *hist1d_beta0 = (TH1D*)temp_beta0->Clone(histname);
  hist1d_beta0->Reset();
  for (int i=0; i<hist1d_beta0->GetXaxis()->GetNbins(); i++)
  {
	hist1d_beta0->SetBinContent(i+1,y_beta0[i]);
	hist1d_beta0->SetBinError(i+1,ey_beta0[i]);
  }

  histname = name+"_hist1d_beta1";
  TH1D *hist1d_beta1 = (TH1D*)temp_beta1->Clone(histname);
  hist1d_beta1->Reset();
  for (int i=0; i<hist1d_beta1->GetXaxis()->GetNbins(); i++)
  {
	hist1d_beta1->SetBinContent(i+1,y_beta1[i]);
	hist1d_beta1->SetBinError(i+1,ey_beta1[i]);
  }

  TString dummyname;
  if (name=="sk_nue"){dummyname = "nue";}
  if (name=="sk_nuebar"){dummyname = "nuebar";}
  if (name=="sk_nue1pi"){dummyname = "nue1pi";}
  if (name=="sk_numu"){dummyname = "numu";}
  if (name=="sk_numubar"){dummyname = "numubar";}

  outfile->cd();
  errorbars->Write(name);
  histname = "spectra2d_"+dummyname;
  spectra->Write(histname);
  histname = "hist_1d_"+dummyname;
  hist1d->Write(histname);
  histname = name+"_beta0";
  errorbars_beta0->Write(histname);
  histname = "spectra2d_"+dummyname+"_beta0";
  spectra_beta0->Write(histname);
  histname = "hist_1d_"+dummyname+"_beta0";
  hist1d_beta0->Write(histname);
  histname = name+"_beta1";
  errorbars_beta1->Write(histname);
  histname = "spectra2d_"+dummyname+"_beta1";
  spectra_beta1->Write(histname);
  histname = "hist_1d_"+dummyname+"_beta1";
  hist1d_beta1->Write(histname);

  // Finally, calculate likelihood ratio between beta=1 and beta=0 (NOTE: THIS WILL ONLY WORK IF MC IS BINNED THE SAME AS THE DATA)
  // Commented out for now
  /*double negLogL_beta0=0, negLogL_beta1=0;
	TFile *f_dat = new TFile("inputs/Run5-6Data_June2015.root","open");
	TH1D *dathist = (TH1D*)f_dat->Get("nue");

	for (int i=0; i<hist1d_beta0->GetXaxis()->GetNbins(); i++)
	{
	double mc = hist1d_beta0->GetBinContent(i+1);
	double dat = dathist->GetBinContent(i+1);
	std::cout << mc << "  " << dat << std::endl;
	if (dat == 0)
	negLogL_beta0 += (mc - dat);
	else
	negLogL_beta0 += (mc - dat + dat*TMath::Log(dat/mc));

	std::cout << negLogL_beta0 << std::endl;
	}

	for (int i=0; i<hist1d_beta1->GetXaxis()->GetNbins(); i++)
	{
	double mc = hist1d_beta1->GetBinContent(i+1);
	double dat = dathist->GetBinContent(i+1);
	std::cout << mc << "  " << dat << std::endl;
	if (dat ==0)
	negLogL_beta1 += (mc - dat);
	else
	negLogL_beta1 += (mc - dat + dat*TMath::Log(dat/mc));

	std::cout << negLogL_beta1 << std::endl;
	}

	std::cout << "negLogL_beta0 = " << negLogL_beta0 << ", negLogL_beta1 = " << negLogL_beta1 << std::endl;
	std::cout << "Likelihood ratio = " << TMath::Exp(negLogL_beta1-negLogL_beta0) << std::endl;
	*/
}

void makeSKspectra_beta(char *infile="sk_nominal.root")
{
  std::cout << std::endl << " ---------------- numu (beta) ----------------" <<std::endl << std::endl;
  makeSKspectra_beta_sample(infile, "sk_numu", 400, 80);
  std::cout << std::endl << " ---------------- nue (beta) ----------------" <<std::endl << std::endl;
  makeSKspectra_beta_sample(infile, "sk_nue", 400, 14);
  std::cout << std::endl << " ---------------- nue1pi (beta) ----------------" <<std::endl << std::endl;
  makeSKspectra_beta_sample(infile, "sk_nue1pi", 400, 14);
  std::cout << std::endl << " ---------------- numubar (beta) ----------------" <<std::endl << std::endl;
  makeSKspectra_beta_sample(infile, "sk_numubar", 400, 10);
  std::cout << std::endl << " ---------------- nuebar (beta) ----------------" <<std::endl << std::endl;
  makeSKspectra_beta_sample(infile, "sk_nuebar", 100, 6);
}

// -------------------------------------------------------------------------------- //



