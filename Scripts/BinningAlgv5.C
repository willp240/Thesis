#include <TH2Poly.h>
#include <TH2D.h>
#include <vector>
#include <TMath.h>
#include <TFile.h>
#include <TTree.h>
#include <THStack.h>
#include <iostream>

  ///////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////
  // Script to set the th2poly binning. First we read in ttree of p-theta for events and fill them in a finely binned histogram. Then we scroll through fixed theta rows, in small steps in MeV. After each step, if there's enough events in the bin (and we the width is above detector resolution) then we keep the bin. If not we keep scrolling until we do
  ///////////////////////////////////////////////////////////////////////////////////



int BinningAlgv5(){

  bool debug=false;
  bool AllPlots=false;
  bool drawErec=false;

  //////////////////////////////////////////////////////////////////////////////////
  //////// Make canvas, open output files
  TCanvas *c = new TCanvas("canv", "canv", 1080, 1080);
  c->SetTopMargin(0.07);
  c->SetBottomMargin(0.11);
  c->SetLeftMargin(0.10);
  c->SetRightMargin(0.18);
  std::string name = "minwidth100";
  std::string pdfFile = (name + ".pdf").c_str();
  c->Print((pdfFile+std::string("[")).c_str());
  gStyle->SetPalette(51);
  TFile * polyFile = new TFile("minwidth100.root","RECREATE");

  int bincount=0;

  //////////////////////////////////////////////////////////////////////////////////
  //////// Open MC and Data files and TTrees
  TFile* fileMC = TFile::Open("../eventTrees/MCeventTree2to8.root");
  TTree* MCtree = (TTree*)fileMC->Get("MCevtTree");
  TFile* fileDATA = TFile::Open("../eventTrees/datEventTree2to8.root");
  TTree* DATAtree = (TTree*)fileDATA->Get("DATevtTree");
  double momMC, thetaMC, weightMC, momDat, thetaDat;
  int sampleMC, sampleDat;
  MCtree->SetBranchAddress("momBranch",&momMC);
  MCtree->SetBranchAddress("thetaBranch",&thetaMC);
  MCtree->SetBranchAddress("weightBranch",&weightMC);
  MCtree->SetBranchAddress("sampleBranch",&sampleMC);
  DATAtree->SetBranchAddress("datmomBranch",&momDat);
  DATAtree->SetBranchAddress("datthetaBranch",&thetaDat);
  DATAtree->SetBranchAddress("datsampleBranch",&sampleDat);
  Long64_t nentriesMC = MCtree->GetEntries();
  Long64_t nentriesDat = DATAtree->GetEntries();
  //X axis of fine binning is 1 bin = 1MeV
  double xBinsFine[30001];
  //Y axis is 1 bin = 0.05 cos theta units
  double yBinsFine[201];
  for (int i =0; i<=30000; i++)
    xBinsFine[i] = i;
  for (int i =0; i<=200; i++){
    yBinsFine[i] = -1 + i*0.01;
  }
  //////////////////////////////////////////////////////////////////////////////////
  ////////// Settings for scanning: resolution, step size, hard code cuts
  const double nsamples= 18;
  int globalbin=0;
  double pHigh=0, pLow=0, tHigh=0, tLow=0, binMCIntegral=0, binDatIntegral=0, prevMCInt=0, prevDatInt=0, MCcut=10, DatCut=1, step=25;
  
  
  ///////////////////////////////////////////////////////////////////////////////////
  //////////// Samples
  //psyche sample enum and names
  double sampleEnum[nsamples] = {3,4,5,19,20,21,43,44,45,49,50,51,55,56,57,61,62,63};
  TString sampleName[nsamples] = {"FGD1_numuCC_0pi", "FGD1_numuCC_1pi", "FGD1_numuCC_other", "FGD2_numuCC_0pi", "FGD2_numuCC_1pi", "FGD2_numuCC_other", "FGD1_anti-numuCC_0pi", "FGD1_anti-numuCC_1pi", "FGD1_anti-numuCC_other", "FGD2_anti-numuCC_0pi", "FGD2_anti-numuCC_1pi", "FGD2_anti-numuCC_other", "FGD1_NuMuBkg_CC0pi_in_AntiNu_Mode","FGD1_NuMuBkg_CC1pi_in_AntiNu_Mode", "FGD1_NuMuBkg_CCother_in_AntiNu_Mode", "FGD2_NuMuBkg_CC0pi_in_AntiNu_Mode", "FGD2_NuMuBkg_CC1pi_in_AntiNu_Mode", "FGD2_NuMuBkg_CCother_in_AntiNu_Mode"};
  //number of fixed theta rows for each sample
  double numRows[nsamples] = {20,17,15,21,15,13,11,6,7,11,6,7,8,7,7,8,7,7};
  double lowXInt[nsamples] = {200, 0, 0, 100, 0, 0, 200, 0, 0, 0, 200, 0, 0, 0, 0, 0, 0, 0};
  double highXInt[nsamples] = {600, 600, 3000, 400, 5000, 30000, 700, 500, 500, 700, 600, 1000, 1000, 1000, 1000, 1000, 1000, 1000};
  
  //hard coded values for rows for each sample
  double yBins[nsamples][25]=  {
    {-1.0,0.27,0.38,0.45,0.5,0.55,0.59,0.63,0.67,0.71,0.74,0.77,0.8,0.83,0.86,0.89,0.92,0.96,1.0}, //FGD1 CC0Pi
    {-1,0.1,0.38,0.48,0.55,0.61,0.66,0.71,0.75,0.79,0.83,0.87,0.91,0.95,0.98,1.0},//FGD1 CC1Pi
    {-1,0.28,0.46,0.52,0.58,0.65,0.7,0.76,0.82,0.86,0.9,0.94,0.97,1.0},//FGD1 CCOther
    {-1,-0.98,-0.79,0.12,0.3,0.39,0.46,0.51,0.56,0.6,0.64,0.68,0.72,0.76,0.8,0.84,0.88,0.92,0.96,1.0}, //FGD2 CC0Pi
    {-1,0.25,0.37,0.47,0.57,0.64,0.71,0.78,0.84,0.89,0.93,0.96,0.98,1.0}, //FGD2 CC1Pi
    {-1,0.4,0.55,0.7,0.78,0.83,0.88,0.92,0.95,0.97,0.99,1.0}, //FGD2 CCOther
    {-1,0.53,0.65,0.73,0.79,0.84,0.88,0.92,0.96,1.0}, //FGD1 nubar CC0Pi
    {-1,0.7,0.84,0.92,1.0}, //FGD1 nubar CC1Pi
    {-1,0.41,0.63,0.78,0.92,1.0}, //FGD1 nubar CCOther
    {-1,0.5,0.64,0.73,0.79,0.84,0.88,0.92,0.96,1.0}, //FGD2 nubar CC0Pi
    {-1,0.7,0.84,0.93,1.0}, //FGD2 nubar CC1Pi
    {-1,0.61,0.79,0.9,0.95,1.0}, //FGD2 nubar CCOther
    {-1,0.5,0.7,0.83,0.9,0.95,1.0}, //FGD1 bkg CC0Pi
    {-1,0.51,0.71,0.83,0.92,1.0}, //FGD1 bkg CC1Pi
    {-1,0.64,0.79,0.88,0.94,1.0}, //FGD1 bkg CCOther
    {-1,0.5,0.7,0.83,0.9,0.95,1.0}, //FGD2 bkg CC0Pi
    {-1,0.57,0.74,0.85,0.93,1.0}, //FGD2 bkg CC1Pi
    {-1,0.63,0.79,0.88,0.95,1.0}  //FGD2 bkg CCOther
  };

  //Loop over all samples
  for(int s = 0; s<18; s++)
    {
      //Gonna make a histogram with fine bins and fill events from tree
      std::cout << sampleName[s] << std::endl;
      TString histname = "veryFineBin"+sampleName[s];

      //Initial integrals=0
      globalbin = 0;
      binMCIntegral = 0;
      binDatIntegral = 0;
      prevMCInt = 0;
      prevDatInt = 0;

      //Fill finely binned mom and theta histogram to get fixed thetas
      TH2D *veryFineBinHisto = new TH2D(histname, histname, 30000, &xBinsFine[0], 200, &yBinsFine[0]);
      for (Long64_t j=0;j<nentriesMC;j++) {
        MCtree->GetEntry(j);
        if(sampleMC==sampleEnum[s]){
          veryFineBinHisto->Fill(momMC,thetaMC,weightMC);
        }
      }
      //Write the finely binned histo
      polyFile->cd();
      veryFineBinHisto->Write();

      //Integral check
      double sampleInt = veryFineBinHisto->Integral(lowXInt[s],highXInt[s],1,200);
      double sampleIntb = veryFineBinHisto->Integral(lowXInt[s],highXInt[s]);
      double intcheck =0;
      std::cout << "sampleInt = " << sampleInt << " /" << numRows[s] << "= " << sampleInt/numRows[s] << " " << sampleIntb << std::endl;

      //dont have bins below 0-200 MeV
      int currentHigh=200, bins=0;
      bins=numRows[s]-1;

      for(int d=0; d<bins; d++)
	std::cout << yBins[s][d] << std::endl;

      histname = "FineBin"+sampleName[s];
      //Fill finely binned MC histogram
      TH2D *FineBinHisto = new TH2D(histname, histname, 30000, &xBinsFine[0], bins-1, &yBins[s][0]);
      for (Long64_t j=0;j<nentriesMC;j++) {
        MCtree->GetEntry(j);
        if(sampleMC==sampleEnum[s]){
          FineBinHisto->Fill(momMC,thetaMC,weightMC);
        }
      }
      //Fill finely binned data histogram
      histname = "FineBinData"+sampleName[s];
      TH2D *FineBinHistoData = new TH2D(histname, histname, 30000, &xBinsFine[0], bins-1, &yBins[s][0]);
      for (Long64_t j=0;j<nentriesDat;j++) {
        DATAtree->GetEntry(j);
        if(sampleDat==sampleEnum[s]){
          FineBinHistoData->Fill(momDat,thetaDat);
        }
      }

      //Make TH2Poly
      std::cout << "Declaring TH2Polys..." << std::endl;
      TH2Poly *hMC = new TH2Poly();
      TH2Poly *hDAT = new TH2Poly();
      histname = sampleName[s];
      std::cout << "Setting titles..." << std::endl;
      hMC->SetName(histname);
      hMC->SetTitle(histname);
      hMC->GetXaxis()->SetTitle("Momentum");
      hMC->GetYaxis()->SetTitle("Theta");
      hMC->GetZaxis()->SetTitle("Events");
      histname = "DATA" + sampleName[s];
      hDAT->SetName(histname);
      hDAT->SetTitle(histname);
      hDAT->GetXaxis()->SetTitle("Momentum");
      hDAT->GetYaxis()->SetTitle("Theta");
      hDAT->GetZaxis()->SetTitle("Events");

      //Vectors saving bin vertices
      vector<vector<double> > xbin;
      vector<vector<double> > ybin;

      //Loop over each row in Y
      std::cout << "Loop over theta rows..." << std::endl;
      for (int tCounter = bins-1; tCounter > 0; tCounter--)
        {
	  std::cout << "Row " << tCounter << ": " << yBins[s][tCounter] << " - " << yBins[s][tCounter-1] << std::endl;
      
	  //Get the current min and max p/theta for this bin
          pHigh = 30000;
          pLow = pHigh-GetMomRes(pHigh);
          tLow = yBins[s][tCounter-1];
          tHigh= yBins[s][tCounter];
          prevMCInt=0;
	  prevDatInt=0;
	  //Loop across this row
	  while(pLow>=0)
            {
		
              binMCIntegral=0;
              binMCIntegral = FineBinHisto->Integral(pLow+1,pHigh,tCounter,tCounter);
	      binDatIntegral=0;
	      binDatIntegral = FineBinHistoData->Integral(pLow+1,pHigh,tCounter,tCounter);
              if(debug)
		std::cout << "pLow = " << FineBinHisto->GetXaxis()->GetBinLowEdge(pLow+1) << ", pHigh = " << FineBinHisto->GetXaxis()->GetBinUpEdge(pHigh) << ", tLow = " << FineBinHisto->GetYaxis()->GetBinLowEdge(tCounter)  << ", tHigh = " << FineBinHisto->GetYaxis()->GetBinUpEdge(tCounter)  << " has " << binMCIntegral << " MC events, " << binDatIntegral << " data events, tCounter " << tCounter << std::endl;

              if(binMCIntegral > MCcut && binDatIntegral > DatCut){//If the bin has > cut events, keep it!
		//Push back bin vertices to vector saving bins
                vector<double> binX;
                binX.push_back(pLow);
                binX.push_back(pHigh);
                binX.push_back(pHigh);
                binX.push_back(pLow);

                vector<double> binY;
                binY.push_back(tLow);
                binY.push_back(tLow);
                binY.push_back(tHigh);
                binY.push_back(tHigh);

                xbin.push_back(binX);
                ybin.push_back(binY);

                pHigh=pLow;//Move along 'one step' in row
		 pLow=pLow-GetMomRes(pHigh);
		//pLow=pLow-step;
                prevMCInt = 0;
		prevDatInt = 0;
                globalbin++;
                if(debug)
		  std::cout << "so keeping. Now pLow = " << pLow << std::endl;
              }//If not > cut then scroll to the next bin
	      else{
                pLow = pLow - GetMomRes(pHigh);
                if(pLow<0)
                  {
		    int checked = 0;
		    double prevpHigh;
		    prevpHigh = xbin.at(globalbin-1).at(1);
                    pLow=0;
		    if(FineBinHisto->Integral(pLow,prevpHigh,tCounter,tCounter) >= 2*MCcut)
		      {
			for(int p = GetMomRes(pLow); p<prevpHigh-GetMomRes(pLow); p++)
                          {
                            if(FineBinHisto->Integral(0,p,tCounter,tCounter)>MCcut&&FineBinHisto->Integral(p,prevpHigh,tCounter,tCounter)>MCcut)
                              {
				if(debug){
				std::cout << " split last bin now (" << prevpHigh << "-" << p <<") with " << FineBinHisto->Integral(p,prevpHigh,tCounter,tCounter)  << " and (" << p << "-0) with " << FineBinHisto->Integral(0,p,tCounter,tCounter) << std::endl;
				std::cout << "Prev integrals: (" << prevpHigh << "-" << pHigh << ") with " << FineBinHisto->Integral(pHigh,prevpHigh,tCounter,tCounter) << " and (" << pHigh << "-0) " <<FineBinHisto->Integral(0,pHigh,tCounter,tCounter) << std::endl;}


                                xbin.erase(xbin.end());
                                ybin.erase(ybin.end());

                                pLow=p;
                                pHigh=prevpHigh;

                                vector<double> binX;
                                binX.push_back(p);
                                binX.push_back(pHigh);
                                binX.push_back(pHigh);
                                binX.push_back(p);

                                vector<double> binY;
                                binY.push_back(tLow);
                                binY.push_back(tLow);
                                binY.push_back(tHigh);
                                binY.push_back(tHigh);

                                xbin.push_back(binX);
                                ybin.push_back(binY);

                                vector<double> binX2;
                                binX2.push_back(0);
                                binX2.push_back(p);
                                binX2.push_back(p);
                                binX2.push_back(0);

                                vector<double> binY2;
                                binY2.push_back(tLow);
                                binY2.push_back(tLow);
                                binY2.push_back(tHigh);
                                binY2.push_back(tHigh);

                                xbin.push_back(binX2);
                                ybin.push_back(binY2);

				globalbin++;
				checked=1;
                                p=prevpHigh;
				pLow=-GetMomRes(pHigh);
                              }
			  }
		      }

		    if(checked==0)
		      {
			if(debug)
			  std::cout << " so setting pLow to " << pLow;
			double currentMCInt = binMCIntegral;
			double currentDatInt = binDatIntegral;
			while(currentMCInt < MCcut || currentDatInt < DatCut)
			  {
			    if(debug)
			      std::cout << " but MC integral " << currentMCInt << " and Data integral " << currentDatInt << " so set pHigh to " << prevpHigh<< std::endl;
			    currentMCInt = FineBinHisto->Integral(pLow,prevpHigh,tCounter,tCounter);
			    currentDatInt = FineBinHistoData->Integral(pLow,prevpHigh,tCounter,tCounter);
			    xbin.erase(xbin.end());
			    ybin.erase(ybin.end());
			    pHigh=prevpHigh;
			    globalbin--;
			  }
			vector<double> binX;
			binX.push_back(pLow);
			binX.push_back(pHigh);
			binX.push_back(pHigh);
			binX.push_back(pLow);
			
			vector<double> binY;
			binY.push_back(tLow);
			binY.push_back(tLow);
			binY.push_back(tHigh);
			binY.push_back(tHigh);
			
			xbin.push_back(binX);
			ybin.push_back(binY);
			pLow=pLow-GetMomRes(pHigh);
			globalbin++;
			if(debug)
			  std::cout << " so setting pLow to " << pLow << " and keeping" << std::endl;
		      }
                  }
		prevMCInt = binMCIntegral;
		prevDatInt = binDatIntegral;

	      }
            }//End while mom >0
        }//End loop over theta rows

      std::cout << "Looped over all rows" << std::endl;

      bincount = bincount + globalbin;
      std::cout << "Bins = " << globalbin << " total = " << bincount << std::endl;

      for(int k=0; k<globalbin; k++)
        {
	  //Loop over all bins and add each one
          if(debug)
	    std::cout << "Adding bin (" << xbin.at(k).at(0) << "," << ybin.at(k).at(0) << ") (" << xbin.at(k).at(1) << "," << ybin.at(k).at(1) << ") (" << xbin.at(k).at(2)<< "," << ybin.at(k).at(2) << ") (" << xbin.at(k).at(3) << "," << ybin.at(k).at(3) << ")" << std::endl;

          double binofx[] = {xbin.at(k).at(0), xbin.at(k).at(1), xbin.at(k).at(2), xbin.at(k).at(3)};
          double binofy[] = {ybin.at(k).at(0), ybin.at(k).at(1), ybin.at(k).at(2), ybin.at(k).at(3)};
          hMC->AddBin(4,binofx,binofy);
          hDAT->AddBin(4,binofx,binofy);
        }

      //Fill MC histo
      for (Long64_t j=0;j<nentriesMC;j++) {
        MCtree->GetEntry(j);
        if(sampleMC ==sampleEnum[s]){
          hMC->Fill(momMC,thetaMC,weightMC);
        }
      }
      //Fill Dat histo
      for (Long64_t j=0;j<nentriesDat;j++) {
	DATAtree->GetEntry(j);
	if(sampleDat ==sampleEnum[s]){
          hDAT->Fill(momDat,thetaDat);
        }
      }

      //Write it all !
      polyFile->cd();
      hMC->Write();
      veryFineBinHisto->Write();

      veryFineBinHisto->Draw("colz");
      if(AllPlots)
	c->Print(pdfFile.c_str());
      veryFineBinHisto->GetXaxis()->SetRangeUser(0,5000);
      veryFineBinHisto->Draw("colz");
      if(AllPlots)
	c->Print(pdfFile.c_str());

      hMC->Draw("colz");
      if(AllPlots)
	c->Print(pdfFile.c_str());
      hMC->GetXaxis()->SetRangeUser(0,5000);
      hMC->Draw("colz");
      //Attempt to draw some nice Erec lines
      if(drawErec){
	auto leg = new TLegend(0.16,0.12,0.36,0.22);	
	TF1 LineErec = GetErecLine(400);
	LineErec->Draw("");
	leg->AddEntry(Line400,"E_{Rec} = 400MeV","l");
	auto LineErec2 = GetErecLine(600);
	LineErec2->SetLineColor(kGreen);
	LineErec2->Draw("same");
	leg->AddEntry(Line600,"E_{Rec} = 600MeV","l");
	auto LineErec3 = GetErecLine(800);
	LineErec3->SetLineColor(kYellow);
	LineErec3->Draw("same");
	leg->AddEntry(Line800,"E_{Rec} = 800MeV","l");
	auto LineErec4 = GetErecLine(1000);
	LineErec4->SetLineColor(kWhite);
	LineErec4->Draw("same");
	leg->AddEntry(Line1000,"E_{Rec} = 1000MeV","l");

	leg->SetTextColor(kWhite);
	leg->Draw("same");
      }

      c->Print(pdfFile.c_str());

      std::cout << "Number of bins: " << hMC->GetNumberOfBins() << std::endl;
      std::cout << "Minimum MC bin content: " << hMC->GetMinimum() << std::endl;
      std::cout << "Maximum MC bin content: " << hMC->GetMaximum() << std::endl;
      std::cout << "Minimum Data bin content: " << hDAT->GetMinimum() << std::endl;
      std::cout << "Maximum Data bin content: " << hDAT->GetMaximum() << std::endl;

      hDAT->Write();
      hDAT->Draw("colz");
      if(AllPlots)
	c->Print(pdfFile.c_str());
      hDAT->GetXaxis()->SetRangeUser(0,5000);
      hDAT->Draw("colz");
      if(AllPlots)
	c->Print(pdfFile.c_str());

      delete hMC;
      delete hDAT;
      delete FineBinHisto;

    }//end loop over samples
  c->Print((pdfFile+std::string("]")).c_str());
}//end tpolyMake

double GetMomRes(double mom){

  return 100.0;
  if(mom < 1000)
    return 50.0;
  else
    return 100.0;

}


TF1 GetErecLine(double Erec){

  const double MLep = 0.10565837*1000;
  const double MNeutron = 0.93956536*1000;
  const double MProton = 0.93827203*1000;
  TString name = Form("Line%.0f",Erec);

  //  TF1 ErecLine("LineErec","1+((sqrt(x*x-[0]*[0])-[1])/(x)) + ((2*sqrt(x*x-[0]*[0])-[0]*[0]-[1]*[1]+[2]*[2])/(2*[3]*x))",1,5000);
  TF1 ErecLine(name,"((sqrt(x*x-[0]*[0])-[1])/x) + (2*[1]*sqrt(x*x-[0]*[0])-[0]*[0]-[1]*[1]+[2]*[2])/(2*[3]*x)",0,5000);
  
  ErecLine->SetParameter(0,MLep);
  ErecLine->SetParameter(1,MNeutron);
  ErecLine->SetParameter(2,MProton);
  ErecLine->SetParameter(3,Erec);
  ErecLine->SetLineColor(kRed);
  ErecLine->SetLineWidth(1);
  return ErecLine;
}
