#include <iostream>
#include <vector>

#include <TH1D.h>
#include <TLine.h>
#include <TStyle.h>

#include "samplePDF/samplePDFND2019.h"
#include "samplePDF/samplePDFND2019GPU.h"
#include "mcmc/mcmc.h"
#include "manager/manager.h"

#include "omp.h"


int main(int argc, char **argv) {

  manager *fitMan = new manager(argv[1]);
  if (fitMan->getGoodConfig() == false) {
    std::cerr << "************* BAD CONFIG ***************" << std::endl;
    throw;
  }

  /// Now we need to see what kind of fit we are doing
  /// Asimov fit    = fit fake data with same parameters than MC
  /// Toy fit       = fit fake data changing the value of nuisance parameters (?)
  /// Fake data fit = fit fake data which are even more different (?)
  /// Real data fit = fit the true ND and SK data (the true fit !)

  bool fakeDataFit = fitMan -> getFakeDataFitFlag();
  bool realDataFit = fitMan -> getRealDataFitFlag();
  bool toyFit      = fitMan -> getToyFitFlag();
  bool asimovFit   = fitMan -> getAsimovFitFlag();

  // The different fits are not compatible between each other so that check that only one type of fit is set
  if (((int)fakeDataFit + (int)realDataFit + (int)toyFit + (int)asimovFit ) != 1) {
    std::cout << "Error: have not chosen one valid type of fit!" << std::endl
      << "Asimov fit = "    << asimovFit    << std::endl
      << "Toy fit = "       << toyFit       << std::endl
      << "Fake data fit = "  << fakeDataFit  << std::endl
      << "Real data fit = " << realDataFit  << std::endl
      << "You have to chose ONLY ONE of them ! " << std::endl  << std::endl;
    throw;
  }

  TString production = fitMan->getNDRuns();

#ifdef CUDA
  samplePDFND2019GPU* sample = new samplePDFND2019GPU(fitMan);
#else
  samplePDFND2019* sample = new samplePDFND2019(fitMan);
#endif

  // Get inputted systematic parameters covariance matrices
  TString xsecCovMatrixFile = fitMan->getXsecCovMatrix();
  TString fluxCovMatrixFile = fitMan->getFluxCovMatrix();
  TString ndDetCovMatrixFile = fitMan->getNDdetCovMatrix();

  covarianceFlux* flux = new covarianceFlux("total_flux_cov", fluxCovMatrixFile);
  covarianceNDDet_2019Poly* det = new covarianceNDDet_2019Poly("nddet_cov", ndDetCovMatrixFile);
  covarianceXsec2015 *xsec = new covarianceXsec2015("xsec_cov", xsecCovMatrixFile);

  // Fill the parameter values
  xsec->setParameters();
  flux->setParameters();
  det->setParameters();

  // Set priors and fix options from file
  // read priors as string from config file
  std::vector<int> NDDetFlatParams = fitMan->getNdDetFlat();
  std::vector<int> NDDetFixParams  = fitMan->getNdDetFix();
  std::vector<int> FluxFlatParams  = fitMan->getFluxFlat();
  std::vector<int> FluxFixParams   = fitMan->getFluxFix();
  std::vector<int> XsecFlatParams  = fitMan->getXsecFlat();
  std::vector<int> XsecFixParams   = fitMan->getXsecFix();

  // Flat flux parameters loop
  if (FluxFlatParams.size() == 1 && FluxFlatParams.at(0) == -1) {
    for (int j = 0; j < flux->getSize(); j++) {
      flux->setEvalLikelihood(j, false);
    }
  } else {
    for (int j = 0; j < FluxFlatParams.size(); j++) {
      flux->setEvalLikelihood(FluxFlatParams.at(j), false);
    }
  }

  // Fixed flux parameters loop
  if (FluxFixParams.size() == 1 && FluxFixParams.at(0) == -1) {
    for (int j = 0; j < flux->getSize(); j++) {
      flux->setParCurrProp(j, flux->getNominal(j));
      flux->toggleFixParameter(j);
    }
  } else {
    for (int j = 0; j < FluxFixParams.size(); j++) {
      flux->setParCurrProp(j, flux->getNominal(j));
      flux->toggleFixParameter(FluxFixParams.at(j));
    }
  }

  std::cout << "------------------------------" << std::endl;

  // Flat ND parameters loop
  if (NDDetFlatParams.size() == 1 && NDDetFlatParams.at(0) == -1) {
    for (int j = 0; j < det->getSize(); j++) {
      det->setEvalLikelihood(j, false);
    }
  } else {
    for (int j = 0; j < NDDetFlatParams.size(); j++) {
      det->setEvalLikelihood(NDDetFlatParams.at(j), false);
    }
  }

  // Fixed ND parameters loop
  if (NDDetFixParams.size() == 1 && NDDetFixParams.at(0) == -1) {
    for (int j = 0; j < det->getSize(); j++) {
      det->setParCurrProp(j, det->getNominal(j));
      det->toggleFixParameter(j);
    }
  } else {
    for (int j = 0; j < NDDetFixParams.size(); j++) {
      det->setParCurrProp(j, det->getNominal(j));
      det->toggleFixParameter(NDDetFixParams.at(j));
    }
  }

  std::cout << "------------------------------" << std::endl;

  // Flat xsec parameters loop
  if (XsecFlatParams.size() == 1 && XsecFlatParams.at(0) == -1) {
    for (int j = 0; j < xsec->getSize(); j++) {
      xsec->setEvalLikelihood(j, false);
    }
  } else {
    for (int j = 0; j < XsecFlatParams.size(); j++) {
      xsec->setEvalLikelihood(XsecFlatParams.at(j), false);
    }
  }

  // Fixed xsec parameters loop
  if (XsecFixParams.size() == 1 && XsecFixParams.at(0) == -1) {
    for (int j = 0; j < xsec->getSize(); j++) {
      xsec->setParCurrProp(j, xsec->getNominal(j));
      xsec->toggleFixParameter(j);
    }
  } else {
    for (int j = 0; j < XsecFixParams.size(); j++) {
      xsec->setParCurrProp(j, xsec->getNominal(j));
      xsec->toggleFixParameter(XsecFixParams.at(j));
    }
  }

  std::cout << "------------------------------" << std::endl;

  // ==========================================================
  //  Set covariance matrices for ND sample
  sample->setXsecCov(xsec);
  sample->setFluxCov(flux);
  sample->setSimpleDetCov(det);

  std::cout << "------------------------------" << std::endl;

  std::cout << "Filling data" << std::endl;
  sample->fillDataFromSamples();
  std::cout << "Filling MC" << std::endl;
  sample->fillReweigtingBins();
#ifdef CUDA
  std::cout << "filling GPU splines" << std::endl;
  sample->fillGPUSplines();
#endif

  std::cout << "ND280 setup complete! " << std::endl;

  // Time reweighting, e.g. GPU vs CPU, splines, event by event...
  TStopwatch reweightClock;
  reweightClock.Start();

  double *fake = 0;
  sample->reweight(fake);

  reweightClock.Stop();
  std::cout << "Single reweight took " << reweightClock.RealTime() << "s to complete" << std::endl;

  std::cout << "Rates after first reweight: " << std::endl;
  sample->printRates();

  // If we're running fake-data fits, we should load up this as data
  if (fakeDataFit) {
    sample->setAsimovFakeData();
  }

  // Now finally get onto the LLH scan stuff
  // Very similar code to MCMC but never start MCMC; just scan over the parameter space

  // only have one sample for ND280 = sample
  std::string outputName = fitMan->getOutputFilename();
  outputName = outputName.substr(0, outputName.find(".root"));
  outputName += "_ND280logL_scan.root";

  TFile *output = new TFile(outputName.c_str(), "RECREATE");
  TDirectory *XSec_LLH = output->mkdir("XSec_LLH");
  TDirectory *Sample_LLH = output->mkdir("Sample_LLH");
  TDirectory *Total_LLH = output->mkdir("Total_LLH");

  // Number of points we do for each LLH scan
  int n_points = 50;
  // We print 5 reweights
  int countwidth = double(n_points)/double(5);

  // Make a vector of covarianceBases to loop over and calculate LLH variations for
  std::vector<covarianceBase*> BaseVector;
  // You should be able to add in any arbitrary covariance base class here
  //  BaseVector.push_back(flux);
  //BaseVector.push_back(det);
  BaseVector.push_back(xsec);
  bool isxsec = false;
  // Loop over the covariance classes
  for (std::vector<covarianceBase*>::iterator it = BaseVector.begin(); it != BaseVector.end(); ++it) {
    // Get the covariance matrix for the current base vector
    TMatrixDSym *cov = (*it)->getCovMatrix();
    //  Get the xsec2015 class when we're dealing with xsec so we can set names properly
    //  Also we need to identify when we have lower and upper bounds so we don't scan in invalid regions
    //  in 2016 this was only applied to xsec params
    covarianceXsec2015 *TempClass = NULL;
    if (std::string((*it)->getName()) == "xsec_cov") {
      isxsec = true;
      TempClass = dynamic_cast<covarianceXsec2015*>(*it);
    } else {
      TempClass = NULL;
    }

    // Scan over all the Q2 normalisation parameters
    for (int i = 10; i < 15; i++) {
      // Get the parameter name
      std::string name = (*it)->getParName(i);
      // For xsec we can get the actual name, hurray for being informative
      if (isxsec) {
        name = TempClass->getParameterName(i);
      }

      // Get the parameter priors and bounds
      double nom = (*it)->getParInit(i);

      // Get the covariance matrix and do the +/- 2 sigma
      int nSigma = 1;
      double lower = (*it)->getNominal(i) - nSigma*sqrt((*cov)(i,i));
      double upper = (*it)->getNominal(i) + nSigma*sqrt((*cov)(i,i));

      if (isxsec) {
        if (lower < TempClass->getParamLowerBound(i)) {
          lower = TempClass->getParamLowerBound(i);
        }
        if (upper > TempClass->getParamUpperBound(i)) {
          upper = TempClass->getParamUpperBound(i);
        }
      }

      std::cout << "Scanning " << name << " with " << n_points << " steps, \nfrom " << lower << " - " << upper << ", nominal = " << nom << std::endl;

      //Now do same for MAQE
      std::string name2 = (*it)->getParName(0);
      // For xsec we can get the actual name, hurray for being informative
      if (isxsec) {
        name2 = TempClass->getParameterName(0);
      }
      // Get the parameter priors and bounds
      double nom2 = (*it)->getParInit(0);

      // Get the covariance matrix and do the +/- 2 sigma
      int nSigma2 = 2;
      double lower2 = (*it)->getNominal(0) - nSigma*sqrt((*cov)(0,0));
      double upper2 = (*it)->getNominal(0) + nSigma*sqrt((*cov)(0,0));

      if (isxsec) {
        if (lower2 < TempClass->getParamLowerBound(0)) {
          lower2 = TempClass->getParamLowerBound(0);
        }
        if (upper2 > TempClass->getParamUpperBound(0)) {
          upper2 = TempClass->getParamUpperBound(0);
        }
      }
      //now do same for maqe
      std::cout << "And " << name2 << " with " << n_points << " steps, \nfrom " << lower2 << " - " << upper2 << ", nominal = " << nom2 << std::endl;
      
      // Make the TH1D
      TH2D *hScan = new TH2D((name+"_"+name2+"_full").c_str(), (name+"_"+name2+"_full").c_str(), n_points, lower, upper, n_points, lower2, upper2);
      hScan->SetTitle(std::string(std::string("2LLH_full, ") + name + "_" + name2 + ";" + name + ";" + name2 + "; -2(ln L_{sample} + ln L_{xsec} + ln L_{flux} + ln L_{det})").c_str());

      TH2D *hScanSam = new TH2D((name+"_"+name2+"_sam").c_str(), (name+"_"+name2+"_sam").c_str(), n_points, lower, upper, n_points, lower2, upper2);
      hScanSam->SetTitle(std::string(std::string("2LLH_sam, ") + name + "_" + name2 + ";" + name + ";" + name2 + "; -2(ln L_{sample})").c_str());

      TH2D *hScanXSec = new TH2D((name+"_"+name2+"_xs").c_str(), (name+"_"+name2+"_xs").c_str(), n_points, lower, upper, n_points, lower2, upper2);
      hScanXSec->SetTitle(std::string(std::string("2LLH_xs, ") + name + "_" + name2 + ";" + name + ";" + name2 + "; -2(ln L_{xsec})").c_str());


      // Scan over the parameter space
      for (int j = 0; j < n_points; j++) {
	for(int k = 0; k < n_points; k++) {

	  if (j % countwidth == 0) {
	    std::cout << j << "/" << n_points << " (" << double(j)/double(n_points) * 100 << "%)" << std::endl;
	  }

	  // Set the parameter
	  (*it)->setParProp(i, hScan->GetXaxis()->GetBinCenter(j+1));
	  (*it)->setParProp(0, hScan->GetYaxis()->GetBinCenter(k+1));
	
	  // Reweight the MC
	  sample->reweight(fake);
	  
	  // Get the -log L likelihoods
	  double samplellh = sample->getLikelihood();
	  double xsecllh = sample->GetXsecCov()->getLikelihood();
	  double fluxllh = sample->GetFluxCov()->getLikelihood();
	  double nd280llh = sample->GetSimpleDetCov()->getLikelihood();
	  double totalllh = samplellh + xsecllh + fluxllh + nd280llh;
	  
	  hScanXSec->SetBinContent(j+1, k+1, 2*xsecllh);
	  hScanSam->SetBinContent(j+1, k+1, 2*samplellh);
	  hScan->SetBinContent(j+1, k+1, 2*totalllh);
	}
      }
      
      XSec_LLH->cd();
      hScanXSec->Write();
      Sample_LLH->cd();
      hScanSam->Write();
      Total_LLH->cd();
      hScan->Write();

      delete hScanXSec;
      delete hScanSam;
      delete hScan;

      hScanXSec = NULL;
      hScanSam = NULL;
      hScan = NULL;

      // Reset the current parameter to its prior for the next parameter scan
      (*it)->setParProp(i, nom);
      (*it)->setParProp(0, nom);
    }
  }

  output->Close();
  std::cout << "Wrote scan to " << output->GetName() << std::endl;
  
  delete output;
  
  return 0;
}
