#include <iostream>
#include <vector>

#include <TH1D.h>
#include <TLine.h>
#include <TStyle.h>

#include "samplePDF/samplePDFND2018.h"
#include "samplePDF/samplePDFND2018GPU.h"
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

  std::vector<double> bestfit = fitMan->getOscParameters();

  TString production = fitMan->getNDRuns();

#ifdef CUDA
  samplePDFND2018GPU* sample = new samplePDFND2018GPU(fitMan);
#else
  samplePDFND2018* sample = new samplePDFND2018(fitMan);
#endif

  // Get inputted systematic parameters covariance matrices
  TString xsecCovMatrixFile = fitMan->getXsecCovMatrix();
  TString fluxCovMatrixFile = fitMan->getFluxCovMatrix();
  TString ndDetCovMatrixFile = fitMan->getNDdetCovMatrix();

  covarianceFlux* flux = new covarianceFlux("total_flux_cov", fluxCovMatrixFile);
  covarianceNDDet_2014Simple* det = new covarianceNDDet_2014Simple("nddet_cov", ndDetCovMatrixFile);
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
  TDirectory *Flux_LLH = output->mkdir("Flux_LLH");
  TDirectory *XSec_LLH = output->mkdir("XSec_LLH");
  TDirectory *ND280_LLH = output->mkdir("ND280_LLH");
  TDirectory *Sample_LLH = output->mkdir("Sample_LLH");
  TDirectory *Total_LLH = output->mkdir("Total_LLH");

  // Number of points we do for each LLH scan
  int n_points = 150;
  // We print 5 reweights
  int countwidth = double(n_points)/double(5);

  // Make a vector of covarianceBases to loop over and calculate LLH variations for
  std::vector<covarianceBase*> BaseVector;
  // You should be able to add in any arbitrary covariance base class here
  //  BaseVector.push_back(flux);
  // BaseVector.push_back(det);
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

    // Scan over all the parameters
    for (int i = 0; i < (*it)->getSize(); i++) {
      // Get the parameter name
      std::string name = (*it)->getParName(i);
      // For xsec we can get the actual name, hurray for being informative
      if (isxsec) {
        name = TempClass->getParameterName(i);
      }

      // Get the parameter priors and bounds
      //double nom = (*it)->getParInit(i);
      double nom = bestfit[i];
      // Get the covariance matrix and do the +/- 2 sigma
      int nSigma = 1;
      double lower = nom - nSigma*sqrt((*cov)(i,i));//(*it)->getNominal(i) - nSigma*sqrt((*cov)(i,i));
      double upper = nom + nSigma*sqrt((*cov)(i,i));//(*it)->getNominal(i) + nSigma*sqrt((*cov)(i,i));

      if (isxsec) {
        if (lower < TempClass->getParamLowerBound(i)) {
          lower = TempClass->getParamLowerBound(i);
        }
        if (upper > TempClass->getParamUpperBound(i)) {
          upper = TempClass->getParamUpperBound(i);
        }
      }

      std::cout << "Scanning " << name << " with " << n_points << " steps, \nfrom " << lower << " - " << upper << ", nominal = " << nom << std::endl;

      // Make the TH1D
      TH1D *hScan = new TH1D((name+"_full").c_str(), (name+"_full").c_str(), n_points, lower, upper);
      hScan->SetTitle(std::string(std::string("2LLH_full, ") + name + ";" + name + "; -2(ln L_{sample} + ln L_{xsec} + ln L_{flux} + ln L_{det})").c_str());

      TH1D *hScanSam = new TH1D((name+"_sam").c_str(), (name+"_sam").c_str(), n_points, lower, upper);
      hScanSam->SetTitle(std::string(std::string("2LLH_sam, ") + name + ";" + name + "; -2(ln L_{sample})").c_str());

      TH1D *hScanXSec = new TH1D((name+"_xs").c_str(), (name+"_xs").c_str(), n_points, lower, upper);
      hScanXSec->SetTitle(std::string(std::string("2LLH_xs, ") + name + ";" + name + "; -2(ln L_{xsec})").c_str());

      TH1D *hScanFlux = new TH1D((name+"_flux").c_str(), (name+"_flux").c_str(), n_points, lower, upper);
      hScanFlux->SetTitle(std::string(std::string("2LLH_flux, ") + name + ";" + name + "; -2(ln L_{flux})").c_str());

      TH1D *hScanND = new TH1D((name+"_nd").c_str(), (name+"_nd").c_str(), n_points, lower, upper);
      hScanND->SetTitle(std::string(std::string("2LLH_nd, ") + name + ";" + name + "; -2(ln L_{ND280})").c_str());

      // Scan over the parameter space
      for (int j = 0; j < n_points; j++) {

        if (j % countwidth == 0) {
          std::cout << j << "/" << n_points << " (" << double(j)/double(n_points) * 100 << "%)" << std::endl;
        }

        // Set the parameter
        (*it)->setParProp(i, hScan->GetBinCenter(j+1));

        // Reweight the MC
        sample->reweight(fake);

        // Get the -log L likelihoods
        double samplellh = sample->getLikelihood();
        double xsecllh = sample->GetXsecCov()->getLikelihood();
        double fluxllh = sample->GetFluxCov()->getLikelihood();
        double nd280llh = sample->GetSimpleDetCov()->getLikelihood();
        double totalllh = samplellh + xsecllh + fluxllh + nd280llh;

        hScanFlux->SetBinContent(j+1, 2*fluxllh);
        hScanXSec->SetBinContent(j+1, 2*xsecllh);
        hScanND->SetBinContent(j+1, 2*nd280llh);
        hScanSam->SetBinContent(j+1, 2*samplellh);
        hScan->SetBinContent(j+1, 2*totalllh);
      }

      Flux_LLH->cd();
      hScanFlux->Write();
      XSec_LLH->cd();
      hScanXSec->Write();
      ND280_LLH->cd();
      hScanND->Write();
      Sample_LLH->cd();
      hScanSam->Write();
      Total_LLH->cd();
      hScan->Write();

      delete hScanFlux;
      delete hScanXSec;
      delete hScanND;
      delete hScanSam;
      delete hScan;

      hScanFlux = NULL;
      hScanXSec = NULL;
      hScanND = NULL;
      hScanSam = NULL;
      hScan = NULL;

      // Reset the current parameter to its prior for the next parameter scan
      (*it)->setParProp(i, nom);
    }
  }

  output->Close();
  std::cout << "Wrote scan to " << output->GetName() << std::endl;

  delete output;


  return 0;

}
