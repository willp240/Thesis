// Exectable for running Markov Chain Monte Carlo with 2019 parameterisation
// Creates ND280 sample, covariance classes for xsec, flux and ND280, and makes the Markov Chain and sets it off
#include <iostream>
#include <vector>
#include <sys/time.h>

#include <TH1D.h>
#include <THStack.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TApplication.h>
#include <TStopwatch.h>

// Include the samplePDFs
#include "samplePDF/samplePDFND2019.h"
#include "samplePDF/samplePDFND2019GPU.h"

// Include the MCMC
#include "mcmc/mcmc.h"

// Include the input manager
#include "manager/manager.h"


int main(int argc, char **argv) {

  if (argc != 2) {
    std::cout << "ND280_MCMC_2019 runs MCMC at ND280" << std::endl;
    std::cout << "Syntax is $: ND280_MCMC_2019 config.cfg" << std::endl;
    std::cout << "Where config.cfg is a valid config file, compatible with the manager class (manager/manager.cpp/h); an example is provided in example_config.cfg" << std::endl;
    exit(-1);
  }

  // Have our lovely fit manager act on argument 1, the config file
  // This prints the inputted parameters after loading
  manager *fitMan = new manager(argv[1]);

  // there's a check inside the manager class that does this; left here for demonstrative purposes
  if (fitMan->getGoodConfig() == false) {
    std::cerr << "Didn't find a good config in input configuration" << std::endl;
    throw;
  }

  /// Create probability distribution functions (PDF)
  /// Look if GPU is set.
  /// SamplePDFND2019 calls init() 
  /// that loads the DataSample for all ND280 runs (MC and data) and add them to the Experiment class
  /// and a create a ToyExperiment to compute systematics later.
#ifdef CUDA
  samplePDFND2019GPU* sample = new samplePDFND2019GPU(fitMan);
#else
  samplePDFND2019* sample = new samplePDFND2019(fitMan);
#endif

  std::cout << "------------------------------" << std::endl;


  // Get inputted systematic parameters covariance matrices
  TString xsecCovMatrixFile = fitMan->getXsecCovMatrix();
  TString fluxCovMatrixFile = fitMan->getFluxCovMatrix();
  TString ndDetCovMatrixFile = fitMan->getNDdetCovMatrix();

  // Setup the covariance matrices
  covarianceXsec2015 *xsec = new covarianceXsec2015("xsec_cov", xsecCovMatrixFile);
  std::cout << "------------------------------" << std::endl;
  covarianceFlux* flux = new covarianceFlux("total_flux_cov", fluxCovMatrixFile);
  std::cout << "------------------------------" << std::endl;
  covarianceNDDet_2019Poly* det = new covarianceNDDet_2019Poly("nddet_cov", ndDetCovMatrixFile);
  std::cout << "------------------------------" << std::endl;

  // Fill the parameter values with their nominal values
  // should _ALWAYS_ be done before overriding with fix or flat
  xsec->setParameters();
  flux->setParameters();
  det->setParameters();

  // ==========================================================
  // Set priors and fix options from file
  // read priors as string from config file
  std::vector<int> NDDetFlatParams = fitMan->getNdDetFlat();
  std::vector<int> NDDetFixParams  = fitMan->getNdDetFix();
  std::vector<int> FluxFlatParams  = fitMan->getFluxFlat();
  std::vector<int> FluxFixParams   = fitMan->getFluxFix();
  std::vector<int> XsecFlatParams  = fitMan->getXsecFlat();
  std::vector<int> XsecFixParams   = fitMan->getXsecFix();

  std::vector<double> Eb = fitMan->getOscParameters();

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
      flux->toggleFixParameter(j);
    }
  } else {
    for (int j = 0; j < FluxFixParams.size(); j++) {
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
      det->toggleFixParameter(j);
    }
  } else {
    for (int j = 0; j < NDDetFixParams.size(); j++) {
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
      xsec->toggleFixParameter(j);
    }
  } else {
    for (int j = 0; j < XsecFixParams.size(); j++) {
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

  xsec->setStepScale(fitMan->getXsecStepScale());
  flux->setStepScale(fitMan->getFluxStepScale());
  det->setStepScale(fitMan->getNdDetStepScale());

  // Now fill the data histograms
  std::cout << "------------------------------" << std::endl;
  TStopwatch clock;
  clock.Start();
  sample->fillDataFromSamples();
  clock.Stop();
  std::cout << "Data fill took " << clock.RealTime() << "s" << std::endl;

  // Now fill the MC histograms
  std::cout << "------------------------------" << std::endl;
  clock.Start();
  sample->fillReweigtingBins();
  clock.Stop();
  std::cout << "MC fill took " << clock.RealTime() << "s" << std::endl;
  std::cout << "------------------------------" << std::endl;

  // If CUDA GPUs are available and enabled in setup, fill the GPU splines
  // Need to call this after fillRewightingBins!
#ifdef CUDA
  clock.Start();
  sample->fillGPUSplines();
  clock.Stop();
  std::cout << "GPU spline fill took " << clock.RealTime() << "s" << std::endl;
#endif

  std::cout << "------------------------------" << std::endl;
  std::cout << "Setup complete! " << std::endl;

  // Do we want oscillation parameters at ND280, I think...
  // Very probably _NOT_ supported!
  double *fake = NULL;

  //WP set eb to whatever we want in config, for making Eb templates
  xsec->setPar(Eb[0],Eb[1]);
  xsec->setPar(Eb[2],Eb[3]);
  xsec->setPar(Eb[4],Eb[5]);
  xsec->setPar(Eb[6],Eb[7]);
  
  // Reweight the Monte-Carlo from nominal to whatever settings we've specified
  clock.Start();
  sample->reweight(fake);
  sample->printRates();
  clock.Stop();
  std::cout << "Single reweight of MC took " << clock.RealTime() << "s to complete" << std::endl;

  // If we're running fake-data fits, we should load up this as data
  bool fakeDataFit = fitMan->getFakeDataFitFlag();
  if (fakeDataFit) {
    std::string FileName = fitMan->getFakeDataFilename();
    if (FileName.empty()) {
      sample->setAsimovFakeData(fitMan->getCustomReWeight());
    } else {
      sample->setAsimovFakeData_FromFile(FileName);
    }
  }

  // Find what we want to do from the config file
  bool SaveNom = fitMan->getSaveNom();
  bool MakeFake = fitMan->getMakeAsimovOnly();
  // If we want to save nominal or just make fake data to file need to save the PDFs
  if (SaveNom || MakeFake) {
    // Save the nominal pdfs after first reweight
    sample->savePDF(std::string(fitMan->getOutputFilename()));
    // If we only want to make fake-data, we can close program now (no need to fit)
    if (MakeFake) {
      std::cout << "Made Asimov fake-data so am now exiting! Bye!" << std::endl;
      return 0;
    }
  }

  // Once we've set the fake-data, move the sample to a random starting position
  bool RandomStart = fitMan->getRandomStart();
  if (RandomStart) {
    sample->RandomStart();
  }

  // Set up Markov Chain (maybe make this whole thing into a class?)
  mcmc *markovChain = new mcmc(fitMan);
  markovChain->setChainLength(fitMan->getNSteps());

  std::cout << "------------------------------" << std::endl;

  // Add samples (only nd280 sample so far)
  markovChain->addSamplePDF(sample);

  bool useXsec = fitMan->getXsecSystOpt();
  bool useFlux = fitMan->getFluxSystOpt();
  bool useNdDet  = fitMan->getNdDetSystOpt();

  // Add the flux if specified
  if (useFlux) {
    markovChain->addSystObj(flux);
  } 

  // Add the xsec if specified
  if (useXsec) {
    markovChain->addSystObj(xsec);
  }

  // Add ND280 if specified
  if (useNdDet) {
    markovChain->addSystObj(det);
  }

  std::cout << "------------------------------" << std::endl;

  // run rabbit run!
  markovChain->runMCMC();

  std::cout << "------------------------------" << std::endl;
  sample->printRates();
  std::cout << "------------------------------" << std::endl;

  return 0;
}
