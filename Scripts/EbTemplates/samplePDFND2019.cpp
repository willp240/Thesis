#include "NIWGReWeight/NIWGFSLepEbMomentumVariation_utils.h"
#include "samplePDFND2019.h"

using namespace niwg::var::FSLepEBVar;                                                                             
using namespace niwg::var;
using namespace niwg;


// ***************************************************************************
// Alternative constructor to make the sample aware of the manager
samplePDFND2019::samplePDFND2019(manager *Manager) : samplePDFBase(Manager->getPOT()) {
// ***************************************************************************

  std::cout << "Creating ND280 2019 instance" << std::endl;
#if USE_SPLINE == USE_TSpline3
  std::cout << "   Using TSpline3" << std::endl;
#elif USE_SPLINE == USE_TSpline3_red
  std::cout << "   Using TSpline3 reduced" << std::endl;
#elif USE_SPLINE == USE_TF1
  std::cout << "   Using TF1" << std::endl;
#elif USE_SPLINE == USE_TF1_red
  std::cout << "   Using TF1 reduced" << std::endl;
#endif

  // Save a pointer to the manager
  FitManager = Manager;
  // What production we're using
  TString production = FitManager->getNDRuns();
  // Currently not set in the manager. Dictates if we run with simple ND280 detector treatment, i.e. treat ND280 systematics as normalisations in pmu cosmu bins for each selection
  bool isSimple = true;
  // Include making distributions by mode?
  modepdf = false;
  // Are we NOT doing rfg (so norfg == false means run with RFG)
  norfg   = false;
  // Are we using the simple detector configuration
  simple  = isSimple;
  // Random start
  HaveIRandomStart = false;
  // Do some checks on the requested production
  if (production == "") {
    std::cerr << "You gave me an empty production!\nExiting!" << std::endl;
    throw;
  }
  if (!production.Contains("P6")) {
    std::cerr << "I only support P6 sorry!" << std::endl;
    throw;
  }

  TString prod_copy(production);
  prod_copy.Remove(prod_copy.First(' '), prod_copy.Length() - prod_copy.First(' '));
  if (prod_copy != "P6" && prod_copy != "P6_S17" && prod_copy != "P6T") {
    std::cerr << "Did not find good good production, exiting" << std::endl;
    std::cerr << "You gave " << prod_copy << ", I need P6 or P6_S17 or P6T" << std::endl;
    throw;
  }
  // Do we really need a new here? Surely setting the string should be fine
  // Warning, something does indeed break, investigate this in the future
  prod = new TString(production);
  // Hard-coded for backward-compatibility...
  nXsecSplines = -1;
  // Hard-code the number of cross-section parameters with xsec splines
  // This then gets checked against the covarianceXsec class' number of spline parameters
  if (prod->Contains("P6T")) {
    nXsecSplines = 19;
  }
  // Pointer to covarianceFlux
  flux      = NULL;
  // Pointer to covarianceNDDet2014_Simple
  detsimple = NULL;
  // Pointer to covarianceXsec2015
  xs2015    = NULL;
  // Psyche experiment class
  exp = NULL;

  // A nice random number to use at our leisure
  rand = new TRandom3(0);
  // Do we care about the MC pdfs by mode
  samplepdfs  = new TObjArray(0);

  // Holds the data pdfs
  datapdfs    = new TObjArray(0);

  // Holds the MC pdfs by mode
  samplemodepdfs = new TObjArray(0);

  // The struct for the ND280 params
  ndobj = new ND280MC();

  // The cross-section model we're using // 2015a = 0, 2015b = 1, 2015c = 2, 2016a = 3, see Structs.h
  xsecmodel = kNIWG_undefined;

  // Set the ND280 test-statistic
  SetTestStatistic(static_cast<TestStatistic>(FitManager->getNDLikelihood()));

  // Initialise the MC with psyche
  LoadPsyche();
  // Then load the samples
  LoadSamples();

  // Enable the mode histograms AFTER addSelection is called
  if (FitManager->getPlotByMode()) {
    enableModeHistograms();
  }

  // If we're using psyche parameters
  if (!simple) {
    det = new covariancePsyche("NDdet", _man.syst().GetCovarianceMatrix());
  }

  if (std::getenv("NIWG") == NULL) {
    std::cerr << "Need NIWG environment variable to run with Eb" << std::endl;
    std::cerr << "EXPORT NIWG to your NIWGReWeight environment" << std::endl;
    throw;
  }

  // Check the preprocessors that have been defined
  // Can currently only run in verbose DEBUG mode in single thread
#ifdef MULTITHREAD
#pragma omp parallel
  {
#if DEBUG > 0
#pragma omp critical
    if (omp_get_num_threads() != 1) {
      std::cerr << "I'm running in DEBUG on multiple threads, this won't work" << std::endl;
      std::cerr << "export OMP_NUM_THREADS=1 and try again please!" << std::endl;
      throw;
    }
#endif
  }
#endif

}

// ***************************************************************************
// Empty destructor
samplePDFND2019::~samplePDFND2019() {
  // ***************************************************************************
  std::cout << "Destroying samplePDFND2019 object" << std::endl;
  delete prod;
  delete rand;
  delete samplepdfs;
  delete datapdfs;
  delete samplemodepdfs;
  delete ndobj;
  delete[] xsecInfo;
  if (!simple) {
    delete exp;
    delete toy;
  }
  if (xs2015) {
    delete xs2015;
  }
  if (flux) {
    delete flux;
  }
  if (detsimple) {
    delete detsimple;
  }
  delete[] xsec_norm_modes;
  delete[] xsec_norm_startbin;
}

// ***************************************************************************
// The initialiser function, only called in the constructor of samplePDFND2019
// Sets up the interface to psyche via class members and loads the data and MC
void samplePDFND2019::LoadPsyche() {
  // ***************************************************************************

  // Quiet psyche
  //QuietPlease();

  // Create the psyche Experiment class which holds our data and MC
  exp = new Experiment("ND280"); 

  // Make and fill the EventSummary even when the selection is not passed.
  if (ND::params().GetParameterI("psycheSteering.Selections.ForceFillEventSummary")) {
    _man.sel().SetForceFillEventSummary(true);
  }

  // Default location for iRODS data and MC
  std::string DataLocation = "MaCh3/P6Data";
  std::string MCLocation = "MaCh3/P6MC";
  if (std::getenv("MACH3_DATA") != NULL) {
    DataLocation = std::string(std::getenv("MACH3_DATA"));
  }
  if (std::getenv("MACH3_MC") != NULL) {
    MCLocation = std::string(std::getenv("MACH3_MC"));
  }

  if (prod->Contains(" 2a")) {
    DataSample* data2a  = new DataSample(file_exists(DataLocation+"/run2aDataSplines.root"));
    DataSample* mc2a    = new DataSample(file_exists(MCLocation+"/run2aMCsplines.root"));
    //    DataSample* sand2a  = new DataSample(file_exists(MCLocation+"/run2aMCsplines.root"));
    SampleGroup run2a("run2a");
    run2a.AddDataSample(data2a);
    run2a.AddMCSample("magnet", mc2a);
    //run2a.AddMCSample("sand", sand2a);
    exp->AddSampleGroup("run2a", run2a); 
  }

  if (prod->Contains(" 2w")) {
    DataSample* data2w  = new DataSample(file_exists(DataLocation+"/run2wDataSplines.root"));
    DataSample* mc2w    = new DataSample(file_exists(MCLocation+"/run2wMCsplines.root"));
    //DataSample* sand2w  = new DataSample(file_exists(MCLocation+"/run2wMCsplines.root"));
    SampleGroup run2w("run2w");
    run2w.AddDataSample(data2w);
    run2w.AddMCSample("magnet", mc2w);
    //run2w.AddMCSample("sand", sand2w);
    exp->AddSampleGroup("run2w", run2w); 
  }

  if (prod->Contains(" 3")) {
    DataSample* data3 = new DataSample(file_exists(DataLocation+"/run3DataSplines.root"));
    DataSample* mc3   = new DataSample(file_exists(MCLocation+"/run3MCsplines.root"));
    //DataSample* sand3 = new DataSample(file_exists(MCLocation+"/run3MCsplines.root"));
    SampleGroup run3("run3");
    run3.AddDataSample(data3);
    run3.AddMCSample("magnet", mc3);
    //run3.AddMCSample("sand", sand3);
    exp->AddSampleGroup("run3", run3);
  }

  if (prod->Contains(" 4a")) {
    DataSample* data4a  = new DataSample(file_exists(DataLocation+"/run4aDataSplines.root"));
    DataSample* mc4a    = new DataSample(file_exists(MCLocation+"/run4aMCsplines.root"));
    //DataSample* sand4a  = new DataSample(file_exists(MCLocation+"/run4aMCsplines.root"));
    SampleGroup run4a("run4a");
    run4a.AddDataSample(data4a);
    run4a.AddMCSample("magnet", mc4a);
    //run4a.AddMCSample("sand", sand4a);
    exp->AddSampleGroup("run4a", run4a);
  }

  if (prod->Contains(" 4w")) {
    DataSample* data4w  = new DataSample(file_exists(DataLocation+"/run4wDataSplines.root"));
    DataSample* mc4w    = new DataSample(file_exists(MCLocation+"/run4wMCsplines.root"));
    //DataSample* sand4w  = new DataSample(file_exists(MCLocation+"/run4wMCsplines.root"));
    SampleGroup run4w("run4w");
    run4w.AddDataSample(data4w);
    run4w.AddMCSample("magnet", mc4w);
    //run4w.AddMCSample("sand", sand4w);
    exp->AddSampleGroup("run4w", run4w); 
  }

  if (prod->Contains(" 5")) {
    DataSample* data5 = new DataSample(file_exists(DataLocation+"/run5DataSplines.root"));
    DataSample* mc5   = new DataSample(file_exists(MCLocation+"/run5MCsplines.root"));
    //DataSample* sand5 = new DataSample(file_exists(MCLocation+"/run5MCsplines.root"));
    SampleGroup run5("run5");
    run5.AddDataSample(data5);
    run5.AddMCSample("magnet", mc5);
    //run5.AddMCSample("sand", sand5);
    exp->AddSampleGroup("run5w", run5);
  }

  if (prod->Contains(" 6")) {
    DataSample* data6  = new DataSample(file_exists(DataLocation+"/run6DataSplines.root"));
    DataSample* mc6    = new DataSample(file_exists(MCLocation+"/run6MCsplines.root"));
    //DataSample* sand6  = new DataSample(file_exists(MCLocation+"/run6MCsplines.root"));
    SampleGroup run6("run6");
    run6.AddDataSample(data6);
    run6.AddMCSample("magnet", mc6);
    //run6.AddMCSample("sand", sand6);
    exp->AddSampleGroup("run6", run6); 
  }

  if (prod->Contains(" 7")) {
    DataSample* data7 = new DataSample(file_exists(DataLocation+"/run7DataSplines.root"));
    DataSample* mc7   = new DataSample(file_exists(MCLocation+"/run7MCsplines.root"));
    //DataSample* sand7 = new DataSample(file_exists(MCLocation+"/run7MCsplines.root"));
    SampleGroup run7("run7");
    run7.AddDataSample(data7);
    run7.AddMCSample("magnet", mc7);
    //run7.AddMCSample("sand", sand7);
    exp->AddSampleGroup("run7", run7); 
  }

  if (prod->Contains(" 8a")) {
    DataSample* data8a = new DataSample(file_exists(DataLocation+"/run8aDataSplines.root"));
    DataSample* mc8a   = new DataSample(file_exists(MCLocation+"/run8aMCsplines.root"));
    //DataSample* sand8a = new DataSample(file_exists(MCLocation+"/run8aMCsplines.root"));
    SampleGroup run8a("run8a");
    run8a.AddDataSample(data8a);
    run8a.AddMCSample("magnet", mc8a);
    //run8a.AddMCSample("sand", sand8a);
    exp->AddSampleGroup("run8a", run8a); 
  }

  if (prod->Contains(" 8w")) {
    DataSample* data8w = new DataSample(file_exists(DataLocation+"/run8wDataSplines.root"));
    DataSample* mc8w   = new DataSample(file_exists(MCLocation+"/run8wMCsplines.root"));
    //DataSample* sand8w = new DataSample(file_exists(MCLocation+"/run8wMCsplines.root"));
    SampleGroup run8w("run8w");
    run8w.AddDataSample(data8w);
    run8w.AddMCSample("magnet", mc8w);
    //run8w.AddMCSample("sand", sand8w);
    exp->AddSampleGroup("run8w", run8w); 
  }

  if (prod->Contains(" 9")) {
    DataSample* data9 = new DataSample(file_exists(DataLocation+"/run9DataSplines.root"));
    DataSample* mc9   = new DataSample(file_exists(MCLocation+"/run9MCsplines.root"));
    //DataSample* sand9 = new DataSample(file_exists(MCLocation+"/run9MCsplines.root"));
    SampleGroup run9("run9");
    run9.AddDataSample(data9);
    run9.AddMCSample("magnet", mc9);
    //run9.AddMCSample("sand", sand9);
    exp->AddSampleGroup("run9", run9);
  }

  // Associate the Experiment class to the psyche AnalysisManager
  _man.SetExperiment(exp);

  // Print the POT for the current selections
  PrintPOT();

  // Psyche supports preloading (load events into memory before looping over)
  // This is faster than processing from .root file but takes a lot of RAM
  // Needs to be set differently depending on if we're getting detector+flux systematics from psyche (== !simple) or if we're using a pre-built covariance (== simple)
  // (false, true) means (!preloadData, preloadMC)
  if (!simple) {
    _man.PreloadEvents();
  } else {
    _man.PreloadEvents(false, false);
  }


  // Comment left for hilariousness:
  // "let's run roughshod over Anselmo! Disable all selections so we can enable them when we want." :D
  // We disable all selections in psyche so we can add whatever we want later on by samplePDFND2019::addSelection
  // Get the selections that are enabled in the AnalysisManager (initialised by psycheSteering/parameters/*.dat
  std::vector<SelectionBase*> selec = _man.sel().GetSelections();
  // Turn them off
  for (std::vector<SelectionBase*>::iterator it = selec.begin(); it != selec.end(); it++) {
    (*it)->Disable();
  }

  // Expand the TObjArrays to be all the psyche samples
  // Could expand only the enabled samples but is OK as long as we initalise unwanted to NULL
  samplepdfs->Expand(SampleId::kNSamples);
  datapdfs->Expand(SampleId::kNSamples);
  // And initialise them to NULL
  // Initialize the arrays that carry the number of dimensions for each psyche sample
  for (int i = 0; i < SampleId::kNSamples; i++) {
    samplepdfs->AddAt(NULL, i);
    datapdfs->AddAt(NULL, i);
    ndims[i] = 0;
  }

  //NowTalk();
}


// ***************************************************************************
// Load up the samples and binning specified by the user
void samplePDFND2019::LoadSamples() {
// ***************************************************************************

  // Can turn on selections from config file
  // neutrino mode:
  // FGD1 CC0pi, FGD1 CC1pi, FGD1 CCOther
  // FGD2 CC0pi, FGD2 CC1pi, FGD2 CCOther
  // anti-neutrino mode:
  // FGD1 CCQE, FGD1 nCCQE, FGD1 numuBCKG CCQE, FGD1 numuBCKG nCCQE
  // FGD2 CCQE, FGD2 nCCQE, FGD2 numuBCKG CCQE, FGD2 numuBCKG nCCQE
  // Probably easiest to loop over a vector of strings

  // ========================================================================
  // Now load up the samples
  std::vector<std::string> Sel_Vec = FitManager->getNDsel();
  // Set the default if empty
  if (Sel_Vec.empty()) {
    std::cout << "******************************" << std::endl;
    std::cout << "No selections specified in config file, adding all!" << std::endl;
    std::cout << "******************************" << std::endl;
    // FGD1 neutrino selections
    Sel_Vec.push_back("FGD1 nu CC0pi");
    Sel_Vec.push_back("FGD1 nu CC1pi");
    Sel_Vec.push_back("FGD1 nu CCOth");

    // FGD2 neutrino selections
    Sel_Vec.push_back("FGD2 nu CC0pi");
    Sel_Vec.push_back("FGD2 nu CC1pi");
    Sel_Vec.push_back("FGD2 nu CCOth");

    /*    // FGD1 anti-neutrino selections
    Sel_Vec.push_back("FGD1 anu CC1Trk");
    Sel_Vec.push_back("FGD1 anu CCNTrk");

    // FGD1 neutrino in anti-neutrino selections
    Sel_Vec.push_back("FGD1 anu nu CC1Trk");
    Sel_Vec.push_back("FGD1 anu nu CCNTrk");

    // FGD1 anti-neutrino selections
    Sel_Vec.push_back("FGD2 anu CC1Trk");
    Sel_Vec.push_back("FGD2 anu CCNTrk");

    // FGD1 neutrino in anti-neutrino selections
    Sel_Vec.push_back("FGD2 anu nu CC1Trk");
    Sel_Vec.push_back("FGD2 anu nu CCNTrk");
    */
    // CC0pi, CC1pi, CCoth anti-neutrino selections
    
    Sel_Vec.push_back("FGD1 anu CC0pi");
    Sel_Vec.push_back("FGD1 anu CC1pi");
    Sel_Vec.push_back("FGD1 anu CCOth");

    Sel_Vec.push_back("FGD1 anu nu CC0pi");
    Sel_Vec.push_back("FGD1 anu nu CC1pi");
    Sel_Vec.push_back("FGD1 anu nu CCOth");

    Sel_Vec.push_back("FGD2 anu CC0pi");
    Sel_Vec.push_back("FGD2 anu CC1pi");
    Sel_Vec.push_back("FGD2 anu CCOth");

    Sel_Vec.push_back("FGD2 anu nu CC0pi");
    Sel_Vec.push_back("FGD2 anu nu CC1pi");
    Sel_Vec.push_back("FGD2 anu nu CCOth");
    

    // Nue
    /*
    Sel_Vec.push_back("FGD1 nu E CC");
    Sel_Vec.push_back("FGD2 nu E CC");
    Sel_Vec.push_back("FGD1 anu E CC");
    Sel_Vec.push_back("FGD2 anu E CC");
    Sel_Vec.push_back("FGD1 anu nu E CC");
    Sel_Vec.push_back("FGD2 anu nu E CC");
    Sel_Vec.push_back("FGD1 nu anu E CC");
    Sel_Vec.push_back("FGD2 nu anu E CC");
    */

    // gamma
    /*
    Sel_Vec.push_back("FGD1 nu gamma");
    Sel_Vec.push_back("FGD2 nu gamma");
    Sel_Vec.push_back("FGD1 anu gamma");
    Sel_Vec.push_back("FGD1 anu gamma");
    */
  }

  // =========================================================
  // Add selections with specified binning
  std::cout << "------------------------------" << std::endl;
  std::cout << "Adding selections at ND280" << std::endl;

  // FGD1 nu CCinc
  if (std::find(Sel_Vec.begin(), Sel_Vec.end(), std::string("FGD1 nu CCinc")) != Sel_Vec.end()) {
    addSelection(SampleId::kFGD1NuMuCC);
  }

  // FGD2 nu CCinc
  if (std::find(Sel_Vec.begin(), Sel_Vec.end(), std::string("FGD2 nu CCinc")) != Sel_Vec.end()) {
    addSelection(SampleId::kFGD2NuMuCC);
  }

  // FGD1 neutrino CC0pi
  if (std::find(Sel_Vec.begin(), Sel_Vec.end(), std::string("FGD1 nu CC0pi")) != Sel_Vec.end()) {
    addSelection(SampleId::kFGD1NuMuCC0Pi);
  }

  // FGD1 neutrino CC1pi
  if (std::find(Sel_Vec.begin(), Sel_Vec.end(), std::string("FGD1 nu CC1pi")) != Sel_Vec.end()) {
    addSelection(SampleId::kFGD1NuMuCC1Pi);
  }

  // FGD1 neutrino CCNpi
  if (std::find(Sel_Vec.begin(), Sel_Vec.end(), std::string("FGD1 nu CCOth")) != Sel_Vec.end()) {
    addSelection(SampleId::kFGD1NuMuCCOther);
  }

  // FGD1 neutrino in anti-neutrino CC1Trk
  if (std::find(Sel_Vec.begin(), Sel_Vec.end(), std::string("FGD1 anu CC1Trk")) != Sel_Vec.end()) {
    addSelection(SampleId::kFGD1AntiNuMuCC1Track);
  }

  // FGD1 neutrino in anti-neutrino CCNTrk
  if (std::find(Sel_Vec.begin(), Sel_Vec.end(), std::string("FGD1 anu CCNTrk")) != Sel_Vec.end()) {
    addSelection(SampleId::kFGD1AntiNuMuCCNTracks);
  }

  // FGD1 anti-neutrino CC0pi
  if (std::find(Sel_Vec.begin(), Sel_Vec.end(), std::string("FGD1 anu CC0pi")) != Sel_Vec.end()) {
    addSelection(SampleId::kFGD1AntiNuMuCC0Pi);
  }

  // FGD1 anti-neutrino CC1pi
  if (std::find(Sel_Vec.begin(), Sel_Vec.end(), std::string("FGD1 anu CC1pi")) != Sel_Vec.end()) {
    addSelection(SampleId::kFGD1AntiNuMuCC1Pi);
  }

  // FGD1 anti-neutrino CCNpi
  if (std::find(Sel_Vec.begin(), Sel_Vec.end(), std::string("FGD1 anu CCOth")) != Sel_Vec.end()) {
    addSelection(SampleId::kFGD1AntiNuMuCCOther);
  }

  // FGD1 neutrino in anti-neutrino CC1Trk
  if (std::find(Sel_Vec.begin(), Sel_Vec.end(), std::string("FGD1 anu nu CC1Trk")) != Sel_Vec.end()) {
    addSelection(SampleId::kFGD1NuMuBkgInAntiNuModeCC1Track);
  }

  // FGD1 neutrino in anti-neutrino CCNTrk
  if (std::find(Sel_Vec.begin(), Sel_Vec.end(), std::string("FGD1 anu nu CCNTrk")) != Sel_Vec.end()) {
    addSelection(SampleId::kFGD1NuMuBkgInAntiNuModeCCNTracks);
  }

  // FGD1 neutrino in anti-neutrino CC0pi
  if (std::find(Sel_Vec.begin(), Sel_Vec.end(), std::string("FGD1 anu nu CC0pi")) != Sel_Vec.end()) {
    addSelection(SampleId::kFGD1NuMuBkgInAntiNuModeCC0Pi);
  }

  // FGD1 neutrino in anti-neutrino CC1pi
  if (std::find(Sel_Vec.begin(), Sel_Vec.end(), std::string("FGD1 anu nu CC1pi")) != Sel_Vec.end()) {
    addSelection(SampleId::kFGD1NuMuBkgInAntiNuModeCC1Pi);
  }

  // FGD1 neutrino in anti-neutrino CCNpi
  if (std::find(Sel_Vec.begin(), Sel_Vec.end(), std::string("FGD1 anu nu CCOth")) != Sel_Vec.end()) {
    addSelection(SampleId::kFGD1NuMuBkgInAntiNuModeCCOther);
  }

  // FGD2 neutrino CC0pi
  if (std::find(Sel_Vec.begin(), Sel_Vec.end(), std::string("FGD2 nu CC0pi")) != Sel_Vec.end()) {
    addSelection(SampleId::kFGD2NuMuCC0Pi);
  }

  // FGD2 neutrino CC1pi
  if (std::find(Sel_Vec.begin(), Sel_Vec.end(), std::string("FGD2 nu CC1pi")) != Sel_Vec.end()) {
    addSelection(SampleId::kFGD2NuMuCC1Pi);
  }

  // FGD2 neutrino CCoth
  if (std::find(Sel_Vec.begin(), Sel_Vec.end(), std::string("FGD2 nu CCOth")) != Sel_Vec.end()) {
    addSelection(SampleId::kFGD2NuMuCCOther);
  }

  // FGD2 anti-neutrino CC1Trk
  if (std::find(Sel_Vec.begin(), Sel_Vec.end(), std::string("FGD2 anu CC1Trk")) != Sel_Vec.end()) {
    addSelection(SampleId::kFGD2AntiNuMuCC1Track);
  }

  // FGD2 anti-neutrino CCNTrk
  if (std::find(Sel_Vec.begin(), Sel_Vec.end(), std::string("FGD2 anu CCNTrk")) != Sel_Vec.end()) {
    addSelection(SampleId::kFGD2AntiNuMuCCNTracks);
  }

  // FGD2 anti-neutrino CC0pi
  if (std::find(Sel_Vec.begin(), Sel_Vec.end(), std::string("FGD2 anu CC0pi")) != Sel_Vec.end()) {
    addSelection(SampleId::kFGD2AntiNuMuCC0Pi);
  }

  // FGD2 anti-neutrino CC1pi
  if (std::find(Sel_Vec.begin(), Sel_Vec.end(), std::string("FGD2 anu CC1pi")) != Sel_Vec.end()) {
    addSelection(SampleId::kFGD2AntiNuMuCC1Pi);
  }

  // FGD2 anti-neutrino CCNpi
  if (std::find(Sel_Vec.begin(), Sel_Vec.end(), std::string("FGD2 anu CCOth")) != Sel_Vec.end()) {
    addSelection(SampleId::kFGD2AntiNuMuCCOther);
  }

  // FGD2 neutrino in anti-neutrino CC1Trk
  if (std::find(Sel_Vec.begin(), Sel_Vec.end(), std::string("FGD2 anu nu CC1Trk")) != Sel_Vec.end()) {
    addSelection(SampleId::kFGD2NuMuBkgInAntiNuModeCC1Track);
  }

  // FGD2 neutrino in anti-neutrino CCNTrk
  if (std::find(Sel_Vec.begin(), Sel_Vec.end(), std::string("FGD2 anu nu CCNTrk")) != Sel_Vec.end()) {
    addSelection(SampleId::kFGD2NuMuBkgInAntiNuModeCCNTracks);
  }

  // FGD2 neutrino in anti-neutrino CC0pi
  if (std::find(Sel_Vec.begin(), Sel_Vec.end(), std::string("FGD2 anu nu CC0pi")) != Sel_Vec.end()) {
    addSelection(SampleId::kFGD2NuMuBkgInAntiNuModeCC0Pi);
  }

  // FGD2 neutrino in anti-neutrino CC1pi
  if (std::find(Sel_Vec.begin(), Sel_Vec.end(), std::string("FGD2 anu nu CC1pi")) != Sel_Vec.end()) {
    addSelection(SampleId::kFGD2NuMuBkgInAntiNuModeCC1Pi);
  }

  // FGD2 neutrino in anti-neutrino CCNpi
  if (std::find(Sel_Vec.begin(), Sel_Vec.end(), std::string("FGD2 anu nu CCOth")) != Sel_Vec.end()) {
    addSelection(SampleId::kFGD2NuMuBkgInAntiNuModeCCOther);
  }

  // Check that at least one sample is enabled
  unsigned int SampleCounter = 0;
  for (int i = 0; i < SampleId::kNSamples; ++i) {
    if (samplepdfs->At(i) != NULL) SampleCounter++;
  }
  if (SampleCounter != Sel_Vec.size()) {
    std::cerr << "You gave " << Sel_Vec.size() << " selections in the config file but psyche only matched " << SampleCounter << " samples!" << std::endl;
    std::cerr << "Config file selections: " << std::endl;
    for (std::vector<std::string>::iterator it = Sel_Vec.begin(); it != Sel_Vec.end(); ++it) {
      std::cerr << "   " << *it << std::endl;
    }
    std::cerr << "Matched selections: " << std::endl;
    for (int i = 0; i < SampleId::kNSamples; ++i) {
      if (samplepdfs->At(i) == NULL) continue;
      std::cerr << "   " << samplepdfs->At(i)->GetName() << std::endl;
    }

    std::cerr << "Look at addSelection in " << __FILE__ << std::endl;
    std::cerr << "I live at " << __FILE__ << ":" << __LINE__ << std::endl;
    throw;
  }

  maxBins=0;
  polybins = new TH2PolyBin**[SampleId::kNSamples];
  for(__int__ i = 0; i<SampleId::kNSamples; i++)
    {
      if (samplepdfs->At(i) == NULL) continue;
      if(maxBins<((TH2Poly*)(samplepdfs->At(i)))->GetNumberOfBins())
	{
	  maxBins=((TH2Poly*)(samplepdfs->At(i)))->GetNumberOfBins();
	}
      polybins[i] = new TH2PolyBin*[maxBins];
      for(__int__ j = 0; j<((TH2Poly*)(samplepdfs->At(i)))->GetNumberOfBins(); j++)
        {
          polybins[i][j] = (TH2PolyBin*)(((TH2Poly*)(samplepdfs->At(i)))->GetBins()->At(j)->Clone());
        }
      for(__int__ j = -9; j<0; j++)//now account for overflow. TH2Polys overflow bins are -1 to -9
        {
          polybins[i][maxBins+j] = (TH2PolyBin*)(((TH2Poly*)(samplepdfs->At(i)))->GetBins()->At(j)->Clone());
        }
    }

  ebtemplates = new TH2Poly****[SampleId::kNSamples];
  for(__int__ i = 0; i<SampleId::kNSamples; i++)
    {
      if (samplepdfs->At(i) == NULL) continue;
      ebtemplates[i] = new TH2Poly***[2];
      for (int j=0; j<2; j++){//loop over target
        ebtemplates[i][j] = new TH2Poly**[4];
        for (int k=0; k<4; k++){//loop over species
	  ebtemplates[i][j][k] = new TH2Poly*[73];
	  for (int l=0; l<73; l++){//loop over species
	    ebtemplates[i][j][k][l] = (TH2Poly*)(samplepdfs->At(i)->Clone());
	  }
	}
      }
    }

  const int var1numubins=73;
  double var1numubinrange[var1numubins+1];
  for(int i=0;i<61; i++){
    var1numubinrange[i] = 3.0/double(60)*i;
  }
  var1numubinrange[61] = 3.25;
  var1numubinrange[62] = 3.5;
  var1numubinrange[63] = 3.75;
  var1numubinrange[64] = 4.0;
  var1numubinrange[65] = 4.5;
  var1numubinrange[66] = 5.0;
  var1numubinrange[67] = 5.5;
  var1numubinrange[68] = 6.0;
  var1numubinrange[69] = 7.0;
  var1numubinrange[70] = 8.0;
  var1numubinrange[71] = 9.0;
  var1numubinrange[72] = 10.0;
  var1numubinrange[73] = 30.0;

      //Erec version
  int var1nuebins=25;
  double var1nuebinrange[var1nuebins+1];
  for(int iBin=0;iBin<var1nuebins+1;iBin++){
    var1nuebinrange[iBin]=1.25/double(var1nuebins)*iBin;
  }

  NumuErec = new TH1D("numuerec","numuerec",var1numubins,var1numubinrange);
  NumubarErec = new TH1D("numubarerec","numubarerec",var1numubins,var1numubinrange);
  NueErec = new TH1D("nueerec","nueerec",var1nuebins,var1nuebinrange);
  NuebarErec = new	TH1D("nuebarerec","nuebarerec",var1nuebins,var1nuebinrange);
} // end LoadSamples

// ***************************************************************************
// Set the datapdfs to the MC-pdfs and treat them as Asimov data
// This might require throwing the starting parameters of the Markov Chain to avoid a stationary chain at the first n steps (n might be as large as 10k!)
void samplePDFND2019::setAsimovFakeData(bool CustomReWeight) {
  // ***************************************************************************

  if (HaveIRandomStart) {
    std::cerr << "samplePDFND2019 has been set to random start and been reweighted to it" << std::endl;
    std::cerr << "And now being asked to set Asimov data to the random start... I think this is wrong!" << std::endl;
    std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
    throw;
  }

  // Check if the relevant covariances have been set
  CheckCovariances();

  // Set Asimov parameters
  // Could use config here instead but will hard-code for now
  if (CustomReWeight) {
    xs2015->setParProp(6, 1.0);
    xs2015->setParProp(7, 1.0);
    double *fake = NULL;
    reweight(fake);
    // Reset the fParCurr and nominal parameters to nominal after the reweight
    std::cout << "Set samplePDFND2019 asimov with cross-section parameters: " << std::endl;
    xs2015->printNominalCurrProp();
    std::cout << "Now resetting..." << std::endl;
    xs2015->setParameters();
  }

  // Loop over all the psyche samples
  for (int i = 0; i < SampleId::kNSamples; ++i) {

    // Move to next sample if we haven't enabled this one
    if (samplepdfs->At(i) == NULL || datapdfs->At(i) == NULL) continue;

    // Can't see why we wouldn't run 2D
    if (ndims[i] != 2) {
      std::cerr << "Can not set Asimov data for other than 2D histograms for " << SampleId::ConvertSample(SampleId::SampleEnum(i)) << std::endl;
      std::cerr << "Simply a matter of implementation; to edit this see " << __FILE__ << ":" << __LINE__ << std::endl;
      throw;
    }

    // Strip out the whitespaces and replace them with underscores in the psyche names
    // Will be used in Clone
    std::string temp = SampleId::ConvertSample(SampleId::SampleEnum(i));
    while (temp.find(" ") != std::string::npos) {
      temp.replace(temp.find(" "), 1, std::string("_"));
    }

    // Just clone the MC for the simple setAsimovFakeData()
    TH2Poly* MCpdf = (TH2Poly*)(getPDF(i))->Clone(temp.c_str());
    // Pass the pointer to addData
    addData(MCpdf, i);
  } // end for loop

  std::cout << "------------------------------------------" << std::endl;
  std::cout << "Have set Asimov fake-data to reweighted MC" << std::endl;
  std::cout << "With cross-section settings: " << std::endl;
  xs2015->printNominalCurrProp();
  //flux->printNominalCurrProp();
  //detsimple->printNominalCurrProp();
  std::cout << "------------------------------------------" << std::endl;

  // Print the post-Asimov rates
  printRates();
} // end setAsimovFakeData()

// ***************************************************************************
// asdf
void samplePDFND2019::setAsimovFakeData_FromFile(std::string &FileName) {
  // ***************************************************************************

  if (HaveIRandomStart) {
    std::cerr << "samplePDFND2019 has been set to random start and been reweighted to it" << std::endl;
    std::cerr << "And now being asked to set Asimov data to the random start... I think this is wrong!" << std::endl;
    std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
    throw;
  }

  // Check if the relevant covariances have been set
  CheckCovariances();

  // Open the TFile
  TFile *AsimovFile = new TFile(FileName.c_str(), "OPEN");
  if (!AsimovFile->IsOpen()) {
    std::cerr << "************" << std::endl;
    std::cerr << "Provided Asimov file " << FileName << " does not exist, exiting" << std::endl;
    std::cerr << "************" << std::endl;
    throw;
  }
  AsimovFile->ls();

  // Loop over all the psyche samples
  for (int i = 0; i < SampleId::kNSamples; ++i) {

    // Move to next sample if we haven't enabled this one
    if (samplepdfs->At(i) == NULL || datapdfs->At(i) == NULL) continue;

    // Can't see why we wouldn't run 2D
    if (ndims[i] != 2) {
      std::cerr << "Can not set Asimov data for other than 2D histograms for " << SampleId::ConvertSample(SampleId::SampleEnum(i)) << std::endl;
      std::cerr << "Simply a matter of implementation; to edit this see " << __FILE__ << ":" << __LINE__ << std::endl;
      throw;
    }

    // Strip out the whitespaces and replace them with underscores in the psyche names
    // Will be used in Clone
    std::string temp = SampleId::ConvertSample(SampleId::SampleEnum(i));
    while (temp.find(" ") != std::string::npos) {
      temp.replace(temp.find(" "), 1, std::string("_"));
    }

    // Clone the MC from the TFile
    TH2Poly* MCpdf = (TH2Poly*)(AsimovFile->Get((std::string("MC_")+temp).c_str())->Clone());

    if (MCpdf == NULL) {
      std::cerr << "Did not file input Asimov in " << AsimovFile << " for " << SampleId::ConvertSample(SampleId::SampleEnum(i)) << std::endl;
      throw;
    }
    MCpdf->SetDirectory(0);

    // Pass the pointer to addData
    addData(MCpdf, i);
  } // end for loop

  std::cout << "------------------------------------------" << std::endl;
  std::cout << "Have set Asimov fake-data to reweighted MC" << std::endl;
  std::cout << "With cross-section settings: " << std::endl;
  xs2015->printNominalCurrProp();
  //flux->printNominalCurrProp();
  //detsimple->printNominalCurrProp();
  std::cout << "------------------------------------------" << std::endl;

  // Print the post-Asimov rates
  printRates();

  delete AsimovFile;
}


// ***************************************************************************
// Similar to samplePDFND2019::setAsimovFakeData()
// Throws the parameters that the fake-data is generated at first so that fake-data is not the nominal parameter set
// For this to work we need the associated xsec, flux and ND280 covariances
void samplePDFND2019::setAsimovFakeDataThrow() {
  // ***************************************************************************

  std::cout << "Setting Asimov fake data to some neat permutation of MC" << std::endl;

  // Check if the relevant covariances have been set
  CheckCovariances();

  // Don't want nominal values so set bool to false
  // Throw xsec
  xs2015->throwNominal(false);
  xs2015->setParameters();
  xs2015->printNominalCurrProp();
  // Throw detector
  detsimple->throwNominal(false);
  detsimple->setParameters();
  // Don't necessarily wan't to print this bad-boy (580 parameters)
  //detsimple->printNominalCurrProp();
  // Throw flux
  flux->throwNominal(false);
  flux->setParameters();
  // Don't necessarily wan't to print this bad-boy (100 parameters)
  //flux->printNominalCurrProp();

  // Now reweight the samplePDF with the new thrown parameters
  double *fake = 0;
  reweight(fake);

  // Loop over all the psyche samples
  for (int i = 0; i < SampleId::kNSamples; ++i) {

    // Move to next sample if we haven't enabled this one
    if (samplepdfs->At(i) == NULL || datapdfs->At(i) == NULL) continue;

    if (ndims[i] != 2) {
      std::cerr << "Can not set Asimov data for other than 2D histograms for " << SampleId::ConvertSample(SampleId::SampleEnum(i)) << std::endl;
      std::cerr << "Simply a matter of implementation; to edit this see " << __FILE__ << ":" << __LINE__ << std::endl;
      throw;
    }

    // Strip out the whitespaces and replace them with underscores in the psyche names
    // Will be used in Clone
    std::string temp = SampleId::ConvertSample(SampleId::SampleEnum(i));
    while (temp.find(" ") != std::string::npos) {
      temp.replace(temp.find(" "), 1, std::string("_"));
    }

    // Just clone the MC for the simple setAsimovFakeData()
    TH2Poly* MCpdf = (TH2Poly*)(getPDF(i))->Clone(temp.c_str());
    // Pass the pointer to addData
    addData(MCpdf, i);

  } // end for loop

  // Let's throw the parameters so that MC != Asimov data
  // Don't want nominal values so set bool to false
  // Throw xsec
  xs2015->throwNominal(false);
  xs2015->setParameters();
  xs2015->printNominalCurrProp();
  // Throw detector
  detsimple->throwNominal(false);
  detsimple->setParameters();
  // Don't necessarily wan't to print this bad-boy
  //detsimple->printNominalCurrProp();
  // Throw flux
  flux->throwNominal(false);
  flux->setParameters();
  // Don't necessarily wan't to print this bad-boy
  //flux->printNominalCurrProp();
  reweight(fake);

  std::cout << "------------------------------------------" << std::endl;
  std::cout << "Have set Asimov fake-data to varied reweighted MC" << std::endl;
  std::cout << "------------------------------------------" << std::endl;

  // Print the post-Asimov rates
  printRates();
} // end setAsimovFakeData()


// ***************************************************************************
// Setup the binning for a given psyche selection
// Not really used but could be called in addSelection...
void samplePDFND2019::SetupBinning(SampleId::SampleEnum Selection, std::vector<double> & BinningX, std::vector<double> & BinningY) {
  // ***************************************************************************

  // LeptonMomentum and LeptonCosTheta
  if (kinvars[Selection][0] == kLeptonMomentum && kinvars[Selection][1] == kLeptonCosTheta) {

    switch (Selection) {
      // FGD1 and FGD2 share binning
      // CC0pi and CCOther have same binning too
      // Also put in the inclusive selection for debugging
    case (SampleId::kFGD1NuMuCC):
    case (SampleId::kFGD2NuMuCC):
    case (SampleId::kFGD1NuMuCC0Pi):
    case (SampleId::kFGD2NuMuCC0Pi):
      {
	const double BinArray_x[] = {0., 200., 300., 400., 450., 500., 550., 600., 650., 700., 750., 800., 850., 900., 950., 1000., 1050., 1100., 1200., 1300., 1400., 1500., 1600., 1700., 1800., 2000., 2500., 3000., 5000., 30000.};
	const double BinArray_y[] = {-1., 0.5, 0.6, 0.7, 0.76, 0.78, 0.8, 0.83, 0.85, 0.88, 0.89, 0.9, 0.91, 0.92, 0.925, 0.93, 0.935, 0.94, 0.945, 0.95, 0.955, 0.96, 0.965, 0.97, 0.975, 0.98, 0.985, 0.99, 0.995, 1.};
	BinningX = MakeVector(BinArray_x);
	BinningY = MakeVector(BinArray_y);
	break;
      }

    // FGD1 and FGD2 share binning
    // Different to CC0pi for CC1pi
    case (SampleId::kFGD1NuMuCC1Pi):
    case (SampleId::kFGD2NuMuCC1Pi):
      {
	const double BinArray_x[] = {0., 300., 350., 400., 500., 600., 650., 700., 750., 800., 900., 1000., 1100., 1200., 1500., 2000., 3000., 5000., 30000.};
	const double BinArray_y[] = {-1., 0.6, 0.7, 0.8, 0.85, 0.88, 0.9, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 0.995, 1.};
	BinningX = MakeVector(BinArray_x);
	BinningY = MakeVector(BinArray_y);
	break;
      }

    case (SampleId::kFGD1NuMuCCOther):
    case (SampleId::kFGD2NuMuCCOther):
      {
	const double BinArray_x[] = {0., 300., 400., 500., 600., 650., 700., 750., 800., 900., 1000., 1100., 1250., 1500., 1750., 2000., 3000., 5000., 30000.};
	const double BinArray_y[] = {-1., 0.6, 0.7, 0.76, 0.8, 0.85, 0.88, 0.89, 0.9, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 0.995, 1.};
	BinningX = MakeVector(BinArray_x);
	BinningY = MakeVector(BinArray_y);
	break;
      }

    // Anti-neutrino binning
    case (SampleId::kFGD1AntiNuMuCC1Track):
    case (SampleId::kFGD2AntiNuMuCC1Track):
      {
	const double BinArray_x[] = {0., 400., 500., 600., 700., 800., 900., 1100., 1400., 2000., 10000.};
	const double BinArray_y[] = {-1., 0.6, 0.7, 0.8, 0.85, 0.88, 0.91, 0.93, 0.95, 0.96, 0.97, 0.98, 0.99, 1.};
	std::vector<double> AnuMom_qe = FitManager->getDetBin_anuMom_qe();
	std::vector<double> AnuTh_qe = FitManager->getDetBin_anuTh_qe();
	if (AnuMom_qe.empty()) {
	  BinningX = MakeVector(BinArray_x);
	} else {
	  std::cout << "Using custom pmu binning" << std::endl;
	  BinningX = AnuMom_qe;
	}
	if (AnuTh_qe.empty()) {
	  BinningY = MakeVector(BinArray_y);
	} else {
	  std::cout << "Using custom costhmu binning" << std::endl;
	  BinningY = AnuTh_qe;
	}
	break;
      }

    case (SampleId::kFGD1AntiNuMuCCNTracks):
    case (SampleId::kFGD2AntiNuMuCCNTracks):
      {
	const double BinArray_x[] = {0., 700., 950., 1200., 1500., 2000., 3000., 10000.};
	const double BinArray_y[] = {-1., 0.75, 0.85, 0.88, 0.91, 0.93, 0.95, 0.96, 0.97, 0.98, 0.99, 1.};
	std::vector<double> AnuMom_nqe = FitManager->getDetBin_anuMom_nqe();
	std::vector<double> AnuTh_nqe = FitManager->getDetBin_anuTh_nqe();
	if (AnuMom_nqe.empty()) {
	  BinningX = MakeVector(BinArray_x);
	} else {
	  std::cout << "Using custom pmu binning" << std::endl;
	  BinningX = AnuMom_nqe;
	}
	if (AnuTh_nqe.empty()) {
	  BinningY = MakeVector(BinArray_y);
	} else {
	  std::cout << "Using custom costhmu binning" << std::endl;
	  BinningY = AnuTh_nqe;
	}
	break;
      }

    case (SampleId::kFGD1NuMuBkgInAntiNuModeCC1Track):
    case (SampleId::kFGD2NuMuBkgInAntiNuModeCC1Track):
      {
	const double BinArray_x[] = {0., 400., 600., 800., 1100., 2000., 10000.};
	const double BinArray_y[] = {-1., 0.7, 0.8, 0.85, 0.90, 0.93, 0.95, 0.96, 0.97, 0.98, 0.99, 1.};
	std::vector<double> Anu_NumuMom_qe = FitManager->getDetBin_numu_anuMom_qe();
	std::vector<double> Anu_NumuTh_qe  = FitManager->getDetBin_numu_anuTh_qe();
	if (Anu_NumuMom_qe.empty()) {
	  BinningX = MakeVector(BinArray_x);
	} else {
	  std::cout << "Using custom pmu binning" << std::endl;
	  BinningX = Anu_NumuMom_qe;
	}
	if (Anu_NumuMom_qe.empty()) {
	  BinningY = MakeVector(BinArray_y);
	} else {
	  std::cout << "Using custom costhmu binning" << std::endl;
	  BinningY = Anu_NumuTh_qe;
	}
	break;
      }

    case (SampleId::kFGD1NuMuBkgInAntiNuModeCCNTracks):
    case (SampleId::kFGD2NuMuBkgInAntiNuModeCCNTracks):
      {
	const double BinArray_x[] = {0., 500., 700., 1000., 1250., 1500., 2000., 3000., 10000.};
	const double BinArray_y[] = {-1., 0.7, 0.8, 0.85, 0.90, 0.93, 0.95, 0.96, 0.97, 0.98, 0.99, 1.};
	std::vector<double> Anu_NumuMom_nqe= FitManager->getDetBin_numu_anuMom_nqe();
	std::vector<double> Anu_NumuTh_nqe = FitManager->getDetBin_numu_anuTh_nqe();
	if (Anu_NumuMom_nqe.empty()) {
	  BinningX = MakeVector(BinArray_x);
	} else {
	  std::cout << "Using custom pmu binning" << std::endl;
	  BinningX = Anu_NumuMom_nqe;
	}
	if (Anu_NumuMom_nqe.empty()) {
	  BinningY = MakeVector(BinArray_y);
	} else {
	  std::cout << "Using custom costhmu binning" << std::endl;
	  BinningY = Anu_NumuTh_nqe;
	}
	break;
      }

    case (SampleId::kFGD1AntiNuMuCC0Pi):
    case (SampleId::kFGD2AntiNuMuCC0Pi):
      {
	const double BinArray_x[] = {0., 300., 400., 500., 550., 600., 650., 700., 750., 800., 900., 1000., 1100., 1200., 1500., 2000., 4000., 30000.};
	const double BinArray_y[] = {-1., 0.6, 0.7, 0.8, 0.85, 0.9, 0.92, 0.93, 0.94, 0.95, 0.96, 0.965, 0.97, 0.975, 0.98, 0.985, 0.99, 0.995, 1.};
	BinningX = MakeVector(BinArray_x);
	BinningY = MakeVector(BinArray_y);
	break;
      }

    case (SampleId::kFGD1AntiNuMuCC1Pi):
    case (SampleId::kFGD2AntiNuMuCC1Pi):
      {
	const double BinArray_x[] = {0., 500., 700., 900., 1300., 2500., 30000.};
	const double BinArray_y[] = {-1, 0.7, 0.8, 0.9, 0.94, 0.96, 0.98, 0.99, 1};
	BinningX = MakeVector(BinArray_x);
	BinningY = MakeVector(BinArray_y);
	break;
      }

    case (SampleId::kFGD1AntiNuMuCCOther):
    case (SampleId::kFGD2AntiNuMuCCOther):
      {
	const double BinArray_x[] = {0., 600., 800., 1000., 1250., 1500., 2000., 4000., 30000.};
	const double BinArray_y[] = {-1., 0.7, 0.8, 0.85, 0.9, 0.93, 0.95, 0.97, 0.98, 0.99, 1.};
	BinningX = MakeVector(BinArray_x);
	BinningY = MakeVector(BinArray_y);
	break;
      }

    case (SampleId::kFGD1NuMuBkgInAntiNuModeCC0Pi):
    case (SampleId::kFGD2NuMuBkgInAntiNuModeCC0Pi):
      {
	const double BinArray_x[] = {0., 300., 500., 700., 800., 900., 1250., 1500., 2000., 4000., 30000.};
	const double BinArray_y[] = {-1., 0.7, 0.8, 0.85, 0.88, 0.9, 0.92, 0.94, 0.96, 0.97, 0.98, 0.99, 1.};
	BinningX = MakeVector(BinArray_x);
	BinningY = MakeVector(BinArray_y);
	break;
      }

    case (SampleId::kFGD1NuMuBkgInAntiNuModeCC1Pi):
    case (SampleId::kFGD2NuMuBkgInAntiNuModeCC1Pi):
      {
	const double BinArray_x[] = {0., 600., 800., 1500., 30000.};
	const double BinArray_y[] = {-1, 0.7, 0.8, 0.86, 0.9, 0.94, 0.96, 0.97, 0.98, 0.99, 1};
	BinningX = MakeVector(BinArray_x);
	BinningY = MakeVector(BinArray_y);
	break;
      }

    case (SampleId::kFGD1NuMuBkgInAntiNuModeCCOther):
    case (SampleId::kFGD2NuMuBkgInAntiNuModeCCOther):
      {
	const double BinArray_x[] = {0., 600., 1000., 1250., 2000., 4000., 30000.};
	const double BinArray_y[] = {-1., 0.7, 0.8, 0.86, 0.9, 0.93, 0.95, 0.97, 0.99, 1.};
	BinningX = MakeVector(BinArray_x);
	BinningY = MakeVector(BinArray_y);
	break;
      }

    default:
      std::cerr << "You gave me a selection I don't have a binning for, whaaa?" << std::endl;
      std::cerr << "Selection: " << Selection << " = " << SampleId::ConvertSample(Selection) << std::endl;
      throw;
      break;
    }

    // PionMomentum vs PionCosTheta
  } else if (kinvars[Selection][0] == kPionMomentum && kinvars[Selection][1] == kPionCosTheta) {

    const double BinArray_x[] = {0., 100., 200., 300., 400., 500., 600., 700., 800., 900., 1000.};
    const double BinArray_y[] = {-1., 0.7, 0.8, 0.85, 0.90, 0.92, 0.94, 0.96, 0.97, 0.98, 0.99, 1.};
    BinningX = MakeVector(BinArray_x);
    BinningY = MakeVector(BinArray_y);

    // For when we want Q2QERec and ErecQE as x and y variables
  } else if (kinvars[Selection][0] == kQ2QE && kinvars[Selection][1] == kErecQE) {

    // Use the same binning for all the samples
    // This is really only for post-fit plots: we don't fit in Q2QErec or ErecQE (yet, thank God)
    const double q2_a[] = {0.0, 0.05, 0.1, 0.15, 0.20, 0.25, 0.30, 0.35, 0.4, 0.45, 0.50, 0.70, 0.90, 1.20, 1.50};
    const double erec_a[] = {0.0, 0.2, 0.4, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0};
    BinningX = MakeVector(q2_a);
    BinningY = MakeVector(erec_a);

    // q0 and q3 in true
  } else if (kinvars[Selection][0] == kq0 && kinvars[Selection][1] == kq3) {

    double q0_a[] = {0.0, 0.05, 0.1, 0.15, 0.20, 0.25, 0.30, 0.35, 0.4, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 1.00, 1.30};
    BinningX = MakeVector(q0_a);
    BinningY = MakeVector(q0_a);

    // Catch
  } else {
    std::cerr << "I have no binning specified for current plotting variables!" << std::endl;
    std::cerr << "Specified variables are " << ND280Kinematic_ToString(kinvars[Selection][0]) << " " << ND280Kinematic_ToString(kinvars[Selection][1]) << std::endl;
    throw;
  }
} // end SetupBinning

// ***************************************************************************
// Add a selection from psyche to samplePDFND2019 and turn it on in psyche
//
// Follows the enum outlined in psyche
// This is currently called by src/testFGD2.cpp and similars
// Could probably make it set-up internally with binning
void samplePDFND2019::addSelection( SampleId::SampleEnum selection, ND280KinematicTypes var1, ND280KinematicTypes var2) {
  // ***************************************************************************

  std::cout << std::setw(20) << "Adding selection: " << SampleId::ConvertSample(selection) << std::endl;

  // Enable the selection in psyche
  _man.sel().EnableSelection(SampleId::ConvertSampleToSelection(selection));

  std::string data = "_data";

  // Set the dimensions for this psyche selection
  ndims[selection]      = 2;

  // Set the kinematic variables for this psyche selection
  kinvars[selection][0] = var1;
  kinvars[selection][1] = var2;
  //kinvars[selection][0] = kq0;
  //kinvars[selection][1] = kq3;

  // Change here to plot in e.g. Q2QE and ErecQE
  // Really make this some setting instead in the future
  //kinvars[selection][0] = kQ2QE;
  //kinvars[selection][1] = kErecQE;

  //kinvars[selection][0] = kPionMomentum;
  //kinvars[selection][1] = kPionCosTheta;

  std::string sampleName = SampleId::ConvertSample(selection);
  std::replace(sampleName.begin(), sampleName.end(), ' ', '_');

  TH2Poly *MC_Sample = NULL;
  TH2Poly *Data_Sample = NULL;

  std::vector<double> xbins;
  std::vector<double> ybins;

  ///WPpoly
  MC_Sample = new TH2Poly();
  Data_Sample = new TH2Poly();

  SetupPoly(sampleName,MC_Sample);
  SetupPoly(sampleName,Data_Sample); //Gets poly binning from input file (defined in config)

  // Add the MC sample
  MC_Sample->GetXaxis()->SetTitle(ND280Kinematic_ToString(kinvars[selection][0]).c_str());
  MC_Sample->GetYaxis()->SetTitle(ND280Kinematic_ToString(kinvars[selection][1]).c_str());
  MC_Sample->GetZaxis()->SetTitle("Events");
  samplepdfs->AddAt(MC_Sample, selection);

  // Add the data sample
  Data_Sample->GetXaxis()->SetTitle(ND280Kinematic_ToString(kinvars[selection][0]).c_str());
  Data_Sample->GetYaxis()->SetTitle(ND280Kinematic_ToString(kinvars[selection][1]).c_str());
  Data_Sample->GetZaxis()->SetTitle("Events");
  datapdfs->AddAt(Data_Sample, selection);

} // end addSelection

// ***************************************************************************
// Function for reading in bins from 
void samplePDFND2019::SetupPoly(std::string sampleName, TH2Poly*& poly) {
// ***************************************************************************

  std::string polyBins =  FitManager->getPolyFile();
  TFile* fp = TFile::Open(polyBins.c_str());
  if(fp==NULL){
    std::cerr << "Input binning file not found. " << std::endl;
    throw;
  }

  poly = (TH2Poly*)fp->Get((sampleName).c_str())->Clone();
  poly->Reset("");

}

// ***************************************************************************
// Function which enables the samplePDFs to be plotted in accordance to interaction mode
void samplePDFND2019::enableModeHistograms() {
  // ***************************************************************************

  if (modepdf == true) {
    std::cout << "Already enabled modepdf but you're trying to do so again" << std::endl;
    std::cout << "Just returning to caller..." << std::endl;
    return;
  }

  modepdf = true;

  // Expand to make all samples
  samplemodepdfs->Expand(SampleId::kNSamples);
  samplemodepdfs->SetOwner(true);

  bool valid = false;

  for (int i = 0; i < SampleId::kNSamples; i++) {
    SampleId::SampleEnum sampleenum  = SampleId::SampleEnum(i);
    if (samplepdfs->At(i) == NULL) continue;
    // Make sure at least one sample is added (unfortunately the calling of enableModeHistograms() needs to happen _AFTER_ addSelection, or else we have no added selections at all
    valid = true;

    modeobjarray[i] = new TObjArray(0);
    modeobjarray[i]->Expand(kMaCh3_nModes+1);
    modeobjarray[i]->SetOwner(true);

    for (int j = 0; j < kMaCh3_nModes+1; j++) {
      TString name = SampleId::ConvertSample(sampleenum).c_str();
      name += j;
      if (ndims[i] == 1) {
        TH1D* tmp = ((TH1D*)samplepdfs->At(i));
        TH1D* clone = (TH1D*)tmp->Clone(name);
        modeobjarray[i]->AddAt(clone, j);
      } else if (ndims[i] == 2) {
        TH2Poly* tmp = ((TH2Poly*)samplepdfs->At(i));//WPpoly 
        TH2Poly* clone = (TH2Poly*)tmp->Clone(name);
        modeobjarray[i]->AddAt(clone, j);
        // Complain if ndims is some weird thing
      } else {
        std::cerr << "ndims[" << i << "] != 1 or 2 in enableModeHistograms()" << std::endl;
        std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
        throw;
      }
    } // End loop over interaction modes
    samplemodepdfs->AddAt(modeobjarray[i], i);
  } // End loop over psyche samples

  if (valid == false) {
    std::cerr << "Tried enabling mode pdfs without having any selections enabled" << std::endl;
    std::cerr << "Did you call enableModeHistograms() before you did addSelection()?" << std::endl;
    std::cerr << "Call addSelection() first, then enableModeHistograms()" << std::endl;
    throw;
  }
}

// ***************************************************************************
// Helper function to print MC rates for the psyche samples
void samplePDFND2019::printMCRates() {
  // ***************************************************************************
  std::cout << std::setw(40) << std::left << "Sample " << std::setw(10) << "MC" << std::endl;
  for (int i = 0; i < SampleId::kNSamples; i++) {
    if (samplepdfs->At(i) != NULL) {
      std::cout << std::setw(40) << std::left << SampleId::ConvertSample(SampleId::SampleEnum(i)) << std::setw(10) << NoOverflowIntegral(((TH2Poly*)(samplepdfs->At(i)))) << std::endl;
    }
  }
}

// ***************************************************************************
// Helper function to print Data rates for the psyche samples
void samplePDFND2019::printDataRates() {
  // ***************************************************************************
  std::cout << std::setw(40) << std::left << "Sample " << std::setw(10) << "Data" << std::endl;
  for (int i = 0; i < SampleId::kNSamples; i++) {
    if (datapdfs->At(i) != NULL) {
      std::cout << std::setw(40) << std::left << SampleId::ConvertSample(SampleId::SampleEnum(i)) << std::setw(10) << NoOverflowIntegral(((TH2Poly*)(datapdfs->At(i)))) << std::endl;
    }
  }
}

// ***************************************************************************
// Helper function to print rates for the psyche samples
void samplePDFND2019::printRates() {
  // ***************************************************************************
  std::cout << std::setw(40) << std::left << "Sample" << std::setw(10) << "Data" << std::setw(10) << "MC" << std::setw(10) << "-LLH" << std::setw(10) << "|" << std::endl;

  double sumData  = 0.0;
  double sumMC    = 0.0;
  for (int i = 0; i < SampleId::kNSamples; i++) {
    if (datapdfs->At(i) != NULL) {
      sumData += NoOverflowIntegral(((TH2Poly*)datapdfs->At(i)));//WPpoly
      sumMC   += NoOverflowIntegral(((TH2Poly*)samplepdfs->At(i)));
      double likelihood = getSampleLikelihood(i);
      std::cout << std::setw(40) << std::left << SampleId::ConvertSample(SampleId::SampleEnum(i)) << std::setw(10) << NoOverflowIntegral(((TH2Poly*)(datapdfs->At(i)))) << std::setw(10) << NoOverflowIntegral(((TH2Poly*)(samplepdfs->At(i)))) << std::setw(10) << likelihood << std::setw(10) << "|" << std::endl;//WPpoly
    }
  }
  double likelihood = getLikelihood();

  std::cout << std::setw(40) << std::left << "Total" << std::setw(10) << sumData << std::setw(10) << sumMC << std::setw(10) << likelihood << std::setw(10) << "|" << std::endl;
}


// ***************************************************************************
// Save the current PDFs (sample distributions) to an output file
void samplePDFND2019::savePDF(std::string inputLoc) {
  // ***************************************************************************

  // Strip out the .root part of saveLoc because we can
  std::string saveLoc = (inputLoc.substr(0, inputLoc.find(".root"))+"_ND280_nom.root");
  TFile *saveRates = new TFile(saveLoc.c_str(), "RECREATE");
  saveRates->cd();

  std::cout << "Saving data and MC distributions in " << saveLoc << std::endl;

  // Loop over the samples
  for (int i = 0; i < SampleId::kNSamples; i++) {
    if (datapdfs->At(i) == NULL) continue;
    // Force entries to be 1
    ((TH2Poly*)datapdfs->At(i))->SetEntries(1);

    // Replace the spaces in the sample name with underscores so we can make it the title of the histogram in the output file (wooew)
    std::string temp = SampleId::ConvertSample(SampleId::SampleEnum(i));
    while (temp.find(" ") != std::string::npos) {
      temp.replace(temp.find(" "), 1, std::string("_"));
    }
    //Get the 1d binning we want. Let's just use SetupBinning to get this as it already exists
    std::vector<double> xbins;
    std::vector<double> ybins;
    SetupBinning(SampleId::SampleEnum(i), xbins, ybins);

    // Write the 2D to file
    std::string TempName = std::string("DATA_")+temp;
    ((TH2Poly*)(datapdfs->At(i)))->Write(TempName.c_str());
    // And the projections onto pmu and cosmu
    // Project without including under/overflow
    //TH1D *ProjectX = ((TH2Poly*)(datapdfs->At(i)))->ProjectionX((TempName+"_x").c_str(), 1, ((TH2Poly*)(datapdfs->At(i)))->GetYaxis()->GetNbins(), "e");
    TH1D *ProjectX = PolyProjectionX(datapdfs->At(i), TempName, xbins);
    ProjectX->SetLineColor(kRed);
    ProjectX->Write((TempName+"_x").c_str());
    //TH1D *ProjectY = ((TH2Poly*)(datapdfs->At(i)))->ProjectionY((TempName+"_y").c_str(), 1, ((TH2Poly*)(datapdfs->At(i)))->GetXaxis()->GetNbins(), "e");
    TH1D *ProjectY = PolyProjectionY(datapdfs->At(i), TempName, ybins);
    ProjectY->SetLineColor(kRed);
    ProjectY->Write((TempName+"_y").c_str());
  }

  // Loop over the samples
  for (int i = 0; i < SampleId::kNSamples; i++) {
    if (samplepdfs->At(i) == NULL) continue;
    // Replace the spaces in the sample name with underscores so we can make it the title of the histogram in the output file (wooew)
    std::string temp = SampleId::ConvertSample(SampleId::SampleEnum(i));
    while (temp.find(" ") != std::string::npos) {
      temp.replace(temp.find(" "), 1, std::string("_"));
    }

    //Get the 1d binning we want. Let's just use SetupBinning to get this as it already exists
    std::vector<double> xbins;
    std::vector<double> ybins;
    SetupBinning(SampleId::SampleEnum(i), xbins, ybins);

    std::string TempName = std::string("MC_")+temp;

    ((TH2Poly*)samplepdfs->At(i))->SetEntries(1);
    // Write the 2D histogram
    ((TH2Poly*)(samplepdfs->At(i)))->Write(TempName.c_str());//WPpoly

    W2Hist[i]->SetEntries(1);
    W2Hist[i]->Write((TempName+"_w2").c_str());

    // Write the 1D histograms
    //TH1D *ProjectX = ((TH2Poly*)(samplepdfs->At(i)))->ProjectionX((TempName+"_x").c_str(), 1, ((TH2Poly*)(samplepdfs->At(i)))->GetYaxis()->GetNbins(), "e");//WPpoly
    TH1D *ProjectX = PolyProjectionX(samplepdfs->At(i), TempName, xbins);
    ProjectX->SetLineColor(kRed);
    ProjectX->Write((TempName+"_x").c_str());
    //TH1D *ProjectY = ((TH2Poly*)(samplepdfs->At(i)))->ProjectionY((TempName+"_y").c_str(), 1, ((TH2Poly*)(samplepdfs->At(i)))->GetXaxis()->GetNbins(), "e");//WPpoly
    TH1D *ProjectY = PolyProjectionY(samplepdfs->At(i), TempName, ybins);
    ProjectY->SetLineColor(kRed);
    ProjectY->Write((TempName+"_y").c_str());

    // Write the modes
    if (modepdf) {
      for (int j = 0; j < kMaCh3_nModes+1; j++) {
        std::string ModeName = std::string("MC_")+temp+"_"+MaCh3mode_ToString(MaCh3_Mode(j));
        TH2Poly* ModePlot = (TH2Poly*)(getPDFMode(i,j)->Clone());
        ModePlot->SetTitle(ModeName.c_str());
        ModePlot->SetLineColor(j);
        ModePlot->GetXaxis()->SetTitle("p_{#mu} (MeV)");
        ModePlot->GetYaxis()->SetTitle("cos #theta_{#mu}");
        ModePlot->GetZaxis()->SetTitle("Events");
        ModePlot->Write(ModeName.c_str());
        delete ModePlot;
      }
    }
    TDirectory *ebtempdir = saveRates->mkdir(("MC_"+temp+"_EbTemplates").c_str());
    ebtempdir->cd();
    for (int j = 0; j < 2; j++) {
      for (int k = 0; k < 4; k++) {
	std::string target;
	std::string species;
	if(j==0) target="C";
	if(j==1) target="O";
	if(k==0) species="numu";
	if(k==1) species="numubar";
	if(k==2) species="nue";
	if(k==3) species="nuebar";
	TDirectory *ebtempsubdir = ebtempdir->mkdir((species+"_"+target).c_str());
	ebtempsubdir->cd();
	for (int l = 0; l < 73; l++) {
	  ebtempsubdir->cd();
	  std::string bin = Form("_Bin_%d", l);
	  std::string TemplateName = std::string("MC_")+temp+"_CCQE_Target_"+target+"_"+species+bin.c_str();
	  TH2Poly* ebtemp = (TH2Poly*)(ebtemplates[i][j][k][l])->Clone();
	  //if(NoOverflowIntegral(ebtemp)>0){
	  //std::cout << temp << " " << target << " " << species << " " << bin.c_str() << " " << NoOverflowIntegral(ebtemp) << std::endl;
	  //}
	  ebtemp->SetTitle(TemplateName.c_str());
	  ebtemp->SetLineColor(j);
	  ebtemp->GetXaxis()->SetTitle("p_{#mu} (MeV)");
	  ebtemp->GetYaxis()->SetTitle("cos #theta_{#mu}");
	  ebtemp->GetZaxis()->SetTitle("Events");
	  ebtemp->Write(TemplateName.c_str());
	  delete ebtemp;
	  ebtempdir->cd();
	}
      }
    }
    saveRates->cd();
  }
  saveRates->Close();
  delete saveRates;
}



// ***************************************************************************
// Fil the data from the samples enabled
void samplePDFND2019::fillDataFromSamples() {
  // ***************************************************************************

  // Reset the datapdfs that we have enabled
  for (int i = 0; i < SampleId::kNSamples; i++) {
    if (datapdfs->At(i) == NULL) continue;
    else ((TH2Poly*)datapdfs->At(i))->Reset("");
  }

  // Make psyche quiet
  QuietPlease(); 

  // Create the Toy experiment with the appropriate format (number of systematics and number of parameters for each systematic)
  // This is essentially the same as the ToyMakerExample in psycheSteering/src but ensures there's no variation
  const unsigned int nExperiments = 1;
  ToyMaker.CreateToyExperiments(nExperiments, _man.syst().GetSystematics());
  toy = ToyMaker.GetToyExperiment(0);

  // Initialise the data TChain (should only do if we're not preloading)
  _man.InitialiseDataTChain();

  // Number of passed events
  int DataCnt = 0;
  int nNotEnabled = 0;

  // Number of data events in the TChain
  const long int nEntries  = _man.GetEntries();
  const unsigned int nData = _man.GetNDataEvents();

  // Create the toy box array (array of PreviousToyBox)
  _man.sel().CreateToyBoxArray(nData);
  // Initialise systematics
  _man.syst().Initialize(_man.sel(), nData);

  // These two are superfluous for data but are needed to call ProcessEvent
  Weight_h totalWeightSyst;
  Weight_h fluxWeightSyst;

  // Number of events
  UInt_t evt = 0;
  // Entry to get from the AnalysisManager (incremented by LoadSuperEvent)
  Long64_t entry = 0;
  // Start running over psyche events
  while (evt < nData && entry < nEntries) {

    // Get the next event in the Experiment
    AnaSuperEventB* sevent = NULL;
    // If we're preloading
    // WARNING, THIS MODIFIES ENTRY TO BE ENTRY += 1
    if (simple) {
      sevent = _man.LoadSuperEvent(entry);   
    } else {
      sevent = _man.GetPreloadedDataSuperEvent(entry);
    }

    // Get the event
    AnaEventB* event = dynamic_cast<AnaEventB*>(sevent->Event);

    // Fill the EventBox
    _man.sel().InitializeEvent(*event);
    _man.syst().InitializeEventSystematics(_man.sel(), *event);

    // Count how many pass selection
    bool passed = _man.ProcessEvent(*toy, *event, totalWeightSyst, fluxWeightSyst);

    // Check if passed
    if (passed) {

      // Get the summary
      AnaEventSummaryB* Summary = dynamic_cast<AnaEventSummaryB*>(event->Summary);
      // Get the Sample Enum
      SampleId::SampleEnum EventSample = Summary->EventSample;

      // Skip if we've turned off this selection
      if (datapdfs->At(EventSample) == NULL) {
        evt++;
        nNotEnabled++;
        continue;
      }

      // Sanity check the sample ID
      if (EventSample == SampleId::kUnassigned) {
        std::cerr << "EventSummary SampleID is unassigned for data evt " << evt << std::endl;
        std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
        throw;
      }

      // Get the x and y variables
      double varx = getKinematicVariable(kinvars[EventSample][0], event);
      double vary = getKinematicVariable(kinvars[EventSample][1], event);
      // 1D binning
      if (ndims[EventSample] == 1) {
        ((TH1*)datapdfs->At(EventSample))->Fill(varx);
        // 2D binning
      } else if (ndims[EventSample] == 2) {
        ((TH2Poly*)datapdfs->At(EventSample))->Fill(varx, vary);
      }
      // Increment this for events that pass the cuts
      DataCnt++;
    }
    // Finalise the event systematics for all event (even those that didn't pass cuts)
    _man.sel().FinalizeEvent(*event);
    _man.syst().FinalizeEventSystematics(*event);
    // Increment the event counter for the while loop
    evt++;
  } // Finish while loop over nevt

  // No longer silence psyche
  NowTalk();

  // Print data counter
  std::cout << DataCnt << " data events passed selection in psyche" << std::endl;
  std::cout << nNotEnabled << " data events didn't pass sample enabled test" << std::endl;
  std::cout << "Looped " << evt << " events and " << nEntries << " in tree" << std::endl;

  // Count the binned and unbinned data
  int BinnedDataCnt = 0;
  int UnbinnedDataCnt = 0;

  std::cout << std::setw(40) << "Sample" << std::setw(10) << "Binned" << std::setw(10) << "Unbinned" << std::endl;
  for (int i = 0; i < SampleId::kNSamples; ++i) {
    if (datapdfs->At(i) == NULL) continue;
    // Binned 
    double Integral = NoOverflowIntegral(((TH2Poly*)datapdfs->At(i)));
    // Include under/overflow
    double IntegralUn = OverflowIntegral((TH2Poly*)datapdfs->At(i));
    
    std::cout << std::setw(40) << SampleId::ConvertSample(SampleId::SampleEnum(i)) << std::setw(10) << Integral << std::setw(10) << IntegralUn << std::endl;
    BinnedDataCnt += Integral;
    UnbinnedDataCnt += IntegralUn;
  }
  std::cout << std::setw(40) << "Total: " << std::setw(10) << BinnedDataCnt << std::setw(10) << UnbinnedDataCnt<< std::endl;

}

// ***************************************************************************
// Change the starting position of the chain by throwing the systematics
void samplePDFND2019::RandomStart() {
  // ***************************************************************************

  HaveIRandomStart = true;

  CheckCovariances();

  std::vector<covarianceBase*> Systematics;
  if (flux) {
    Systematics.push_back(flux);
  }
  if (xs2015) {
    Systematics.push_back(xs2015);
  }
  if (detsimple) {
    Systematics.push_back(detsimple);
  }

  // Reconfigure the systematics randomly
  for (std::vector<covarianceBase*>::iterator it = Systematics.begin(); it != Systematics.end(); it++) {
    (*it)->RandomConfiguration();
  }

  // Then reweight the simulation to the random configuration
  std::cout << "Reweighting the sample for a random start..." << std::endl;
  double *fake = NULL;
  reweight(fake);

  // Print the rates
  std::cout << "Printing rates after the random start:" << std::endl;
  printRates();
}

// ***************************************************************************
// Reweight the Monte-Carlo at ND280 for one step of the Markov Chain
void samplePDFND2019::reweight(double* oscpar) {
  // ***************************************************************************

  // Reset the histograms before reweight
  ResetHistograms();
  // If we want to run in ND280 "simple" mode; i.e. throw the flux and detector parameters from the covariances
  // We have to set the simple variable
  if (simple) {
    ReWeight_NoPsyche();
    return;
  }
  if (xsecmodel >= kNIWG_2020a) {
    std::cerr << "Can't run non-simple reweight with NIWG 2016a or later" << std::endl;
    std::cerr << "This is simply a lack of implementation rather than physical impossibility!" << std::endl;
    std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
    throw;
  }

  // Here is the "not-so-simple" ND280 reconfiguration, which interfaces to psyche to get flux and detector weights on an event by event basis
  // Has not been updated for psyche v3
  // Check e.g. samplePDFND2014 for how it was done before
  // Removed for code bloat
}

// ***************************************************************************
void samplePDFND2019::ReWeight_NoPsyche() {
  // ***************************************************************************

  if(detsimple) {
    detsimple->transformSpace();
  }

  // Prepare the weights
  // e.g. find the relevant spline segments, pass stuff to GPU
  PrepareWeights();

  // Call entirely different routine if we're running with openMP
#ifdef MULTITHREAD
  ReWeight_NoPsyche_MP();
#else

  if (xsecmodel >= kNIWG_2020a) {

    int nEvents = ndobj->nEvents;
    // Loop over all the events
    // Fill their weights for this reconfiguration
    for (int i = 0; i < nEvents; i++) {

      // If sample isn't enabled
      if (samplepdfs->At(ndobj->event_sample[i]) == NULL) continue;

      double beamw  = 1.;
      double detw   = 1.;
      double xsecw  = 1.;

      // Get the flux weight
      if(flux) {
        beamw = flux->calcReWeight(ndobj->flux_e_bin[i]);
      }

      // Get the detector weight
      if(detsimple) {
        detw = detsimple->calcReWeight(ndobj->det_bin[i]);
      }

      // Get the cross-section weights
      if (xs2015) {
        xsecw = CalcXsecWeight(i);
      }

      // Fill the MC
      double weight = ndobj->flux_w[i]*beamw*detw*xsecw;
      ((TH1*)samplepdfs->At(ndobj->event_sample[i]))->AddBinContent(ndobj->hist_bin[i], weight);

      // This one is needed for Barlow-Beeston calculation
      // Could use TH2D::GetSumw2() instead but that wouldn't work for the multi-threaded case, so use this one for validation, and check against samplepdfs
      W2Hist[ndobj->event_sample[i]]->AddBinContent(ndobj->hist_bin[i], weight*weight);

      if (modepdf) {
        ((TH1*)((TObjArray*)samplemodepdfs->At(ndobj->event_sample[i]))->At(xsecInfo[i].mode))->AddBinContent(ndobj->hist_bin[i],ndobj->flux_w[i]*beamw*detw*xsecw);
      }

    } // end for loop

  }
#endif // end the else in openMP
}

#ifdef MULTITHREAD
// ***************************************************************************
// routine for ReWeight_NoPsyche on openMP
void samplePDFND2019::ReWeight_NoPsyche_MP() {
  // ***************************************************************************

  double ebParamValCnu = xs2015->calcReWeight(Eb_C_nu_pos);
  double ebParamValCnubar = xs2015->calcReWeight(Eb_C_nubar_pos);
  double ebParamValOnu = xs2015->calcReWeight(Eb_O_nu_pos);
  double ebParamValOnubar = xs2015->calcReWeight(Eb_O_nubar_pos);

  //WP:Set EB dials
  varier_eb->SetSystematic(niwg::rew::kEbTwkDial_Vary_C12_numu, ebParamValCnu);
  varier_eb->SetSystematic(niwg::rew::kEbTwkDial_Vary_C12_numub, ebParamValCnubar);
  varier_eb->SetSystematic(niwg::rew::kEbTwkDial_Vary_C12_nue, ebParamValCnu);
  varier_eb->SetSystematic(niwg::rew::kEbTwkDial_Vary_C12_nueb, ebParamValCnubar);
  varier_eb->SetSystematic(niwg::rew::kEbTwkDial_Vary_O16_numu, ebParamValOnu);
  varier_eb->SetSystematic(niwg::rew::kEbTwkDial_Vary_O16_numub, ebParamValOnubar);
  varier_eb->SetSystematic(niwg::rew::kEbTwkDial_Vary_O16_nue, ebParamValOnu);
  varier_eb->SetSystematic(niwg::rew::kEbTwkDial_Vary_O16_nueb, ebParamValOnubar);
  varier_eb->Reconfigure();

  int nEvents = ndobj->nEvents;
  // First find how many samples samplePDF contains; set in samplePDF::reweight
  // Each event has: a hist_bin[i], mom, theta, event_sample, flux_e_bin, det_bin, xsecInfo or xsecInfo, flux_w

  // The master array
  double **samplePDF_array = new double*[nSamples]();
  // The master array of weights^2, needed for Barlow Beeston
  double **samplePDF_w2_array = new double*[nSamples]();
  for (__int__ i = 0; i < nSamples; i++) {
    samplePDF_array[i] = new double[maxBins]();
    samplePDF_w2_array[i] = new double[maxBins]();
    for (__int__ j = 0; j < maxBins; ++j) {
      samplePDF_array[i][j] = 0.0;
      samplePDF_w2_array[i][j] = 0.0;
    }
  }

  // The master mode array
  // Has an extra index to denote the mode
  double ***samplePDF_mode_array = NULL;

  // Allocate only if modepdf
  if (modepdf) {
    samplePDF_mode_array = new double**[nSamples]();
    for (__int__ i = 0; i < nSamples; i++) {
      samplePDF_mode_array[i] = new double*[kMaCh3_nModes+1]();
      for (__int__ j = 0; j < kMaCh3_nModes+1; j++) {
        samplePDF_mode_array[i][j] = new double[maxBins]();
        for (__int__ k = 0; k < maxBins; ++k) {
          samplePDF_mode_array[i][j][k] = 0.0;
        }
      }
    }
  }

  //For making Eb templates
  double *****ebtemplates_array = NULL;
  ebtemplates_array = new double****[nSamples]();
    for (__int__ i = 0; i < nSamples; i++) {
      ebtemplates_array[i] = new double***[2]();
      for (__int__ j = 0; j < 2; j++) {
        ebtemplates_array[i][j] = new double**[4]();
        for (__int__ k = 0; k < 4; ++k) {
          ebtemplates_array[i][j][k] = new double*[73]();
	  for(__int__ l = 0; l < 73; ++l) {
	    ebtemplates_array[i][j][k][l] = new double[maxBins]();
	    for(__int__ m=0; m < maxBins; ++m) {
	      ebtemplates_array[i][j][k][l][m] = 0.0;
	    }
	  }
	}
      }
    }
  // The per-thread array
  double **samplePDF_array_private = NULL;
  // And the mode
  double ***samplePDF_mode_array_private = NULL;
  //eb template
  double *****ebtemplates_private_array = NULL;
  // The per-thread array of weight^2
  double **samplePDF_w2_array_private = NULL;

  // Declare the omp parallel region
  // The parallel region needs to stretch beyond the for loop!
#pragma omp parallel private(samplePDF_array_private, samplePDF_mode_array_private, samplePDF_w2_array_private)
  {

    // private to each thread
    samplePDF_array_private = new double*[nSamples]();
    samplePDF_w2_array_private = new double*[nSamples]();
    for (__int__ i = 0; i < nSamples; i++) {
      samplePDF_array_private[i] = new double[maxBins]();
      samplePDF_w2_array_private[i] = new double[maxBins]();
      for(__int__ j = 0; j < maxBins; ++j) {
        samplePDF_array_private[i][j] = 0.0;
        samplePDF_w2_array_private[i][j] = 0.0;
      }
    }

    // Do the same for the mode arrays
    if (modepdf) {
      samplePDF_mode_array_private = new double**[nSamples]();
      for (__int__ i = 0; i < nSamples; i++) {
        samplePDF_mode_array_private[i] = new double*[kMaCh3_nModes+1]();
        for (__int__ j = 0; j < kMaCh3_nModes+1; j++) {
          samplePDF_mode_array_private[i][j] = new double[maxBins]();
          for (__int__ k = 0; k < maxBins; ++k) {
            samplePDF_mode_array_private[i][j][k] = 0.0;
          }
        }
      }

      // Set the samplePDF_mode_array_private to zero, even if we're aren't running modepdf
      // This is purely to avoid compiler warning
    } else {
      samplePDF_mode_array_private = new double**[1]();
      samplePDF_mode_array_private[0] = new double*[1]();
      samplePDF_mode_array_private[0][0] = new double[1]();
      samplePDF_mode_array_private[0][0][0] = 0.0;
    }

  //For making Eb templates
  ebtemplates_private_array = new double****[nSamples]();
    for (__int__ i = 0; i < nSamples; i++) {
      ebtemplates_private_array[i] = new double***[2]();
      for (__int__ j = 0; j < 2; j++) {
        ebtemplates_private_array[i][j] = new double**[4]();
        for (__int__ k = 0; k < 4; ++k) {
          ebtemplates_private_array[i][j][k] = new double*[73]();
          for(__int__ l = 0; l < 73; ++l) {
            ebtemplates_private_array[i][j][k][l] = new double[maxBins]();
	    for(__int__ m = 0; m < maxBins; ++m) {
              ebtemplates_private_array[i][j][k][l][m] = 0.0;
            }
          }
      	}
      }
    }

    // The master for loop: this is where we gain all performace by multi-threading
#pragma omp for
    for (int i = 0; i < nEvents; i++) {

      if (samplepdfs->At(ndobj->event_sample[i]) == NULL) continue;

      // Only do Eb shift for 12C and 16O, when mode is CCQE and Enu is < 4 GeV
      if((xsecInfo[i].target == kTarget_O || xsecInfo[i].target == kTarget_C) && (xsecInfo[i].mode == kMaCh3_CCQE) && (xsecInfo[i].Enu < 4) && ndobj->hist_bin[i] > 0)
      {
	//WP:Applying EB dial
	ndobj->mom_adj[i] = ndobj->mom[i] + varier_eb->GetMomentumVariation(ndobj->eb_dial[i],ndobj->eb_bin[i]);
	
	//Need to check what bin we're now in after Eb shift. Could just do FindBin but that's sloooow, especially for Polys. The shifts are of order a few MeV, and bin widths are ~50MeV so vast majority of shifts should not change the bin or change by +/- 1 in x direction. So we check those first then. Bear in mind this only works if the adjacent bins are numbered consecutively. i.e the bin immediately to the left is one global index lower, and the one to the right is one higher. For TPolys, this is dependent on the order you 'AddBin'-ed when you set up the initial poly. Also this will work best if bins are constant in theta and all rectangular

	if(ndobj->mom_adj[i]<0){
	  ndobj->mom_adj[i]=0.0001;
	}

	if(polybins[ndobj->event_sample[i]][ndobj->hist_bin[i]-1]->IsInside(ndobj->mom_adj[i], ndobj->theta[i]) )
	  ;
	else if(ndobj->hist_bin[i] < ((TH2Poly*)(samplepdfs->At(ndobj->event_sample[i])))->GetNumberOfBins() && polybins[ndobj->event_sample[i]][ndobj->hist_bin[i]]->IsInside(ndobj->mom_adj[i], ndobj->theta[i]) )
	  {
	    ndobj->hist_bin[i]++;
	  }
	else if (ndobj->hist_bin[i] > 1 && polybins[ndobj->event_sample[i]][ndobj->hist_bin[i]-2]->IsInside(ndobj->mom_adj[i], ndobj->theta[i]))
	  {
	    ndobj->hist_bin[i]--;
	  }
	else
	  {
	  ndobj->hist_bin[i] = ((TH2Poly*)samplepdfs->At(ndobj->event_sample[i]))->FindBin(ndobj->mom_adj[i], ndobj->theta[i]);
	  }
	//WPPoly Fix->NotFix
      }
      
      double beamw   = 1.;
      double detw    = 1.;
      double xsecw   = 1.;

      // Get the flux weight
      if (flux) {
        beamw = flux->calcReWeight(ndobj->flux_e_bin[i]);
      }
      // Get the detector weight
      if (detsimple) {
	//a bit hacky: if the event is in the overflow, we assign the detector bin as the last one in that row. BANFF are just not applying a weight though, but it would be messy to not assign a bin at all. So now we're still assigning that last bin in the row, but now applying the weight if we're above 30000. This only affects the overflow, and is a very small effect anyway. See covariance/covarianceNDDet_2014Simple:getBin() for how the bin is assignment is done
	if(ndobj->mom[i]<30000 && ndobj->mom[i]>0){
	  detw = detsimple->calcReWeight(ndobj->det_bin[i]);
	}
      }

      // Get the cross-section weights
      if (xs2015) {
        xsecw = CalcXsecWeight(i);
      }

      // Is simply the total weight to apply to one event
      double WeightTot = ndobj->flux_w[i]*beamw*detw*xsecw;

      // Fill each threads' private array up
      if(ndobj->hist_bin[i]>=0){
	samplePDF_array_private[ndobj->event_sample[i]][ndobj->hist_bin[i]] += WeightTot;
        // And the w2
        samplePDF_w2_array_private[ndobj->event_sample[i]][ndobj->hist_bin[i]] += WeightTot*WeightTot;
      } else {
	samplePDF_array_private[ndobj->event_sample[i]][maxBins+ndobj->hist_bin[i]] += WeightTot;
        // And the w2
        samplePDF_w2_array_private[ndobj->event_sample[i]][maxBins+ndobj->hist_bin[i]] += WeightTot*WeightTot;
      }

      int targ=0, spec=0, erecbin=0;
      if(xsecInfo[i].target==12) targ=0;
      if(xsecInfo[i].target==16) targ=1;
      if(xsecInfo[i].species==12) spec=0;
      if(xsecInfo[i].species==-12) spec=1;
      if(xsecInfo[i].species==14) spec=2;
      if(xsecInfo[i].species==-14) spec=3;
      if(spec < 2) erecbin = NumuErec->FindBin(xsecInfo[i].Enu);
      else erecbin = NueErec->FindBin(xsecInfo[i].Enu);

      if((xsecInfo[i].target == kTarget_O || xsecInfo[i].target == kTarget_C) && (xsecInfo[i].mode == kMaCh3_CCQE) && (xsecInfo[i].Enu < 4) && ndobj->hist_bin[i] > 0){
	ebtemplates_private_array[ndobj->event_sample[i]][targ][spec][erecbin][ndobj->hist_bin[i]]+=WeightTot;
      }
      
      // Fill each threads' private array up
      if (modepdf) {
	if(ndobj->hist_bin[i]>=0){
	  samplePDF_mode_array_private[ndobj->event_sample[i]][xsecInfo[i].mode][ndobj->hist_bin[i]] += WeightTot;
        } else {
          samplePDF_mode_array_private[ndobj->event_sample[i]][xsecInfo[i].mode][maxBins+ndobj->hist_bin[i]] += WeightTot;
	}
      }

#if DEBUG > 0
      // Get limits of plot
      //This is a bit mad, SetBinContent does different things for the overflow than to the rest of the bins for a TH2Poly
      //The gift that keeps on giving...
      if(ndobj->hist_bin[i]<0){
        XsecOnly[ndobj->event_sample[i]]->SetBinContent(ndobj->hist_bin[i], POTWeights[i]*xsecw);
        NDCovOnly[ndobj->event_sample[i]]->SetBinContent(ndobj->hist_bin[i], POTWeights[i]*detw);
        BeamCovOnly[ndobj->event_sample[i]]->SetBinContent(ndobj->hist_bin[i], POTWeights[i]*beamw);
        AllOnly[ndobj->event_sample[i]]->SetBinContent(ndobj->hist_bin[i], ndobj->flux_w[i]*beamw*detw*xsecw);
      } else {
	XsecOnly[ndobj->event_sample[i]]->SetBinContent(ndobj->hist_bin[i], POTWeights[i]*xsecw+XsecOnly[ndobj->event_sample[i]]->GetBinContent(ndobj->hist_bin[i]));
	NDCovOnly[ndobj->event_sample[i]]->SetBinContent(ndobj->hist_bin[i], POTWeights[i]*detw+NDCovOnly[ndobj->event_sample[i]]->GetBinContent(ndobj->hist_bin[i]));
	BeamCovOnly[ndobj->event_sample[i]]->SetBinContent(ndobj->hist_bin[i], POTWeights[i]*beamw+BeamCovOnly[ndobj->event_sample[i]]->GetBinContent(ndobj->hist_bin[i]));
	AllOnly[ndobj->event_sample[i]]->SetBinContent(ndobj->hist_bin[i], ndobj->flux_w[i]*beamw*detw*xsecw+AllOnly[ndobj->event_sample[i]]->GetBinContent(ndobj->hist_bin[i]));
      }
      // Check for negative or large weights!
      bool foundweird = false;
      if (ndobj->flux_w[i] < 0 || ndobj->flux_w[i] > __LARGE_WEIGHT__) {
        std::cerr << "NEGATIVE WEIGHT flux weight = " << ndobj->flux_w[i] << std::endl;
        std::cerr << "Event # " << i << std::endl;
        foundweird = true;
      }
      if (beamw < 0 || beamw > __LARGE_WEIGHT__) {
        std::cerr << "NEGATIVE WEIGHT beam weight = " << beamw << std::endl;
        std::cerr << "Event # " << i << std::endl;
        foundweird = true;
      }
      if (xsecw < 0 || xsecw > __LARGE_WEIGHT__) {
        std::cerr << "NEGATIVE WEIGHT xsecw = " << xsecw << std::endl;
        std::cerr << "Event # " << i << std::endl;
        foundweird = true;
      }
      if (detw < 0 || detw > __LARGE_WEIGHT__) {
        std::cerr << "NEGATIVE WEIGHT det weight = " << detw << std::endl;
        std::cerr << "Event # " << i << std::endl;
        foundweird = true;
      }
      if (std::isnan(WeightTot)) {
        std::cerr << "Found a nan total weight!" << std::endl;
        foundweird = true;
      }
#endif

      // Print the structs if we're interested in event by event breakdowns
#if DEBUG == 2
      PrintStructs(xsecw, beamw, detw, i);
#endif

      // If we found something weird print everything we've got!
#if DEBUG > 0
      if (foundweird) {
        DumpSplines(xsecw, beamw, detw, i);
        throw;
      }
#endif

    } // end for loop and openMP for, but still keep the omp region!
    // Now we can write the individual arrays from each thread to the main array
    // This can probably be done by directly using samplepdfs->At(i)...
#pragma omp parallel shared(samplePDF_array, samplePDF_w2_array)
    {
      for (__int__ i = 0; i < nSamples; i++) {
        for (__int__ j = 0; j < maxBins; j++) {
#pragma omp atomic
          samplePDF_array[i][j] += samplePDF_array_private[i][j];
          samplePDF_w2_array[i][j] += samplePDF_w2_array_private[i][j];
        }
      }
    } // end the shared region for samplePDF_array

    // Do the same for the mode histograms
#pragma omp parallel shared(samplePDF_mode_array)
    {
      if (modepdf) {
        for (__int__ i = 0; i < nSamples; i++) {
          for (__int__ j = 0; j < kMaCh3_nModes+1; j++) {
            for (__int__ k = 0; k < maxBins; k++) {
              // Do the same summation for the mode arrays
#pragma omp atomic
              samplePDF_mode_array[i][j][k] += samplePDF_mode_array_private[i][j][k];
            }
          }
        }
      } // end the shared region for samplePDF_mode_array
    }

        // Do the same for the ebtemplates
#pragma omp parallel shared(ebtemplates_array)
    {
        for (__int__ i = 0; i < nSamples; i++) {
          for (__int__ j = 0; j < 2; j++) {
            for (__int__ k = 0; k < 4; k++) {
	      for (__int__ l = 0; l < 73; l++) {
		for (__int__ m = 0; m < maxBins; m++) {
		  // Do the same summation for the eb templates
#pragma omp atomic
		  ebtemplates_array[i][j][k][l][m] += ebtemplates_private_array[i][j][k][l][m];
		}
	      }
	    }
	  }
	} // end the shared region for eb templates
    }
    

    
    // Delete each thread's private array (still in OMP block)
    for (__int__ i = 0; i < nSamples; i++) {
      delete[] samplePDF_array_private[i];
      delete[] samplePDF_w2_array_private[i];
    }
    delete[] samplePDF_array_private;
    delete[] samplePDF_w2_array_private;

    // Delete each thread's private mode array (we're still in a OMP block)
    // This have been assigned even if we aren't running with modepdf
    if (modepdf) {
      for (__int__ i = 0; i < nSamples; i++) {
        for (__int__ j = 0; j < kMaCh3_nModes+1; j++) {
          delete[] samplePDF_mode_array_private[i][j];
        }
        delete[] samplePDF_mode_array_private[i];
      }
      delete[] samplePDF_mode_array_private;
    } else {
      delete[] samplePDF_mode_array_private[0][0];
      delete[] samplePDF_mode_array_private[0];
      delete[] samplePDF_mode_array_private;
    }

    for (__int__ i = 0; i < nSamples; i++) {
        for (__int__ j = 0; j < 2; j++) {
	  for (__int__ k = 0; k < 4; k++) {
	    for (__int__ l =0; l < 73; l++) {
	      delete[] ebtemplates_private_array[i][j][k][l];
	    }
	    delete[] ebtemplates_private_array[i][j][k];
	  }
	  delete[] ebtemplates_private_array[i][j];
	}
	delete[] ebtemplates_private_array[i];
    }
    delete[] ebtemplates_private_array;
      
    
  } // end the #pragma omp parallel region

  // Write the collected array samplePDF_array to the samplepdf object (TH2D or TH1D)
  for (__int__ i = 0; i < nSamples; i++) {
    if (samplepdfs->At(i) == NULL) continue;
    // Check if this is the first time we fill this sample (only want to fill once)
    bool firsttime = (NoOverflowIntegral(W2Hist[i]) == 0);
    for (__int__ j = 0; j < maxBins-9; j++) {
      ((TH2Poly*)samplepdfs->At(i))->SetBinContent(j, samplePDF_array[i][j]+((TH2Poly*)samplepdfs->At(i))->GetBinContent(j));
      // Fill only the first time, Hack!
      if (firsttime) W2Hist[i]->SetBinContent(j, samplePDF_w2_array[i][j]+W2Hist[i]->GetBinContent(j));
    }
    for (__int__ j =-9; j<0; j++) {
      ((TH2Poly*)samplepdfs->At(i))->SetBinContent(j, samplePDF_array[i][maxBins+j]+((TH2Poly*)samplepdfs->At(i))->GetBinContent(j));
      // Fill only the first time, Hack!
      if (firsttime) W2Hist[i]->SetBinContent(j, samplePDF_w2_array[i][maxBins+j]+W2Hist[i]->GetBinContent(j));
    }
  }

  if (modepdf) {
    for (__int__ i = 0; i < nSamples; i++) {
      if (samplepdfs->At(i) == NULL) continue;
      for (__int__ j = 0; j < kMaCh3_nModes+1; j++) {
        for (__int__ k = 0; k < maxBins-9; k++) {
          ((TH2Poly*)((TObjArray*)(samplemodepdfs->At(i)))->At(j))->SetBinContent(k, samplePDF_mode_array[i][j][k]+((TH2Poly*)((TObjArray*)(samplemodepdfs->At(i)))->At(j))->GetBinContent(k));
        }
        for (__int__ k =-9; k<0; k++) {
	  ((TH2Poly*)((TObjArray*)(samplemodepdfs->At(i)))->At(j))->SetBinContent(k, samplePDF_mode_array[i][j][maxBins+k]+((TH2Poly*)((TObjArray*)(samplemodepdfs->At(i)))->At(j))->GetBinContent(k));
	}
      }
    }
  }


for (__int__ i = 0; i < nSamples; i++) {
      if (samplepdfs->At(i) == NULL) continue;
      for (__int__ j = 0; j < 2; j++) {
        for (__int__ k = 0; k < 4; k++) {
	  for(__int__ l = 0; l < 73; l++) {
	    for(__int__ m = 0; m < maxBins; m++) {
	      if(ebtemplates_array[i][j][k][l][m] > 0){
		ebtemplates[i][j][k][l]->SetBinContent(m,ebtemplates_array[i][j][k][l][m] + ebtemplates[i][j][k][l]->GetBinContent(m));
	      }
	    }
	  }
	}
      }
 }


  // Delete the master array
  for (__int__ i = 0; i < nSamples; i++) {
    delete[] samplePDF_array[i];
    delete[] samplePDF_w2_array[i];
  }
  delete[] samplePDF_array;
  delete[] samplePDF_w2_array;

  // Delete the master mode array
  if (modepdf) {
    for (__int__ i = 0; i < nSamples; i++) {
      for (__int__ j = 0; j < kMaCh3_nModes+1; j++) {
        delete[] samplePDF_mode_array[i][j];
      }
      delete[] samplePDF_mode_array[i];
    }
    delete[] samplePDF_mode_array;
  }

  for (__int__ i = 0; i < nSamples; i++) {
        for (__int__ j = 0; j < 2; j++) {
          for (__int__ k = 0; k < 4; k++) {
            for (__int__ l =0; l<73; l++) {
              delete[] ebtemplates_array[i][j][k][l];
            }
            delete[] ebtemplates_array[i][j][k];
          }
          delete[] ebtemplates_array[i][j];
        }
        delete[] ebtemplates_array[i];
    }
    delete[] ebtemplates_array;

  
#if DEBUG > 0
  std::cout << "**************PRINTING DETAILED SAMPLE INFORMATION AFTER LOOP************" << std::endl;
  std::string FileName;
  if (FitManager != NULL) {
    FileName = FitManager->getOutputFilename();
    // Strip out the ending .root extension
    FileName = FileName.substr(0, FileName.find(".root"))+"_diagnostics.root";
  } else {
    FileName = "Diagnostics.root";
  }
  // Get the TFile
  TFile *DiagFile = new TFile(FileName.c_str(), "RECREATE");
  double MCOnlyBin = 0.0;
  double MCOnlyUn = 0.0;
  double POTOnlyBin = 0.0;
  double POTOnlyUn = 0.0;
  double FluxOnlyBin = 0.0;
  double FluxOnlyUn = 0.0;
  double DetOnlyBin = 0.0;
  double DetOnlyUn = 0.0;
  double XsecOnlyBin = 0.0;
  double XsecOnlyUn = 0.0;
  double NDCovOnlyBin = 0.0;
  double NDCovOnlyUn = 0.0;
  double BeamCovOnlyBin = 0.0;
  double BeamCovOnlyUn = 0.0;
  double AllOnlyBin = 0.0;
  double AllOnlyUn = 0.0;
  // Loop over the samples
  for (int i = 0; i < SampleId::kNSamples; ++i) {
    if (samplepdfs->At(i) == NULL || NoOverflowIntegral(((TH2Poly*)(samplepdfs->At(i)))) == 0) continue;//WPpoly
    std::cout << SampleId::ConvertSample(SampleId::SampleEnum(i)) << ":" << std::endl;
    std::cout << "  MC only:             " << NoOverflowIntegral(MCOnly[i]) << std::endl;
    MCOnlyBin += NoOverflowIntegral(MCOnly[i]) ;
    std::cout << "  MC only nobin:       " << OverflowIntegral(MCOnly[i]) << std::endl;
    MCOnlyUn += OverflowIntegral(MCOnly[i]) ;

    std::cout << "  POT only:            " << NoOverflowIntegral(POTOnly[i]) << std::endl;
    POTOnlyBin += NoOverflowIntegral(POTOnly[i]) ;
    std::cout << "  POT only nobin:      " << OverflowIntegral(POTOnly[i]) << std::endl;
    POTOnlyUn += OverflowIntegral(POTOnly[i]) ;

    std::cout << "  Flux only:           " << NoOverflowIntegral(FluxOnly[i]) << std::endl;
    FluxOnlyBin += NoOverflowIntegral(FluxOnly[i]) ;
    std::cout << "  Flux only nobin:     " << OverflowIntegral(FluxOnly[i]) << std::endl;
    FluxOnlyUn += OverflowIntegral(FluxOnly[i]) ;

    std::cout << "  Xsec only:           " << NoOverflowIntegral(XsecOnly[i]) << std::endl;
    XsecOnlyBin += NoOverflowIntegral(XsecOnly[i]) ;
    std::cout << "  Xsec only nobin:     " << OverflowIntegral(XsecOnly[i]) << std::endl;
    XsecOnlyUn += OverflowIntegral(XsecOnly[i]) ;

    std::cout << "  Det only:            " << NoOverflowIntegral(DetOnly[i]) << std::endl;
    DetOnlyBin += NoOverflowIntegral(DetOnly[i]) ;
    std::cout << "  Det only nobin:      " << OverflowIntegral(DetOnly[i]) << std::endl;
    DetOnlyUn += OverflowIntegral(DetOnly[i]) ;

    std::cout << "  ND Cov only:         " << NoOverflowIntegral(NDCovOnly[i]) << std::endl;
    NDCovOnlyBin += NoOverflowIntegral(NDCovOnly[i]) ;
    std::cout << "  ND Cov nobin:        " << OverflowIntegral(NDCovOnly[i]) << std::endl;
    NDCovOnlyUn += OverflowIntegral(NDCovOnly[i]) ;

    std::cout << "  Beam Cov only:       " << NoOverflowIntegral(BeamCovOnly[i]) << std::endl;
    BeamCovOnlyBin += NoOverflowIntegral(BeamCovOnly[i]) ;
    std::cout << "  Beam Cov nobin:      " << OverflowIntegral(BeamCovOnly[i]) << std::endl;
    BeamCovOnlyUn += OverflowIntegral(BeamCovOnly[i]) ;

    std::cout << "  All comb:            " << NoOverflowIntegral(AllOnly[i]) << std::endl;
    AllOnlyBin += NoOverflowIntegral(AllOnly[i]) ;
    std::cout << "  All comb nobin:      " << OverflowIntegral(AllOnly[i]) << std::endl;
    AllOnlyUn += OverflowIntegral(AllOnly[i]) ;

    // Write
    DiagFile->cd();

    // Do this or else root doesn't draw these (it thinks entries == 0 so doesn't draw anything). This is because we're using AddBinContent rather than Fill
    MCOnly[i]->SetEntries(1);
    POTOnly[i]->SetEntries(1);
    FluxOnly[i]->SetEntries(1);
    XsecOnly[i]->SetEntries(1);
    DetOnly[i]->SetEntries(1);
    NDCovOnly[i]->SetEntries(1);
    BeamCovOnly[i]->SetEntries(1);
    AllOnly[i]->SetEntries(1);
    NuMuPmu[i]->SetEntries(1);

    MCOnly[i]->Write();
    POTOnly[i]->Write();
    FluxOnly[i]->Write();
    XsecOnly[i]->Write();
    DetOnly[i]->Write();
    NDCovOnly[i]->Write();
    BeamCovOnly[i]->Write();
    AllOnly[i]->Write();
    NuMuPmu[i]->Write();
  }
  std::cout << "Totals: " << std::endl;
  std::cout << "  MCOnly bin:      " << MCOnlyBin << std::endl;
  std::cout << "  MCOnly un:       " << MCOnlyUn << std::endl;
  std::cout << "  POTOnly bin:     " << POTOnlyBin << std::endl;
  std::cout << "  POTOnly un:      " << POTOnlyUn << std::endl;
  std::cout << "  FluxOnly bin:    " << FluxOnlyBin << std::endl;
  std::cout << "  FluxOnly un:     " << FluxOnlyUn << std::endl;
  std::cout << "  XsecOnly bin:    " << XsecOnlyBin << std::endl;
  std::cout << "  XsecOnly un:     " << XsecOnlyUn << std::endl;
  std::cout << "  DetOnly bin:     " << DetOnlyBin << std::endl;
  std::cout << "  DetOnly un:      " << DetOnlyUn << std::endl;
  std::cout << "  NDCovOnly bin:   " << NDCovOnlyBin << std::endl;
  std::cout << "  NDCovOnly un:    " << NDCovOnlyUn << std::endl;
  std::cout << "  BeamCovOnly bin: " << BeamCovOnlyBin << std::endl;
  std::cout << "  BeamCovOnly un:  " << BeamCovOnlyUn << std::endl;
  std::cout << "  AllOnly bin:     " << AllOnlyBin << std::endl;
  std::cout << "  AllOnly un:      " << AllOnlyUn << std::endl;

  DiagFile->Write();
  DiagFile->Close();
  delete DiagFile;
  std::cout << "With xsec systematics: " << std::endl;
  xs2015->printNominalCurrProp();
#endif

} // end function
#endif


// *******************************************
// Calculate the cross-section weight
double samplePDFND2019::CalcXsecWeight(const int i) {
  // *******************************************

  double xsecw = 1.;

  // Should apply the "once only" weight
  // These are weights that apply to the nominal MC only once
  xsecw *= xsecInfo[i].weight;
#if DEBUG > 0
  if (xsecw < 0 || xsecw > __LARGE_WEIGHT__) {
    std::cerr << "NEGATIVE XSEC ONE-TIME WEIGHT" << std::endl;
    std::cerr << "one-time weight = " << xsecInfo[i].weight << std::endl;
    std::cerr << "Setting to zero!" << std::endl;
    xs2015->printNominalCurrProp();
  }
#if DEBUG == 2
  std::cout << "===============\nPRINTING EVENT BY EVENT WEIGHTS\n===============" << std::endl;
  std::cout << "EVENT " << i << std::endl;
  std::cout << "One-time weight = " << xsecw << std::endl;
#endif
#endif

  // Get the spline weights
  // These are different for CPU and GPU code because we don't call TSpline3->Eval for GPU
  xsecw *= CalcXsecWeight_Spline(i);
#if DEBUG > 0
  if (xsecw < 0 || xsecw > __LARGE_WEIGHT__) {
    std::cerr << "NEGATIVE XSEC SPLINE WEIGHT" << std::endl;
    std::cerr << "spline weight = " << CalcXsecWeight_Spline(i) << std::endl;
    std::cerr << "Setting to zero!" << std::endl;
    xs2015->printNominalCurrProp();
  }
#if DEBUG == 2
  std::cout << "Spline weight = " <<  CalcXsecWeight_Spline(i) << std::endl;
#endif

#endif

  // Get the normalisatiSon weights
  // These are the same for CPU and GPU code
  xsecw *= CalcXsecWeight_Norm(i);
#if DEBUG > 0
  if (xsecw < 0 || xsecw > __LARGE_WEIGHT__) {
    std::cerr << "NEGATIVE XSEC NORM WEIGHT" << std::endl;
    std::cerr << "norm weight = " << CalcXsecWeight_Norm(i) << std::endl;
    std::cerr << "Setting to zero!" << std::endl;
    xs2015->printNominalCurrProp();
  }
#endif

  // Do a check on the cross-section weight
  // If negative (maybe we've gone outside the spline interpolation), set to zero
  if (xsecw < 0) {
    xsecw = 0;
    // If large, cap at 100
  } else if (xsecw > __LARGE_WEIGHT__) {
    xsecw = __LARGE_WEIGHT__;
  }

#if DEBUG == 2
  std::cout << "Total xsec weight = " << xsecw << std::endl;
#endif

  return xsecw;
}

#if USE_SPLINE > USE_TSpline3_red
// ***************************************************************************
// Calculate the TF1 weight for one event i
// Only good for 2019a parameterisation
double samplePDFND2019::CalcXsecWeight_Spline(const int i) {
  // ***************************************************************************
  // The returned cross-section weight
  double xsecw = 1.0;
  // Loop over the spline parameters and get their responses
  for (int id = 0; id < nSplineParams; id++) {
    // Check that the TF1 exists
    if (xsecInfo[i].GetFunc(id)) {
      // Get the global index in cross-section parameterisation for this spline parameter
      int GlobalIndex = splineParsIndex[id];
      // Then get the variation we want to reweight to
      double xvar = xs2015->calcReWeightFrac(GlobalIndex);
      // Calculate xsecw by evaluating the TF1
      xsecw *= xsecInfo[i].Eval(id, xvar);
    }
  } // End the for loop
  return xsecw;
}

#else
// ***************************************************************************
// Calculate the spline weight for one event i
// Only good for 2017a parameterisation
double samplePDFND2019::CalcXsecWeight_Spline(const int i) {
  // ***************************************************************************
  double xsecw = 1.0;
  // Loop over the spline parameters and get their responses
  for (int spline = 0; spline < nSplineParams; spline++) {
    xsecw *= FastSplineEval((xsecInfo[i].GetFunc(spline)), spline);
  } // End the for loop
  return xsecw;
}
#endif

// ***************************************************************************
// Calculate the normalisation weight for one event
// Only good for 2019a parameterisation
double samplePDFND2019::CalcXsecWeight_Norm(const int i) {
  // ***************************************************************************

  double xsecw = 1.0;

  // Do the mode based norm rw
  if (xsecInfo[i].normbin >= 0) {
    xsecw *= xs2015->calcReWeight(xsecInfo[i].normbin);
#if DEBUG > 0
    EventsPerNorm[xsecInfo[i].normbin]++;
#endif
  }
#if DEBUG == 2
  std::cout << "Normbin weight = " << xsecw << std::endl;
#endif

// For NIWG 2019, the base model changed from RFG->SF, so BeRPA is no longer relevant. To keep some freedom in Q2, new normalisation parameters were introduced
  if (xsecInfo[i].mode == kMaCh3_CCQE && (xsecInfo[i].target == kTarget_C || xsecInfo[i].target == kTarget_O)) {
    if(xsecInfo[i].Q2>0.0 && xsecInfo[i].Q2<0.05){
      xsecw *= xs2015->calcReWeight(Q2_norm_0_pos);
#if DEBUG > 0
      EventsPerNorm[Q2_norm_0_pos]++;
#endif
    } else if(xsecInfo[i].Q2>0.05 && xsecInfo[i].Q2<0.1){
      xsecw *= xs2015->calcReWeight(Q2_norm_1_pos);
#if DEBUG > 0
      EventsPerNorm[Q2_norm_1_pos]++;
#endif
    } else if(xsecInfo[i].Q2>0.1 && xsecInfo[i].Q2<0.15){
      xsecw *= xs2015->calcReWeight(Q2_norm_2_pos);
#if DEBUG > 0
      EventsPerNorm[Q2_norm_2_pos]++;
#endif
    } else if(xsecInfo[i].Q2>0.15 && xsecInfo[i].Q2<0.2){
      xsecw *= xs2015->calcReWeight(Q2_norm_3_pos);
#if DEBUG > 0
      EventsPerNorm[Q2_norm_3_pos]++;
#endif
    } else if(xsecInfo[i].Q2>0.2 && xsecInfo[i].Q2<0.25){
      xsecw *= xs2015->calcReWeight(Q2_norm_4_pos);
#if DEBUG > 0
      EventsPerNorm[Q2_norm_4_pos]++;
#endif
    }

  }

  // For NIWG 2019a the 2p2h norm is applied to both C and O interactions, separately for neutrino and anti-neutrino
  // O interactions furthermore get a "CtoO" normalisation, which is the same for neutrino and anti-neutrino interactions
  // Apply the 2p2h norm neutrino for 2p2h neutrino events on oxygen
  if (xsecInfo[i].normbin == mec_norm_CtoO_pos) {
    // If neutrino, apply neutrino 2p2h
    if (xsecInfo[i].species > 0) {
      xsecw *= xs2015->calcReWeight(mec_norm_nu_pos);
#if DEBUG > 0
      EventsPerNorm[mec_norm_nu_pos]++;
#endif
      // Or apply the 2p2h norm anti-neutrino for 2p2h anti-neutrino events on oxygen
    } else {
      xsecw *= xs2015->calcReWeight(mec_norm_anu_pos);
#if DEBUG > 0
      EventsPerNorm[mec_norm_anu_pos]++;
#endif
    }
  }

  // And electron neutrino CC interactions should additionally get the electron neutrino normalisation applied
  if (xsecInfo[i].species == kNue && (xsecInfo[i].mode == kMaCh3_2p2h || xsecInfo[i].mode == kMaCh3_CCMisc || xsecInfo[i].mode < kMaCh3_NC1pi0)) {
    xsecw *= xs2015->calcReWeight(nue_norm_nu_pos);
#if DEBUG > 0
    EventsPerNorm[nue_norm_nu_pos]++;
#endif
    // And electron anti-neutrino CC interactions should additionally get the electron anti-neutrino normalisation applied
  } else if (xsecInfo[i].species == kNue_bar && (xsecInfo[i].mode == kMaCh3_2p2h || xsecInfo[i].mode == kMaCh3_CCMisc || xsecInfo[i].mode < kMaCh3_NC1pi0)) {
    xsecw *= xs2015->calcReWeight(nue_norm_anu_pos);
#if DEBUG > 0
    EventsPerNorm[nue_norm_anu_pos]++;
#endif
  }

  // WP: 2018 CC normalisations
  // Apply only to CC events, (< kMaCh3_NC1pi0 || == kMaCh3_2p2h)
  // With energy range 0.4 to 0.6 GeV
  // And target C or O
  bool isCC = (xsecInfo[i].mode < kMaCh3_NC1pi0 || xsecInfo[i].mode == kMaCh3_2p2h || xsecInfo[i].mode == kMaCh3_CCMisc) 
    && (xsecInfo[i].target == kTarget_C || xsecInfo[i].target == kTarget_O)
    && (0.4 <= xsecInfo[i].Enu && xsecInfo[i].Enu <= 0.6);
  // Then get if it's neutrino or anti-neutrino
  bool isCCnu = (isCC && xsecInfo[i].species > 0);
  bool isCCnub = (isCC && xsecInfo[i].species < 0);
  if (isCCnu) {
    xsecw *= xs2015->calcReWeight(CC_nu_pos);
#if DEBUG > 0
    EventsPerNorm[CC_nu_pos]++;
#endif
  } else if (isCCnub) {
    xsecw *= xs2015->calcReWeight(CC_nubar_pos);
#if DEBUG > 0
    EventsPerNorm[CC_nubar_pos]++;
#endif
  }

  // Applies to CC Misc events only 
  if(xsecInfo[i].mode == kMaCh3_CCMisc && (xsecInfo[i].target == kTarget_C || xsecInfo[i].target == kTarget_O) && (xsecInfo[i].species == kNumu || xsecInfo[i].species == kNumubar || xsecInfo[i].species == kNue || xsecInfo[i].species == kNuebar) ) {
    xsecw *= xs2015->calcReWeight(CC_misc_norm_pos);
#if DEBUG > 0
    EventsPerNorm[CC_misc_norm_pos]++;
#endif
  }
  
  // CC DIS and Multi Pi normalisations
  if((xsecInfo[i].mode == kMaCh3_CCMpi || xsecInfo[i].mode == kMaCh3_CCDIS) && (xsecInfo[i].target == kTarget_C || xsecInfo[i].target == kTarget_O) && (xsecInfo[i].species == kNumu || xsecInfo[i].species == kNumubar || xsecInfo[i].species == kNue || xsecInfo[i].species == kNuebar)) {
    if (xsecInfo[i].species > 0) {
      xsecw *= xs2015->calcReWeight(CC_dis_multipi_norm_nu_pos);
#if DEBUG > 0
      EventsPerNorm[CC_dis_multipi_norm_nu_pos]++;
#endif
    }
    else if (xsecInfo[i].species < 0) {
      xsecw *= xs2015->calcReWeight(CC_dis_multipi_norm_nubar_pos);
#if DEBUG > 0
      EventsPerNorm[CC_dis_multipi_norm_nubar_pos]++;
#endif
    }
  }

  return xsecw;
}

// ***************************************************************************
// Get memory, which is probably silly
int samplePDFND2019::getValue(){ //Note: this value is in KB!
  // ***************************************************************************
  FILE* file = fopen("/proc/self/status", "r");
  int result = -1;
  char line[128];

  while (fgets(line, 128, file) != NULL){
    if (strncmp(line, "VmSize:", 7) == 0){
      result = parseLine(line);
      break;
    }
  }
  fclose(file);
  return result;
}
// ***************************************************************************


// ***************************************************************************
// Get memory, which is probably silly
int samplePDFND2019::parseLine(char* line){
  // ***************************************************************************
  int i = strlen(line);
  while (*line < '0' || *line > '9') line++;
  line[i-3] = '\0';
  i = atoi(line);
  return i;
}


// ***************************************************************************
// Here we set the cross-section covariance matrix for the samplePDF
// Essentially involves: Checking the input root file for the covarianceXsec2015 class
//                       Finding what parameters are normalisation parameters
//                       Finding what parameters are spline parameters (actually doesn't happen here in the code, see instead fillReweightingBins()
//                       Once all the normalisation parameters are found, set their modes (e.g. CC coherent normalisation should only apply to coherent interaction modes)
//                       Most of this is done automatically now, other than setting what xsec model it is by looking at the size (could probably be made even better by string comparison on the input root file name!) and the modes (could also probably be done by string comparison, e.g. if (string.find("CCQE") != std::string::npos) norm_mode = kMaCh3_CCQE)
void samplePDFND2019::setXsecCov(covarianceXsec2015 * const xs, bool norms) {
  // ***************************************************************************

  // Set the xs2015 var
  xs2015 = xs;

  // Get the NIWG model
  xsecmodel = int(xs->getNIWGmodel());
  if (xsecmodel != kNIWG_2020a) {
    std::cerr << "DoItSmart input covariance matrix only supported for NIWG 2019a so far" << std::endl;
    std::cerr << "Can't continue without new implementation" << std::endl;
    std::cout << "Model: " << NIWGmodel_ToString((NIWGmodel_enum)xsecmodel) << std::endl;
    std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
    throw;
  }

  // Get the number of normalisation parameters
  nxsec_norm_modes = xs->getNumND280NormParams();

  // Get the names from xs
  xsec_norm_names = xs2015->getND280NormParsNames();

  // Index for norm parameters in the total number of cross-section parameters (e.g. startbin 2 for 2p2h C means the 3 parameter is a normalisation parameter for 2p2h on Carbon)
  xsec_norm_startbin  = new int[nxsec_norm_modes];

  // Set up which modes the normalisation parameters should apply to (this is only filled for DoItSmart)
  xsec_norm_modes     = new int[nxsec_norm_modes];

  // Get the map between the normalisation parameters index (out of the total), their name, what mode they should apply to, what target, and what PDG
  normParsModes  = xs2015->getND280NormParsModes();
  normParsIndex  = xs2015->getND280NormParsIndex();
  normParsTarget = xs2015->getND280NormParsTarget();
  normParsPDG    = xs2015->getND280NormParsPDG();

  // Have a old-school iterator which counts up to nxsec_norm_modes (i.e. iterators our pointer xsec_norm_modes and others)
  for (int i = 0; i < nxsec_norm_modes; ++i) {
    xsec_norm_startbin[i] = normParsIndex.at(i);
    xsec_norm_modes[i]    = normParsModes.at(i);
  }

  // Number of spline parameters
  nSplineParams   = xs2015->getNumND280SplineParams();
  splineParsIndex = xs2015->getND280SplineParsIndex();
  splineParsNames = xs2015->getND280SplineParsNames();

  // Number of unique splines (neutrino and anti-neutrino can share spline name but have a different parameter in covarianceXsec2015)
  nSplineParamsUniq   = xs2015->getNumSplineParamsUniq();
  splineParsUniqIndex = xs2015->getSplineParsUniqIndex();
  splineParsUniqNames = xs2015->getSplineParsUniqNames();

  // Search through which parameters exist in splinePars and don't exist in splineParsUniq
  nSplineParamsShare    = xs2015->getNumSplineParamsShare();
  splineParsShareIndex  = xs2015->getSplineParsShareIndex();
  splineParsShareToUniq = xs2015->getSplineParsShareToUniq();
  splineParsShareNames  = xs2015->getSplineParsShareNames();

  // Find the functional parameters (Q2 and Eb for 2019)
  nFuncParams = xs2015->getNumND280FuncParams();
  funcParsNames = xs2015->getND280FuncParsNames();
  funcParsIndex = xs2015->getND280FuncParsIndex();

  // Keep the bit-field first entry
  for (int i = 0; i < xs2015->getNumParams(); i++) {
    xsecBitField.push_back(xs2015->getXSecParamID(i,1));
  }

  // Output the normalisation parameters as a sanity check!
  std::cout << "Normalisation parameters: " << std::endl;
  for (int i = 0; i < nxsec_norm_modes; ++i) {
    std::cout << std::setw(5) << i << "|" << std::setw(5) << xsec_norm_startbin[i] << "|" << std::setw(30) << xsec_norm_names[i] << "|" <<  std::setw(5) << xsec_norm_modes[i] << std::endl;
  }

  std::cout << "Total spline parameters: " << std::endl;
  for (int i = 0; i < nSplineParams; ++i) {
    std::cout << std::setw(5) << i << "|" << std::setw(5) << splineParsIndex.at(i) << "|" << std::setw(30) << splineParsNames.at(i) << std::endl;
  }

  std::cout << "Unique spline parameters: " << std::endl;
  for (int i = 0; i < nSplineParamsUniq; ++i) {
    std::cout << std::setw(5) << i << "|" << std::setw(5) << splineParsUniqIndex.at(i) << "|" << std::setw(30) << splineParsUniqNames.at(i) << std::endl;
  }

  std::cout << "Shared spline parameters: " << std::endl;
  for (int i = 0; i < nSplineParamsShare; ++i) {
    std::cout << std::setw(5) << i << "|" << std::setw(5) << splineParsShareIndex.at(i) << "|" << std::setw(30) << splineParsShareNames.at(i) << std::endl;
  }

#if USE_SPLINE < USE_TF1
  // Initialise
  SplineInfoArray = new FastSplineInfo[nSplineParams];
  for (int i = 0; i < nSplineParams; ++i) {
    SplineInfoArray[i].nPts = -999;
    SplineInfoArray[i].xPts = NULL;
    SplineInfoArray[i].CurrSegment = -999;
  }
#endif

  // Clear cross-section normalisation parameters
  // Put the normalisation parameters into the struct, defined in Structs.h
  xsec_norms.clear();
  xsec_norms.reserve(nxsec_norm_modes);
  for (int i = 0; i < nxsec_norm_modes; i++) {
    XsecNorms2 tmp_norm;

    TAxis *temp = NULL;
    tmp_norm.ebins = temp;
    tmp_norm.name       = xsec_norm_names[i];
    tmp_norm.mode       = xsec_norm_modes[i];
    tmp_norm.pdgs       = normParsPDG[i];
    tmp_norm.targets    = normParsTarget[i];
    tmp_norm.startbin   = xsec_norm_startbin[i];
    xsec_norms.push_back(tmp_norm);
  }

  // Check that we've actually set a cross-section model
  if (xsecmodel == kNIWG_undefined) {
    std::cerr << "Unable to identify the NIWG model from the input covariance matrix!" << std::endl;
    std::cerr << "This is probably the sizes not matching the hard-coded ID" << std::endl;
    std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
    throw;
  }

  std::cout << "Finished setXsecCov in samplePDFND2019, xsecmodel = " << NIWGmodel_ToString(NIWGmodel_enum(xsecmodel)) << std::endl;


} // End setting of setXsecCov for normalisation parameters


// ***************************************************************************
// Get kinematic variables through psyche
double samplePDFND2019::getKinematicVariable(ND280KinematicTypes variable, AnaEventB* event) {
  // ***************************************************************************

  // Beware that some of these already exist in the trees that have gone through processing in T2KRW
  // As of Nov 2016 the sample_sum TTree contains: (see T2KRW/app/genWeightsFromNRooTracker_BANFF_*.cxx)
  //
  // Enu
  // Pmu
  // CosThetaMu
  // TruePmu
  // TrueCosThetamu
  // Q2
  // Q2QE
  // NeutrinoCode

  // The lepton candidate
  AnaEventSummaryB* Summary = dynamic_cast<AnaEventSummaryB*>(event->Summary);
  // The sample which this event falls under
  SampleId::SampleEnum index = Summary->EventSample;
  // Check we want this index
  if (samplepdfs->At(index) == NULL || datapdfs->At(index) == NULL) {
    std::cerr << "getKinematicVariable was passed sample " << SampleId::SampleEnum(index) << "(" << index << ") without sample being enabled" << std::endl;
    return 0.0;
  }
  AnaParticleMomB* LeptonCandidate = dynamic_cast<AnaParticleMomB*>(Summary->LeptonCandidate[index]);

  // Check the variable
  switch (variable) {

    // Reconsturcted lepton momentum from psyche
    case (kLeptonMomentum):
      return LeptonCandidate->Momentum;
      break;

      // Reconstructed lepton cos theta from psyche
    case (kLeptonCosTheta):
      return LeptonCandidate->DirectionStart[2];
      break;

      // Reconstructed Lepton theta from psyche
    case (kLeptonTheta):
      return acos(LeptonCandidate->DirectionStart[2])*180/M_PI;
      break;

      // Neutrino energy from psyche
    case (kNeutrinoEnergy):
      return Summary->TrueVertex[index]->NuEnergy;
      break;

      // Reconstructed Enu assuming CCQE
    case (kErecQE):
      {
        double mn = 939.6; //neutron mass
        double mm = 105.65; //muon mass
        double eb = 25; //binding energy C (27 for O, but good enough)
        double pm = LeptonCandidate[index].Momentum;
        double Em = sqrt(pm*pm+mm*mm);
        double cm = LeptonCandidate[index].DirectionStart[2];

        return ((mn-eb)*Em-(-2*mn*eb+eb*eb+mm*mm)/2.0)/(mn-eb-Em+pm*cm)/1.E3; //translate to GeV
      }
      break;

      // Reconstructed Q2 assuming CCQE
    case (kQ2QE):
      {
        double mn= 939.6; //neutron mass  
        double mm= 105.65; //muon mass    
        double eb= 25; //binding energy C (27 for O, but good enough)       
        double pm= LeptonCandidate->Momentum;
        double Em= sqrt(pm*pm+mm*mm);
        double cm= LeptonCandidate->DirectionStart[2];

        double erec = ((mn-eb)*Em-(-2*mn*eb+eb*eb+mm*mm)/2.0)/(mn-eb-Em+pm*cm);

        return (2*erec*(Em-pm*cm)-mm*mm)/1.E6; //translate to GeV^2
      }
      break;

      // True q0
    case (kq0):
      {
        // The true vertex
        AnaTrueVertexB *vtx = Summary->TrueVertex[index];
        if (vtx == NULL) {
          std::cerr << "Could not find TrueVertex object, sorry" << std::endl;
          std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
          return 0.0;
        }

        // The neutrino energy
        double Enu = vtx->NuEnergy/1.E3;
        // The neutrino PDG
        int NuPDG = vtx->NuPDG;
        // The charged lepton PDG
        int LepPDG = 0;
        if (NuPDG < 0) {
          LepPDG = NuPDG+1;
        } else {
          LepPDG = NuPDG-1;
        }

        // Get the highest momentum muon track
        TLorentzVector Muon = MaCh3Utils::GetHMParticle(Summary, LepPDG);
        // This is true muon momentum
        double Emu = Muon.E();

        double q0 = fabs(Enu - Emu);
        return q0;
      }
      break;

      // True q3
    case (kq3):
      {
        AnaTrueVertexB *vtx = Summary->TrueVertex[index];
        if (vtx == NULL) {
          std::cerr << "Could not find TrueVertex object, sorry" << std::endl;
          std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
          return 0.0;
        }

        // The neutrino energy (GeV)
        double Enu = vtx->NuEnergy/(1.E3);
        // The neutrino PDG
        int NuPDG = vtx->NuPDG;
        // The charged lepton PDG
        int LepPDG = 0;
        if (NuPDG < 0) {
          LepPDG = NuPDG+1;
        } else {
          LepPDG = NuPDG-1;
        }

        // Get the highest momentum muon track (GeV)
        TLorentzVector Muon = MaCh3Utils::GetHMParticle(Summary, LepPDG);

        // Make the Enu TLorentzVector with true neutrino direction
        TLorentzVector EnuVect(MaCh3Utils::ND280NuDir[0], MaCh3Utils::ND280NuDir[1], MaCh3Utils::ND280NuDir[2] , 1);
        EnuVect *= Enu;

        // Get pmu
        double Pmu = Muon.Vect().Mag();
        double CosMu = cos(EnuVect.Vect().Angle(Muon.Vect()));

        double q3 = sqrt(Enu*Enu + Pmu*Pmu - 2*Enu*Pmu*CosMu);
        return q3;
      }
      break;

    default:
      std::cerr << "Couldn't find variable " << variable << " exiting!" << std::endl;
      throw;

  } // end switch

}

// ***************************************************************************
// Get individual sample likelihood
double samplePDFND2019::getSampleLikelihood(int isample) {
  // ***************************************************************************

  double negLogLsample = 0.;

  // Skip disabled samples
  if (datapdfs->At(isample) == NULL || samplepdfs->At(isample) == NULL) return 0.0;

  // 1D loop
  if(ndims[isample]==1) {

    for (int i = 1; i <= ((TH1*)datapdfs->At(isample))->GetNbinsX(); i++) {

      double mc = ((TH1*)samplepdfs->At(isample))->GetBinContent(i);
      double dat = ((TH1*)datapdfs->At(isample))->GetBinContent(i);

      if(mc > 0 && dat > 0) {
        negLogLsample += (mc - dat + dat * TMath::Log(dat/mc));
      } else if(mc > 0 && dat == 0) {
        negLogLsample += mc;
      }

    }

    // 2D loop
  } else if (ndims[isample]==2) {

    for(int j = 0; j <= ((TH2Poly*)datapdfs->At(isample))->GetNumberOfBins(); j++) {
      double mc = ((TH2Poly*)samplepdfs->At(isample))->GetBinContent(j);
      double dat = ((TH2Poly*)datapdfs->At(isample))->GetBinContent(j);//WPpoly
      double w2 = W2Hist[isample]->GetBinContent(j);
      negLogLsample += getTestStatLLH(dat, mc, w2);
    }
  }

  return negLogLsample;
}

// ***************************************************************************
// Get the likelihood
double samplePDFND2019::getLikelihood() {
  // ***************************************************************************

  double negLogL = 0.;

  int isample = 0;
  int i       = 0;
  int j       = 0;

  // This only brings speed increase of 2 for 8 processors but might as well leave in
  // If DEBUG > 0 we also verify that the multi-threaded version equates to single thread (see below)
#ifdef MULTITHREAD
#pragma omp parallel for reduction(+:negLogL) private(isample, i, j)
#endif
  for (isample = 0; isample < SampleId::kNSamples; isample++) {

    double negLogLsample = 0.;

    if (datapdfs->At(isample) == NULL || samplepdfs->At(isample) == NULL) continue;

    // 1D loop
    if(ndims[isample]==1) {

      for (i = 1; i < ((TH1*)datapdfs->At(isample))->GetNbinsX(); ++i) {
        double mc = ((TH1*)samplepdfs->At(isample))->GetBinContent(i+1);
        double dat = ((TH1*)datapdfs->At(isample))->GetBinContent(i+1);
        if(mc > 0 && dat > 0) {
          negLogLsample += (mc - dat + dat * TMath::Log(dat/mc));
        } else if(mc > 0 && dat == 0) {
          negLogLsample += mc;
        }

      }

      // 2D loop
    } else if (ndims[isample]==2) {

      for(j = 0; j <= ((TH2Poly*)datapdfs->At(isample))->GetNumberOfBins(); j++) {
	double mc = ((TH2Poly*)samplepdfs->At(isample))->GetBinContent(j);
	double dat = ((TH2Poly*)datapdfs->At(isample))->GetBinContent(j);//WPpoly
        double w2 = W2Hist[isample]->GetBinContent(j);
        negLogLsample += getTestStatLLH(dat, mc, w2);
      }
    }
    negLogL += negLogLsample;
  } // end sample for loop


  // *************************************************************
  // THREAD SAFE VERIFICATION
  // TESTS IF MULTI-THREADED MODE GIVES DIFFERENT LOGL TO SINGLE-THREADED MODE
  // *************************************************************
#ifdef DEBUG
  // Test single-thread version
  double negLogLsingle = 0.;
  for (isample = 0; isample < SampleId::kNSamples; isample++) {
    if (datapdfs->At(isample) == NULL || samplepdfs->At(isample) == NULL) continue;
    // 1D loop
    if(ndims[isample]==1) {
      for (i = 1; i <= ((TH1*)datapdfs->At(isample))->GetNbinsX(); i++) {
	double mc = ((TH1*)samplepdfs->At(isample))->GetBinContent(i);
	double dat = ((TH1*)datapdfs->At(isample))->GetBinContent(i);
        if(mc > 0 && dat > 0) {
          negLogLsingle += (mc - dat + dat * TMath::Log(dat/mc));
        } else if(mc > 0 && dat == 0) {
          negLogLsingle += mc;
        }
      }
      // 2D loop
    } else if (ndims[isample]==2) {
      for(j = 0; j <= ((TH2Poly*)datapdfs->At(isample))->GetNumberOfBins(); j++) {
	double mc = ((TH2Poly*)samplepdfs->At(isample))->GetBinContent(j);
	double dat = ((TH2Poly*)datapdfs->At(isample))->GetBinContent(j);//WPpoly
        // Can use normal Poisson, Barlow-Beeston or "IceCube" LLH
        double w2 = W2Hist[isample]->GetBinContent(j);
        negLogLsingle += getTestStatLLH(dat, mc, w2);
      }
    }
  } // end sample for loop

  // Arbitrarily defined 1.E-6, but should be sufficient
  if (fabs(negLogL - negLogLsingle) > 1.E-6) {
    std::cout << "negLogL_MP - negLogL_SP > 1.E-6!" << std::endl;
    std::cout << fabs(negLogL - negLogLsingle) << std::endl;
  }
  // *************************************************************
  // END OF THREAD SAFE VERIFICATION
  // *************************************************************
#endif

  return negLogL;
}

// *************************
// Calculate the Barlow-Beeston likelhood contribution from MC statistics
// Assumes the beta scaling parameters are Gaussian distributed
// Follows arXiv:1103.0354 section 5 and equation 8, 9, 10, 11 on page 4/5
// Essentially solves equation 11
// data is data, mc is mc, w2 is Sum(w_{i}^2) (sum of weights squared), which is sigma^2_{MC stats}
double samplePDFND2019::getTestStatLLH(double data, double mc, double w2) {
  // *************************

  // Need some MC
  if (mc == 0) return 0.0;

  // The MC used in the likeliihood calculation
  // Is allowed to be changed by Barlow Beeston beta parameters
  double newmc = mc;

  // Not full Barlow-Beeston or what is referred to as "light": we're not introducing any more parameters
  // Assume the MC has a Gaussian distribution around generated
  // As in https://arxiv.org/abs/1103.0354 eq 10, 11

  // The penalty from MC statistics using Barlow-Beeston
  double penalty = 0;
  if (fTestStatistic == kBarlowBeeston) {
    // Barlow-Beeston uses fractional uncertainty on MC, so sqrt(sum[w^2])/mc
    double fractional = sqrt(w2)/mc;
    // -b/2a in quadratic equation
    double temp = mc*fractional*fractional-1;
    // b^2 - 4ac in quadratic equation
    double temp2 = temp*temp + 4*data*fractional*fractional;
    if (temp2 < 0) {
      std::cerr << "Negative square root in Barlow Beeston coefficient calculation!" << std::endl;
      std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
      throw;
    }
    // Solve for the positive beta
    double beta = (-1*temp+sqrt(temp2))/2.;
    newmc = mc*beta;
    // And penalise the movement in beta relative the mc uncertainty
    if (fractional > 0) penalty = (beta-1)*(beta-1)/(2*fractional*fractional);
    else penalty = 0;
  }

  // Calculate the new Poisson likelihood
  // For Barlow-Beeston newmc is modified, so can only calculate Poisson likelihood after Barlow-Beeston
  // For the Poisson likelihood, this is just the usual calculation
  // For IceCube likelihood, we calculate it later
  double stat = 0;
  // All likelihood calculations may use the bare Poisson likelihood, so calculate here
  if (data == 0) stat = newmc;
  else if (newmc > 0) stat = newmc-data+data*TMath::Log(data/newmc);

  // Also try the IceCube likelihood
  // It does not modify the MC content
  // https://arxiv.org/abs/1901.04645
  // Argüelles, C.A., Schneider, A. & Yuan, T. J. High Energ. Phys. (2019) 2019: 30. https://doi.org/10.1007/JHEP06(2019)030
  // We essentially construct eq 3.16 and take the logarithm
  // in eq 3.16, mu is MC, sigma2 is w2, k is data
  if (fTestStatistic == kIceCube) {
    // If there for some reason is 0 mc uncertainty, return the Poisson LLH
    if (w2 == 0) return stat;

    // Reset the penalties if there is mc uncertainty
    stat = 0.0;
    penalty = 0.0;
    // Auxillary variables
    long double b = mc/w2;
    long double a = mc*b+1;
    long double k = data;
    // Use C99's implementation of log of gamma function to not be C++11 dependent
    stat = -1*(a * logl(b) + lgammal(k+a) - lgammal(k+(long double)1) - ((k+a)*log1pl(b)) - lgammal(a));
  }

  // Return the statistical contribution and penalty
  return stat+penalty;
}


// *************************
// Add data from a vector of doubles
// Only OK for TH1D
void samplePDFND2019::addData(std::vector<double> &dat, int index) {
  // *************************

  if (datapdfs->At(index) == NULL) {
    std::cerr << "Failed to replace data because selection " << SampleId::ConvertSample(SampleId::SampleEnum(index)) << "(" << index << ")" <<  " is not enabled" << std::endl;
    std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
    throw;
  }

  // Check our data is TH1D
  if (ndims[index] != 1) {
    std::cerr << "Failed to replace data because datapdf is not TH1D" << std::endl;
    std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
    throw;
  }

  // Reset the current data
  ((TH1*)datapdfs->At(index))->Reset();

  // Replace the data with new data
  for (int i = 0; i < int(dat.size()); i++) {
    ((TH1*)datapdfs->At(index))->Fill(dat[i]);
  }

  return;
}

// *************************
// Add data from a vector of vector of doubles
void samplePDFND2019::addData(std::vector< vector <double> > &dat, int index) {
  // *************************

  if (datapdfs->At(index) == NULL) {
    std::cerr << "Failed to replace data because selection " << SampleId::ConvertSample(SampleId::SampleEnum(index)) << "(" << index << ")" <<  " is not enabled" << std::endl;
    std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
    throw;
  }

  // Check our data is TH2D
  if (ndims[index] != 2) {
    std::cerr << "Failed to replace data because datapdf is not TH1D" << std::endl;
    std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
    throw;
  }

  ((TH2Poly*)datapdfs->At(index))->Reset("");//An option is required here in the reset function for TH2Polys so have to have ""

  for(int i=0; i<int(dat.size()); i++) {
    ((TH2Poly*)datapdfs->At(index))->Fill(dat[0][i], dat[1][i]);
  }

  return;

}

// *************************
// Add data from TH1D
void samplePDFND2019::addData(TH1D* binneddata, int index){
  // *************************

  if (datapdfs->At(index)) {
    datapdfs->RemoveAt(index);
    datapdfs->AddAt(binneddata,index);
  } else {
    datapdfs->AddAt(binneddata,index);
  }

  return;

}

// *************************
// Add data from TH2D
void samplePDFND2019::addData(TH2Poly* binneddata, int index){
  // *************************

  if(datapdfs->At(index)) {
    datapdfs->RemoveAt(index);
    datapdfs->AddAt(binneddata,index);
  } else {
    datapdfs->AddAt(binneddata,index);
  }

  return;
}

// *************************
// Get the event rate
double samplePDFND2019::getEventRate(int index) {
  // *************************

  if(samplepdfs->At(index)) {
    return NoOverflowIntegral(((TH2Poly*)samplepdfs->At(index))); 
  } else {
    return 0;
  }
}


// ***************************************************************************
// Make the TChain for getting spline information once we've looped over relevant psyche events
TChain* samplePDFND2019::MakeXsecChain() {
  // ***************************************************************************

  // The TTree we want inside each ROOT file is the sample_sum tree, created by looping the psyche events in T2KRW
  TChain* chain = new TChain("sample_sum");

  //Get the file paths as in AnalysisManager
  std::map< std::string, SampleGroup >& sampleGroups = exp->GetSampleGroups();
  std::map< std::string, SampleGroup >::iterator sit;

  // Loop over the sample groups and add files to the TChain
  for (sit = sampleGroups.begin(); sit != sampleGroups.end(); sit++) {

    SampleGroup& sampleGroup = sit->second;
    std::map< std::string, DataSample*>& mcSamples = sampleGroup.GetMCSamples();
    std::map< std::string, DataSample*>::iterator it2;

    for (it2 = mcSamples.begin(); it2 != mcSamples.end(); it2++) {
      DataSample* sample = it2->second;      
      std::cout << "Making TChain from: " << sample->GetFilePath().c_str() << std::endl;
      chain->AddFile(sample->GetFilePath().c_str());
    }      
  }

  // Check that there are entries in the chain
  if (chain->GetEntries() <= 0) {
    std::cerr << "Found no entries in the MC you provided..." << std::endl;
    std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
    throw;
  }

  std::cout << "Total entries in TChain: " << chain->GetEntries() << std::endl;

  // Disable all branches and turn on when we need them
  chain->SetBranchStatus("*", false);

  return chain;
}

#if USE_SPLINE > USE_TSpline3_red
// ***************************************************************************
// TGraph** because each xsecgraph points to a TObjArray and each collection of xsecgraph has nSplines entries
//
// #########################
// WARNING THIS NEEDS TO MATCH THE ENUMERATION IN chain->SetBranchAddress IN samplePDFND2019::fillReweightingBins()
void samplePDFND2019::SetSplines(TGraph** &xsecgraph, const int EventNumber) {
  // ***************************************************************************

  xsecInfo[EventNumber].SetSplineNumber(nSplineParams);

  if (xsecmodel != kNIWG_2020a) {
    std::cerr << "Sorry, I only support NIWG 2019a currently" << std::endl;
    std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
    throw;
  }

#ifdef DEBUG_TF1
  std::stringstream ss;
  ss << EventNumber;
  dumpfile->cd();
#endif

  // Now loop over and set the splines that we have
  // This is definitely the set of unique parameters
  // When we do TSpline3->Eval(paramVal) we don't want paramVal to be the same though!
  for (int j = 0; j < nSplineParams; ++j) {

    // The TF1 object we build from fitting the TGraph
    TF1 *Fitter = NULL;
#if USE_SPLINE == USE_TF1_red
    TF1_red *tf1_red = NULL;
#endif

    // For a valid spline we require it not be null and have more than one point
    if (xsecgraph[j] && xsecgraph[j]->GetN() > 1) {

      // For 2p2h shape C and O we can't fit a polynomial: try a linear combination of two linear functions around 0
      if (j == 3 || j == 4) {
        Fitter = new TF1(xsecgraph[j]->GetName(), "(x<=0)*(1+[0]*x)+(x>0)*([1]*x+1)", xsecgraph[j]->GetX()[0], xsecgraph[j]->GetX()[xsecgraph[j]->GetN()-1]);
      // Fit 5hd order polynomial for all other parameters
      } else {
        Fitter = new TF1(xsecgraph[j]->GetName(), "1+[0]*x+[1]*x*x+[2]*x*x*x+[3]*x*x*x*x+[4]*x*x*x*x*x", xsecgraph[j]->GetX()[0], xsecgraph[j]->GetX()[xsecgraph[j]->GetN()-1]);
      }
      // Fit the TF1 to the graph
      xsecgraph[j]->Fit(Fitter, "Q0");

#if USE_SPLINE == USE_TF1_red
      // Make the reduced TF1 if we want
      tf1_red = new TF1_red(Fitter, j);
#endif

      // Debug the TF1 by comparing it to the TSpline3
#ifdef DEBUG_TF1
      // Make the comparison spline
      TSpline3 *spline = new TSpline3(xsecgraph[j]->GetName(), xsecgraph[j]);

      // The total xrange
      double xrange = xsecgraph[j]->GetX()[xsecgraph[j]->GetN()-1] - xsecgraph[j]->GetX()[0];
      bool checkit = false;
      double point = -999;
      // Scan 100 points and look for differences in weights
      const int npoints = 100;
      for (int z = 0; z < npoints; ++z) {
        double val = xsecgraph[j]->GetX()[0] + z*xrange/double(npoints);
        double fitted = Fitter->Eval(val);
        double splined = spline->Eval(val);
        if (fabs(fitted/splined-1) > 0.05) {
          checkit = true;
          point = val;
          break;
        }
      }

      // If we've found a dodgy event write to file
      if (checkit == true) {
        std::string name = xsecgraph[j]->GetName();
        name = ss.str()+"_"+name;
        xsecgraph[j]->SetNameTitle((name+"_gr").c_str(), (name+"_gr").c_str());
        Fitter->SetNameTitle((name+"_tf").c_str(), (name+"_tf").c_str());
        spline->SetNameTitle((name+"_sp").c_str(), (name+"_sp").c_str());

        TCanvas *temp = new TCanvas("canv","canv", 800,800);
        temp->SetName((name+"_canv").c_str());
        temp->SetTitle((name+"_canv").c_str());
        xsecgraph[j]->SetLineColor(kBlack);
        xsecgraph[j]->SetLineWidth(2);
        xsecgraph[j]->SetLineStyle(kSolid);
        xsecgraph[j]->SetFillStyle(0);
        xsecgraph[j]->SetFillColor(0);
        Fitter->SetLineColor(kRed);
        Fitter->SetLineStyle(kDashed);
        Fitter->SetLineWidth(2);
        spline->SetLineColor(kOrange);
        spline->SetLineStyle(kSolid);
        spline->SetLineWidth(2);

        temp->cd();
        xsecgraph[j]->Draw();
        TLine *line = new TLine(point, 0, point, 20);
        line->SetLineColor(kRed);
        line->SetLineWidth(2);
        xsecgraph[j]->GetXaxis()->SetTitle("Variation (rel nominal)");
        xsecgraph[j]->GetYaxis()->SetTitle("Weight");
        spline->Draw("same");
        Fitter->Draw("same");
        line->Draw("same");
        TLegend *leg = new TLegend(0.65, 0.85, 1.0, 1.0);
        leg->SetHeader(name.c_str());
        leg->AddEntry(xsecgraph[j], "Graph");
        leg->AddEntry(spline, "TSpline3");
        leg->AddEntry(Fitter, "TF1");
        leg->Draw("same");

        dumpfile->cd();
        xsecgraph[j]->Write();
        Fitter->Write();
        spline->Write();
        temp->Write();

        delete temp;
        delete leg;
        delete line;
      }
      delete spline;
#endif

      // For events without the jth TGraph, set to NULL
    } else {
      //Fitter = dummy;
      Fitter = NULL;
#if USE_SPLINE == USE_TF1_red
      tf1_red = NULL;
#endif
    }

    // Now TF1* Fitter is either a TF1 or a pointer to NULL for the jth parameter
    // Save it to the struct
#if USE_SPLINE == USE_TF1_red
    xsecInfo[EventNumber].SetFunc(j, tf1_red);
#elif USE_SPLINE == USE_TF1
    xsecInfo[EventNumber].SetFunc(j, Fitter);
#endif
  } // end for j loop

  // Delete the TGraph pointers
  for (int i = 0; i < nSplineParams; ++i) {
    delete xsecgraph[i];
  }
} // end function

#else
// ***************************************************************************
// TGraph** because each xsecgraph points to a TObjArray and each collection of xsecgraph has nSplines entries
//
// #########################
// WARNING THIS NEEDS TO MATCH THE ENUMERATION IN chain->SetBranchAddress IN samplePDFND2019::fillReweightingBins()
void samplePDFND2019::SetSplines(TGraph** &xsecgraph, const int EventNumber) {
  // ***************************************************************************

  xsecInfo[EventNumber].SetSplineNumber(nSplineParams);

  if (xsecmodel != kNIWG_2020a) {
    std::cerr << "Sorry, I only support NIWG 2019a currently" << std::endl;
    std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
    throw;
  }

  // Now loop over and set the splines that we have
  // This is definitely the set of unique parameters
  // When we do TSpline3->Eval(paramVal) we don't want paramVal to be the same though!
  for (int j = 0; j < nSplineParams; ++j) {

    // Here's the TSpline3
    TSpline3* spline = NULL;
#if USE_SPLINE == USE_TSpline3_red
    TSpline3_red *spline_red = NULL;
#endif

    // For a valid spline we require it not be null and have more than one point
    if (xsecgraph[j] && xsecgraph[j]->GetN() > 1) {

      // Create the TSpline3* from the TGraph* and build the coefficients
      spline = new TSpline3(xsecgraph[j]->GetName(), xsecgraph[j]);
      spline->SetNameTitle(xsecgraph[j]->GetName(), xsecgraph[j]->GetName());
#if USE_SPLINE == USE_TSpline3_red
      // Make the reduced TSpline3 format and delete the old spline
      spline_red = new TSpline3_red(spline, j);
#endif

      // Fill the SplineInfoArray entries with information on each splinified parameter
      if (SplineInfoArray[j].xPts == NULL) {
        // Fill the number of points
#if USE_SPLINE == USE_TSpline3_red
        SplineInfoArray[j].nPts = spline_red->GetNp();
#else
        SplineInfoArray[j].nPts = spline->GetNp();
#endif

        // Fill the x points
        SplineInfoArray[j].xPts = new double[SplineInfoArray[j].nPts];
        for (int k = 0; k < SplineInfoArray[j].nPts; ++k) {
          double xtemp = -999.99;
          double ytemp = -999.99;
#if USE_SPLINE == USE_TSpline3_red
          spline_red->GetKnot(k, xtemp, ytemp);
#else
          spline->GetKnot(k, xtemp, ytemp);
#endif
          SplineInfoArray[j].xPts[k] = xtemp;
        }
      }
      // For events without the jth spline, set to NULL
    } else {
      spline = NULL;
#if USE_SPLINE == USE_TSpline3_red
      spline_red = NULL;
#endif
    }
    // Save the TSpline3* into the struct
#if USE_SPLINE == USE_TSpline3_red
    xsecInfo[EventNumber].SetFunc(j, spline_red);
#else
    xsecInfo[EventNumber].SetFunc(j, spline);
#endif
  } // end for j loop

  // Delete the TGraph pointers
  for (int i = 0; i < nSplineParams; ++i) {
    delete xsecgraph[i];
  }

} // end function

// ***************************************************************************
// TGraph** because each xsecgraph points to a TObjArray and each collection of xsecgraph has nSplines entries
//
// #########################
// WARNING THIS NEEDS TO MATCH THE ENUMERATION IN chain->SetBranchAddress IN samplePDFND2019::fillReweightingBins()
// Reduce number of points in TGraph to 3
void samplePDFND2019::SetSplines_Reduced(TGraph** &xsecgraph, const int EventNumber) {
  // ***************************************************************************

  xsecInfo[EventNumber].SetSplineNumber(nSplineParams);
  if (xsecmodel != kNIWG_2020a) {
    std::cerr << "Sorry, I only support NIWG 2019a currently" << std::endl;
    std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
    throw;
  }

  // First look at the TGraphs
  for (int i = 0; i < nSplineParams; ++i) {

    // Here's the reduced graph
    TGraph* ReducedGraph = NULL;
    // And the TSpline3 we build from that reduced graph
    TSpline3* Spline = NULL;
#if USE_SPLINE == USE_TSpline3_red
    TSpline3_red *spline_red = NULL;
#endif

    // If this criteria is true then the spline response is not flat
    if (xsecgraph[i] && xsecgraph[i]->GetN() > 1) {
      // Get the TGraph for this spline parameter
      TGraph *TempGraph = xsecgraph[i];
      // Now look at the x coordinates
      double *xaxis = TempGraph->GetX();
      // And the responses in y
      double *yresp = TempGraph->GetY();
      // The number of points in the original graph
      int nPoints = TempGraph->GetN();
      // The reduced x and y response (set to be half of original)
      double *x_red = new double[nPoints/2+1];
      double *y_red = new double[nPoints/2+1];
      if (nPoints % 2 != 1) {
        std::cerr << "Reduced spline points method only works for odd number of points, sorry" << std::endl;
        std::cerr << "Change the code in " << __FILE__ << ":" << __LINE__ << std::endl;
        std::cerr << "Or revert to using _NOT_ reduced spline code at ND280 (see SetSplines function in " << __FILE__ << ")" << std::endl;
        throw;
      }
      // Skip every second point
      for (int j = 0; j < (nPoints/2)+1; ++j) {
        // If we're on the last point in the reduced we want to pick out the last point in the total TGraph
        x_red[j] = xaxis[j*2];
        y_red[j] = yresp[j*2];
      }
      // Now make the new TGraphs from these x and y
      ReducedGraph = new TGraph(nPoints/2+1, x_red, y_red);
      ReducedGraph->SetName(xsecgraph[i]->GetName());

      // Now build the TSpline3
      Spline = new TSpline3(ReducedGraph->GetName(), ReducedGraph);
      Spline->SetNameTitle(ReducedGraph->GetName(), ReducedGraph->GetName());
      // Make the reduced TSpline3 format and delete the old spline
#if USE_SPLINE == USE_TSpline3_red
      spline_red = new TSpline3_red(Spline, i);
#endif

      // Fill the SplineInfoArray entries with information on each splinified parameter
      if (SplineInfoArray[i].xPts == NULL) {
        // Fill the number of points
#if USE_SPLINE == USE_TSpline3_red
        SplineInfoArray[i].nPts = spline_red->GetNp();
#else
        SplineInfoArray[i].nPts = Spline->GetNp();
#endif
        // Fill the x points
        SplineInfoArray[i].xPts = new double[SplineInfoArray[i].nPts];
        for (int k = 0; k < SplineInfoArray[i].nPts; ++k) {
          double xtemp = -999.99;
          double ytemp = -999.99;
#if USE_SPLINE == USE_TSpline3_red
          spline_red->GetKnot(k, xtemp, ytemp);
#else
          Spline->GetKnot(k, xtemp, ytemp);
#endif
          SplineInfoArray[i].xPts[k] = xtemp;
        }
      }

      // Delete the temporary memory
      delete[] x_red;
      delete[] y_red;
      delete ReducedGraph;
    } else {
      ReducedGraph = NULL;
      Spline = NULL;
    }
    // Finall save the pointer to the reduced spline in the struct
#if USE_SPLINE == USE_TSpline3_red
    xsecInfo[EventNumber].SetFunc(i, spline_red);
#else
    xsecInfo[EventNumber].SetFunc(i, Spline);
#endif
  } // end loop over xsec parameter

  // Delete the TGraph pointers
  for (int i = 0; i < nSplineParams; ++i) {
    delete xsecgraph[i];
  }
} // end function
#endif

// ***************************************************************************
// Interface with psyche; set up our experiment, loop over the events (data and MC), associate cross-section splines with events
// Needs to have addSelection called!
void samplePDFND2019::fillReweigtingBins() {
  // ***************************************************************************

  // Check that the covariances are filled
  CheckCovariances();

  // Start a little clock for timing
  TStopwatch clock;
  clock.Start();

  // General thing that happens here:
  // 1) Loop over events in the Highland flat tree
  // 2) See if it passes selection by interfacing to psyche and getting systematics
  //    2a) psyche systematics for ND280 and flux are ONLY used if we have detsimple = 0 (eats a tonne of RAM, 60GB-ish and is very slow, so is not the default)
  // 3) Save the characteristics of the events (e.g. momentum, weights, splines) in Structs so that we can loop over them later with the reweight call

  // nXsecSplines is set based on the TString prod
  // We have also set a covariance matrix which we can count the number of spline parameters
  nXsecSplinesFromCov = xs2015->getNumSplineParamsUniq();
  if (nXsecSplinesFromCov!= nXsecSplines) {
    std::cerr << "nXsecSplinesFromCov != nXsecSplines, THIS WILL CAUSE PAIN AND SUFFERING!" << std::endl;
    std::cerr << "nXsecSplines is hard-coded in samplePDFND2019 constructor to force compliance, nXsecSplinesFromCov is read from the input covariance matrix" << std::endl;
    std::cerr << "nXsecSplinesFromCov = " << nXsecSplinesFromCov << std::endl;
    std::cerr << "nXsecSplines = "  << nXsecSplines << std::endl;
    std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
    throw;
  }

  // Here are all the things we need when we read in from the TChain and use SetBranchAddress
  // 18 is the number of cross-section parameters that come in splines
  // Is 17 sometimes, ha
  // For NIWG2015... we need two extra parameters in the splines because of SCCV and SCCA splines (we aren't actually used!)
  if (xsecmodel < kNIWG_2017a) {
    nXsecSplines += 2;
  }

  QuietPlease(); 
  // Initialise the psyche MC TChain
  _man.InitialiseMCTChain();
  NowTalk();

  // Pointer to the sample_sum TChain we're using (created through looping over psyche events in T2KRW/app/genWeights*BANFF*.cxx
  // Warning, uses new!
  TChain* chain = NULL;

  //These are arrays as previously for each event we had a spline (and set of variables, enu, q2, fluxweight etc which came from the spline file) per sample. For the 2020 analysis, splines were made only for 1 sample per event, but still with the objects all being vectors. So we just grab the 0th element of them. 
  TGraph** xsecgraph = new TGraph*[nXsecSplines];
  for(int i = 0; i < nXsecSplines; i++) {
    xsecgraph[i] = NULL;
  }
  TObjArray** grapharrays = new TObjArray*[nXsecSplines];
  for(int i = 0; i < nXsecSplines; i++) {
    grapharrays[i] = NULL;
  }

  // Number of samples in psyche
  const unsigned int nPsycheSamples = SampleId::kNSamples;

  // The event-level variables, sort of
  double enu_a          [nPsycheSamples];
  double q2_a           [nPsycheSamples];

  // Different weights read in from the flat tree (used in SetBranchAddress and subsequent GetEntry)
  double fluxweight_a   [nPsycheSamples];
  double detweight_tree = -999.99;

  // Interaction specifics
  int mode_a    [nPsycheSamples];
  int species_a [nPsycheSamples];
  int truthv_a  [nPsycheSamples];
  int target_a  [nPsycheSamples];

  // Get the cross-section splines that are saved for the events through T2KRW
  if (xs2015) {

    // Make the xsec chain (calls a new)
    // Also disables all branches
    chain = MakeXsecChain();

    if (xsecmodel == kNIWG_2020a) {

      // Loop over the spline parameters as determined by the xsec covariance 2015 class
      // Do this to set up the TGraph points that get fed to our TSpline3
      for (int i = 0; i < nSplineParamsUniq; i++) {

        // The name of the current spline parameter (enables matching input covariance with input spline)
        // This way we enforce an enum on the structure; namely it follows whatever the input covariance matrix structure follows for parameters that have id > 0
        std::string paramName = splineParsUniqNames.at(i);
        int detOn = xsecBitField.at(splineParsUniqIndex.at(i));

        if (!detOn&1) {
          std::cerr << "Found a spline parameter which is set to be off, is this correct?" << std::endl;
          std::cerr << "On spline parameter " << i << " which has global " << splineParsUniqIndex.at(i) << " (" << xs2015->getParameterName(splineParsUniqIndex.at(i)) << ")" << std::endl;
          std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
          throw;
        }

        // Be super safe here and check the names once
        // Hard-coded names here, be careful when implementations change!

        // The CCQE parameters
        if          (paramName == "MAQE") {
          chain->SetBranchStatus( "MAQEGraph", true);
          chain->SetBranchAddress("MAQEGraph",        &(grapharrays[i]));

          // 2p2h parameters
        } else if   (paramName == "2p2h_shape_C") {
          chain->SetBranchStatus( "PDD_CGraph", true);
          chain->SetBranchAddress("PDD_CGraph",       &(grapharrays[i]));
        } else if   (paramName == "2p2h_shape_O") {
          chain->SetBranchStatus( "PDD_OGraph", true);
          chain->SetBranchAddress("PDD_OGraph",       &(grapharrays[i]));

	  // 2p2h E Dep parameters
	} else if     (paramName == "2p2h_Edep_lowEnu") {
	  chain->SetBranchStatus( "MEC_lowEnuGraph", true);
	  chain->SetBranchAddress("MEC_lowEnuGraph",  &(grapharrays[i]));
	} else if     (paramName == "2p2h_Edep_highEnu") {
	  chain->SetBranchStatus( "MEC_highEnuGraph", true);
	  chain->SetBranchAddress("MEC_highEnuGraph",  &(grapharrays[i]));
	} else if     (paramName == "2p2h_Edep_lowEnubar") {
          chain->SetBranchStatus( "MEC_lowEnubarGraph", true);
          chain->SetBranchAddress("MEC_lowEnubarGraph",  &(grapharrays[i]));
	} else if     (paramName == "2p2h_Edep_highEnubar") {
          chain->SetBranchStatus( "MEC_highEnubarGraph", true);
          chain->SetBranchAddress("MEC_highEnubarGraph",  &(grapharrays[i]));
	
          // Single pion parameters
        } else if     (paramName == "CA5") {
          chain->SetBranchStatus( "CA5Graph", true);
          chain->SetBranchAddress("CA5Graph",         &(grapharrays[i]));
        } else if   (paramName == "MARES") {
          chain->SetBranchStatus( "MARESGraph", true);
          chain->SetBranchAddress("MARESGraph",       &(grapharrays[i]));
        } else if   (paramName == "ISO_BKG_LowPPi") {
          chain->SetBranchStatus( "BgSclRes_lowPPiGraph", true);
          chain->SetBranchAddress("BgSclRes_lowPPiGraph",    &(grapharrays[i]));
        } else if   (paramName == "ISO_BKG") {
	chain->SetBranchStatus( "BgSclResGraph", true);
	chain->SetBranchAddress("BgSclResGraph",    &(grapharrays[i]));

          // Multi-pion and DIS parameter
        } else if     (paramName == "CC_BY_DIS") {
          chain->SetBranchStatus( "DIS_BY_corrGraph", true);
          chain->SetBranchAddress("DIS_BY_corrGraph",   &(grapharrays[i]));
        } else if     (paramName == "CC_BY_MPi") {
	  chain->SetBranchStatus( "MultiPi_BYGraph", true);
	  chain->SetBranchAddress("MultiPi_BYGraph",   &(grapharrays[i]));
        } else if     (paramName == "CC_AGKY_Mult") {
	  chain->SetBranchStatus( "MultiPi_Xsec_AGKYGraph", true);
 	  chain->SetBranchAddress("MultiPi_Xsec_AGKYGraph",   &(grapharrays[i]));

	  // FSI parameters
        } else if   (paramName == "FEFQE") {
          chain->SetBranchStatus( "FSI_INEL_LOGraph", true);
          chain->SetBranchAddress("FSI_INEL_LOGraph", &(grapharrays[i]));
        } else if   (paramName == "FEFQEH") {
          chain->SetBranchStatus( "FSI_INEL_HIGraph", true);
          chain->SetBranchAddress("FSI_INEL_HIGraph", &(grapharrays[i]));
        } else if   (paramName == "FEFINEL") {
          chain->SetBranchStatus( "FSI_PI_PRODGraph", true);
          chain->SetBranchAddress("FSI_PI_PRODGraph", &(grapharrays[i]));
        } else if   (paramName == "FEFABS") {
          chain->SetBranchStatus( "FSI_PI_ABSGraph", true);
          chain->SetBranchAddress("FSI_PI_ABSGraph",  &(grapharrays[i]));
        } else if   (paramName == "FEFCX") {
          chain->SetBranchStatus( "FSI_CEX_LOGraph", true);
          chain->SetBranchAddress("FSI_CEX_LOGraph",  &(grapharrays[i]));
        } //else if   (paramName == "FSI_CEX_HI") {
        //chain->SetBranchStatus( "FSI_CEX_HIGraph", true);
        //chain->SetBranchAddress("FSI_CEX_HIGraph",  &(grapharrays[i]));

        // Any other parameter should indicate that we don't have a spline for it, which is a serious error!
        //}
        else {
          std::cerr << "Trying to match an input spline parameter with a spline that doesn't seem to exist in the input tree" << std::endl;
          std::cerr << "This matching is hard-coded though, so it might because it isn't updated!!!" << std::endl;
          std::cerr << "Parameter name = " << paramName << std::endl;
          std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
          throw;
        }
      } // End the for loop over the unique spline parameters
    } else {
      std::cerr << "Not supporting the xsecmodel you've input" << std::endl;
      std::cerr << "xsecmodel = " << NIWGmodel_enum(xsecmodel) << " = " << NIWGmodel_ToString(NIWGmodel_enum(xsecmodel)) << std::endl;
      std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
      throw;
    }

    // Now set up the other branches in the file
    // Interaction specifics
    chain->SetBranchStatus( "Enu", true);
    chain->SetBranchAddress("Enu",      &enu_a);
    chain->SetBranchStatus( "Q2", true);
    chain->SetBranchAddress("Q2",       &q2_a);
    chain->SetBranchStatus( "TgtMat", true);
    chain->SetBranchAddress("TgtMat",   &target_a);

    // NEUT specifics
    // Interaction code
    chain->SetBranchStatus( "ReactionCode", true);
    chain->SetBranchAddress("ReactionCode",     &(mode_a));
    // Incoming neutrino PDG
    chain->SetBranchStatus( "NeutrinoCode", true);
    chain->SetBranchAddress("NeutrinoCode",     &(species_a));


    // Psyche specifics
    chain->SetBranchStatus( "TruthVtx", true);
    chain->SetBranchAddress("TruthVtx",   &truthv_a);
    // Flux weight for each event
    chain->SetBranchStatus( "FluxWeight", true);
    chain->SetBranchAddress("FluxWeight", &fluxweight_a);
    // Detector weight for each event
    chain->SetBranchStatus( "DetNomWeight", true);
    chain->SetBranchAddress("DetNomWeight", &detweight_tree);
  } // end if xsec2015

  W2Hist.reserve(SampleId::kNSamples);
  for (int i = 0; i < SampleId::kNSamples; ++i) {
    if (samplepdfs->At(i) == NULL) {
      W2Hist.push_back(NULL);
      continue;
    } else {
      // Holds raw number of MC events
      W2Hist.push_back(((TH2Poly*)(samplepdfs->At(i)->Clone())));
      W2Hist[i]->SetDirectory(0);
      W2Hist[i]->Reset("");
      W2Hist[i]->Fill(0.0, 0.0, 0.0);
      W2Hist[i]->SetTitle((std::string(W2Hist[i]->GetTitle())+"_W2Hist").c_str());
      W2Hist[i]->SetName(W2Hist[i]->GetTitle());
    }
  }


  // Initialise the debugging TH1D
#if DEBUG > 0
  nEventsInTree = _man.GetNMCEvents();
  MCOnly    .reserve(SampleId::kNSamples);
  POTOnly   .reserve(SampleId::kNSamples);
  FluxOnly  .reserve(SampleId::kNSamples);
  XsecOnly  .reserve(SampleId::kNSamples);
  DetOnly   .reserve(SampleId::kNSamples);
  NDCovOnly .reserve(SampleId::kNSamples);
  BeamCovOnly.reserve(SampleId::kNSamples);
  AllOnly   .reserve(SampleId::kNSamples);
  NuMuPmu   .reserve(SampleId::kNSamples);

  for (int i = 0; i < SampleId::kNSamples; ++i) {
    if (samplepdfs->At(i) == NULL) continue;
    // Holds raw number of MC events
    MCOnly[i] = (((TH2Poly*)(samplepdfs->At(i)->Clone())));
    MCOnly[i]->SetDirectory(0);
    MCOnly[i]->Reset("");
    MCOnly[i]->SetTitle((std::string(MCOnly[i]->GetTitle())+"_MC").c_str());
    MCOnly[i]->SetName(MCOnly[i]->GetTitle());
    // Holds POT scaled number of MC events
    POTOnly[i] = (((TH2Poly*)(samplepdfs->At(i)->Clone())));
    POTOnly[i]->SetDirectory(0);
    POTOnly[i]->Reset("");
    POTOnly[i]->SetTitle((std::string(POTOnly[i]->GetTitle())+"_POT").c_str());
    POTOnly[i]->SetName(POTOnly[i]->GetTitle());
    // Holds Flux scaled number of MC events
    FluxOnly[i] = (((TH2Poly*)(samplepdfs->At(i)->Clone())));
    FluxOnly[i]->Reset("");
    FluxOnly[i]->SetDirectory(0);
    FluxOnly[i]->SetTitle((std::string(FluxOnly[i]->GetTitle())+"_FLUX").c_str());
    FluxOnly[i]->SetName(FluxOnly[i]->GetTitle());
    // Holds xsec+POT scaled number of MC events
    XsecOnly[i] = (((TH2Poly*)(samplepdfs->At(i)->Clone())));
    XsecOnly[i]->Reset("");
    XsecOnly[i]->SetDirectory(0);
    XsecOnly[i]->SetTitle((std::string(XsecOnly[i]->GetTitle())+"_XSEC").c_str());
    XsecOnly[i]->SetName(XsecOnly[i]->GetTitle());
    // Holds det scaled number of MC events
    DetOnly[i] = (((TH2Poly*)(samplepdfs->At(i)->Clone())));
    DetOnly[i]->Reset("");
    DetOnly[i]->SetDirectory(0);
    DetOnly[i]->SetTitle((std::string(DetOnly[i]->GetTitle())+"_DET").c_str());
    DetOnly[i]->SetName(DetOnly[i]->GetTitle());
    // Holds ND280 cov+POT scaled number of MC events
    NDCovOnly[i] = (((TH2Poly*)(samplepdfs->At(i)->Clone())));
    NDCovOnly[i]->Reset("");
    NDCovOnly[i]->SetDirectory(0);
    NDCovOnly[i]->SetTitle((std::string(NDCovOnly[i]->GetTitle())+"_NDCOV").c_str());
    NDCovOnly[i]->SetName(NDCovOnly[i]->GetTitle());
    // Holds beam cov+POT scaled number of MC events
    BeamCovOnly[i] = (((TH2Poly*)(samplepdfs->At(i)->Clone())));
    BeamCovOnly[i]->Reset("");
    BeamCovOnly[i]->SetDirectory(0);
    BeamCovOnly[i]->SetTitle((std::string(BeamCovOnly[i]->GetTitle())+"_BEAMCOV").c_str());
    BeamCovOnly[i]->SetName(BeamCovOnly[i]->GetTitle());
    // Holds raw number of MC events
    AllOnly[i] = (((TH2Poly*)(samplepdfs->At(i)->Clone())));
    AllOnly[i]->Reset("");
    AllOnly[i]->SetDirectory(0);
    AllOnly[i]->SetTitle((std::string(AllOnly[i]->GetTitle())+"_ALL").c_str());
    AllOnly[i]->SetName(AllOnly[i]->GetTitle());

    // Holds scaled event rates matching numu
    std::string Title = ((TH2Poly*)(samplepdfs->At(i)))->GetTitle();
    Title += "_NuMu";
    NuMuPmu[i] = new TH1D(Title.c_str(), Title.c_str(), 20, 0.0, 5.0);
    NuMuPmu[i]->SetTitle(Title.c_str());
    NuMuPmu[i]->SetName(Title.c_str());
    NuMuPmu[i]->GetXaxis()->SetTitle("p_{#mu} (GeV/c)");
    NuMuPmu[i]->GetYaxis()->SetTitle("N_{events} (MC*POTw*DETw*FLUXw)");
  }
#endif

  // Important that we've run addSelection by here; so that selections are enabled in psyche
  // Check addSelection has been called by looping over samplepdfs making sure they aren't all NULL
  bool FoundSample = false;
  for (int i = 0; i < SampleId::kNSamples; ++i) {
    if (samplepdfs->At(i) != NULL) {
      FoundSample = true;
      break;
    }
  }
  if (!FoundSample) {
    std::cerr << "I seem to have no samples turned on but you're asking me to initialize psyche systematics. No can do" << std::endl;
    std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
    throw;
  }

  // Number of events in the psyche tree
  unsigned int nMC = _man.GetNMCEvents();

  // Create the array of PreviousToyBox
  _man.sel().CreateToyBoxArray(nMC);
  _man.syst().Initialize(_man.sel(), nMC);

  // Counter of passed events
  int nPassed = 0;
  int nNoPassed = 0;

  // Holds our weights (psyche typedef)
  Weight_h detWeightSyst;
  Weight_h fluxWeightSyst;

  // Diagnostics for counting negative and overflow weights
  int n10Weights = 0;
  int nZeroWeights = 0;
  // Count the number of negative detector weights from input tree
  int NegativeWeights = 0;

  // Now store each class of weights: combine to later give what BANFF provide for cross-validation
  // e.g. data events, raw number of MC events, POT scaled MC evt, POT+det scaled, etc

  // evt counts number of actual events
  UInt_t evt = 0;
  // entry counts number of entries in tree
  Long64_t entry = 0;
  Long64_t nEntries = _man.GetEntries();

  // How often we loop
  int countwidth = double(nMC)/double(10)+1;

  // Find the normalisation bins for current cross-section parameterisation
  FindNormBins();

  double DialNominals[] = {0,0,0,0,0,0,0,0};
  size_t NDials = 8;
  varier_eb =  new NIWGFSLepEbMomentumVariation_continuous(Dials, DialNominals, NDials);


#ifdef DEBUG_TF1
  dumpfile = new TFile((std::string(FitManager->getOutputFilename())+"_tf1_dump.root").c_str(), "recreate");
#endif

  // Reserve memory for cross-section processing
  ReserveMemory(nMC);
#if DEBUG > 0
  EventsPerNorm.reserve(xs2015->getNumParams());
  EventsPerMode.reserve(kMaCh3_nModes);
  for (int i = 0; i < xs2015->getNumParams(); ++i) EventsPerNorm[i] = 0;
  for (int i = 0; i < kMaCh3_nModes; ++i) EventsPerMode[i] = 0;
#endif

  std::cout << "Looping over MC events..." << std::endl;
  // Loop over the events
  while (evt < nMC && entry < nEntries) {

    // Be a little verbose
    if (evt % countwidth == 0) {
      std::cout << "   Event " << evt << "/" << nMC << " (" << int(double(evt)/double((nMC))*100) << "%)" << std::endl;
    }

    // Get the next event in the Experiment class
    AnaSuperEventB* sevent = NULL;
    if (!simple) {
      sevent = _man.GetPreloadedMCSuperEvent(entry);
    } else {
      sevent = _man.LoadSuperEvent(entry);
    }

    // Get the event class
    AnaEventB* event = dynamic_cast<AnaEventB*>(sevent->Event);
    // Fill the psyche EventBox
    _man.sel().InitializeEvent(*event);
    // Initialize the SystBox for variation systematics
    _man.syst().InitializeEventSystematics(_man.sel(), *event);

    // Check if this event passed the processing and return the detector and flux weights
    bool passed = _man.ProcessEvent((*toy), *event, detWeightSyst, fluxWeightSyst);
    if (!passed) {
      _man.syst().FinalizeEventSystematics(*event);
      _man.sel().FinalizeEvent(*event);
      nNoPassed++;
      evt++;
      continue;
    }

    // Get the event weights from psyche
    double POTWeight = sevent->POTWeight;
    double PileUpWeight = event->Weight;
    double DetWeight = detWeightSyst.Correction;
    double FluxWeight = fluxWeightSyst.Correction;

    // Save the event summary pointer
    AnaEventSummaryB* Summary = dynamic_cast<AnaEventSummaryB*>(event->Summary);
    // Get the event sample that this event falls under
    SampleId::SampleEnum SampleID = Summary->EventSample;

    // Check that the sample ID is valid
    if (SampleID == SampleId::kUnassigned) {
      std::cerr << "Error, the sample ID is unassigned for this event through psyche" << std::endl;
      std::cerr << "Indicates something has gone awry with psyche processing!" << std::endl;
      std::cerr << "...or there's a bug in psyche, enjoy!" << std::endl;
      std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
      throw;
    }

    // Check if we have turned off this sample, and check that it passed selection
    if (samplepdfs->At(SampleID) == NULL) {
      // Even when we have failed the sampleID test we have to Finalise the event to make sure the next event is processed OK
      _man.syst().FinalizeEventSystematics(*event);
      _man.sel().FinalizeEvent(*event);
      evt++;
      continue;
    }

    // Check that the truth vertex for this sample is actually filled
    if (Summary->TrueVertex[SampleID] == NULL) {
      std::cerr << "No true vertex assocaited with event " << entry << " and psyche SampleID " << SampleId::ConvertSample(SampleId::SampleEnum(SampleID)) << " (" << SampleID << ")" << std::endl;
      std::cerr << "Skipping this event but this is problematic!" << std::endl;
      std::cerr << "Indicates psyche has some bug in it or you're feeding unreasonable parameters" << std::endl;
      _man.syst().FinalizeEventSystematics(*event);
      _man.sel().FinalizeEvent(*event);
      evt++;
      continue;
    }

    // **************************************************
    // Now we are sure that this event passes our selection
    // No more selection happens after this
    // **************************************************

    // We're on event "evt" in the ROOT file, so get this entry in the TChain
    // We're on "nPassed" event that have passed selection
    // Now get the filled stuff in the TChain (e.g. splines, Q2, species)
    chain->GetEntry(evt);

    // Now first perform checks that we have the same event in the TChain as psyche
    // Check the truth vertex
    int PsycheVertex = Summary->TrueVertex[SampleID]->ID;
    int SavedVertex = truthv_a[0];
    if (PsycheVertex != SavedVertex) {
      std::cerr << "Mismatched truth vertex from psyche and input ROOT file" << std::endl;
      std::cerr << "This _STRONGLY_ indicates you're running with different psycheSteering.parameters.dat in psyche/psycheSteering/*/parameters/ than the input splines were created with" << std::endl;
      std::cerr << "Entry in loop: " << entry << std::endl;
      std::cerr << "Entry in tree: " << evt << std::endl;
      std::cerr << "PsycheVertex: " << PsycheVertex << std::endl;
      std::cerr << "SavedVertex:  " << SavedVertex << std::endl;

      std::cerr << "Sample Number: " << SampleId::ConvertSample(SampleId::SampleEnum(SampleID)) << " (" << SampleID << ")" << std::endl;
      std::cerr << "Neutrino code: " << species_a[0] << std::endl;
      std::cerr << "Reaction code: " << mode_a[0] << std::endl;
      std::cerr << "Enu:           " << enu_a[0] << std::endl;
      std::cerr << "Q2:            " << q2_a[0] << std::endl;

      for (int z = 0; z < SampleId::kNSamples; ++z) {
        if (truthv_a[z] != -999) {
          std::cerr << "Found non-999 saved truth vertex in truthv_a[" << z << "] = " << truthv_a[0] << std::endl;
          std::cerr << "Sample Number: " << SampleId::ConvertSample(SampleId::SampleEnum(z)) << "(" << z << ")" << std::endl;
          std::cerr << "Neutrino code: " << species_a[0] << std::endl;
          std::cerr << "Reaction code: " << mode_a[0] << std::endl;
          std::cerr << "Enu:           " << enu_a[0] << std::endl;
          std::cerr << "Q2:            " << q2_a[0] << std::endl;
        }
      }
      throw;
    }

    // Get the target
    int target = target_a[0];
    // Returns lepton momentum
    double pLep = getKinematicVariable(kinvars[SampleID][0], event);
    // Apply the momentum shift to the lepton
    if (kinvars[SampleID][0] == kLeptonMomentum) {
      pLep += ApplyCoulombShift(event, target, abs(mode_a[0]), kLeptonMomentum);
#if DEBUG == 2
      std::cout << "===================== EVENT " << nPassed << " =====================" << std::endl;
      std::cout << "*** Shifted reconstructed lepton momentum by " << ApplyCoulombShift(event, target, abs(mode_a[0]), kLeptonMomentum) << std::endl;
      std::cout << "    Enu:     " << enu_a[0] << std::endl;
      std::cout << "    Q2:      " << q2_a[0] << std::endl;
      std::cout << "    pLep:    " << pLep << std::endl;
      std::cout << "    Target:  " << target_a[0] << std::endl;
      std::cout << "    Species: " << species_a[0] << std::endl;
      std::cout << "    Mode:    " << mode_a[0] << std::endl;
#endif
    }
    // Get the cos theta of the outgoing lepton
    double CosThLep = getKinematicVariable(kinvars[SampleID][1], event);

    // Push back vitals into the struct
    ndobj->event_sample.push_back(SampleID);
    ndobj->mom.push_back(pLep);
    ndobj->theta.push_back(CosThLep);
    ndobj->mom_adj.push_back(pLep);
    //WP: Eb dial:

    // If 12C or 16O target, get Eb dial
    if((target_a[0]==kTarget_C || target_a[0]==kTarget_O))
    {
      ndobj->eb_dial.push_back(GetDial(species_a[0], target_a[0]));
      ndobj->eb_bin.push_back(varier_eb->GetInputCacheBin(GetDial(species_a[0], target_a[0]),enu_a[0],acos(CosThLep)));
    }else{
      ndobj->eb_dial.push_back(GetDial(12,12));//just a place holder so eb_dial vector is still alligned with number of events
      ndobj->eb_bin.push_back(0);
    }

    // Check the read detector weight from the file is < 10
    if (detweight_tree > 10) {
      n10Weights++;
      std::cerr << "   ***WARNING, detector weight from psyche in input file = " << detweight_tree << " --- Capping at 10!" << std::endl;
      std::cerr << "      Selection:    " << SampleId::ConvertSample(SampleId::SampleEnum(SampleID)) << std::endl;
      std::cerr << "      Momentum:     " << pLep << std::endl;
      std::cerr << "      CosTheta:     " << CosThLep << std::endl;
      std::cerr << "      POTWeight:    " << POTWeight << std::endl;
      std::cerr << "      PileUpWeight: " << PileUpWeight << std::endl;
      std::cerr << "      DetWeight:    " << DetWeight << std::endl;
      std::cerr << "      DetWeight_a:  " << detweight_tree << std::endl;
      std::cerr << "      FluxWeight:   " << FluxWeight << std::endl;
      std::cerr << "      FluxWeight_a: " << fluxweight_a[0] << std::endl;
      detweight_tree = 10;
    }

    // Check the read detector weight from the file is > 0
    if (detweight_tree < 0) {
      nZeroWeights++;
      std::cerr << "   ***WARNING, detector weight from psyche in input file = " << detweight_tree << " --- Setting to 0!" << std::endl;
      std::cerr << "      Selection:    " << SampleId::ConvertSample(SampleId::SampleEnum(SampleID)) << std::endl;
      std::cerr << "      Momentum:     " << pLep << std::endl;
      std::cerr << "      CosTheta:     " << CosThLep << std::endl;
      std::cerr << "      POTWeight:    " << POTWeight << std::endl;
      std::cerr << "      PileUpWeight: " << PileUpWeight << std::endl;
      std::cerr << "      DetWeight:    " << DetWeight << std::endl;
      std::cerr << "      DetWeight_a:  " << detweight_tree << std::endl;
      std::cerr << "      FluxWeight:   " << FluxWeight << std::endl;
      std::cerr << "      FluxWeight_a: " << fluxweight_a[0] << std::endl;
      detweight_tree = 0;
    }

    // The flux weight from flux group
    // flux_w is the total weight to apply to an event from psyche and flux group
    // Here we multiply in the central value from the flux group, provided in the input TChain
    // And the central value detector weight
    ndobj->flux_w.push_back(POTWeight*PileUpWeight*fluxweight_a[0]*detweight_tree);

    // Need the POT weight for each event for the diagnostics TH1Ds
#if DEBUG > 0
    POTWeights.push_back(POTWeight);
#endif

    // Check for negative weights
    if (ndobj->flux_w.back() < 0.0) {
      NegativeWeights++;
      std::cerr << "***Found negative total nominal weight to apply to an event!" << std::endl;
      std::cerr << std::left << std::setw(5) << " " << std::setw(20) << "Event no:" << nPassed << std::endl;
      std::cerr << std::left << std::setw(5) << " " << std::setw(20) << "Flux weight (input from flattree):  " << fluxweight_a[0] << std::endl;
      std::cerr << std::left << std::setw(5) << " " << std::setw(20) << "ND280 weight (input from flattree): " << detweight_tree << std::endl;
      std::cerr << std::left << std::setw(5) << " " << std::setw(20) << "Overall weight: " << ndobj->flux_w[nPassed] << std::endl;
      ndobj->flux_w.back() = 0.0;
    }

    // Push back the q0 q3
    //ndobj->q0.push_back(getKinematicVariablkq0, event));
    //ndobj->q3.push_back(getKinematicVariable(kq3, event));

    // Find the global histogram bin number for this event
    // Do this last in case we have to momentum shift the muon
    // Essentially here so we avoid doing a Fill for every reconfigure
    // 1D
    if (ndims[SampleID] == 1) {
      ndobj->hist_bin.push_back(((TH1*)samplepdfs->At(SampleID))->FindFixBin(pLep));
      // 2D
    } else if (ndims[SampleID] == 2) {
      ndobj->hist_bin.push_back((((TH2Poly*)samplepdfs->At(SampleID))->FindBin(pLep, CosThLep)));
    } else {
      std::cerr << "ndims doesn't appear to be set yet!" << std::endl;
      std::cerr << "EventSample = " << SampleID << " ndims = " << ndims[SampleID] << std::endl;
      throw;
    }

    // Fill the debugging trees
#if DEBUG > 0
    //This is a bit mad, SetBinContent does different things for the overflow than to the rest of the bins for a TH2Poly
    if(ndobj->hist_bin.back()<0){
      MCOnly[SampleID]->SetBinContent(ndobj->hist_bin.back(), 1.0);
      POTOnly[SampleID]->SetBinContent(ndobj->hist_bin.back(), POTWeight);
      FluxOnly[SampleID]->SetBinContent(ndobj->hist_bin.back(), POTWeight*fluxweight_a[0]);
      DetOnly[SampleID]->SetBinContent(ndobj->hist_bin.back(), POTWeight*detweight_tree);
    } else { 
      MCOnly[SampleID]->SetBinContent(ndobj->hist_bin.back(), 1.0+MCOnly[SampleID]->GetBinContent(ndobj->hist_bin.back()));
      POTOnly[SampleID]->SetBinContent(ndobj->hist_bin.back(), POTWeight+POTOnly[SampleID]->GetBinContent(ndobj->hist_bin.back()));
      FluxOnly[SampleID]->SetBinContent(ndobj->hist_bin.back(), POTWeight*fluxweight_a[0]+FluxOnly[SampleID]->GetBinContent(ndobj->hist_bin.back()));
      DetOnly[SampleID]->SetBinContent(ndobj->hist_bin.back(), POTWeight*detweight_tree+DetOnly[SampleID]->GetBinContent(ndobj->hist_bin.back()));
    }
    NuMuPmu[SampleID]->Fill(pLep/1000., POTWeight*detweight_tree*fluxweight_a[0]*PileUpWeight);
#endif

    // If we're using flux, find the fluxbin which the event corresponds to
    if (flux) {
      // Get the neutrino energy
      double E = Summary->TrueVertex[SampleID]->NuEnergy;
      // Get the PDG of incoming neutrino
      int nupdg = Summary->TrueVertex[SampleID]->NuPDG;
      // Check if run period is RHC
      bool isRHC = anaUtils::IsRHC(event->EventInfo.Run);
      // Find the bin which this event sits in
      ndobj->flux_e_bin.push_back(flux->findBin(nupdg, kND280, E/1000., isRHC));
    }

    // Find what ND280 detector systematic parameter applies to this event
    // If we're using simple detector (not using psyche for event by event detector systematics), find the global ND280 systematics bin which the event sits in
    // Do this last in case we have to momentum shift the muon
    if (detsimple) {
      // If 1D
      if (ndims[SampleID] == 1) {
        // detector systematics binning is in lepton momentum
        ndobj->det_bin.push_back(detsimple->getBin(SampleID, pLep));
        // If 2D
      } else if (ndims[SampleID] == 2) {
        ndobj->det_bin.push_back(detsimple->getBin(SampleID, pLep, CosThLep));
      } else {
        std::cerr << "ndims doesn't appear to be set yet!" << std::endl;
        std::cerr << "EventSample = " << SampleId::ConvertSample(SampleId::SampleEnum(SampleID)) << " ndims = " << ndims[SampleID] << std::endl;
        throw;
      }
    }

    // Set up the 2015 cross-section parameterisation
    if (xs2015) {
      // Get the NEUT mode
      int NEUTmode = TMath::Abs(mode_a[0]);
      // Convert to a MaCh3 mode
      xsecInfo[nPassed].mode = NEUTMode_ToMaCh3Mode(NEUTmode);

      // Now set the xsecgraphs for the struct
      for (int k = 0; k < nXsecSplines; k++) {
        xsecgraph[k] = (TGraph*)(grapharrays[k]->At(0)->Clone());
      }

      // Interaction specifics 
      // Need to recalculate Q2 unfortunately because of the lepton momentum shift!
      xsecInfo[nPassed].Q2      = q2_a[0];

      // Need to shift the Q2 because of the Coulomb shift
      // Because we don't have neutrino direction in psyche (only average), we have to use this Q2 to recalculate the new one
      //NIWG actually say Q2 shouldn't be updated. The shift is meant to characterize the effect of the coulomb field as the lepton leaves the nuclear environment, so doesn't affect Q2
      //      xsecInfo[nPassed].Q2 += ApplyCoulombShift(event, target, NEUTmode, kQ2, xsecInfo[nPassed].Q2);

      xsecInfo[nPassed].Enu     = enu_a[0];
      // Write the target
      xsecInfo[nPassed].target  = target_a[0];
      // Species of the interaction
      xsecInfo[nPassed].species = species_a[0];
      // Now sort out the normalisation parameters
      // __This requires the target, species and mode to be set__
      SetEventNormBin(nPassed);

      // Splines are loaded in here
      // for GPU, this is overloaded
      // setSpline _AFTER_ xsecInfo[i].species
      SetSplines(xsecgraph, nPassed);
      // The reduced SetSplines function cuts the TGraph to have half the points
      // Is only left here for documentation purpose pretty much!
      //SetSplines_Reduced(xsecgraph, nPassed);

#if DEBUG > 0
      EventsPerMode[xsecInfo[nPassed].mode]++;
#endif
      // Set the total event weight from stuff like RFG, RPA and coherent
      // This would be parameters that we have reweighted the MC instead of creating an entirely new production

      //For kNIWG2019 we don't need either of these. SF is base model and the cc coherent model this tunes to is the nominal in neut 5.4.0
      xsecInfo[nPassed].weight = 1;//cccohweight_a[SampleID];//rfgweight_a[SampleID] * cccohweight_a[SampleID];

      // Double the weight for NC1gamma events because NEUT was wrong by 1/2 at the time of production (I think... Kendall knows more about this)
      if (xsecInfo[nPassed].mode == kMaCh3_NC1gam) {
        xsecInfo[nPassed].weight *= 2.0;
      }

      /*
      // Patrick Stowell and Kevin McFarland's MINERvA q0 q3 correction
      // Needs to ndobj->q3 to be saved!!
      if (xsecInfo[i].mode == kMaCh3_2p2h) {
      xsecInfo[i].weight *= MINERvA_q0q3_corr(ndobj->q0[i], ndobj->q3[i], false);
      }

      if (xsecInfo[i].mode == kMaCh3_CCQE) {
      xsecInfo[i].weight *= MINERvA_q0q3_corr(ndobj->q0[i], ndobj->q3[i], true);
      }
      */
    }

#if DEBUG == 2
    std::cout << "*** Applying weights to event in fillReweightingBins" << std::endl;
    std::cout << "    Event number:          " << evt << std::endl;
    std::cout << "    Entry number:          " << entry << std::endl;
    std::cout << "    Passed number:         " << nPassed << std::endl;
    xsecInfo[nPassed].Print();
    ndobj->Print(nPassed);
    std::cout << "    POTWeight (psyche):    " << POTWeight << std::endl;
    std::cout << "    PileUpWeight (psyche): " << PileUpWeight << std::endl;
    std::cout << "    DetWeight (psyche):    " << DetWeight << std::endl;
    std::cout << "    fluxweight (tree):     " << fluxweight_a[0] << std::endl;
    std::cout << "    detweight (tree):      " << detweight_tree << std::endl;
#endif

    // Finalise the event systematics and the event
    _man.syst().FinalizeEventSystematics(*event);
    _man.sel().FinalizeEvent(*event);
    // Increment the counters
    evt++;
    // If we reach it this far we've passed the event
    nPassed++;
  } // end while loop on psyche


  // Set the number of events simply to the size of the vector
  ndobj->nEvents = ndobj->mom.size();
  // Shrink the memory to size
  ShrinkMemory(ndobj->nEvents);

  clock.Stop();
  std::cout << "------------------------------------------------------" << std::endl;
  std::cout << "Finished looping psyche events (took " << clock.RealTime() << "s)" << std::endl;
  std::cout << "There are " << ndobj->nEvents << " which passed cuts" << std::endl;
  std::cout << "Allocated " << ndobj->nEvents*(sizeof(__float__)*3 + sizeof(__int__)*4)/1.E6 << " MB to the detector structs" << std::endl;
  std::cout << "          " << ndobj->nEvents*(sizeof(__int__)*4 + sizeof(__float__)*3)/1.E6 << " MB to the xsec structs (excluding new TSpline3 or TF1)" << std::endl;

  std::cout << "  Found " << nZeroWeights << " events with NEGATIVE detector weight from psyche!" << std::endl;
  std::cout << "  Found " << n10Weights << " events with >10 detector weight from psyche!" << std::endl;

#ifdef DEBUG_TF1
  std::cout << "Bye" << std::endl;
  dumpfile->Close();
  delete dumpfile;
  throw;
#endif

  // Now that we have everything set-up, do a nice summary of the samples which we can use later for openMP
#ifdef MULTITHREAD
  maxBins   = 0;
  nSamples  = 0;

  // Find the maximum of hist_bin and event_sample
  // Don't multi-thread this either!!!
  for (unsigned int i = 0; i < ndobj->nEvents; i++) {
    // find maximum hist bin
    if (ndobj->hist_bin[i] > maxBins) {
      maxBins = ndobj->hist_bin[i];
    }

    // find maximum nSamples
    if (ndobj->event_sample[i] > nSamples) {
      nSamples = ndobj->event_sample[i];
    }
  }

  // Take the overflow into account
  nSamples += 1;
  //For th2polys, the overflow bins are -1 to -9, depending on which direction your overflowing. For multithreading, we first add to a private array and then set this as the bin content at the end of the loop over events where array index = bin number. You can't have an array with negative indices so we'll use the final 9 bins as overflow. The last array element will equate to bin -1, the penultimate to bin -2 and so.
  maxBins += 10;
#endif

  clock.Stop();
  std::cout << "Finished cross-section loading: " << ndobj->nEvents*(sizeof(short int)*4+sizeof(double)*4)/1.E6 << " MB + splines allocated" << std::endl;
  std::cout << "Cross-section loading took " << clock.RealTime() << "s" << std::endl;

  if (NegativeWeights > 0) {
    std::cerr << "**********************" << std::endl;
    std::cerr << "Found some _NEGATIVE_WEIGHTS_ for events!" << std::endl;
    std::cerr << "Total number: " << NegativeWeights << std::endl;
    std::cerr << "**********************" << std::endl;
  }

#if DEBUG > 0
  std::cout << "Events by parameter: " << std::endl;
  for (int i = 0; i < xs2015->getNumParams(); ++i) {
    if (EventsPerNorm[i] == 0) continue;
    std::cout << "   " << std::setw(20) << xs2015->getParameterName(i) << EventsPerNorm[i] << std::endl;
  }
  std::cout << "Events by mode: " << std::endl;
  for (int i = 0; i < kMaCh3_nModes; ++i) {
    std::cout << "   " << std::setw(20) << MaCh3mode_ToString(MaCh3_Mode(i)) << EventsPerMode[i] << std::endl;
  }
#endif

  // Delete the allocated memory
  // Delete the input TChain for the splines
  delete chain;
  // Delete the psyche experiment and toys if we're not using psyche for event by event reweighting
  if (simple) {
    // This delete ruins everything because ROOT's TFile delete seems to be broken, argh
    /*
    // Maybe call Experiment destructor too?
    // AddDataSample allocates with new also and does not seem to delete
    // Delete AnalysisManager also? I think we can delete all the psyche stuff here if we use isSimple
    // Loop through the added MC sampleGroups and delete them
    std::map<std::string, SampleGroup> SampleGroups = exp->GetSampleGroups();
    std::map< std::string, SampleGroup >::iterator sit;
    for (sit = SampleGroups.begin(); sit != SampleGroups.end(); sit++) {
    SampleGroup& sampleGroup = sit->second;
    std::map< std::string, DataSample*>& mcSamples = sampleGroup.GetMCSamples();
    std::map< std::string, DataSample*>::iterator it2;
    // Delete MC
    for (it2 = mcSamples.begin(); it2 != mcSamples.end(); it2++) {
    DataSample* sample = it2->second;
    std::cout << "Deleting psyche MC: " << sample->GetFilePath().c_str() << std::endl;
    delete sample;
    }

    // Delete Data
    std::vector<DataSample*>& dataSamples = sampleGroup.GetDataSamples();
    std::vector<DataSample*>::iterator it3;
    for (it3 = dataSamples.begin(); it3 != dataSamples.end(); it3++) {
    DataSample* sample = (*it3);
    std::cout << "Deleting psyche Data: " << sample->GetFilePath().c_str() << std::endl;
    delete sample;
    }
    }
    */
    delete exp;
    delete toy;
  }

  delete[] xsecgraph;
  delete[] grapharrays;
} // End fillReweightingBins() function


// ***************************************************************************
// Function to reserve memory for ND objects
void samplePDFND2019::ReserveMemory(int nEvents) {
  // ***************************************************************************

  // Number of events which passed the above selection
  ndobj->nEvents      = nEvents;
  // Make nice floats (these are saved throughout execution in the ndobj struct, so are floats to save RAM)
  ndobj->mom          .reserve(nEvents);
  ndobj->theta        .reserve(nEvents);
  ndobj->flux_w       .reserve(nEvents);
  ndobj->hist_bin     .reserve(nEvents);
  ndobj->event_sample .reserve(nEvents);
  //ndobj->q0           .reserve(nEvents);
  //ndobj->q3           .reserve(nEvents);

  if (flux) {
    ndobj->flux_e_bin .reserve(nEvents);
  } 

  if (detsimple) {
    ndobj->det_bin    .reserve(nEvents);
  }

  if (xs2015) {
    xsecInfo = new xsec2018<__SPLINE_TYPE__*>[nEvents];
    // Print size of xsec2018 variables
#if DEBUG > 0
    std::cout << "Size of mode (int) in xsec2018: " << sizeof(xsecInfo[0].mode) << std::endl;
    std::cout << "Size of Q2 (float) in xsec2018: " << sizeof(xsecInfo[0].Q2) << std::endl;
#endif
  }

}

// ************************
// Shrink the structs to size
void samplePDFND2019::ShrinkMemory(int nEvents) {
  // ************************
  std::vector<__float__>(ndobj->mom).swap(ndobj->mom);
  std::vector<__float__>(ndobj->theta).swap(ndobj->theta);
  std::vector<__float__>(ndobj->flux_w).swap(ndobj->flux_w);
  //std::vector<double>(ndobj->q0).swap(ndobj->q0);
  //std::vector<double>(ndobj->q3).swap(ndobj->q3);

  std::vector<__int__>(ndobj->event_sample).swap(ndobj->event_sample);
  std::vector<__int__>(ndobj->flux_e_bin).swap(ndobj->flux_e_bin);
  std::vector<__int__>(ndobj->det_bin).swap(ndobj->det_bin);
  std::vector<__int__>(ndobj->hist_bin).swap(ndobj->hist_bin);
}

// ***************************************************************************
// Helper function to find normalisation parameter bins
// This is to avoid hard-coding the normalisations and instead use the input to determine where each normalisation parameter lives in the global parameter index
// e.g. 2p2h normalisation for C/O lives in position x (which we could just hard-code); I prefer to use a search instead in case parameterisations change. It's really up to the end-user!
void samplePDFND2019::FindNormBins() {
  // ***************************************************************************

  // Before we loop over the events we sometimes have dials that are split by carbon or oxygen (e.g. MEC C, MEC O)
  // Let's find what those parameters are in our xsec parameterisation
  // This is sort of silly because we also have nxsec_norm_modes and the other private members right here in the class...
  if (xsecmodel == kNIWG_2020a) {

    // Initalise all the positions to -999
    mec_norm_nu_pos               = -999;
    mec_norm_anu_pos              = -999;

    mec_norm_CtoO_pos             = -999;

    nue_norm_nu_pos               = -999;
    nue_norm_anu_pos              = -999;

    coh_c_pos                     = -999;
    coh_o_pos                     = -999;

    Q2_norm_0_pos                 = -999;
    Q2_norm_1_pos                 = -999;
    Q2_norm_2_pos                 = -999;
    Q2_norm_3_pos                 = -999;
    Q2_norm_4_pos                 = -999;

    Eb_C_nu_pos                   = -999;
    Eb_C_nubar_pos                = -999;
    Eb_O_nu_pos                   = -999;
    Eb_O_nubar_pos                = -999;

    CC_misc_norm_pos              = -999;

    CC_dis_multipi_norm_nu_pos    = -999;
    CC_dis_multipi_norm_nubar_pos = -999;

    CC_nu_pos                     = -999;
    CC_nubar_pos                  = -999;

    // Loop over the Q2 and Eb parameters to find what their indices are
    int berpa_it = 0; // still called berpa_it but now loops over Eb+Q2 instead
    for (std::vector<int>::iterator it = funcParsIndex.begin(); it != funcParsIndex.end(); ++it, ++berpa_it) {

      std::string name = funcParsNames.at(berpa_it);

      if (name == "Q2_norm_0") {
        Q2_norm_0_pos = *it;
      } else if (name == "Q2_norm_1") {
        Q2_norm_1_pos = *it;
      } else if (name == "Q2_norm_2") {
        Q2_norm_2_pos = *it;
      } else if (name == "Q2_norm_3") {
        Q2_norm_3_pos = *it;
      } else if (name == "Q2_norm_4") {
        Q2_norm_4_pos = *it;
      } else if (name == "EB_dial_C_nu") {
        Eb_C_nu_pos = *it;
      } else if (name == "EB_dial_C_nubar") {
	Eb_C_nubar_pos = *it;
      } else if (name == "EB_dial_O_nu") {
        Eb_O_nu_pos = *it;
      } else if (name == "EB_dial_O_nubar") {
        Eb_O_nubar_pos = *it;
      } else {
        std::cerr << "Found a functional parameter which wasn't Q2 or Eb in samplePDFND2019" << std::endl;
        throw;
      }
    }

    // Make a normal iterator here too (used for indexing)
    int stdIt = 0;
    // Loop over the cross-section mode normalisation parameters
    // These are loaded up in setXsecCov
    for (std::vector<XsecNorms2>::iterator it = xsec_norms.begin(); it != xsec_norms.end(); ++it, ++stdIt) {

      // Get the mode this applies to
      int mode = normParsModes.at(stdIt);
      // And the index
      int pos = normParsIndex.at(stdIt);

      // What target(s) does this parameter apply to 
      std::vector<int> targetVec = normParsTarget.at(stdIt);
      // What PDG(s) does this parameter apply to
      std::vector<int> pdgVec = normParsPDG.at(stdIt);

      // Some parameters will have more than one target, those are not what we are looking for here
      // target = 12 or 16, or -1 (both)
      int target = -999;
      // If size is one, just read the target as the first entry
      if (targetVec.size() == 1) {
        target = targetVec.at(0);
        // The other possibility is that this parameter should apply to multiple targets
      } else {
        for (std::vector<int>::iterator jt = targetVec.begin(); jt != targetVec.end(); ++jt) {
          if (*jt != -999) {
            target = 0;
          } else {
            target = *jt;
          }
        } // end for targetVec loop
      }

      // Now we have target == 0 for both C and O, target == 12 for only C, target == 16 for only O

      // If pdg vector is only 1 long it's likely to be nue and nuebar
      int pdg = -999;
      if (pdgVec.size() == 1) {
        pdg = pdgVec.at(0);
      } else {

        // -1 = anti-neutrino, 0 = both, 1 = neutrino
        // Loop over and find if this is a neutrino only
        for (std::vector<int>::iterator jt = pdgVec.begin(); jt != pdgVec.end(); ++jt) {
          // We're only dealing with neutrino
          if (*jt > 0 && (pdg == -999 || pdg == 1)) {
            pdg = 1;
          } else if (*jt < 0 && (pdg == -999 || pdg == -1)) {
            pdg = -1;
          } else {
            pdg = 0;
          }
        }
      }

      //WP: 2018 normalisations for CC events
      if (mode == -3) {
        if (pdg == 1) {
          CC_nu_pos = pos;
        } else if (pdg == -1) {
          CC_nubar_pos = pos;
        }
      }	

      // Now we have anti-netrino = -1, both = 0, neutrino = 1

      // Do the switch on mode
      // We look for 2p2h, CC coh, nue
      else if (mode == kMaCh3_2p2h) {
        // 2p2h C and O (target == 0 means both Carbon and Oxygen)
        if (target == 0) {
          // 2p2h C numu (pdg == 1 means nu)
          if (pdg == 1) {
            mec_norm_nu_pos = pos;
            // 2p2h C numubar
          } else if (pdg == -1) {
            mec_norm_anu_pos = pos;
          }
          // 2p2h C to O (only oxygen target!)
        } else if (target == kTarget_O) {
          mec_norm_CtoO_pos = pos;
          // Something bad happened (2p2h but not oxygen or carbon)
        } else {
          std::cerr << "Running a 2p2h interaction but it isn't written as being on carbon or oxygen, neutrino" << std::endl;
          std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
          throw;
        }

        // Look for CC coherent
      } else if (mode == kMaCh3_CCcoh) {
        if (pdg != 0) {
          std::cerr << "Something went wrong: coherent parameters should be for both neutrino and anti-neutrino" << std::endl;
          throw;
        }
        // Target is carbon
        if (target == kTarget_C) {
          coh_c_pos = pos;
          // Target is oxygen
        } else if (target == kTarget_O) {
          coh_o_pos = pos;
          // Something bad happened
        } else {
          std::cerr << "Running a coherent interaction but it isn't written as being on carbon or oxygen" << std::endl;
          std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
          throw;
        }

        // The nue and nuebar
      } else if (mode == -1) {

        if (target != 0) {
          std::cerr << "Something went wrong: nue_numu and nuebar_numubar shouldn't care about the target, but I've found something which looks like it does!" << std::endl;
          throw;
        }

        // Nue scaling
        if (pdg == kNue) {
          nue_norm_nu_pos = pos;
          // Nue bar scaling
        } else if (pdg == kNue_bar) {
          nue_norm_anu_pos = pos;
          // Something bad happened
        } else {
          std::cerr << "Running a coherent interaction but it isn't written as being on carbon or oxygen" << std::endl;
          std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
          throw;
        }
      } else if(mode == kMaCh3_CCMisc) {
        CC_misc_norm_pos = pos;
      }

      else if(mode == -4) {
        if (pdg == 1) {
          CC_dis_multipi_norm_nu_pos = pos;
        }
        else if (pdg == -1) {
          CC_dis_multipi_norm_nubar_pos= pos;
        }
      } // End if on mode
    } // End the for loop

    // Check that all necessary params have been set
    if (mec_norm_nu_pos == -999) {
      std::cerr << "mec norm nu pos == -999" << std::endl;
      throw;

    } else if (mec_norm_anu_pos == -999) {
      std::cerr << "mec norm anu pos == -999" << std::endl;
      throw;

    } else if (mec_norm_CtoO_pos == -999) {
      std::cerr << "mec norm CtoO pos == -999" << std::endl;
      throw;

    } else if (coh_c_pos == -999) {
      std::cerr << "coh c pos == -999" << std::endl;
      throw;

    } else if (coh_o_pos == -999) {
      std::cerr << "coh c pos == -999" << std::endl;
      throw;

    } else if (nue_norm_nu_pos == -999) {
      std::cerr << "nue norm nu pos == -999" << std::endl;
      throw;

    } else if (nue_norm_anu_pos == -999) {
      std::cerr << "nue norm anu pos == -999" << std::endl;
      throw;

    } else if (Q2_norm_0_pos == -999) {
      std::cerr << "Q2 norm 0 pos == -999" << std::endl;
      throw;
    } else if (Q2_norm_1_pos == -999) {
      std::cerr << "Q2 norm 1 pos == -999" << std::endl;
      throw;
    } else if (Q2_norm_2_pos == -999) {
      std::cerr << "Q2 norm20 pos == -999" << std::endl;
      throw;
    } else if (Q2_norm_3_pos == -999) {
      std::cerr << "Q2 norm 3 pos == -999" << std::endl;
      throw;
    } else if (Q2_norm_4_pos == -999) {
      std::cerr << "Q2 norm 4 pos == -999" << std::endl;
      throw;
    } else if (CC_nu_pos == -999) {
      std::cerr << "CC nu pos == -999" << std::endl;
      throw;

    } else if (CC_nubar_pos == -999) {
      std::cerr << "CC nubar pos == -999" << std::endl;
      throw;

    } else if (Eb_C_nu_pos == -999) {
      std::cerr << "Eb C nu pos == -999" << std::endl;
      throw;

    } else if (Eb_C_nubar_pos == -999) {
      std::cerr << "Eb C nubar pos == -999" << std::endl;
      throw;

    } else if (Eb_O_nu_pos == -999) {
      std::cerr << "Eb O nu pos == -999" << std::endl;
      throw;

    } else if (Eb_O_nubar_pos == -999) {
      std::cerr << "Eb O nubar pos == -999" << std::endl;
      throw;

    } else if (CC_misc_norm_pos == -999) {
      std::cerr << "CC Misc pos == -999" << std::endl;
      throw;

    } else if (CC_dis_multipi_norm_nu_pos == -999) {
      std::cerr << "CC DIS & Multi Pi Norm Nu pos == -999" << std::endl;
      throw;

    } else if (CC_dis_multipi_norm_nubar_pos == -999) {
      std::cerr << "CC DIS & Multi Pi Norm Nubar pos == -999" << std::endl;
      throw;
    }

  } else {
    std::cerr << "Sorry I currently only support NIWG 2019a" << std::endl;
    std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
    throw;
  } // End if on xsecmodel

#ifdef DEBUG
  std::cout << "*********************************************" << std::endl;
  std::cout << "Found all the normalisation parameters in samplePDFND2019 as: " << std::endl;
  std::cout << std::left << std::setw(35) << "2p2h norm C and O neutrino" << mec_norm_nu_pos << std::endl;
  std::cout << std::left << std::setw(35) << "2p2h norm C and O anti-neutrino" << mec_norm_anu_pos << std::endl;
  std::cout << std::left << std::setw(35) << "2p2h norm CtoO neutrino" << mec_norm_CtoO_pos << std::endl;
  std::cout << std::left << std::setw(35) << "CC Coherent on C" << coh_c_pos << std::endl;
  std::cout << std::left << std::setw(35) << "CC Coherent on O" << coh_o_pos << std::endl;
  std::cout << std::left << std::setw(35) << "nue_numu norm" << nue_norm_nu_pos << std::endl;
  std::cout << std::left << std::setw(35) << "nuebar_numubar norm" << nue_norm_anu_pos << std::endl;
  std::cout << std::left << std::setw(35) << "Q2 norm 0" << Q2_norm_0_pos << std::endl;
  std::cout << std::left << std::setw(35) << "Q2 norm 1" << Q2_norm_1_pos << std::endl;
  std::cout << std::left << std::setw(35) << "Q2 norm 2" << Q2_norm_2_pos << std::endl;
  std::cout << std::left << std::setw(35) << "Q2 norm 3" << Q2_norm_3_pos << std::endl;
  std::cout << std::left << std::setw(35) << "Q2 norm 4" << Q2_norm_4_pos << std::endl;
  std::cout << std::left << std::setw(35) << "Eb C nu" << Eb_C_nu_pos << std::endl;
  std::cout << std::left << std::setw(35) << "Eb C nubar" << Eb_C_nubar_pos << std::endl;
  std::cout << std::left << std::setw(35) << "Eb O nu" << Eb_O_nu_pos << std::endl;
  std::cout << std::left << std::setw(35) << "Eb O nubar" << Eb_C_nubar_pos << std::endl;
  std::cout << std::left << std::setw(35) << "CC nu" << CC_nu_pos << std::endl;
  std::cout << std::left << std::setw(35) << "CC nubar" << CC_nubar_pos << std::endl;
  std::cout << std::left << std::setw(35) << "CC misc" << CC_misc_norm_pos << std::endl;
  std::cout << std::left << std::setw(35) << "CC dis & multi pi nu" << CC_dis_multipi_norm_nu_pos << std::endl;
  std::cout << std::left << std::setw(35) << "CC dis & multi pi nubar" << CC_dis_multipi_norm_nubar_pos << std::endl;
  std::cout << "*********************************************" << std::endl;
#endif
}

// ***********************************
// Set the normbin for each event, needs energy for 2018 CC normalisations
void samplePDFND2019::SetEventNormBin(const int EventNumber) {
  // ***********************************

  // If we're doing things the smart way
  if (xsecmodel != kNIWG_2020a) {
    std::cerr << "Only support NIWG 2020a" << std::endl;
    std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
  }

  xsecInfo[EventNumber].normbin = -1;

  // Loop the normalisation parameters (read in on setXsecCov)
  for (std::vector<XsecNorms2>::iterator it = xsec_norms.begin(); it != xsec_norms.end(); ++it) {

    // Look at the mode for this event and match it to the parameter
    if ((*it).mode != xsecInfo[EventNumber].mode) {
      continue;
    }

    // Loop over target material and match
    if ((*it).targets.size() > 0) {
      bool TargetMatch = false;
      for (std::vector<int>::iterator jt = (*it).targets.begin(); jt != (*it).targets.end(); ++jt) {
        if ((*jt) != xsecInfo[EventNumber].target) {
          continue;
        } else {
          TargetMatch = true;
          break;
        }
      }
      if (!TargetMatch) continue;
    }

    // Loop over neutrino pdg and match
    if ((*it).pdgs.size() > 0) {
      bool PDGMatch = false;
      for (std::vector<int>::iterator jt = (*it).pdgs.begin(); jt != (*it).pdgs.end(); ++jt) {
        if ((*jt) != xsecInfo[EventNumber].species) {
          continue;
        } else {
          PDGMatch = true;
          break;
        }
      }
      if (!PDGMatch) continue;
    }

    // If we have separate dials for carbon and oxygen, we need to set these here
    // e.g. if the event mode is 2p2h and the target is oxygen, set the normalisation bin to that equivalent
    if (xsecInfo[EventNumber].mode == kMaCh3_2p2h) {
      // If 2p2h on Carbon
      if (xsecInfo[EventNumber].target == kTarget_C) {
        // Neutrino interaction
        if (xsecInfo[EventNumber].species > 0) {
          xsecInfo[EventNumber].normbin = mec_norm_nu_pos;
          // Anti-neutrino interaction
        } else if (xsecInfo[EventNumber].species < 0) {
          xsecInfo[EventNumber].normbin = mec_norm_anu_pos;
        } // Also have electron neutrinos here so don't have the else
        // If 2p2h on Oxygen
      } else if (xsecInfo[EventNumber].target == kTarget_O) {
        xsecInfo[EventNumber].normbin = mec_norm_CtoO_pos;
      } // Interactions happen on other elements too, e.g. 28, so don't include an else here
      // Do the Coherent modes, which are separated by carbon and oxygen in 2016a
    } else if (xsecInfo[EventNumber].mode == kMaCh3_CCcoh) {
      // The Carbon target
      if (xsecInfo[EventNumber].target == kTarget_C) {
        xsecInfo[EventNumber].normbin = coh_c_pos;
        // The Oxygen target
      } else if (xsecInfo[EventNumber].target == kTarget_O) {
        xsecInfo[EventNumber].normbin = coh_o_pos;
      } // Interactions happen on other elements too, e.g. 50, so don't include an else here

      // If none of the above special cases, just set the normbin to the bin we've earlier read in
      // It essentially just means that we have a normalisation parameter which applies to every type of interaction target and species
    } else {
      int detOn = xsecBitField[(*it).startbin];
      if (int(detOn&1) == 0) {
        continue;
      }
      xsecInfo[EventNumber].normbin = (*it).startbin;
    }
  }
}

// ******************************************
// Apply the Coulomb shift in 2018 analyses
// -3.6 MeV for negative leptons on C
// -4.3 MeV for negative leptons on O
// +2.6 MeV for positive leptons on C
// +3.3 MeV for positive leptons on O
// Sadly we have to pass the target because psyche doesn't seem to save it
// Warning, ONLY APPLY TO CC INTERACTIONS
// InitialValue is only used when kKinType == kQ2 since we can't calculate Q2 directly from psyche we need to find the incoming neutrino angle. Defaults to 0.0
double samplePDFND2019::ApplyCoulombShift(AnaEventB *event, int target, int reaction, ND280KinematicTypes kKinType, double InitialValue) {
  // ******************************************

  // Only want CC events
  if (reaction >= kNEUT_NC1pi01n) return 0.0;

  // The event summary
  AnaEventSummaryB* Summary = dynamic_cast<AnaEventSummaryB*>(event->Summary);
  // Get the selection this event falls into
  SampleId::SampleEnum index = Summary->EventSample;

  // Make sure the lepton candidate has true particle ID of lepton (so we only shift true leptons that are reconstructed as leptons, not shift true pions that get reconstructed as lepton)
  // Get the lepton candidate
  AnaParticleB* LeptonCandidate = dynamic_cast<AnaParticleB*>(Summary->LeptonCandidate[index]);
  // Get the true particle giving rise to our lepton candidate
  AnaTrueParticleB* TrueParticle = dynamic_cast<AnaTrueParticleB*>(LeptonCandidate->GetTrueParticle());
  // Get the true particle's PDG
  int TruePDG = TrueParticle->PDG;

  // Get the neutrino PDG
  AnaTrueVertexB *vtx = Summary->TrueVertex[index];
  int NuPDG = vtx->NuPDG;

  // The amount we shift by
  double shift = 0.0;

  // Check the lepton candidate matches the neutrino PDG
  // If anti-neutrino interaction, check that NuPDG 
  if (NuPDG < 0) {
    if (NuPDG != TruePDG-1) return shift;
  } else {
    if (NuPDG != TruePDG+1) return shift;
  }
  // Check PID of true particle making this lepton candidate
  // This is now done above as before. The following line is only left for historical comparison
  //if (abs(TruePDG) != 11 && abs(TruePDG) != 13 && abs(TruePDG) != 15) return 0.0;

  // Shift the lepton momentum
  if (kKinType == kLeptonMomentum) {
    // Neutrino
    if (TruePDG > 0) {
      // If target was C
      if (target == kTarget_C) {
        shift = -3.6;
        // If target was O
      } else if (target == kTarget_O) {
        shift = -4.3;
      }
      // Anti-neutrino
    } else if (TruePDG < 0) {
      // If target was C
      if (target == kTarget_C) {
        shift = 2.6;
      } else if (target == kTarget_O) {
        shift = 3.3;
      }
    }
    // Return the shiften lepton momentum
    return shift;
    // Shift the Q2 if requested
  } else if (kKinType == kQ2) {
    // Get the true vertex
    AnaTrueVertexB *vtx = Summary->TrueVertex[index];
    if (vtx == NULL) {
      std::cerr << "Could not find TrueVertex object, sorry" << std::endl;
      std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
      return 0.0;
    }
    double EnuTrue = vtx->NuEnergy/1.E3;

    // Then we get the TLorentzVector of the highest momentum muon
    double ELepTrue = ((MaCh3Utils::GetHMParticle(Summary, TruePDG)).E());
    double PLepTrue = ((MaCh3Utils::GetHMParticle(Summary, TruePDG)).Vect().Mag());
    double MLep = sqrt(ELepTrue*ELepTrue - PLepTrue*PLepTrue);

    // Now we can get the ThNuMu true
    double CosThNuMu = (2.0*EnuTrue*ELepTrue - MLep*MLep - InitialValue)/(2.0*EnuTrue*PLepTrue);

    // Now shift the muon momentum (returned shift is in MeV)
    PLepTrue += (ApplyCoulombShift(event, target, reaction, kLeptonMomentum))/1.E3;
    // And update the lepton momentum
    ELepTrue = sqrt(PLepTrue*PLepTrue + MLep*MLep);

    // Calculate the new Q2
    double Q2Upd = -(MLep*MLep) + 2.0*EnuTrue*(ELepTrue - PLepTrue*CosThNuMu);

    // And finally return the shift
    shift = Q2Upd - InitialValue;

    if (std::isnan(shift)) {
      std::cerr << "Found nan shift" << std::endl;
      std::cerr << "TruePDG:  " << TruePDG << std::endl;
      std::cerr << "PLepTrue: " << PLepTrue << std::endl;
      std::cerr << "ELepTrue: " << ELepTrue << std::endl;
      std::cerr << "MLep:     " << MLep << std::endl;
      std::cerr << "CosThNuMu: " << CosThNuMu << std::endl;
      std::cerr << "EnuTrue:  " << EnuTrue << std::endl;
      std::cout << "Q2Upd:    " << Q2Upd << std::endl;
      throw;
    }
    return shift;
  } else {
    std::cerr << "Can't Coulomb shift the requested variable " << ND280Kinematic_ToString(kKinType) << " (" << kKinType << ")" << std::endl;
    std::cerr << "Will throw error" << std::endl;
    std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
    throw;
  }

  return shift;
}

// **************************************************************************
// Recalculate Q^2 after Eb shift. Takes in shifted lepton momentum, lepton angle, and true neutrino energy
double samplePDFND2019::CalculateQ2(double PLep, double PUpd, double EnuTrue, double InitialQ2){
  // ***************************************************************************

  const double MLep = 0.10565837;

  // Caluclate muon energy
  double ELep = sqrt((MLep*MLep)+(PLep*PLep));

  double CosTh = (2*EnuTrue*ELep - MLep*MLep - InitialQ2)/(2*EnuTrue*PLep);

  ELep = sqrt((MLep*MLep)+(PUpd*PUpd));

  // Calculate the new Q2
  double Q2Upd = -(MLep*MLep) + 2.0*EnuTrue*(ELep - PUpd*CosTh);

  return Q2Upd - InitialQ2;
}


// **************************************************************************
// Recalculate Enu after Eb shift. Takes in shifted lepton momentum, lepton angle, and binding energy change, and if nu/anu
double samplePDFND2019::CalculateEnu(double PLep, double costh, double Eb, bool neutrino){
  // ***************************************************************************

  double mNeff = 0.93956536 - Eb / 1000.;
  double mNoth = 0.93827203;

  if (!neutrino) {
    mNeff = 0.93827203 - Eb / 1000.;
    mNoth = 0.93956536;
  }

  double mLep = 0.10565837;
  double eLep = sqrt(PLep * PLep + mLep * mLep);

  double Enu = (2 * mNeff * eLep - mLep * mLep + mNoth * mNoth - mNeff * mNeff) /(2 * (mNeff - eLep + PLep * costh));

  return Enu;

}

// ***************************************************************************
// Helper function to check if the covariances have been set
void samplePDFND2019::CheckCovariances() {
  // ***************************************************************************

  // Check if we're doing smart processing with a 2015 parameterisation
  // This shouldn't happen because (to my knowledge) there is no such root file yet
  if (xsecmodel < kNIWG_2017a) {
    std::cerr << "Trying to use smart read-in with 2015 NIWG parameterisation, which is probably wrong" << std::endl;
    std::cerr << "So far only >= 2016a has been processed like this" << std::endl;
    throw;
  }

  if (xs2015 == NULL) {
    std::cerr << "xs2015 covariance not set, can't setup samplePDFND2019 in " << __FILE__ << ":" << __LINE__ << std::endl;
    throw;
  }

  if (flux == NULL) {
    std::cerr << "flux covariance not set, can't setup samplePDFND2019 in " << __FILE__ << ":" << __LINE__ << std::endl;
    throw;
  }

  if (detsimple == NULL) {
    std::cerr << "detector covariance not set, can't setup samplePDFND2019 in " << __FILE__ << ":" << __LINE__ << std::endl;
    throw;
  }

}

#if DEBUG > 0
// ***************************************************************************
// Dump the spline draw for a given event to an output file
// This is useful to check if the splines are bad
void samplePDFND2019::DumpSplines(double xsecw, double beamw, double detw, int i) {
  // ***************************************************************************

  std::cerr << "===============\nFOUND STRANGE EVENT IN SAMPLE\n===============" << std::endl;

  std::cerr << "Current xs2015 configuration:" << std::endl;
  xs2015->printNominalCurrProp();

  PrintStructs(xsecw, beamw, detw, i);

  // Make a string describing the run
  std::string OutputName = std::string((*prod).Data());
  while (OutputName.find(" ") != std::string::npos) {
    OutputName.replace(OutputName.find(" "), 1, std::string("_"));
  }
  std::stringstream ss;
  ss << i;
  OutputName += "_event_" + ss.str() + ".root";

#if USE_SPLINE != USE_TF1_red
  // Make a TFile
  TFile *file = new TFile(OutputName.c_str(), "recreate");
  TCanvas *canv = new TCanvas("canv", "canv", 1080, 1080);
  file->cd();
  canv->cd();
  canv->Draw();
  // Loop over the spline parameters and write them to file
  for (int j = 0; j < nSplineParams; ++j) {
    if (xsecInfo[i].GetFunc(j) == NULL) continue;
    xsecInfo[i].GetFunc(j)->Draw("LP*");
    canv->Write(xsecInfo[i].GetFunc(j)->GetTitle());
  }

  file->Write();
  file->Close();
  delete file;
  delete canv;
  std::cerr << "Wrote collected spline response to " << OutputName << std::endl;
#endif
}
#endif

#ifdef DEBUG
// ***************************************************************************
// Prints the relevant struct information and weights for one event
void samplePDFND2019::PrintStructs(double xsecw, double beamw, double detw, const int i) {
  // ***************************************************************************
  std::cout << "===============\nPRINTING EVENT BY EVENT WEIGHTS\n===============" << std::endl;
  std::cout << "EVENT " << i << std::endl;
  xsecInfo[i].Print();
  ndobj->Print(i);
  std::cout << std::setw(20) << "beamw         = " << beamw << std::endl;
  std::cout << std::setw(20) << "detw     = "   << detw << std::endl;
  std::cout << std::setw(20) << "xsecw    = "  << xsecw << std::endl;
  std::cout << std::setw(20) << "totalw (xsec*det*beam*flux)  = " << xsecw*detw*beamw*ndobj->flux_w[i] << std::endl;
}
#endif

// ***************************************************************************
// Silence cout and cerr. Last is risky but psyche persists on spamming both
void samplePDFND2019::QuietPlease() {
  // ***************************************************************************
#if DEBUG > 0
  return;
#else
  buf = std::cout.rdbuf();
  errbuf = std::cerr.rdbuf();
  std::cout.rdbuf( NULL );
  std::cerr.rdbuf( NULL );
#endif
}

// ***************************************************************************
// Reset cout and cerr
void samplePDFND2019::NowTalk() {
  // ***************************************************************************
#if DEBUG > 0
  return;
#else
  std::cout.rdbuf(buf);
  std::cerr.rdbuf(errbuf);
#endif
}


#if USE_SPLINE < USE_TF1
// *************************
// Only need to do the binary search once per parameter, not once per event!
// Takes down the number of binary searches from 1.2M to 17, ha!
// For ROOT version see root/hist/hist/src/TSpline3.cxx TSpline3::FindX(double)
void samplePDFND2019::FindSplineSegment() {
  // *************************

  if (xsecmodel != kNIWG_2020a) {
    std::cerr << "ERROR" << std::endl;
    std::cerr << "FindSplineSegment() only supported in NIWG 2019a" << std::endl;
    std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
    throw;
  }

  // Loop over the splines
  for (int i = 0; i < nSplineParams; i++) {

    const int nPoints = SplineInfoArray[i].nPts;
    const double* xArray = SplineInfoArray[i].xPts;

    if (nPoints == -999 || xArray == NULL) {
      std::cerr << "ERROR" << std::endl;
      std::cerr << "SplineInfoArray[" << i << "] isn't set yet" << std::endl;
      std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
      continue;
      //      throw;
    }

    // Get the variation for this reconfigure for the ith parameter
    const int GlobalIndex = splineParsIndex[i];
    const double xvar = xs2015->calcReWeightFrac(GlobalIndex);

#if DEBUG > 0
    std::cout << "Finding spline segment for xsec parameter " << splineParsNames[i] << " at " << xvar << std::endl;
#endif

    // The segment we're interested in (klow in ROOT code)
    int segment = 0;
    int kHigh = nPoints-1;

    // If the variation is below the lowest saved spline point
    if (xvar <= xArray[0]) {
      segment = 0;
      // If the variation is above the highest saved spline point
    } else if (xvar >= xArray[nPoints-1]) {
      // Yes, the -2 is indeed correct, see TSpline.cxx:814 and //see: https://savannah.cern.ch/bugs/?71651
      segment = kHigh;
      // If the variation is between the maximum and minimum, perform a binary search
    } else {
      // The top point we've got
      int kHalf = 0;
      // While there is still a difference in the points (we haven't yet found the segment)
      // This is a binary search, incrementing segment and decrementing kHalf until we've found the segment
      while (kHigh - segment > 1) {
        // Increment the half-step 
        kHalf = (segment + kHigh)/2;
        // If our variation is above the kHalf, set the segment to kHalf
        if (xvar > xArray[kHalf]) {
          segment = kHalf;
          // Else move kHigh down
        } else {
          kHigh = kHalf;
        }
      } // End the while: we've now done our binary search
    } // End the else: we've now found our point

    if (segment >= nPoints-1 && nPoints > 1) segment = nPoints-2;

    // Save the segment for the ith parameter
    SplineInfoArray[i].CurrSegment = segment;

#ifdef DEBUG
    if (SplineInfoArray[i].xPts[segment] > xvar && segment != 0) {
      std::cerr << "Found a segment which is _ABOVE_ the variation!" << std::endl;
      std::cerr << "IT SHOULD ALWAYS BE BELOW! (except when segment 0)" << std::endl;
      std::cerr << splineParsNames[i] << std::endl;

      std::cerr << "Found segment   = " << segment << std::endl;
      std::cerr << "Doing variation = " << xvar << std::endl;
      std::cerr << "x in spline     = " << SplineInfoArray[i].xPts[segment] << std::endl;
      for (int j = 0; j < SplineInfoArray[j].nPts; ++j) {
        std::cerr << "    " << j << " = " << SplineInfoArray[i].xPts[j] << std::endl;
      }
      std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
      throw;
    }
#endif
  }
}
#endif

// *************************
// Prepare the weights before the reweight loop starts
void samplePDFND2019::PrepareWeights() {
  // *************************
  // Find the relevant spline segments for these parameter variations
#if USE_SPLINE < USE_TF1
  FindSplineSegment();
#endif
}

// *************************
void samplePDFND2019::GetKinVars(ND280KinematicTypes &TypeX, ND280KinematicTypes &TypeY) {
  // *************************
  ND280KinematicTypes Type1 = kLeptonMomentum;
  ND280KinematicTypes Type2 = kLeptonCosTheta;

  for (int i = 0; i < SampleId::kNSamples; ++i) {
    if (samplepdfs->At(i) == NULL) continue;
    Type1 = kinvars[i][0];
    Type2 = kinvars[i][1];
  }
  TypeX = Type1;
  TypeY = Type2;
}

// ************************
// Print POT from psyche
// Psyche also has similar functions in psycheCore/*/src/Experiment.hxx/cxx
void samplePDFND2019::PrintPOT() {
  // ************************

  if (exp == NULL) {
    std::cerr << "PsycheCore Experiment instance is NULL and you're trying to count the POT" << std::endl;
    std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
    throw;
  }

  // Count the POT
  float TotalPOTData = 0.0;
  float TotalPOTMC = 0.0;
  float TotalPOTSand = 0.0;
  std::cout << "================" << std::endl;
  // Sample groups
  std::map<std::string, SampleGroup> SampleGroups = exp->GetSampleGroups();
  for (std::map<std::string, SampleGroup>::iterator sit = SampleGroups.begin(); sit != SampleGroups.end(); ++sit) {
    // Get the current sample group
    SampleGroup TempGroup = (*sit).second;

    // Get the data
    std::vector<DataSample*> DataSamples = TempGroup.GetDataSamples();
    std::cout << "-------------------------" << std::endl;
    std::cout << "Sample group " << sit->first << std::endl;
    std::cout << "Data:" << std::endl;
    for (std::vector<DataSample*>::iterator it = DataSamples.begin(); it != DataSamples.end(); ++it) {
      std::cout << (*it)->GetFilePath() << " POT:" << std::endl;
      (*it)->DumpPOT();
    }

    // Get the MC
    std::map<std::string, DataSample*> MCSampleGroups = TempGroup.GetMCSamples();
    std::cout << "MC:" << std::endl;
    for (std::map<std::string, DataSample*>::iterator it = MCSampleGroups.begin(); it != MCSampleGroups.end(); ++it) {
      std::cout << (*it).first << " from " << (*it).second->GetFilePath() << " POT:" << std::endl;
      (*it).second->DumpPOT();
    }
  }

  std::cout << "Total POT Data: " << TotalPOTData << std::endl;
  std::cout << "Total POT MC  : " << TotalPOTMC << std::endl;
  std::cout << "Total POT Sand: " << TotalPOTSand << std::endl;
  std::cout << "================" << std::endl;
}

// **************************************************
// Helper function to reset the data and MC histograms
void samplePDFND2019::ResetHistograms() {
  // **************************************************

  // Loop over the samples and reset the and fill with zeros
  // Don't openMP this; no signficiant gain
  for (int i = 0; i < SampleId::kNSamples; i++) {	 
    if (samplepdfs->At(i) == NULL) continue;
    // Reset the number of samples that are turned on

    if (ndims[i] == 1) {
      ((TH1*)samplepdfs->At(i))->Reset();
      ((TH1*)samplepdfs->At(i))->Fill(0.0, 0.0);
    } else if (ndims[i] == 2) {
      ((TH2Poly*)samplepdfs->At(i))->Reset("");
      ((TH2Poly*)samplepdfs->At(i))->Fill(0.0, 0.0, 0.0);
      // Don't reset the MC error histogram
      //W2Hist[i]->Reset("");
      //W2Hist[i]->Fill(0.0, 0.0, 0.0)
    }
    // If we want to plot according to mode
    if (modepdf) {
      for (int j = 0; j < kMaCh3_nModes+1; j++) {
        // If 1D
        if (ndims[i] == 1) {
          ((TH1*)((TObjArray*)samplemodepdfs->At(i))->At(j))->Reset();
          ((TH1*)((TObjArray*)samplemodepdfs->At(i))->At(j))->Fill(0.0, 0.0);
          // If 2D
        } else if (ndims[i] == 2) {
          ((TH2Poly*)((TObjArray*)samplemodepdfs->At(i))->At(j))->Reset("");
          ((TH2Poly*)((TObjArray*)samplemodepdfs->At(i))->At(j))->Fill(0.0, 0.0, 0.0);
        }
      } // end the mode loop
    } // end the if modepdf
#if DEBUG > 0
    // Only reset the histograms that get refilled on event-by-event basis
    // i.e. don't need to reset POT scaled because they're fixed once we've looped once
    XsecOnly[i]->Reset("");
    NDCovOnly[i]->Reset("");
    BeamCovOnly[i]->Reset("");
    AllOnly[i]->Reset("");
#endif
  } // end loop over samples
} // end function
