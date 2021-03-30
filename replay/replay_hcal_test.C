R__ADD_INCLUDE_PATH($SBS/include)
R__ADD_LIBRARY_PATH($SBS/lib64)
R__ADD_LIBRARY_PATH($SBS/lib)
R__LOAD_LIBRARY(libsbs.so)

#if !defined(__CLING__) || defined(__ROOTCLING__)
#include <iostream>

#include "TSystem.h"
#include "TList.h"
#include "THaRun.h"
#include "THaEvent.h"
#include "THaAnalyzer.h"
#include "THaApparatus.h"

//#include "SBSGEMStand.h"
//#include "SBSBigBite.h"
#include "SBSEArm.h"
#include "SBSHCal.h"
#endif

void replay_hcal_test(Int_t runnum = 931, Int_t lastEvent = -1){

  gSystem->Load("../libsbs.so");
  SBSHCal *hcal = new SBSHCal("hcal","HCAL");

  SBSEArm *harm = new SBSEArm("sbs","Hadron Arm with HCal");
  harm->AddDetector(hcal);

  //
  //  Steering script for Hall A analyzer demo
  //

  // Set up the equipment to be analyzed.

  // add the two spectrometers with the "standard" configuration
  // (VDC planes, S1, and S2)
  // Collect information about a easily modified random set of channels
  // (see DB_DIR/*/db_D.dat)
  /*
     THaApparatus* DECDAT = new THaDecData("D","Misc. Decoder Data");
     gHaApps->Add( DECDAT );
     */


  // Set up the analyzer - we use the standard one,
  // but this could be an experiment-specific one as well.
  // The Analyzer controls the reading of the data, executes
  // tests/cuts, loops over Apparatus's and PhysicsModules,
  // and executes the output routines.
  THaAnalyzer* analyzer = new THaAnalyzer;

  gHaApps->Add(harm);

  // A simple event class to be output to the resulting tree.
  // Creating your own descendant of THaEvent is one way of
  // defining and controlling the output.
  THaEvent* event = new THaEvent;

  // Define the run(s) that we want to analyze.
  // We just set up one, but this could be many.
  //  THaRun* run = new THaRun( "prod12_4100V_TrigRate25_4.dat" );
  THaRun* run = new THaRun(TString::Format("data/fadc_%d.dat.0",runnum) );
  run->SetLastEvent(lastEvent);

  run->SetDataRequired(0);
  run->SetDate(TDatime());

  analyzer->SetVerbosity(0);
  analyzer->SetOdefFile("output_hcal_test.def");

  // Define the analysis parameters
  analyzer->SetEvent( event );
  analyzer->SetOutFile( TString::Format("rootfiles/fadc_%d.root",runnum));
  // File to record cuts accounting information
  analyzer->SetSummaryFile("sbs_hcal_test.log"); // optional

  //analyzer->SetCompressionLevel(0); // turn off compression
  analyzer->Process(run);     // start the actual analysis
}
