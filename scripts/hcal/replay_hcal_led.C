#include "TSystem.h"
#include "TList.h"
#include "THaRun.h"
#include "THaEvent.h"
#include "THaAnalyzer.h"
#include "THaApparatus.h"

//#include "SBSGEMStand.h"
//#include "SBSBigBite.h"
R__ADD_INCLUDE_PATH(../../)
R__LOAD_LIBRARY(../../build/libsbs)
#include "SBSEArm.h"
#include "SBSHCal.h"

void replay_hcal_led(Int_t runnum = 288, Int_t lastEvent = -1){

  gSystem->Load("../libsbs.so");
  SBSHCal *hcal = new SBSHCal("hcal","HCAL");
  hcal->SetWithTDC(false);
  //SBSCalorimeter *hcal = new SBSCalorimeter("hcal","HCAL");
  hcal->SetWithADCSamples(true);

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
  //THaRun* run = new THaRun(TString::Format("/home/daq/data/fadcAERO_%d.dat.0",runnum) );
  THaRun* run = new THaRun(TString::Format("%s/fadcAERO_%d.dat.0",getenv("HCAL_DATA"),runnum));
  run->SetLastEvent(lastEvent);

  run->SetDataRequired(0);
  run->SetDate(TDatime());

  analyzer->SetVerbosity(2);
  analyzer->SetOdefFile("output_hcal_led.def");
  analyzer->SetMarkInterval(500);

  // Define the analysis parameters
  analyzer->SetEvent( event );
  //analyzer->SetOutFile( TString::Format("rootfiles/fadcAERO_%d.root",runnum));
  analyzer->SetOutFile( TString::Format("%s/fadcAERO_%d.root",getenv("HCAL_ROOTFILES"),runnum));
  // File to record cuts accounting information
  analyzer->SetSummaryFile(TString::Format("sbs_hcal_led_%d.log",runnum)); // optional

  //analyzer->SetCompressionLevel(0); // turn off compression
  analyzer->Process(run);     // start the actual analysis
}
