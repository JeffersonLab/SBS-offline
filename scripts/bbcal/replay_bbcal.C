#include <iostream>

#include "TSystem.h"
#include "TString.h"
#include "TFile.h"

#include "THaShower.h"
#include "THaEvData.h"
#include "THaRun.h"
#include "THaAnalyzer.h"
#include "THaVarList.h"

#include "SBSBigBite.h"
#include "SBSBBShower.h"

// Simple example replay script
//
// Ole Hansen, 11 April 2016
void replay_bbcosmics(int run_number = 124, uint nev = -1)
{
  //load SBS-offline
  gSystem->Load("libsbs.so");
  //--- Define the experimental configuration, i.e. spectrometers, detectors ---

  //THaHRS* bb = new THaHRS("R", "Right HRS" );
  //hrs->AddDetector( new THaVDC("vdc", "Vertical Drift Chambers") );
  // gHaApps->Add(hrs);
  
  SBSBigBite* bigbite = new SBSBigBite("bb", "BigBite spectrometer" );
  //bigbite->AddDetector( new THaShower("ps", "BigBite preshower") );
  bigbite->AddDetector( new SBSBBShower("ps", "BigBite preshower") );
  bigbite->AddDetector( new SBSBBShower("sh", "BigBite shower") );
  gHaApps->Add(bigbite);
  
  // Ideal beam (perfect normal incidence and centering)
  // THaIdealBeam* ib = new THaIdealBeam("IB", "Ideal beam");
  // gHaApps->Add(ib);

  //--- Set up the run we want to replay ---

  // This often requires a bit of coding to search directories, test
  // for non-existent files, etc.

  TString experiment = "bbcal";

  // Get the run number
  TString data_dir = gSystem->Getenv("DATA_DIR");
  if( data_dir.IsNull() ) 
    data_dir = ".";
  TString run_file = data_dir + "/test_" + experiment + Form("_%d.dat",run_number);
  if( gSystem->AccessPathName(run_file) ) {
   Error("replay.C", "Input file does not exist: %s", run_file.Data() );
   exit(-1);
  }

  THaRun* run = new THaRun( run_file, "BB cosmics" );
  run->SetDataRequired(7);//for the time being
  
  cout << "Number of events to replay (-1=all)? ";
  if( nev > 0 )
    run->SetLastEvent(nev);
  
  //--- Set up any physics calculations we want to do ---

  // Extract the reconstructed target quantities of the golden track
  // Not really a physics calculation, but a convenience function.
  // It effectively converts the L.tr.* variables, which are arrays, 
  // to scalers L.gold.*

  //gHaPhysics->Add( new THaGoldenTrack( "R.gold", "RHRS golden track", "R" ));

  // Single-arm electron kinematics for the one spectrometer we have set up.
  // We assume a carbon-12 target (12 AMU)
  //gHaPhysics->Add( new THaPrimaryKine( "R.ekine", "RHRS electron kinematics",
  //"R", 0.511e-3, 12*0.9315 ));

  // Vertex position calculated from RHRS golden track and ideal beam
  // (will poor resolution if raster is on)
  //gHaPhysics->Add( new THaReactionPoint( "R.vx", "Vertex R", "R", "IB" ));

  //--- Define what we want the analyzer to do ---
  // The only mandatory items are the output definition and output file names
  
  THaAnalyzer* analyzer = new THaAnalyzer;

  TString out_dir = gSystem->Getenv("OUT_DIR");
  if( out_dir.IsNull() )
    out_dir = ".";
  TString out_file = out_dir + "/" + experiment + Form("_%d.root", run_number);

  analyzer->SetOutFile( out_file );
  
  analyzer->SetCutFile( "replay.cdef" );
  analyzer->SetOdefFile( "replay.odef" );

  analyzer->SetVerbosity(2);  // write cut summary to stdout
  analyzer->EnableBenchmarks();

  //--- Analyze the run ---
  // Here, one could replay more than one run with
  // a succession of Process calls. The results would all go into the
  // same ROOT output file

  run->Print();
  

  analyzer->Process(run);

  // Clean up

  analyzer->Close();
  delete analyzer;
  //gHaCuts->Clear();
  gHaVars->Clear();
  gHaPhysics->Delete();
  gHaApps->Delete();

  // Open the ROOT file so that a user doing interactive analysis can 
  // immediately look at the results. Of course, don't do this in batch jobs!
  // To close the file later, simply type "delete rootfile" or just exit

  //TFile* rootfile = new TFile( out_file, "READ" );
}
