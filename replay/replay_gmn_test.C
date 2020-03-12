#include <iostream>

#include "TSystem.h"
#include "TString.h"
#include "TFile.h"
#include "TList.h"
#include "TObject.h"

#include "THaEvData.h"
#include "THaRun.h"
#include "THaAnalyzer.h"
#include "THaVarList.h"
#include "THaInterface.h"

#include "SBSBigBite.h"
#include "SBSBBShower.h"
#include "SBSEArm.h"
#include "SBSHCal.h"

#include "SBSSimDecoder.h"

void replay_gmn_test(const char* filebase = "simdig_gmn13.5_100evts_NB", uint nev = -1)
{
  SBSBigBite* bigbite = new SBSBigBite("bb", "BigBite spectrometer" );
  bigbite->AddDetector( new SBSBBShower("ps", "BigBite preshower") );
  bigbite->AddDetector( new SBSBBShower("sh", "BigBite shower") );
  gHaApps->Add(bigbite);
  
  SBSEArm *harm = new SBSEArm("sbs","Hadron Arm with HCal");
  harm->AddDetector( new SBSHCal("hcal","HCAL") );
  gHaApps->Add(harm);

  THaAnalyzer* analyzer = new THaAnalyzer;
  
  THaInterface::SetDecoder( SBSSimDecoder::Class() );
  
  TString run_file = Form("digitized/%s.root", filebase);
  if( gSystem->AccessPathName(run_file) ) {
    Error("replay.C", "Input file does not exist: %s", run_file.Data() );
    exit(-1);
  }
  
  THaRunBase *run = new SBSSimFile(run_file.Data());
  run->SetFirstEvent(0);

  cout << "Number of events to replay (-1=all)? ";
  if( nev > 0 )
    run->SetLastEvent(nev);
  
  run->SetDataRequired(0);
  run->SetDate(TDatime());
  
  TString out_dir = gSystem->Getenv("OUT_DIR");
  if( out_dir.IsNull() )
    out_dir = ".";
  TString out_file = out_dir + "/" + Form("replayed_%s.root", filebase);
  
  
  analyzer->SetOutFile( out_file.Data() );
  // File to record cuts accounting information
  analyzer->SetSummaryFile("sbs_hcal_test.log"); // optional

  // Change the cratemap to point to the sim one
  analyzer->SetCrateMapFileName("sbssim_cratemap");
  
  analyzer->SetCutFile( "replay_gmn.cdef" );
  analyzer->SetOdefFile( "replay_gmn.odef" );
  
  analyzer->SetVerbosity(2);  // write cut summary to stdout
  analyzer->EnableBenchmarks();
  
  run->Print();
  
  analyzer->Process(run);

  // Clean up

  analyzer->Close();
  delete analyzer;
  //gHaCuts->Clear();
  gHaVars->Clear();
  gHaPhysics->Delete();
  gHaApps->Delete();
 
  
}
