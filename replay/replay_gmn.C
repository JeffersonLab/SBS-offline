R__ADD_INCLUDE_PATH($SBS/include)
R__ADD_LIBRARY_PATH($SBS/lib64)
R__ADD_LIBRARY_PATH($SBS/lib)
R__LOAD_LIBRARY(libsbs.so)

#if !defined(__CLING__) || defined(__ROOTCLING__)
#include <iostream>

#include "TSystem.h"
#include "THaGlobals.h"
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
#include "SBSBBTotalShower.h"
#include "SBSGRINCH.h"
#include "SBSEArm.h"
#include "SBSHCal.h"
#include "SBSGEMStand.h"
#include "SBSTimingHodoscope.h"

#include "SBSSimDecoder.h"
#endif

void replay_gmn_test(const char* filebase, uint nev = -1, TString experiment="gmn")
{
  SBSBigBite* bigbite = new SBSBigBite("bb", "BigBite spectrometer" );
  //bigbite->AddDetector( new SBSBBShower("ps", "BigBite preshower") );
  //bigbite->AddDetector( new SBSBBShower("sh", "BigBite shower") );
  bigbite->AddDetector( new SBSBBTotalShower("ts", "sh", "ps", "BigBite shower") );
  bigbite->AddDetector( new SBSGRINCH("grinch", "GRINCH PID") );
  bigbite->AddDetector( new SBSTimingHodoscope("hodo", "timing hodo") );
  bigbite->AddDetector( new SBSGEMSpectrometerTracker("gem", "GEM tracker") );
  gHaApps->Add(bigbite);
  
  SBSEArm *harm = new SBSEArm("sbs","Hadron Arm with HCal");
  harm->AddDetector( new SBSHCal("hcal","HCAL") );
  gHaApps->Add(harm);

  //bigbite->SetDebug(2);
  //harm->SetDebug(2);

  THaAnalyzer* analyzer = new THaAnalyzer;
  
  THaInterface::SetDecoder( SBSSimDecoder::Class() );
  
  TString run_file = Form("%s.root", filebase);
  if(std::getenv("DATA_DIR")){
    run_file = Form("%s/%s.root", std::string(std::getenv("DATA_DIR")).c_str(), filebase);
  }
  
  if( gSystem->AccessPathName(run_file) ) {
    Error("replay.C", "Input file does not exist: %s", run_file.Data() );
    exit(-1);
  }
  
  THaRunBase *run = new SBSSimFile(run_file.Data(), "gmn", "");
  run->SetFirstEvent(0);

  cout << "Number of events to replay (-1=all)? ";
  //if( nev > 0 )
  //run->SetFirstEvent(110);
  run->SetLastEvent(nev);
  
  run->SetDataRequired(0);
  run->SetDate(TDatime());
  
  TString out_dir = gSystem->Getenv("OUT_DIR");
  if( out_dir.IsNull() )
    out_dir = ".";
  TString out_file = out_dir + "/" + Form("replayed_%s.root", filebase);
  
  analyzer->SetOutFile( out_file.Data() );
  cout << "output file " << out_file.Data() << " set up " << endl; 
  // File to record cuts accounting information
  analyzer->SetSummaryFile("sbs_hcal_test.log"); // optional

  // Change the cratemap to point to the sim one
  analyzer->SetCrateMapFileName("sbssim_cratemap");

  cout << "sim crate map setup " << endl;
  
  analyzer->SetCutFile( "replay_gmn.cdef" );
  analyzer->SetOdefFile( "replay_gmn.odef" );
  
  cout << "cut file and out file processed " << endl;
  
  analyzer->SetVerbosity(2);  // write cut summary to stdout
  analyzer->EnableBenchmarks();
  
  run->Print();
  
  cout << "about to process " << endl;
  analyzer->Process(run);

  // Clean up

  analyzer->Close();
  delete analyzer;
  //gHaCuts->Clear();
  gHaVars->Clear();
  gHaPhysics->Delete();
  gHaApps->Delete();
 
  
}

int main(int argc, char *argv[])
{
  new THaInterface( "The Hall A analyzer", &argc, argv, 0, 0, 1 );
  //const 
  string filebase; 
  uint nev = -1;
  if(argc<2 || argc>3){
    cout << "Usage: replay_gmn filebase(char*) nev(uint)" << endl;
    return -1;
  }
  filebase = argv[1];
  if(argc==3) nev = atoi(argv[2]);

  replay_gmn_test(filebase.c_str(), nev);
  return 0;
}
