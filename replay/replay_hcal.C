#if !defined(__CLING__) || defined(__ROOTCLING__)
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
#include "SBSEArm.h"
#include "SBSBBShower.h"
#include "SBSHCal.h"
#include "SBSTimingHodoscope.h"
#endif

// Simple example replay script
//
// Ole Hansen, 11 April 2016
void replay_hcal(int run_number = 124, uint nev = -1, uint nseg = 0)
{
  //load SBS-offline
  gSystem->Load("libsbs.so");
  //--- Define the experimental configuration, i.e. spectrometers, detectors ---
  
  //THaHRS* bb = new THaHRS("R", "Right HRS" );
  //hrs->AddDetector( new THaVDC("vdc", "Vertical Drift Chambers") );
  //gHaApps->Add(hrs);
  
  //Ideal beam (perfect normal incidence and centering)
  //THaIdealBeam* ib = new THaIdealBeam("IB", "Ideal beam");
  //gHaApps->Add(ib);

  SBSEArm *harm = new SBSEArm("sbs","Hadron Arm with HCal");
  SBSHCal* hcal =  new SBSHCal("hcal","HCAL");
  hcal->SetStoreRawHits(kTRUE);
  harm->AddDetector(hcal);

  SBSGenericDetector* trig= new SBSGenericDetector("trig","HCal trigs");
  trig->SetModeADC(SBSModeADC::kWaveform);
  trig->SetStoreRawHits(kTRUE);
  harm->AddDetector( trig );
  
  gHaApps->Add(harm);
  
  //--- Set up the run we want to replay ---
  
  // This often requires a bit of coding to search directories, test
  // for non-existent files, etc.
  //Variables for searching for split data files.
  //int seg = 0;
  bool seg_ok = true;
  TString exp = "hcal";
  // Create file name patterns.
  string firstname = "hcal_trigtest_%d";
  //string firstname = "e1209019_trigtest_%d";
  //string firstname = "hcal_adc_tdc_%d";
  //string firstname = "hcal_%d";

  THaAnalyzer* analyzer = new THaAnalyzer;
  
  while(seg_ok) {
    string endname = Form(".evio.%d",nseg);
    //string endname = Form(".dat.%d",nseg);
    //string endname = Form(".evio");
    string combined(string(firstname)+endname);
    const char* RunFileNamePattern = combined.c_str();
    vector<TString> pathList;
    pathList.push_back(".");
    //string data_path = "/adaqfs/home/a-onl/sbs/data/";
    //pathList.push_back(Form("%s/data","/adaqfs/home/a-onl/sbs"));
    string data_path = "/volatile/halla/sbs/seeds/evio/";
    pathList.push_back(Form("%s/evio","/volatile/halla/sbs/seeds"));
    //pathList.push_back(Form("%s/data","/adaqfs/home/a-onl/skbarcus"));
    
    // Check if segment exits
    string run_name = Form(RunFileNamePattern, run_number);
    string seg_path_str = data_path+run_name;
    //Convert the std::string to const char * pointer expected by gSystem->AccessPathName.
    const char * seg_path = seg_path_str.c_str();
    std::cout<<seg_path<<std::endl;
    //if( gSystem->AccessPathName("/adaqfs/home/a-onl/sbs/data/hcal_trigtest_183.evio.0")) {
    
    if( gSystem->AccessPathName(seg_path)) {
      seg_ok = false;
      std::cout << "Segment " << nseg << " not found. Exiting" << std::endl;
      continue;
    }
    else{
      std::cout<<"Found "<<Form(RunFileNamePattern,run_number)<<"."<<std::endl;
    }
    
    THaRun* run = new THaRun( pathList, Form(RunFileNamePattern, run_number) );
    
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
    
    //THaAnalyzer* analyzer = new THaAnalyzer;
    
    //TString out_dir = gSystem->Getenv("OUT_DIR");
    TString out_dir ="/volatile/halla/sbs/seeds/rootfiles";
    if( out_dir.IsNull() )  out_dir = ".";
    TString out_file = out_dir + "/" + exp + Form("_%d_%d.root", run_number,nev);
    
    analyzer->SetOutFile( out_file );
    
    analyzer->SetCutFile( "replay_hcal.cdef" );
    analyzer->SetOdefFile( "replay_hcal.odef" );
    
    analyzer->SetVerbosity(2);  // write cut summary to stdout
    analyzer->EnableBenchmarks();
    
    //--- Analyze the run ---
    // Here, one could replay more than one run with
    // a succession of Process calls. The results would all go into the
    // same ROOT output file
    
    run->Print();
    
    
    analyzer->Process(run);
    nseg++;
    }

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
