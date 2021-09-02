#include "TSystem.h"
#include "TList.h"
#include "THaRun.h"
#include "THaEvent.h"
#include "THaAnalyzer.h"
#include "THaApparatus.h"
#include "TString.h"
#include "TClonesArray.h"

#include "SBSGEMSpectrometerTracker.h"
#include "SBSBigBite.h"
//#include "SBSGEMStand.h"
//#include "SBSBigBite.h"

void replay_BBGEM( int runnum=220, int firstsegment=0, int maxsegments=1, long firstevent=0, long nevents=-1 ){

  //  gSystem->Load("libsbs.so");

  SBSBigBite   *bb = new SBSBigBite("bb", "Generic apparatus");
  //SBSGEMStand *gems = new SBSGEMStand("gems", "Collection of GEMs in stand");
  SBSGEMSpectrometerTracker *bbgem = new SBSGEMSpectrometerTracker("gem", "BigBite Hall A GEM data");
    
  bb->AddDetector(bbgem);

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
  
  gHaApps->Add(bb);

  // A simple event class to be output to the resulting tree.
  // Creating your own descendant of THaEvent is one way of
  // defining and controlling the output.
  THaEvent* event = new THaEvent;

  TString prefix = gSystem->Getenv("DATA_DIR");
  
  bool segmentexists = true;
  int segment=firstsegment; 

  TClonesArray *filelist = new TClonesArray("THaRun",10);

  

  //This loop adds all file segments found to the list of THaRuns to process:
  while( segmentexists && segment-firstsegment < maxsegments ){

    TString codafilename;
    codafilename.Form( "%s/bbgem_%d.evio.%d", prefix.Data(), runnum, segment );

    if( gSystem->AccessPathName( codafilename.Data() ) ){
      segmentexists = false;
    } else if( segment == 0 ){
      new( (*filelist)[segment] ) THaRun( codafilename.Data() );
    } else {
      THaRun *rtemp = ( (THaRun*) (*filelist)[segment-1] ); //make otherwise identical copy of previous run in all respects except coda file name:
      new( (*filelist)[segment] ) THaRun( *rtemp );
      ( (THaRun*) (*filelist)[segment] )->SetFilename( codafilename.Data() );
    }
    segment++;
  }

  prefix = gSystem->Getenv("OUT_DIR");

  TString outfilename;
  outfilename.Form( "%s/bbgem_replayed_%d.root", prefix.Data(), runnum );

  // Define the run(s) that we want to analyze.
  // We just set up one, but this could be many.
  //  THaRun* run = new THaRun( "prod12_4100V_TrigRate25_4.dat" );
  //THaRun* run = new THaRun( "5GEM_sample.dat" );
  //THaRun* run = new THaRun( "/Users/puckett/WORK/GEM_ALIGNMENT/RAWDATA/gem_cleanroom_2811.evio.31" );
  //THaRun* run = new THaRun( codafilename.Data() );
  //THaRun* run = new THaRun( "/Users/puckett/WORK/GEM_ALIGNMENT/RAWDATA/gem_cleanroom_2805.evio.0" );

  

  analyzer->SetVerbosity(0);

  analyzer->EnableBenchmarks();
  
  // Define the analysis parameters
  analyzer->SetEvent( event );
  analyzer->SetOutFile( outfilename.Data() );
  // File to record cuts accounting information
  analyzer->SetSummaryFile("summary_example.log"); // optional

  analyzer->SetOdefFile( "replay_BB_TEDF.odef" );
  
  //analyzer->SetCompressionLevel(0); // turn off compression

  filelist->Compress();

  for( int iseg=0; iseg<filelist->GetEntries(); iseg++ ){
    THaRun *run = ( (THaRun*) (*filelist)[iseg] );
    if( nevents > 0 ) run->SetLastEvent(nevents); //not sure if this will work as we want it to for multiple file segments chained together

    run->SetFirstEvent( firstevent );
    
    run->SetDataRequired(0);
    run->SetDate(TDatime());

    analyzer->Process(run);     // start the actual analysis
  }
}
