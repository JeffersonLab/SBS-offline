#include <iostream>
#include <string>

#include "THaApparatus.h"
#include "THaTrackingDetector.h"
#include "THaRunBase.h"
#include "THaCrateMap.h"
#include "THaAnalysisObject.h"

#include "THaSpectrometer.h"
#include "SBSGEMPolarimeterTracker.h"
#include "SBSGEMModule.h"
#include "THaTrack.h"
#include "TClonesArray.h"
#include "Helper.h"

using namespace Podd;

SBSGEMPolarimeterTracker::SBSGEMPolarimeterTracker( const char* name, const char* desc, THaApparatus* app ):
  THaNonTrackingDetector(name,desc,app), SBSGEMTrackerBase() {

  fModules.clear();
  fModulesInitialized = false;
  
  fIsMC = false;//by default!
  //fCrateMap = 0;
  fPedestalMode = false;
  fSubtractPedBeforeCommonMode = false;
  fPedMode_DBoverride = false; //Only if the user script invokes the SetPedestalMode method do we override the database value:

  fNegSignalStudy = false;
  
  fIsSpectrometerTracker = false; //used by tracker base
  fIsPolarimeterTracker = true;

  fMultiTrackSearch = false; //Use the "best" track from front tracker, one track search only;

  fUseOpticsConstraint = false;
  fUseForwardOpticsConstraint = false;
  fUseSlopeConstraint = false;

  // fUseFrontTrackConstraint = false;
  // fFrontTrackInitialized = false;

  //fFrontTracks = new TClonesArray("THaTrack",1);
}

SBSGEMPolarimeterTracker::~SBSGEMPolarimeterTracker(){
  //fFrontTracks->Clear("C");
  //delete fFrontTracks;
}


THaAnalysisObject::EStatus SBSGEMPolarimeterTracker::Init( const TDatime& date ){
  //assert( fCrateMap == 0 );
 
  // Why THaTrackingDetector::Init() here? THaTrackingDetector doesn't implement its own Init() method
  // Is that only because this class inherits THaTrackingDetector directly? Is there any advantage to this
  // over calling THaAnalysisObject::Init directly?
  THaAnalysisObject::EStatus status = THaNonTrackingDetector::Init(date);
  //Note: hopefully this triggers the calling of ReadDatabase for the tracker itself!

  if( status == kOK ){

    //    int modcounter=0;

    for( auto& module: fModules ) {
      status = module->Init(date); //This triggers calling of ReadDatabase for each module (I hope)!
      if( status != kOK ){
	return status;
      }
    }

    CompleteInitialization();

    //I'm not sure we really actually want any of these lines of code for the polarimeter context:
    // if( !fFrontTrackInitialized ){
    
    //   new( (*fFrontTracks)[0] ) THaTrack();
    //   fFrontTrackInitialized = true;
    // }
    
  } else {
    return kInitError;
  }

  
  return kOK;
}

Int_t SBSGEMPolarimeterTracker::ReadDatabase( const TDatime& date ){
  std::cout << "[Reading SBSGEMPolarimeterTracker database]" << std::endl;

  fIsInit = kFALSE;

  FILE* file = OpenFile( date );
  if( !file ) return kFileError;

  std::string modconfig;
  //std::vector<Int_t> *cmap = new std::vector<Int_t>;
  //it appears that cmap is not actually used in any way. TBD

  //copying commented structure of DBRequest here for reference (see VarDef.h):
  // struct DBRequest {
  //   const char*      name;     // Key name
  //   void*            var;      // Pointer to data (default to Double*)
  //   VarType          type;     // (opt) data type (see VarType.h, default Double_t)
  //   UInt_t           nelem;    // (opt) number of array elements (0/1 = 1 or auto)
  //   Bool_t           optional; // (opt) If true, missing key is ok
  //   Int_t            search;   // (opt) Search for key along name tree
  //   const char*      descript; // (opt) Key description (if 0, same as name)
  // };
  //std::string pedfilename = ""; //Optional: load pedestals from separate file

  // std::cout << "before LoadDB, pedestalmode = " << fPedestalMode << ", fSubtractPedBeforeCommonMode = " << fSubtractPedBeforeCommonMode << " fPedModeDBoverride = "
  // 	    << fPedMode_DBoverride << std::endl;
  
  int pedestalmode_flag = fPedestalMode ? pow(-1,fSubtractPedBeforeCommonMode) : 0;
  int doefficiency_flag = fMakeEfficiencyPlots ? 1 : 0;
  //int onlinezerosuppressflag = fOnlineZeroSuppression ? 1 : 0;
  int useconstraintflag = fUseConstraint ? 1 : 0; //use constraint on track search region from other detectors in the parent THaSpectrometer (or other)
  //int usefronttrackconstraintflag = fUseFrontTrackConstraint ? 1 : 0; 
  
  int mc_flag = fIsMC ? 1 : 0;
  int fasttrack_flag = fTryFastTrack ? 1 : 0;
  //int useforwardopticsconstraint = fUseForwardOpticsConstraint ? 1 : 0;
  int negsignalstudy_flag = fNegSignalStudy ? 1 : 0;
  int usetrigtime = fUseTrigTime ? 1 : 0;

  int multitracksearch = fMultiTrackSearch ? 1 : 0;

  int nontrackmode = fNonTrackingMode ? 1 : 0;
  
  //  std::vector<int> mingoodhits; 
  //std::vector<double> chi2cut_space;
  //std::vector<double> chi2cut_hits;
  
  DBRequest request[] = {
    { "modules",  &modconfig, kString, 0, 0, 1 }, //read the list of modules:
    { "pedfile",  &fpedfilename, kString, 0, 1 },
    { "cmfile",  &fcmfilename, kString, 0, 1 },
    { "is_mc",        &mc_flag,    kInt, 0, 1, 1 }, //NOTE: is_mc can also be defined via the constructor in the replay script
    { "minhitsontrack", &fMinHitsOnTrack, kInt, 0, 1},
    { "maxhitcombos", &fMaxHitCombinations, kInt, 0, 1},
    { "maxhitcombos_inner", &fMaxHitCombinations_InnerLayers, kInt, 0, 1},
    { "maxhitcombos_total", &fMaxHitCombinations_Total, kDouble, 0, 1},
    { "tryfasttrack", &fasttrack_flag, kInt, 0, 1 },
    { "gridbinwidthx", &fGridBinWidthX, kDouble, 0, 1},
    { "gridbinwidthy", &fGridBinWidthY, kDouble, 0, 1},
    { "gridedgetolerancex", &fGridEdgeToleranceX, kDouble, 0, 1},
    { "gridedgetolerancey", &fGridEdgeToleranceY, kDouble, 0, 1},
    { "trackchi2cut", &fTrackChi2Cut, kDoubleV, 0, 1},
    { "useconstraint", &useconstraintflag, kInt, 0, 1},
    { "sigmahitpos", &fSigma_hitpos, kDouble, 0, 1},
    { "pedestalmode", &pedestalmode_flag, kInt, 0, 1, 1},
    { "do_neg_signal_study", &negsignalstudy_flag, kUInt, 0, 1, 1}, //(optional, search): toggle doing negative signal analysis
    { "do_efficiencies", &doefficiency_flag, kInt, 0, 1, 1},
    { "dump_geometry_info", &fDumpGeometryInfo, kInt, 0, 1, 1},
    { "efficiency_bin_width_1D", &fBinSize_efficiency1D, kDouble, 0, 1, 1 },
    { "efficiency_bin_width_2D", &fBinSize_efficiency2D, kDouble, 0, 1, 1 },
    { "usetrigtime", &usetrigtime, kInt, 0, 1, 1 },
    { "trigtime_crate", &fCrate_RefTime, kUInt, 0, 1, 1 },
    { "trigtime_slot", &fSlot_RefTime, kUInt, 0, 1, 1 },
    { "trigtime_chan", &fChan_RefTime, kUInt, 0, 1, 1 },
    { "trigtime_t0", &fRefTime_Offset, kDouble, 0, 1, 1},
    { "trigtime_calib", &fRefTime_CAL, kDouble, 0, 1, 1},
    { "use_enhanced_chi2", &fUseEnhancedChi2, kInt, 0, 1, 1},
    { "trackchi2cut_hitquality", &fTrackChi2CutHitQuality, kDoubleV, 0, 1, 1},
    { "minhighqualityhitsontrack", &fMinHighQualityHitsOnTrack, kIntV, 0, 1, 1},
    { "sigmatrackt0", &fSigmaTrackT0, kDouble, 0, 1, 1 },
    { "cuttrackt0", &fCutTrackT0, kDouble, 0, 1, 1 },
    { "multitracksearch", &multitracksearch, kInt, 0, 1, 1},
    { "nontrackingmode", &nontrackmode, kInt, 0, 1, 1},
    {0}
  };

  
  
  Int_t err = LoadDB( file, date, request, fPrefix );
  fclose(file);
  if( err )
    return err;

  fMakeEfficiencyPlots = (doefficiency_flag != 0 );
  if( !fPedMode_DBoverride ) {
    //std::cout << "pedestalmode_flag = " << pedestalmode_flag << std::endl;
    
    fPedestalMode = ( pedestalmode_flag != 0 );
    fSubtractPedBeforeCommonMode = ( pedestalmode_flag < 0 );
  }

  if( !fNonTrackingMode_DBoverride ){
    fNonTrackingMode = ( nontrackmode != 0 );
  }
  
  fNegSignalStudy = negsignalstudy_flag != 0;

  fIsMC = (mc_flag != 0);
  fTryFastTrack = (fasttrack_flag != 0);
  
  //fOnlineZeroSuppression = (onlinezerosuppressflag != 0);
  fUseConstraint = (useconstraintflag != 0);

  //fUseFrontTrackConstraint = (usefronttrackconstraintflag != 0);
  fMinHitsOnTrack = std::max(3,fMinHitsOnTrack);

  if( fPedestalMode ){ //then we will just dump raw data to the tree and/or histograms:
    //fZeroSuppress = false;
    fMakeEfficiencyPlots = false; //If in pedestal mode, we don't want efficiency plots
  }

  fUseTrigTime = (usetrigtime != 0);

  fMultiTrackSearch = (multitracksearch != 0 ); //controls whether multiple track searches are even in principle possible
  
  //vsplit is a Podd function that "tokenizes" a string into a vector<string> by whitespace:
  std::vector<std::string> modules = vsplit(modconfig);
  if( modules.empty()) {
    Error("", "[SBSGEMPolarimeterTracker::ReadDatabase] No modules defined");
  }

  int modcounter = 0;

  

  //If a previous module configuration exists, check if anything changed.
  //We will still re-initialize the modules (for now), but we just need to decide whether to re-allocate the modules:

  if( fModulesInitialized ){
    //first check if the size changed:
    if( fModules.size() != modules.size() ) { //If the size of the configuration changed, we know that we have to re-initialize everything!
      DeleteContainer(fModules);
      fModules.resize( modules.size() );
      fModulesInitialized = false;
    } else { //size stayed the same, but we need to check if any of the names changed:

      for( auto it = modules.begin(); it != modules.end(); ++it ) {
	if( fModules[modcounter] ){

	  std::string modname = fModules[modcounter]->GetName();
	  if( modname != *it ){ //name of one or more modules changed, assume we need to re-init everything:
            DeleteContainer(fModules);
	    fModules.resize( modules.size() );
	    fModulesInitialized = false;
	    break;
	  }
	} else { //one or more modules not allocated, re-init everything:
          DeleteContainer(fModules);
	  fModules.resize( modules.size() );
	  fModulesInitialized = false;
	  break;
	}

	modcounter++;
      }
    }
  } else {
    fModules.resize( modules.size() );
  }

  modcounter = 0;

  

  std::cout << "Modules initialized = " << fModulesInitialized << ", number of modules = " << fModules.size() << std::endl;
  //we need to implement some checks to make sure the modules are not already allocated:
  for( const auto& module: modules ) {

    if( !fModulesInitialized ) {
      std::cout << "Initializing module " << module << "... ";
      fModules[modcounter] = new SBSGEMModule(module.c_str(), module.c_str(), this);
      std::cout << " done." << std::endl;
    }
    fModules[modcounter]->fModule = modcounter; //just a dummy index in the module array
    modcounter++;
  }

  if( !fModulesInitialized ) fModulesInitialized = true;

  //Define the number of modules from the "modules" string:
  fNmodules = fModules.size();
  //number of layers should 

  std::cout << "fNmodules = " << fNmodules << std::endl;
  
  //Actually, the loading of the geometry was moved to SBSGEMModule::ReadDatabase(), since this information is specified on a module-by-module basis:
  // Int_t err = ReadGeometry( file, date );
  // if( err ) {
  //     fclose(file);
  //     return err;
  // }
    

  fIsInit = kTRUE;
    
  return kOK;
}


Int_t SBSGEMPolarimeterTracker::Begin( THaRunBase* run ) {
  for( auto& module: fModules ) {
    module->Begin(run);
  }

  TString appname = GetApparatus()->GetName();
  appname += "_";
  TString detname = GetName();

  detname.Prepend(appname);
  
  InitEfficiencyHistos(detname.Data()); //create efficiency histograms (see SBSGEMTrackerBase)
  
  
  return 0;
}

void SBSGEMPolarimeterTracker::Clear( Option_t *opt ){
  
  THaNonTrackingDetector::Clear(opt);

  SBSGEMTrackerBase::Clear();

  //fTrigTime = 0.0;
  
  for( auto& module: fModules ) {
    module->Clear(opt);
  }

  fTrackTheta.clear();
  fTrackPhi.clear();
  fTrackSClose.clear();
  fTrackZClose.clear();

  fFrontTrackIsSet = false;
  
}

Int_t SBSGEMPolarimeterTracker::Decode(const THaEvData& evdata ){
  //return 0;
  //std::cout << "[SBSGEMPolarimeterTracker::Decode], decoding all modules, event ID = " << evdata.GetEvNum() <<  "...";

  //attempt to decode trigger time. No error-checking is done here on the channel info so you'd better provide this correctly in
  //the DB:

  fTrigTime = 0.0;

  //  bool TriggerTimeIsDecoded = false;
  
  if( fUseTrigTime ){
    int ntrigtdchits = evdata.GetNumHits( fCrate_RefTime, fSlot_RefTime, fChan_RefTime );

    int besthit=-1;

    double mintdiff = 1000.0;
    double besttime = 0.0;
    
    //    std::cout << "[SBSGEMPolarimeterTracker::Decode]: ntrigtdchits = " << ntrigtdchits << std::endl;
    
    for( int ihit=0; ihit<ntrigtdchits; ihit++ ){
      UInt_t rawtdc = evdata.GetData( fCrate_RefTime, fSlot_RefTime, fChan_RefTime, ihit );
      UInt_t edge = evdata.GetRawData( fCrate_RefTime, fSlot_RefTime, fChan_RefTime, ihit );

      double rawtime = rawtdc * fRefTime_CAL;
      
      if( edge == 0 ){ //make sure this is a leading-edge TDC
	if(besthit < 0 || fabs( rawtime - fRefTime_Offset ) < mintdiff ){
	  mintdiff = fabs( rawtime - fRefTime_Offset );
	  besthit = ihit;
	  besttime = rawtime - fRefTime_Offset;
	}
      }
    }

    if( besthit >= 0 ){
      fTrigTime = besttime;
      //  TriggerTimeIsDecoded = true;
      //std::cout << "Decoded trigger time, ttrig, offset = " << besttime << ", " << fRefTime_Offset << std::endl; 
    }
  }
  
  //Triggers decoding of each module:

  //Int_t stripcounter = 0;
  for( auto& module: fModules ) {

    module->SetTriggerTime( fTrigTime );
    
    module->Decode(evdata);

    //std::cout << "Decoding module " << (*it)->GetName() << ", nstrips fired = " << (*it)->fNstrips_hit << std::endl;


    //stripcounter += module->fNstrips_hit;
    //std::cout << "done..." << std::endl;
  }

  //std::cout << "done, fNstrips_hit = " << stripcounter << std::endl;
  
  return 0;
}


Int_t SBSGEMPolarimeterTracker::End( THaRunBase* run ){

  UInt_t runnum = run->GetNumber(); 

  TString fname_neg_events;
  
  fname_neg_events.Form("GEM_neg_on_track_run%d.dat",runnum);

  if(fNegSignalStudy)
    PrintNegEvents(fname_neg_events.Data());


  //To automate the printing out of pedestals for database and DAQ, 
  if( fPedestalMode ){
    TString fname_dbped, fname_daqped, fname_dbcm, fname_daqcm;
    
    
    TString specname = GetApparatus()->GetName();
    TString detname = GetName();
    fname_dbped.Form( "db_ped_%s_%s_run%d.dat", specname.Data(), detname.Data(), runnum );
    fname_dbcm.Form( "db_cmr_%s_%s_run%d.dat", specname.Data(), detname.Data(), runnum );
    fname_daqped.Form( "daq_ped_%s_%s_run%d.dat", specname.Data(), detname.Data(), runnum );
    fname_daqcm.Form( "daq_cmr_%s_%s_run%d.dat", specname.Data(), detname.Data(), runnum );
    std::cout<<"\n\n\n"<<specname<<" "<<detname<<"\n\n\n";
    
    //fpedfile_dbase.open( fname_dbase.Data() );
    fCMfile_dbase.open( fname_dbcm.Data() );
    fpedfile_daq.open( fname_daqped.Data() );
    fCMfile_daq.open( fname_daqcm.Data() );
    
    
    TString sdate = run->GetDate().AsString();
    sdate.Prepend( "#" );
    
    TString message;

    message.Form( "# Copy file into sbs-onl@sbsvtp#:~/cfg/pedestals for online pedestal subtraction" );
    fCMfile_daq << sdate << std::endl;
    fCMfile_daq << message << std::endl;
    fCMfile_daq << "# format = crate, slot, mpd, adc_ch, CM min, CM max"
    		  << std::endl;

    message.Form( "# Copy this file into $DB_DIR/gemped to use these pedestals for analysis");
    fCMfile_dbase << sdate << std::endl;
    fCMfile_dbase << message << std::endl;
    fCMfile_dbase << "# format = crate, slot, mpd, adc_ch, CM mean, CM RMS"
		  << std::endl;
    
    message.Form( "# Copy file into sbs-onl@sbsvtp#:~/cfg/pedestals for online pedestal subtraction" );
    
    fpedfile_daq << sdate << std::endl;
    fpedfile_daq <<  message << std::endl;
    fpedfile_daq << "# format = APV        crate       slot       mpd_id       adc_ch followed by " << std::endl
     		 << "# APV channel number      pedestal mean      pedestal rms " << std::endl;

    //message.Form( "# This file defines the common-mode range for the online zero-suppression for the GEM DAQ. Copy its contents into (location TBD) to set these values for detector %s.%s", specname.Data(), detname.Data() );
    
    // fpedfile_cmr << sdate << std::endl;
    // fpedfile_cmr << message << std::endl;
    // fpedfile_cmr << "# format = crate     slot      mpd_id     adc_ch      commonmode min      commonmode max" << std::endl;
  }

  for( auto& module: fModules ) {
    if( fPedestalMode ) { module->PrintPedestals(fCMfile_dbase, fpedfile_daq, fCMfile_daq); }
    module->End(run);
  }

  if( fMakeEfficiencyPlots ){
    CalcEfficiency();
  
    hdidhit_x_layer->Compress();
    hdidhit_y_layer->Compress();
    hdidhit_xy_layer->Compress();

    hshouldhit_x_layer->Compress();
    hshouldhit_y_layer->Compress();
    hshouldhit_xy_layer->Compress();

    hefficiency_x_layer->Compress();
    hefficiency_y_layer->Compress();
    hefficiency_xy_layer->Compress();

    hdidnothit_x_layer->Compress();
    hdidnothit_y_layer->Compress();

    hdidhit_fullreadout_x_layer->Compress();
    hdidhit_fullreadout_y_layer->Compress();

    hneghit_x_layer->Compress();
    hneghit_y_layer->Compress();

    hneghit1D_x_layer->Compress();
    hneghit1D_y_layer->Compress();

    hneghit_good_x_layer->Compress();
    hneghit_good_y_layer->Compress();

    hneghit_good1D_x_layer->Compress();
    hneghit_good1D_y_layer->Compress();

    hdidhit_x_layer->Write(0,kOverwrite);
    hdidhit_y_layer->Write(0,kOverwrite);
    hdidhit_xy_layer->Write(0,kOverwrite);

    hshouldhit_x_layer->Write(0,kOverwrite);
    hshouldhit_y_layer->Write(0,kOverwrite);
    hshouldhit_xy_layer->Write(0,kOverwrite);

    hefficiency_x_layer->Write(0,kOverwrite);
    hefficiency_y_layer->Write(0,kOverwrite);
    hefficiency_xy_layer->Write(0,kOverwrite);

    hdidnothit_x_layer->Write(0,kOverwrite);
    hdidnothit_y_layer->Write(0,kOverwrite);

    hdidhit_fullreadout_x_layer->Write(0,kOverwrite);
    hdidhit_fullreadout_y_layer->Write(0,kOverwrite);

    hneghit_x_layer->Write(0,kOverwrite);
    hneghit_y_layer->Write(0,kOverwrite);

    hneghit1D_x_layer->Write(0,kOverwrite);
    hneghit1D_y_layer->Write(0,kOverwrite);

    hneghit_good_x_layer->Write(0,kOverwrite);
    hneghit_good_y_layer->Write(0,kOverwrite);

    hneghit_good1D_x_layer->Write(0,kOverwrite);
    hneghit_good1D_y_layer->Write(0,kOverwrite);
  }

  if( fDumpGeometryInfo ){ //Print out geometry info for alignment:

    //Maybe this should go to $OUT_DIR:
    TString fnametemp;
    fnametemp.Form( "GEM_alignment_info_%s_%s_run%d.txt",
		    GetApparatus()->GetName(),
		    GetName(), runnum );

    PrintGeometry( fnametemp.Data() );

    
  }
  
  return 0;
}

void SBSGEMPolarimeterTracker::Print(const Option_t* opt) const {
  std::cout << "GEM Stand " << fName << " with " << fModules.size() << " planes defined:" << std::endl;
  /*
    for (std::vector<SBSGEMModule *>::iterator it = fModules.begin() ; it != fModules.end(); ++it){
    std::cout << "\t"
    (*it)->Print(opt);
    }
  */
  for( const auto& module: fModules ) {
    module->Print(opt);
  }
}


void SBSGEMPolarimeterTracker::SetDebug( Int_t level ) {
  THaNonTrackingDetector::SetDebug(level);
  for( auto& module: fModules ) {
    module->SetDebug(level);
  }
}

Int_t SBSGEMPolarimeterTracker::DefineVariables( EMode mode ){
  if( mode == kDefine and fIsSetup ) return kOK;
  fIsSetup = ( mode == kDefine );
  
  RVarDef vars[] = {
    { "trigtime", "trigger time (ns)", "fTrigTime" },
    { "track.ntrack", "number of tracks found", "fNtracks_found" },
    { "track.nhits", "number of hits on track", "fNhitsOnTrack" },
    { "track.x", "Track X (TRANSPORT)", "fXtrack" }, //might be redundant with spectrometer variables, but probably needed for "non-tracking" version
    { "track.y", "Track Y (TRANSPORT)", "fYtrack" },
    { "track.xp", "Track dx/dz (TRANSPORT)", "fXptrack" },
    { "track.yp", "Track dy/dz (TRANSPORT)", "fYptrack" },
    { "track.chi2ndf", "Track Chi2/ndf", "fChi2Track" },
    { "track.chi2ndf_hitquality", "Track Chi2/ndf for hit ADC and time correlations", "fChi2TrackHitQuality" },
    { "track.besttrack", "Index of 'best' track", "fBestTrackIndex" },
    { "track.ngoodhits", "Number of high quality hits on track", "fNgoodhitsOnTrack" },
    { "track.t0", "Track t0 time (weighted average of hit times relative to expected, ns)", "fT0track" },
    { "track.theta", "Track polar theta wrt front track", "fTrackTheta" },
    { "track.phi", "Track polar phi wrt front track", "fTrackPhi" },
    { "track.sclose", "Track distance of closest approach wrt front track", "fTrackSClose" },
    { "track.zclose", "Track point of closest approach wrt front track", "fTrackZClose" },
    { "hit.ngoodhits", "Total number of hits on all found tracks", "fNgoodhits" },
    { "hit.trackindex", "Index of track containing this hit", "fHitTrackIndex" },
    { "hit.module", "Module index of this hit", "fHitModule" },
    { "hit.layer", "Layer index of this hit", "fHitLayer" },
    { "hit.nstripu", "number of U strips on this hit", "fHitNstripsU" },
    { "hit.nstripv", "number of V strips on this hit", "fHitNstripsV" },
    { "hit.ustripmax", "index of u strip with max ADC in this hit", "fHitUstripMax" },
    { "hit.vstripmax", "index of v strip with max ADC in this hit", "fHitVstripMax" },
    { "hit.ustriplo", "index of minimum u strip in this hit", "fHitUstripLo" },
    { "hit.vstriplo", "index of minimum v strip in this hit", "fHitVstripLo" },
    { "hit.ustriphi", "index of maximum u strip in this hit", "fHitUstripHi" },
    { "hit.vstriphi", "index of maximum v strip in this hit", "fHitVstripHi" },
    { "hit.u", "reconstructed hit position along u", "fHitUlocal" },
    { "hit.v", "reconstructed hit position along v", "fHitVlocal" },
    { "hit.xlocal", "reconstructed local x position of hit (internal module coordinates)", "fHitXlocal" },
    { "hit.ylocal", "reconstructed local y position of hit (internal module coordinates)", "fHitYlocal" },
    { "hit.xglobal", "reconstructed global x position of hit", "fHitXglobal" },
    { "hit.yglobal", "reconstructed global y position of hit", "fHitYglobal" },
    { "hit.zglobal", "reconstructed global z position of hit", "fHitZglobal" },
    { "hit.umoment", "U cluster moment (consult source code or A. Puckett for definition)", "fHitUmoment" },
    { "hit.vmoment", "V cluster moment (consult source code or A. Puckett for definition)", "fHitVmoment" },
    { "hit.usigma", "U cluster rms", "fHitUsigma" },
    { "hit.vsigma", "V cluster rms", "fHitVsigma" },
    { "hit.residu", "u hit residual with fitted track (inclusive method)", "fHitResidU" },
    { "hit.residv", "v hit residual with fitted track (inclusive method)", "fHitResidV" },
    { "hit.eresidu", "u hit residual with fitted track (exclusive method)", "fHitEResidU" },
    { "hit.eresidv", "v hit residual with fitted track (exclusive method)", "fHitEResidV" },
    { "hit.ADCU", "cluster ADC sum, U strips", "fHitUADC" },
    { "hit.ADCV", "cluster ADC sum, V strips", "fHitVADC" },
    { "hit.ADCU_deconv", "cluster deconvoluted ADC sum, U strips", "fHitUADCclust_deconv" },
    { "hit.ADCV_deconv", "cluster deconvoluted ADC sum, V strips", "fHitVADCclust_deconv" },
    { "hit.ADCmaxsampUclust_deconv", "max deconv. cluster-summed ADC sample, U strips", "fHitUADCclust_maxsamp_deconv" },
    { "hit.ADCmaxsampVclust_deconv", "max deconv. cluster-summed ADC sample, V strips", "fHitVADCclust_maxsamp_deconv" },
    { "hit.ADCmaxcomboUclust_deconv", "max deconv. cluster-summed two-sample combo, U strips", "fHitUADCclust_maxcombo_deconv" },
    { "hit.ADCmaxcomboVclust_deconv", "max deconv. cluster-summed two-sample comboe, V strips", "fHitVADCclust_maxcombo_deconv" },
    { "hit.ADCavg", "cluster ADC average", "fHitADCavg" },
    { "hit.ADCavg_deconv", "cluster ADC average deconvoluted (ADCU+ADCV)/2", "fHitADCavg_deconv" },
    { "hit.ADCmaxstripU", "ADC sum of max U strip", "fHitUADCmaxstrip" },
    { "hit.ADCmaxstripV", "ADC sum of max V strip", "fHitVADCmaxstrip" },
    { "hit.ADCmaxsampU", "max sample of max U strip", "fHitUADCmaxsample" },
    { "hit.ADCmaxsampV", "max sample of max V strip", "fHitVADCmaxsample" },
    { "hit.ADCmaxsampUclust", "max U cluster-summed ADC time sample", "fHitUADCmaxclustsample" },
    { "hit.ADCmaxsampVclust", "max V cluster-summed ADC time sample", "fHitVADCmaxclustsample" },
    { "hit.DeconvADCmaxstripU", "deconv ADC sum of max U strip", "fHitUADCmaxstrip_deconv" },
    { "hit.DeconvADCmaxstripV", "deconv ADC sum of max V strip", "fHitVADCmaxstrip_deconv" },
    { "hit.DeconvADCmaxsampU", "deconv max sample of max U strip", "fHitUADCmaxsample_deconv" },
    { "hit.DeconvADCmaxsampV", "deconv max sample of max V strip", "fHitVADCmaxsample_deconv" },
    { "hit.DeconvADCmaxcomboU", "deconv max two-sample combo of max U strip", "fHitUADCmaxcombo_deconv" },
    { "hit.DeconvADCmaxcomboV", "deconv max two-sample combo of max V strip", "fHitVADCmaxcombo_deconv" },
    { "hit.ADCasym", "Hit ADC asymmetry: (ADCU - ADCV)/(ADCU + ADCV)", "fHitADCasym" },
    { "hit.ADCasym_deconv", "Hit ADC asymmetry (deconv): (ADCU-ADCV)/(ADCU+ADCV)", "fHitADCasym_deconv" },
    { "hit.Utime", "cluster timing based on U strips", "fHitUTime" },
    { "hit.Vtime", "cluster timing based on V strips", "fHitVTime" },
    { "hit.UtimeDeconv", "U strip deconvoluted cluster time", "fHitUTimeDeconv" },
    { "hit.VtimeDeconv", "V strip deconvoluted cluster time", "fHitVTimeDeconv" },
    { "hit.UtimeFit", "cluster timing based on U strips", "fHitUTimeFit" },
    { "hit.VtimeFit", "cluster timing based on V strips", "fHitVTimeFit" },
    { "hit.UtimeMaxStrip", "cluster timing based on U strips", "fHitUTimeMaxStrip" },
    { "hit.VtimeMaxStrip", "cluster timing based on V strips", "fHitVTimeMaxStrip" },
    { "hit.UtimeMaxStripDeconv", "cluster timing based on U strips", "fHitUTimeMaxStripDeconv" },
    { "hit.VtimeMaxStripDeconv", "cluster timing based on V strips", "fHitVTimeMaxStripDeconv" },
    { "hit.UtimeMaxStripFit", "Strip fitted t0 for max strip in cluster", "fHitUTimeMaxStripFit" },
    { "hit.VtimeMaxStripFit", "Strip fitted t0 for max strip in cluster", "fHitVTimeMaxStripFit" },
    { "hit.deltat", "cluster U time - V time", "fHitDeltaT" },
    { "hit.Tavg", "hit T average", "fHitTavg" },
    { "hit.deltat_deconv", "deconvoluted cluster U time - Vtime", "fHitDeltaTDeconv" },
    { "hit.Tavg_deconv", "deconvoluted average U,V cluster time", "fHitTavgDeconv" },
    { "hit.deltat_fit", "cluster U - V fit time", "fHitDeltaTFit" },
    { "hit.Tavg_fit", "cluster (U+V)/2 fit time", "fHitTavgFit" },
    { "hit.isampmaxUclust", "peak time sample in cluster-summed U ADC samples", "fHitIsampMaxUclust" },
    { "hit.isampmaxVclust", "peak time sample in cluster-summed V ADC samples", "fHitIsampMaxVclust" },
    { "hit.isampmaxUclustDeconv", "peak time sample max. deconv. U cluster sum", "fHitIsampMaxUclustDeconv" },
    { "hit.isampmaxVclustDeconv", "peak time sample max. deconv. V cluster sum", "fHitIsampMaxVclustDeconv" },
    { "hit.isampmaxUstrip", "peak time sample in max U strip", "fHitIsampMaxUstrip" },
    { "hit.isampmaxVstrip", "peak time sample in max V strip", "fHitIsampMaxVstrip" },
    { "hit.isampmaxUstripDeconv", "peak time sample in max U strip", "fHitIsampMaxUstripDeconv" },
    { "hit.isampmaxVstripDeconv", "peak time sample in max V strip", "fHitIsampMaxVstripDeconv" },
    { "hit.icombomaxUstripDeconv", "peak time sample combo in max U strip", "fHitIcomboMaxUstripDeconv" },
    { "hit.icombomaxVstripDeconv", "peak time sample combo in max V strip", "fHitIcomboMaxVstripDeconv" },
    { "hit.icombomaxUclustDeconv", "max two-sample combo deconvoluted U cluster", "fHitIcomboMaxUclustDeconv" },
    { "hit.icombomaxVclustDeconv", "max two-sample combo deconvoluted V cluster", "fHitIcomboMaxVclustDeconv" },
    { "hit.ccor_clust", "correlation coefficient between cluster-summed U and V samples", "fHitCorrCoeffClust" },
    { "hit.ccor_strip", "correlation coefficient between U and V samples on strips with max ADC", "fHitCorrCoeffMaxStrip" },
    { "hit.ccor_clust_deconv", "correlation coefficient between cluster-summed U and V deconvoluted samples", "fHitCorrCoeffClustDeconv" },
    { "hit.ccor_strip_deconv", "correlation coefficient between U and V max strip deconvoluted samples", "fHitCorrCoeffMaxStripDeconv" },
    { "hit.ENABLE_CM_U", "Enable CM flag for max U strip in this hit", "fHitU_ENABLE_CM" },
    { "hit.ENABLE_CM_V", "Enable CM flag for max V strip in this hit", "fHitV_ENABLE_CM" },
    { "hit.CM_GOOD_U", "Enable CM flag for max U strip in this hit", "fHitU_CM_GOOD" },
    { "hit.CM_GOOD_V", "Enable CM flag for max V strip in this hit", "fHitV_CM_GOOD" },
    { "hit.BUILD_ALL_SAMPLES_U", "Enable CM flag for max U strip in this hit", "fHitU_BUILD_ALL_SAMPLES" },
    { "hit.BUILD_ALL_SAMPLES_V", "Enable CM flag for max V strip in this hit", "fHitV_BUILD_ALL_SAMPLES" },
    { "hit.ADCfrac0_Umax", "Max U strip ADC0/ADCsum", "fHitADCfrac0_MaxUstrip" },
    { "hit.ADCfrac1_Umax", "Max U strip ADC1/ADCsum", "fHitADCfrac1_MaxUstrip" },
    { "hit.ADCfrac2_Umax", "Max U strip ADC2/ADCsum", "fHitADCfrac2_MaxUstrip" },
    { "hit.ADCfrac3_Umax", "Max U strip ADC3/ADCsum", "fHitADCfrac3_MaxUstrip" },
    { "hit.ADCfrac4_Umax", "Max U strip ADC4/ADCsum", "fHitADCfrac4_MaxUstrip" },
    { "hit.ADCfrac5_Umax", "Max U strip ADC5/ADCsum", "fHitADCfrac5_MaxUstrip" },
    { "hit.ADCfrac0_Vmax", "Max V strip ADC0/ADCsum", "fHitADCfrac0_MaxVstrip" },
    { "hit.ADCfrac1_Vmax", "Max V strip ADC1/ADCsum", "fHitADCfrac1_MaxVstrip" },
    { "hit.ADCfrac2_Vmax", "Max V strip ADC2/ADCsum", "fHitADCfrac2_MaxVstrip" },
    { "hit.ADCfrac3_Vmax", "Max V strip ADC3/ADCsum", "fHitADCfrac3_MaxVstrip" },
    { "hit.ADCfrac4_Vmax", "Max V strip ADC4/ADCsum", "fHitADCfrac4_MaxVstrip" },
    { "hit.ADCfrac5_Vmax", "Max V strip ADC5/ADCsum", "fHitADCfrac5_MaxVstrip" },
    { "hit.DeconvADC0_Umax", "Max U strip Deconv ADC0", "fHitDeconvADC0_MaxUstrip" },
    { "hit.DeconvADC1_Umax", "Max U strip Deconv ADC1", "fHitDeconvADC1_MaxUstrip" },
    { "hit.DeconvADC2_Umax", "Max U strip deconv ADC2", "fHitDeconvADC2_MaxUstrip" },
    { "hit.DeconvADC3_Umax", "Max U strip deconv ADC3", "fHitDeconvADC3_MaxUstrip" },
    { "hit.DeconvADC4_Umax", "Max U strip deconv ADC4", "fHitDeconvADC4_MaxUstrip" },
    { "hit.DeconvADC5_Umax", "Max U strip deconv ADC5", "fHitDeconvADC5_MaxUstrip" },
    { "hit.DeconvADC0_Vmax", "Max V strip deconv ADC0", "fHitDeconvADC0_MaxVstrip" },
    { "hit.DeconvADC1_Vmax", "Max V strip deconv ADC1", "fHitDeconvADC1_MaxVstrip" },
    { "hit.DeconvADC2_Vmax", "Max V strip deconv ADC2", "fHitDeconvADC2_MaxVstrip" },
    { "hit.DeconvADC3_Vmax", "Max V strip deconv ADC3", "fHitDeconvADC3_MaxVstrip" },
    { "hit.DeconvADC4_Vmax", "Max V strip deconv ADC4", "fHitDeconvADC4_MaxVstrip" },
    { "hit.DeconvADC5_Vmax", "Max V strip deconv ADC5", "fHitDeconvADC5_MaxVstrip" },
    { "hit.TSchi2_Umax", "Max U strip TS chi2", "fHitTSchi2MaxUstrip" },
    { "hit.TSchi2_Vmax", "Max V strip TS chi2", "fHitTSchi2MaxVstrip" },
    { "hit.TSprob_Umax", "Max U strip TS prob", "fHitTSprobMaxUstrip" },
    { "hit.TSprob_Vmax", "Max V strip TS prob", "fHitTSprobMaxVstrip" },
    { "hit.Tavg_corr", "Corrected hit time (ns)", "fHitTavgCorrected" },
    { "hit.Ugain","Applied gain factor U", "fHitUgain" },
    { "hit.Vgain","Applied gain factor V", "fHitVgain" },
    { "hit.crate_U", "VTP ROC ID of max U strip", "fHitCrate_U" },
    { "hit.mpd_U", "MPD/fiber number of max U strip", "fHitMPD_U" },
    { "hit.adc_id_U", "ADC channel number of APV containing max U strip", "fHitADCID_U" },
    { "hit.crate_V", "VTP ROC ID of max V strip", "fHitCrate_V" },
    { "hit.mpd_V", "MPD/fiber number of max V strip", "fHitMPD_V" },
    { "hit.adc_id_V", "ADC channel number of APV containing max V strip", "fHitADCID_V" },
    { "nlayershit", "number of layers with any strip fired", "fNlayers_hit" },
    { "nlayershitu", "number of layers with any U strip fired", "fNlayers_hitU" },
    { "nlayershitv", "number of layers with any V strip fired", "fNlayers_hitV" },
    { "nlayershituv", "number of layers with at least one 2D hit", "fNlayers_hitUV" },
    { "nstripsu_layer", "total number of U strips fired by layer", "fNstripsU_layer" },
    { "nstripsv_layer", "total number of V strips fired by layer", "fNstripsV_layer" },
    { "nstripsu_layer_neg", "total number of negative U strips fired by layer", "fNstripsU_layer_neg" },
    { "nstripsv_layer_neg", "total number of negative V strips fired by layer", "fNstripsV_layer_neg" },
    { "nstripsu_layer_neg_miss", "total U strips near track by layer with hits not found", "fNstripsU_layer_neg_miss" },
    { "nstripsv_layer_neg_miss", "total V strips near track by layer with hits notfound", "fNstripsV_layer_neg_miss" },
    { "nstripsu_layer_neg_hit", "total U strips near track by layer with hits found", "fNstripsU_layer_neg_hit" },
    { "nstripsv_layer_neg_hit", "total V strips near track by layer with hits found", "fNstripsV_layer_neg_hit" },
    { "nclustu_layer", "total number of U clusters by layer", "fNclustU_layer" },
    { "nclustv_layer", "total number of V clusters by layer", "fNclustV_layer" },
    { "nclustu_layer_neg", "total number of negative U clusters by layer", "fNclustU_layer_neg" },
    { "nclustv_layer_neg", "total number of negative V clusters by layer", "fNclustV_layer_neg" },
    { "n2Dhit_layer", "total_number of 2D hits by layer", "fN2Dhit_layer" },
    //{ "nclustu_layer_miss", "total number of U clusters by layer in a module missing hits", "fNclustU_layer_miss" },
    //{ "nclustv_layer_miss", "total number of V clusters by layer in a module missing hits", "fNclustV_layer_miss" },
    { nullptr }
  };
  DefineVarsFromList( vars, mode );

  return 0;
}


Int_t SBSGEMPolarimeterTracker::CoarseProcess( TClonesArray& tracks ){
  
  if( !fUseConstraint && !fPedestalMode ){
    //std::cout << "SBSGEMPolarimeterTracker::CoarseTrack...";
    //If no external constraints on the track search region are being used/defined, we do the track-finding in CoarseTrack (before processing all the THaNonTrackingDetectors in the parent spectrometer):
    //std::cout << "calling find_tracks..." << std::endl;
    if( !ftracking_done ) find_tracks();

    if( fFrontTrackIsSet && fNtracks_found > 0 ){ //there is currently no scenario is SBSEArm/GEN-RP context in which this condition could evaluate to true
      CalcScatteringParameters();
    }
    
    //The following lines aren't applicable for PolarimeterTracker:
    // for( int itrack=0; itrack<fNtracks_found; itrack++ ){
    //   THaTrack *Track = AddTrack( tracks, fXtrack[itrack], fYtrack[itrack], fXptrack[itrack], fYptrack[itrack] );

    //   int ndf = 2*fNhitsOnTrack[itrack]-4;
    //   double chi2 = fChi2Track[itrack]*ndf;
      
    //   Track->SetChi2( chi2, ndf );

    //   int index = tracks.GetLast();
    //   Track->SetIndex( index );
    // }

    //std::cout << "done. found " << fNtracks_found << " tracks" << std::endl;
    
  }
  
  return 0;
}
Int_t SBSGEMPolarimeterTracker::FineProcess( TClonesArray& tracks ){


  // In the polarimeter context, when FineProcess gets invoked, the TClonesArray &tracks refers to the tracks reconstructed by any
  // "spectrometer trackers" or "front" trackers upstream of the polarimeter tracker. However, in the GEN-RP context, we
  // also need to worry about the "non-track" front track constraint.
  // In principle, this method can get called more than once.
  // In particular, the SBSEArm class used by GEN-RP will call this method twice, once
  // in SBSEArm::Track() and again in SBSEArm::Reconstruct() which just invokes THaSpectrometer::Reconstruct()!
  // The first time, there will be no front track set and the tracking will not have been done.
  // So the first call to this method from SBSEArm::Track() will cause tracking to happen, which will set
  // ftracking_done to true.
  // The second call to this method from SBSEArm::Reconstruct() will happen AFTER the "front track" parameters are set,
  // based on the results from the front tracking. This will cause CalcScatteringParameters() to be invoked.
  
  if( fUseConstraint && !fPedestalMode ){ //
    //std::cout << "SBSGEMPolarimeterTracker::FineTrack..."; 
    //Calls SBSGEMTrackerBase::find_tracks(), which takes no arguments:
    //std::cout << "calling find_tracks" << std::endl;

    //Is using a constraint, only attempt tracking if constraints have actually been initialized with info from external detectors:
    fConstraintInitialized = fConstraintPoint_Front_IsInitialized && fConstraintPoint_Back_IsInitialized &&
      fConstraintWidth_Front_IsInitialized && fConstraintWidth_Back_IsInitialized;

    if( fConstraintInitialized ){

      //In the polarimeter context, we want to loop on all the tracks in the tracks array that gets passed as argument to to FineProcess as well as any front "pseudotracks" based on a projected incident neutron trajectory 

      if( !ftracking_done ) find_tracks();

    }
    //std::cout << "done. found " << fNtracks_found << " tracks" << std::endl;
    
  }

  // regardless of constraint and whether tracking was done here or in CoarseProcess, if front track was set and
  // we have at least one track, calculate the scattering parameters:
  if( fFrontTrackIsSet && fNtracks_found > 0 ){
    CalcScatteringParameters();
  }
  
  return 0;
}

void SBSGEMPolarimeterTracker::CalcScatteringParameters(){
  //Just in case there is anything leftover here from a previous event:
  fTrackTheta.clear();
  fTrackPhi.clear();
  fTrackSClose.clear();
  fTrackZClose.clear();
  if( !fFrontTrackIsSet ){ //just assign kBig to each track
    for( int itr=0; itr<fNtracks_found; itr++ ){
      fTrackTheta.push_back( 1.e20 );
      fTrackPhi.push_back( 1.e20 );
      fTrackSClose.push_back( 1.e20 );
      fTrackZClose.push_back( 1.e20 );
    }
  } else {
    //Calculate theta, phi, SClose, ZClose (the formulas for these parameters can found in e.g.,
    // Andrew Puckett's Ph.D. thesis). It is ASSUMED that the back tracker coordinate system
    // is the same as the front, if it isn't, then the results of these calculations
    // won't be trustworthy:
    TVector3 FrontDir(fFrontTrackXp,fFrontTrackYp,1.0);
    FrontDir = FrontDir.Unit();
    TVector3 FrontPos( fFrontTrackX, fFrontTrackY, 0.0 );
    for( int itr=0; itr<fNtracks_found; itr++ ){
      TVector3 BackDir( fXptrack[itr], fYptrack[itr], 1.0 );
      BackDir = BackDir.Unit();
      TVector3 BackPos( fXtrack[itr], fYtrack[itr], 0.0 );

      // Theta and phi calculations are simple. Theta is the simplest. Since
      // both direction vectors are unit vectors, we simply take the arccosine of back dot front:
      double theta = acos(BackDir.Dot(FrontDir));
      //For the phi direction, we need to calculate the basis vectors of the comoving coordinate
      //system:
      // The standard convention for azimuthal scattering angles is that we choose the y axis
      // of the comoving coordinate system to be
      // perpendicular to the z axis and parallel to the yz plane of TRANSPORT coordinates, and
      // then we choose the x axis to be perpendicular to both:
      TVector3 xaxis(1,0,0);
      TVector3 yaxis = (FrontDir.Cross(xaxis)).Unit(); //y axis is perpendicular to both z and the x axis of transport coordinates.
      xaxis = (yaxis.Cross(FrontDir)).Unit(); //Final x axis is chosen by requiring orthogonality to both y and z and making sure the coordinate system is right-handed
      // now take the projections of the scattered (back) track along y and x;
      // The phi angle is given by tan(phi) = by/bx
      double phi = atan2( BackDir.Dot( yaxis ), BackDir.Dot( xaxis ) );
      fTrackTheta.push_back( theta ); //in radians
      fTrackPhi.push_back( phi ); //in radians, -PI < phi < PI
      //Next up: distance and point of closest approach. Refer to Andrew Puckett Ph.D. thesis for the
      // formulas (equations 4.26, 4.29, 4.30, etc)
      double a = 1.0 + pow( fFrontTrackXp, 2 ) + pow( fFrontTrackYp, 2 );
      double c = 1.0 + pow( fXptrack[itr], 2 ) + pow( fYptrack[itr], 2 );
      double b = 1.0 + fFrontTrackXp*fXptrack[itr] + fFrontTrackYp*fYptrack[itr];

      double det = a*c-b*b;
      double d1 = -fFrontTrackXp*(fFrontTrackX-fXtrack[itr])-fFrontTrackYp*(fFrontTrackY-fYtrack[itr]);
      double d2 = fXptrack[itr]*(fFrontTrackX-fXtrack[itr])+fYptrack[itr]*(fFrontTrackY-fYtrack[itr]);

      double z1 = (c*d1 + b*d2)/det;
      double z2 = (b*d1 + a*d2)/det;

      fTrackZClose.push_back( 0.5*(z1+z2) );

      double x1close = fFrontTrackX + z1*fFrontTrackXp;
      double x2close = fXtrack[itr]+z2*fXptrack[itr];
      double y1close = fFrontTrackY + z1*fFrontTrackYp;
      double y2close = fYtrack[itr]+z2*fYptrack[itr];
      double s2 = pow(x1close-x2close,2)+pow(y1close-y2close,2)+pow(z1-z2,2);
      fTrackSClose.push_back( sqrt(s2) );
      
    }
  }
}

void SBSGEMPolarimeterTracker::SetFrontTrack( TVector3 pos, TVector3 dir){
  fFrontTrackXp = dir.X()/dir.Z();
  fFrontTrackYp = dir.Y()/dir.Z();
  fFrontTrackX = pos.X() - fFrontTrackXp * pos.Z();
  fFrontTrackY = pos.Y() - fFrontTrackYp * pos.Z();
  fFrontTrackIsSet = true;
}

void SBSGEMPolarimeterTracker::SetFrontTrack( double x, double y, double theta, double phi ){
  fFrontTrackX = x;
  fFrontTrackY = y;
  fFrontTrackXp = theta;
  fFrontTrackYp = phi;
  fFrontTrackIsSet = true;
}

ClassImp(SBSGEMPolarimeterTracker)
