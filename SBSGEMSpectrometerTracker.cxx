#include <iostream>
#include <string>

#include "THaApparatus.h"
#include "THaTrackingDetector.h"
#include "THaRunBase.h"
#include "THaCrateMap.h"
#include "THaAnalysisObject.h"

#include "THaSpectrometer.h"
#include "SBSGEMSpectrometerTracker.h"
#include "SBSGEMModule.h"
#include "THaTrack.h"
#include "TClonesArray.h"
#include "Helper.h"

using namespace Podd;

SBSGEMSpectrometerTracker::SBSGEMSpectrometerTracker( const char* name, const char* desc, THaApparatus* app ):
  THaTrackingDetector(name,desc,app), SBSGEMTrackerBase() {

  fModules.clear();
  fModulesInitialized = false;
  
  fIsMC = false;//by default!
  //fCrateMap = 0;
  fPedestalMode = false;
  fSubtractPedBeforeCommonMode = false;
  fPedMode_DBoverride = false; //Only if the user script invokes the SetPedestalMode method do we override the database value:
  
  fIsSpectrometerTracker = true; //used by tracker base
  fUseOpticsConstraint = false;
  
  fTestTrackInitialized = false;

  fTestTracks = new TClonesArray("THaTrack",1);
}

SBSGEMSpectrometerTracker::~SBSGEMSpectrometerTracker(){
  fTestTracks->Clear("C");
  delete fTestTracks;
}


THaAnalysisObject::EStatus SBSGEMSpectrometerTracker::Init( const TDatime& date ){
  //assert( fCrateMap == 0 );
 
  // Why THaTrackingDetector::Init() here? THaTrackingDetector doesn't implement its own Init() method
  // Is that only because this class inherits THaTrackingDetector directly? Is there any advantage to this
  // over calling THaAnalysisObject::Init directly?
  THaAnalysisObject::EStatus status = THaTrackingDetector::Init(date);
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

    if( !fTestTrackInitialized ){
    
      new( (*fTestTracks)[0] ) THaTrack();
      fTestTrackInitialized = true;
    }
    
  } else {
    return kInitError;
  }

 
  return kOK;
}

Int_t SBSGEMSpectrometerTracker::ReadDatabase( const TDatime& date ){
  std::cout << "[Reading SBSGEMSpectrometerTracker database]" << std::endl;

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
  int useconstraintflag = fUseConstraint ? 1 : 0; //use constraint on track search region from other detectors in the parent THaSpectrometer (or other
  int mc_flag = fIsMC ? 1 : 0;
  int fasttrack_flag = fTryFastTrack ? 1 : 0;
  int useopticsconstraint = fUseOpticsConstraint ? 1 : 0;
  int useslopeconstraint = fUseSlopeConstraint ? 1 : 0;
  
  DBRequest request[] = {
    { "modules",  &modconfig, kString, 0, 0, 1 }, //read the list of modules:
    { "pedfile",  &fpedfilename, kString, 0, 1 },
    { "is_mc",        &mc_flag,    kInt, 0, 1, 1 }, //NOTE: is_mc can also be defined via the constructor in the replay script
    { "minhitsontrack", &fMinHitsOnTrack, kInt, 0, 1},
    { "maxhitcombos", &fMaxHitCombinations, kInt, 0, 1},
    { "maxhitcombos_inner", &fMaxHitCombinations_InnerLayers, kInt, 0, 1},
    { "maxhitcombos_total", &fMaxHitCombinations_Total, kLong, 0, 1},
    { "tryfasttrack", &fasttrack_flag, kInt, 0, 1 },
    { "gridbinwidthx", &fGridBinWidthX, kDouble, 0, 1},
    { "gridbinwidthy", &fGridBinWidthY, kDouble, 0, 1},
    { "gridedgetolerancex", &fGridEdgeToleranceX, kDouble, 0, 1},
    { "gridedgetolerancey", &fGridEdgeToleranceY, kDouble, 0, 1},
    { "trackchi2cut", &fTrackChi2Cut, kDouble, 0, 1},
    { "useconstraint", &useconstraintflag, kInt, 0, 1},
    { "useopticsconstraint", &useopticsconstraint, kInt, 0, 1},
    { "sigmahitpos", &fSigma_hitpos, kDouble, 0, 1},
    { "pedestalmode", &pedestalmode_flag, kInt, 0, 1, 1},
    { "do_efficiencies", &doefficiency_flag, kInt, 0, 1, 1},
    { "dump_geometry_info", &fDumpGeometryInfo, kInt, 0, 1, 1},
    { "efficiency_bin_width_1D", &fBinSize_efficiency1D, kDouble, 0, 1, 1 },
    { "efficiency_bin_width_2D", &fBinSize_efficiency2D, kDouble, 0, 1, 1 },
    { "xptar_min", &fxptarmin_track, kDouble, 0, 1, 1},
    { "xptar_max", &fxptarmax_track, kDouble, 0, 1, 1},
    { "yptar_min", &fyptarmin_track, kDouble, 0, 1, 1},
    { "yptar_max", &fyptarmax_track, kDouble, 0, 1, 1},
    { "ytar_min", &fytarmin_track, kDouble, 0, 1, 1},
    { "ytar_max", &fytarmax_track, kDouble, 0, 1, 1},
    { "pmin", &fPmin_track, kDouble, 0, 1, 1},
    { "pmax", &fPmax_track, kDouble, 0, 1, 1},
    { "useslopeconstraint", &useslopeconstraint, kInt, 0, 1, 1 },
    { "xpfp_min", &fxpfpmin, kDouble, 0, 1, 1 },
    { "xpfp_max", &fxpfpmax, kDouble, 0, 1, 1 },
    { "ypfp_min", &fypfpmin, kDouble, 0, 1, 1 },
    { "ypfp_max", &fypfpmax, kDouble, 0, 1, 1 },
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
  
  fIsMC = (mc_flag != 0);
  fTryFastTrack = (fasttrack_flag != 0);

  fUseSlopeConstraint = (useslopeconstraint != 0 );
  
  //std::cout << "pedestal file name = " << fpedfilename << std::endl;
  
  fUseOpticsConstraint = (useopticsconstraint != 0 );

  // std::cout << "pedestal mode flag = " << pedestalmode_flag << std::endl;
  // std::cout << "do efficiency flag = " << doefficiency_flag << std::endl;
  // std::cout << "pedestal mode, efficiency plots = " << fPedestalMode << ", " << fMakeEfficiencyPlots << std::endl;
  
  //fOnlineZeroSuppression = (onlinezerosuppressflag != 0);
  fUseConstraint = (useconstraintflag != 0);
  fMinHitsOnTrack = std::max(3,fMinHitsOnTrack);

  if( fPedestalMode ){ //then we will just dump raw data to the tree and/or histograms:
    //fZeroSuppress = false;
    fMakeEfficiencyPlots = false; //If in pedestal mode, we don't want efficiency plots
  }
  
  //vsplit is a Podd function that "tokenizes" a string into a vector<string> by whitespace:
  std::vector<std::string> modules = vsplit(modconfig);
  if( modules.empty()) {
    Error("", "[SBSGEMSpectrometerTracker::ReadDatabase] No modules defined");
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


Int_t SBSGEMSpectrometerTracker::Begin( THaRunBase* run ) {
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

void SBSGEMSpectrometerTracker::Clear( Option_t *opt ){
  THaTrackingDetector::Clear(opt);
  SBSGEMTrackerBase::Clear();

  for( auto& module: fModules ) {
    module->Clear(opt);
  }
}

Int_t SBSGEMSpectrometerTracker::Decode(const THaEvData& evdata ){
  //return 0;
  //std::cout << "[SBSGEMSpectrometerTracker::Decode], decoding all modules, event ID = " << evdata.GetEvNum() <<  "...";

  //Triggers decoding of each module:

  Int_t stripcounter = 0;
  for( auto& module: fModules ) {

    module->Decode(evdata);

    //std::cout << "Decoding module " << (*it)->GetName() << ", nstrips fired = " << (*it)->fNstrips_hit << std::endl;


    stripcounter += module->fNstrips_hit;
    //std::cout << "done..." << std::endl;
  }

  //std::cout << "done, fNstrips_hit = " << stripcounter << std::endl;
  
  return 0;
}


Int_t SBSGEMSpectrometerTracker::End( THaRunBase* run ){

  UInt_t runnum = run->GetNumber(); 
  
  //To automate the printing out of pedestals for database and DAQ, 
  if( fPedestalMode ){
    TString fname_dbase, fname_daq, fname_cmr, fname_dbcm;
    
    
    TString specname = GetApparatus()->GetName();
    TString detname = GetName();
    fname_dbase.Form( "db_ped_%s_%s_run%d.dat", specname.Data(), detname.Data(), runnum );
    fname_dbcm.Form( "db_cmr_%s_%s_run%d.dat", specname.Data(), detname.Data(), runnum );
    fname_daq.Form( "daq_ped_%s_%s_run%d.dat", specname.Data(), detname.Data(), runnum );
    fname_cmr.Form( "CommonModeRange_%s_%s_run%d.txt", specname.Data(), detname.Data(), runnum );
    
    fpedfile_dbase.open( fname_dbase.Data() );
    fCMfile_dbase.open( fname_dbcm.Data() );
    fpedfile_daq.open( fname_daq.Data() );
    fpedfile_cmr.open( fname_cmr.Data() );
    
    TString sdate = run->GetDate().AsString();
    sdate.Prepend( "#" );
    
    TString message;

    message.Form( "# Copy the contents of this file into $DB_DIR/db_%s.%s.dat to use these pedestals for analysis", specname.Data(), detname.Data() );
    fCMfile_dbase << sdate << std::endl;
    fCMfile_dbase << message << std::endl;
    fCMfile_dbase << "# format = detname.commonmode_mean(U,V) and detname.commonmode_rms(U,V) = common-mode mean and RMS for U, V strips by APV card in order of position"
		  << std::endl;
    
    message.Form( "# Copy the contents of this file into $DB_DIR/db_%s.%s.dat to use these pedestals for analysis", specname.Data(), detname.Data() );
    
    fpedfile_dbase << sdate << std::endl;
    fpedfile_dbase << message << std::endl;
    fpedfile_dbase << "#format = detname.ped(u,v) and detname.rms(u,v) = pedestal mean and rms by (u,v) strip number in order of position" << std::endl;
    
    message.Form("# copy the contents of this file into (location TBD) to update CODA thresholds for detector %s.%s", specname.Data(), detname.Data() );
    
    // fpedfile_daq << sdate << std::endl;
    // fpedfile_daq <<  message << std::endl;
    // fpedfile_daq << "# format = APV        crate       slot       mpd_id       adc_ch followed by " << std::endl
    // 		 << "# APV channel number      pedestal mean      pedestal rms (for average over time samples)" << std::endl;

    message.Form( "# This file defines the common-mode range for the online zero-suppression for the GEM DAQ. Copy its contents into (location TBD) to set these values for detector %s.%s", specname.Data(), detname.Data() );
    
    // fpedfile_cmr << sdate << std::endl;
    // fpedfile_cmr << message << std::endl;
    // fpedfile_cmr << "# format = crate     slot      mpd_id     adc_ch      commonmode min      commonmode max" << std::endl;
  }

  for( auto& module: fModules ) {
    if( fPedestalMode ) { module->PrintPedestals(fpedfile_dbase, fCMfile_dbase, fpedfile_daq, fpedfile_cmr); }
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

    hdidhit_x_layer->Write(0,kOverwrite);
    hdidhit_y_layer->Write(0,kOverwrite);
    hdidhit_xy_layer->Write(0,kOverwrite);

    hshouldhit_x_layer->Write(0,kOverwrite);
    hshouldhit_y_layer->Write(0,kOverwrite);
    hshouldhit_xy_layer->Write(0,kOverwrite);

    hefficiency_x_layer->Write(0,kOverwrite);
    hefficiency_y_layer->Write(0,kOverwrite);
    hefficiency_xy_layer->Write(0,kOverwrite);
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

void SBSGEMSpectrometerTracker::Print(const Option_t* opt) const {
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


void SBSGEMSpectrometerTracker::SetDebug( Int_t level ) {
  THaTrackingDetector::SetDebug(level);
  for( auto& module: fModules ) {
    module->SetDebug(level);
  }
}

Int_t SBSGEMSpectrometerTracker::DefineVariables( EMode mode ){
  if( mode == kDefine and fIsSetup ) return kOK;
  fIsSetup = ( mode == kDefine );
  
  RVarDef vars[] = {
    { "track.ntrack", "number of tracks found", "fNtracks_found" },
    { "track.nhits", "number of hits on track", "fNhitsOnTrack" },
    { "track.x", "Track X (TRANSPORT)", "fXtrack" }, //might be redundant with spectrometer variables, but probably needed for "non-tracking" version
    { "track.y", "Track Y (TRANSPORT)", "fYtrack" },
    { "track.xp", "Track dx/dz (TRANSPORT)", "fXptrack" },
    { "track.yp", "Track dy/dz (TRANSPORT)", "fYptrack" },
    { "track.chi2ndf", "Track Chi2/ndf", "fChi2Track" },
    { "track.besttrack", "Index of 'best' track", "fBestTrackIndex" },
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
    { "hit.ADCavg", "cluster ADC average", "fHitADCavg" },
    { "hit.ADCmaxstripU", "ADC sum of max U strip", "fHitUADCmaxstrip" },
    { "hit.ADCmaxstripV", "ADC sum of max V strip", "fHitVADCmaxstrip" },
    { "hit.ADCmaxsampU", "max sample of max U strip", "fHitUADCmaxsample" },
    { "hit.ADCmaxsampV", "max sample of max V strip", "fHitVADCmaxsample" },
    { "hit.ADCmaxsampUclust", "max U cluster-summed ADC time sample", "fHitUADCmaxclustsample" },
    { "hit.ADCmaxsampVclust", "max V cluster-summed ADC time sample", "fHitVADCmaxclustsample" },
    { "hit.ADCasym", "Hit ADC asymmetry: (ADCU - ADCV)/(ADCU + ADCV)", "fHitADCasym" },
    { "hit.Utime", "cluster timing based on U strips", "fHitUTime" },
    { "hit.Vtime", "cluster timing based on V strips", "fHitVTime" },
    { "hit.UtimeMaxStrip", "cluster timing based on U strips", "fHitUTimeMaxStrip" },
    { "hit.VtimeMaxStrip", "cluster timing based on V strips", "fHitVTimeMaxStrip" },
    { "hit.deltat", "cluster U time - V time", "fHitDeltaT" },
    { "hit.Tavg", "hit T average", "fHitTavg" },
    { "hit.isampmaxUclust", "peak time sample in cluster-summed U ADC samples", "fHitIsampMaxUclust" },
    { "hit.isampmaxVclust", "peak time sample in cluster-summed V ADC samples", "fHitIsampMaxVclust" },
    { "hit.isampmaxUstrip", "peak time sample in max U strip", "fHitIsampMaxUstrip" },
    { "hit.isampmaxVstrip", "peak time sample in max V strip", "fHitIsampMaxVstrip" },
    { "hit.ccor_clust", "correlation coefficient between cluster-summed U and V samples", "fHitCorrCoeffClust" },
    { "hit.ccor_strip", "correlation coefficient between U and V samples on strips with max ADC", "fHitCorrCoeffMaxStrip" },
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
    { "nlayershit", "number of layers with any strip fired", "fNlayers_hit" },
    { "nlayershitu", "number of layers with any U strip fired", "fNlayers_hitU" },
    { "nlayershitv", "number of layers with any V strip fired", "fNlayers_hitV" },
    { "nlayershituv", "number of layers with at least one 2D hit", "fNlayers_hitUV" },
    { "nstripsu_layer", "total number of U strips fired by layer", "fNstripsU_layer" },
    { "nstripsv_layer", "total number of V strips fired by layer", "fNstripsV_layer" },
    { "nclustu_layer", "total number of U clusters by layer", "fNclustU_layer" },
    { "nclustv_layer", "total number of V clusters by layer", "fNclustV_layer" },
    { "n2Dhit_layer", "total_number of 2D hits by layer", "fN2Dhit_layer" },
    { nullptr }
  };
  DefineVarsFromList( vars, mode );

  return 0;
}


Int_t SBSGEMSpectrometerTracker::CoarseTrack( TClonesArray& tracks ){
  
  
  if( !fUseConstraint && !fPedestalMode ){
    //std::cout << "SBSGEMSpectrometerTracker::CoarseTrack...";
    //If no external constraints on the track search region are being used/defined, we do the track-finding in CoarseTrack (before processing all the THaNonTrackingDetectors in the parent spectrometer):
    //std::cout << "calling find_tracks..." << std::endl;
    find_tracks();

    for( int itrack=0; itrack<fNtracks_found; itrack++ ){
      THaTrack *Track = AddTrack( tracks, fXtrack[itrack], fYtrack[itrack], fXptrack[itrack], fYptrack[itrack] );

      int ndf = 2*fNhitsOnTrack[itrack]-4;
      double chi2 = fChi2Track[itrack]*ndf;
      
      Track->SetChi2( chi2, ndf );

      int index = tracks.GetLast();
      Track->SetIndex( index );
    }

    //std::cout << "done. found " << fNtracks_found << " tracks" << std::endl;
    
  }
  
  return 0;
}
Int_t SBSGEMSpectrometerTracker::FineTrack( TClonesArray& tracks ){

  
  if( fUseConstraint && !fPedestalMode ){ //
    //std::cout << "SBSGEMSpectrometerTracker::FineTrack..."; 
    //Calls SBSGEMTrackerBase::find_tracks(), which takes no arguments:
    //std::cout << "calling find_tracks" << std::endl;
    find_tracks();

    //We don't necessarily know 
    
    for( int itrack=0; itrack<fNtracks_found; itrack++ ){
      //AddTrack returns a pointer to the created THaTrack:
      THaTrack *Track = AddTrack( tracks, fXtrack[itrack], fYtrack[itrack], fXptrack[itrack], fYptrack[itrack] );	// Then we can set additional properties of the track using the returned pointer:

      int ndf = 2*fNhitsOnTrack[itrack]-4;
      double chi2 = fChi2Track[itrack]*ndf;
      
      Track->SetChi2( chi2, ndf );

      int index = tracks.GetLast();
      Track->SetIndex( index );
      
    }

    //std::cout << "done. found " << fNtracks_found << " tracks" << std::endl;
    
  }
  
  return 0;
}

bool SBSGEMSpectrometerTracker::PassedOpticsConstraint( TVector3 track_origin, TVector3 track_direction ){
  
  // std::cout << "[SBSGEMSpectrometerTracker::PassedOpticsConstraint]: Checking target reconstruction" 
  // 	    << std::endl;

  double xptemp = track_direction.X()/track_direction.Z();
  double yptemp = track_direction.Y()/track_direction.Z();

  //Project back to z = 0 if appropriate:
  double xtemp = track_origin.X() - xptemp * track_origin.Z(); 
  double ytemp = track_origin.Y() - yptemp * track_origin.Z();
  
  ( (THaTrack*) (*fTestTracks)[0] )->Set( xtemp, ytemp, xptemp, yptemp );

  THaSpectrometer *spec = static_cast<THaSpectrometer *>( GetApparatus() );

  if( spec ){
    spec->FindVertices( *fTestTracks );
  } else {
    return false;
  }

  THaTrack *trtemp = ( (THaTrack*) (*fTestTracks)[0] );

  if( trtemp->HasTarget() && trtemp->HasVertex() ){
    double Ptemp = trtemp->GetP();
    double xptartemp = trtemp->GetTTheta();
    double yptartemp = trtemp->GetTPhi();
    double ytartemp = trtemp->GetTY();

    return ( fxptarmin_track <= xptartemp && xptartemp <= fxptarmax_track &&
	     fyptarmin_track <= yptartemp && yptartemp <= fyptarmax_track &&
	     fytarmin_track <= ytartemp && ytartemp <= fytarmax_track &&
	     fPmin_track <= Ptemp && Ptemp <= fPmax_track );
  } else {
    return false;
  }
}



ClassImp(SBSGEMSpectrometerTracker)
