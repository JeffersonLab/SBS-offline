#include <iostream>
#include <string>

#include "THaApparatus.h"
#include "THaTrackingDetector.h"
#include "THaRunBase.h"
#include "THaCrateMap.h"
#include "THaAnalysisObject.h"

#include "SBSGEMSpectrometerTracker.h"
#include "SBSGEMModule.h"
#include "THaTrack.h"
#include "TClonesArray.h"

using namespace Podd;

SBSGEMSpectrometerTracker::SBSGEMSpectrometerTracker( const char* name, const char* desc, THaApparatus* app ):
  THaTrackingDetector(name,desc,app), SBSGEMTrackerBase() {

  fModules.clear();
  fIsMC = false;//by default!
  //fCrateMap = 0;
  fPedestalMode = false;
}

SBSGEMSpectrometerTracker::~SBSGEMSpectrometerTracker(){
  return;
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
    
    for (std::vector<SBSGEMModule *>::iterator it = fModules.begin() ; it != fModules.end(); ++it){
      status = (*it)->Init(date); //This triggers calling of ReadDatabase for each module (I hope)!
      if( status != kOK ){
	return status;
      }
    }

    CompleteInitialization();
    
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
  
  int pedestalmode_flag = 0;
  int doefficiency_flag = 1;
  int onlinezerosuppressflag = 0;
  int useconstraintflag = 0; //use constraint on track search region from other detectors in the parent THaSpectrometer (or other
  DBRequest request[] = {
    { "modules",  &modconfig, kString, 0, 0, 1 }, //read the list of modules:
    { "pedfile",  &fpedfilename, kString, 0, 1 },
    { "is_mc",        &fIsMC,    kInt, 0, 1, 1 }, //NOTE: is_mc can also be defined via the constructor in the replay script
    { "onlinezerosuppression", &onlinezerosuppressflag, kInt, 0, 1},
    { "minhitsontrack", &fMinHitsOnTrack, kInt, 0, 1},
    { "maxhitcombos", &fMaxHitCombinations, kInt, 0, 1},
    { "gridbinwidthx", &fGridBinWidthX, kDouble, 0, 1},
    { "gridbinwidthy", &fGridBinWidthY, kDouble, 0, 1},
    { "gridedgetolerancex", &fGridEdgeToleranceX, kDouble, 0, 1},
    { "gridedgetolerancey", &fGridEdgeToleranceY, kDouble, 0, 1},
    { "trackchi2cut", &fTrackChi2Cut, kDouble, 0, 1},
    { "useconstraint", &useconstraintflag, kInt, 0, 1},
    { "sigmahitpos", &fSigma_hitpos, kDouble, 0, 1},
    { "pedestalmode", &pedestalmode_flag, kInt, 0, 1, 1},
    { "do_efficiencies", &doefficiency_flag, kInt, 0, 1, 1},
    {0}
  };

  
  
  Int_t status = kInitError;
  LoadDB( file, date, request, fPrefix );
  fclose(file);

  fMakeEfficiencyPlots = (doefficiency_flag != 0 );
  fPedestalMode = ( pedestalmode_flag != 0 );

  std::cout << "pedestal file name = " << fpedfilename << std::endl;
  
  // std::cout << "pedestal mode flag = " << pedestalmode_flag << std::endl;
  // std::cout << "do efficiency flag = " << doefficiency_flag << std::endl;
  // std::cout << "pedestal mode, efficiency plots = " << fPedestalMode << ", " << fMakeEfficiencyPlots << std::endl;
  
  fOnlineZeroSuppression = (onlinezerosuppressflag != 0);
  fUseConstraint = (useconstraintflag != 0);
  fMinHitsOnTrack = std::max(3,fMinHitsOnTrack);

  if( fPedestalMode ){ //then we will just dump raw data to the tree:
    fZeroSuppress = false;
    fMakeEfficiencyPlots = false; //If in pedestal mode, we don't want efficiency plots
  }
  
  //vsplit is a Podd function that "tokenizes" a string into a vector<string> by whitespace:
  std::vector<std::string> modules = vsplit(modconfig);
  if( modules.empty()) {
    Error("", "[SBSGEMSpectrometerTracker::ReadDatabase] No modules defined");
  }

  int modcounter = 0;

  for (std::vector<std::string>::iterator it = modules.begin() ; it != modules.end(); ++it){
    fModules.push_back(new SBSGEMModule( (*it).c_str(), (*it).c_str(), this) );

    modcounter = fModules.size()-1;
    fModules[modcounter]->fModule = modcounter; //just a dummy index in the module array
  }

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
    
  status = kOK;

  if( status != kOK )
    return status;

  fIsInit = kTRUE;
    
  return kOK;
}


Int_t SBSGEMSpectrometerTracker::Begin( THaRunBase* run ){
  for (std::vector<SBSGEMModule *>::iterator it = fModules.begin() ; it != fModules.end(); ++it){
    (*it)->Begin(run);
  }

  InitEfficiencyHistos(GetName()); //create efficiency histograms (see SBSGEMTrackerBase)
  
  
  return 0;
}

void SBSGEMSpectrometerTracker::Clear( Option_t *opt ){
  for (std::vector<SBSGEMModule *>::iterator it = fModules.begin() ; it != fModules.end(); ++it){
    (*it)->Clear(opt);
  }

  return;
}

Int_t SBSGEMSpectrometerTracker::Decode(const THaEvData& evdata ){
  //return 0;
  std::cout << "[SBSGEMSpectrometerTracker::Decode], decoding all modules, event ID = " << evdata.GetEvNum() <<  std::endl;

  //Triggers decoding of each module:
  
  for (std::vector<SBSGEMModule *>::iterator it = fModules.begin() ; it != fModules.end(); ++it){
    (*it)->Decode(evdata);
  }
  
  return 0;
}


Int_t SBSGEMSpectrometerTracker::End( THaRunBase* run ){

  //To automate the printing out of pedestals for database and DAQ, 
  if( fPedestalMode ){
    TString fname_dbase, fname_daq, fname_cmr;
    
    UInt_t runnum = run->GetNumber(); 
    TString specname = GetApparatus()->GetName();
    TString detname = GetName();
    fname_dbase.Form( "db_ped_%s_%s_run%d.dat", specname.Data(), detname.Data(), runnum );
    fname_daq.Form( "daq_ped_%s_%s_run%d.dat", specname.Data(), detname.Data(), runnum );
    fname_cmr.Form( "CommonModeRange_%s_%s_run%d.txt", specname.Data(), detname.Data(), runnum );
    
    fpedfile_dbase.open( fname_dbase.Data() );
    fpedfile_daq.open( fname_daq.Data() );
    fpedfile_cmr.open( fname_cmr.Data() );
    
    TString sdate = run->GetDate().AsString();
    sdate.Prepend( "#" );
    
    TString message;
    
    message.Form( "# Copy the contents of this file into $DB_DIR/db_%s.%s.dat to use these pedestals for analysis", specname.Data(), detname.Data() );
    
    fpedfile_dbase << sdate << std::endl;
    fpedfile_dbase << message << std::endl;
    fpedfile_dbase << "#format = detname.ped(u,v) and detname.rms(u,v) = pedestal mean and rms by (u,v) strip number in order of position" << std::endl;
    
    message.Form("# copy the contents of this file into (location TBD) to update CODA thresholds for detector %s.%s", specname.Data(), detname.Data() );
    
    fpedfile_daq << sdate << std::endl;
    fpedfile_daq <<  message << std::endl;
    fpedfile_daq << "# format = APV        crate       slot       mpd_id       adc_ch followed by " << std::endl
		 << "# APV channel number      pedestal mean      pedestal rms (for average over time samples)" << std::endl;

    message.Form( "# This file defines the common-mode range for the online zero-suppression for the GEM DAQ. Copy its contents into (location TBD) to set these values for detector %s.%s", specname.Data(), detname.Data() );
    
    fpedfile_cmr << sdate << std::endl;
    fpedfile_cmr << message << std::endl;
    fpedfile_cmr << "# format = crate     slot      mpd_id     adc_ch      commonmode min      commonmode max" << std::endl;
  }
  
  for (std::vector<SBSGEMModule *>::iterator it = fModules.begin() ; it != fModules.end(); ++it){
    if( fPedestalMode ) { (*it)->PrintPedestals( fpedfile_dbase, fpedfile_daq, fpedfile_cmr ); }
    (*it)->End(run);
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

    hdidhit_x_layer->Write();
    hdidhit_y_layer->Write();
    hdidhit_xy_layer->Write();

    hshouldhit_x_layer->Write();
    hshouldhit_y_layer->Write();
    hshouldhit_xy_layer->Write();

    hefficiency_x_layer->Write();
    hefficiency_y_layer->Write();
    hefficiency_xy_layer->Write();
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
  for( unsigned int i = 0; i < fModules.size(); i++ ){
    fModules[i]->Print(opt);
  }

  return;
}


void SBSGEMSpectrometerTracker::SetDebug( Int_t level ){
  THaTrackingDetector::SetDebug( level );
  for (std::vector<SBSGEMModule *>::iterator it = fModules.begin() ; it != fModules.end(); ++it){
    (*it)->SetDebug(level);
  }

  return;
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
    { "hit.ADCasym", "Hit ADC asymmetry: (ADCU - ADCV)/(ADCU + ADCV)", "fHitADCasym" },
    { "hit.Utime", "cluster timing based on U strips", "fHitUTime" },
    { "hit.Vtime", "cluster timing based on V strips", "fHitVTime" },
    { "hit.deltat", "cluster U time - V time", "fHitDeltaT" },
    { "hit.ccor_clust", "correlation coefficient between cluster-summed U and V samples", "fHitCorrCoeffClust" },
    { "hit.ccor_strip", "correlation coefficient between U and V samples on strips with max ADC", "fHitCorrCoeffMaxStrip" },
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

  //std::cout << "SBSGEMSpectrometerTracker::CoarseTrack" << std::endl;
  if( !fUseConstraint && !fPedestalMode ){
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

    //std::cout << "found " << fNtracks_found << " tracks" << std::endl;
    
  }
  
  return 0;
}
Int_t SBSGEMSpectrometerTracker::FineTrack( TClonesArray& tracks ){

  //std::cout << "SBSGEMSpectrometerTracker::FineTrack" << std::endl;
  if( fUseConstraint && !fPedestalMode ){ //

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
    
  }
  
  return 0;
}

ClassImp(SBSGEMSpectrometerTracker)
