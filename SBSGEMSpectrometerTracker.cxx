#include <iostream>
#include <string>

#include "THaTrackingDetector.h"
#include "THaRunBase.h"
#include "THaCrateMap.h"
#include "THaAnalysisObject.h"

#include "SBSGEMSpectrometerTracker.h"
#include "SBSGEMModule.h"

SBSGEMSpectrometerTracker::SBSGEMSpectrometerTracker( const char* name, const char* desc, THaApparatus* app ):
  THaTrackingDetector(name,desc,app) : SBSGEMTrackerBase() {

        fModules.clear();
	fIsMC = false;//by default!
        fCrateMap = 0;	
}

SBSGEMSpectrometerTracker::~SBSGEMSpectrometerTracker(){
    return;
}


THaAnalysisObject::EStatus SBSGEMSpectrometerTracker::Init( const TDatime& date ){
    assert( fCrateMap == 0 );
 
    // Why THaTrackingDetector::Init() here? THaTrackingDetector doesn't implement its own Init() method
    //Does this trigger the invocation of THaTrackingDetector's dedicated "ReadDatabase" method, perhaps? Even that doesn't really make any sense
    THaAnalysisObject::EStatus status = THaTrackingDetector::Init(date);
    //Note: 

    if( status == kOK ){
        for (std::vector<SBSGEMModule *>::iterator it = fModules.begin() ; it != fModules.end(); ++it){
	  status = (*it)->Init(date); //This triggers calling of ReadDatabase for each module!
            if( status != kOK ){
                return status;
            }
        }
    } else {
        return kInitError;
    }

    return kOK;
}

Int_t 

Int_t SBSGEMSpectrometerTracker::ReadDatabase( const TDatime& date ){
    std::cout << "[Reading SBSGEMSpectrometerTracker database]" << std::endl;

    fIsInit = kFALSE;

    FILE* file = OpenFile( date );
    if( !file ) return kFileError;

    

    std::string modconfig;
    std::vector<Int_t> *cmap = new std::vector<Int_t>;
    //it appears that cmap is not actually used yet in any way. TBD

    DBRequest request[] = {
      { "modules",       &modconfig,       kString   }, //read the list of modules:
      { "is_mc",             &fIsMC,             kInt, 0, 1}, //NOTE: is_mc can also be defined via the constructor in the replay script
      {0}
    };

   
    
    Int_t status = kInitError;
    err = LoadDB( file, date, request, fPrefix );
    fclose(file);

    //vsplit is a Podd function that "tokenizes" a string into a vector<string> by whitespace:
    std::vector<std::string> modules = vsplit(modconfig);
    if( modules.empty()) {
            Error("", "[SBSGEMSpectrometerTracker::ReadDatabase] No modules defined");
    }

    for (std::vector<std::string>::iterator it = planes.begin() ; it != planes.end(); ++it){
      fModules.push_back(new SBSGEMModule( (*it).c_str(), (*it).c_str(), this, fIsMC));
    }

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
    for (std::vector<SBSGEMPlane *>::iterator it = fModules.begin() ; it != fModules.end(); ++it){
        (*it)->Begin(run);
    }

    return 0;
}

void SBSGEMSpectrometerTracker::Clear( Option_t *opt ){
    for (std::vector<SBSGEMPlane *>::iterator it = fModules.begin() ; it != fModules.end(); ++it){
        (*it)->Clear(opt);
    }

    return;
}

Int_t SBSGEMSpectrometerTracker::Decode(const THaEvData& evdata ){
  //return 0;
  //std::cout << "[SBSGEMSpectrometerTracker::Decode]" << std::endl;

  //Triggers decoding of each module:
  
  for (std::vector<SBSGEMPlane *>::iterator it = fModules.begin() ; it != fModules.end(); ++it){
    (*it)->Decode(evdata);
  }
  
  return 0;
}


Int_t SBSGEMSpectrometerTracker::End( THaRunBase* run ){
  for (std::vector<SBSGEMPlane *>::iterator it = fModules.begin() ; it != fModules.end(); ++it){
    (*it)->End(run);
  }


  return 0;
}

void SBSGEMSpectrometerTracker::Print(const Option_t* opt) const {
  std::cout << "GEM Stand " << fName << " with " << fModules.size() << " planes defined:" << std::endl;
  /*
    for (std::vector<SBSGEMPlane *>::iterator it = fModules.begin() ; it != fModules.end(); ++it){
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
  for (std::vector<SBSGEMPlane *>::iterator it = fModules.begin() ; it != fModules.end(); ++it){
    (*it)->SetDebug(level);
  }

  return;
}

Int_t SBSGEMSpectrometerTracker::DefineVariables( EMode mode ){
  if( mode == kDefine and fIsSetup ) return kOK;
  fIsSetup = ( mode == kDefine );
  RVarDef vars[] = {
    //        { "trkstat", "Track reconstruction status",  "fTrkStat" },
    { 0 },
  };
  DefineVarsFromList( vars, mode );

  return 0;
}


Int_t SBSGEMSpectrometerTracker::CoarseTrack( TClonesArray& tracks ){
  return 0;
}
Int_t SBSGEMSpectrometerTracker::FineTrack( TClonesArray& tracks ){
  return 0;
}

