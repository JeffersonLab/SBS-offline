#include <iostream>
#include <string>

#include "THaTrackingDetector.h"
#include "THaRunBase.h"
#include "THaCrateMap.h"
#include "THaAnalysisObject.h"

#include "SBSGEMStand.h"
#include "SBSGEMPlane.h"
#include "Textvars.h"

using namespace Podd;

SBSGEMStand::SBSGEMStand( const char* name, const char* desc, THaApparatus* app ):
    THaTrackingDetector(name,desc,app) {

        fPlanes.clear();
	fIsMC = false;//by default!
        fCrateMap = 0;	
}

SBSGEMStand::~SBSGEMStand(){
    return;
}


THaAnalysisObject::EStatus SBSGEMStand::Init( const TDatime& date ){
    assert( fCrateMap == 0 );
 
    // Why THaTrackingDetector::Init() here? THaTrackingDetector doesn't implement its own Init() method
    THaAnalysisObject::EStatus status = THaTrackingDetector::Init(date);

    if( status == kOK ){
        for (std::vector<SBSGEMPlane *>::iterator it = fPlanes.begin() ; it != fPlanes.end(); ++it){
            status = (*it)->Init(date);
            if( status != kOK ){
                return status;
            }
        }
    } else {
        return kInitError;
    }

    return kOK;
}


Int_t SBSGEMStand::ReadDatabase( const TDatime& date ){
    std::cout << "[Reading SBSGEMStand database]" << std::endl;

    fIsInit = kFALSE;

    FILE* file = OpenFile( date );
    if( !file ) return kFileError;

    //As far as I can tell, this doesn't do anything yet (AJRP):
    Int_t err = ReadGeometry( file, date );
    if( err ) {
        fclose(file);
        return err;
    }

    std::string planeconfig;
    std::vector<Int_t> *cmap = new std::vector<Int_t>;
    //it appears that cmap is not actually used yet in any way. TBD

    DBRequest request[] = {
        { "planeconfig",       &planeconfig,       kString   },
        { "cratemap",          cmap,               kIntV     },
        { "is_mc",             &fIsMC,             kInt, 0, 1},
        {0}
    };

    err = LoadDB( file, date, request, fPrefix );
    fclose(file);
    if(err)
      return err;

    //vsplit is a Podd function that "tokenizes" a string into a vector<string> by whitespace:
    std::vector<std::string> planes = vsplit(planeconfig);
    if( planes.empty()) {
            Error("", "[SBSGEMStand::ReadDatabase] No planes defined");
    }

    for (const auto& name : planes ) {
      fPlanes.push_back(new SBSGEMPlane( name.c_str(), name.c_str(), this, fIsMC));
    }

    fIsInit = kTRUE;
    
    return kOK;
}


Int_t SBSGEMStand::Begin( THaRunBase* run ){
    for (std::vector<SBSGEMPlane *>::iterator it = fPlanes.begin() ; it != fPlanes.end(); ++it){
        (*it)->Begin(run);
    }

    return 0;
}

void SBSGEMStand::Clear( Option_t *opt ){
  THaTrackingDetector::Clear(opt);
  for (auto* plane : fPlanes ){
        plane->Clear(opt);
    }

    return;
}

Int_t SBSGEMStand::Decode(const THaEvData& evdata ){
  //return 0;
  //std::cout << "[SBSGEMStand::Decode]" << std::endl;

    for (auto* plane : fPlanes ) {
      plane->Decode(evdata);
    }

    return 0;
}


Int_t SBSGEMStand::End( THaRunBase* run ){
    for (auto* plane : fPlanes ) {
        plane->End(run);
    }


    return 0;
}

void SBSGEMStand::Print(const Option_t* opt) const {
    std::cout << "GEM Stand " << fName << " with " << fPlanes.size() << " planes defined:" << std::endl;
    /*
    for (std::vector<SBSGEMPlane *>::iterator it = fPlanes.begin() ; it != fPlanes.end(); ++it){
        std::cout << "\t"
        (*it)->Print(opt);
    }
    */
    for( auto* plane : fPlanes.size() ){
        plane->Print(opt);
    }

    return;
 }


void SBSGEMStand::SetDebug( Int_t level ){
      THaTrackingDetector::SetDebug( level );
    for (auto* plane : fPlanes ) {
        plane->SetDebug(level);
    }

    return;
}

Int_t SBSGEMStand::DefineVariables( EMode mode ){
    if( mode == kDefine and fIsSetup ) return kOK;
    fIsSetup = ( mode == kDefine );
    RVarDef vars[] = {
//        { "trkstat", "Track reconstruction status",  "fTrkStat" },
        { 0 },
    };
    DefineVarsFromList( vars, mode );

    return 0;
}


Int_t SBSGEMStand::CoarseTrack( TClonesArray& tracks ){
    return 0;
}
Int_t SBSGEMStand::FineTrack( TClonesArray& tracks ){
    return 0;
}

