#include "SBSVTP.h"
#include "VTPModule.h"
#include <iostream>

using namespace std;

ClassImp(SBSVTP);

SBSVTP::SBSVTP( const char* name, const char* description, 
  THaApparatus* apparatus) :
  THaNonTrackingDetector(name, description, apparatus)
{
  // constructor
  fVTPErrorFlag = 0;
}

//______________________________________________________________
SBSVTP::~SBSVTP()
{
  DefineVariables( kDelete );
}

//______________________________________________________________
void SBSVTP::Clear( Option_t* opt )
{
  fVTPErrorFlag = 0;
  fVTPClusters.clear();
}

//______________________________________________________________
Int_t SBSVTP::ReadDatabase( const TDatime& date )
{
  const char*  here = "ReadDatabase";
  vector<Int_t> detmap;

  FILE* file = OpenFile( date );
  if( !file ) return kFileError;

  Int_t err = kOK;
  DBRequest config_request[] = {
    { "detmap", &detmap, kIntV },
    { 0 }
  };
  err = LoadDB( file, date, config_request, fPrefix );

  //  UInt_t flags = THaDetMap::kFillLogicalChannel | THaDetMap::kFillModel;
  UInt_t flags = THaDetMap::kSkipLogicalChannel;
  if( !err && FillDetMap(detmap, flags, here) <= 0 ) {
    err = kInitError;
  }

  fclose(file);
  if(err)
    return err;

  return kOK;  
}

//______________________________________________________________
Int_t SBSVTP::DefineVariables( EMode mode )
{

  if( mode == kDefine && fIsSetup ) return kOK;
  fIsSetup = (mode == kDefine );

  RVarDef vars[] = {
    {"errflag",   "VTP Error Flag",                     "fVTPErrorFlag"},
    {"detid",     "VTP detector ID",                    "fVTPClusters.fDet"},
    {"clus.x",    "VTP clusters x coord",               "fVTPClusters.fX"},
    {"clus.y",    "VTP clusters y coord",               "fVTPClusters.fY"},
    {"clus.e",    "VTP clusters energy",                "fVTPClusters.fE"},
    {"clus.time", "VTP clusters time",                  "fVTPClusters.fTime"},
    {"clus.size", "VTP clusters size (number of hits)", "fVTPClusters.fSize"},
    {0}
  };

  return DefineVarsFromList( vars, mode );
}

//______________________________________________________________
Int_t SBSVTP::Decode( const THaEvData& evdata )
{

  Int_t fNVTPClusters = 0;
  Int_t Nvtpfound = 0;
  for(UInt_t i = 0; i < fDetMap->GetSize(); i++) {
    THaDetMap::Module* d = fDetMap->GetModule(i);
    // cout << "SBSVTP::Decode crate, slot: " << d->crate << " " << d->slot << endl;
    Decoder::VTPModule* vtp = dynamic_cast<Decoder::VTPModule*>(evdata.GetModule(d->crate, d->slot));

    if(vtp) {
      Nvtpfound++;
      if(evdata.GetEvNum() != vtp->GetTriggerNum()) {
        fVTPErrorFlag = 1;
      }
      else {
	size_t ncluster = vtp->GetClusterX().size();
	if(Nvtpfound == 1) {
	  // Detector ID HCAL or ECAL
	  fVTPClusters.fDet = vtp->GetDetectorID();
	  // Cluster information
	  fVTPClusters.fX = vtp->GetClusterX();
	  fVTPClusters.fY = vtp->GetClusterY();
	  fVTPClusters.fE = vtp->GetClusterEnergy();
	  fVTPClusters.fTime = vtp->GetClusterTime();
	  fVTPClusters.fSize = vtp->GetClusterSize();
	}
	else {
	  fVTPClusters.fDet = vtp->GetDetectorID();
	  for(auto& x : vtp->GetClusterX()) fVTPClusters.fX.emplace_back(x);
	  for(auto& x : vtp->GetClusterY()) fVTPClusters.fY.emplace_back(x);
	  for(auto& x : vtp->GetClusterEnergy()) fVTPClusters.fE.emplace_back(x);
	  for(auto& x : vtp->GetClusterTime()) fVTPClusters.fTime.emplace_back(x);
	  for(auto& x : vtp->GetClusterSize()) fVTPClusters.fSize.emplace_back(x);
	}

        fNVTPClusters += ncluster;

      }// EvNum and TrigNum match
    }// found vtp module
  }

  return fNVTPClusters;
}

//______________________________________________________________
Int_t SBSVTP::CoarseProcess( TClonesArray& tracks )
{
  return 0;
}

//______________________________________________________________
Int_t SBSVTP::FineProcess( TClonesArray& tracks )
{
  return 0;
}
