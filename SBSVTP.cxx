#include "SBSVTP.h"
#include <iostream>

using namespace std;

ClassImp(SBSVTP);

SBSVTP::SBSVTP( const char* name, const char* description, 
  THaApparatus* apparatus) :
  THaDetector(name, description, apparatus)
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
void SBSVTP::DefineVariables( EMode mode )
{
  if( mode == kDefine & fIsSetup ) return kOK;
  fIsSetup = (mode == kDefine );

  RVarDef vars[] = {
    {"vtp.detid",     "VTP detector ID",                    "fVTPClusters.fDet"},
    {"vtp.clus.x",    "VTP clusters x coord",               "fVTPClusters.fX"},
    {"vtp.clus.y",    "VTP clusters y coord",               "fVTPClusters.fY"},
    {"vtp.clus.e",    "VTP clusters energy",                "fVTPClusters.fE"},
    {"vtp.clus.time", "VTP clusters time",                  "fVTPClusters.fTime"},
    {"vtp.clus.size", "VTP clusters size (number of hits)", "fVTPClusters.fSize"},
    {0}
  };

  return DefineVarsFromList( vars, mode );
}

//______________________________________________________________
Int_t SBSVTP::Decode( const THaEvData& evdata )
{
  Int_t fNVTPClusters = 0;
  for(UInt_t i = 0; i < fDetMap->GetSize(); i++) {
    THaDetMap::Module* d = fDetMap->GetModule(i);
    cout << "SBSVTP::Decode crate, slot: " << d->crate << " " << d->slot << endl;
    Decoder::VTPModule* vtp = dynamic_cast<Decoder::VTPModule*>(evdata.GetModule(d->crate, d->slot));

    if(vtp) {
      if(evdataa.GetEvNum() != vtp->GetTriggerNum()) {
        fVTPErrorFlag = 1;
      }
      else {
        // Detector ID HCAL or ECAL
        fVTPClusters.fDet = vtp->GetDetectorID();

        // Cluster information
        fVTPClusters.fX.insert(fVTPClusters.fX.end(), vtp->GetClusterX().begin(), vtp->GetClusterX().end());
        fVTPClusters.fY.insert(fVTPClusters.fY.end(), vtp->GetClusterY().begin(), vtp->GetClusterY().end());
        fVTPClusters.fE.insert(fVTPClusters.fE.end(), vtp->GetClusterE().begin(), vtp->GetClusterE().end());
        fVTPClusters.fTime.insert(fVTPClusters.fTime.end(), vtp->GetClusterTime().begin(), vtp->GetClusterTime().end());
        fVTPClusters.fSize.insert(fVTPClusters.fSize.end(), vtp->GetClusterSize().begin(), vtp->GetClusterSize().end());

        fNVTPClusters++;
      }
    }// found vtp module
  }

  return fNVTPClusters;
}
