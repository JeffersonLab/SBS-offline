#include "SBSBPM.h"
static const UInt_t NCHAN = 4;
using namespace std;
//_____________________________________________________________________________
SBSBPM::SBSBPM( const char* name, const char* description,
				  THaApparatus* apparatus ) :
  THaBeamDet(name,description,apparatus),
  fRawSignal(NCHAN),fPedestals(NCHAN),fCorSignal(NCHAN),fRotPos(NCHAN/2),
  fRot2HCSPos(NCHAN/2,NCHAN/2)
{
  // Constructor
}
//_____________________________________________________________________________
SBSBPM::~SBSBPM()
{
  // Destructor. Remove variables from global list.

  if( fIsSetup )
    RemoveVariables();
}
//_____________________________________________________________________________
Int_t SBSBPM::ReadDatabase( const TDatime& date )
{
  // ReadDatabase:  if detectors cant be added to detmap
  //                or entry for bpms is missing           -> kInitError
  //                otherwise                              -> kOk

  const char* const here = "ReadDatabase";

  vector<Int_t> detmap;
  Double_t pedestals[NCHAN], rotations[NCHAN], offsets[2];

  // std::cout << GetDBFileName() << std::endl;

  FILE* file = OpenFile( date );
  if( !file )
    return kInitError;

  fOrigin.SetXYZ(0,0,0);
  Int_t err = ReadGeometry( file, date );

  if( !err ) {
    // Read configuration parameters
    DBRequest config_request[] = {
      { "detmap",    &detmap,  kIntV },
      { 0 }
    };
    err = LoadDB( file, date, config_request, fPrefix );
  }

  UInt_t flags = THaDetMap::kFillLogicalChannel | THaDetMap::kFillModel;
  if( !err && FillDetMap(detmap, flags, here) <= 0 ) {
    err = kInitError;  // Error already printed by FillDetMap
  }

  if( !err ) {
    memset( pedestals, 0, sizeof(pedestals) );
    memset( rotations, 0, sizeof(rotations) );
    memset( offsets  , 0, sizeof( offsets ) );
    DBRequest calib_request[] = {
      { "calib_rot",   &fCalibRot },
      { "pedestals",   pedestals, kDouble, NCHAN, 1 },
      { "rotmatrix",   rotations, kDouble, NCHAN, 1 },
      { "offsets"  ,   offsets,   kDouble, 2    , 1 },
      { 0 }
    };
    err = LoadDB( file, date, calib_request, fPrefix );
  }
  fclose(file);
  if( err )
    return err;

  fOffset(0) = offsets[0];
  fOffset(1) = offsets[1];

  fPedestals.SetElements( pedestals );
  fRot2HCSPos(0,0) = rotations[0];
  fRot2HCSPos(0,1) = rotations[1];
  fRot2HCSPos(1,0) = rotations[2];
  fRot2HCSPos(1,1) = rotations[3];

  // printf(Form("OFFSETS for BPMs in file %d: %f, %f\n",date.GetDate(),offsets[0],offsets[1]));
  return kOK;
}
//_____________________________________________________________________________
Int_t SBSBPM::DefineVariables( EMode mode )
{
  // Initialize global variables and lookup table for decoder


  if( mode == kDefine && fIsSetup ) return kOK;
  fIsSetup = ( mode == kDefine );

  // Register variables in global list

  RVarDef vars[] = {
    { "rawcur.1", "current in antenna 1", "GetRawSignal0()"},
    { "rawcur.2", "current in antenna 2", "GetRawSignal1()"},
    { "rawcur.3", "current in antenna 3", "GetRawSignal2()"},
    { "rawcur.4", "current in antenna 4", "GetRawSignal3()"},
    { "x", "reconstructed x-position", "fPosition.fX"},
    { "y", "reconstructed y-position", "fPosition.fY"},
    { "z", "reconstructed z-position", "fPosition.fZ"},
    { "rotpos1", "position in bpm system","GetRotPosX()"},
    { "rotpos2", "position in bpm system","GetRotPosY()"},
    { 0 }
  };
    

  return DefineVarsFromList( vars, mode );

}
//_____________________________________________________________________________
void SBSBPM::Clear( Option_t* opt )
{
  // Reset per-event data.
  THaBeamDet::Clear(opt);
  fPosition.SetXYZ(0.,0.,-10000.);
  fDirection.SetXYZ(0.,0.,1.);
  fNfired=0;
  for( UInt_t k=0; k<NCHAN; ++k ) {
    fRawSignal(k)=-1;
    fCorSignal(k)=-1;
  }
}
//_____________________________________________________________________________
Int_t SBSBPM::Decode( const THaEvData& evdata )
{
  // loops over all modules defined in the detector map
  // copies raw data into local variables
  // performs pedestal subtraction

  const char* const here = "Decode()";

  // to eliminate compilation warnings 
  Int_t detMapSize = fDetMap->GetSize(); 
  Int_t numCh=0,data=0;
  UInt_t k=0,chan=0; 


  for (Int_t i = 0; i < detMapSize; i++ ){
    THaDetMap::Module* d = fDetMap->GetModule( i );
    numCh = evdata.GetNumChan(d->crate,d->slot);
    for (Int_t j=0; j< numCh; j++) {
      chan = evdata.GetNextChan( d->crate, d->slot, j);
      if ((chan>=d->lo)&&(chan<=d->hi)) {
	data = evdata.GetData( Decoder::kPulseIntegral,d->crate, d->slot, chan, 0 );
	k = d->first + ((d->reverse) ? d->hi - chan : chan - d->lo) -1;
	if ((k<NCHAN)&&(fRawSignal(k)==-1)) {
	  fRawSignal(k)= data;
	  fNfired++;
	}
	else {
	  Warning( Here(here), "Illegal detector channel: %d", k );
	}
      }
    }
  }

  char msg[200]; 
  sprintf(msg,"====> DEBUG: Nfired = %d, NCHAN = %d",(int)fNfired,(int)NCHAN); 

  if (fNfired!=NCHAN) {
    Warning( Here(here), "******* Number of fired Channels out of range. "
	     "Setting beam position to nominal values");
    // THaDetMap::Print() format is: crate, slot, lo, hi, first, model-type, refchan, refindex, resolution, plane, signal
    std::cout << msg << std::endl; 
    std::cout << "SBSBPM::Decode() calling THaDetMap::Print() " << std::endl;
    fDetMap->Print(); 
    std::cout << "SBSBPM: detMapSize = " << detMapSize << std::endl; 
  }
  else {
    fCorSignal=fRawSignal;
    fCorSignal-=fPedestals;
  }

  return 0;
}
//____________________________________________________
Int_t SBSBPM::Process( )
{
 
  // called by the beam apparaturs 
  // uses the pedestal substracted signals from the antennas
  // to get the position in the bpm coordinate system 
  // and uses the transformation matrix defined in the database
  // to transform it into the HCS
  // directions are not calculated, they are always set parallel to z

  Double_t ap, am;

  for( UInt_t k=0; k<NCHAN; k+=2 ) {
    if( fCorSignal(k)+fCorSignal(k+1)!=0.0 ) {
      ap=fCorSignal(k);
      am=fCorSignal(k+1);
      fRotPos(k/2)=fCalibRot*(ap-am)/(ap+am);
    }
    else {
      fRotPos(k/2)=0.0;
    }
  }

  TVectorD dum(fRotPos);

  dum*=fRot2HCSPos;

  // FIXME
  //Double_t targetXYZ[3] = {0.0,0.0,-115.0};
  //printf(Form("Using target position != fOrigin -> = (%f,%f,%f)",targetXYZ[0],targetXYZ[1],targetXYZ[2]));
 
  fPosition.SetXYZ(
		   //dum(0)+targetXYZ[0]+fOffset(0),
		   //dum(1)+targetXYZ[1]+fOffset(1),
		   //targetXYZ[2]
		   dum(0)+fOrigin(0)+fOffset(0),
		   dum(1)+fOrigin(1)+fOffset(1),
		   fOrigin(2)
		   );

  return 0 ;
}
//______________________________________________________________________________
ClassImp(SBSBPM)
