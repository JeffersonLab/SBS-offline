//////////////////////////////////////////////////////////////////////////
//
// Bare-bones SBB BigBite spectrometer class 
// Used for testing.
//
//////////////////////////////////////////////////////////////////////////

#include "SBSGEPEArm.h"
#include "TList.h"
#include "SBSHCal.h"
#include "SBSECal.h"
#include "SBSGEMSpectrometerTracker.h"
#include "SBSGEMPolarimeterTracker.h"
#include "THaTrack.h"
#include "SBSRasteredBeam.h"
#include "THaTrackingDetector.h"
#include "TClass.h"

using namespace std;

ClassImp(SBSGEPEArm)


//_____________________________________________________________________________
SBSGEPEArm::SBSGEPEArm( const char* name, const char* description ) :
THaSpectrometer( name, description )
{
  // Constructor. Defines standard detectors

  //fOpticsOrder = -1; //default to zero for optics order. Not sure if this causes a seg fault, but safer.

  SetPID( false );

  
  fECALdist = 4.7; // default 4.7 m (highest-Q2 GEP setting)
  
}

//_____________________________________________________________________________
SBSGEPEArm::~SBSGEPEArm()
{
  // Destructor
}

Int_t SBSGEPEArm::ReadRunDatabase( const TDatime &date ){
  Int_t err = THaSpectrometer::ReadRunDatabase( date );
  if( err ) return err;
  
  FILE* file = OpenRunDBFile( date );
  if( !file ) return kFileError;
  
  //Require magdist:
  const DBRequest req[] = {
    { "ecaldist", &fECALdist, kDouble, 0, 0, 1 },
    { nullptr }
  };
  err = LoadDB( file, date, req );
  fclose(file);
  if( err )
    return kInitError;
  
  //fOpticsOrigin.SetXYZ( 0.0, 0.0, fMagDist + 2.025 );
  
  return kOK;
}

void SBSGEPEArm::Clear( Option_t *opt )
{
  THaSpectrometer::Clear(opt);
  // fFrontConstraintX.clear();
  // fFrontConstraintY.clear();
  // fFrontConstraintZ.clear();
  // fBackConstraintX.clear();
  // fBackConstraintY.clear();
  // fBackConstraintZ.clear();
}


Int_t SBSGEPEArm::ReadDatabase( const TDatime& date )
{

  FILE* file = OpenFile( date );
  if( !file ){
    std::cerr << "SBSGEPEArm::ReadDatabase(): database not found!"<< std::endl;
    return kFileError;
  }
      
  fIsInit = true;
  
  return kOK;
}
  
Int_t SBSGEPEArm::DefineVariables( EMode mode ){
  THaSpectrometer::DefineVariables(mode);
  
  if( mode == kDefine and fIsSetup ) return kOK;
  fIsSetup = ( mode == kDefine );

  RVarDef ecalanglevars[] = {
    { "ECALth_n", "xECAL/ECALdist", "fECALtheta_n" },
    { "ECALph_n", "yECAL/ECALdist", "fECALphi_n" },
    { "ECALdir_x", "x component of ECAL unit vector", "fECALdir_x" },
    { "ECALdir_y", "y component of ECAL unit vector", "fECALdir_y" },
    { "ECALdir_z", "z component of ECAL unit vector", "fECALdir_z" },
    { nullptr }
  };
  DefineVarsFromList( ecalanglevars, mode );
    
  
  return 0;
}

//_____________________________________________________________________________
Int_t SBSGEPEArm::FindVertices( TClonesArray &tracks )
{
  // Eventually this will just project a straight line of the ECAL+CDET "tracks" back to the target;
  // Perhaps once we have the vertex information from tracking, we can do a refinement of the angle reconstruction
  // in Reconstruct();
  
  return 0;
}


//_____________________________________________________________________________
Int_t SBSGEPEArm::TrackCalc()
{
  // Additioal track calculations

  // TODO

  return 0;
}

//_____________________________________________________________________________
Int_t SBSGEPEArm::CoarseTrack()
{
  
  THaSpectrometer::CoarseTrack();
  //This routine will probably never do anything
  // TODO
  //std::cout << " call SBSBigBite::CoarseTrack" << std::endl;
  //std::cout << "done" << std::endl;
  return 0;
}

Int_t SBSGEPEArm::CoarseReconstruct()
{

  // std::cout << "SBSArm::CoarseReconstruct(): polarimeter mode = "
  // 	    << fPolarimeterMode << std::endl;
  
  fECALtheta_n = kBig;
  fECALphi_n = kBig;

  fECALdir_x = kBig;
  fECALdir_y = kBig;
  fECALdir_z = kBig;

  //The following line will cause CoarseProcess to be invoked for ECAL and CDET:
  THaSpectrometer::CoarseReconstruct();

  // Need to add some lines to get the clusters here;
  // For the initial "online" code for testing purposes,
  // we'll just grab the "best" and/or highest-energy cluster and calculate its global position in the Hall
  // coordinate system:

  SBSECal *ECal = nullptr;
  TIter next( fNonTrackingDetectors );

  //Grab a pointer to ECAL:
  
  while( auto* theNonTrackDetector =
	 static_cast<THaNonTrackingDetector*>( next() )) {
    if(theNonTrackDetector->InheritsFrom("SBSECal")){
      ECal = static_cast<SBSECal*>(theNonTrackDetector);
    }
  }

  if( ECal != nullptr && ECal->GetNclust() > 0 ){
    double xclust = ECal->GetX();
    double yclust = ECal->GetY();
    double eclust = ECal->GetE();

    fECALtheta_n = xclust/fECALdist;
    fECALphi_n = yclust/fECALdist;

    TVector3 ECALdir_global;

    TransportToLab( 1.0, fECALtheta_n, fECALphi_n, ECALdir_global );

    fECALdir_x = ECALdir_global.X();
    fECALdir_y = ECALdir_global.Y();
    fECALdir_z = ECALdir_global.Z();

    int itrack = fTracks->GetLast() + 1;

    THaTrack *Ttemp = new( (*fTracks)[itrack] ) THaTrack( xclust, yclust, fECALtheta_n, fECALphi_n );

    //We'll set "Target" the same as "FP" unless and until we
    //develop any corrections for e.g., fringe field curvature
    Ttemp->SetTarget( xclust, yclust, fECALtheta_n, fECALphi_n );
    Ttemp->SetMomentum( eclust );
    Ttemp->SetEnergy( eclust );
    Ttemp->SetPvect( eclust * ECALdir_global );

    //TO-DO: Add CDET
    
  }
  
  return 0;
}

//_____________________________________________________________________________
Int_t SBSGEPEArm::Track()
{

  return 0;  
}

//_____________________________________________________________________________
Int_t SBSGEPEArm::Reconstruct()
{
  // Fine Reconstruction of particles in spectrometer
  //std::cout << "SBSBigBite::Reconstruct()..." << std::endl;

  THaSpectrometer::Reconstruct();

  // This routine is also unlikely to do anything unless we write ANOTHER interstage module to update the electron scattering
  // angle information based on tracking results (for the vertex)
  
  return 0;
  
}

//_____________________________________________________________________________

Int_t SBSGEPEArm::CalcPID(){
  //for now, do nothing, return electron by default
  return 11;
}
//_______________________


//_____________________________________________________________________________

