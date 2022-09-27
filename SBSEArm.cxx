//////////////////////////////////////////////////////////////////////////
//
// Bare-bones SBB BigBite spectrometer class 
// Used for testing.
//
//////////////////////////////////////////////////////////////////////////

#include "SBSEArm.h"
#include "TList.h"
#include "SBSHCal.h"
#include "SBSGEMSpectrometerTracker.h"

using namespace std;

ClassImp(SBSEArm)


//_____________________________________________________________________________
SBSEArm::SBSEArm( const char* name, const char* description ) :
  THaSpectrometer( name, description )
{
  // Constructor. Defines standard detectors


  fFrontConstraintWidthX = 1.5;
  fFrontConstraintWidthY = 1.5;
  fBackConstraintWidthX = 1.5;
  fBackConstraintWidthY = 1.5;
  fFrontConstraintX0 = 0.0;
  fFrontConstraintY0 = 0.0;
  fBackConstraintX0 = 0.0;
  fBackConstraintY0 = 0.0;
}


Int_t SBSEArm::DefineVariables( EMode mode ){
  THaSpectrometer::DefineVariables(mode);
  
  if( mode == kDefine and fIsSetup ) return kOK;
  fIsSetup = ( mode == kDefine );
  
  
  
  RVarDef constraintvars[] = {
    { "x_fcp", "front track constraint x", "fFrontConstraintX" },
    { "y_fcp", "front track constraint y", "fFrontConstraintY" },
    { "z_fcp", "front track constraint z", "fFrontConstraintZ" },
    { "x_bcp", "back track constraint x", "fBackConstraintX" },
    { "y_bcp", "back track constraint y", "fBackConstraintY" },
    { "z_bcp", "back track constraing z", "fBackConstraintZ" },
    { nullptr }
  };
  DefineVarsFromList( constraintvars, mode );
  
  return 0;
}

//_____________________________________________________________________________
SBSEArm::~SBSEArm()
{
  // Destructor
}


void SBSEArm::Clear( Option_t *opt )
{
  THaSpectrometer::Clear(opt);
  fFrontConstraintX.clear();
  fFrontConstraintY.clear();
  fFrontConstraintZ.clear();
  fBackConstraintX.clear();
  fBackConstraintY.clear();
  fBackConstraintZ.clear();
}


Int_t SBSEArm::ReadDatabase( const TDatime& date )
{

  FILE* file = OpenFile( date );
  if( !file ){
    std::cerr << "SBSEArm::ReadDatabase(): database not found!"<< std::endl;
    return kFileError;
  }
    
  std::vector<Double_t> firstgem_offset;
  
   
  const DBRequest request[] = {
    { "frontconstraintwidth_x", &fFrontConstraintWidthX, kDouble, 0, 1, 0},
    { "frontconstraintwidth_y", &fFrontConstraintWidthY, kDouble, 0, 1, 0},
    { "backconstraintwidth_x", &fBackConstraintWidthX, kDouble, 0, 1, 0},
    { "backconstraintwidth_y", &fBackConstraintWidthY, kDouble, 0, 1, 0},
    { "frontconstraint_x0", &fFrontConstraintX0, kDouble, 0, 1, 0},
    { "frontconstraint_y0", &fFrontConstraintY0, kDouble, 0, 1, 0},
    { "backconstraint_x0", &fBackConstraintX0, kDouble, 0, 1, 0},
    { "backconstraint_y0", &fBackConstraintY0, kDouble, 0, 1, 0},
    { "gemorigin_xyz",    &firstgem_offset, kDoubleV,  0, 1, 1},
    {0}
  };

  Int_t status = LoadDB( file, date, request, fPrefix, 1 ); //The "1" after fPrefix means search up the tree
  fclose(file);
  if( status != 0 ){
    return status;
  }

  fGEMorigin.SetXYZ( firstgem_offset[0],
		     firstgem_offset[1],
		     firstgem_offset[2] );
  


  
  return kOK;
}



//_____________________________________________________________________________
  Int_t SBSEArm::FindVertices( TClonesArray& /* tracks */ )
{
  // Reconstruct target coordinates for all tracks found.

  // TODO

  return 0;
}

//_____________________________________________________________________________
Int_t SBSEArm::TrackCalc()
{
  // Additioal track calculations

  // TODO

  return 0;
}

Int_t SBSEArm::CoarseReconstruct()
{

  THaSpectrometer::CoarseReconstruct(); 

  Double_t x_fcp = 0, y_fcp = 0, z_fcp = 0;
  Double_t x_bcp = 0, y_bcp = 0, z_bcp = 0;
 

  TIter next( fNonTrackingDetectors );
  while( auto* theNonTrackDetector =
 	static_cast<THaNonTrackingDetector*>( next() )) {
   if(theNonTrackDetector->InheritsFrom("SBSHCal")){

     SBSHCal* HCal = reinterpret_cast<SBSHCal*>(theNonTrackDetector);
     
     if(HCal->GetNclust() == 0) return 0;

     std::vector<SBSCalorimeterCluster*> HCalClusters = HCal->GetClusters();

     int i_max = 0;
     double E_max = 0;
     
     for(unsigned int i = 0; i<HCalClusters.size(); i++){
       
       if(HCalClusters[i]->GetE() > E_max){
	 i_max = i;
	 E_max = HCalClusters[i]->GetE();
       }
     }
          

     x_bcp = HCalClusters[i_max]->GetX() + HCal->GetOrigin().X();
     y_bcp = HCalClusters[i_max]->GetY() + HCal->GetOrigin().Y();
     z_bcp = HCal->GetOrigin().Z();
          
     x_fcp = fGEMorigin.X();
     y_fcp = fGEMorigin.Y();
     z_fcp = fGEMorigin.Z();


     fFrontConstraintX.push_back(x_fcp);
     fFrontConstraintY.push_back(y_fcp);
     fFrontConstraintZ.push_back(z_fcp);
     fBackConstraintX.push_back(x_bcp);
     fBackConstraintY.push_back(y_bcp);
     fBackConstraintZ.push_back(z_bcp);

     TIter next2( fTrackingDetectors );
     while( auto* theTrackDetector =
	    static_cast<THaTrackingDetector*>( next2() )) {
       if(theTrackDetector->InheritsFrom("SBSGEMSpectrometerTracker")){
	 SBSGEMSpectrometerTracker* SBSGEM = reinterpret_cast<SBSGEMSpectrometerTracker*>(theTrackDetector);
	 //std::cout << "setting constraints for tracks" << std::endl;
	 SBSGEM->SetFrontConstraintPoint(x_fcp + fFrontConstraintX0, y_fcp + fFrontConstraintY0, z_fcp);
	 SBSGEM->SetBackConstraintPoint(x_bcp + fBackConstraintX0, y_bcp + fBackConstraintY0, z_bcp);
	 SBSGEM->SetFrontConstraintWidth(fFrontConstraintWidthX, 
	 				   fFrontConstraintWidthY);
	 SBSGEM->SetBackConstraintWidth(fBackConstraintWidthX, 
	 			       fBackConstraintWidthY);
	 
       }//End inherits from SBSGEMSpectrometerTracker
     }//End over tracking detectors
     
   }//End inherits from SBSHCal class
 } //End loop over not tracking detectors

 return 0;
}


//_____________________________________________________________________________

