//////////////////////////////////////////////////////////////////////////
//
// Bare-bones SBB BigBite spectrometer class 
// Used for testing.
//
//////////////////////////////////////////////////////////////////////////

#include "SBSBigBite.h"
#include "THaTrack.h"
#include "TList.h"
#include "SBSBBShower.h"
#include "SBSBBTotalShower.h"
#include "SBSGEMSpectrometerTracker.h"
#include "THaTrackingDetector.h"

using namespace std;

ClassImp(SBSBigBite)

//_____________________________________________________________________________
SBSBigBite::SBSBigBite( const char* name, const char* description ) :
  THaSpectrometer( name, description )
{
  // Constructor. Defines standard detectors
  //The standard BigBite detector package in the 12 GeV/SBS era will include:
  // pre-shower + shower calorimeters (inherit from THaNonTrackingDetector OR THaPidDetector)
  // Timing hodoscope (inherit from THaNonTrackingDetector)
  // GRINCH (inherit from THaPidDetector)
  // GEMs (five-layer) (inherit from THaTrackingDetector)

}

//_____________________________________________________________________________
SBSBigBite::~SBSBigBite()
{
  // Destructor
}

//_____________________________________________________________________________
Int_t SBSBigBite::ReadDatabase( const TDatime& date )
{
  // Hack from THaVDC::ReadDatabase()
  const char* const here = "SBSBigBite::ReadDatabase()";
  
  FILE* file = OpenFile( date );
  if( !file ){
    std::cerr << here << "(): database not found!"<< std::endl;
    return kFileError;
  }
  
  UInt_t nparams;
  std::vector<Double_t> optics_param;
  const DBRequest request[] = {
    { "optics_order",    &fOpticsOrder, kUInt,  0, 0, 1},
    { "optics_nelem",     &nparams,      kUInt,   0, 0, 1},
    { "optics_parameters", &optics_param, kDoubleV, 0, 0, 1},
    { "frontconstraintwidth_x", &fFrontConstraintWidthX, kDouble, 0, 0, 0},
    { "frontconstraintwidth_y", &fFrontConstraintWidthY, kDouble, 0, 0, 0},
    { "backconstraintwidth_x", &fBackConstraintWidthX, kDouble, 0, 0, 0},
    { "backconstraintwidth_y", &fBackConstraintWidthY, kDouble, 0, 0, 0},
    {0}
  };
  
  Int_t status = LoadDB( file, date, request, fPrefix, 1 ); //The "1" after fPrefix means search up the tree
  
  if( status != 0 ){
    fclose(file);
    return status;
  }
  
  int n_elem = TMath::FloorNint(optics_param.size()/nparams);
  
  if(n_elem<4){
    std::cerr << "Warning: mismatch between " << optics_param.size()/n_elem
	      << " optics parameters provided and " << nparams
	      << " optics parameters expected!" << std::endl
	      << " Fix database! " << endl;
    return kInitError;
  }
  
  //int o_i, o_j, o_k, o_l, o_m;// shall we use those???
  fb_xptar.resize(nparams);
  fb_yptar.resize(nparams);
  fb_ytar.resize(nparams);
  fb_pinv.resize(nparams);
  
  for(int i=0; i<nparams; i++){
    fb_xptar[i] = optics_param[n_elem*i];
    fb_yptar[i] = optics_param[n_elem*i+1];
    fb_ytar[i] = optics_param[n_elem*i+2];
    fb_pinv[i] = optics_param[n_elem*i+3];
  }
  
  
  ///////////THEN in the code (but where)
  	  // int ipar = 0;
	  // for(int i=0; i<=order; i++){
	  //   for(int j=0; j<=order-i; j++){
	  //     for(int k=0; k<=order-i-j; k++){
	  // 	for(int l=0; l<=order-i-j-k; l++){
	  // 	  for(int m=0; m<=order-i-j-k-l; m++){
	  // 	    double term = pow(xfp_fit,m)*pow(yfp_fit,l)*pow(xpfp_fit,k)*pow(ypfp_fit,j)*pow(xtar,i);
	  // 	    xptar_fit += b_xptar(ipar)*term;
	  // 	    yptar_fit += b_yptar(ipar)*term;
	  // 	    ytar_fit += b_ytar(ipar)*term;
	  // 	    pthetabend_fit += b_pinv(ipar)*term;
	  // 	    //pinv_fit += b_pinv(ipar)*term;
	  // 	    // cout << ipar << " " << term << " " 
	  // 	    //      << b_xptar(ipar) << " " << b_yptar(ipar) << " " 
	  // 	    //      << b_ytar(ipar) << " " << b_pinv(ipar) << endl;
		      
	  // 	    ipar++;
	  // 	  }
	  // 	}
	  //     }
	  //   }
	  // }
  ////////////
  
  fIsInit = true;
  return kOK;
}


//_____________________________________________________________________________
Int_t SBSBigBite::CoarseTrack()
{
  // Coarse track Reconstruction
  THaSpectrometer::CoarseTrack();
  // TODO
  //std::cout << " call SBSBigBite::CoarseTrack" << std::endl;
  return 0;
}

//_____________________________________________________________________________
Int_t SBSBigBite::CoarseReconstruct()
{
  // Coarse Reconstruction of particles in spectrometer
  THaSpectrometer::CoarseReconstruct(); 
  // TODO
  // fetch the clusters from SBSBBShower detectors
  // FOR NOW: fetch the highest clusters from SBSBBShower detectors
  double x_fcp = 0, y_fcp = 0, z_fcp = 0; //, wx_fcp = 0, wy_fcp = 0;
  double x_bcp = 0, y_bcp = 0, z_bcp = 0; //, wx_bcp = 0, wy_bcp = 0;
  double sumweights_x = 0, sumweights_y = 0;
  double Etot = 0;
  int npts = 0;
  TIter next( fNonTrackingDetectors );
  while( auto* theNonTrackDetector =
	 static_cast<THaNonTrackingDetector*>( next() )) {
    //if(theNonTrackDetector->InheritsFrom("SBSBBShower")){
    if(theNonTrackDetector->InheritsFrom("SBSCalorimeter")){
      SBSBBTotalShower* BBTotalShower = reinterpret_cast<SBSBBTotalShower*>(theNonTrackDetector);
      //BBShower->EresMax();
      // gather here the info useful for
      if(BBTotalShower->GetShower()->GetNclust()){
	//cout << BBTotalShower->GetShower()->GetName() << " " << BBTotalShower->GetShower()->GetX() << " " << BBTotalShower->GetShower()->GetY() << " " << BBTotalShower->GetShower()->GetOrigin().Z() << endl;
	
	Etot+= BBTotalShower->GetShower()->GetE();
	x_bcp+= BBTotalShower->GetShower()->GetX()/(BBTotalShower->GetShower()->SizeRow()/sqrt(12));
	y_bcp+= BBTotalShower->GetShower()->GetY()/(BBTotalShower->GetShower()->SizeCol()/sqrt(12));
	z_bcp+= BBTotalShower->GetShower()->GetOrigin().Z();
	npts++;
	sumweights_x+=1./(BBTotalShower->GetShower()->SizeRow()/sqrt(12));
	sumweights_y+=1./(BBTotalShower->GetShower()->SizeCol()/sqrt(12));
	//wx_bcp+=BBTotalShower->GetShower()->SizeRow()/sqrt(12);
	//wy_bcp+=BBTotalShower->GetShower()->SizeCol()/sqrt(12);
      }
      
      if(BBTotalShower->GetPreShower()->GetNclust()){
	//cout << BBTotalShower->GetPreShower()->GetName() << " " << BBTotalShower->GetPreShower()->GetX() << " " << BBTotalShower->GetPreShower()->GetY() << " " << BBTotalShower->GetPreShower()->GetOrigin().Z() << endl;
	
	Etot+= BBTotalShower->GetPreShower()->GetE();
	x_bcp+= BBTotalShower->GetPreShower()->GetX()/(BBTotalShower->GetShower()->SizeRow()/sqrt(12));
	y_bcp+= BBTotalShower->GetPreShower()->GetY()/(BBTotalShower->GetShower()->SizeCol()/sqrt(12));
	z_bcp+= BBTotalShower->GetPreShower()->GetOrigin().Z();
	npts++;
	sumweights_x+=1./(BBTotalShower->GetPreShower()->SizeRow()/sqrt(12));
	sumweights_y+=1./(BBTotalShower->GetPreShower()->SizeCol()/sqrt(12));
	//wx_bcp+=BBTotalShower->GetPreShower()->SizeRow()/sqrt(12);
	//wy_bcp+=BBTotalShower->GetPreShower()->SizeCol()/sqrt(12);
      }
      
    }
    
  }
  if(npts){
    x_bcp/=npts;
    y_bcp/=npts;
    z_bcp/=npts;
    
    //wx_bcp/=npts;
    //wy_bcp/=npts;
    
    // std::cout << "Back constraint point x, y, z: " 
    // 	      << x_bcp << ", " << y_bcp << ", "<< z_bcp 
    //   //<< "; width x, y: " << wx_bcp << ", " << wy_bcp 
    // 	      << endl;
    
    // apply first order optics???
    // Yes, with the electron energy
    //TODO: replace hard-coded coefficients with optics coefficients
    double dx = (x_bcp*(0.522*Etot-0.121)+0.1729*Etot-0.278)/(Etot*2.224-0.249);
    double dy = y_bcp*0.251;
    
    z_fcp = 0;
    x_fcp = x_bcp+dx*(z_fcp-z_bcp);
    y_fcp = y_bcp+dy*(z_fcp-z_bcp);
    
    //wx_fcp = wx_bcp;
    //wy_fcp = wy_bcp;
    
    // std::cout << "Front constraint point x, y, z: " 
    // 	      << x_fcp << ", " << y_fcp << ", "<< z_fcp 
    //   //<< "; width x, y: " << wx_fcp << ", " << wy_fcp 
    // 	      << endl;
    
    TIter next2( fTrackingDetectors );
    while( auto* theTrackDetector =
	   static_cast<THaTrackingDetector*>( next2() )) {
      if(theTrackDetector->InheritsFrom("SBSGEMSpectrometerTracker")){
	SBSGEMSpectrometerTracker* BBGEM = reinterpret_cast<SBSGEMSpectrometerTracker*>(theTrackDetector);
	//std::cout << "setting constraints for tracks" << std::endl;
	BBGEM->SetFrontConstraintPoint(x_fcp, y_fcp, z_fcp);
	BBGEM->SetBackConstraintPoint(x_bcp, y_bcp, z_bcp);
	BBGEM->SetFrontConstraintWidth(fFrontConstraintWidthX, 
				       fFrontConstraintWidthY);
	//(wx_fcp, wy_fcp);
	BBGEM->SetBackConstraintWidth(fBackConstraintWidthX, 
				      fBackConstraintWidthY);
	//(wx_bcp, wy_bcp);
	/*
	BBGEM->SetFrontConstraintPoint(TVector3(x_fcp, y_fcp, z_fcp));
	BBGEM->SetBackConstraintPoint(TVector3(x_bcp, y_bcp, z_bcp));
	BBGEM->SetFrontConstraintWidth(TVector2(fFrontConstraintWidthX, fFrontConstraintWidthY));
	BBGEM->SetBackConstraintWidth(TVector2(fBackConstraintWidthX, fBackConstraintWidthY));
	*/
      }
    }

  }
  //std::cout << " call SBSBigBite::CoarseReconstruct" << std::endl;
  //THaSpectrometer::CoarseReconstruct();
  return 0;
}

//_____________________________________________________________________________
Int_t SBSBigBite::Track()
{
  // Fine track Reconstruction
  THaSpectrometer::Track();
  // TODO

  return 0;
  
}

//_____________________________________________________________________________
Int_t SBSBigBite::Reconstruct()
{
  // Fine Reconstruction of particles in spectrometer
  THaSpectrometer::Reconstruct();
  // TODO

  return 0;
  
}

//_____________________________________________________________________________
  Int_t SBSBigBite::FindVertices( TClonesArray& /* tracks */ )
{
  // Reconstruct target coordinates for all tracks found.

  // TODO

  return 0;
}

//_____________________________________________________________________________
Int_t SBSBigBite::TrackCalc()
{
  // Additioal track calculations

  // TODO

  return 0;
}
