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
#include "SBSGEMTrackerBase.h"
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
  const char* const here = "ReadDatabase";
  
  FILE* file = OpenFile( date );
  if( !file ){
    std::cerr << "SBSBigBite::" << here << "(): database not found!"<< std::endl;
    return kFileError;
  }
    
  // Read TRANSPORT matrices
  fXptarMatrixElems.clear();
  fYptarMatrixElems.clear();
  fYtarMatrixElems.clear();
  fPinvMatrixElems.clear();
  fXtarMatrixElems.clear();


  /*
    // Maybe we want to do that instead???? 
    // Don't know, I guess we need to have an example of 
    // optics use in the analyzer to see what is the most convenient
  cout << "reading optics" << endl;
  
  int order = 2;
  ifstream opticsfile("BBoptics.txt");
  if(!opticsfile.is_open()){
    cout << "No optics file, exit..." << endl;
    exit(-1);
  }
  int nparams = 0;
  opticsfile >> nparams;
  int o_i, o_j, o_k, o_l, o_m;
  TVectorD b_xptar(nparams), b_yptar(nparams), b_ytar(nparams), b_pinv(nparams);// xtar "fixed"
  for(int i=0; i<nparams; i++){
    opticsfile >> b_xptar(i);
    opticsfile >> b_yptar(i);
    opticsfile >> b_ytar(i);
    opticsfile >> b_pinv(i);
    opticsfile >> o_m >> o_l >> o_k >> o_j >> o_i;
    if(opticsfile.eof() && i<nparams-1){
      cout << "Optics file shorter than expected (probably not corresponding order), exit..." << endl;
      exit(-1);
    }
  }
  
  ///////////THEN in the code 
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
  
  */
  
  /*
  // Read fOrigin and fSize (currently unused)
  Int_t err = ReadGeometry( file, date );
  if( err ) {
    fclose(file);
    return err;
  }

  //FIXME: move to HRS
  fTMatrixElems.clear();
  fDMatrixElems.clear();
  fPMatrixElems.clear();
  fPTAMatrixElems.clear();
  fYMatrixElems.clear();
  fYTAMatrixElems.clear();
  fLMatrixElems.clear();

  fFPMatrixElems.clear();
  fFPMatrixElems.resize(3);

  map<string,MEdef_t> matrix_map;

  // TRANSPORT to focal-plane tensors
  matrix_map["t"]   = MEdef_t( 3, &fFPMatrixElems, true, 0 );
  matrix_map["y"]   = MEdef_t( 3, &fFPMatrixElems, true, 1 );
  matrix_map["p"]   = MEdef_t( 3, &fFPMatrixElems, true, 2 );
  // Standard focal-plane to target matrix elements (D=delta, T=theta, Y=y, P=phi)
  matrix_map["D"]   = MEdef_t( 3, &fDMatrixElems );
  matrix_map["T"]   = MEdef_t( 3, &fTMatrixElems );
  matrix_map["Y"]   = MEdef_t( 3, &fYMatrixElems );
  matrix_map["P"]   = MEdef_t( 3, &fPMatrixElems );
  // Additional matrix elements describing the dependence of y-target and
  // phi-target on the /absolute value/ of theta, found necessary in optimizing
  // the septum magnet optics (R. Feuerbach, March 1, 2005)
  matrix_map["YTA"] = MEdef_t( 4, &fYTAMatrixElems );
  matrix_map["PTA"] = MEdef_t( 4, &fPTAMatrixElems );
  // Matrix for calculating pathlength from z=0 (target) to focal plane (meters)
  // (R. Feuerbach, October 16, 2003)
  matrix_map["L"]   = MEdef_t( 4, &fLMatrixElems );

  string MEstring, TCmodule;
  DBRequest request1[] = {
    { "matrixelem",  &MEstring, kString },
    { "time_cor",    &TCmodule, kString, 0, true },
    { nullptr }
  };
  err = LoadDB( file, date, request1, fPrefix );
  if( err ) {
    fclose(file);
    return err;
  }
  if( MEstring.empty() ) {
    Error( Here(here), "No matrix elements defined. Set \"maxtrixelem\" in database." );
    fclose(file);
    return kInitError;
  }
  // Parse the matrix elements
  err = ParseMatrixElements( MEstring, matrix_map, fPrefix );
  if( err ) {
    fclose(file);
    return err;
  }
  MEstring.clear();

  // Ensure that we have all three focal plane matrix elements, else we cannot
  // do anything sensible with the tracks
  if( fFPMatrixElems[T000].order == 0 ) {
    Error( Here(here), "Missing FP matrix element t000. Fix database." );
    err = kInitError;
  }
  if( fFPMatrixElems[Y000].order == 0 ) {
    Error( Here(here), "Missing FP matrix element y000. Fix database." );
    err = kInitError;
  }
  if( fFPMatrixElems[P000].order == 0 ) {
    Error( Here(here), "Missing FP matrix element p000. Fix database." );
    err = kInitError;
  }
  if( err ) {
    fclose(file);
    return err;
  }

  // If given, find the module for calculating an event-by-event
  // time offset correction
  if( !TCmodule.empty() ) {
    fTimeCorrectionModule = dynamic_cast<Podd::TimeCorrectionModule*>
      (FindModule(TCmodule.c_str(), "Podd::TimeCorrectionModule", false));
    if( !fTimeCorrectionModule ) {
       Warning( Here(here), "Time correction module \"%s\" not found. "
            "Event-by-event time offsets will NOT be used!\nCheck \"time_cor\" database key",
            TCmodule.c_str() );
    }
  }

  // Compute derived geometry quantities
  fTan_vdc  = fFPMatrixElems[T000].poly[0];
  fVDCAngle = TMath::ATan(fTan_vdc);
  fSin_vdc  = TMath::Sin(fVDCAngle);
  fCos_vdc  = TMath::Cos(fVDCAngle);

  // Define the VDC coordinate axes in the "TRANSPORT system" (z along particle
  // direction at central momentum)
  DefineAxes(fVDCAngle);

  // Read configuration parameters
  fErrorCutoff = 1e9;
  fNumIter = 1;
  fCoordType = kRotatingTransport;
  Int_t disable_tracking = 0, disable_finetrack = 0, only_fastest_hit = 1;
  Int_t do_tdc_hardcut = 1, do_tdc_softcut = 0, ignore_negdrift = 0;
#ifdef MCDATA
  Int_t mc_data = 0;
#endif
  string coord_type;

  DBRequest request[] = {
    { "max_matcherr",      &fErrorCutoff,      kDouble, 0, true },
    { "num_iter",          &fNumIter,          kInt,    0, true },
    { "coord_type",        &coord_type,        kString, 0, true },
    { "disable_tracking",  &disable_tracking,  kInt,    0, true },
    { "disable_finetrack", &disable_finetrack, kInt,    0, true },
    { "only_fastest_hit",  &only_fastest_hit,  kInt,    0, true },
    { "do_tdc_hardcut",    &do_tdc_hardcut,    kInt,    0, true },
    { "do_tdc_softcut",    &do_tdc_softcut,    kInt,    0, true },
    { "ignore_negdrift",   &ignore_negdrift,   kInt,    0, true },
#ifdef MCDATA
    { "MCdata",            &mc_data,           kInt,    0, true },
#endif
    { nullptr }
  };

  err = LoadDB( file, date, request, fPrefix );
  fclose(file);
  if( err )
    return err;

  // Analysis control flags
  SetBit( kOnlyFastest,     only_fastest_hit );
  SetBit( kHardTDCcut,      do_tdc_hardcut );
  SetBit( kSoftTDCcut,      do_tdc_softcut );
  SetBit( kIgnoreNegDrift,  ignore_negdrift );
#ifdef MCDATA
  SetBit( kMCdata,          mc_data );
#endif
  SetBit( kDecodeOnly,      disable_tracking );
  SetBit( kCoarseOnly,      !disable_tracking && disable_finetrack );

  if( !coord_type.empty() ) {
    if( THaString::CmpNoCase(coord_type, "Transport") == 0 )
      fCoordType = kTransport;
    else if( THaString::CmpNoCase(coord_type, "RotatingTransport") == 0 )
      fCoordType = kRotatingTransport;
    else {
      Error( Here(here), "Invalid coordinate type coord_type = %s. "
             "Must be \"Transport\" or \"RotatingTransport\". Fix database.",
             coord_type.c_str() );
      return kInitError;
    }
  }

  // Sanity checks of parameters
  if( fErrorCutoff < 0.0 ) {
    Warning( Here(here), "Negative max_matcherr = %6.2lf makes no sense, "
             "taking absolute.", fErrorCutoff );
    fErrorCutoff = -fErrorCutoff;
  } else if( fErrorCutoff == 0.0 ) {
    Error( Here(here), "Illegal parameter max_matcherr = 0.0. Must be > 0. "
           "Fix database." );
    return kInitError;
  }
  if( fNumIter < 0) {
    Warning( Here(here), "Negative num_iter = %d makes no sense, "
             "taking absolute.", fNumIter );
    fNumIter = -fNumIter;
  } else if( fNumIter > 10 ) {
    Error( Here(here), "Illegal parameter num_iter = %d. Must be <= 10. "
           "Fix database.", fNumIter );
    return kInitError;
  }

  if( fDebug > 0 ) {
#ifdef MCDATA
    Info( Here(here), "VDC flags fastest/hardcut/softcut/noneg/mcdata/"
          "decode/coarse = %d/%d/%d/%d/%d/%d/%d", TestBit(kOnlyFastest),
          TestBit(kHardTDCcut), TestBit(kSoftTDCcut), TestBit(kIgnoreNegDrift),
          TestBit(kMCdata), TestBit(kDecodeOnly), TestBit(kCoarseOnly) );
#else
    Info( Here(here), "VDC flags fastest/hardcut/softcut/noneg/"
          "decode/coarse = %d/%d/%d/%d/%d/%d", TestBit(kOnlyFastest),
          TestBit(kHardTDCcut), TestBit(kSoftTDCcut), TestBit(kIgnoreNegDrift),
          TestBit(kDecodeOnly), TestBit(kCoarseOnly) );
#endif
  }

  // figure out the track length from the origin to the s1 plane
  // since we take the VDC to be the origin of the coordinate
  // space, this is actually pretty simple
  const THaDetector* s1 = nullptr;
  if( GetApparatus() )
    // TODO: need? if so, change to HRS reference detector
    s1 = GetApparatus()->GetDetector("s1");
  if(s1 == nullptr)
    fCentralDist = 0;
  else
    fCentralDist = s1->GetOrigin().Z();

  CalcMatrix(1.,fLMatrixElems); // tensor without explicit polynomial in x_fp
  */
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
  double x_fcp = 0, y_fcp = 0, z_fcp = 0, wx_fcp = 0, wy_fcp = 0;
  double x_bcp = 0, y_bcp = 0, z_bcp = 0, wx_bcp = 0, wy_bcp = 0;
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
	x_bcp+= BBTotalShower->GetShower()->GetX();
	y_bcp+= BBTotalShower->GetShower()->GetY();
	z_bcp+= BBTotalShower->GetShower()->GetOrigin().Z();
	npts++;
	
	wx_bcp+=BBTotalShower->GetShower()->SizeRow()/sqrt(12);
	wy_bcp+=BBTotalShower->GetShower()->SizeCol()/sqrt(12);
      }
      
      if(BBTotalShower->GetShower()->GetNclust()){
	//cout << BBTotalShower->GetPreShower()->GetName() << " " << BBTotalShower->GetPreShower()->GetX() << " " << BBTotalShower->GetPreShower()->GetY() << " " << BBTotalShower->GetPreShower()->GetOrigin().Z() << endl;
	
	Etot+= BBTotalShower->GetPreShower()->GetE();
	x_bcp+= BBTotalShower->GetPreShower()->GetX();
	y_bcp+= BBTotalShower->GetPreShower()->GetY();
	z_bcp+= BBTotalShower->GetPreShower()->GetOrigin().Z();
	npts++;
	
	wx_bcp+=BBTotalShower->GetPreShower()->SizeRow()/sqrt(12);
	wy_bcp+=BBTotalShower->GetPreShower()->SizeCol()/sqrt(12);
      }
      
    }
    
  }
  if(npts){
    x_bcp/=npts;
    y_bcp/=npts;
    z_bcp/=npts;
    
    wx_bcp/=npts;
    wy_bcp/=npts;
    
    std::cout << "Back constraint point x, y, z: " 
     	      << x_bcp << ", " << y_bcp << ", "<< z_bcp 
     	      << "; width x, y: " << wx_bcp << ", " << wy_bcp << endl;
    
    // apply first order optics???
    // Yes, with the electron energy
    //TODO: replace hard-coded coefficients with optics coefficients
    double dx = (x_bcp*(0.522*Etot-0.121)+0.1729*Etot-0.278)/(Etot*2.224-0.249);
    double dy = y_bcp*0.251;
    
    z_fcp = 0;
    x_fcp = x_bcp+dx*(z_fcp-z_bcp);
    y_fcp = y_bcp+dy*(z_fcp-z_bcp);
    
    wx_bcp = wx_fcp;
    wy_bcp = wy_fcp;
    
    std::cout << "Back constraint point x, y, z: " 
     	      << x_fcp << ", " << y_fcp << ", "<< z_fcp 
     	      << "; width x, y: " << wx_fcp << ", " << wy_fcp << endl;
    
    TIter next2( fTrackingDetectors );
    while( auto* theTrackDetector =
	   static_cast<THaTrackingDetector*>( next2() )) {
      if(theTrackDetector->InheritsFrom("SBSGEMTrackerBase")){
	SBSGEMTrackerBase* BBGEM = reinterpret_cast<SBSGEMTrackerBase*>(theTrackDetector);
	BBGEM->SetFrontConstraintPoint(x_fcp, y_fcp, z_fcp);
	BBGEM->SetBackConstraintPoint(x_bcp, y_bcp, z_bcp);
	BBGEM->SetFrontConstraintWidth(wx_fcp, wy_fcp);
	BBGEM->SetBackConstraintWidth(wx_bcp, wy_bcp);
	/*
	BBGEM->SetFrontConstraintPoint(TVector3(x_fcp, y_fcp, z_fcp));
	BBGEM->SetBackConstraintPoint(TVector3(x_bcp, y_bcp, z_bcp));
	BBGEM->SetFrontConstraintWidth(TVector2(wx_fcp, wy_fcp));
	BBGEM->SetBackConstraintWidth(TVector2(wx_bcp, wy_bcp));
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
