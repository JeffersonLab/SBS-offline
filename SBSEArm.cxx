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
#include "THaTrack.h"
#include "SBSRasteredBeam.h"
#include "THaTrackingDetector.h"
#include "TClass.h"

using namespace std;

ClassImp(SBSEArm)


//_____________________________________________________________________________
SBSEArm::SBSEArm( const char* name, const char* description ) :
THaSpectrometer( name, description )
{
  // Constructor. Defines standard detectors

  fOpticsOrder = -1; //default to zero for optics order. Not sure if this causes a seg fault, but safer.

  SetPID( false );

  fFrontConstraintWidthX = 1.5;
  fFrontConstraintWidthY = 1.5;
  fBackConstraintWidthX = 1.5;
  fBackConstraintWidthY = 1.5;
  fFrontConstraintX0 = 0.0;
  fFrontConstraintY0 = 0.0;
  fBackConstraintX0 = 0.0;
  fBackConstraintY0 = 0.0;
  fGEMorigin.SetXYZ(0,0,0);

  fMagDist = 2.25; //This is a required parameter from ReadRunDatabase
  fHCALdist = 17.0; //default 17 m (GEN setting). But mandatory in readrundb
  fBdL = 1.58; // T*m (optional: calculate avg. proton deflection for pcentral)
  
  fOpticsAngle = 0.0;
  fOpticsOrigin.SetXYZ( 0.0, 0.0, fMagDist + 2.025 ); //In g4sbs, the front tracker first plane starts (by default) at 2.025 m downstream of the SBS magnet front face
  InitOpticsAxes( fOpticsAngle );

  fGEMtheta = fOpticsAngle;
  fGEMphi   = 0.0*TMath::DegToRad();
  fGEMorigin = fOpticsOrigin; 

  InitGEMAxes( fGEMtheta, fGEMphi, fGEMorigin );

  fPrecon_flag = 0;

  fFrontConstraintX.clear();
  fFrontConstraintY.clear();
  fFrontConstraintZ.clear();
  fBackConstraintX.clear();
  fBackConstraintY.clear();
  fBackConstraintZ.clear();

  fb_xptar.clear();
  fb_yptar.clear();
  fb_ytar.clear();
  fb_pinv.clear();

  f_oi.clear();
  f_oj.clear();
  f_ok.clear();
  f_ol.clear();
  f_om.clear();

  fb_xfp.clear();
  fb_yfp.clear();
  fb_xpfp.clear();
  fb_ypfp.clear();

  f_foi.clear();
  f_foj.clear();
  f_fok.clear();
  f_fol.clear();
  f_fom.clear();

  

}

//_____________________________________________________________________________
SBSEArm::~SBSEArm()
{
  // Destructor
}

Int_t SBSEArm::ReadRunDatabase( const TDatime &date ){
  Int_t err = THaSpectrometer::ReadRunDatabase( date );
  if( err ) return err;
  
  FILE* file = OpenRunDBFile( date );
  if( !file ) return kFileError;
  
  //Require magdist:
  const DBRequest req[] = {
    { "magdist", &fMagDist, kDouble, 0, 0, 1 },
    { "hcaldist", &fHCALdist, kDouble, 0, 0, 1 },
    { nullptr }
  };
  err = LoadDB( file, date, req );
  fclose(file);
  if( err )
    return kInitError;
  
  fOpticsOrigin.SetXYZ( 0.0, 0.0, fMagDist + 2.025 );
  
  return kOK;
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
  std::vector<Double_t> optics_origin;
  std::vector<Double_t> optics_param;
  
  double gemthetadeg = fGEMtheta * TMath::RadToDeg();
  double gemphideg   = fGEMphi * TMath::RadToDeg();
  double opticsthetadeg = fOpticsAngle * TMath::RadToDeg(); 
  
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
    { "gemtheta", &gemthetadeg, kDouble, 0, 1, 1},
    { "gemphi", &gemphideg, kDouble, 0, 1, 1},
    { "opticstheta", &opticsthetadeg, kDouble, 0, 1, 1},
    { "optics_origin", &optics_origin, kDoubleV, 0, 1, 1},
    { "optics_order",    &fOpticsOrder, kInt,  0, 1, 1},
    { "optics_parameters", &optics_param, kDoubleV, 0, 1, 1},
    { "preconflag", &fPrecon_flag, kUInt, 0, 1, 1 },
    { "BdL", &fBdL, kDouble, 0, 1, 1 },
    {0}
  };

  Int_t status = LoadDB( file, date, request, fPrefix, 1 ); //The "1" after fPrefix means search up the tree
  fclose(file);
  if( status != 0 ){
    return status;
  }

  fOpticsAngle = opticsthetadeg * TMath::DegToRad();
  if( optics_origin.size() == 3 ){ //database overrides default values:
    fOpticsOrigin.SetXYZ( optics_origin[0],
			  optics_origin[1],
			  optics_origin[2] );
  }

  InitOpticsAxes( fOpticsAngle );

  fGEMtheta = gemthetadeg * TMath::DegToRad();
  fGEMphi = gemphideg * TMath::DegToRad();
  
  if( firstgem_offset.size() == 3 ){
    fGEMorigin.SetXYZ( firstgem_offset[0],
		       firstgem_offset[1],
		       firstgem_offset[2] );
  }

  InitGEMAxes( fGEMtheta, fGEMphi );

  //Optics model initialization (copy of BigBite for now):
  if(fOpticsOrder>=0){
    unsigned int nterms = 0;
    
    for(int i = 0; i<=fOpticsOrder; i++){ //x 
      for( int j=0; j<=fOpticsOrder-i; j++){ //y
	for( int k=0; k<=fOpticsOrder-i-j; k++){
	  for( int l=0; l<=fOpticsOrder-i-j-k; l++){
	    for( int m=0; m<=fOpticsOrder-i-j-k-l; m++ ){
	      nterms++;
	    }
	  }
	}
      }
    }
    cout << nterms << " lines of parameters expected for optics of order " << fOpticsOrder << endl;
    
    //int n_elem = TMath::FloorNint(optics_param.size()/nparam);
    
    //we expect 9 parameters per line: four coefficients plus five exponents:

    unsigned int nparams = 9*nterms;

    if(nparams!=optics_param.size()){
      std::cerr << "Warning: mismatch between " << optics_param.size()
		<< " optics parameters provided and " << nparams*9
		<< " optics parameters expected!" << std::endl;
      std::cerr << " Fix database! " << std::endl;
      return kInitError;
    }
    
    //int o_i, o_j, o_k, o_l, o_m;// shall we use those???
    fb_xptar.resize(nterms);
    fb_yptar.resize(nterms);
    fb_ytar.resize(nterms);
    fb_pinv.resize(nterms);
    f_oi.resize(nterms);
    f_oj.resize(nterms);
    f_ok.resize(nterms);
    f_ol.resize(nterms);
    f_om.resize(nterms);
    
    for(unsigned int i=0; i<nterms; i++){
      fb_xptar[i] = optics_param[9*i];
      fb_yptar[i] = optics_param[9*i+1];
      fb_ytar[i] = optics_param[9*i+2];
      fb_pinv[i] = optics_param[9*i+3];
      f_om[i] = int(optics_param[9*i+4]);
      f_ol[i] = int(optics_param[9*i+5]);
      f_ok[i] = int(optics_param[9*i+6]);
      f_oj[i] = int(optics_param[9*i+7]);
      f_oi[i] = int(optics_param[9*i+8]);
    }
  }
  
  fIsInit = true;
  
  return kOK;
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
    { "z_bcp", "back track constraint z", "fBackConstraintZ" },
    { nullptr }
  };
  DefineVarsFromList( constraintvars, mode );

  RVarDef hcalanglevars[] = {
    { "HCALth_n", "xHCAL/HCALdist", "fHCALtheta_n" },
    { "HCALph_n", "yHCAL/HCALdist", "fHCALphi_n" },
    { "HCALdir_x", "x component of HCAL unit vector", "fHCALdir_x" },
    { "HCALdir_y", "y component of HCAL unit vector", "fHCALdir_y" },
    { "HCALdir_z", "z component of HCAL unit vector", "fHCALdir_z" },
    { nullptr }
  };
  DefineVarsFromList( hcalanglevars, mode );
    
  
  return 0;
}

//_____________________________________________________________________________
Int_t SBSEArm::FindVertices( TClonesArray &tracks )
{
  // Reconstruct target coordinates for all tracks found.

  // TODO
  //std::cout << "SBSBigBite::FindVertices()...";
  // Reconstruct target coordinates for all tracks found.
  Int_t n_trk = tracks.GetLast()+1;
  for( Int_t t = 0; t < n_trk; t++ ) {
    auto* theTrack = static_cast<THaTrack*>( tracks.At(t) );
    CalcOpticsCoords(theTrack);
    
    if(fOpticsOrder>=0)CalcTargetCoords(theTrack);
  }

  if( GetNTracks() > 0 ) {
    // Select first track in the array. If there is more than one track
    // and track sorting is enabled, then this is the best fit track
    // (smallest chi2/ndof).  Otherwise, it is the track with the best
    // geometrical match (smallest residuals) between the U/V clusters
    // in the upper and lower VDCs (old behavior).
    // 
    // Chi2/dof is a well-defined quantity, and the track selected in this
    // way is immediately physically meaningful. The geometrical match
    // criterion is mathematically less well defined and not usually used
    // in track reconstruction. Hence, chi2 sorting is preferable, albeit
    // obviously slower.
    
    fGoldenTrack = static_cast<THaTrack*>( fTracks->At(0) );
    fTrkIfo      = *fGoldenTrack;
    fTrk         = fGoldenTrack;
  } else {
    fGoldenTrack = nullptr;  
  }
  
  return 0;
}

void SBSEArm::CalcOpticsCoords( THaTrack* track )
{
  Double_t x_fp, y_fp, xp_fp, yp_fp;
  //Double_t px, py, pz;// NB: not the actual momentum!
  
  x_fp = track->GetX();
  y_fp = track->GetY();
  xp_fp = track->GetTheta();
  yp_fp = track->GetPhi();

  TVector3 TrackPosLocal_GEM( x_fp, y_fp, 0.0 );

  //std::cout << "calculating optics coordinates: (xfp,yfp,xpfp,ypfp)=(" << x_fp << ", " << y_fp << ", " << xp_fp << ", " << yp_fp << ")" << std::endl;

  TVector3 TrackPosGlobal_GEM = fGEMorigin + TrackPosLocal_GEM.X() * fGEMxaxis_global + TrackPosLocal_GEM.Y() * fGEMyaxis_global + TrackPosLocal_GEM.Z() * fGEMzaxis_global;
  
  //std::cout << "Track pos global = " << endl;
  //TrackPosGlobal_GEM.Print(); 

  TVector3 TrackDirLocal_GEM( xp_fp, yp_fp, 1.0 );
  TrackDirLocal_GEM = TrackDirLocal_GEM.Unit();

  //  std::cout << "Track direction local = " << endl;

  // TrackDirLocal_GEM.Print();

  TVector3 TrackDirGlobal_GEM = TrackDirLocal_GEM.X() * fGEMxaxis_global + TrackDirLocal_GEM.Y() * fGEMyaxis_global + TrackDirLocal_GEM.Z() * fGEMzaxis_global;
  TrackDirGlobal_GEM = TrackDirGlobal_GEM.Unit(); //Likely unnecessary, but safer (I guess)
  
  //  std::cout << "Track direction global = " << endl;
  //TrackDirGlobal_GEM.Print();

  //Now project track to the z = 0 plane of the ideal optics system:
  // recall the formula to intersect a ray with a plane:
  // (x + s * trackdir - x0) dot planedir = 0.0

  double sintersect = (fOpticsOrigin - TrackPosGlobal_GEM).Dot(fOpticsZaxis_global)/ (TrackDirGlobal_GEM.Dot( fOpticsZaxis_global ) );

  TVector3 TrackIntersect_global = TrackPosGlobal_GEM + sintersect * TrackDirGlobal_GEM;
  
  //  std::cout << "Track intersection point, global = " << endl;
  //TrackIntersect_global.Print();

  //rather than modifying the x, y, theta, phi directly, let's use the RX, RY, RTheta, and RPhi coordinates:
  //TVector3 TrackIntersect_ = TrackIntersect_global - fOpticsOrigin;

  double xoptics = (TrackIntersect_global - fOpticsOrigin).Dot( fOpticsXaxis_global );
  double yoptics = (TrackIntersect_global - fOpticsOrigin).Dot( fOpticsYaxis_global );

  double xpoptics = TrackDirGlobal_GEM.Dot( fOpticsXaxis_global )/TrackDirGlobal_GEM.Dot( fOpticsZaxis_global );
  double ypoptics = TrackDirGlobal_GEM.Dot( fOpticsYaxis_global )/TrackDirGlobal_GEM.Dot( fOpticsZaxis_global );

  //std::cout << "GEM origin = " << std::endl;
  //fGEMorigin.Print();
  //std::cout << "Optics origin = " << std::endl;
  //fOpticsOrigin.Print();

  //std::cout << "GEM z axis global = " << std::endl;
  //fGEMzaxis_global.Print();
  //std::cout << "Optics z axis global = " << std::endl;
  //fOpticsZaxis_global.Print();
  
  //std::cout << "GEM (x,y,xp,yp) = " << x_fp << ", " << y_fp << ", " << xp_fp << ", " << yp_fp << std::endl;
  // std::cout << "Optics (x,y,xp,yp) = " << xoptics << ", " << yoptics << ", " << xpoptics << ", " << ypoptics << endl;
  
  track->SetR( xoptics, yoptics, xpoptics, ypoptics );
  
  //The following line is no longer necessary as we are using the "Rotated TRANSPORT coordinates" to store the track parameters in ideal optics system
  //track->Set(x_fp, y_fp, xp_fp, yp_fp);
  
  
}

void SBSEArm::CalcTargetCoords( THaTrack* track )
{
  //std::cout << "SBSBigBite::CalcTargetCoords()...";
  
  //const double //make it configurable
  const double th_sbs = GetThetaGeo();//retrieve the actual angle

  //th_sbs < 0 for beam right... this makes SBS_zaxis along -x  
  
  TVector3 SBS_zaxis( sin(th_sbs), 0, cos(th_sbs) );
  TVector3 SBS_xaxis(0,-1,0);
  TVector3 SBS_yaxis = (SBS_zaxis.Cross(SBS_xaxis)).Unit();

  TVector3 spec_xaxis_fp,spec_yaxis_fp, spec_zaxis_fp;
  TVector3 spec_xaxis_tgt,spec_yaxis_tgt, spec_zaxis_tgt;
  
  spec_xaxis_tgt = SBS_xaxis;
  spec_yaxis_tgt = SBS_yaxis;
  spec_zaxis_tgt = SBS_zaxis;
  
  spec_zaxis_fp = SBS_zaxis;
  spec_yaxis_fp = SBS_yaxis;
  spec_zaxis_fp.Rotate(-fOpticsAngle, spec_yaxis_fp);
  spec_xaxis_fp = spec_yaxis_fp.Cross(spec_zaxis_fp).Unit();
 
  Double_t x_fp, y_fp, xp_fp, yp_fp;
  //if( fCoordType == kTransport ) {

  if( track->HasRot() ){
    //    std::cout << "using rotated track coordinates for optics: " << endl;
    x_fp = track->GetRX();
    y_fp = track->GetRY();
    xp_fp = track->GetRTheta();
    yp_fp = track->GetRPhi();
  } else {
    //std::cout << "using non-rotated track coordinates for optics: " << endl;
    x_fp = track->GetX();
    y_fp = track->GetY();
    xp_fp = track->GetTheta();
    yp_fp = track->GetPhi();
  }
  //}
  //cout << x_fp << " " << y_fp << " " << xp_fp << " " << yp_fp << endl;

  double /*vx, vy, vz, */px, py, pz;
  double p_fit, xptar_fit, yptar_fit, ytar_fit, xtar;
  xtar = 0.0;
  double thetabend_fit;
  double pthetabend_fit;
  double vz_fit;

  xptar_fit = 0.0;
  yptar_fit = 0.0;
  ytar_fit = 0.0;
  pthetabend_fit = 0.0;
  
  int ipar = 0;
  for(int i=0; i<=fOpticsOrder; i++){
    for(int j=0; j<=fOpticsOrder-i; j++){
      for(int k=0; k<=fOpticsOrder-i-j; k++){
	for(int l=0; l<=fOpticsOrder-i-j-k; l++){
	  for(int m=0; m<=fOpticsOrder-i-j-k-l; m++){
	    double term = pow(x_fp,m)*pow(y_fp,l)*pow(xp_fp,k)*pow(yp_fp,j)*pow(xtar,i);
	    xptar_fit += fb_xptar[ipar]*term;
	    yptar_fit += fb_yptar[ipar]*term;
	    ytar_fit += fb_ytar[ipar]*term;
	    pthetabend_fit += fb_pinv[ipar]*term;
	    //pinv_fit += b_pinv(ipar)*term;
	    // cout << ipar << " " << term << " " 
	    //      << b_xptar(ipar) << " " << b_yptar(ipar) << " " 
	    //      << b_ytar(ipar) << " " << b_pinv(ipar) << endl;
	    ipar++;
	  }
	}
      }
    }
  }

  TVector3 phat_tgt_fit(xptar_fit, yptar_fit, 1.0 );
  phat_tgt_fit = phat_tgt_fit.Unit();

  TVector3 phat_fp(xp_fp, yp_fp, 1.0 );
  phat_fp = phat_fp.Unit();

  TVector3 phat_fp_rot = phat_fp.X() * fOpticsXaxis_global + 
    phat_fp.Y() * fOpticsYaxis_global + 
    phat_fp.Z() * fOpticsZaxis_global; 
  
  thetabend_fit = acos( phat_fp_rot.Dot( phat_tgt_fit ) );

  // TVector3 phat_tgt_fit_global = phat_tgt_fit.X() * spec_xaxis_tgt +
  //   phat_tgt_fit.Y() * spec_yaxis_tgt +
  //   phat_tgt_fit.Z() * spec_zaxis_tgt;
  
  // TVector3 phat_fp_fit(xp_fp, yp_fp, 1.0 );
  // phat_fp_fit = phat_fp_fit.Unit();
  
  // TVector3 phat_fp_fit_global = phat_fp_fit.X() * spec_xaxis_fp +
  //   phat_fp_fit.Y() * spec_yaxis_fp +
  //   phat_fp_fit.Z() * spec_zaxis_fp;
  
  //thetabend_fit = acos( phat_fp_fit_global.Dot( phat_tgt_fit_global ) );

  //if( fPrecon_flag != 1 ){
  p_fit = pthetabend_fit/thetabend_fit;
    //} else {
    //double delta = pthetabend_fit;
    //double p_firstorder = fA_pth1 * ( 1.0 + (fB_pth1 + fC_pth1*fMagDist)*xptar_fit ) / thetabend_fit;
    //p_fit = p_firstorder * (1.0 + delta);
    //}

  //For SBS, which is on beam right, we have ytar = vz * sin(|theta|) - yptar * vz * cos(theta)
  // i.e., ytar = vz * ( sin(|theta|) - yptar * cos(theta) )
  // --> vz = ytar/(sin(|theta|) - yptar * cos(theta) )
  // But this is the same thing as -ytar/(sin(-|theta|) + yptar * cos(theta))
  // So ASSUMING th_sbs < 0 for beam right, the formula below can be used unchanged: 
  vz_fit = -ytar_fit / (sin(th_sbs) + cos(th_sbs)*yptar_fit);
  
  pz = p_fit*sqrt( 1.0/(xptar_fit*xptar_fit+yptar_fit*yptar_fit+1.) );
  px = xptar_fit * pz;
  py = yptar_fit * pz;

  TVector3 pvect_SBS = TVector3(px, py, pz);

  //In SBS transport coordinates, py is toward small angles, px is vertically down (toward the floor), pz is along spectrometer central axis. To translate to HALL coordinates, we have:
  // x to beam left, y vertically up, and z along the beam direction:

  //since th_sbs < 0 for beam right, the formula below can be used unchanged (I think):
  px = +pvect_SBS.Z()*sin(th_sbs)+pvect_SBS.Y()*cos(th_sbs);
  py = -pvect_SBS.X();
  pz = pvect_SBS.Z()*cos(th_sbs)-pvect_SBS.Y()*sin(th_sbs);

  //We should move this to before the optics calculations in case we want to actually correct the optics for the beam position:
  double ybeam=0.0,xbeam=0.0;
  
  //retrieve beam position, if available, to calculate xtar.
  TIter aiter(gHaApps);
  THaApparatus* app = 0;
  while( (app=(THaApparatus*)aiter()) ){
    if(app->InheritsFrom("SBSRasteredBeam")){
      SBSRasteredBeam* RasterBeam = reinterpret_cast<SBSRasteredBeam*>(app);
      //double xbeam = RasterBeam->GetPosition().X();
      ybeam = RasterBeam->GetPosition().Y()/1000.0; //if this is given in mm, we need to convert to meters (also for BB)
      xbeam = RasterBeam->GetPosition().X()/1000.0;
      xtar = - ybeam - cos(GetThetaGeo()) * vz_fit * xptar_fit;
    }
    //cout << var->GetName() << endl;
  }
  //  f_xtg_exp.push_back(xtar);
  
  track->SetTarget(xtar, ytar_fit, xptar_fit, yptar_fit);
  track->SetMomentum(p_fit);
  track->SetPvect(TVector3(px, py, pz));
  track->SetVertex(TVector3(xbeam, ybeam, vz_fit));

  //cout << px << " " << py << " " << pz << "   " << vz_fit << endl;
  //cout << track->GetLabPx() << " " << track->GetLabPy() << " " << track->GetLabPz() 
  //   << "   " << track->GetVertexZ() << endl;
  
  //std::cout << "Done." << std::endl;
}

//_____________________________________________________________________________
Int_t SBSEArm::TrackCalc()
{
  // Additioal track calculations

  // TODO

  return 0;
}

//_____________________________________________________________________________
Int_t SBSEArm::CoarseTrack()
{
  // Coarse track Reconstruction
  // std::cout << " SBSEArm::CoarseTrack()..." << std::endl;

  // if( !fTrackingDetectors ) {
  //   cerr << "fTrackingDetectors == NULL?" << endl;
  //   exit(1);
  // }
  // cout << fTrackingDetectors->IsA()->GetName() << endl;
  
  // fTrackingDetectors->Print();
  
  //  std::cout << " fTrackingDetectors->Print() invoked..." << std::endl;

  THaSpectrometer::CoarseTrack();
  // TODO
  //std::cout << " call SBSBigBite::CoarseTrack" << std::endl;
  //std::cout << "done" << std::endl;
  return 0;
}

Int_t SBSEArm::CoarseReconstruct()
{

  fHCALtheta_n = kBig;
  fHCALphi_n = kBig;

  fHCALdir_x = kBig;
  fHCALdir_y = kBig;
  fHCALdir_z = kBig;

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
          
      fHCALtheta_n = HCalClusters[i_max]->GetX()/fHCALdist;
      fHCALphi_n = HCalClusters[i_max]->GetY()/fHCALdist;

      TVector3 HCALdir_global; 

      TransportToLab( 1.0, fHCALtheta_n, fHCALphi_n, HCALdir_global );

      fHCALdir_x = HCALdir_global.X();
      fHCALdir_y = HCALdir_global.Y();
      fHCALdir_z = HCALdir_global.Z();

      //x_fcp = fGEMorigin.X();
      //y_fcp = fGEMorigin.Y();
      //z_fcp = fGEMorigin.Z();

      x_fcp = 0.0;
      y_fcp = 0.0;
      z_fcp = 0.0;

      fFrontConstraintX.push_back( x_fcp );
      fFrontConstraintY.push_back( y_fcp );
      fFrontConstraintZ.push_back( z_fcp );

      fBackConstraintX.push_back( x_bcp );
      fBackConstraintY.push_back( y_bcp );
      fBackConstraintZ.push_back( z_bcp );


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
Int_t SBSEArm::Track()
{
  // Fine track Reconstruction
  THaSpectrometer::Track();
  // TODO

  return 0;
  
}

//_____________________________________________________________________________
Int_t SBSEArm::Reconstruct()
{
  // Fine Reconstruction of particles in spectrometer
  //std::cout << "SBSBigBite::Reconstruct()..." << std::endl;
  
  THaSpectrometer::Reconstruct();
  // TODO

  //std::cout << "Done..." << std::endl;
  
  return 0;
  
}

//_____________________________________________________________________________

Int_t SBSEArm::CalcPID(){
  //for now, do nothing
  return 0;
}
//_______________________

void SBSEArm::InitOpticsAxes(double BendAngle, const TVector3 &Origin ){
  fOpticsOrigin = Origin;
  fOpticsYaxis_global.SetXYZ(0,1,0);
  fOpticsZaxis_global.SetXYZ(-sin(BendAngle), 0, cos(BendAngle) );
  fOpticsXaxis_global.SetXYZ(cos(BendAngle), 0, sin(BendAngle) );
}

void SBSEArm::InitOpticsAxes(double BendAngle ){
  // fOpticsOrigin = Origin;
  fOpticsYaxis_global.SetXYZ(0,1,0);
  fOpticsZaxis_global.SetXYZ(-sin(BendAngle), 0, cos(BendAngle) );
  fOpticsXaxis_global.SetXYZ(cos(BendAngle), 0, sin(BendAngle) );
}

void SBSEArm::InitGEMAxes( double theta, double phi, const TVector3 &Origin ){
  fGEMorigin = Origin;
  fGEMzaxis_global.SetXYZ( sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta) );
  fGEMxaxis_global = (fOpticsYaxis_global.Cross( fGEMzaxis_global) ).Unit(); // check to make sure this is consistent with definition in the zero-field alignment code
  fGEMyaxis_global = (fGEMzaxis_global.Cross(fGEMxaxis_global)).Unit();
}

void SBSEArm::InitGEMAxes( double theta, double phi ){
  //fGEMorigin = Origin;
  fGEMzaxis_global.SetXYZ( sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta) );
  fGEMxaxis_global = (fOpticsYaxis_global.Cross( fGEMzaxis_global) ).Unit(); // check to make sure this is consistent with definition in the zero-field alignment code
  fGEMyaxis_global = (fGEMzaxis_global.Cross(fGEMxaxis_global)).Unit();
}

//_____________________________________________________________________________

