//////////////////////////////////////////////////////////////////////////
//
// Bare-bones SBB BigBite spectrometer class 
// Used for testing.
//
//////////////////////////////////////////////////////////////////////////

#include "SBSBigBite.h"
#include "THaTrack.h"
#include "THaPIDinfo.h"
#include "TList.h"
#include "SBSBBShower.h"
#include "SBSBBTotalShower.h"
#include "SBSTimingHodoscope.h"
#include "SBSGRINCH.h"
#include "SBSGEMSpectrometerTracker.h"
#include "THaTrackingDetector.h"
#include "TH2D.h"

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
  /*
  h1_yVx_bcp = new TH2D("h1_yVx_bcp", ";x_{bcp} (m);y_{bcp} (m)", 300, -1.5, 1.5, 100, -0.5, 0.5);
  h1_x_fcpVbcp = new TH2D("h1_x_fcpVbcp", ";x_{bcp} (m);x_{fcp} (m)", 300, -1.5, 1.5, 300, -1.5, 1.5);
  h1_yVx_fcp = new TH2D("h1_yVx_fcp", ";x_{fcp} (m);y_{fcp} (m)", 300, -1.5, 1.5, 100, -0.5, 0.5);
  h1_y_fcpVbcp = new TH2D("h1_y_fcpVbcp", ";y_{bcp} (m);y_{fcp} (m)", 100, -0.5, 0.5, 100, -0.5, 0.5);
  h1_dyVdx = new TH2D("h1_dyVdx",";dx/dz;dy/dz", 100, -0.5, 0.5, 50, -0.25, 0.25);
  */
}

//_____________________________________________________________________________
SBSBigBite::~SBSBigBite()
{
  // Destructor
}

void SBSBigBite::Clear( Option_t *opt )
{
  THaSpectrometer::Clear(opt);
  fEpsEtotRatio.clear();
  fFrontConstraintX.clear();
  fFrontConstraintY.clear();
  fBackConstraintX.clear();
  fBackConstraintY.clear();
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
    { "trackgrinchcorr_const", &fTrackGrinchClusCorr_0, kDouble, 0, 1, 0},
    { "trackgrinchcorr_slope", &fTrackGrinchClusCorr_1, kDouble, 0, 1, 0},
    { "trackgrinchcorr_sigma", &fTrackGrinchClusCorr_Sigma, kDouble, 0, 1, 0},
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
  
  fIsInit = true;
  return kOK;
}

Int_t SBSBigBite::DefineVariables( EMode mode ){
  THaSpectrometer::DefineVariables(mode);
  
  if( mode == kDefine and fIsSetup ) return kOK;
  fIsSetup = ( mode == kDefine );
  RVarDef vars[] = {
    { "x_fcp", "front track constraint x", "fFrontConstraintX" },
    { "y_fcp", "front track constraint y", "fFrontConstraintY" },
    { "x_bcp", "back track constraint x", "fBackConstraintX" },
    { "y_bcp", "back track constraint y", "fBackConstraintY" },
    { nullptr }
  };
  DefineVarsFromList( vars, mode );

  return 0;
}

//_____________________________________________________________________________
// Int_t SBSBigBite::End( THaRunBase* run )
// {
//   h1_yVx_bcp->Write();
//   h1_x_fcpVbcp->Write();
//   h1_yVx_fcp->Write();
//   h1_y_fcpVbcp->Write();
//   h1_dyVdx->Write();
// } 

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
  double EpsEtotRatio = 0;
  TIter next( fNonTrackingDetectors );
  while( auto* theNonTrackDetector =
	 static_cast<THaNonTrackingDetector*>( next() )) {
    //if(theNonTrackDetector->InheritsFrom("SBSBBShower")){
    if(theNonTrackDetector->InheritsFrom("SBSCalorimeter")){
      SBSBBTotalShower* BBTotalShower = reinterpret_cast<SBSBBTotalShower*>(theNonTrackDetector);
      //BBShower->EresMax();
      // gather here the info useful for
      if(BBTotalShower->GetShower()->GetNclust()){
	//cout << BBTotalShower->GetShower()->GetName() << " " << BBTotalShower->GetShower()->GetX() << " " << BBTotalShower->GetShower()->GetY() << " " << BBTotalShower->GetShower()->GetOrigin().Z() << " " << 1./(BBTotalShower->GetShower()->SizeRow()/sqrt(12)) << " " << 1./(BBTotalShower->GetShower()->SizeCol()/sqrt(12)) << endl;
	
	Etot+= BBTotalShower->GetShower()->GetECorrected();
	x_bcp+= -BBTotalShower->GetShower()->GetX()/(BBTotalShower->GetShower()->SizeRow()/sqrt(12));
	//y_bcp+= BBTotalShower->GetShower()->GetY()/(BBTotalShower->GetShower()->SizeCol()/sqrt(12));
	y_bcp = BBTotalShower->GetShower()->GetY();
	z_bcp+= BBTotalShower->GetShower()->GetOrigin().Z();
	npts++;
	sumweights_x+=1./(BBTotalShower->GetShower()->SizeRow()/sqrt(12));
	sumweights_y+=1.;//1./(BBTotalShower->GetShower()->SizeCol()/sqrt(12));
	//wx_bcp+=BBTotalShower->GetShower()->SizeRow()/sqrt(12);
	//wy_bcp+=BBTotalShower->GetShower()->SizeCol()/sqrt(12);
      }
      if(BBTotalShower->GetPreShower()->GetNclust()){
	//cout << BBTotalShower->GetPreShower()->GetName() << " " << BBTotalShower->GetPreShower()->GetX() << " " << BBTotalShower->GetPreShower()->GetY() << " " << BBTotalShower->GetPreShower()->GetOrigin().Z() << " " << 1./(BBTotalShower->GetPreShower()->SizeRow()/sqrt(12)) << " " << 1./(BBTotalShower->GetPreShower()->SizeCol()/sqrt(12)) << endl;
	
	Etot+= BBTotalShower->GetPreShower()->GetECorrected();
	EpsEtotRatio = BBTotalShower->GetPreShower()->GetECorrected()/Etot;
	fEpsEtotRatio.push_back(EpsEtotRatio);
	fEtot.push_back(Etot);
	//x_bcp+= -BBTotalShower->GetPreShower()->GetX()/(BBTotalShower->GetShower()->SizeRow()/sqrt(12));
	//y_bcp+= BBTotalShower->GetPreShower()->GetY()/(BBTotalShower->GetShower()->SizeCol()/sqrt(12));
	//z_bcp+= BBTotalShower->GetPreShower()->GetOrigin().Z();
	//npts++;
	//sumweights_x+=1./(BBTotalShower->GetPreShower()->SizeRow()/sqrt(12));
	//sumweights_y+=1./(BBTotalShower->GetPreShower()->SizeCol()/sqrt(12));
	//wx_bcp+=BBTotalShower->GetPreShower()->SizeRow()/sqrt(12);
	//wy_bcp+=BBTotalShower->GetPreShower()->SizeCol()/sqrt(12);
      }
    }
    
  }
  if(npts){
    x_bcp/=sumweights_x;
    y_bcp/=sumweights_y;
    //z_bcp/=npts;
    
    //wx_bcp/=npts;
    //wy_bcp/=npts;
    
    // std::cout << "Back constraint point x, y, z: " 
    //  	      << x_bcp << ", " << y_bcp << ", "<< z_bcp 
    //   //<< "; width x, y: " << wx_bcp << ", " << wy_bcp 
    //  	      << endl;
    
    // apply first order optics???
    // Yes, with the electron energy
    //TODO: replace hard-coded coefficients with optics coefficients
    // relationship 
    // x_5 = xfp+z_5*xpfp
    // thetabend = 10deg+xptar-xpfp
    // p ~= Ecalo
    // p*thetabend = 0.277+0.122*xfp-0.063*xpfp = Ecalo*(10deg+xptar-xpfp)
    // xptar = 0.523 * xfp - 0.414 * xpfp
    // 0.277+0.122*xfp-0.063*xpfp = Ecalo*(10deg+(0.523 * xfp -0.414 * xpfp)-xpfp)
    // 0.277+0.122*(x_5-z_5*xpfp)-0.063*xpfp = 
    //   Ecalo*(10deg + 0.523*(x_5-z_5*xpfp) + (-0.414-1)*xpfp) =>OK
    // 
    // 0.277 = fb_pinv_00000 = fb_pinv[0] = M_{p0}
    // 0.122 = fb_pinv_00001 = fb_pinv[1] = M_{px}
    // -0.063 = fb_pinv_00100 = fb_pinv[6] = M_{px'}
    // 0.523 = fb_xptar_00001 = fb_xptar[1] = M_{x'x}
    // -0.414 = fb_xptar_00100 = fb_xptar[6] = M_{x'x'}

    // fb_pinv_00000+fb_pinv_00001*(x_bcp-z_bcp*xpfp)+fb_pinv_00100*xpfp = 
    //   Ecalo*(10deg+fb_xptar_00001*(x_bcp-z_bcp*xpfp)+(fb_xptar_00100-1)*xpfp)

    // +fb_pinv[0]
    // +fb_pinv[1]*x_bcp
    // -fb_pinv[1]*z_bcp*xpfp
    // +fb_pinv[6]*xpfp
    //  =
    // +Etot*10.*TMath::DegToRad()
    // +Etot*fb_xptar[1]*x_bcp
    // -Etot*fb_xptar[1]*z_bcp*xpfp
    // +Etot*(fb_xptar[6]-1)*xpfp
    
    // -fb_pinv[1]*z_bcp*xpfp
    // +fb_pinv[6]*xpfp
    // -Etot*(fb_xptar[6]-1)*xpfp
    // +Etot*fb_xptar[1]*z_bcp*xpfp
    //  =
    // +Etot*10.*TMath::DegToRad()
    // +Etot*fb_xptar[1]*x_bcp
    // -fb_pinv[0]
    // -fb_pinv[1]*x_bcp
        
    double dx = (Etot*10.*TMath::DegToRad() -fb_pinv[0] + x_bcp * (Etot*fb_xptar[1]-fb_pinv[1])) /
      (-fb_pinv[1]*z_bcp+fb_pinv[6]+Etot*(fb_xptar[1]*z_bcp+1-fb_xptar[6]));
    double dy = y_bcp*0.251;//y_bcp*fb_yptar[3]/(1+fb_yptar[3]*z_bcp-fb_yptar[10]);
    
    //cout << "(x_bcp*(" << fb_xptar[1] << "*Etot-" << fb_pinv[1] << ")+" 
    //<< 10.*TMath::DegToRad() << "*Etot-" << fb_pinv[0] << ")/" << endl 
    //<< " (Etot*" << (fb_xptar[1]*z_bcp+1-fb_xptar[6]) << "+" << -fb_pinv[1]*z_bcp+fb_pinv[6] << ")" << endl;
    //cout << fb_yptar[3]/(1+fb_yptar[3]*z_bcp-fb_yptar[10]) << endl;
    
    //(x_bcp*(0.522891*Etot-0.121773)+0.174533*Etot-0.276919)/
    //(Etot*2.43719+-0.301327)
    
    //double dx_2 = (x_bcp*(0.522*Etot-0.121)+0.1729*Etot-0.278)/(Etot*2.224-0.249);
    //double dy_2 = y_bcp*0.251;
    
    z_fcp = 0;
    x_fcp = x_bcp+dx*(z_fcp-z_bcp);
    y_fcp = y_bcp+dy*(z_fcp-z_bcp);
    
    //cout << x_fcp-(x_bcp+dx_2*(z_fcp-z_bcp)) << " " << y_fcp-(y_bcp+dy_2*(z_fcp-z_bcp)) << endl;
    
    /*
    h1_yVx_bcp->Fill(x_bcp, y_bcp);
    h1_x_fcpVbcp->Fill(x_bcp, x_fcp);
    h1_yVx_fcp->Fill(x_fcp, y_fcp);
    h1_y_fcpVbcp->Fill(y_bcp, y_fcp);
    h1_dyVdx->Fill(dx, dy);
    */
    
    fFrontConstraintX.push_back(x_fcp);
    fFrontConstraintY.push_back(y_fcp);
    fBackConstraintX.push_back(x_bcp);
    fBackConstraintY.push_back(y_bcp);

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
  Int_t SBSBigBite::FindVertices( TClonesArray& tracks )
{
  
  // Reconstruct target coordinates for all tracks found.
  Int_t n_trk = tracks.GetLast()+1;
  for( Int_t t = 0; t < n_trk; t++ ) {
    auto* theTrack = static_cast<THaTrack*>( tracks.At(t) );
    CalcTargetCoords(theTrack);
  }  

  
  return 0;
}

void SBSBigBite::CalcTargetCoords( THaTrack* track )
{
  const double tracker_pitch_angle = 10.0*TMath::DegToRad();
  const double th_bb = 33.0*TMath::DegToRad();//temporary
 
  TVector3 BB_zaxis( sin(th_bb), 0, cos(th_bb) );
  TVector3 BB_xaxis(0,-1,0);
  TVector3 BB_yaxis = (BB_zaxis.Cross(BB_xaxis)).Unit();

  TVector3 spec_xaxis_fp,spec_yaxis_fp, spec_zaxis_fp;
  TVector3 spec_xaxis_tgt,spec_yaxis_tgt, spec_zaxis_tgt;
  
  spec_xaxis_tgt = BB_xaxis;
  spec_yaxis_tgt = BB_yaxis;
  spec_zaxis_tgt = BB_zaxis;
  
  spec_zaxis_fp = BB_zaxis;
  spec_yaxis_fp = BB_yaxis;
  spec_zaxis_fp.Rotate(-tracker_pitch_angle, spec_yaxis_fp);
  spec_xaxis_fp = spec_yaxis_fp.Cross(spec_zaxis_fp).Unit();
 
  Double_t x_fp, y_fp, xp_fp, yp_fp;
  //if( fCoordType == kTransport ) {
  x_fp = track->GetX();
  y_fp = track->GetY();
  xp_fp = track->GetTheta();
  yp_fp = track->GetPhi();
  //}

  double vx, vy, vz, px, py, pz;
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
  
  TVector3 phat_tgt_fit_global = phat_tgt_fit.X() * spec_xaxis_tgt +
    phat_tgt_fit.Y() * spec_yaxis_tgt +
    phat_tgt_fit.Z() * spec_zaxis_tgt;
  
  TVector3 phat_fp_fit(xp_fp, yp_fp, 1.0 );
  phat_fp_fit = phat_fp_fit.Unit();
  
  TVector3 phat_fp_fit_global = phat_fp_fit.X() * spec_xaxis_fp +
    phat_fp_fit.Y() * spec_yaxis_fp +
    phat_fp_fit.Z() * spec_zaxis_fp;
  
  thetabend_fit = acos( phat_fp_fit_global.Dot( phat_tgt_fit_global ) );
  
  p_fit = pthetabend_fit/thetabend_fit;
  vz_fit = -ytar_fit / (sin(th_bb) + cos(th_bb)*yptar_fit);
  
  pz = p_fit*sqrt( 1.0/(xptar_fit*xptar_fit+yptar_fit*yptar_fit+1) );
  px = xptar_fit * pz;
  py = yptar_fit * pz;
  
  track->SetTarget(xtar, ytar_fit, xptar_fit, yptar_fit);
  track->SetMomentum(p_fit);
  track->SetPvect(TVector3(px, py, pz));
  track->SetVertex(TVector3(0, 0, vz_fit));
  
}

//_____________________________________________________________________________
Int_t SBSBigBite::TrackCalc()
{
  // Additioal track calculations
  // Timing here???
  // PID info here???
  // TODO
  
  for( Int_t t = 0; t < fTracks->GetLast()+1; t++ ) {
    auto* theTrack = static_cast<THaTrack*>( fTracks->At(t) );
    THaPIDinfo* PIDinfo = new THaPIDinfo(2, 2);
    theTrack->SetPIDinfo(PIDinfo);
    CalcTimingPID(theTrack);
  }
  
  return 0;
}

//_____________________________________________________________________________
void SBSBigBite::CalcTimingPID(THaTrack* the_track)
{
  // Additioal track calculations
  // Timing here???
  // PID info here???
  // TODO
  
  //particles: 0: electron, 1: pion
  //detectors: 0: ps/sh
  THaPIDinfo* pidinfo = new THaPIDinfo(3, 2);
  pidinfo->SetDefaultPriors();
  
  the_track->SetPIDinfo(pidinfo);
  the_track->GetPIDinfo()->SetProb(0, 0, eproba_pssh(fEpsEtotRatio[the_track->GetIndex()]));
  the_track->GetPIDinfo()->SetProb(1, 0, eproba_gemcal(fEtot[the_track->GetIndex()]/the_track->GetP()));
  
  the_track->GetPIDinfo()->SetProb(0, 1, piproba_pssh(fEpsEtotRatio[the_track->GetIndex()]));
  the_track->GetPIDinfo()->SetProb(1, 1, piproba_gemcal(fEtot[the_track->GetIndex()]/the_track->GetP()));
  
  TIter next( fNonTrackingDetectors );
  while( auto* theNonTrackDetector =
	 static_cast<THaNonTrackingDetector*>( next() )) {
    //if(theNonTrackDetector->InheritsFrom("SBSBBShower")){
    // match a hodoscope cluster to a track:
    if(theNonTrackDetector->InheritsFrom("SBSTimingHodoscope")){
      SBSTimingHodoscope* TH = reinterpret_cast<SBSTimingHodoscope*>(theNonTrackDetector);
      
      //x, y of track at z = Z_hodoscope
      double x_track = the_track->GetX()+
	the_track->GetTheta()*TH->GetOrigin().Z();
      double y_track = the_track->GetY()+
	the_track->GetPhi()*TH->GetOrigin().Z();
      // Not sure what to use for the hodoscope. 
      // Perhaps we'd have to complete the class with a cluster
    }
    
    //bool match = false;
    int NGRINCHPMTs_match = 0;
    // match a GRINCH cluster to a track:
    if(theNonTrackDetector->InheritsFrom("SBSGRINCH")){
      SBSGRINCH* GRINCH = reinterpret_cast<SBSGRINCH*>(theNonTrackDetector);

      //x, y of track at z = Z_GRINCH
      double x_track = the_track->GetX()+
	the_track->GetTheta()*GRINCH->GetOrigin().Z();
      //double y_track = the_track->GetY()+
      //the_track->GetPhi()*GRINCH->GetOrigin().Z();

      for(int i = 0; i<GRINCH->GetNumClusters(); i++){
	SBSGRINCH_Cluster* gc_clus = GRINCH->GetCluster(i);
	
	if( fabs(x_track-gc_clus->GetXcenter()*fTrackGrinchClusCorr_1-fTrackGrinchClusCorr_0)<fTrackGrinchClusCorr_Sigma)NGRINCHPMTs_match = gc_clus->GetNHits();
	
	//if(y_track)
      }
    }

    the_track->GetPIDinfo()->SetProb(2, 0, eproba_grinch(NGRINCHPMTs_match));
    the_track->GetPIDinfo()->SetProb(2, 1, piproba_grinch(NGRINCHPMTs_match, the_track->GetP()));
    
  }
}


Double_t SBSBigBite::eproba_pssh(Double_t eps_etot_ratio)
{
  return 0;
}

Double_t SBSBigBite::eproba_gemcal(Double_t etot_p_ratio)
{
  return 0;
}

Double_t SBSBigBite::eproba_grinch(Int_t npmt)
{
  return 0;
}

Double_t SBSBigBite::piproba_pssh(Double_t eps_etot_ratio)
{
  return 0;
}

Double_t SBSBigBite::piproba_gemcal(Double_t etot_p_ratio)
{
  return 0;
}

Double_t SBSBigBite::piproba_grinch(Int_t npmt, Double_t p)
{
  return 0;
}
