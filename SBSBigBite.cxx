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
#include "SBSRasteredBeam.h"
#include "THaTrackingDetector.h"
#include "TH2D.h"

using namespace std;

ClassImp(SBSBigBite)

//_____________________________________________________________________________
SBSBigBite::SBSBigBite( const char* name, const char* description ) :
  THaSpectrometer( name, description )
{
  SetMultiTracks(false);
  SetTrSorting(false);
  fTrackerPitchAngle = 10.0;
  fDetectorStackYaw = 0;
  fDetectorStackRoll = 0;
  fOpticsOrder = -1;
  fFrontConstraintWidthX = 1.5;
  fFrontConstraintWidthY = 1.5;
  fBackConstraintWidthX = 1.5;
  fBackConstraintWidthY = 1.5;
  fTrackGrinchClusCorr_0 = 0.0;
  fTrackGrinchClusCorr_1 = 0.0;
  fTrackGrinchClusCorr_Sigma = 1.5;
  fPtheta_00000 = 0.0;
  fPtheta_10000 = 0.0;
  fPtheta_00100 = 0.0;
  fXptar_10000 = 0.0;
  fXptar_00100 = 0.0;
  fYtar_01000 = 0.0;
  fYtar_00010 = 0.0;
  fECaloFudgeFactor = 1.0;
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
  f_xtg_exp.clear();
  fEpsEtotRatio.clear();
  fEtot.clear();
  fEtotPratio.clear();
  fEpsEtotRatio.clear();
  fFrontConstraintX.clear();
  fFrontConstraintY.clear();
  fBackConstraintX.clear();
  fBackConstraintY.clear();
  fProbaE.clear();
  fProbaPi.clear();
}


//_____________________________________________________________________________
Int_t SBSBigBite::ReadDatabase( const TDatime& date )
{
  // Hack from THaVDC::ReadDatabase()
  const char* const here = "SBSBigBite::ReadDatabase()";
  
  //THaSpectrometer::ReadRunDatabase();
  
  FILE* file = OpenFile( date );
  if( !file ){
    std::cerr << here << "(): database not found!"<< std::endl;
    return kFileError;
  }

  int pidflag = fPID ? 1 : 0;
  
  std::vector<Double_t> firstgem_offset;
  
  std::vector<Double_t> optics_param;
  std::vector<Double_t> pssh_pidproba;
  std::vector<Double_t> pcal_pidproba;
  std::vector<Double_t> grinch_pidproba;
  const DBRequest request[] = {
    { "tracker_pitch_angle",    &fTrackerPitchAngle, kDouble,  0, 1, 1},
    { "tracker_yaw_angle",    &fDetectorStackYaw, kDouble,  0, 1, 1},
    { "tracker_roll_angle",    &fDetectorStackRoll, kDouble,  0, 1, 1},
    { "firstgemoffset_xyz",    &firstgem_offset, kDoubleV,  0, 1, 1},
    { "optics_order",    &fOpticsOrder, kUInt,  0, 1, 1},
    { "optics_parameters", &optics_param, kDoubleV, 0, 1, 1},
    { "ecalo_fudgefactor", &fECaloFudgeFactor, kDouble, 0, 1, 1},
    { "do_pid",    &pidflag, kInt,  0, 1, 1},
    { "frontconstraintwidth_x", &fFrontConstraintWidthX, kDouble, 0, 1, 0},
    { "frontconstraintwidth_y", &fFrontConstraintWidthY, kDouble, 0, 1, 0},
    { "backconstraintwidth_x", &fBackConstraintWidthX, kDouble, 0, 1, 0},
    { "backconstraintwidth_y", &fBackConstraintWidthY, kDouble, 0, 1, 0},
    { "trackgrinchcorr_const", &fTrackGrinchClusCorr_0, kDouble, 0, 1, 0},
    { "trackgrinchcorr_slope", &fTrackGrinchClusCorr_1, kDouble, 0, 1, 0},
    { "trackgrinchcorr_sigma", &fTrackGrinchClusCorr_Sigma, kDouble, 0, 1, 0},
    { "psshPIDprobatable",    &pssh_pidproba, kDoubleV,  0, 1, 0},
    { "pcalPIDprobatable",    &pcal_pidproba, kDoubleV,  0, 1, 0},
    { "grinchPIDpbins",    &fP_table, kDoubleV,  0, 1, 0},
    { "grinchPIDprobatable",    &grinch_pidproba, kDoubleV,  0, 1, 0},
    {0}
  };
    
  Int_t status = LoadDB( file, date, request, fPrefix, 1 ); //The "1" after fPrefix means search up the tree
  fclose(file);
  if( status != 0 ){
    return status;
  }
  
  if(fECaloFudgeFactor!=1.0)cout << "Setting the calorimeter energy fudge factor to " << fECaloFudgeFactor << endl;
  
  fPID = ( pidflag != 0 );

  fDetectorStackPitch = (fTrackerPitchAngle-10.0)*TMath::DegToRad();
  fTrackerPitchAngle*= TMath::DegToRad();
  fDetectorStackYaw*= TMath::DegToRad();
  fDetectorStackRoll*= TMath::DegToRad();
  
  double fDetectorStackThSph;
  double fDetectorStackPhSph;

  //can we maybe neglect roll for the moment? I think we can.
  
  GeoToSph(fDetectorStackYaw, fDetectorStackPitch, 
	   fDetectorStackThSph, fDetectorStackPhSph);
  
  Double_t st = TMath::Sin(fDetectorStackThSph);
  Double_t ct = TMath::Cos(fDetectorStackThSph);
  Double_t sp = TMath::Sin(fDetectorStackPhSph);
  Double_t cp = TMath::Cos(fDetectorStackPhSph);

  // Compute the rotation from Optics (ideal) to Detector (actual) and vice versa.
  // If the pitch is positive i.e. the actual tracker pitch angle is larger than 
  // the nominal 10deg, the transformation from the detector coordinates 
  // to the optics coordinate should transfomr the track dx into a more negative dx
  // The transformation effectively does that if 
  // fDetectorStackPitch>0 when fTrackerPitchAngle>10.
  // But, for consistency reasons, I would like 
  // fDetectorStackPitch<0 when fTrackerPitchAngle>10.
  // Since I'm too dumb (and lack time) to rewrite the transformation, 
  // I prefer to evaluate the transformation with fDetectorStackPitch>0 
  // and *then* transform fDetectorStackPitch into -fDetectorStackPitch...
  Double_t norm = TMath::Sqrt(ct*ct + st*st*cp*cp);
  //TVector3 nx( st*st*sp*cp/norm, -norm, st*ct*sp/norm );
  //TVector3 ny( ct/norm,          0.0,   -st*cp/norm   );
  //TVector3 nz( st*cp,            st*sp, ct            );
  TVector3 nx( norm,   st*st*sp*cp/norm, st*ct*sp/norm );
  TVector3 ny( 0.0,    ct/norm,          -st*cp/norm   );
  TVector3 nz( -st*sp, st*cp,            ct            );
  fDet2OptRot.SetToIdentity().RotateAxes( nx, ny, nz );
  fOpt2DetRot = fDet2OptRot.Inverse();
  
  fDetectorStackPitch = -fDetectorStackPitch;
  
  cout << firstgem_offset.size() << endl;
  if(firstgem_offset.size()==3){
    fFirstGEMLayerOffset = TVector3(firstgem_offset[0],
				    firstgem_offset[1],
				    firstgem_offset[2]);
  }else{
    fFirstGEMLayerOffset = TVector3(0, 0, 0);
    std::cout << "First GEM offset not set or set improperly, setting to 0, 0, 0" << std::endl;
  }
  
  //do we have non tracking detectors
  bool nontrackdet = false;
  TIter next( fNonTrackingDetectors );
  while( auto* theNonTrackDetector =
	 static_cast<THaNonTrackingDetector*>( next() )) {
    cout << theNonTrackDetector->GetName() << endl;
    nontrackdet = true;
  }//if we do not find any non tracking detectors, we force fPID to be false.
  if(!nontrackdet)fPID = false;

  if(fOpticsOrder>=0){
    int nterms = 0;
    
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

    int nparams = 9*nterms;

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
    
    for(int i=0; i<nterms; i++){
      fb_xptar[i] = optics_param[9*i];
      fb_yptar[i] = optics_param[9*i+1];
      fb_ytar[i] = optics_param[9*i+2];
      fb_pinv[i] = optics_param[9*i+3];
      f_om[i] = int(optics_param[9*i+4]);
      f_ol[i] = int(optics_param[9*i+5]);
      f_ok[i] = int(optics_param[9*i+6]);
      f_oj[i] = int(optics_param[9*i+7]);
      f_oi[i] = int(optics_param[9*i+8]);
      if(f_om[i]==0 && f_ol[i]==0 && f_ok[i]==0 && f_oj[i]==0 && f_oi[i]==0){
	fPtheta_00000 = fb_pinv[i];
      }
      if(f_om[i]==1 && f_ol[i]==0 && f_ok[i]==0 && f_oj[i]==0 && f_oi[i]==0){
	fPtheta_10000 = fb_pinv[i];
	fXptar_10000 = fb_xptar[i];
      }
      if(f_om[i]==0 && f_ol[i]==0 && f_ok[i]==1 && f_oj[i]==0 && f_oi[i]==0){
	fPtheta_00100 = fb_pinv[i];
	fXptar_00100 = fb_xptar[i];
      }
      if(f_om[i]==0 && f_ol[i]==1 && f_ok[i]==0 && f_oj[i]==0 && f_oi[i]==0)
	fYtar_01000 = fb_ytar[i];
      if(f_om[i]==0 && f_ol[i]==0 && f_ok[i]==0 && f_oj[i]==1 && f_oi[i]==0)
	fYtar_00010 = fb_ytar[i];
    }
  }
  
  //PID stuff
  fEpsEtotRatio_table.clear();
  fProba_e_PSSH_table.clear();
  fProba_pi_PSSH_table.clear();
  
  if(!pssh_pidproba.empty()){
    int npts = pssh_pidproba.size()/3;
    fEpsEtotRatio_table.resize(npts);
    fProba_e_PSSH_table.resize(npts);
    fProba_pi_PSSH_table.resize(npts);
    
    for(int i = 0; i<npts; i++){
      fEpsEtotRatio_table[i] = pssh_pidproba[3*i];
      fProba_e_PSSH_table[i] = pssh_pidproba[3*i+1];
      fProba_pi_PSSH_table[i] = pssh_pidproba[3*i+2];
    }
  }
  
  //PID stuff
  fEtotPratio_table.clear();
  fProba_e_PCAL_table.clear();
  fProba_pi_PCAL_table.clear();
  
  if(!pcal_pidproba.empty()){
    int npts = pcal_pidproba.size()/3;
    fEtotPratio_table.resize(npts);
    fProba_e_PCAL_table.resize(npts);
    fProba_pi_PCAL_table.resize(npts);
    
    for(int i = 0; i<npts; i++){
      fEtotPratio_table[i] = pcal_pidproba[3*i];
      fProba_e_PCAL_table[i] = pcal_pidproba[3*i+1];
      fProba_pi_PCAL_table[i] = pcal_pidproba[3*i+2];
    }
  }
  
  fNGRINCHPMTs_table.clear();
  fProba_e_GRINCH_table.clear();
  fProba_pi_GRINCH_table.clear();
  
  if(!grinch_pidproba.empty() && fP_table.size()){
    int n_ppts = 2+fP_table.size();
    int npts = grinch_pidproba.size()/(n_ppts);
    fNGRINCHPMTs_table.resize(npts);
    fProba_e_GRINCH_table.resize(npts);
    fProba_pi_GRINCH_table.resize(fP_table.size());
    for(int j = 0; j<fP_table.size(); j++){
      fProba_pi_GRINCH_table[j].resize(npts);
    }
    for(int i = 0; i<npts; i++){
      fNGRINCHPMTs_table[i] = grinch_pidproba[n_ppts*i];
      fProba_e_GRINCH_table[i] = grinch_pidproba[n_ppts*i+1];
      for(int j = 0; j<fP_table.size(); j++){
	fProba_pi_GRINCH_table[j][i] = grinch_pidproba[n_ppts*i+2+j];
      }
    }
  }
  
  fIsInit = true;
  return kOK;
}

Int_t SBSBigBite::DefineVariables( EMode mode ){
  THaSpectrometer::DefineVariables(mode);
  
  if( mode == kDefine and fIsSetup ) return kOK;
  fIsSetup = ( mode == kDefine );
  
  
  // removing all that stuff since apparently I'm not able to code properly...

  RVarDef beamtrackvars[] = {
    { "tr.tg_x", "track constraint y", "f_xtg_exp" },
    { nullptr }
  };
  DefineVarsFromList( beamtrackvars, mode );
  
  RVarDef constraintvars[] = {
    { "x_fcp", "front track constraint x", "fFrontConstraintX" },
    { "y_fcp", "front track constraint y", "fFrontConstraintY" },
    { "x_bcp", "back track constraint x", "fBackConstraintX" },
    { "y_bcp", "back track constraint y", "fBackConstraintY" },
    { nullptr }
  };
  DefineVarsFromList( constraintvars, mode );
  
  RVarDef pidvars[] = {
    { "eps_over_etot", "electron probability", "fEpsEtotRatio" },
    { "etot_over_p", "electron probability", "fEtotPratio" },
    { "prob_e", "electron probability", "fProbaE" },
    { "prob_pi", "pion probability", "fProbaPi" },
    { nullptr }
  };
  DefineVarsFromList( pidvars, mode );
  
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
  //std::cout << " SBSBigBite::CoarseTrack()...";
  THaSpectrometer::CoarseTrack();
  // TODO
  //std::cout << " call SBSBigBite::CoarseTrack" << std::endl;
  //std::cout << "done" << std::endl;
  return 0;
}

//_____________________________________________________________________________
Int_t SBSBigBite::CoarseReconstruct()
{

  //std::cout << "SBSBigBite::CoarseReconstruct()..."; 
  // Coarse Reconstruction of particles in spectrometer
  THaSpectrometer::CoarseReconstruct(); 
  // TODO
  // fetch the clusters from SBSBBShower detectors
  // FOR NOW: fetch the highest clusters from SBSBBShower detectors
  double x_fcp = 0, y_fcp = 0, z_fcp = 0;
  double x_bcp = 0, y_bcp = 0, z_bcp = 0;
  double sumweights_x = 0, sumweights_y = 0;
  double Etot = 0;
  //npts is incremented only if there are clusters in the preshower and shower
  int npts = 0;
  double EpsEtotRatio = 0;
  TIter next( fNonTrackingDetectors );
  while( auto* theNonTrackDetector =
	 static_cast<THaNonTrackingDetector*>( next() )) {
    //if(theNonTrackDetector->InheritsFrom("SBSBBShower")){
    //if(theNonTrackDetector->InheritsFrom("SBSCalorimeter")){
    if(theNonTrackDetector->InheritsFrom("SBSBBTotalShower")){//More explicit
      //cout << "found SBSBBTtotalShower" << endl;
      SBSBBTotalShower* BBTotalShower = reinterpret_cast<SBSBBTotalShower*>(theNonTrackDetector);
      //BBShower->EresMax();
      // explicitely return and 0 if there is no cluster in the calorimeter;
      if(BBTotalShower->GetShower()->GetNclust()==0 || BBTotalShower->GetPreShower()->GetNclust()==0){
	//cout << BBTotalShower->GetPreShower()->GetNclust() << " clusters in preshower, " << BBTotalShower->GetShower()->GetNclust() << " clusters in shower, kill here" << endl;
	return 0;
      }
      //cout << "shower cluster E " << BBTotalShower->GetShower()->GetECorrected() 
      //   << " X = " << BBTotalShower->GetShower()->GetX() 
      //   << " Y = " << BBTotalShower->GetShower()->GetY()
      //   << endl;
      //cout << "preshower cluster E " << BBTotalShower->GetPreShower()->GetECorrected() 
      //   << " X = " << BBTotalShower->GetPreShower()->GetX() 
      //   << " Y = " << BBTotalShower->GetPreShower()->GetY()
      //   << endl;
      if(GetMultiTracks()){
	std::vector<SBSCalorimeterCluster*> ShowerClusters = BBTotalShower->GetShower()->GetClusters();
	std::vector<SBSCalorimeterCluster*> PreShowerClusters = BBTotalShower->GetPreShower()->GetClusters();
	z_bcp = BBTotalShower->GetShower()->GetOrigin().Z();
	sumweights_x = 1./(BBTotalShower->GetShower()->SizeRow()/sqrt(12));
	sumweights_y = 1.;
	//cout << BBTotalShower->GetShower()->GetECorrected() << endl;
	for(int i = 0; i<ShowerClusters.size(); i++){
	  Etot = ShowerClusters[i]->GetE();
	  npts = 1;
	  x_bcp = ShowerClusters[i]->GetX()/(BBTotalShower->GetShower()->SizeRow()/sqrt(12));
	  y_bcp = ShowerClusters[i]->GetY();
	  
	  if(BBTotalShower->PSMatchClusIdx(i)<PreShowerClusters.size()){
	    Etot+= PreShowerClusters[BBTotalShower->PSMatchClusIdx(i)]->GetE();
	    fEpsEtotRatio.push_back(EpsEtotRatio);
	    fEtot.push_back(Etot);

	    x_bcp/=sumweights_x;
	    y_bcp/=sumweights_y;
	    
	    double dx = (Etot*fECaloFudgeFactor*fTrackerPitchAngle - fPtheta_00000 + x_bcp * (Etot*fECaloFudgeFactor*fXptar_10000-fPtheta_10000)) /
	      (-fPtheta_10000*z_bcp+fPtheta_00100+Etot*(fXptar_10000*z_bcp+1-fXptar_00100));
	    double dy = y_bcp*fb_ytar[3]/(fb_ytar[3]*z_bcp-fb_ytar[10]);
	    
	    z_fcp = 0;
	    x_fcp = x_bcp+dx*(z_fcp-z_bcp);
	    y_fcp = y_bcp+dy*(z_fcp-z_bcp);
	    
	    fFrontConstraintX.push_back(x_fcp);
	    fFrontConstraintY.push_back(y_fcp);
	    fBackConstraintX.push_back(x_bcp);
	    fBackConstraintY.push_back(y_bcp);
	    
	    //now what???
	  }
	}//end loop on Shower clusters
      }else{//end if(GetMultiTracks())
	//if(BBTotalShower->GetShower()->GetNclust()){
	//cout << BBTotalShower->GetShower()->GetName() << " " << BBTotalShower->GetShower()->GetX() << " " << BBTotalShower->GetShower()->GetY() << " " << BBTotalShower->GetShower()->GetOrigin().Z() << " " << 1./(BBTotalShower->GetShower()->SizeRow()/sqrt(12)) << " " << 1./(BBTotalShower->GetShower()->SizeCol()/sqrt(12)) << endl;
	// TODO: so far we use only the "main" cluster... 
	//       we might want to check the others...
	//y_bcp+= BBTotalShower->GetShower()->GetY()/(BBTotalShower->GetShower()->SizeCol()/sqrt(12));
	Etot+= BBTotalShower->GetShower()->GetECorrected();
	x_bcp+= BBTotalShower->GetShower()->GetX()/pow(BBTotalShower->GetShower()->SizeRow()/sqrt(12.),2);
	y_bcp = BBTotalShower->GetShower()->GetY();
	z_bcp+= BBTotalShower->GetShower()->GetOrigin().Z();
	npts++;
	sumweights_x+=1./pow(BBTotalShower->GetShower()->SizeRow()/sqrt(12),2);
	sumweights_y+=1.;//1./(BBTotalShower->GetShower()->SizeCol()/sqrt(12));
	//}
	//if(BBTotalShower->GetPreShower()->GetNclust()){
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
	//}
      
	//if we're here we've found the totalshower
	//if(npts){
	x_bcp/=sumweights_x;
	y_bcp/=sumweights_y;
	//z_bcp/=npts;
	
	//std::cout << "Back constraint point x, y, z: " 
	//	  << x_bcp << ", " << y_bcp << ", "<< z_bcp 
	//	  << endl;
	
        // to account for the angle and position offsets of the detector stack: 
	// simplest approximate way to do it: 
	// we have x_bcp^det 
	// x_bcp^opt ~ x_bcp^det+fFirstGEMLayerOffset.X()+fDetectorStackPitch*z_bcp
	// => calculate dx^opt(x_bcp^opt)
	// => calculate x_fcp^opt(dx^opt)
	// => x_fcp^det = x_fcp^opt+fFirstGEMLayerOffset.X()
	// similarly with y_bcp^det 
	// y_bcp^opt ~ y_bcp^det+fFirstGEMLayerOffset.Y()+fDetectorStackYaw*z_bcp
	// => calculate dy^opt(y_bcp^opt)
	// => calculate x_fcp^opt(dy^opt)
	// => y_fcp^det = y_fcp^opt-fFirstGEMLayerOffset.Y()
	// Of course all of the above would hold only for angles less than a few degrees.
	
	//transformation in optics coordinate
	//commenting these out for now: AJRP
	//x_bcp+=fFirstGEMLayerOffset.X()+fDetectorStackPitch*z_bcp;
	//y_bcp+=fFirstGEMLayerOffset.Y()+fDetectorStackYaw*z_bcp;
	
	// Use 10.0 degrees instead of fTrackerPitchAngle 
	// because we are now in the "ideal" system.
	//double xp_bcp = (Etot*fECaloFudgeFactor*fTrackerPitchAngle - fPtheta_00000 + x_bcp * (Etot*fXptar_10000-fPtheta_10000)) /
	//  (-fPtheta_10000*z_bcp+fPtheta_00100+Etot*fECaloFudgeFactor*(fXptar_10000*z_bcp+1-fXptar_00100));
	double yp_bcp = y_bcp*fYtar_01000/(fYtar_01000*z_bcp-fYtar_00010);

	double Efudge = Etot * fECaloFudgeFactor;
	
	double xp_bcp = ( fPtheta_00000 + fPtheta_10000 * x_bcp - Efudge * ( fTrackerPitchAngle + fXptar_10000 * x_bcp ) ) /
	  ( Efudge * (fXptar_00100 - 1.0 - fXptar_10000*z_bcp) + fPtheta_10000 * z_bcp - fPtheta_00100 );
	
	// The dy equation is correct under the assumption ytarget = 0: can we refine?
	// double dx = (Etot*10.*TMath::DegToRad() -fb_pinv[0] + x_bcp * (Etot*fb_xptar[1]-fb_pinv[1])) /
	// (-fb_pinv[1]*z_bcp+fb_pinv[6]+Etot*(fb_xptar[1]*z_bcp+1-fb_xptar[6]));
	// double dy = y_bcp*fb_ytar[3]/(fb_ytar[3]*z_bcp-fb_ytar[10]);

	//The x' equation is:
	// Ecalo = p
	// ECALO*thetabend = fPtheta_00000 + fPtheta_10000 * xfp
	// ECALO*( 10 deg. + xptar - xpfp ) = ptheta0 + pthetax * xfp
	// xptar = xptar_0 + xptar_x * xfp + xptar_xp * xpfp
	// xbcp = xfp + xpfp * zbcp
	// xfp = xbcp - xpfp * zbcp
	// ECALO * ( 10 deg. + xptar - xpfp ) = pth0 + Mpthx * (xbcp - xpfp * zbcp );
	// ECALO * ( 10 deg. + xp0 + Mxpx * (xbcp - xpfp * zbcp) + Mxpxp * xpfp - xpfp ) = pth0 + Mpthx * (xbcp-xpfp*zbcp)

	// ECAL * (10 deg. + xp0 + Mxpx * xbcp) + xpfp * ECALO * ( Mxpxp - Mxpx * zbcp - 1 ) = Mpthx * xbcp - xpfp*Mpthx*zbcp
	// --> xpfp*[ ECALO * ( Mxpxp - Mxpx*zbcp - 1 ) - Mpthx*zbcp ] = 
	
	//cout << "(x_bcp*(" << fXptar_10000 << "*Etot-" << fPtheta_10000 << ")+" 
	//   << 10.*TMath::DegToRad() << "*Etot-" << fPtheta_00000 << ")/" << endl 
	//   << " (Etot*" << (fXptar_10000*z_bcp+1-fXptar_00100) << "+" << -fPtheta_10000*z_bcp+fPtheta_00100 << ")" << endl;
	//cout << fYtar_01000/(fYtar_01000*z_bcp-fYtar_00010) << endl;
	
	z_fcp = 0;
	x_fcp = x_bcp+xp_bcp*(z_fcp-z_bcp);
	y_fcp = y_bcp+yp_bcp*(z_fcp-z_bcp);

	//commenting this out for now
	//x_bcp+=-fFirstGEMLayerOffset.X();
	//y_bcp+=-fFirstGEMLayerOffset.Y();
	
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
	
	// std::cout << "Front constraint point x, y, z: " 
	// 	      << x_fcp << ", " << y_fcp << ", "<< z_fcp 
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
	    BBGEM->SetBackConstraintWidth(fBackConstraintWidthX, 
					  fBackConstraintWidthY);
	    /*
	      BBGEM->SetFrontConstraintPoint(TVector3(x_fcp, y_fcp, z_fcp));
	      BBGEM->SetBackConstraintPoint(TVector3(x_bcp, y_bcp, z_bcp));
	      BBGEM->SetFrontConstraintWidth(TVector2(fFrontConstraintWidthX, fFrontConstraintWidthY));
	      BBGEM->SetBackConstraintWidth(TVector2(fBackConstraintWidthX, fBackConstraintWidthY));
	    */
	  }
	}
	//}//end if(npts>0);
      }//end else of if(multitracks)
      
    }//end if(inheritsfrom(SBSCalorimeter))
    
  }//end loop on non tracking detectors
  
  //std::cout << " call SBSBigBite::CoarseReconstruct" << std::endl;
  //THaSpectrometer::CoarseReconstruct();

  //std::cout << "done." << std::endl; 
  
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
  //std::cout << "SBSBigBite::Reconstruct()..." << std::endl;
  
  THaSpectrometer::Reconstruct();
  // TODO

  //std::cout << "Done..." << std::endl;
  
  return 0;
  
}

//_____________________________________________________________________________
Int_t SBSBigBite::FindVertices( TClonesArray& tracks )
{

  //std::cout << "SBSBigBite::FindVertices()...";
  // Reconstruct target coordinates for all tracks found.
  Int_t n_trk = tracks.GetLast()+1;
  for( Int_t t = 0; t < n_trk; t++ ) {
    auto* theTrack = static_cast<THaTrack*>( tracks.At(t) );
    CalcOpticsCoords(theTrack);

    if(fOpticsOrder>=0)CalcTargetCoords(theTrack);
  }
  
  //sort on other criteria than chi^2
  if( GetTrSorting() ) {
    fTracks->Sort();
    // Reassign track indexes. Sorting may have changed the order
    for( int i = 0; i < fTracks->GetLast()+1; i++ ) {
      auto* theTrack = static_cast<THaTrack*>( fTracks->At(i) );
      assert( theTrack );
      theTrack->SetIndex(i);
    }
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
  } else
    fGoldenTrack = nullptr;  

  //std::cout << "done." << std::endl;
  
  return 0;
}


void SBSBigBite::CalcOpticsCoords( THaTrack* track )
{
  Double_t x_fp, y_fp, xp_fp, yp_fp;
  Double_t px, py, pz;// NB: not the actual momentum!
  
  x_fp = track->GetX()+fFirstGEMLayerOffset.X();
  y_fp = track->GetY()+fFirstGEMLayerOffset.Y();
  xp_fp = track->GetTheta();
  yp_fp = track->GetPhi();
  
  // cout << x_fp << " " << y_fp << " " << xp_fp << " " << yp_fp << endl;
  
  // x_fp = track->GetX()+fFirstGEMLayerOffset.X();
  // y_fp = track->GetY()+fFirstGEMLayerOffset.Y();
  
  // cout << x_fp << " " << y_fp << " " 
  //  << xp_fp+fDetectorStackPitch << " " << yp_fp+fDetectorStackYaw << endl;
  // provided involved angles are small maybe we could almost do 
  // xp_fp+= -fDetectorStackPitch
  // yp_fp+= fDetectorStackYaw
  // and forgo transformation calculation altogether? 
  // defo should use that as a sanity check
  pz = sqrt( 1.0/(xp_fp*xp_fp+yp_fp*yp_fp+1.0) );
  px = xp_fp * pz;
  py = yp_fp * pz;
  
  TVector3 p_trk(px, py, pz);
  p_trk.Transform(fDet2OptRot);
  
  xp_fp = p_trk.X()/p_trk.Z();
  yp_fp = p_trk.Y()/p_trk.Z();
  
  // cout << x_fp << " " << y_fp << " " << xp_fp << " " << yp_fp << endl;
  
  track->Set(x_fp, y_fp, xp_fp, yp_fp);
  
  
}

void SBSBigBite::CalcTargetCoords( THaTrack* track )
{
  //std::cout << "SBSBigBite::CalcTargetCoords()...";
  
  //const double //make it configurable
  const double th_bb = GetThetaGeo();//retrieve the actual angle
  
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
  spec_zaxis_fp.Rotate(-fTrackerPitchAngle, spec_yaxis_fp);
  spec_xaxis_fp = spec_yaxis_fp.Cross(spec_zaxis_fp).Unit();
 
  Double_t x_fp, y_fp, xp_fp, yp_fp;
  //if( fCoordType == kTransport ) {
  x_fp = track->GetX();
  y_fp = track->GetY();
  xp_fp = track->GetTheta();
  yp_fp = track->GetPhi();
  //}
  //cout << x_fp << " " << y_fp << " " << xp_fp << " " << yp_fp << endl;

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

  TVector3 pvect_BB = TVector3(px, py, pz);
  
  px = +pvect_BB.Z()*sin(th_bb)+pvect_BB.Y()*cos(th_bb);
  py = -pvect_BB.X();
  pz = pvect_BB.Z()*cos(th_bb)-pvect_BB.Y()*sin(th_bb);
  
  //retrieve beam position, if available, to calculate xtar.
  TIter aiter(gHaApps);
  THaApparatus* app = 0;
  while( (app=(THaApparatus*)aiter()) ){
    if(app->InheritsFrom("SBSRasteredBeam")){
      SBSRasteredBeam* RasterBeam = reinterpret_cast<SBSRasteredBeam*>(app);
      //double xbeam = RasterBeam->GetPosition().X();
      double ybeam = RasterBeam->GetPosition().Y();
      
      xtar = - ybeam - cos(GetThetaGeo()) * vz_fit * xptar_fit;
    }
    //cout << var->GetName() << endl;
  }
  f_xtg_exp.push_back(xtar);
  
  track->SetTarget(xtar, ytar_fit, xptar_fit, yptar_fit);
  track->SetMomentum(p_fit);
  track->SetPvect(TVector3(px, py, pz));
  track->SetVertex(TVector3(0, 0, vz_fit));

  //cout << px << " " << py << " " << pz << "   " << vz_fit << endl;
  //cout << track->GetLabPx() << " " << track->GetLabPy() << " " << track->GetLabPz() 
  //   << "   " << track->GetVertexZ() << endl;
  
  //std::cout << "Done." << std::endl;
}

//_____________________________________________________________________________
Int_t SBSBigBite::TrackCalc()
{
  // Additioal track calculations
  // Timing calculation goes here
  // PID calculation goes here!
  for( Int_t t = 0; t < fTracks->GetLast()+1; t++ ) {
    auto* theTrack = static_cast<THaTrack*>( fTracks->At(t) );
    CalcTrackTiming(theTrack);
  }
  
  return 0; 
}

//_____________________________________________________________________________
Int_t SBSBigBite::CalcPID()
{
  // PID calculation goes here!
  for( Int_t t = 0; t < fTracks->GetLast()+1; t++ ) {
    auto* theTrack = static_cast<THaTrack*>( fTracks->At(t) );
    CalcTrackPID(theTrack);
  }
  
}

//_____________________________________________________________________________
void SBSBigBite::CalcTrackTiming(THaTrack* the_track)
{
  TIter next( fNonTrackingDetectors );
  while( auto* theNonTrackDetector =
	 static_cast<THaNonTrackingDetector*>( next() )) {
    // match a hodoscope cluster to a track:
    //the hodoscope has to be found for anything to be done.
    if(theNonTrackDetector->InheritsFrom("SBSTimingHodoscope")){
      SBSTimingHodoscope* TH = reinterpret_cast<SBSTimingHodoscope*>(theNonTrackDetector);
      
      //x, y of track at z = Z_hodoscope
      double x_track = the_track->GetX(TH->GetOrigin().Z());
      double y_track = the_track->GetY(TH->GetOrigin().Z());
      //cout << x_track << " " << y_track << endl;
      // Not sure what to use for the hodoscope. 
      // Perhaps we'd have to complete the class with a cluster
      
      for(int i=0; i<TH->GetNClusters(); i++){
	SBSTimingHodoscopeCluster* clus = TH->GetCluster(i);
	if(clus->GetXmean()-clus->GetSize()*TH->SizeCol()/2<x_track && 
	   x_track<clus->GetXmean()+clus->GetSize()*TH->SizeCol()/2){
	  //std::cout << clus->GetSize() << " " << clus->GetTmean() << " " << clus->GetXmean() << " " << clus->GetYmean() << " " << clus->GetToTmean() << std::endl;
	  the_track->SetTime(clus->GetTmean());
	}  
      }
    }//end if inherits from hodoscope
    
    if(theNonTrackDetector->InheritsFrom("SBSBBTotalShower")){
      SBSBBTotalShower* BBTotalShower = reinterpret_cast<SBSBBTotalShower*>(theNonTrackDetector);
      double Z_cst =  BBTotalShower->GetShower()->GetOrigin().Z();
      //check that the track we consider is consistent with the calorimeter constraint 
      int i_match = -1;
      for(int i = 0; i<fEtot.size(); i++){
	/*
	cout << "back X: " << fBackConstraintX[i]-fBackConstraintWidthX 
	     << " <? " << the_track->GetX(Z_cst) << " <? " 
	     << fBackConstraintX[i]+fBackConstraintWidthX << endl;
	cout << "back Y: " << fBackConstraintY[i]-fBackConstraintWidthY 
	     << " <? " << the_track->GetY(Z_cst) << " <? " 
	     << fBackConstraintY[i]+fBackConstraintWidthY << endl;
	cout << "front X: " << fFrontConstraintX[i]-fFrontConstraintWidthX 
	     << " <? " << the_track->GetX() << " <? " 
	     << fFrontConstraintX[i]+fFrontConstraintWidthX << endl;
	cout << "front Y: " << fFrontConstraintY[i]-fFrontConstraintWidthY 
	     << " <? " << the_track->GetY() << " <? " 
	     << fFrontConstraintY[i]+fFrontConstraintWidthY << endl;
	*/
	/*
	if(fBackConstraintX[i]-fBackConstraintWidthX > the_track->GetX(Z_cst) ||
	   the_track->GetX(Z_cst) > fBackConstraintX[i]+fBackConstraintWidthX)
	  cout << "back X: " << fBackConstraintX[i]-fBackConstraintWidthX 
	       << " <? " << the_track->GetX(Z_cst) << " <? " 
	       << fBackConstraintX[i]+fBackConstraintWidthX << endl;
	
	if(fBackConstraintY[i]-fBackConstraintWidthY > the_track->GetY(Z_cst) ||
	   the_track->GetY(Z_cst) > fBackConstraintY[i]+fBackConstraintWidthY)
	  cout << "back Y: " << fBackConstraintY[i]-fBackConstraintWidthY 
	       << " <? " << the_track->GetY(Z_cst) << " <? " 
	       << fBackConstraintY[i]+fBackConstraintWidthY << endl;
	
	if(fFrontConstraintX[i]-fFrontConstraintWidthX > the_track->GetX() ||
	   the_track->GetX() > fFrontConstraintX[i]+fFrontConstraintWidthX)
	  cout << "front X: " << fFrontConstraintX[i]-fFrontConstraintWidthX 
	       << " <? " << the_track->GetX(Z_cst) << " <? " 
	       << fFrontConstraintX[i]+fFrontConstraintWidthX << endl;
	
	if(fFrontConstraintY[i]-fFrontConstraintWidthY > the_track->GetY() ||
	   the_track->GetY() > fFrontConstraintY[i]+fFrontConstraintWidthY)
	  cout << "front Y: " << fFrontConstraintY[i]-fFrontConstraintWidthY 
	       << " <? " << the_track->GetY(Z_cst) << " <? " 
	       << fFrontConstraintY[i]+fFrontConstraintWidthY << endl;
	*/
	
	if(fBackConstraintX[i]-fBackConstraintWidthX < the_track->GetX(Z_cst) &&  
	   the_track->GetX(Z_cst) < fBackConstraintX[i]+fBackConstraintWidthX && 
	   fBackConstraintY[i]-fBackConstraintWidthY < the_track->GetY(Z_cst) &&  
	   the_track->GetY(Z_cst) < fBackConstraintY[i]+fBackConstraintWidthY ){
	  i_match = i;
	}
      }
      if(i_match<0)continue;
      
      fEtotPratio.push_back(fEtot[i_match]/the_track->GetP());
    }//end if(inheritsfrom(SBSBBTotalShower))

    
  }
  
}

//_____________________________________________________________________________
void SBSBigBite::CalcTrackPID(THaTrack* the_track)
{
  if(fEpsEtotRatio.size()==0 || fEtot.size()==0 || 
     fFrontConstraintX.size()==0 || fFrontConstraintY.size()==0 ||
     fBackConstraintX.size()==0 || fBackConstraintY.size()==0)return;
  
  //particles: 0: electron, 1: pion
  //detectors: 0: ps/sh
  THaPIDinfo* pidinfo = new THaPIDinfo(2, 2);
  pidinfo->SetDefaultPriors();
  
  the_track->SetPIDinfo(pidinfo);
  
  double eproba, piproba;
  
  TIter next( fNonTrackingDetectors );
  while( auto* theNonTrackDetector =
	 static_cast<THaNonTrackingDetector*>( next() )) {
    //first, Calorimeter PID
    //the calorimeter has to be found for anything to be done.
    if(theNonTrackDetector->InheritsFrom("SBSBBTotalShower")){
      SBSBBTotalShower* BBTotalShower = reinterpret_cast<SBSBBTotalShower*>(theNonTrackDetector);
      double Z_cst =  BBTotalShower->GetShower()->GetOrigin().Z();
      //check that the track we consider is consistent with the calorimeter constraint 
      int i_match = -1;
      for(int i = 0; i<fEtot.size(); i++){
	if(fBackConstraintX[i]-fBackConstraintWidthX < the_track->GetX(Z_cst) &&  
	   the_track->GetX(Z_cst) < fBackConstraintX[i]+fBackConstraintWidthX && 
	   fBackConstraintY[i]-fBackConstraintWidthY < the_track->GetY(Z_cst) &&  
	   the_track->GetY(Z_cst) < fBackConstraintY[i]+fBackConstraintWidthY ){
	  i_match = i;
	}
      }
      if(i_match<0)continue;
      
      double pr1, pr2;
      //fEtotPratio.push_back(fEtot[i_match]/the_track->GetP());
      proba_pssh(fEpsEtotRatio[i_match], eproba, piproba);
      proba_pcal(fEtot[i_match]/the_track->GetP(), pr1, pr2);
      eproba*=pr1;
      piproba*=pr2;
      the_track->GetPIDinfo()->SetProb(0, 0, eproba);
      the_track->GetPIDinfo()->SetProb(0, 1, piproba);
    }//end if(inheritsfrom(SBSBBTotalShower))
    
    // then, GRINCH PID
    // match a GRINCH cluster to a track:
    // again, the grinch has to be found for anything to be done.
    if(theNonTrackDetector->InheritsFrom("SBSGRINCH")){
      SBSGRINCH* GRINCH = reinterpret_cast<SBSGRINCH*>(theNonTrackDetector);
      
      //x, y of track at z = Z_GRINCH
      double x_track = the_track->GetX(GRINCH->GetZ());
      //double y_track = the_track->GetY()+
      //the_track->GetPhi()*GRINCH->GetOrigin().Z();
      
      //cout << "x, y track = " << x_track << ", " 
      //<< the_track->GetY()+the_track->GetPhi()*GRINCH->GetZ() << endl;
      
      //cout << GRINCH->GetNumClusters() << " GRINCH clusters " << endl;
      int NGRINCHPMTs_match = 0;
      
      for(int i = 0; i<GRINCH->GetNumClusters(); i++){
	SBSGRINCH_Cluster* gc_clus = GRINCH->GetCluster(i);
	//cout << "N hits = " << gc_clus->GetNHits() << endl;
	//cout << "x,y grinch " << gc_clus->GetXcenter() << ", "
	//   << gc_clus->GetYcenter() << endl;
	
	//cout << gc_clus->GetXcenter()*fTrackGrinchClusCorr_1+fTrackGrinchClusCorr_0 << " " <<  x_track << endl;
	//cout << fabs(x_track-(gc_clus->GetXcenter()*fTrackGrinchClusCorr_1+fTrackGrinchClusCorr_0)) << " " << fTrackGrinchClusCorr_Sigma << endl;
	
	if( fabs(x_track-gc_clus->GetXcenter()*fTrackGrinchClusCorr_1-fTrackGrinchClusCorr_0)<fTrackGrinchClusCorr_Sigma)NGRINCHPMTs_match = gc_clus->GetNHits();
      }
      
      proba_grinch(NGRINCHPMTs_match, the_track->GetP(), eproba, piproba);
      the_track->GetPIDinfo()->SetProb(1, 0, eproba);
      the_track->GetPIDinfo()->SetProb(1, 1, piproba);
    }
    
    the_track->GetPIDinfo()->CombinePID();
    
    // cout << " Eps/Etot = " << fEpsEtotRatio[the_track->GetIndex()] 
    // 	 << " Etot/p = " << fEtot[the_track->GetIndex()]/the_track->GetP()
    // 	 << " N GRINCH PMTs = " << NGRINCHPMTs_match 
    // 	 << ", P = " << the_track->GetP() << endl;
    // cout << " => combined track PID: electron " 
    // 	 << the_track->GetPIDinfo()->GetProb(0, 0) << " "
    // 	 << the_track->GetPIDinfo()->GetProb(1, 0) << " "
    // 	 << the_track->GetPIDinfo()->GetCombinedProb(0) 
    // 	 << " pion " << the_track->GetPIDinfo()->GetProb(0, 1) << " "
    // 	 << the_track->GetPIDinfo()->GetProb(1, 1) << " "
    // 	 << the_track->GetPIDinfo()->GetCombinedProb(1) << endl;
    
    fProbaE.push_back(the_track->GetPIDinfo()->GetCombinedProb(0));
    fProbaPi.push_back(the_track->GetPIDinfo()->GetCombinedProb(1));
    
  }
}


Int_t SBSBigBite::proba_pssh(Double_t eps_etot_ratio, 
			     Double_t& proba_e, Double_t& proba_pi)
{
  if(fProba_e_PSSH_table.size()==0)return -1;
  proba_e = fProba_e_PSSH_table[fProba_e_PSSH_table.size()-1];
  proba_pi = fProba_pi_PSSH_table[fProba_pi_PSSH_table.size()-1];
  for(size_t i = 0; i<fEpsEtotRatio_table.size()-1; i++){
    if(fEpsEtotRatio_table[i]<=eps_etot_ratio && 
       eps_etot_ratio<fEpsEtotRatio_table[i+1]){
      proba_e = fProba_e_PSSH_table[i]+
	(fProba_e_PSSH_table[i+1]-fProba_e_PSSH_table[i])/
	(fEpsEtotRatio_table[i+1]-fEpsEtotRatio_table[i])*
	(eps_etot_ratio-fEpsEtotRatio_table[i]);
      proba_pi = fProba_pi_PSSH_table[i]+
	(fProba_pi_PSSH_table[i+1]-fProba_pi_PSSH_table[i])/
	(fEpsEtotRatio_table[i+1]-fEpsEtotRatio_table[i])*
	(eps_etot_ratio-fEpsEtotRatio_table[i]);
    }
  }
  return 0;
}

Int_t SBSBigBite::proba_pcal(Double_t etot_p_ratio, 
			     Double_t& proba_e, Double_t& proba_pi)
{
  if(fEtotPratio_table.size()==0)return -1;
  proba_e = fProba_e_PCAL_table[fProba_e_PCAL_table.size()-1];
  proba_pi = fProba_e_PCAL_table[fProba_pi_PCAL_table.size()-1];
  for(size_t i = 0; i<fEtotPratio_table.size()-1; i++){
    if(fEtotPratio_table[i]<=etot_p_ratio && etot_p_ratio<fEtotPratio_table[i+1]){
      proba_e = fProba_e_PCAL_table[i]+
	(fProba_e_PCAL_table[i+1]-fProba_e_PCAL_table[i])/
	(fEtotPratio_table[i+1]-fEtotPratio_table[i])*
	(etot_p_ratio-fEtotPratio_table[i]);
      proba_pi = fProba_pi_PCAL_table[i]+
	(fProba_pi_PCAL_table[i+1]-fProba_pi_PCAL_table[i])/
	(fEtotPratio_table[i+1]-fEtotPratio_table[i])*
	(etot_p_ratio-fEtotPratio_table[i]);
    }
  }
  return 0;
}

Int_t SBSBigBite::proba_grinch(Int_t npmt, Double_t p,
			       Double_t& proba_e, Double_t& proba_pi)
{
  if(fProba_e_GRINCH_table.size()==0)return -1;
  int j = fP_table.size()-1;
  if(j==-1)return -1;
  for(size_t i = 0; i<fP_table.size()-1;i++){
    if(fP_table[i]<p && p<fP_table[i])j = i;
  }
  proba_e = fProba_e_GRINCH_table[fProba_e_GRINCH_table.size()-1];
  proba_pi = fProba_e_GRINCH_table[fProba_pi_GRINCH_table.size()-1];
  for(size_t i = 0; i<fNGRINCHPMTs_table.size()-1; i++){
    if(fNGRINCHPMTs_table[i]<=npmt && npmt<fNGRINCHPMTs_table[i+1]){
      proba_e = fProba_e_GRINCH_table[i]+
	(fProba_e_GRINCH_table[i+1]-fProba_e_GRINCH_table[i])/
	(fNGRINCHPMTs_table[i+1]-fNGRINCHPMTs_table[i])*
	(npmt-fNGRINCHPMTs_table[i]);
      proba_pi = fProba_pi_GRINCH_table[j][i]+
	(fProba_pi_GRINCH_table[j][i+1]-fProba_pi_GRINCH_table[j][i])/
	(fNGRINCHPMTs_table[i+1]-fNGRINCHPMTs_table[i])*
	(npmt-fNGRINCHPMTs_table[i]);
    }
  }
  return 0;
}

//_____________________________________________________________________________
Bool_t SBSBigBite::SetTrSorting( Bool_t set )
{
  Bool_t oldset = TestBit(kSortTracks);
  SetBit( kSortTracks, set );
  return oldset;
}

//_____________________________________________________________________________
Bool_t SBSBigBite::GetTrSorting() const
{
  return TestBit(kSortTracks);
}

//_____________________________________________________________________________
Bool_t SBSBigBite::SetMultiTracks( Bool_t set )
{
  Bool_t oldset = TestBit(kMultiTracks);
  SetBit( kMultiTracks, set );
  return oldset;
}

//_____________________________________________________________________________
Bool_t SBSBigBite::GetMultiTracks() const
{
  return TestBit(kMultiTracks);
}
