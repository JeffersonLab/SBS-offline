//*-- Author :    Andrew Puckett 2025-03-23:

//////////////////////////////////////////////////////////////////////////
//
// SBSGEPRegionOfInterestModule
//
// Grab  ECAL cluster position (and possibly other information)
// Store it  
//
//////////////////////////////////////////////////////////////////////////

#include "SBSGEPRegionOfInterestModule.h"
#include "InterStageModule.h"
#include "THaGlobals.h"
#include "SBSGEPEArm.h" //For the electron arm:
#include "SBSEArm.h" //For the proton arm
#include "SBSECal.h"
#include "SBSHCal.h"
#include "SBSGEMSpectrometerTracker.h"
#include "SBSGEMPolarimeterTracker.h"
#include "TClonesArray.h"
#include "THaTrack.h"
#include "TMath.h"
#include "TList.h"
//_____________________________________________________________________________
SBSGEPRegionOfInterestModule::SBSGEPRegionOfInterestModule( const char *name, const char *description, Int_t stage ) : InterStageModule(name,description,stage){
  //Constructor; for now, does nothing other than instantiate
  //Default Earm and Parm names:
  fEarmName = "earm";
  fParmName = "sbs";
  fEarmDetName = "ecal";
  fParmDetName = "gemFT";
  fParmDetNamePol = "gemFPP";
  fParmDetNameCalo = "hcal";
  
  fTestTracks = new TClonesArray("THaTrack",1);

  fTargZ0 = 0.0;
  
  fDataValid = false; 
}
//_____________________________________________________________________________
SBSGEPRegionOfInterestModule::~SBSGEPRegionOfInterestModule()
{
  //Destructor; for now, does nothing except call THaAnalysisObject::RemoveVariables();
  RemoveVariables();

  delete fTestTracks;
};

//_____________________________________________________________________________
//Clear method: Invoke standard InterStageModule::Clear():
void SBSGEPRegionOfInterestModule::Clear( Option_t *opt )
{
  InterStageModule::Clear(opt);
  //Clear out any other event-level variables here:
  fECALclusterpos_global.SetXYZ(kBig,kBig,kBig);

  fECAL_energy = kBig;
  
  fetheta_central = kBig;
  fephi_central = kBig;
  fEprime_central = kBig;
  fptheta_central = kBig;
  fpphi_central = kBig;
  fPp_central = kBig;
  fECAL_energy = kBig;

  fTestTracks->Clear("C");
}

//_____________________________________________________________________________
Int_t SBSGEPRegionOfInterestModule::DefineVariables( THaAnalysisObject::EMode mode )
{
  // Define/delete event-by-event global variables

  Int_t ret = InterStageModule::DefineVariables(mode);// exports fDataValid etc.
  if( ret )
    return ret;

  RVarDef vars[] = {
    { "xECAL_global", "Global ECAL cluster x (m)", "fECALclusterpos_global.X()" },
    { "yECAL_global", "Global ECAL cluster y (m)", "fECALclusterpos_global.Y()" },
    { "zECAL_global", "Global ECAL cluster z (m)", "fECALclusterpos_global.Z()" },
    { "ECAL_energy", "ECAL best cluster energy (GeV)", "fECAL_energy" },
    { "etheta",  "electron polar angle (rad)", "fetheta_central" },
    { "ephi",  "electron azimuthal angle (rad)", "fephi_central" },
    { "Eprime", "electron expected energy (GeV)", "fEprime_central" },
    { "ptheta", "proton expected polar angle (rad)", "fptheta_central" },
    { "pphi", "proton expected azimuthal angle (rad)", "fpphi_central" },
    { "pp", "proton expected momentum (GeV/c)", "fPp_central" },
    { "xfp0", "predicted X at fp (assuming point target at origin)", "fxfp_central" },
    { "yfp0", "predicted Y at fp (assuming point target at origin)", "fyfp_central" },
    { "xpfp0", "predicted X' at fp (assuming point target at origin)", "fxpfp_central" },
    { "ypfp0", "predicted Y' at fp (assuming point target at origin)", "fypfp_central" },
    { nullptr }
  };

  return DefineVarsFromList( vars, mode );
}
//_____________________________________________________________________________
Int_t SBSGEPRegionOfInterestModule::ReadRunDatabase( const TDatime &date ){
  //Load beam energy:

  FILE* file = OpenRunDBFile( date );
  if( !file ) return kFileError;

  double ebeamtemp;
  
  const DBRequest req[] = {
    { "ebeam", &ebeamtemp, kDouble, 0, 0, 1 },
    { nullptr }
  };
  Int_t err = LoadDB( file, date, req );
  fclose(file);
  if( err )
    return kInitError;

  //We're neglecting electron mass here (for the purposes of this module it won't matter):
  fBeam4Vect.SetPxPyPzE( 0.0, 0.0, ebeamtemp, ebeamtemp);
  
  return kOK; 
}

Int_t SBSGEPRegionOfInterestModule::ReadDatabase( const TDatime &date ){
  //Load vertex z bin definitions, earm name and parm name:
  FILE* file = OpenFile( date );
  if( !file ){
    std::cerr << "SBSGEPRegionOfInterestModule::ReadDatabase(): database not found!"<< std::endl;
    return kFileError;
  }

  //load vertex z bins for scanning target length:
  const DBRequest request[] = {
    { "nbins_zvertex", &fNbinsVertexZ, kInt, 0, 1, 1 },
    { "zvertex_min", &fVertexZmin, kDouble, 0, 1, 1 },
    { "zvertex_max", &fVertexZmax, kDouble, 0, 1, 1 },
    { "earm_name", &fEarmName, kString, 0, 1, 1 },
    { "parm_name", &fParmName, kString, 0, 1, 1 },
    { "edet_name", &fEarmDetName, kString, 0, 1, 1 },
    { "pdet_name", &fParmDetName, kString, 0, 1, 1 },
    { "pdetpol_name", &fParmDetNamePol, kString, 0, 1, 1 },
    { "pdetcalo_name", &fParmDetNameCalo, kString, 0, 1, 1 },
    { "z0targ", &fTargZ0, kDouble, 0, 1, 1 },
    { nullptr }
  };
  
  Int_t status = LoadDB( file, date, request, fPrefix, 1 ); //The "1" after fPrefix means search up the tree
  fclose(file);
  if( status != 0 ){
    return status;
  }

  fIsInit = true;

  return kOK;
  
}

//_____________________________________________________________________________
Int_t SBSGEPRegionOfInterestModule::Process( const THaEvData &evdata ){
  //Okay here we go: we've written the code needed to start writing the code.

  THaApparatus *app = 0;

  bool gotEarm = false;
  bool gotParm = false;
  bool gotEdet = false;
  bool gotPdet = false;
  bool gotPdetPol = false;
  bool gotPdetCalo = false;
  
  TIter aiter(gHaApps);

  SBSGEPEArm *Earm = nullptr;
  SBSEArm *Parm = nullptr;

  SBSECal *Edet = nullptr;
  SBSGEMSpectrometerTracker *Pdet = nullptr;
  SBSGEMPolarimeterTracker *PdetPol = nullptr;

  SBSHCal *PdetCalo = nullptr;
  
  while( (app = (THaApparatus*) aiter()) ){
    std::string appname = app->GetName();
    if( app->InheritsFrom("SBSGEPEArm") ){
      if( appname == fEarmName ){
	Earm = dynamic_cast<SBSGEPEArm*>(app);
	gotEarm = true;

	Edet = dynamic_cast<SBSECal*>(Earm->GetDetector(fEarmDetName.c_str()));

	if( Edet ) gotEdet = true;
      }
    }
    if( app->InheritsFrom("SBSEArm") ){
      if( appname == fParmName ){
	Parm = dynamic_cast<SBSEArm*>(app);
	gotParm = true;

	Pdet = dynamic_cast<SBSGEMSpectrometerTracker*>(Parm->GetDetector(fParmDetName.c_str()));
	if( Pdet ) gotPdet = true;

	PdetPol = dynamic_cast<SBSGEMPolarimeterTracker*>(Parm->GetDetector(fParmDetNamePol.c_str()));

	if( PdetPol ) gotPdetPol = true;
	
	PdetCalo = dynamic_cast<SBSHCal*>(Parm->GetDetector(fParmDetNameCalo.c_str()));
	if( PdetCalo ) gotPdetCalo = true;
	
      }
    }
  }

  if( !gotParm || !gotEarm || !gotEdet || !gotPdet || !gotPdetPol || !gotPdetCalo ){
    std::cout << "Error: missing Earm and/or Parm and/or Edet and/or Pdet and/or PdetPol and/or PdetCalo! (gotEarm, gotParm, gotEdet, gotPdet, gotPdetPol, gotPdetCalo)=(" << gotEarm << ", " << gotParm << ", "
	      << gotEdet << ", " << gotPdet << ", " << gotPdetPol
	      << ", " << gotPdetCalo << ")" << std::endl;
    fDataValid = false;
    return 0;
  }

  //If we reached this point, then we can grab E arm and P arm info:

  //Always clear out proton arm front tracker constraint points before evaluating ROI:
  
  Pdet->ClearConstraints();
  PdetPol->ClearConstraints();

  //First grab HCAL info:
  if( PdetCalo->GetNclust() <= 0 ) return 0;

  auto HCalClusters = PdetCalo->GetClusters();

  int ibest_hcal = PdetCalo->GetBestClusterIndex();

  double xHCAL = HCalClusters[ibest_hcal]->GetX() + PdetCalo->GetOrigin().X();
  double yHCAL = HCalClusters[ibest_hcal]->GetY() + PdetCalo->GetOrigin().Y();
  double zHCAL = PdetCalo->GetOrigin().Z();
  
  double ThetaEarm = Earm->GetThetaGeo(); //E arm is ordinarily on beam left, so this angle SHOULD be positive
  double ThetaParm = Parm->GetThetaGeo(); //P arm is ordinarily on beam right, so this angle SHOULD be negative

  // The following lines assume ThetaEarm > 0 for beam left:
  TVector3 Earm_zaxis( sin(ThetaEarm), 0.0, cos(ThetaEarm) );
  TVector3 Earm_xaxis(0,-1,0); //TRANSPORT system; +x = down
  TVector3 Earm_yaxis = Earm_zaxis.Cross( Earm_xaxis ).Unit();
  
  //The following lines assume ThetaParm < 0 for beam right:
  TVector3 Parm_zaxis( sin(ThetaParm), 0.0, cos(ThetaParm) );
  TVector3 Parm_xaxis( 0, -1, 0 ); //TRANSPORT system; +x = down
  TVector3 Parm_yaxis = Parm_zaxis.Cross( Parm_xaxis ).Unit();

  //Set HCAL global cluster position. Not yet clear whether and/or how we will use this:
  fHCALclusterpos_global = Parm->GetHCALdist() * Parm_zaxis +
    HCalClusters[ibest_hcal]->GetX() * Parm_xaxis +
    HCalClusters[ibest_hcal]->GetY() * Parm_yaxis; 
  
  //Grab the E arm "track" 
  if( Earm->GetNTracks() >= 1 ){
    TClonesArray *EarmTracks = Earm->GetTracks();

    THaTrack *EarmTrack = ( (THaTrack*) (*EarmTracks)[0] );

    double xclust = EarmTrack->GetX();
    double yclust = EarmTrack->GetY();
    double ECALdist = Earm->GetECalDist();

    fECAL_energy = EarmTrack->GetEnergy();
    
    //TVector3 ECALpos(xclust,yclust,ECALdist);
    TVector3 ECALpos_global = xclust * Earm_xaxis + yclust * Earm_yaxis + ECALdist * Earm_zaxis;

    fECALclusterpos_global = ECALpos_global;

    Pdet->SetECALpos( ECALpos_global );
    
    TVector3 vertex_central(0,0,fTargZ0);
    
    // Central ECAL direction:
    TVector3 ECALdir_global = (ECALpos_global - vertex_central).Unit();

    double ebeam = fBeam4Vect.E();
    double Mp = fmass_proton_GeV;

    Pdet->SetBeamE( ebeam );
    
    fetheta_central = ECALdir_global.Theta();
    fephi_central = ECALdir_global.Phi();
    fEprime_central = ebeam/(1.0+ebeam/Mp*(1.0-cos(fetheta_central)));

    double Q2 = 2.0*ebeam*fEprime_central*(1.0-cos(fetheta_central));
    double tau = Q2/(4.0*Mp*Mp);
    fPp_central = sqrt(Q2*(1.0+tau)); // = sqrt(nu^2 + 2M nu)
    fptheta_central = acos( (ebeam-fEprime_central*cos(fetheta_central))/fPp_central );
    fpphi_central = fephi_central + TMath::Pi();

    // std::cout << "SBSGEMRegionOfInterestModule: multi tracks enabled = "
    // 	      << Pdet->MultiTracksEnabled() << std::endl;
    //Get constraint point offsets for centering:

    //Calculate "central" expected proton track regardless of whether we're using the z-vertex binning:

    // This calculation assumes a point target at the origin:
    TVector3 pnhat_central( sin(fptheta_central)*cos(fpphi_central),sin(fptheta_central)*sin(fpphi_central),cos(fptheta_central));
    TVector3 ProtonMomentum = fPp_central * pnhat_central;
    double raytemp[6];
    
    TVector3 vdummy(0,0,fTargZ0);
    TVector3 dummy;
    Parm->LabToTransport( vdummy, ProtonMomentum, dummy, raytemp );
    
    double xptar = raytemp[1];
    double yptar = raytemp[3];
    double xtar = raytemp[0];
    double ytar = raytemp[2];
    
    int itrack = fTestTracks->GetLast()+1;
    THaTrack *Ttemp = new( (*fTestTracks)[itrack] ) THaTrack();
    
    Ttemp->SetTarget( xtar, ytar, xptar, yptar );
    Ttemp->SetMomentum( fPp_central );
    
    Parm->CalcFpCoords( Ttemp );
    
    Ttemp->Set( Ttemp->GetDX(), Ttemp->GetDY(), Ttemp->GetDTheta(), Ttemp->GetDPhi() );
    
    double xfp = Ttemp->GetX();
    double yfp = Ttemp->GetY();
    double xpfp = Ttemp->GetTheta();
    double ypfp = Ttemp->GetPhi();
    
    //set output variables:
    fxfp_central = xfp;
    fyfp_central = yfp;
    fxpfp_central = xpfp;
    fypfp_central = ypfp;
    
    double x0fcp = Parm->GetFrontConstraintX0(0);
    double y0fcp = Parm->GetFrontConstraintY0(0);
    double x0bcp = Parm->GetBackConstraintX0(0);
    double y0bcp = Parm->GetBackConstraintY0(0);
    
    if( Pdet->MultiTracksEnabled() ){ //use multiple constraint points:

      //First clear out any existing constraint points for FT:
      
      TVector3 vertex;

      double zbinwidth = (fVertexZmax - fVertexZmin)/double(fNbinsVertexZ);
      //      std::cout << "GEP region of interest vertex scan:" << std::endl;
      for( int ibin=0; ibin<fNbinsVertexZ; ibin++ ){
	
	//Set vertex assumption:
	vertex.SetXYZ(0.0,0.0,fVertexZmin + (ibin+0.5)*zbinwidth);

	// std::cout << "(ibin, zvertex)=(" << ibin << ", " << vertex.Z() << ")"
	// 	  << std::endl;
	
	//Now calculate electron scattering angle and the rest of e and p kinematic variables:
	TVector3 enhat = (ECALpos_global-vertex).Unit();

	double etheta = enhat.Theta();
	double ephi = enhat.Phi();
	double eprime = ebeam/(1.0+ebeam/Mp*(1.0-cos(etheta)));
	Q2 = 2.0*ebeam*eprime*(1.0-cos(etheta));
	tau = Q2/(4.0*Mp*Mp);

	double pp = sqrt(Q2*(1.0+tau));
	double ptheta = acos( (ebeam-eprime*cos(etheta))/pp );
	double pphi = ephi + TMath::Pi();

	//Proton direction (unit vector):
	TVector3 pnhat(sin(ptheta)*cos(pphi),sin(ptheta)*sin(pphi),cos(ptheta));

	ProtonMomentum = pp*pnhat;
	//Next we need to calculate this in SBS (Parm) transport coordinates:

	//	double raytemp[6];

	TVector3 tvert; 
	
	Parm->LabToTransport( vertex, ProtonMomentum, tvert, raytemp );
	//	TVector3 pnhat_SBS( pnhat.Dot( Parm_xaxis ), pnhat.Dot( Parm_yaxis ), pnhat.Dot( Parm_zaxis ) );
	//double xptar_p = pnhat_SBS.X()/pnhat_SBS.Z();
	//double yptar_p = pnhat_SBS.Y()/pnhat_SBS.Z();

	xptar = raytemp[1];
	yptar = raytemp[3];
	xtar = raytemp[0];
	ytar = raytemp[2];

	//Now with these quantities calculated, we are able to use the forward optics matrix to predict the FP track:
	
	itrack = fTestTracks->GetLast() + 1;

	//Initialize with the default constructor (no arguments);
	Ttemp = new( (*fTestTracks)[itrack] ) THaTrack();

	Ttemp->SetTarget( xtar, ytar, xptar, yptar );
	Ttemp->SetMomentum( pp );

	//Track to focal plane:
	Parm->CalcFpCoords( Ttemp );

	//Although it doesn't matter that much, set the "regular" FP coordinates based on the "detector" coordinates:
	Ttemp->Set( Ttemp->GetDX(), Ttemp->GetDY(), Ttemp->GetDTheta(), Ttemp->GetDPhi() );
	//After the line above, the "det" coordinates are the same as the "regular" coordinates.
	//We could, of course, just grab the "det" coordinates directly:
	
	xfp = Ttemp->GetX();
	yfp = Ttemp->GetY();
	xpfp = Ttemp->GetTheta();
	ypfp = Ttemp->GetPhi();

	//Now we add front and back constraint points based on these calculated track parameters.
	// What Z value should we assume for the back constraint point? For the front it's easy.
	// For the back 

	double zfront = 0.0;
	double zback = Pdet->GetZmaxLayer() + 0.05; //the 0.05 here (5 cm past the back GEM layer) is arbitrary... may need adjustment
	
	Pdet->SetFrontConstraintPoint( xfp + x0fcp, yfp + y0fcp, 0.0 );
	Pdet->SetBackConstraintPoint( xfp + xpfp * zback + x0bcp, yfp + ypfp*zback + y0bcp, zback );
	
	//Also add constraint points to diagnostic outputs:
	Parm->AddFrontConstraintPoint( xfp + x0fcp, yfp + y0fcp, 0.0 );
	Parm->AddBackConstraintPoint( xfp + xpfp * zback + x0bcp, yfp + ypfp*zback + y0bcp, zback );
	
	fDataValid = true;
	//Now, IF the multi-track search is enabled (and if the code is written correctly), this should be sufficient
      } //end loop over z vertex bins 
    } else { // end if multi-tracks enabled)
      
      double zfront = 0.0;
      double zback = Pdet->GetZmaxLayer() + 0.05;
      
      Pdet->SetFrontConstraintPoint(fxfp_central + x0fcp, fyfp_central + y0fcp, 0.0 );
      Pdet->SetBackConstraintPoint( fxfp_central+fxpfp_central*zback + x0bcp, fyfp_central+fypfp_central*zback + y0bcp, zback );

      fDataValid = true;
      
    }
  }

  
  return 0;
}

ClassImp(SBSGEPRegionOfInterestModule);
