//*-- Author :    Andrew Puckett 2025-04-02:

//////////////////////////////////////////////////////////////////////////
//
// SBSGEPHeepCoinModule
//
// Calculate Exclusivity cut variables for "online" analysis
//
//////////////////////////////////////////////////////////////////////////

#include "SBSGEPHeepCoinModule.h"
#include "THaPhysicsModule.h"
#include "THaTrackingModule.h"
#include "SBSGEPEArm.h"
#include "SBSEArm.h"
#include "TClonesArray.h"
#include "THaTrack.h"
#include "SBSHCal.h"
#include "SBSECal.h"
#include "TList.h"
//_____________________________________________________________________________
SBSGEPHeepCoinModule::SBSGEPHeepCoinModule( const char *name, const char *desc, const char *espectro, const char *pspectro) : THaPhysicsModule(name, desc), fEarmName(espectro), fParmName(pspectro), fPspectro(nullptr), fEspectro(nullptr){

  //  fEarmNameespectro;
  //fParmName(psectro);
  fProtonMass = 0.93827208816; //PDG value as of 4/8/2025
  fTarget4vect.SetPxPyPzE( 0, 0, 0, fProtonMass );
}

//_____________________________________________________________________________
SBSGEPHeepCoinModule::~SBSGEPHeepCoinModule(){
  //standard destructor
  RemoveVariables();
}

void SBSGEPHeepCoinModule::Clear( Option_t *opt ){
  THaPhysicsModule::Clear(opt);

  fvertex.SetXYZ(kBig,kBig,kBig);

  fetheta = fephi = fPtheta = fPphi = fEcalo = fPp = kBig;
  fEprime_eth = fPp_eth = fPth_eth = fPph_eph = kBig;
  fPp_pth = fEprime_pth = feth_pth = feph_pph = kBig;
  fQ2_pp = fQ2_eth = fQ2_p4vect = fQ2_e4vect = kBig;

  fEprime_pp = feth_pp = kBig;
  
  fQ2_pth = kBig;
  fepsilon_eth = fepsilon_pth = fepsilon_pp = fepsilon_p4vect = kBig;
  
  fProton4vect.SetPxPyPzE( kBig, kBig, kBig, kBig );
  fElectron4vect.SetPxPyPzE( kBig, kBig, kBig, kBig );

  fDp_pth = fDp_eth = fdphi = facoplanarity = kBig;

  fdxECAL = fdyECAL = fdxECAL_4vect = fdyECAL_4vect = kBig; 

  fKinFact_eth = fKinFact_pth = fKinFact_pp = fKinFact_p4vect = kBig;
  
  fdeltat_ADC = fdeltat_TDC = kBig;
}

//_____________________________________________________________________________
Int_t SBSGEPHeepCoinModule::DefineVariables( EMode mode ){
  RVarDef vars[] = {
    { "datavalid", "data valid? T/F", "fDataValid" },
    { "etheta", "electron polar angle", "fetheta" },
    { "ephi", "electron azimuthal angle", "fephi" },
    { "ptheta", "proton polar angle", "fPtheta" },
    { "pphi", "proton azimuthal angle", "fPphi" },
    { "ecalo", "electron energy (GeV, from ECAL)", "fEcalo" },
    { "pp", "proton momentum (GeV/c)", "fPp" },
    { "eprime_eth", "electron energy from electron angle", "fEprime_eth" },
    { "pp_eth", "proton momentum from electron angle", "fPp_eth" },
    { "pth_eth", "proton angle from electron angle", "fPth_eth" },
    { "pp_pth", "proton momentum from proton angle", "fPp_pth" },
    { "eprime_pth", "electron energy from proton angle", "fEprime_pth" },
    { "eth_pth", "electron angle from proton angle", "feth_pth" },
    { "eprime_pp", "electron energy from proton momentum", "fEprime_pp" },
    { "eth_pp", "electron angle from proton momentum", "feth_pp" },
    { "Q2_pp", "Q^2 from proton momentum", "fQ2_pp" },
    { "Q2_eth", "Q^2 from electron angle", "fQ2_eth" },
    { "Q2_pth", "Q^2 from proton angle", "fQ2_pth" },
    { "Q2_p4vect", "Q^2 from proton four-vector", "fQ2_p4vect" },
    { "Q2_e4vect", "Q^2 from electron four-vector", "fQ2_e4vect" },
    { "eps_pp", "epsilon from proton momentum", "fepsilon_pp" },
    { "eps_eth", "epsilon from electron angle", "fepsilon_eth" },
    { "eps_pth", "epsilon from proton angle", "fepsilon_pth" },
    { "eps_4vect", "epsilon from proton four-vector", "fepsilon_p4vect" },
    { "K_eth", "Kinematic factor for mu GE/GM from electron angle", "fKinFact_eth" },
    { "K_pth", "Kinematic factor for mu GE/GM from proton angle", "fKinFact_pth" },
    { "K_pp", "Kinematic factor for mu GE/GM from proton momentum", "fKinFact_pp" },
    { "K_p4vect", "Kinematic factor for mu GE/GM from proton four-vector", "fKinFact_p4vect" },
    { "dpp", "pproton/pexpect(thetaproton)-1", "fDp_pth" },
    { "dpe", "pproton/pexpect(thetaelectron)-1", "fDp_eth" },
    { "dphi", "(phi_e - phi_p + PI)", "fdphi" },
    { "acoplanarity", "acoplanarity (angle between electron and proton planes)", "facoplanarity" },
    { "inelasticity", "inelasticity (Pbeam + Ptarg - Pproton).m2()", "finelasticity_proton" },
    { "dxECAL", "x ECAL - x expect (angles-only)", "fdxECAL" },
    { "dyECAL", "y ECAL - y expect (angles-only)", "fdyECAL" },
    { "dxECAL4vect", "x ECAL - x expect (four-vector)", "fdxECAL_4vect" },
    { "dyECAL4vect", "y ECAL - y expect (four-vector)", "fdyECAL_4vect" },
    { "dt_ADC", "tECAL-tHCAL (ADC time of seed block)", "fdeltat_ADC" },
    { "dt_TDC", "tCDET-tHCAL (TDC time of seed hit(s))", "fdeltat_TDC" }, 
    { nullptr }
  };

  return DefineVarsFromList( vars, mode );
}

//_____________________________________________________________________________

Int_t SBSGEPHeepCoinModule::ReadRunDatabase( const TDatime &date ){
  //Get beam energy (and perhaps spectrometer angles?) from run db. Actually only beam energy is really needed. Or IS IT? 
  //Load beam energy:

  FILE* file = OpenRunDBFile( date );
  if( !file ) return kFileError;

  double ebeamtemp = febeam;
  
  const DBRequest req[] = {
    { "ebeam", &ebeamtemp, kDouble, 0, 0, 1 },
    { nullptr }
  };
  Int_t err = LoadDB( file, date, req );
  fclose(file);
  if( err )
    return kInitError;

  //We're neglecting electron mass here (for the purposes of this module it won't matter):
  fBeam4vect.SetPxPyPzE( 0.0, 0.0, ebeamtemp, ebeamtemp);

  febeam = ebeamtemp;
  
  return kOK; 
}

//_____________________________________________________________________________

Int_t SBSGEPHeepCoinModule::ReadDatabase( const TDatime &date ){
  FILE* file = OpenFile( date );
  if( !file ){
    std::cerr << "SBSGEPRegionOfInterestModule::ReadDatabase(): database not found!"<< std::endl;
    return kFileError;
  }

  std::string earmname = fEarmName.Data();
  std::string parmname = fParmName.Data();
  
  //Load spectrometer names (I THINK that's all we need for the moment):
  const DBRequest request[] = {
    { "earm_name", &earmname, kString, 0, 1, 1 },
    { "parm_name", &parmname, kString, 0, 1, 1 },
    { nullptr }
  };

  Int_t status = LoadDB( file, date, request, fPrefix, 1 ); //The "1" after fPrefix means search up the tree
  fclose(file);
  if( status != 0 ){
    return status;
  }

  fEarmName = earmname.c_str();
  fParmName = parmname.c_str(); 
  
  fIsInit = true;

  return kOK;
  
}
//_____________________________________________________________________________
Int_t SBSGEPHeepCoinModule::Begin( THaRunBase *run ){ //At this point all the apparatuses have been initialized:
  //bool gotEarm = false;
  //bool gotParm = false;
  // these arne't needed!
 
  THaApparatus *app = 0;

  TIter aiter( gHaApps );

  //We create this "Begin" method to initialize spectrometer-local x,y,z axes, so they are only calculated ONCE (since they are
  // certainly constant on a per-run basis)
  while( (app = (THaApparatus*) aiter() ) ){
    TString appname = app->GetName();
    if( app->InheritsFrom( "SBSGEPEArm") && appname == fEarmName ){
      fEspectro = dynamic_cast<SBSGEPEArm*>(app);
      //  gotEarm = true;

      double theta = fEspectro->GetThetaGeo();
      fEarm_zaxis.SetXYZ( sin(theta), 0, cos(theta) );
      fEarm_xaxis.SetXYZ(  0, -1, 0 );
      fEarm_yaxis = fEarm_zaxis.Cross(fEarm_xaxis).Unit();
    }

    if( app->InheritsFrom( "SBSEArm" ) && appname == fParmName ){
      fPspectro = dynamic_cast<SBSEArm*>(app);
      //  gotParm = true;

      double theta = fPspectro->GetThetaGeo(); //should be negative for beam right:
      fParm_zaxis.SetXYZ( sin(theta), 0, cos(theta) );
      fParm_xaxis.SetXYZ( 0, -1, 0 );
      fParm_yaxis = fParm_zaxis.Cross(fParm_xaxis).Unit();
    }
  }

  return kOK;
};
//_____________________________________________________________________________
Int_t SBSGEPHeepCoinModule::Process( const THaEvData &evdata ){
  
  if( fEspectro == nullptr || fPspectro == nullptr ){
    fDataValid = false;
    return 0;
  }

  double ECALDIST = fEspectro->GetECalDist();
  
  if( fPspectro->GetNTracks() >= 1 && fEspectro->GetNTracks() >= 1 ){

    double &E = febeam; //short-hand;
    double &Mp = fProtonMass;
    
    THaTrack *EarmTrack = fEspectro->GetGoldenTrack();
    THaTrack *ParmTrack = fPspectro->GetGoldenTrack();

    //Now for the Earm track, at this point of the code, the "xfp/yfp" variables contain the cluster positions:
    fvertex.SetXYZ(0,0,ParmTrack->GetVertexZ());

    double xclust = EarmTrack->GetX();
    double yclust = EarmTrack->GetY();
    fEcalo = EarmTrack->GetEnergy();
    
    TVector3 clpos_local(xclust,yclust,ECALDIST);

    TVector3 clpos_global;
    //How this works:
    // th = xclust/ECALDIST
    // ph = yclust/ECALDIST
    // p = clpos_local.Mag() = sqrt( xclus^2 + yclus^2 + ECALDIST^2 ) = distance from origin
    // clpos_global will contain the result of "fToLabRot * p * nhat,
    // where nhat is the unit vector along clpos_local, when we pass the following arguments to "TransportToLab": 
    fEspectro->TransportToLab( clpos_local.Mag(), xclust/ECALDIST, yclust/ECALDIST, clpos_global );

    //Now we calculate updated tracks:
    TVector3 enhat_global = (clpos_global - fvertex).Unit(); 

    //This completes the electron scattering angle reconstruction in the lab frame (will later need to be adjusted to include CDET info):
    fetheta = enhat_global.Theta();
    fephi = enhat_global.Phi();
    
    //We should also update the E arm "track" to reflect the new angles:
    TVector3 vdummy;
    double raytemp[6];
    fEspectro->LabToTransport( fvertex, fEcalo * enhat_global, vdummy, raytemp );

    double extar = raytemp[0];
    double exptar = raytemp[1];
    double eytar = raytemp[2];
    double eyptar = raytemp[3];

    EarmTrack->SetTarget( extar, eytar, exptar, eyptar );
    EarmTrack->Set( xclust, yclust, exptar, eyptar );
    EarmTrack->SetPvect( fEcalo * enhat_global );
    EarmTrack->SetVertex( fvertex );

    fElectron4vect.SetPxPyPzE( EarmTrack->GetLabPx(), EarmTrack->GetLabPy(), EarmTrack->GetLabPz(), EarmTrack->GetEnergy() );
    
    //Now grab proton arm info:

    fPp = ParmTrack->GetP();
    fPtheta = acos( ParmTrack->GetLabPz()/fPp );
    fPphi = atan2( ParmTrack->GetLabPy(), ParmTrack->GetLabPx() ); //This will be in -PI < phi < PI. Let's ADD 2pi to negative phi angles:
    if( fPphi < 0 ) fPphi += 2.0*TMath::Pi(); //This way ephi - pphi will be centered at -PI and (ephi - pphi + PI) will be centered at zero

    //Now we can start calculating some exclusivity cut variables:

    fEprime_eth = E/(1.0+E/Mp*(1.-cos(fetheta)));
    fQ2_eth = 2.0*E*fEprime_eth*(1.-cos(fetheta));

    double tau_eth = fQ2_eth/(4.0*Mp*Mp);
    double nu_eth = E-fEprime_eth;
    fPp_eth = sqrt(fQ2_eth*(1.0+tau_eth));
    
    fPth_eth = acos( (E-fEprime_eth*cos(fetheta))/fPp_eth );
    fPph_eph = fephi + TMath::Pi();
    
    fPp_pth = 2.0*Mp*E*(Mp+E)*cos(fPtheta)/(pow(Mp,2)+2.0*Mp*E + pow(E*sin(fPtheta),2));

    //now going the other way around:
    // Proton kinetic energy from proton angle:
    double nu_pth = sqrt( pow( fPp_pth, 2 ) + pow(Mp,2) ) - Mp;
    fQ2_pth = 2.0*Mp*nu_pth;

    double tau_pth = fQ2_pth/(4.0*Mp*Mp);
    
    fEprime_pth = E - nu_pth;
    feth_pth = acos( 1.0 - fQ2_pth/(2.0*E*fEprime_pth) );
    feph_pph = fPphi - TMath::Pi();

    double nu_pp = sqrt(pow(Mp,2)+pow(fPp,2))-Mp;
    fQ2_pp = 2.0*Mp*nu_pp;
    fEprime_pp = E - nu_pp;
    
    feth_pp = acos( 1.0 - fQ2_pp/(2.0*E*fEprime_pp) );
    
    double tau_pp = fQ2_pp/pow(2.0*Mp,2);
    
    //now calculate epsilon values:
    fepsilon_eth = 1.0/(1.0 + 2.0*(1.0+tau_eth)*pow(tan(fetheta/2.0),2));
    fepsilon_pth = 1.0/(1.0 + 2.0*(1.0+tau_pth)*pow(tan(feth_pth/2.0),2));
    fepsilon_pp  = 1.0/(1.0 + 2.0*(1.0+tau_pp)*pow(tan(feth_pp/2.0),2));

    //Kin fact = sqrt(tau*(1+eps)/(2*eps))
    fKinFact_eth = sqrt(tau_eth*(1.0+fepsilon_eth)/(2.0*fepsilon_eth));
    fKinFact_pth = sqrt(tau_pth*(1.0+fepsilon_pth)/(2.0*fepsilon_pth));
    fKinFact_pp = sqrt(tau_pp*(1.0+fepsilon_pp)/(2.0*fepsilon_pp));
    
    fProton4vect.SetPxPyPzE( ParmTrack->GetLabPx(),
			     ParmTrack->GetLabPy(),
			     ParmTrack->GetLabPz(),
			     sqrt(pow(fPp,2)+pow(Mp,2)) );

    fQ2_e4vect = -(fBeam4vect-fElectron4vect).M2();
    // Pproton = Ptarget + q --> q = Pproton-Ptarget
    fQ2_p4vect = -(fProton4vect-fTarget4vect).M2();

    fDp_pth = fPp/fPp_pth-1.0;
    fDp_eth = fPp/fPp_eth-1.0;
    fdphi = fephi - fPphi + TMath::Pi();

    TVector3 enhat_react = (fBeam4vect.Vect().Unit().Cross( enhat_global )).Unit(); //normal vector to electron scattering plane, defined by direction vectors only, k x k' --> points vertically up for electrons on beam left
    TVector3 pnhat_react = (fProton4vect.Vect().Unit().Cross( fBeam4vect.Vect().Unit() ) ).Unit(); //normal vector to proton scattering plane, defined by direction vectors only, q x k --> points vertically up for protons on beam right

    facoplanarity = acos( enhat_react.Dot( pnhat_react ) ); //angle between lepton and proton scattering planes

    //This is the electron 4-vector predicted from the proton 4-vector:
    // beam + target = e + p --> e = beam + 
    TLorentzVector e4vect_p = fBeam4vect + fTarget4vect - fProton4vect; 
    finelasticity_proton = e4vect_p.M2(); //missing mass of beam + target - proton (should give electron mass squared ~= zero)

    //Now we just need to calculate epsilon and kinfact from proton 4-vector: how should we define?
    double tau_p4vect = fQ2_p4vect/pow(2.0*Mp,2);
    double etheta_p4vect = e4vect_p.Vect().Theta();
    fepsilon_p4vect = 1.0/(1.0+2.0*(1.0+tau_p4vect)*pow(tan(etheta_p4vect/2.0),2));
    fKinFact_p4vect = sqrt(tau_p4vect*(1.0+fepsilon_p4vect)/(2.0*fepsilon_p4vect));

    //Last thing to calculate is the "dx, dy" variables (and also the time coincidence stuff... eventually)
    // Angles only method uses proton scattering angle reconstruction:
    TVector3 evect_pth(sin(feth_pth)*cos(feph_pph),sin(feth_pth)*sin(feph_pph),cos(feth_pth));

    fEspectro->LabToTransport( fvertex, fEprime_pth*evect_pth, vdummy, raytemp );

    
    //Ray temp contains coordinates in target system at zspec = 0:
    fdxECAL = EarmTrack->GetX() - (raytemp[0]+raytemp[1]*ECALDIST);
    fdyECAL = EarmTrack->GetY() - (raytemp[2]+raytemp[3]*ECALDIST);

    //Now, use the proton 4-vector:

    fEspectro->LabToTransport( fvertex, e4vect_p.Vect(), vdummy, raytemp );

    fdxECAL_4vect = EarmTrack->GetX() - (raytemp[0]+raytemp[1]*ECALDIST);
    fdyECAL_4vect = EarmTrack->GetY() - (raytemp[2]+raytemp[3]*ECALDIST);

    //I guess we stick with the convention that our "diffs" are earm - parm:
    fdeltat_ADC = fEspectro->GetAtimeECAL() - fPspectro->GetAtimeHCAL(); 
    fdeltat_TDC = fEspectro->GetTDCtimeCDET() - fPspectro->GetTDCtimeHCAL();

    fDataValid = true;
    
  } //And with that, I think we're done!

  
  return 0; 
};

ClassImp(SBSGEPHeepCoinModule);
