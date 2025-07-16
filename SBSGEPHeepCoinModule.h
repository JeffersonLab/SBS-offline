#ifndef SBSGEPHeepCoinModule_h_
#define SBSGEPHeepCoinModule_h_

///////////////////////////////////////////////////////////////////////////////////////////////////////
//
// SBSGEPHeepCoinModule: Utility physics module to calculate GEP elastic ep kinematic correlation stuff
//
///////////////////////////////////////////////////////////////////////////////////////////////////////

#include "THaPhysicsModule.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TString.h"

//We're going to store the pointers to Earm and Parm directly here:
class SBSEArm; class SBSGEPEArm;

class SBSGEPHeepCoinModule : public THaPhysicsModule {
  
public:
  SBSGEPHeepCoinModule( const char *name, const char *description, const char *espectro,
			const char *pspectro );

  virtual ~SBSGEPHeepCoinModule();

  virtual void Clear( Option_t *opt="" );
  virtual Int_t Process( const THaEvData& );

  

protected:

  virtual Int_t DefineVariables( EMode mode = kDefine );
  virtual Int_t ReadRunDatabase( const TDatime& date ); //need this at minimum to define beam energy
  //Not yet clear if a "ReadDatabase" method is actually needed
  virtual Int_t ReadDatabase( const TDatime &date );
  virtual Int_t Begin( THaRunBase *r );
  
  //Reconstructed particle scattering angles:
  TVector3 fvertex;  //reconstructed scattering vertex:

  Double_t febeam; //grab from run db as before
  
  Double_t fetheta, fephi, fPtheta, fPphi; 
  Double_t fEcalo, fPp; // reconstructed energy of electron and momentum of proton

  //Various derived quantities:
  // Electron energy and proton scattering angles and momentum from electron angles and beam energy: 
  Double_t fEprime_eth, fPp_eth, fPth_eth, fPph_eph;
  Double_t fpthtar_e, fpphtar_e;
  
  // When using the proton kinematics to predict electron kinematics, we pretty much want to use
  // two methods:
  // 1. angles-only
  // 2. four-vector
  // For some variables we could use a "momentum-only" prediction, but given the relatively poor momentum resolution of SBS,
  // these are less useful.

  // Proton momentum predicted from proton scattering angle:
  Double_t fPp_pth;

  // Predicted electron kinematics from proton scattering angle ONLY:  
  Double_t fEprime_pth, feth_pth, feph_pph;
  Double_t feth_pp, fEprime_pp;
  
  Double_t fQ2_pp;
  Double_t fQ2_eth;
  Double_t fQ2_pth; 
  Double_t fQ2_p4vect;
  Double_t fQ2_e4vect;

  Double_t fepsilon_eth;
  Double_t fepsilon_pth;
  Double_t fepsilon_pp;
  Double_t fepsilon_p4vect;

  Double_t fKinFact_eth;
  Double_t fKinFact_pth;
  Double_t fKinFact_pp;
  Double_t fKinFact_p4vect;
  
  TLorentzVector fProton4vect;
  TLorentzVector fElectron4vect; //Warning! Includes calorimeter energy resolution! 
  TLorentzVector fBeam4vect; 
  TLorentzVector fTarget4vect; 
  //Exclusivity cut variables:

  Double_t fDp_pth, fDp_eth; // equivalents of "pmissp" and "pmisse" from GEP-III analysis, ("central" momentum definition TBD).
  Double_t fdphi; //
  Double_t facoplanarity;

  Double_t finelasticity_proton; //(beam + target - proton)^2 (four vectors)
  
  // Basically since we only have tracking for the proton, we want to calculate dxECAL, dyECAL using two methods:
  // Using the proton four-vector, and "angles only":

  Double_t fdxECAL, fdyECAL; //"angles only" method:
  Double_t fdxECAL_4vect, fdyECAL_4vect;

  Double_t fdeltat_ADC; //Time difference between ECAL and HCAL cluster ADC times
  Double_t fdeltat_TDC; //CDET - HCAL TDC (to be implemented)
  
  //We may choose to add others later...
  TString fEarmName; //electron spectrometer name (default = "earm");
  TString fParmName; //proton spectrometer name (default = "sbs");

  //Pointers to spectrometer objects (technically the parent classes)
  SBSEArm *fPspectro; 
  SBSGEPEArm *fEspectro;

  Double_t fProtonMass;

  TVector3 fEarm_zaxis, fEarm_xaxis, fEarm_yaxis;
  TVector3 fParm_zaxis, fParm_xaxis, fParm_yaxis;
    
  ClassDef(SBSGEPHeepCoinModule,0)
};

#endif
