#ifndef SBSGEPRegionOfInterestModule_h_
#define SBSGEPRegionOfInterestModule_h_

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// This is adapted from Podd_TimeCorrectionModule example:
//
// Its role is to grab cluster position (and possibly other variables) from SBSGEPEArm
// after the CoarseReconstruct stage for all detectors.
// It then populates a list of front and back constraint points and widths for the GEP
// front tracker
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "InterStageModule.h"
#include "TVector3.h"
#include "TLorentzVector.h"

class TClonesArray;
class THaTrack;
//class InterStageModule;

using namespace Podd;

class SBSGEPRegionOfInterestModule : public InterStageModule {
public:
  SBSGEPRegionOfInterestModule( const char *name, const char *description, Int_t stage );
  virtual ~SBSGEPRegionOfInterestModule();

  virtual void Clear( Option_t *opt="" );
  virtual Int_t Process( const THaEvData & );

  //I don't think we need a custom Init method here:
  //virtual EStatus Init( const TDatime& date );

  Double_t GetXfpCentral() const { return fxfp_central; }
  Double_t GetYfpCentral() const { return fyfp_central; }
  Double_t GetThfpCentral() const { return fxpfp_central; }
  Double_t GetPhfpCentral() const { return fypfp_central; }
  
protected:

  virtual Int_t  DefineVariables( EMode mode = kDefine );
  virtual Int_t  ReadDatabase( const TDatime& date );
  virtual Int_t  ReadRunDatabase( const TDatime& date ); //This is to load beam energy (redundant, I know, but whatever)

  //no need for this yet
  //  virtual Int_t   Begin( THaRunBase* r=0 );
  
  //constant (per-run) parameters:
  //-----------------------------------------------------------------------------------------------------------------------
  //This should be loaded from gHaRun->GetParameters()->GetBeamE(); alternatively we could just load it from the database?  
  TLorentzVector fBeam4Vect; //Where is the most convenient place to get the beam energy from? --> Run database 

  const double fmass_proton_GeV = 0.93827208816;
  //Define z vertex bins for scanning target extent (these should be loaded from regular DB):
  Int_t fNbinsVertexZ;
  Double_t fVertexZmin;
  Double_t fVertexZmax;

  Double_t fTargZ0;
  
  // Names of Earm and Parm: read from DB; I don't have a strong preference for
  // how to store these; might as well use std::string 
  std::string fEarmName;
  std::string fParmName;

  std::string fEarmDetName;
  std::string fParmDetName;

  std::string fParmDetNamePol;
  
  //We might as well store spectrometer 3-vectors here, or would that be redundant with the ones in the spectrometer classes? 

  //variable (per-event) parameters:
  //------------------------------------------------------------------------------------------------------------------------------
  TVector3 fECALclusterpos_global; //ECAL cluster position in "global" Hall A Coordinates (+x to beam left, +y up, +z along beam) 

  Double_t fECAL_energy;
  
  //Electron and proton final-state kinematics from ECAL cluster pos for point-target assumption:
  Double_t fetheta_central;
  Double_t fephi_central;
  Double_t fEprime_central;
  Double_t fptheta_central;
  Double_t fpphi_central;
  Double_t fPp_central;

  //Add variables to define the "central" elastically scattered proton ray (assuming point target at the origin):
  Double_t fxfp_central;
  Double_t fyfp_central;
  Double_t fxpfp_central;
  Double_t fypfp_central;
  
  TClonesArray *fTestTracks;
  
  // We may want to add some CDET-related info here once we understand what's going on there:
  
  
  ClassDef(SBSGEPRegionOfInterestModule,0)
};

#endif
