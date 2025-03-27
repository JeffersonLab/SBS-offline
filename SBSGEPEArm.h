#ifndef SBSGEPEArm_h
#define SBSGEPEArm_h

#include "THaSpectrometer.h"

class TList;
class THaTrack;

class SBSGEPEArm : public THaSpectrometer {

public:
  SBSGEPEArm( const char *name, const char *description );
  virtual ~SBSGEPEArm();


  virtual void  Clear( Option_t* opt="");
  
  virtual Int_t FindVertices( TClonesArray& tracks );
  virtual Int_t TrackCalc();
  
  virtual Int_t CoarseReconstruct();
  virtual Int_t	CoarseTrack();
  virtual Int_t	Reconstruct();
  virtual Int_t	Track();
  virtual Int_t CalcPID();

  inline Double_t GetECalDist() const { return fECALdist; };
  
  // void SetPolarimeterMode( Bool_t ispol );
  
protected:
  virtual Int_t ReadDatabase( const TDatime& date );
  virtual Int_t ReadRunDatabase( const TDatime& date );
  virtual Int_t DefineVariables( EMode mode = kDefine );
  
  //Reconstructed angles in transport coordinates:
  Double_t fECALtheta_n; //xECAL/ECALdist
  Double_t fECALphi_n; //yECAL/ECALdist

  //Reconstructed angles in global coordinates:
  Double_t fECALdir_x;
  Double_t fECALdir_y;
  Double_t fECALdir_z;

  //Distance to ECAL
  Double_t fECALdist; //add to run database (this only changes when kinematics change). But should it be optional or mandatory? 
  

  
  ClassDef(SBSGEPEArm,0) // BigBite spectrometer

};


#endif//ROOT_TreeSearch_SBSGEPEArm

