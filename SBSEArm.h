#ifndef SBSEarm_h
#define SBSEarm_h

#include "THaSpectrometer.h"

class TList;
class THaTrack;

class SBSEArm : public THaSpectrometer {

public:
  SBSEArm( const char *name, const char *description );
  virtual ~SBSEArm();


  virtual void  Clear( Option_t* opt="");
  
  virtual Int_t FindVertices( TClonesArray& tracks );
  virtual Int_t TrackCalc();
  
  virtual Int_t CoarseReconstruct();
  virtual Int_t	CoarseTrack();
  virtual Int_t	Reconstruct();
  virtual Int_t	Track();
  virtual Int_t CalcPID();
  
protected:
  virtual Int_t ReadDatabase( const TDatime& date );
  virtual Int_t ReadRunDatabase( const TDatime& date );
  virtual Int_t DefineVariables( EMode mode = kDefine );

  void CalcOpticsCoords( THaTrack* the_track );//calculate optics coords from det coords
  void CalcTargetCoords( THaTrack* the_track );//calculate target coords from optics coords
  
  void InitOpticsAxes(double, const TVector3 & );
  void InitOpticsAxes(double); //version with only bend angle argument
  void InitGEMAxes(double, double, const TVector3 & );
  void InitGEMAxes(double, double); //version with only angle arguments:

  Double_t fFrontConstraintWidthX;
  Double_t fFrontConstraintWidthY;
  Double_t fBackConstraintWidthX;
  Double_t fBackConstraintWidthY;
  Double_t fFrontConstraintX0;
  Double_t fFrontConstraintY0;
  Double_t fBackConstraintX0; 
  Double_t fBackConstraintY0;

  //for output only... Vectors instead?
  std::vector<double> fFrontConstraintX;
  std::vector<double> fFrontConstraintY;
  std::vector<double> fFrontConstraintZ;
  std::vector<double> fBackConstraintX;
  std::vector<double> fBackConstraintY;
  std::vector<double> fBackConstraintZ;

  Double_t fHCALtheta_n; //xHCAL/HCALdist
  Double_t fHCALphi_n; //yHCAL/HCALdist

  Double_t fHCALdir_x;
  Double_t fHCALdir_y;
  Double_t fHCALdir_z;

  TVector3 fGEMorigin;  //Absolute position of GEM origin relative to target center, in TARGET transport coordinates
  Double_t fGEMtheta; //Polar angle of GEM stack Z axis relative to SBS Z axis
  Double_t fGEMphi; //Azimuthal angle of GEM stack Z axis relative to SBS Z axis

  Double_t fMagDist; //mandatory parameter from run database
  Double_t fHCALdist; //add to run database (this only changes when kinematics change). But should it be optional or mandatory? 
  
  Double_t fBdL; //define BdL (assumed units = T*m)

  TVector3 fGEMxaxis_global;
  TVector3 fGEMyaxis_global;
  TVector3 fGEMzaxis_global;

  TVector3 fOpticsOrigin; //Give origin of ideal optics system
  double fOpticsAngle; //Ideal central bend angle of GEM wrt BigBite
  TVector3 fOpticsXaxis_global;
  TVector3 fOpticsYaxis_global;
  TVector3 fOpticsZaxis_global;
  
  //TRotation fOpt2DetRot;// transformation from optics (ideal) to detector (actual)
  //TRotation fDet2OptRot;// transformation from detector (actual) to optics (ideal)

  UInt_t fPrecon_flag; //Indicate which momentum reconstruction formalism we are using:

  int fOpticsOrder;
  std::vector<double> fb_xptar;
  std::vector<double> fb_yptar;
  std::vector<double> fb_ytar;
  std::vector<double> fb_pinv;
  //AJRP: changed the exponents to integers here for speed:
  std::vector<int> f_oi;
  std::vector<int> f_oj;
  std::vector<int> f_ok;
  std::vector<int> f_ol;
  std::vector<int> f_om;

  //Also include (optional) forward optics model to aid in false track rejection. 
  int fForwardOpticsOrder;
  std::vector<double> fb_xfp;
  std::vector<double> fb_yfp;
  std::vector<double> fb_xpfp;
  std::vector<double> fb_ypfp;
  //AJRP: changed the exponents to integers here for speed:
  std::vector<int> f_foi;
  std::vector<int> f_foj;
  std::vector<int> f_fok;
  std::vector<int> f_fol;
  std::vector<int> f_fom;
  
  ClassDef(SBSEArm,0) // BigBite spectrometer

};


#endif//ROOT_TreeSearch_SBSEArm

