#ifndef SBSEarm_h
#define SBSEarm_h

#include "THaSpectrometer.h"

class TList;
class THaTrack;
class TVector3;
class TLorentzVector;
  
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

  //Override/extend THaAnalysisObject::Begin method:
  virtual Int_t   Begin( THaRunBase* r=0 );
  
  void SetPolarimeterMode( Bool_t ispol );

  //Change these methods to public access so other classes can use them:
  void CalcOpticsCoords( THaTrack* the_track );//calculate optics coords from det coords
  void CalcTargetCoords( THaTrack* the_track );//calculate target coords from optics coords  
  void CalcFpCoords( THaTrack* the_track ); //Calculate fp coords from

  //Add public methods to add diagnostic output variables for constraint points:
  // NOTE! These methods do not modify constraint points for SBSGEMTrackerBase or its derived classes,
  // that has to be done using the tracker methods. 
  void AddFrontConstraintPoint( TVector3 cpoint );
  void AddFrontConstraintPoint( double x, double y, double z );
  void AddBackConstraintPoint( TVector3 cpoint );
  void AddBackConstraintPoint( double x, double y, double z );

  Double_t GetFrontConstraintWidthX(int icp=0);
  Double_t GetFrontConstraintWidthY(int icp=0);
  Double_t GetBackConstraintWidthX(int icp=0);
  Double_t GetBackConstraintWidthY(int icp=0);
  Double_t GetFrontConstraintX0(int icp=0);
  Double_t GetFrontConstraintY0(int icp=0);
  Double_t GetBackConstraintX0(int icp=0);
  Double_t GetBackConstraintY0(int icp=0);

  Bool_t IsPolarimeter() const { return fPolarimeterMode; }
  
  
  
protected:
  virtual Int_t ReadDatabase( const TDatime& date );
  virtual Int_t ReadRunDatabase( const TDatime& date );
  virtual Int_t DefineVariables( EMode mode = kDefine );

  //target coords using forward optics model
  
  void InitOpticsAxes(double, const TVector3 & );
  void InitOpticsAxes(double); //version with only bend angle argument
  void InitGEMAxes(double, double, const TVector3 & );
  void InitGEMAxes(double, double); //version with only angle arguments:
  void InitGEMAxes(double, double, double, const TVector3 &); //version that takes three angles, consistent with more correct alignment procedure
  void InitGEMAxes(double, double, double);
  
  void CheckConstraintOffsetsAndWidths();
  
  //We have to make these vectors to accommodate the polarimeter mode; separate offsets and widths for front and back trackers:
  
  std::vector<Double_t> fFrontConstraintWidthX;
  std::vector<Double_t> fFrontConstraintWidthY;
  std::vector<Double_t> fBackConstraintWidthX;
  std::vector<Double_t> fBackConstraintWidthY;
  std::vector<Double_t> fFrontConstraintX0;
  std::vector<Double_t> fFrontConstraintY0;
  std::vector<Double_t> fBackConstraintX0; 
  std::vector<Double_t> fBackConstraintY0;

  //Idea is that slope of the track along x and y is roughly linearly correlated with the position of the back constraint:

  bool fUseDynamicConstraint; //The "dynamic constraint" sets the front constraint point automatically based on the back constraint point; it is useful for applying effective loose constraints based on spectrometer optics (correlation between x and theta, y and phi, etc)
  double fDynamicConstraintSlopeX;
  double fDynamicConstraintOffsetX;
  //double fDynamicWidthX; 
  double fDynamicConstraintSlopeY;
  double fDynamicConstraintOffsetY; 
  //double fDynamicWidthY;
  
  
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
  //Double_t fGEMtheta; //Polar angle of GEM stack Z axis relative to SBS Z axis
  //Double_t fGEMphi; //Azimuthal angle of GEM stack Z axis relative to SBS Z axis

  //X,Y,Z rotation angles (yaw,pitch,roll, resp.):
  Double_t fGEMax;
  Double_t fGEMay;
  Double_t fGEMaz;
  
  Double_t fMagDist; //mandatory parameter from run database
  Double_t fHCALdist; //add to run database (this only changes when kinematics change). But should it be optional or mandatory? 
  
  Double_t fBdL; //define BdL (assumed units = T*m)

  // TRotation fGEM_Rpos; //Total rotation to apply to space coordinates;
  // TRotation fGEM_Rdir; //Total rotation to apply to track direction
  
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

  Bool_t fUseBeamPosInOptics;
  
  int fOpticsOrder;
  int fOpticsNterms;
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

  Bool_t fPolarimeterMode; //Use polarimeter mode
  Bool_t fPolarimeterMode_DBoverride; //flag to override DB value
  
  Double_t fAnalyzerZ0; //Z of midpoint of analyzer. 
  Double_t fAnalyzerThick; //Total thickness of analyzer

  int fNbinsZBackTrackerConstraint; // we'll divide the analyzer thickness for the back polarimeter tracker into bins for the constraint definition
  
  //Also include (optional) forward optics model to aid in false track rejection. 
  int fForwardOpticsOrder;
  int fForwardOpticsNterms;
  
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

  Bool_t fGEPtrackingMode; //Flag to turn on GEP tracking mode. 

  
  
  ClassDef(SBSEArm,0) // BigBite spectrometer

};


#endif//ROOT_TreeSearch_SBSEArm

