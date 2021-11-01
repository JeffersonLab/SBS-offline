#ifndef ROOT_TreeSearch_SBSBigBite
#define ROOT_TreeSearch_SBSBigBite

#include "THaSpectrometer.h"

class TList;
class THaTrack;
class TH2D;

class SBSBigBite : public THaSpectrometer {

public:
  SBSBigBite( const char *name, const char *description );
  virtual ~SBSBigBite();
    
  virtual void             Clear( Option_t* opt="");
  
  virtual Int_t	CoarseReconstruct();
  virtual Int_t	CoarseTrack();
  virtual Int_t	Reconstruct();
  virtual Int_t	Track();
  virtual Int_t CalcPID();
   
  virtual Int_t FindVertices( TClonesArray& tracks );
  virtual Int_t TrackCalc();
    
  //copied from THaHRS...
  Bool_t GetTrSorting() const;
  Bool_t SetTrSorting( Bool_t set = false );

  Bool_t GetMultiTracks() const;
  Bool_t SetMultiTracks( Bool_t set = false );
    
  //virtual Int_t   Begin( THaRunBase* r=0 );
  //virtual Int_t   End( THaRunBase* r=0 );

protected:
  virtual Int_t ReadDatabase( const TDatime& date );
  virtual Int_t ReadRunDatabase( const TDatime& date );
  virtual Int_t DefineVariables( EMode mode = kDefine );
  virtual void  DefinePidParticles();

  void CalcOpticsCoords( THaTrack* the_track );//calculate optics coords from det coords
  void CalcTargetCoords( THaTrack* the_track );//calculate target coords from optics coords
  void CalcTrackTiming( THaTrack* the_track );
  void CalcTrackPID( THaTrack* the_track );
    
  Int_t proba_pssh(Double_t eps_etot_ratio, 
		   Double_t& proba_e, Double_t& proba_pi);
  Int_t proba_pcal(Double_t etot_p_ratio, 
		   Double_t& proba_e, Double_t& proba_pi);
  Int_t proba_grinch(Int_t npmt, Double_t p, 
		     Double_t& proba_e, Double_t& proba_pi);

  void InitOpticsAxes(double, const TVector3 & );
  void InitOpticsAxes(double); //version with only bend angle argument
  void InitGEMAxes(double, double, const TVector3 & );
  void InitGEMAxes(double, double); //version with only angle arguments:
  
  // My current understanding is that fPointingOffset designates the point 
  // towards the beamline that the "mouth" of the spectrometer is pointing to...
  // which we might also need to define. 
  // Hence, we might want to define an additional set of parameters:
  // detector stack pitch, yaw, roll + actual position of first GEM tracker?
  // These angles are wrt the "ideal" central ray
  // double fDetectorStackPitch;
  // double fDetectorStackYaw;
  // double fDetectorStackRoll;
  
  // AJRP: reworking these parameters to be more consistent with the zero-field alignment code:
  Double_t fMagDist; //We could use "colldist" from THaSpectrometer for this, but that would be needlessly confusing
  
  Double_t fGEMtheta;   //polar angle of GEM z axis wrt TARGET transport coordinates
  Double_t fGEMphi;     //azimuthal angle of GEM z axis wrt TARGET transport coordinates
  TVector3 fGEMorigin;  //Absolute position of GEM origin relative to target center, in TARGET transport coordinates

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
    
  Double_t fPtheta_00000;
  Double_t fPtheta_10000;
  Double_t fPtheta_00100;
  Double_t fXptar_10000;
  Double_t fXptar_00100;
  Double_t fYtar_01000;
  Double_t fYtar_00010;

  std::vector<double> f_xtg_exp;


  Double_t fFrontConstraintWidthX;
  Double_t fFrontConstraintWidthY;
  Double_t fBackConstraintWidthX;
  Double_t fBackConstraintWidthY;
  //Let's add some handy-dandy offsets for the front point only,
  //so we can center the peaks at zero, tighten up the windows, 
  //and get tracking to run faster until things are better calibrated:
  Double_t fFrontConstraintX0;
  Double_t fFrontConstraintY0; 

    
  //for output only... Vectors instead?
  std::vector<double> fFrontConstraintX;
  std::vector<double> fFrontConstraintY;
  std::vector<double> fBackConstraintX;
  std::vector<double> fBackConstraintY;
    
  std::vector<double> fEpsEtotRatio;
  std::vector<double> fEtot;
  std::vector<double> fEtotPratio;
    
  double fTrackGrinchClusCorr_0;
  double fTrackGrinchClusCorr_1;
  double fTrackGrinchClusCorr_Sigma;
    
  std::vector<double> fEpsEtotRatio_table;
  std::vector<double> fProba_e_PSSH_table;
  std::vector<double> fProba_pi_PSSH_table;
    
  std::vector<double> fEtotPratio_table;
  std::vector<double> fProba_e_PCAL_table;
  std::vector<double> fProba_pi_PCAL_table;
    
  std::vector<double> fP_table;
  std::vector<double> fNGRINCHPMTs_table;
  std::vector<double> fProba_e_GRINCH_table;
  std::vector<std::vector<double>> fProba_pi_GRINCH_table;
    
  std::vector<double> fProbaE;
  std::vector<double> fProbaPi;

  //This is now redundant with fOpticsAngle
  //  double fTrackerPitchAngle;
    
  TH2D* h1_yVx_bcp;
  TH2D* h1_x_fcpVbcp;
  TH2D* h1_yVx_fcp;
  TH2D* h1_y_fcpVbcp;
  TH2D* h1_dyVdx;
  
    
  double fECaloFudgeFactor;// poor man's solution to apply the calorimeter constraint 
    
  enum {
    kMultiTracks  = BIT(13), // Tracks are to be sorted by chi2
    kSortTracks   = BIT(14), // Tracks are to be sorted by chi2
    kAutoStdDets  = BIT(15)  // Auto-create standard detectors if no "vdc"
  };
    
  ClassDef(SBSBigBite,0) // BigBite spectrometer
};


#endif//ROOT_TreeSearch_SBSBigBite

