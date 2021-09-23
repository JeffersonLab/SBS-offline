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
    virtual Int_t       CalcPID();
   
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
    virtual Int_t DefineVariables( EMode mode = kDefine );
    
    void CalcTargetCoords( THaTrack* the_track );
    void CalcTrackTiming( THaTrack* the_track );
    void CalcTrackPID( THaTrack* the_track );
    
    Int_t proba_pssh(Double_t eps_etot_ratio, 
		     Double_t& proba_e, Double_t& proba_pi);
    Int_t proba_pcal(Double_t etot_p_ratio, 
		     Double_t& proba_e, Double_t& proba_pi);
    Int_t proba_grinch(Int_t npmt, Double_t p, 
		       Double_t& proba_e, Double_t& proba_pi);
    
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

    Double_t fFrontConstraintWidthX;
    Double_t fFrontConstraintWidthY;
    Double_t fBackConstraintWidthX;
    Double_t fBackConstraintWidthY;

    //for output only... Vectors instead?
    std::vector<double> fFrontConstraintX;
    std::vector<double> fFrontConstraintY;
    std::vector<double> fBackConstraintX;
    std::vector<double> fBackConstraintY;
    
    std::vector<double> fEpsEtotRatio;
    std::vector<double> fEtot;
    
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
    
    double fTrackerPitchAngle;
    /*
    TH2D* h1_yVx_bcp;
    TH2D* h1_x_fcpVbcp;
    TH2D* h1_yVx_fcp;
    TH2D* h1_y_fcpVbcp;
    TH2D* h1_dyVdx;
    */

    enum {
      kMultiTracks  = BIT(13), // Tracks are to be sorted by chi2
      kSortTracks   = BIT(14), // Tracks are to be sorted by chi2
      kAutoStdDets  = BIT(15)  // Auto-create standard detectors if no "vdc"
    };
    
    ClassDef(SBSBigBite,0) // BigBite spectrometer
};


#endif//ROOT_TreeSearch_SBSBigBite

