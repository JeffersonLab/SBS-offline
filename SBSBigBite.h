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
    
    virtual Int_t FindVertices( TClonesArray& tracks );
    virtual Int_t TrackCalc();
    
    //virtual Int_t   Begin( THaRunBase* r=0 );
    //virtual Int_t   End( THaRunBase* r=0 );

    protected:
    virtual Int_t ReadDatabase( const TDatime& date );
    virtual Int_t DefineVariables( EMode mode = kDefine );
    
    void CalcTargetCoords( THaTrack* the_track );
    void CalcTimingPID( THaTrack* the_track );
    
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
   
    
    /*
    TH2D* h1_yVx_bcp;
    TH2D* h1_x_fcpVbcp;
    TH2D* h1_yVx_fcp;
    TH2D* h1_y_fcpVbcp;
    TH2D* h1_dyVdx;
    */
    
    ClassDef(SBSBigBite,0) // BigBite spectrometer
};


#endif//ROOT_TreeSearch_SBSBigBite

