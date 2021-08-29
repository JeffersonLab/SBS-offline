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

    virtual Int_t	CoarseReconstruct();
    virtual Int_t	CoarseTrack();
    virtual Int_t	Reconstruct();
    virtual Int_t	Track();
    
    virtual Int_t FindVertices( TClonesArray& tracks );
    virtual Int_t TrackCalc();
    
    /*
    // Class for storing matrix element data
    class THaMatrixElement {
    public:
    THaMatrixElement() : iszero(true), order(0), v(0)
	{ pw.reserve(5); poly.reserve(kPORDER); }
      bool match( const THaMatrixElement& rhs ) const
      { assert(pw.size() == rhs.pw.size()); return ( pw == rhs.pw ); }
      void clear()
      { iszero = true; pw.clear(); order = 0; v = 0.0; poly.clear(); }
      
      bool iszero;             // whether the element is zero
      std::vector<int> pw;     // exponents of matrix element
      //   e.g. D100 = { 1, 0, 0 }
      int  order;
      double v;                // its computed value
      std::vector<double> poly;// the associated polynomial
    };
    enum { kPORDER = 4 };

    std::vector<THaMatrixElement> fXptarMatrixElems;
    std::vector<THaMatrixElement> fYptarMatrixElems;
    std::vector<THaMatrixElement> fYtarMatrixElems;
    std::vector<THaMatrixElement> fPinvMatrixElems;
    std::vector<THaMatrixElement> fXtarMatrixElems;
    */
    //virtual Int_t   Begin( THaRunBase* r=0 );
    //virtual Int_t   End( THaRunBase* r=0 );

    protected:
    virtual Int_t ReadDatabase( const TDatime& date );
    virtual Int_t DefineVariables( EMode mode = kDefine );
    
    void CalcTargetCoords( THaTrack* the_track );
    
    int fOpticsOrder;
    std::vector<double> fb_xptar;
    std::vector<double> fb_yptar;
    std::vector<double> fb_ytar;
    std::vector<double> fb_pinv;
    
    Double_t fFrontConstraintWidthX;
    Double_t fFrontConstraintWidthY;
    Double_t fBackConstraintWidthX;
    Double_t fBackConstraintWidthY;

    //for output only
    Double_t fFrontConstraintX;
    Double_t fFrontConstraintY;
    Double_t fBackConstraintX;
    Double_t fBackConstraintY;
    
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

