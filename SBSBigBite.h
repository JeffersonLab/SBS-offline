#ifndef ROOT_TreeSearch_SBSBigBite
#define ROOT_TreeSearch_SBSBigBite

#include "THaSpectrometer.h"

class TList;
class THaTrack;

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
    */
    protected:
    virtual Int_t ReadDatabase( const TDatime& date );
    /*
    std::vector<THaMatrixElement> fXptarMatrixElems;
    std::vector<THaMatrixElement> fYptarMatrixElems;
    std::vector<THaMatrixElement> fYtarMatrixElems;
    std::vector<THaMatrixElement> fPinvMatrixElems;
    std::vector<THaMatrixElement> fXtarMatrixElems;
    */
    
    ClassDef(SBSBigBite,0) // BigBite spectrometer
};


#endif//ROOT_TreeSearch_SBSBigBite

