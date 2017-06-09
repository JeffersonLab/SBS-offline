#ifndef ROOT_THaBBTotalShower
#define ROOT_THaBBTotalShower

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// THaBBTotalShower                                                          //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

//------------------------------------------------------//
//
//	Debug Definitions
//	place this section below any other head files
//
//------------------------------------------------------//
#ifdef DEBUG_LEVEL
#   undef DEBUG_LEVEL
#endif

//	DEBUG_LEVEL;	
//	=0	or not define: no debug, full speed 
//	>=1	enable debug extra warning (suggested setting)
//	>=2	above + enable debug assert
//	>=3	above + enable debug extra info
//	>=4	above + massive info (in a for or while)
#define DEBUG_LEVEL   2

#include    "DebugDef.h"
//------------------------------------------------------//

#include "THaPidDetector.h"

class THaBBShower;

class THaBBTotalShower : public THaPidDetector {

public:
    THaBBTotalShower( const char* name, const char* description = "",
        THaApparatus* a = NULL );
    THaBBTotalShower( const char* name, const char* shower_name,
        const char* preshower_name, const char* description = "",
        THaApparatus* a = NULL );
    virtual ~THaBBTotalShower();

    virtual Int_t      Decode( const THaEvData& );
    virtual Int_t      CoarseProcess( TClonesArray& tracks );
    virtual Int_t      FineProcess( TClonesArray& tracks );
    Float_t    GetE(Int_t i) const           { return fE[i]; }
    Float_t    GetY(Int_t i) const           { return fY[i]; }
    Float_t    GetX(Int_t i) const           { return fX[i]; }
    Int_t      GetID(Int_t i) const          { return fID[i]; }
    THaBBShower* GetShower() const      { return fShower; }
    THaBBShower* GetPreShower() const   { return fPreShower; }
    virtual EStatus    Init( const TDatime& run_time );
    virtual void       SetApparatus( THaApparatus* );
    void               LoadMCHitAt( Double_t x, Double_t y, Double_t E );
protected:

    // Maximum number of clusters
    static const Int_t kMaxNClust = 16;

    // Subdetectors
    THaBBShower* fShower;      // Shower subdetector
    THaBBShower* fPreShower;   // Preshower subdetector

    // Parameters
    Float_t    fMaxDx;       // Maximum dx between shower and preshower centers
    Float_t    fMaxDy;       // Maximum dx between shower and preshower centers

    // Per event data
    Int_t       fNclust;      // Number of clusters
    Float_t*    fE;           //[fNClust] Total shower energy
    Float_t*    fX;           //[fNClust] Total shower X
    Float_t*    fY;           //[fNClust] Total shower Y
    Int_t*      fID;          //[fNClust] ID of Presh and Shower coincidence

    void           ClearEvent();
    virtual Int_t  ReadDatabase( const TDatime& date );
    virtual Int_t  DefineVariables( EMode mode = kDefine );

private:
    void           Setup( const char* name,  const char* desc, 
        const char* shnam, const char* psnam,
        THaApparatus* app, bool mode );

    ClassDef(THaBBTotalShower,0)    //A total shower detector (shower plus preshower)
};

///////////////////////////////////////////////////////////////////////////////

#endif
