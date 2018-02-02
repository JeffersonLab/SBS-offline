#ifndef ROOT_TreeSearch_SBSBigBite
#define ROOT_TreeSearch_SBSBigBite

#include "THaSpectrometer.h"

/*
  Plug in "SBSTrackInfo as defined in TreeSearch::SBSSpec" if needed
 */

class SBSBigBite : public THaSpectrometer {

    public:
    SBSBigBite( const char *name, const char *description, 
		UInt_t nsectors = 1, Bool_t cer = true, Bool_t sci = true, Bool_t cal = true);
    virtual ~SBSBigBite();

    virtual void      Clear( Option_t* opt="" );
    virtual EStatus   Init( const TDatime& run_time );
    virtual Int_t     FindVertices( TClonesArray& tracks );
    virtual Int_t     TrackCalc();
    
    //TClonesArray*     GetTrackInfo() { return fSBSTrackInfo; }
    
    //protected:
    
    //TClonesArray*     fSBSTrackInfo;   // SBS-specific per-track info

    //virtual Int_t     DefineVariables( EMode mode = kDefine );
    //virtual Int_t     ReadRunDatabase( const TDatime& date );

    protected:
    ClassDef(SBSBigBite,0) // BigBite spectrometer
};


#endif//ROOT_TreeSearch_SBSBigBite
