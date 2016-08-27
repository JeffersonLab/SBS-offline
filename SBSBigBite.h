#ifndef ROOT_TreeSearch_SBSBigBite
#define ROOT_TreeSearch_SBSBigBite

#include "THaSpectrometer.h"

class SBSBigBite : public THaSpectrometer {

    public:
    SBSBigBite( const char *name, const char *description );
    virtual ~SBSBigBite();

    virtual Int_t FindVertices( TClonesArray& tracks );
    virtual Int_t TrackCalc();

    protected:
    ClassDef(SBSBigBite,0) // BigBite spectrometer
};


#endif//ROOT_TreeSearch_SBSBigBite

