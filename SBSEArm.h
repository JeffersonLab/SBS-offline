#ifndef ROOT_TreeSearch_SBSEArm
#define ROOT_TreeSearch_SBSEArm

#include "THaSpectrometer.h"

class SBSEArm : public THaSpectrometer {

    public:
    SBSEArm( const char *name, const char *description );
    virtual ~SBSEArm();

    virtual Int_t FindVertices( TClonesArray& tracks );
    virtual Int_t TrackCalc();

    protected:
    ClassDef(SBSEArm,0) // BigBite spectrometer
};


#endif//ROOT_TreeSearch_SBSEArm

