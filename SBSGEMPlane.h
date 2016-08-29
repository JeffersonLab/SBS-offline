#include "THaSubDetector.h"

class THaDetectorBase;
class THaEvData;
class THaRunBase;

struct mpdmap_t {
    UInt_t crate;
    UInt_t slot;
    UInt_t mpd_id;
    UInt_t gem_id;
    UInt_t adc_id;
    UInt_t i2c;
    UInt_t pos;
    UInt_t invert;
};

class SBSGEMPlane : public THaSubDetector {
    public:

        SBSGEMPlane( const char *name, const char *description = "",
                THaDetectorBase* parent = 0 );

        virtual ~SBSGEMPlane();

        virtual void    Clear( Option_t* opt="" );
        virtual Int_t   Decode( const THaEvData& );
        virtual void    Print( Option_t* opt="" ) const;

        virtual Int_t   ReadDatabase(const TDatime& );

        virtual Int_t   Begin( THaRunBase* r=0 );
        virtual Int_t   End( THaRunBase* r=0 );

    private:
        std::vector<mpdmap_t>    fMPDmap;
        std::vector<Int_t>    fChanMapData;


        ClassDef(SBSGEMPlane,0)

};
