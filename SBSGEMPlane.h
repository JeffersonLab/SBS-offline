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

#define N_APV25_CHAN    128
#define N_MPD_TIME_SAMP 6
#define MPDMAP_ROW_SIZE 8


class SBSGEMPlane : public THaSubDetector {
    public:

        SBSGEMPlane( const char *name, const char *description = "",
		     THaDetectorBase* parent = 0, bool ismc = false );

        virtual ~SBSGEMPlane();

        virtual void    Clear( Option_t* opt="" );
        virtual Int_t   Decode( const THaEvData& );
        virtual void    Print( Option_t* opt="" ) const;

        virtual Int_t   ReadDatabase(const TDatime& );
        virtual Int_t   DefineVariables( EMode mode );

        virtual Int_t   Begin( THaRunBase* r=0 );
        virtual Int_t   End( THaRunBase* r=0 );

    private:
        std::vector<mpdmap_t>    fMPDmap;
        std::vector<Int_t>       fChanMapData;

        Double_t fZeroSuppressRMS;
        Bool_t fZeroSuppress;

        // Output variables
        Int_t  fNch;
        Int_t  fNchEff;
        Int_t *fStrip; // [fNch]
        Int_t *fadc[N_MPD_TIME_SAMP]; 
        // Being obnoxious so we match the stand alone more closely
        Int_t *fadc0; // [fNch]
        Int_t *fadc1; // [fNch]
        Int_t *fadc2; // [fNch]
        Int_t *fadc3; // [fNch]
        Int_t *fadc4; // [fNch]
        Int_t *fadc5; // [fNch]
	std::vector<Int_t> fStrip_; 
	std::vector<Int_t> fadc_0; 
	std::vector<Int_t> fadc_1;
	std::vector<Int_t> fadc_2;
	std::vector<Int_t> fadc_3;
	std::vector<Int_t> fadc_4;
	std::vector<Int_t> fadc_5;
        Double_t *fPedestal; 
        Double_t *fRMS; 
	
	bool fIsMC;
	
        ClassDef(SBSGEMPlane,0)

};




