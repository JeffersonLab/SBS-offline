#ifndef MPDModuleVMEv4_
#define MPDModuleVMEv4_

/////////////////////////////////////////////////////////////////////
//
//   MPDModule
//   This is the MPD INFN (GEM readout for SBS) module decoder
//   Note, it must inherit from VmeModule.
//   Feel free to copy this and make appropriate changes.
//
//   other steps
//   1. Register (add "DoRegister" call and include the header)
//   note: the number (4444) that is registered must appear in db_cratemap.dat
//   2. Add to namespace Decoder.h
//   3. Add to Makefile
//   4. Add to haDecode_LinkDef.h
//   5. Add line(s) for [crate,slot] in db_cratemap.dat
//
//   if the pre-compiler flag "LIKEV792" is defined, the decoding is
//   sort of like a V792 ... as an example.
//
/////////////////////////////////////////////////////////////////////

//#define LIKEV792x 0

// (jc2) What was this used for? Seems to conflict with a value defined
// in Caen775. So I've commented out.
//#define NTDCCHAN   32
#define MAXHIT    2048

#include "VmeModule.h"

namespace Decoder {

  class MPDModuleVMEv4 : public VmeModule {

  public:

    MPDModuleVMEv4() = default;
    MPDModuleVMEv4(Int_t crate, Int_t slot);
    virtual ~MPDModuleVMEv4();

    using VmeModule::GetData;
    using VmeModule::LoadSlot;
    using VmeModule::Init;

    virtual UInt_t GetData( UInt_t adc, UInt_t sample, UInt_t chan) const;
    virtual void Init();
    virtual void Clear(const Option_t *opt="");
    virtual Int_t Decode(const UInt_t *p); // { return 0; };
    
    /*
    void Config(Int_t mode, Int_t sampleperiod, Int_t nsample, Int_t nadc=16, Int_t nch=128) {
      fAcqMode = mode;
      fSamplePeriod = sampleperiod;
      fNumSample = nsample;
      fNumADC=nadc;
      fNumChan=nch;
      fData.clear();
      fFrameHeader.clear();
      fFrameTrailer.clear();
      for (Int_t i=0;i<fNumADC*fNumSample*fNumChan;i++) {
	fData.push_back(0);
      }
      for (Int_t i=0;i<fNumADC*fNumSample;i++) { 
	fFrameHeader.push_back(0);
	fFrameTrailer.push_back(0);
      }
      IsInit = kTRUE;
      //     CheckSetMode();
    }
    */

    
//#ifdef LIKEV792x
    // Loads slot data.  if you don't define this, the base class's method is used
    virtual UInt_t LoadSlot( THaSlotData *sldat, const UInt_t *evbuffer, UInt_t pos, UInt_t len);
//#endif

    //void CommonModeSubtraction();
    
  private:

    //bool fOnlineZeroSuppression; //if true, assumes that raw-data are already zero-suppressed and baseline-subtracted
    
    // configuration parameters
    Int_t fAcqMode; // normal, zero suppression, histogram, synch ...
    Int_t fSamplePeriod; // 25 ns, 75 ns ...
    UInt_t fNumSample; // number of sample / event
    
    Int_t fNumADC; // number of ADC fifos (number of front end cards served by the MPD)
    
    // current indices
    Int_t fIdxA; // ADC
    Int_t fIdxS; // Sample
    Int_t fIdxC; // Channel

    Int_t fIdxMPD; // MPD ID
  
    Int_t fCountS; // Sample Counter from electronics
    Int_t fCountW; // Word 

    Int_t fNumHits;

    //New members/methods to parse MPD4 VME format used by UVA GEM cosmic test stand, Jan. 2021:
    UInt_t fBlockHeader; //Default = 0x0
    UInt_t fAPVHeader;   //Default = 0x4
    
  
    std::vector<Int_t> fFrameHeader;  // Frame Header
    std::vector<Int_t> fFrameTrailer;  // Frame Trailer

    static TypeIter_t fgThisType;

    // linearization of the indeces 
    // inline Int_t as2i(Int_t adc, Int_t sample) const {
    //   return adc*fNumSample + sample;
    // };

    // inline UInt_t asc2i(UInt_t adc, UInt_t sample, UInt_t chan) const {
    //   return adc*fNumSample*fNumChan + sample*fNumChan + chan;
    // };
    
    ClassDef(MPDModuleVMEv4,0)  //  INFN MPD Module 

  };

}

#endif
