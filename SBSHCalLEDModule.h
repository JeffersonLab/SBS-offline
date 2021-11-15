#ifndef HCalLED_
#define HCalLED_

/////////////////////////////////////////////////////////////////////
//
//   HCalLED
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

  class HCalLED : public VmeModule {

  public:

    HCalLED() {};
    HCalLED(Int_t crate, Int_t slot);
    virtual ~HCalLED();

    using VmeModule::GetData;
    using VmeModule::LoadSlot;
    using VmeModule::Init;

    //virtual Int_t GetData(Int_t adc, Int_t sample, Int_t chan) const;
    virtual void Init();
    //virtual void Clear(const Option_t *opt);
    virtual Int_t Decode(const UInt_t *p); // { return 0; };
    
    virtual UInt_t LoadSlot(THaSlotData *sldat,  const UInt_t *evbuffer, UInt_t pos, UInt_t len);
//#endif

  private:
    static TypeIter_t fgThisType;
    
    ClassDef(HCalLED,0)

  };

}

#endif
