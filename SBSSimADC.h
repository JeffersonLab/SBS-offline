#ifndef SBSSimADC_
#define SBSSimADC_

/////////////////////////////////////////////////////////////////////
//
//   SBSSimADC
//   G4SBS Simplified Simulated ADC Module
//
//   Right now, all it supports is an ADC sum and individual ADC samples
//
/////////////////////////////////////////////////////////////////////

#include "PipeliningModule.h"
#include "stdint.h"
#include <vector>
#include "SBSSimDataDecoder.h"

namespace Decoder {

  class SBSSimADC : public PipeliningModule {   // Inheritance

  public:

    SBSSimADC();                         // Default constructor
    SBSSimADC(Int_t crate, Int_t slot);  // Constructor
    virtual ~SBSSimADC();                // Virtual constructor

    // Use parent class functions
    using PipeliningModule::Init;

    virtual void Clear(const Option_t *opt="");
    virtual void Init();
    virtual void CheckDecoderStatus() const;
    virtual void ClearDataVectors();
    virtual UInt_t LoadSlot(THaSlotData *sldat, const UInt_t* evbuffer, const UInt_t *pstop);
    virtual UInt_t LoadSlot( THaSlotData *sldat, const UInt_t* evbuffer, UInt_t pos, UInt_t len);
    
    // We don't need these functions for simulated data, but must be defined
    // so that this won't be an abstract class
    UInt_t LoadNextEvBuffer(THaSlotData*) {return 0;};// needs return something for compilation
    UInt_t LoadThisBlock(THaSlotData*, const std::vector<UInt_t >&) {return 0;};// needs return something for compilation
    Int_t Decode(const UInt_t *) { return 0; }; // use DecodeOneWord instead

    /*
    struct fadc_data_struct {
      std::vector<uint32_t> samples;
      std::vector<uint32_t> integrals;
      void clear() {
        samples.clear();
        integrals.clear();
      }
    };  // fadc_data_struct
    */

  private:
    static const size_t NADCCHAN = 2048; // Max ADC channels
    static TypeIter_t fgThisType;
    static TypeIter_t fgType1;
    static TypeIter_t fgType2;
    static TypeIter_t fgType3;
    //static std::vector<TypeIter_t> fgTypeVec;//Let's try something...
    std::vector<SimEncoder::sadc_data> sadc_data;

    ClassDef(SBSSimADC,0)  //  Generic SimADC module

  };  // SBSSimADC class

}  // Decoder namespace

#endif
