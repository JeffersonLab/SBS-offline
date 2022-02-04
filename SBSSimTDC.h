#ifndef SBSSimTDC_
#define SBSSimTDC_

/////////////////////////////////////////////////////////////////////
//
//   SBSSimTDC
//   G4SBS Simplified Simulated ADC Module
//
//   Right now, all it supports is an ADC sum and individual ADC samples
//
/////////////////////////////////////////////////////////////////////

#include "PipeliningModule.h"
#include "stdint.h"
#include <vector>

namespace Decoder {

  class SBSSimTDC : public PipeliningModule {   // Inheritance

  public:

    SBSSimTDC();                         // Default constructor
    SBSSimTDC(Int_t crate, Int_t slot);  // Constructor
    virtual ~SBSSimTDC();                // Virtual constructor

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



    struct tdc_data_struct {
      std::vector<uint32_t> lead_time;
      std::vector<uint32_t> trail_time;
    };  // tdc_data_struct

  private:
    static const size_t NTDCCHAN = 129; // Max TDC channels, *including ref chan*
    static TypeIter_t fgThisType;
    static TypeIter_t fgType0;
    static TypeIter_t fgType1;
    static TypeIter_t fgType2;
    static TypeIter_t fgType3;
    std::vector<tdc_data_struct> tdc_data;

    ClassDef(SBSSimTDC,0)  //  Generic SimTDC module

  };  // SBSSimTDC class

}  // Decoder namespace

#endif
