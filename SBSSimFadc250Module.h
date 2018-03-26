#ifndef SBSSimFadc250Module_
#define SBSSimFadc250Module_

/////////////////////////////////////////////////////////////////////
//
//   SBSSimFadc250Module
//   G4SBS Simplified Simulated JLab FADC 250 Module
//   Right now it only has the basic features needed to get SBSHCal
//   Nothing else in the Fadc250Module has been copied or implemented
//
/////////////////////////////////////////////////////////////////////

#include "Fadc250Module.h"
#include "stdint.h"
#include <vector>

namespace Decoder {

  class SBSSimFadc250Module : public Fadc250Module {   // Inheritance

  public:

    SBSSimFadc250Module();                         // Default constructor
    SBSSimFadc250Module(Int_t crate, Int_t slot);  // Constructor
    virtual ~SBSSimFadc250Module();                // Virtual constructor

    using Module::GetData;
    using Module::LoadSlot;

    void SetPulseSamplesVector(Int_t chan, std::vector<uint32_t> pulses);

    virtual void Clear(const Option_t *opt="");
    virtual void Init();
    virtual void CheckDecoderStatus() const;
    virtual std::vector<uint32_t> GetPulseSamplesVector(Int_t chan) const;
    virtual Int_t GetPulseIntegralData(Int_t chan, Int_t ievent) const;
    virtual Int_t GetFadcMode() const { return 8; };
    virtual Int_t GetMode() const { return GetFadcMode(); };
    virtual Int_t GetNumFadcEvents(Int_t chan) const;
    virtual Int_t GetNumFadcSamples(Int_t chan, Int_t ievent) const;
    virtual Bool_t HasCapability(Decoder::EModuleType type);
    virtual void ClearDataVectors();
    virtual Int_t LoadSlot(THaSlotData *sldat, const UInt_t* evbuffer, const UInt_t *pstop);
    virtual Int_t LoadSlot(THaSlotData *sldat, const UInt_t* evbuffer, Int_t pos, Int_t len);
//    virtual Int_t DecodeOneWord(UInt_t pdat);
// We dont use the Decode() but if you dont define it the class is abstract and wont be instantiated
    Int_t Decode(const UInt_t *pdat) { return 0; }; // use DecodeOneWord instead
//    virtual Bool_t IsMultiFunction();
//    Int_t LoadNextEvBuffer(THaSlotData *sldat);
//    Int_t GetData(EModuleType mtype, Int_t chan, Int_t ievent) const;
//    Int_t GetNumEvents(EModuleType mtype, Int_t ichan) const;
//    Int_t GetNumEvents() const { return GetNumEvents(0); } ;
//    Int_t GetNumEvents(Int_t ichan) const { return GetNumFadcEvents(ichan); } ;
//    Int_t GetNumSamples(Int_t ichan) const { return GetNumFadcSamples(ichan, 0);};


    struct fadc_data_struct {
      std::vector<uint32_t> samples;
      std::vector<uint32_t> integrals;
      void clear() {
        samples.clear();
        integrals.clear();
      }
    };  // fadc_data_struct

  private:
    static const size_t NADCCHAN = 16;
    static TypeIter_t fgThisType;
    std::vector<fadc_data_struct> fadc_data;
    ClassDef(SBSSimFadc250Module,0)  //  JLab FADC 250 Module

  };  // SBSSimFadc250Module class

}  // Decoder namespace

#endif
