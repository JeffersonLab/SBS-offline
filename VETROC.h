#ifndef Podd_VETROCModule_h_
#define Podd_VETROCModule_h_

/////////////////////////////////////////////////////////////////////
//
//   VETROCModule
//
/////////////////////////////////////////////////////////////////////


#include "VmeModule.h"
#include <cstring>  // for memset
#include <vector>

namespace Decoder {

  class VETROCModule : public VmeModule {

  public:

    VETROCModule() : slot_data(nullptr) {}
    VETROCModule(Int_t crate, Int_t slot);
    virtual ~VETROCModule() = default;

    using VmeModule::GetData;
    using VmeModule::Init;
    using VmeModule::GetOpt;

    virtual void  Init();
    virtual void  Clear(Option_t *opt="");
    virtual Int_t Decode(const UInt_t *p);
    virtual UInt_t GetData( UInt_t chan, UInt_t hit) const;
    virtual UInt_t GetOpt( UInt_t chan, UInt_t hit) const;

    // Loads slot data.  if you don't define this, the base class's method is used
    virtual UInt_t LoadSlot( THaSlotData *sldat, const UInt_t *evbuffer, const UInt_t *pstop );
    // Loads slot data for bank structures
    virtual UInt_t LoadSlot( THaSlotData *sldat, const UInt_t *evbuffer, UInt_t pos, UInt_t len);

  private:

    std::vector<UInt_t> fNumHits;
    std::vector<UInt_t> fTdcData;  // Raw data
    std::vector<UInt_t> fTdcOpt;  // Edge flag =0 Leading edge, = 1 Trailing edge

    THaSlotData *slot_data;  // Need to fix if multi-threading becomes available
   
    class tdcData {
    public:
      tdcData() :
      glb_hdr_evno(0), glb_hdr_slno(0), evh_trig_num(0), chan(0), raw(0), opt(0),trig_time_l(0),trig_time_h(0), trig_time(0),
        status(0) {}
      void clear() { memset(this, 0, sizeof(tdcData)); }
      UInt_t glb_hdr_evno, glb_hdr_slno, evh_trig_num;
      UInt_t chan, raw , opt ,trig_time_l,trig_time_h, trig_time;
      Int_t status;
    } tdc_data;

    static TypeIter_t fgThisType;
    ClassDef(VETROCModule,0)  //  VETROC of a module; make your replacements

      };

}

#endif
