#ifndef SBSDecodeF1TDCModule_
#define SBSDecodeF1TDCModule_

/////////////////////////////////////////////////////////////////////
//
//   SBSDecodeF1TDCModule
//   JLab F1 TDC Module
//
/////////////////////////////////////////////////////////////////////

#include "VmeModule.h"

namespace Decoder {

class SBSDecodeF1TDCModule : public VmeModule {

public:

   SBSDecodeF1TDCModule() {};
   SBSDecodeF1TDCModule(Int_t crate, Int_t slot);
   virtual ~SBSDecodeF1TDCModule();

   using Module::GetData;
   using Module::LoadSlot;

   enum EResolution { ILO = 0, IHI = 1 };

   virtual void Init();
   virtual void Clear(const Option_t *opt="");
   virtual Bool_t IsSlot(UInt_t rdata);
   virtual Int_t GetData(Int_t chan, Int_t hit) const;

   void SetResolution(Int_t which=0) {
     fResol = IHI;
     if (which==0) fResol=ILO;
     return;
   }
   EResolution GetResolution() const { return fResol; };
   Bool_t IsHiResolution() const { return (fResol==IHI); };

   Int_t GetNumHits(Int_t chan);// const { return fNumHits; };
   Int_t Decode(const UInt_t *p) { return 0; };

   // For multiple slots
   Int_t GetNumSlots() const { return nF1; };

   // Loads slot data for bank structures
   virtual Int_t LoadSlot(THaSlotData *sldat, const UInt_t *evbuffer, Int_t pos, Int_t len);
// Loads sldat and increments ptr to evbuffer
   Int_t LoadSlot(THaSlotData *sldat,  const UInt_t* evbuffer, const UInt_t *pstop );

private:


   //Int_t *fNumHits;
   //Int_t *fTdcData;  // Raw data (either samples or pulse integrals)
   std::vector<Int_t> fNumHits;
   std::vector<Int_t> fTdcData;  // Raw data (either samples or pulse integrals)
   EResolution fResol;
   Bool_t IsInit;
   Int_t slotmask, chanmask, datamask;

   // For multiple slots
   Int_t nF1;
   Int_t *F1slots;

   static TypeIter_t fgThisType;
   ClassDef(SBSDecodeF1TDCModule,0)  //  JLab F1 TDC Module, test version

};

}

#endif
