#ifndef ROOT_SBSAdcHit
#define ROOT_SBSAdcHit

/////////////////////////////////////////////////////////////////////////////
//                                                                         //
// SBSAdcHit                                                               //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

#include "TObject.h"
#include "TRef.h"
#include "SBSScintPMT.h"

class SBSAdcHit : public TObject {
  
 public:
  explicit SBSAdcHit(SBSScintPMT* pmt=nullptr, Int_t rawampl=0);
  virtual ~SBSAdcHit() = default;

  virtual void Clear(Option_t *s="");

  SBSScintPMT* GetPMT() const { return (SBSScintPMT*)fPMT.GetObject(); }
  Int_t GetRawAmpl() const {return fRawAmpl;}
  Int_t GetAmplPedCor() const {return fAmplPedCor;}
  Double_t GetAmpl() const {return fAmpl;}
  Int_t GetBarNum() const {return fBarNum;}
  Int_t GetSide() const {return fSide;}

  void CorrectHit();
  void SetRawAmpl(Int_t rawampl) { fRawAmpl = rawampl;}
  void SetBarNum(Int_t barnum) {fBarNum=barnum;}
  void SetSide(Int_t side) {fSide=side;}

  Bool_t IsSortable() const { return kTRUE; }

  Int_t  Compare(const TObject* obj) const;

 protected:

  TRef fPMT;  // reference to "real" PMT, kept in detector arrays
  Int_t fRawAmpl;
  Int_t fAmplPedCor;
  Double_t fAmpl;
  Int_t fBarNum;
  Int_t fSide;

  ClassDef(SBSAdcHit,1)  // ADC and calibrated amplitude hit
};
/////////////////////////////////////////////////////////////////
#endif
