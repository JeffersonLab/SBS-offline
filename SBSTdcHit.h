#ifndef ROOT_SBSTdcHit
#define ROOT_SBSTdcHit

/////////////////////////////////////////////////////////////////////////////
//                                                                         //
// SBSTdcHit                                                             //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

#include "TObject.h"
#include "TRef.h"
#include "SBSScintPMT.h"
#include <cstdio>

class SBSTdcHit : public TObject {
  
 public:
  SBSTdcHit(const SBSScintPMT* pmt=NULL, Int_t rawtime=0, Double_t ext_offset=0);
  virtual ~SBSTdcHit() {}

  SBSScintPMT* GetPMT() const { return (SBSScintPMT*)fPMT.GetObject(); }

  Int_t    GetRawTime() const {return fRawTime;}
  Double_t GetTime() const {return fTime;}
  Int_t    GetBarNum() const {return fBarNum;}
  Int_t    GetSide() const {return fSide;}

  void UpdateTime(Double_t ext_offset=0.);
  void SetRawTime(Int_t rawtime) { fRawTime = rawtime;}

  Bool_t IsSortable() const { return kTRUE; }

  Int_t  Compare(const TObject* obj) const;

  void   Clear(Option_t *s="");
  
 protected:

  TRef fPMT;  // reference to "real" PMT, kept in detector arrays
  Int_t fRawTime;
  Double_t fTime;
  Int_t fBarNum;
  Int_t fSide;

 public:
  ClassDef(SBSTdcHit,1)  // TDC and real-time for a Hit. Per PMT
};
/////////////////////////////////////////////////////////////////
#endif
