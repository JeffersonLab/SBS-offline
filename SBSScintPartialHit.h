#ifndef ROOT_SBSScintPartialHit
#define ROOT_SBSScintPartialHit

/////////////////////////////////////////////////////////////////////////////
//                                                                         //
// SBSScintPartialHit                                                             //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

#include "TObject.h"
#include "TRef.h"
#include "SBSScintBar.h"
#include <cstdio>

class SBSScintPartialHit : public TObject {
  
 public:
  
  SBSScintPartialHit( SBSScintBar* bar = NULL, Int_t barnum=-7, Int_t CaseNum=0, Double_t lt=0.0, Double_t lt_raw=0.0,
		 Double_t rt=0.0, Double_t rt_raw=0.0, Double_t la=0.0,Double_t la_raw=0.0,
		 Double_t ra=0.0 , Double_t ra_raw=0.0 );
    
  virtual ~SBSScintPartialHit() {;}

  SBSScintBar* GetScintBar() const { return (SBSScintBar*)fScBar.GetObject(); } 
  Int_t GetBarNum() {return fBarNum;}
  Int_t GetCaseNum() {return fCaseNum;}
  Double_t GetLt() {return fLt;}
  Double_t GetRt() {return fRt;}
  Double_t GetLa() {return fLa;}
  Double_t GetRa() {return fRa;}
  Double_t GetLt_raw() {return fLt_raw;}
  Double_t GetRt_raw() {return fRt_raw;}
  Double_t GetLa_raw() {return fLa_raw;}
  Double_t GetRa_raw() {return fRa_raw;}
 
  void SetScintBar(SBSScintBar* bar) {fScBar=bar;}
  void SetBarNum(Int_t barnum) {fBarNum=barnum;}
  void SetCaseNum(Int_t Val) { fCaseNum = Val;}
  void SetLt(Double_t Val) { fLt = Val;}
  void SetRt(Double_t Val) { fRt = Val;}
  void SetLa(Double_t Val) { fLa = Val;}
  void SetRa(Double_t Val) { fRa = Val;}
  void SetLt_raw(Double_t Val) { fLt_raw = Val;}
  void SetRt_raw(Double_t Val) { fRt_raw = Val;}
  void SetLa_raw(Double_t Val) { fLa_raw = Val;}
  void SetRa_raw(Double_t Val) { fRa_raw = Val;}

 protected:

  TRef fScBar;             // reference to the bar
  Int_t fBarNum;
  Int_t fCaseNum;
  Double_t fLt;
  Double_t fLt_raw;
  Double_t fRt;
  Double_t fRt_raw;
  Double_t fLa;
  Double_t fLa_raw;
  Double_t fRa;
  Double_t fRa_raw;

 public:
  ClassDef(SBSScintPartialHit,1) // Partial(not complete L/R A/T) scintillator hit
};
/////////////////////////////////////////////////////////////////
#endif
