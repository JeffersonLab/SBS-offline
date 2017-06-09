#ifndef ROOT_SBSScintPMT
#define ROOT_SBSScintPMT


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// SBSScintPMT                                                               //
//                                                                           //
// Class to represent a PMT on the neutron bars                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

//class SBSScintPMT;
class SBSScintBar;

#include "TObject.h"
#include "TRef.h"

enum ESide { kLeft = 0, kRight = 1 };

class SBSScintPMT : public TObject {

 public: 

  SBSScintPMT( 
	  Double_t gain=1.0, 
	  Int_t ped=0, 
	  Double_t res=1.0, 
	  Double_t off=0.0,
	  Double_t walk=0.0, 
	  SBSScintBar* bar=0, 
	  Int_t barnum=0, 
	  Int_t side=0,
      Int_t lowlim=0, 
	  Int_t uplim=65536, 
	  Double_t wraparound=0,
	  Double_t wexp=-0.5
	  );
    
  virtual ~SBSScintPMT(void);

  SBSScintBar* GetScintBar() const { return (SBSScintBar*)fScBar.GetObject(); }

  void SetGain(Double_t gain) {fGain=gain;}
  void SetPed(Int_t ped) {fPed=ped;}
  void SetTDCRes(Double_t res) {fTDCRes=res;}
  void SetTOffset(Double_t off) {fTOffset=off;}
  void SetTimeWalk(Double_t walk) {fTimeWalkPar=walk;}
  void SetTimeWExp(Double_t wexp) {fTimeWalkExp=wexp;}
  void SetScintBar(SBSScintBar* bar);
  void SetBarNum(Int_t barnum) {fBarNum=barnum;}
  void SetSide(Int_t side) {fSide=side;}
  void SetRawLowLim(Int_t lowlim) {fRawLowLim=lowlim;}
  void SetRawUpLim(Int_t uplim) {fRawUpLim=uplim;}
  void SetRawWrapAround(Double_t wraparound) {fRawWrapAround=wraparound;}

  Double_t GetGain()     const     {return fGain;}
  Int_t    GetPed()      const     {return fPed;}
  Double_t GetTDCRes()   const     {return fTDCRes;}
  Double_t GetTOffset()  const     {return fTOffset;}
  Double_t GetTimeWalk() const     {return fTimeWalkPar;}
  Double_t GetTimeWExp() const     {return fTimeWalkExp;}
  Int_t    GetBarNum()   const     {return fBarNum;}
  Int_t    GetSide()     const     {return fSide;}

  Int_t GetRawLowLim()   const     {return fRawLowLim;}
  Int_t GetRawUpLim()    const     {return fRawUpLim;}
  Double_t GetRawWrapAround() const {return fRawWrapAround;}

 protected: 
  Double_t fGain;
  Int_t fPed;
  Double_t fTDCRes;
  Double_t fTOffset;
  Double_t fTimeWalkPar;
  Double_t fTimeWalkExp;
  TRef   fScBar;       // to the original bar
  Int_t fBarNum;
  Int_t fSide;
  Int_t fRawLowLim;
  Int_t fRawUpLim;
  Double_t fRawWrapAround;

 public:
  ClassDef(SBSScintPMT,2)  // Class to represent a PMT on the neutron bars

};
////////////////////////////////////////////////////////////////////////////////

#endif
