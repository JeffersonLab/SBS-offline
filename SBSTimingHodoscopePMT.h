/* #ifndef ROOT_SBSTimingHodoscopePMT */
/* #define ROOT_SBSTimingHodoscopePMT */
#ifndef SBSTimingHodoscopePMT_h
#define SBSTimingHodoscopePMT_h


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// SBSTimingHodoscopePMT                                                     //
//                                                                           //
// Class to represent a PMT on the timing hodoscope                          //
// rmontgom@jlab.org drafting jul 2021                                       //
///////////////////////////////////////////////////////////////////////////////

#include "TObject.h"
#include "SBSElement.h"

class SBSTimingHodoscopePMT{

public:
  SBSTimingHodoscopePMT( SBSElement* element, Double_t walkpar0, Double_t walkpar1,Int_t barnum, Int_t side, Int_t id );
  virtual ~SBSTimingHodoscopePMT();

  // getter functions
  Double_t               GetTimeWalkPar0() const { return fTimeWalkPar0;}
  Double_t               GetTimeWalkPar1() const { return fTimeWalkPar1;}
  Int_t                  GetBarNum()       const { return fBarNum;}
  Int_t                  GetSide()         const { return fSide;}
  Int_t                  GetId()           const { return fId;}
  SBSElement*            GetPMTElement()         { return fPMTElement;}

  // setter functions
  void SetTimeWalkPar0(Double_t walkpar0) {fTimeWalkPar0=walkpar0;}
  void SetTimeWalkPar1(Double_t walkpar1) {fTimeWalkPar1=walkpar1;}
  void SetBarNum(Int_t barnum) {fBarNum=barnum;}
  void SetSide(Int_t side) {fSide=side;}
  void SetId(Int_t id) {fId=id;}
  void ClearEvent();

 /* private:  */
 protected:
  SBSElement* fPMTElement;
  Double_t fTimeWalkPar0;
  Double_t fTimeWalkPar1;
  Int_t fBarNum;
  Int_t fSide;
  Int_t fId; // this is the global element id in fElements
/* const Int_t nTDC */

 /* public: */
  ClassDef(SBSTimingHodoscopePMT,5)  // Class to represent a PMT on the timing hodoscope

};
////////////////////////////////////////////////////////////////////////////////

#endif
