 /* #ifndef ROOT_SBSTimingHodoscopePMT  */
 /* #define ROOT_SBSTimingHodoscopePMT  */
#ifndef SBSTimingHodoscopePMT_h
#define SBSTimingHodoscopePMT_h


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// SBSTimingHodoscopePMT                                                     //
//                                                                           //
// Class to represent a PMT on the timing hodoscope                          //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

/* #include "TObject.h" */
#include "SBSElement.h"

class SBSTimingHodoscopePMT{

public:
  SBSTimingHodoscopePMT( SBSElement* element, Float_t walkpar0, Float_t walkpar1,Int_t barnum, Int_t side, Int_t id );
  virtual ~SBSTimingHodoscopePMT();

  // getter functions
  Float_t                GetTimeWalkPar0() const { return fTimeWalkParameter0;}
  Float_t                GetTimeWalkPar1() const { return fTimeWalkParameter1;}
  Int_t                  GetBarNum()       const { return fBarNum;}
  Int_t                  GetSide()         const { return fSide;}
  Int_t                  GetId()           const { return fId;}
  SBSElement*            GetPMTElement()         { return fPMTElement;}

  // setter functions
  void SetTimeWalkPar0(Float_t walkpar0) {fTimeWalkParameter0=walkpar0;}
  void SetTimeWalkPar1(Float_t walkpar1) {fTimeWalkParameter1=walkpar1;}
  void SetBarNum(Int_t barnum) {fBarNum=barnum;}
  void SetSide(Int_t side) {fSide=side;}
  void SetId(Int_t id) {fId=id;}
  void ClearEvent();

 /* private:  */
 protected:
  SBSElement* fPMTElement;
  Float_t fTimeWalkParameter0;
  Float_t fTimeWalkParameter1;
  Int_t fBarNum;
  Int_t fSide;
  Int_t fId; // this is the global element id in fElements
/* const Int_t nTDC */

 /* public: */
  ClassDef(SBSTimingHodoscopePMT,5)  // Class to represent a PMT on the timing hodoscope

};
////////////////////////////////////////////////////////////////////////////////

#endif
