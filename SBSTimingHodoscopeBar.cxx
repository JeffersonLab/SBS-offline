///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// SBSTimingHodoscopeBar                                                     //
//                                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "SBSTimingHodoscopeBar.h"

ClassImp(SBSTimingHodoscopeBar)

// constructor 
SBSTimingHodoscopeBar::SBSTimingHodoscopeBar( Int_t barnum,
					      SBSTimingHodoscopePMT* leftpmt,
					      SBSTimingHodoscopePMT* rightpmt,
					      Int_t baroff) :
fBarNum(barnum),fBarOff(baroff),fLPMT(leftpmt), fRPMT(rightpmt)
{
}

//____________________________________________________________________

SBSTimingHodoscopeBar::~SBSTimingHodoscopeBar() { 
  ClearEvent();
}
//_____________________________________________________________
void SBSTimingHodoscopeBar::ClearEvent() {
  fBarNum = 0;
  fBarOff = 0;
}
//____________________________________________________________________





///////////////////////////////////////////////////////////////////////////////
