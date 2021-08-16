///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// SBSTimingHodoscopeBar                                                     //
//                                                                           //
// rmontgom@jlab.org drafting jul 2021                                       //
///////////////////////////////////////////////////////////////////////////////

#include "SBSTimingHodoscopeBar.h"

ClassImp(SBSTimingHodoscopeBar)

// constructor 
SBSTimingHodoscopeBar::SBSTimingHodoscopeBar( Int_t barnum,
					      SBSTimingHodoscopePMT* leftpmt,
					      SBSTimingHodoscopePMT* rightpmt,
					      Int_t baroff) :
fBarNum(barnum),fLPMT(leftpmt),fRPMT(rightpmt), fBarOff(baroff)
{
}

//____________________________________________________________________

SBSTimingHodoscopeBar::~SBSTimingHodoscopeBar() { 
  ClearEvent();
}
//_____________________________________________________________
void SBSTimingHodoscopeBar::ClearEvent() {
  fBarNum = 0;
  delete fLPMT;
  delete fRPMT;
  fBarOff = 0;
}
//____________________________________________________________________





///////////////////////////////////////////////////////////////////////////////
