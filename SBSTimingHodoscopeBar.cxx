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

SBSTimingHodoscopeBar::~SBSTimingHodoscopeBar() = default;

//_____________________________________________________________
void SBSTimingHodoscopeBar::Clear( Option_t* ) {
  fBarNum = 0;
  fBarOff = 0;
}
//____________________________________________________________________





///////////////////////////////////////////////////////////////////////////////
