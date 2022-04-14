///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// SBSTimingHodoscopePMT                                                     //
//                                                                           //
// Class to represent a PMT on the timing hodoscope                          //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "SBSTimingHodoscopePMT.h"

ClassImp(SBSTimingHodoscopePMT);

//_____________________________________________________________________________
SBSTimingHodoscopePMT::SBSTimingHodoscopePMT(
  SBSElement* element, Double_t walkpar0, Double_t walkpar1, Int_t barnum,
  Int_t side, Int_t id )
  : fPMTElement(element),
    fTimeWalkParameter0(walkpar0),
    fTimeWalkParameter1(walkpar1),
    fBarNum(barnum),
    fSide(side),
    fId(id) {
}

//_____________________________________________________________________________
SBSTimingHodoscopePMT::~SBSTimingHodoscopePMT() = default;

//_____________________________________________________________
void SBSTimingHodoscopePMT::Clear( Option_t* )
{
  // fPMTElement = 0;
  fTimeWalkParameter0 = 0.0;
  fTimeWalkParameter1 = 0.0;
  fBarNum = 0;
  fSide = 0;
  fId = 0;
}
///////////////////////////////////////////////////////////////////////////////
