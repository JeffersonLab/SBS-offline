///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// SBSTimingHodoscopePMT                                                     //
//                                                                           //
// Class to represent a PMT on the timing hodoscope                          //
// rmontgom@jlab.org drafting jul 2021                                       //
///////////////////////////////////////////////////////////////////////////////

#include "SBSTimingHodoscopePMT.h"

ClassImp(SBSTimingHodoscopePMT);

//_____________________________________________________________________________
SBSTimingHodoscopePMT::SBSTimingHodoscopePMT( SBSElement* element, Float_t walkpar0, Float_t walkpar1, Int_t barnum, Int_t side, Int_t id) : fPMTElement(element), fTimeWalkParameter0(walkpar0), fTimeWalkParameter1(walkpar1), fBarNum(barnum), fSide(side), fId(id){
}

//_____________________________________________________________________________
SBSTimingHodoscopePMT::~SBSTimingHodoscopePMT() {
  ClearEvent();
}
//_____________________________________________________________
void SBSTimingHodoscopePMT::ClearEvent() {
  // fPMTElement = 0;
  delete fPMTElement;
  fTimeWalkParameter0 = 0.0;
  fTimeWalkParameter1 = 0.0;
  fBarNum = 0;
  fSide = 0;
  fId = 0;
}
///////////////////////////////////////////////////////////////////////////////
