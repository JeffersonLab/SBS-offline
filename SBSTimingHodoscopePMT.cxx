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
SBSTimingHodoscopePMT::SBSTimingHodoscopePMT( SBSElement* element, Double_t walkpar0, Double_t walkpar1, Int_t barnum, Int_t side, Int_t id, bool tdcflag, bool adcflag) : fPMTElement(element), fTimeWalkPar0(walkpar0), fTimeWalkPar1(walkpar1), fBarNum(barnum), fSide(side), fId(id), fTdcFlag(tdcflag), fAdcFlag(adcflag){
}

//_____________________________________________________________________________
SBSTimingHodoscopePMT::~SBSTimingHodoscopePMT() {
  ClearEvent();
}
//_____________________________________________________________
void SBSTimingHodoscopePMT::ClearEvent() {
  // fPMTElement = 0;
  delete fPMTElement;
  fTimeWalkPar0 = 0.0;
  fTimeWalkPar1 = 0.0;
  fBarNum = 0;
  fSide = 0;
  fId = 0;
  fTdcFlag = false;
  fAdcFlag = false;
}
///////////////////////////////////////////////////////////////////////////////
