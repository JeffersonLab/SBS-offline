///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// SBSScintPMT                                                               //
//                                                                           //
// Class to represent a PMT on the neutron bars                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "SBSScintPMT.h"
#include "SBSScintBar.h"

//_____________________________________________________________________________
SBSScintPMT::SBSScintPMT( Double_t gain, Int_t ped, Double_t res, Double_t off,
			  Double_t walk, SBSScintBar* bar, Int_t barnum,
			  Int_t side,
			  Int_t lowlim, Int_t uplim, Double_t wraparound,
			  Double_t wexp ) :
  fGain(gain), fPed(ped), fTDCRes(res), fTOffset(off), fTimeWalkPar(walk),
  fTimeWalkExp(wexp), fScBar(bar), fBarNum(barnum), fSide(side),
  fRawLowLim(lowlim), fRawUpLim(uplim), fRawWrapAround(wraparound) { ; }

//_____________________________________________________________________________
SBSScintPMT::~SBSScintPMT(void) { ; }

//_____________________________________________________________________________
void SBSScintPMT::SetScintBar(SBSScintBar* bar)
{
  fScBar=bar;
}

//_____________________________________________________________________________
ClassImp(SBSScintPMT)

///////////////////////////////////////////////////////////////////////////////
