///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// SBSAdcHit                                                                 //
//                                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "SBSAdcHit.h"
#include "SBSScintPMT.h"

using namespace std;


//_____________________________________________________________________________
SBSAdcHit::SBSAdcHit(SBSScintPMT* pmt, Int_t rawampl):fPMT(pmt),
						      fRawAmpl(rawampl) 
{
  // Construct and correct the Adc hit for calibrations.
  if (pmt) {
    fBarNum = pmt->GetBarNum();
    fSide = pmt->GetSide();
  } else {
    fBarNum = 0;
    fSide = 0;
  }
  CorrectHit();
}

//_____________________________________________________________________________
void SBSAdcHit::CorrectHit() {
  // Apply the associated PMT's pedestal and gain correction factors.

  SBSScintPMT* pmt = GetPMT();
  if (pmt) {
    fAmplPedCor = fRawAmpl - pmt->GetPed();
    fAmpl = fAmplPedCor * pmt->GetGain();
  } else {
    fAmplPedCor = 0;
    fAmpl = 0;
  }
}

//_____________________________________________________________________________
Int_t SBSAdcHit::Compare(const TObject *obj) const {
  // sort SBSAdcHits according to side and bar number.
  //   obj should be a SBSAdcHit, and we will assume as much
  //
  // Returns -1 if this < obj, 0 if this==obj, and 1 if this>obj
  const SBSAdcHit *h = static_cast<const SBSAdcHit*>(obj);

  // Side (Left-right) comparison
  if (fSide < h->fSide) return -1;
  if (fSide > h->fSide) return  1;

  // Bar-number comparison
  if (fBarNum < h->fBarNum) return -1;
  if (fBarNum > h->fBarNum) return  1;

  // finally, Amplitude comparison, highest Amplitude "wins" (first)
  if (fAmpl < h->fAmpl) return  1;
  if (fAmpl > h->fAmpl) return -1;
  
  return 0;
}

//_____________________________________________________________________________
void SBSAdcHit::Clear(Option_t *) {
  // clear the data inside SBSAdcHit, so it does not have to be deleted each
  // time

  fPMT = nullptr;
  fRawAmpl = fAmplPedCor = fSide = 0;
  fBarNum = -1;
  fAmpl = 0;
}

///////////////////////////////////////////////////////////////////////////////

ClassImp(SBSAdcHit)
