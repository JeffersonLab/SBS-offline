///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// SBSScintPartialHit                                                             //
//                                                                           //
// Class to represent a Partial hit on the neutron bars                      //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "SBSScintPartialHit.h"

ClassImp(SBSScintPartialHit)

SBSScintPartialHit::SBSScintPartialHit( SBSScintBar* bar, Int_t barnum, Int_t CaseNum,
			      Double_t lt, Double_t lt_raw,
			      Double_t rt, Double_t rt_raw,
			      Double_t la, Double_t la_raw,
			      Double_t ra, Double_t ra_raw ):
  fScBar(bar),fBarNum(barnum), fCaseNum(CaseNum), fLt(lt), fLt_raw(lt_raw), fRt(rt), fRt_raw(rt_raw),
  fLa(la), fLa_raw(la_raw),  fRa(ra), fRa_raw(ra_raw)
{;}


///////////////////////////////////////////////////////////////////////////////
