///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// SBSScintBar                                                               //
//                                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "SBSScintPMT.h"
#include "SBSScintBar.h"

// constructor 

SBSScintBar::SBSScintBar( Double_t x,Double_t y,Double_t z,
			  Double_t wx, Double_t wy, Double_t wz,
			  Double_t c, Double_t att,
			  Double_t lgain, Int_t lped, Double_t lres,
			  Double_t loff, Double_t lwalk,
			  Double_t rgain, Int_t rped, Double_t rres,
			  Double_t roff, Double_t rwalk, Int_t barnum,
			  Int_t llowlim, Int_t luplim, Double_t lwrapa, 
			  Int_t rlowlim, Int_t ruplim, Double_t rwrapa ) :
  fLPMT(lgain, lped, lres, loff, lwalk, this, barnum, 0, llowlim, luplim, lwrapa ),
  fRPMT(rgain, rped, rres, roff, rwalk, this, barnum, 1, rlowlim, ruplim, rwrapa ),
  fXPosPlane(x), fYPosPlane(y), fZPosPlane(z), fXWidth(wx), fYWidth(wy),  fZWidth(wz), 
  fc(c), fatt(att), fBarType(kSpecial), fBarNum(barnum), fBarNum_nd(0)
{
  // that's it!
}


//____________________________________________________________________




//______________________________________________________________________

SBSScintBar::~SBSScintBar()
{
  
}




ClassImp(SBSScintBar)


///////////////////////////////////////////////////////////////////////////////
