///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// SBSScintHit                                                               //
//                                                                           //
// Class representing a single hit for a scint                               //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "SBSScintHit.h"


//_______________________________________________________________________________
SBSScintHit::SBSScintHit(  const SBSScintBar* bar,Int_t planenum, Int_t barnum,
			   Double_t ypos,
			   Double_t Tof, Double_t HitEnergy, Double_t Tdiff)
{
  Clear();

  fScBar = const_cast<SBSScintBar*>(bar);
  fPlaneNum = planenum;
  fBarNum = barnum;
  fHitYPos = ypos;
  fHitTOF = Tof;
  fHitEdep = HitEnergy;
  fTdiff = Tdiff;
  
  if (bar) {
    fHitXPos = bar->GetXPos();
    fHitZPos = bar->GetZPos();
    fBarNum_nd = bar->GetBarNum_nd();
  } else {
    fHitXPos = -1.e5;
    fHitZPos = -1.e5;
  }
}

//_______________________________________________________________________________
SBSScintHit::SBSScintHit(const SBSScintHit* pScHit)
{
  Clear();
  CopyScintHit(pScHit);
}

//_______________________________________________________________________________
SBSScintHit::SBSScintHit(const SBSScintHit* pScHit, Int_t clusternum)
{
  Clear();
  fClusterNum=clusternum;
  CopyScintHit(pScHit);
}

//_______________________________________________________________________________
SBSScintHit::SBSScintHit(const SBSScintHit* pScHit, Int_t planenum, Int_t Bar_nd)
{
  Clear();
  fPlaneNum=planenum;
  fBarNum_nd=Bar_nd;
  CopyScintHit(pScHit);
}

//___________________________________________________________________________
void SBSScintHit::Clear(Option_t *) {
  fScBar = nullptr;
  fPlaneNum = -1;
  fBarNum = 0;
  fBarNum_nd = -1;
  fOrder = 0;
  fClusterNum = -1;
  
  fHitXPos = fHitYPos = fHitZPos = fHitTOF = fHitEdep = fTdiff = 0;
}
  
//___________________________________________________________________________
SBSScintHit::~SBSScintHit()
{
  fScBar = 0;
}

//_____________________________________________________________________________
Int_t SBSScintHit::CopyScintHit(const SBSScintHit* pScHit)
{
  if(pScHit){
    fScBar = pScHit->fScBar;
    fBarNum = pScHit->fBarNum;
    if(fBarNum_nd<0) fBarNum_nd = pScHit->fBarNum_nd;
    fHitXPos = pScHit->fHitXPos;
    fHitYPos = pScHit->fHitYPos;
    fHitZPos = pScHit->fHitZPos;
    fHitTOF = pScHit->fHitTOF;
    fHitEdep = pScHit->fHitEdep;
    fTdiff = pScHit->fTdiff;
    fOrder = pScHit->fOrder;
    if(fClusterNum<0) fClusterNum = pScHit->fClusterNum;
    if(fPlaneNum<0) fPlaneNum = pScHit->fPlaneNum;
  }

  return 0;
}

//_____________________________________________________________________________
Int_t SBSScintHit::Compare(const TObject* obj) const
{
  // Compare two hits via plane number, and then bar number,
  // and finally time
  // Returns -1 if this < obj, 0 if this==obj, and 1 if this>obj
  const SBSScintHit *h = static_cast<const SBSScintHit*>(obj);

  // Side (Left-right) comparison
  if (fPlaneNum < h->fPlaneNum) return -1;
  if (fPlaneNum > h->fPlaneNum) return  1;

  // Bar-number comparison
  if (fBarNum < h->fBarNum) return -1;
  if (fBarNum > h->fBarNum) return  1;
  
  // finally, time comparison, earlier wins
  if (fHitTOF < h->fHitTOF) return -1;
  if (fHitTOF > h->fHitTOF) return  1;

  return 0;
}

//////////////////////////////////////////////////////////////////////////////

ClassImp(SBSScintHit)
