///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// SBSTimingHodoscopeCluster                                                 //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include "SBSTimingHodoscopeCluster.h"
#include <iostream>
#include <DataType.h>

ClassImp(SBSTimingHodoscopeCluster)   // Generic shower cluster class

using namespace std;


/*
//_____________________________________________________________
SBSTimingHodoscopeCluster::SBSTimingHodoscopeCluster(Int_t nmaxblk, SBSTimingHodoscope* bar) 
: fNMaxElements(nmaxblk), fMaxElement(0)
{
   fElements.clear();
   fElements.reserve(fNMaxElements);
   fElements.push_back(bar);
   fXmean = bar->GetX();
   fYmean = bar->GetY();
   fTmean = bar->GetT();
   fMult = 1;
}
*/

//_____________________________________________________________
SBSTimingHodoscopeCluster::SBSTimingHodoscopeCluster(Int_t nmaxblk) 
: fNMaxElements(nmaxblk)
{
  fElements.clear();
  fElements.reserve(fNMaxElements);
  fXmean = 0;
  fYmean = 0;
  fTmean = 0;
  fMult = 0;
}


//_____________________________________________________________
SBSTimingHodoscopeCluster::SBSTimingHodoscopeCluster() {
  fXmean = kBig;
  fYmean = kBig;
  fTmean = kBig;
  fMult = 0;
}

//_____________________________________________________________
SBSTimingHodoscopeCluster::~SBSTimingHodoscopeCluster()
{ 
}

//_____________________________________________________________
void SBSTimingHodoscopeCluster::AddElement(SBSTimingHodoscopeBar* bar) {
  if (fMult<fNMaxElements) {
    fElements.push_back(bar);
    fMult = fElements.size();
    /*
    block->SetStat(1);
    fX = (fX*fE + block->GetX()*block->GetE()) / (fE+block->GetE());
    fY = (fY*fE + block->GetY()*block->GetE()) / (fE+block->GetE());
    fE += block->GetE();
        if(block->GetE() > fEblk) {
          fEblk = block->GetE();
          fRow = block->GetRow();
          fCol = block->GetCol();
        }
        // Keep a pointer to the element with the highest energy
        if(!fMaxElement) {
          fMaxElement = block;
        } else if ( block->GetE() > fMaxElement->GetE() ) {
          fMaxElement = block;
        }
    */
  }
}

//_____________________________________________________________
void SBSTimingHodoscopeCluster::ClearEvent() {
  fMult=0;fXmean=fYmean=fTmean=0.;
  fElements.clear();
}

//_____________________________________________________________
SBSTimingHodoscopeBar* SBSTimingHodoscopeCluster::GetElement(UInt_t i)
{
  SBSTimingHodoscopeBar* bar=0;
  if(i < fElements.size()) bar = fElements[i];
  return bar;
}



ClassImp(SBSTimingHodoscopeCluster)
