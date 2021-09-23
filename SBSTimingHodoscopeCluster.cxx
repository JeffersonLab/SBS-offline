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


//_____________________________________________________________
SBSTimingHodoscopeCluster::SBSTimingHodoscopeCluster(Int_t nmaxblk, SBSTimingHodoscopeBar* bar) 
: fNMaxElements(nmaxblk)
{
   fElements.clear();
   fElements.reserve(fNMaxElements);
   fElements.push_back(bar);
   fXmean = bar->GetElementPos();
   fYmean = bar->GetHitPos();
   fTmean = bar->GetMeanTime();
   fMult = 1;
}

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
    fXmean = (fXmean*(fMult-1)+bar->GetElementPos())/fMult;
    fYmean = (fYmean*(fMult-1)+bar->GetHitPos())/fMult;
    fTmean = (fTmean*(fMult-1)+bar->GetMeanTime())/fMult;
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
