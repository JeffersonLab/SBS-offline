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
  : fNMaxElements(nmaxblk), fMaxElement(bar)
{
   fElements.clear();
   fElements.reserve(fNMaxElements);
   fElements.push_back(bar);
   fXmean = bar->GetElementPos();
   fYmean = bar->GetHitPos();
   fTmean = bar->GetMeanTime();
   fToTmean = bar->GetMeanToT();
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
  fToTmean = 0;
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
Bool_t SBSTimingHodoscopeCluster::AddElement(SBSTimingHodoscopeBar* bar) {
  //returns *true* if the element can be added to the cluster
  //returns *false* otherwise
  if (fMult<fNMaxElements) {
    //Add a criterion to reduce the time
    fElements.push_back(bar);
    fMult = fElements.size();
    //recompute the avergae of time, positions, and time over threshold
    //Should we use a weighted average? using time over threshold.
    //between ToT and amplitude:
    double sumtotmean_prev = fToTmean*(fMult-1);
    fToTmean = (fToTmean*(fMult-1)+bar->GetMeanToT());
    fXmean = (fXmean*sumtotmean_prev+bar->GetElementPos()*bar->GetMeanToT())/fToTmean;
    fYmean = (fYmean*sumtotmean_prev+bar->GetHitPos()*bar->GetMeanToT())/fToTmean;
    fTmean = (fTmean*sumtotmean_prev+bar->GetMeanTime()*bar->GetMeanToT())/fToTmean;
    fToTmean/= fMult;
    //fXmean = (fXmean*(fMult-1)+bar->GetElementPos())/fMult;
    //fYmean = (fYmean*(fMult-1)+bar->GetHitPos())/fMult;
    //fTmean = (fTmean*(fMult-1)+bar->GetMeanTime())/fMult;
    //we might want to see later for ToT weighted quantities
    //replace max element if the new element is the maximum one
    if(bar->GetMeanToT()>fMaxElement->GetMeanToT()){
      fMaxElement = bar;
    }
    return true;
  }else{
    return false;
  }
}

//_____________________________________________________________
void SBSTimingHodoscopeCluster::Clear( Option_t* ) {
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
