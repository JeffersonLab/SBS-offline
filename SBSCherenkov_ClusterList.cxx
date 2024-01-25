//*-- Author :    Ole Hansen   26 March 2001; completed by Guido Urciuoli;
//*-- last change: 9 January 2004.

//////////////////////////////////////////////////////////////////////////
//
// SBSCherenkov_ClusterList
//
//////////////////////////////////////////////////////////////////////////


#include "SBSCherenkov_ClusterList.h"
#include "TClonesArray.h"
#include "TMath.h"
#include <DataType.h>

#include <fstream>
#include <cstdio>
#include <cstring>
#include <cmath>

using namespace std;

//_____________________________________________________________________________
SBSCherenkov_Hit::SBSCherenkov_Hit(): 
  fPMTNum(-1), fRow(-1), fCol(-1), 
  //fTDC(-1), fToT(-1), fADC(-1),
  fX(kBig), fY(kBig), fTime(kBig), fAmp(kBig)
  //fFlag(0), fVeto(0) 
{
  fClustIndex = -1;
  fTrackIndex = -1;
} 

//_____________________________________________________________________________
SBSCherenkov_Hit::SBSCherenkov_Hit( Int_t pmtnum, Int_t i, Int_t j, 
				    Double_t x, Double_t y, Double_t t, Double_t a ):
  fPMTNum(pmtnum), fRow(i), fCol(j), 
  //fTDC(TDC), fToT(ToT), fADC(ADC), 
  fX(x), fY(y), fTime(t), fAmp(a)
  //fFlag(0), fVeto(0), tdcr_set(false), tdcf_set(false) 
{
  fClustIndex = -1;
  fTrackIndex = -1;
}

//_____________________________________________________________________________
Int_t SBSCherenkov_Hit::Compare( const TObject* theOtherHit ) const
{
  if (fTime < static_cast<const SBSCherenkov_Hit*>( theOtherHit )->fTime ||
      fAmp < static_cast<const SBSCherenkov_Hit*>( theOtherHit )->fAmp)
    return -1;
  if (fTime > static_cast<const SBSCherenkov_Hit*>( theOtherHit )->fTime ||
      fAmp > static_cast<const SBSCherenkov_Hit*>( theOtherHit )->fAmp)
    return +1;
  else
    return 0;
}

void SBSCherenkov_Hit::Clear( Option_t *opt ){ 
  fPMTNum = -1;
  fRow = -1;
  fCol = -1;
  fClustIndex = -1;
  fTrackIndex = -1;
  fX = kBig;
  fY = kBig;
  fTime = kBig;
  fAmp = kBig;
}

/*
//_____________________________________________________________________________
void SBSCherenkov_Hit::Show(FILE * fout1, FILE* fout2) 
{
  // FIXME comment on fprintf(fout1...) should be changed accordingly 
  // when ntuple will be created 
  fprintf(fout2," PMT num ");
  fprintf(fout2,"%4d",fPMTNum);
  fprintf(fout2," Row, Col ");
  fprintf(fout2,"%4d,%4d",fRow,fCol);
  fprintf(fout2," ; X, Y ");
  fprintf(fout2,"%4f,%4f",fX,fY);
  fprintf(fout2,"; fTDC, fToT ");
  fprintf(fout2,"%4d,%4d",fTDC, fToT);
  fprintf(fout2,"; fADC = ");
  fprintf(fout2,"%4d",fADC);
  fprintf(fout2,"\n");
  //  Show(fout1);
}

void SBSCherenkov_Hit::Show(FILE * fout1) 
{
  // FIXME comment on fprintf(fout1...) should be changed accordingly 
  // when ntuple will be created
  fprintf(fout1," %4d",fPMTNum);
  fprintf(fout1," %4d %4d",fRow,fCol);
  fprintf(fout1,"% 4f %4f",fX,fY);
  fprintf(fout1," %4d %4d \n",fTime, fAmp);
  //fprintf(fout1," %4d \n",fADC);
}
*/

//=============================================================================
// SBSCherenkov_Cluster
//=============================================================================

//_____________________________________________________________________________
SBSCherenkov_Cluster::SBSCherenkov_Cluster() : // f(0)
  fXcenter(0), fYcenter(0),
  fXcenter_w(0), fYcenter_w(0), fCharge(0),
  fMeanTime(0), fMeanAmp(0),
  fTimeRMS(0), fAmpRMS(0),
  fTrackMatch(false), fTrack(0)
{
  fTrackIndex = -1;
  fMirrorIndex = -1;

  fTrackMatch_dx = kBig;
  fTrackMatch_dy = kBig; 
  
  fHitList = new TList(); 
}

//_____________________________________________________________________________
SBSCherenkov_Cluster::SBSCherenkov_Cluster( const SBSCherenkov_Cluster& rhs ) : // f(rhs.f)
  TObject(rhs), fXcenter(rhs.fXcenter), fYcenter(rhs.fYcenter),
  fXcenter_w(rhs.fXcenter_w), fYcenter_w(rhs.fYcenter_w), fCharge(rhs.fCharge), 
  fMeanTime(rhs.fMeanTime), fMeanAmp(rhs.fMeanAmp),
  fTimeRMS(rhs.fTimeRMS), fAmpRMS(rhs.fAmpRMS),
  fTrackMatch(rhs.fTrackMatch), fTrack(rhs.fTrack), fTrackIndex(rhs.fTrackIndex),
  fMirrorIndex(rhs.fMirrorIndex), fTrackMatch_dx(rhs.fTrackMatch_dx), fTrackMatch_dy(rhs.fTrackMatch_dy)
{
  fHitList = new TList();
  if( rhs.fHitList && (rhs.fHitList->GetSize() > 0 )) {
    TIter next( rhs.fHitList );
    while( SBSCherenkov_Hit* pHit = static_cast<SBSCherenkov_Hit*>( next() ))
      fHitList->AddLast(pHit);
  }
}

//_____________________________________________________________________________
SBSCherenkov_Cluster& SBSCherenkov_Cluster::operator=( const SBSCherenkov_Cluster& rhs ) // f = rhs.f;
{
  // Assignment operator
  if( this != &rhs ) {
    fXcenter = rhs.fXcenter;
    fYcenter = rhs.fYcenter;
    fXcenter_w = rhs.fXcenter_w;
    fYcenter_w = rhs.fYcenter_w;
    fCharge = rhs.fCharge;
    fMeanTime = rhs.fMeanTime;
    fMeanAmp = rhs.fMeanAmp;
    fTimeRMS = rhs.fTimeRMS;
    fAmpRMS = rhs.fAmpRMS;
    fTrackMatch = rhs.fTrackMatch;
    fTrack = rhs.fTrack;

    fTrackIndex = rhs.fTrackIndex;
    fMirrorIndex = rhs.fMirrorIndex;

    fTrackMatch_dx = rhs.fTrackMatch_dx;
    fTrackMatch_dy = rhs.fTrackMatch_dy;
    
    if( !fHitList )
      fHitList = new TList;
    else
      fHitList->Clear("nodelete");
    if( rhs.fHitList && (rhs.fHitList->GetSize() > 0 )) {
      TIter next( rhs.fHitList ); 
      while( SBSCherenkov_Hit* pHit = static_cast<SBSCherenkov_Hit*>( next() ))
	fHitList->AddLast(pHit);
    }
  }
  return *this;
}

//_____________________________________________________________________________
void SBSCherenkov_Cluster::MergeCluster( const SBSCherenkov_Cluster& rhs )
{//adds the cluster in argument to the cluster which the method is applied to 
  if( !fHitList ) fHitList = new TList;
  Int_t list1size = fHitList->GetSize();
  Int_t list2size = 0;
  
  if(list1size==0){
    *this = rhs;
    //return *this;
  }
  
  if( rhs.fHitList && (rhs.fHitList->GetSize() > 0 )) {
    list2size = rhs.fHitList->GetSize();
    TIter next( rhs.fHitList ); 
    while( SBSCherenkov_Hit* pHit = static_cast<SBSCherenkov_Hit*>( next() ))
      fHitList->AddLast(pHit);
  }
  
  if(list2size==0)return;// return *this;
  
  fXcenter = (fXcenter*((Double_t)(list1size))+ rhs.fXcenter*((Double_t)(list2size)))/
    ((Double_t)(list1size+list2size));
  fYcenter = (fYcenter*((Double_t)(list1size))+ rhs.fYcenter*((Double_t)(list2size)))/
    ((Double_t)(list1size+list2size));
  
  fXcenter_w = (fXcenter_w*fCharge+rhs.fXcenter_w*rhs.fCharge)/(fCharge+rhs.fCharge);
  fYcenter_w = (fYcenter_w*fCharge+rhs.fYcenter_w*rhs.fCharge)/(fCharge+rhs.fCharge);
  
  fCharge += rhs.fCharge;
  
  fMeanTime = (fMeanTime*((Double_t)list1size)+rhs.fMeanTime*((Double_t)list2size))/
    ((Double_t)(list1size+list2size));
  fMeanAmp = (fMeanAmp*((Double_t)list1size)+rhs.fMeanAmp*((Double_t)list2size))/
    ((Double_t)(list1size+list2size));
  
  fTimeRMS = sqrt( (pow(fTimeRMS, 2)*((Double_t)list1size) + pow(rhs.fTimeRMS, 2)*((Double_t)list2size) )/((Double_t)(list1size+list2size)) );
  fAmpRMS = sqrt( (pow(fAmpRMS, 2)*((Double_t)list1size) + pow(rhs.fAmpRMS, 2)*((Double_t)list2size) )/((Double_t)(list1size+list2size)) );
  //return *this;
}

//_____________________________________________________________________________
void SBSCherenkov_Cluster::Clear( Option_t* opt ) // f = 0;
{
  //Reset the cluster to an empty state.

  //Don't delete the hits, just clear the internal list.
  if( fHitList ) 
    fHitList->Clear("nodelete");

  // Full clear
  if( opt && opt[0] == 'F' ) {
    fXcenter = 0;
    fYcenter = 0;
    fXcenter_w = 0;
    fYcenter_w = 0;
    fCharge = 0;
    fMeanTime = 0;
    fMeanAmp = 0;
    fTimeRMS = 0;
    fAmpRMS = 0;
    fTrackMatch = false;
    fTrack = 0;
    fTrackIndex=-1;
    fMirrorIndex = -1;
    fTrackMatch_dx = kBig;
    fTrackMatch_dy = kBig;
  } else {
    // Fast clear for clearing TClonesArrays of clusters
    // Needs to be deleted since a TList allocates memory
    // FIXME: performance issue?
    delete fHitList;
    fHitList = 0;
  }  
}

//_____________________________________________________________________________
void SBSCherenkov_Cluster::Insert( SBSCherenkov_Hit* theHit )
{
  //Add a hit to the cluster
    
  if( !fHitList ) fHitList = new TList;
  fHitList->AddLast( theHit );
  
  Int_t listnewsize = fHitList->GetSize();
  fXcenter = (fXcenter*((Double_t)(listnewsize-1))+theHit->GetX())/((Double_t)listnewsize);
  fYcenter = (fYcenter*((Double_t)(listnewsize-1))+theHit->GetY())/((Double_t)listnewsize);
  
  fXcenter_w = fXcenter_w*fCharge;
  fYcenter_w = fYcenter_w*fCharge;
  fCharge+= theHit->GetAmp();
  fXcenter_w+= theHit->GetAmp()*theHit->GetX();
  fYcenter_w+= theHit->GetAmp()*theHit->GetY();
  fXcenter_w = fXcenter_w/fCharge;
  fYcenter_w = fYcenter_w/fCharge;
  
  fMeanTime = (fMeanTime*((Double_t)(listnewsize-1))+theHit->GetTime())/((Double_t)listnewsize);
  fMeanAmp = (fMeanAmp*((Double_t)(listnewsize-1))+theHit->GetAmp())/((Double_t)listnewsize);
  fTimeRMS = sqrt((pow(fTimeRMS, 2)*((Double_t)(listnewsize-1))+ pow(theHit->GetTime(), 2))/
  			((Double_t)listnewsize));
  fAmpRMS = sqrt((pow(fAmpRMS, 2)*((Double_t)(listnewsize-1))+ pow(theHit->GetAmp(), 2))/
  			 ((Double_t)listnewsize));
}

//_____________________________________________________________________________
void SBSCherenkov_Cluster::Remove( SBSCherenkov_Hit* theHit )
{
  
  if( !fHitList ) return;//if list does not exist, nothing to do
  if(fHitList->IndexOf(theHit)<0) return;//if hit not in list, nothing to do
    
  Int_t listnewsize = fHitList->GetSize();
  fXcenter = (fXcenter*((Double_t)(listnewsize+1))-theHit->GetX())/((Double_t)listnewsize);
  fYcenter = (fYcenter*((Double_t)(listnewsize+1))-theHit->GetY())/((Double_t)listnewsize);
  
  fXcenter_w = fXcenter_w*fCharge;
  fYcenter_w = fYcenter_w*fCharge;
  fCharge-= theHit->GetAmp();
  fXcenter_w-= theHit->GetAmp()*theHit->GetX();
  fYcenter_w-= theHit->GetAmp()*theHit->GetY();
  fXcenter_w = fXcenter_w/fCharge;
  fYcenter_w = fYcenter_w/fCharge;
  
  fMeanTime = (fMeanTime*((Double_t)(listnewsize+1))-theHit->GetTime())/((Double_t)listnewsize);
  fMeanAmp = (fMeanAmp*((Double_t)(listnewsize+1))-theHit->GetAmp())/((Double_t)listnewsize);
  fTimeRMS = sqrt((pow(fTimeRMS, 2)*((Double_t)(listnewsize+1))-pow(theHit->GetTime(), 2))/
  			((Double_t)listnewsize));
  fAmpRMS = sqrt((pow(fAmpRMS, 2)*((Double_t)(listnewsize+1))-pow(theHit->GetAmp(), 2))/
  			 ((Double_t)listnewsize));
}

//_____________________________________________________________________________
Bool_t SBSCherenkov_Cluster::IsNeighbor(const SBSCherenkov_Hit* theHit, Double_t par)
{
  //cout << "SBSCherenkov_Cluster::IsNeighbor: " << endl;
  //cout << fHitList << endl;
  Double_t dx,dy,dist2;
  if( !fHitList )//{
    return 0;
  //}else{cout << " size ? " << fHitList->GetLast()+1 << endl;}
  TIter next( fHitList );

  //if(theHit){
  //cout << " theHit pos " << theHit->GetX() << " " << theHit->GetY() << endl;
  //}else{cout << "SBSCherenkov_Cluster::IsNeighbor : theHit is 0" << endl; }
  while( SBSCherenkov_Hit* pHit = static_cast<SBSCherenkov_Hit*>( next() )) {
    if(pHit){
      //cout << " pHit pos: " <<endl;
      //cout << pHit->GetX() << endl; 
      dx   = theHit->GetX() - pHit->GetX();
      //cout<< " " << pHit->GetY() << endl;
      dy   = theHit->GetY() - pHit->GetY();
      dist2 = dx*dx + dy*dy;
      assert(dist2>=0);
      if( sqrt(dist2)<par )
	return true;
    }
  }
  return false;
}

//_____________________________________________________________________________
ClassImp(SBSCherenkov_Hit)
ClassImp(SBSCherenkov_Cluster)


