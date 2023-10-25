#ifndef ROOT_SBSCherenkov_ClusterList
#define ROOT_SBSCherenkov_ClusterList

//****************************************************************
// The header file for the "ClusterRICH.C"
//****************************************************************

#include "TObject.h"
#include "TList.h"

class THaTrack;
class TClonesArray;

#include <iostream>
#include <stdio.h>

// Data class contains the data of a single PMT
// (i,j) coordinate of a pad; (x,y) coordinated of a pad; ADC value
// and a Veto to point it cannot be a local minimum inside the cluster
// useful to tell overlapping clusters apart). 
// Above all it contains a module to compare adc values. This allows 
// the ordering according the adc values itself of the components of a cluster.

class SBSCherenkov_Hit : public TObject {

 public:
  SBSCherenkov_Hit(); 
  SBSCherenkov_Hit( Int_t pmtnum, Int_t i, Int_t j, //Int_t TDC, Int_t ToT
		    Double_t x, Double_t y, Double_t t, Double_t a );
  virtual ~SBSCherenkov_Hit() {}
  
  virtual void    Clear( Option_t* opt);

  //void       Show(FILE * fout1);
  //void       Show(FILE * fout1, FILE * fout2);
  
  Int_t      GetPMTNum()   const {return fPMTNum;}
  Int_t      GetRow()      const {return fRow;}
  Int_t      GetCol()      const {return fCol;}

  Int_t      GetClustIndex() const { return fClustIndex; }
  Int_t      GetTrackIndex() const { return fTrackIndex; }
  //Int_t      GetTDC()      const {return fTDC;}
  //Int_t      GetToT()      const {return fToT;}
  //Int_t      GetADC()      const {return fADC;}
  Double_t    GetX()        const {return fX;}
  Double_t    GetY()        const {return fY;}
  Double_t    GetTime()     const {return fTime;}
  Double_t    GetAmp()      const {return fAmp;}
  //Int_t      GetFlag()     const {return fFlag;}
  //Int_t      GetVeto()     const {return fVeto;}
  
  void       SetPMTNum( Int_t pmtnum ) {fPMTNum = pmtnum;}
  void       SetRow( Int_t i )         {fRow = i;}
  void       SetCol( Int_t j )         {fCol = j;}

  void SetClustIndex( Int_t icl ){ fClustIndex = icl; }
  void SetTrackIndex( Int_t itr ){ fTrackIndex = itr; }
  //void       SetTDC( Int_t TDC )       {fTDC = TDC;}
  //void       SetToT( Int_t ToT )       {fToT = ToT;}
  //void       SetADC( Int_t ADC )       {fADC = ADC;}
  void       SetX( Double_t x )         {fX = x;}
  void       SetY( Double_t y )         {fY = y;}
  void       SetTime( Double_t t )      {fTime = t;}
  void       SetAmp( Double_t a )       {fAmp = a;}
  //void       SetFlag( Int_t Flag )     {fFlag = Flag;}
  //void       SetVeto( Int_t Veto )     {fVeto = Veto;}

  virtual Int_t   Compare( const TObject* ) const;
  virtual Bool_t  IsSortable() const { return kTRUE; }
  
  //Bool_t     TDC_risSet() {return tdcr_set;}
  //Bool_t     TDC_fisSet() {return tdcf_set;}

private:
  Int_t     fPMTNum; // Hit PMT number
  Int_t     fRow;      // Hit row number in PMT matrix
  Int_t     fCol;      // Hit column number in PMT matrix
  
  Int_t     fClustIndex; //index of cluster to which this hit belongs
  Int_t     fTrackIndex; //index of track to which this hit is associated
  
  //Int_t     fTDC;  // Hit rise time TDC
  //Int_t     fToT;  // Hit time over threshold
  //Int_t     fADC;    // Hit ADC (if available)
  Double_t   fX;      // Hit X position in PMT matrix
  Double_t   fY;      // Hit Y position in PMT matrix
  Double_t   fTime;      // Time from TDC
  Double_t   fAmp;      // Amplitude from ADC or ToT
  //Int_t     fFlag;   // ?
  //Int_t     fVeto;   // ?
  
  //Bool_t tdcr_set;// lead TDC is set
  //Bool_t tdcf_set;// fall TDC is set
  
  ClassDef(SBSCherenkov_Hit,0)   //A hit in the RICH
};


// --------------------------------------------------------------

// Cluster: List of Hits of elements that make up a cluster in GRINCH. 
// Specifically, in the GRINCH, a cluster is a list of adjacent PMT hits;
// 

class SBSCherenkov_Cluster : public TObject {
 public:
  SBSCherenkov_Cluster();
  SBSCherenkov_Cluster( const SBSCherenkov_Cluster& rhs );
  SBSCherenkov_Cluster& operator=( const SBSCherenkov_Cluster& rhs );

  ~SBSCherenkov_Cluster() { delete fHitList; }
  
  void    Clear( Option_t* opt="" );
  void    Insert( SBSCherenkov_Hit* theHit );
  void    Remove( SBSCherenkov_Hit* theHit );
  
  Bool_t  IsNeighbor( const SBSCherenkov_Hit* theHit,  Double_t par);
  
  Int_t   GetNHits()   const { return fHitList->GetSize(); }
  TList*  GetHitList()      { return fHitList; }
  
  Double_t GetXcenter() const { return fXcenter; }
  Double_t GetYcenter() const  { return fYcenter; }
  Double_t GetXcenter_w() const { return fXcenter_w; }
  Double_t GetYcenter_w() const  { return fYcenter_w; }
  Double_t GetCharge()  const     { return fCharge; }
  
  Double_t GetMeanTime() const { return fMeanTime; }
  Double_t GetMeanAmp() const { return fMeanAmp; }
  Double_t GetTimeRMS() const  { return fTimeRMS; }
  Double_t GetAmpRMS() const  { return fAmpRMS; }

  Int_t GetTrackIndex() const { return fTrackIndex; }
  Int_t GetMirrorIndex() const { return fMirrorIndex; }

  Double_t GetDX() const { return fTrackMatch_dx; }
  Double_t GetDY() const { return fTrackMatch_dy; }
  
  THaTrack* GetTrack() const            { return fTrack; }

  void    SetXcenter(Double_t xcenter)    { fXcenter = xcenter; }
  void    SetYcenter(Double_t ycenter)     { fYcenter = ycenter; }
  void    SetXcenter_w(Double_t xcenter_w)  { fXcenter_w = xcenter_w; }
  void    SetYcenter_w(Double_t ycenter_w)   { fYcenter_w = ycenter_w; }
  void    SetCharge(Double_t charge)          { fCharge = charge; }
  
  void    SetMeanTime(Double_t meantdc) { fMeanTime = meantdc; }
  void    SetMeanAmp(Double_t meantot) { fMeanAmp = meantot; }
  void    SetTimeRMS(Double_t tdcrms)   { fTimeRMS = tdcrms; }
  void    SetAmpRMS(Double_t totrms)   { fAmpRMS = totrms; }
  
  void    SetTrack( THaTrack* track ) { fTrack = track; }
  void    SetTrackIndex( Int_t itr ){ fTrackIndex = itr; }
  void    SetMirrorIndex( Int_t imirr ){ fMirrorIndex = imirr; }

  void    SetDX( Double_t dx ){ fTrackMatch_dx = dx; }
  void    SetDY( Double_t dy ){ fTrackMatch_dy = dy; }

  void    MergeCluster( const SBSCherenkov_Cluster& rhs );
  
 private:
  TList*     fHitList;   //List of hits belonging to this cluster
  
  Double_t    fXcenter;   // X mean of all hits in the list
  Double_t    fYcenter;   // Y mean of all hits in the list
  Double_t    fXcenter_w; // Weighted X mean : (Sum of x*adc)/(sum adc) of all hits in the list
  Double_t    fYcenter_w; // Weighted Y mean : (Sum of y*adc)/(sum adc) of all hits in the list
  Double_t    fCharge;    // Sum of ADCs
  
  Double_t    fMeanTime;   // Mean rising time of all hits in the list
  Double_t    fMeanAmp;   // Mean time over threshold of all hits in the list
  Double_t    fTimeRMS;    // Rising time RMS of all hits in the list
  Double_t    fAmpRMS;    // time-over-threshold RMS of all hits in the list

  Bool_t     fTrackMatch;// true if a track can be matched to the cluster
  THaTrack*  fTrack;     // Track best associated with this cluster
  
  Int_t      fTrackIndex; //ROOT-tree friendly variable to store track index with which this cluster is associated
  Int_t      fMirrorIndex; //Store mirror index for track-matched clusters

  Double_t   fTrackMatch_dx; //x GRINCH - x predicted from track;
  Double_t   fTrackMatch_dy; //y GRINCH - y predicted from track;
  
  ClassDef(SBSCherenkov_Cluster,0)  //A cluster of hits in the GRINCH
};

#endif


















