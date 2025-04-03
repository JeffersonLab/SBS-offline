#ifndef SBSCDet_Hit_h
#define SBSCDet_Hit_h

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

class SBSCDet_Hit : public TObject {

 public:
  SBSCDet_Hit(); 
  SBSCDet_Hit( Int_t pmtnum, Int_t i, Int_t j, Int_t k,//Int_t TDC, Int_t ToT
		    Double_t x, Double_t y, Double_t z, 
		    Double_t le, Double_t te, Double_t tot );
  virtual ~SBSCDet_Hit() {}
  
  virtual void    Clear( Option_t* opt);

  //void       Show(FILE * fout1);
  //void       Show(FILE * fout1, FILE * fout2);
  
  Int_t      GetPMTNum()   const {return fPMTNum;}
  Int_t      GetRow()      const {return fRow;}
  Int_t      GetCol()      const {return fCol;}
  Int_t      GetLayer()    const {return fLayer;}

  Int_t      GetTrackIndex() const { return fTrackIndex; }
  Double_t    GetX()        const {return fX;}
  Double_t    GetY()        const {return fY;}
  Double_t    GetZ()        const {return fZ;}
  Double_t    GetTDC_LE()     const {return fTDC_LE;}
  Double_t    GetTDC_TE()     const {return fTDC_TE;}
  Double_t    GetToT()     const {return fToT;}
  //Int_t      GetFlag()     const {return fFlag;}
  //Int_t      GetVeto()     const {return fVeto;}
  
  void       SetPMTNum( Int_t pmtnum ) {fPMTNum = pmtnum;}
  void       SetRow( Int_t i )         {fRow = i;}
  void       SetCol( Int_t j )         {fCol = j;}
  void       SetLayer( Int_t k )       {fLayer = k;}


  //void SetClustIndex( Int_t icl ){ fClustIndex = icl; }
  void SetTrackIndex( Int_t itr ){ fTrackIndex = itr; }
  void       SetX( Double_t x )         {fX = x;}
  void       SetY( Double_t y )         {fY = y;}
  void       SetZ( Double_t z )         {fZ = z;}
  void       SetTDC_LE( Int_t le )       {fTDC_LE = le;}
  void       SetTDC_TE( Int_t te )       {fTDC_TE = te;}
  void       SetToT( Int_t tot )       {fToT = tot;}
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
  Int_t     fLayer;      // Hit column number in PMT matrix
  
  //Int_t     fClustIndex; //index of cluster to which this hit belongs
  Int_t     fTrackIndex; //index of track to which this hit is associated
  
  Double_t   fX;      // Hit X position in PMT matrix
  Double_t   fY;      // Hit Y position in PMT matrix
  Double_t   fZ;      // Hit Z position in PMT matrix
  Double_t   fTDC_LE;      // LE Time from TDC
  Double_t   fTDC_TE;      // TE Time from TDC
  Double_t   fToT;      // Time over Threshold from TDC
  //Int_t     fFlag;   // ?
  //Int_t     fVeto;   // ?
  
  //Bool_t tdcr_set;// lead TDC is set
  //Bool_t tdcf_set;// fall TDC is set
  
  ClassDef(SBSCDet_Hit,0)   //A hit in the RICH
};

#endif






