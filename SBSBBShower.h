#ifndef SBSBBShower_h
#define SBSBBShower_h

////////////////////////////////////////////////////////////////////////////////
//
// SBSBBShower
//
// A sub-class of SBSCalorimeter for the Pre-shower and Shower in BigBite
//
////////////////////////////////////////////////////////////////////////////////

#include "SBSCalorimeter.h"

class SBSBBShower : public SBSCalorimeter {
public:
  SBSBBShower( const char* name, const char* description = "",
      THaApparatus* a = NULL);
  virtual ~SBSBBShower();

  // Standard apparatus re-implemented functions
  virtual Int_t  CoarseProcess(TClonesArray& tracks);
  virtual Int_t  FineProcess(TClonesArray& tracks);
  virtual Int_t  ReadDatabase( const TDatime& date );
  virtual Int_t  DefineVariables( EMode mode = kDefine );
  virtual void ClearEvent();

  Int_t GetRowMax() {return GetRow();}
  Int_t GetColMax() {return GetCol();}

  
  Float_t EresMax() {return fEres;}
  Float_t XresMax() {return fXres;}
  Float_t YresMax() {return fYres;}
  
  Float_t Eres(int i) {return fE_cl_res[i];}
  Float_t Xres(int i) {return fX_cl_res[i];}
  Float_t Yres(int i) {return fY_cl_res[i];}
    
  //two methods to set search region.
  void SetSearchRegion(int rowmin, int rowmax, int colmin, int colmax);

  void       LoadMCHitAt( Double_t x, Double_t y, Double_t E );
protected:
  bool fSearchRegion;
  Int_t fSearchRowmin;
  Int_t fSearchRowmax;
  Int_t fSearchColmin;
  Int_t fSearchColmax;

  bool fMultClus;// allow multiple clustering
  bool fMCdata;// easy way to enable/disable the use of MC data.

  // Cluster parameters Other parameters
  Float_t fClusRadius; // radius for cluster search around the max 
  Int_t fClusBlockRadX;
  Int_t fClusBlockRadY;

  Float_t   fEres;        // [fNclus] Energy resolution of main cluster
  Float_t   fXres;        // [fNclus] x position (m) of main cluster
  Float_t   fYres;        // [fNclus] y position (m) of main cluster

  std::vector<Float_t>   fE_cl_res;        // [fNclus] Energy resolution of clusters
  std::vector<Float_t>   fX_cl_res;        // [fNclus] x position (m) of clusters
  std::vector<Float_t>   fY_cl_res;        // [fNclus] y position (m) of clusters 
  Float_t**  fEblk_cl;     // [fNclublk] Block energies composing main cluster

  // Useful derived quantities for internal usage.
  Double_t fThrADC;

  ClassDef(SBSBBShower,1)     // BigBite pre and shower class
};
#endif // SBSBBShower_h
