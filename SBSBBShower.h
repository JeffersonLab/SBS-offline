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

  explicit SBSBBShower( const char* name, const char* description = "",
                        THaApparatus* a = nullptr);

  virtual ~SBSBBShower();

  // Standard apparatus re-implemented functions
  virtual Int_t  CoarseProcess(TClonesArray& tracks);
  virtual Int_t  FineProcess(TClonesArray& tracks);
  virtual Int_t  ReadDatabase( const TDatime& date );
  virtual Int_t  DefineVariables( EMode mode = kDefine );
  virtual void   MakeCluster(Int_t nblk_size,SBSElement* blk);
  virtual void   MakeCluster(Int_t nblk_size);
  virtual void   AddToCluster(Int_t nc,SBSElement* blk);

  virtual void MakeMainCluster(Int_t icl=0); 
  virtual void Clear( Option_t* opt="" );

  //virtual Int_t SelectBestCluster(); //Potential override/extension of SBSCalorimeter::SelectBestCluster

  Int_t GetRowMax() {return GetRow();}
  Int_t GetColMax() {return GetCol();}

  
  Double_t EresMax() const {return fEres;}
  Double_t XresMax() const {return fXres;}
  Double_t YresMax() const {return fYres;}
  
  Double_t Eres(int i) const {return fE_cl_res[i];}
  Double_t Xres(int i) const {return fX_cl_res[i];}
  Double_t Yres(int i) const {return fY_cl_res[i];}

  SBSElement* GetElement(UInt_t i);
    
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
  Double_t fClusRadius; // radius for cluster search around the max 
  Int_t fClusBlockRadX;
  Int_t fClusBlockRadY;

  Double_t   fEres;        // [fNclus] Energy resolution of main cluster
  Double_t   fXres;        // [fNclus] x position (m) of main cluster
  Double_t   fYres;        // [fNclus] y position (m) of main cluster

  std::vector<Double_t>   fE_cl_res;        // [fNclus] Energy resolution of clusters
  std::vector<Double_t>   fX_cl_res;        // [fNclus] x position (m) of clusters
  std::vector<Double_t>   fY_cl_res;        // [fNclus] y position (m) of clusters 
  Double_t**  fEblk_cl;     // [fNclublk] Block energies composing main cluster

  // Useful derived quantities for internal usage.
  Double_t fThrADC;

  ClassDef(SBSBBShower,1)     // BigBite pre and shower class
};
#endif // SBSBBShower_h
