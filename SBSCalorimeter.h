#ifndef SBSCalorimeter_h
#define SBSCalorimeter_h

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// SBSCalorimeter                                                            //
//                                                                           //
// A generic calorimeter class which could contain the following types       //
// of data:                                                                  //
//  - Single-valued ADC                                                      //
//  - Single-valued TDC                                                      //
//  - Multi-valued ADC                                                       //
//  - Multi-valued ADC + Single-valued TDC                                   //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "SBSGenericDetector.h"
//#include "THaShower.h"
#include "SBSCalorimeterCluster.h"
#include "TRotation.h"
#include "TVector3.h"
#include "THaDetMap.h"


class SBSCalorimeter : public SBSGenericDetector {
  //class SBSCalorimeter : public THaShower {

public:
  SBSCalorimeter( const char* name, const char* description = "",
      THaApparatus* a = NULL);
  virtual ~SBSCalorimeter();

  virtual void ClearEvent();
  virtual void ClearOutputVariables();

  // Standard apparatus re-implemented functions
  virtual Int_t      CoarseProcess(TClonesArray& tracks);
  virtual Int_t      FineProcess(TClonesArray& tracks);

protected:

  virtual Int_t  ReadDatabase( const TDatime& date );
  virtual Int_t  DefineVariables( EMode mode = kDefine );

  // Configuration
  Int_t  fNclublk;      ///< Max number of blocks composing a cluster
  Int_t  fNclubr;       ///< Max number of row-blocks composing a cluster
  Int_t  fNclubc;       ///< Max number of col-blocks composing a cluster

  // Mapping (see also fDetMap)
  UShort_t   fChanMapStart; ///< Starting number for block number (i.e. 0 or 1)
  std::vector<std::vector<UShort_t> > fChanMap; //< Maps modules in THaDetMap to calorimeter block number

  // Output variables
  // Cluster output
  std::vector<Float_t> fEclus;    //< [] Energy (MeV) of clusters
  std::vector<Float_t> fEclus_c;  //< [] Corrected energy (MeV) of clusters
  std::vector<Float_t> fEclusBlk; //< [] Energy (MeV) of block with highest E in clusters
  std::vector<Float_t> fEclusBlk_c;//< [] Corrected energy (MeV) of block with highest E in clusters
  std::vector<Float_t> fXclus;    //< [] x position (m) of clusters
  std::vector<Float_t> fYclus;    //< [] y position (m) of clusters
  std::vector<Float_t> fNblkclus; //< [] number of blocks in clusters
  std::vector<Int_t> fRowblkclus; //< [] row of blocks in clusters
  std::vector<Int_t> fColblkclus; //< [] col of blocks in clusters
  // Single valued variables for the cluster with highest energy
  Float_t fE;                     //< Energy (MeV) of cluster with highest energy
  Float_t fE_c;                   //< Corrected energy (MeV) of cluster with highest energy
  Float_t fEblk;                  //< Energy (MeV) of largest block in cluster with highest energy
  Float_t fEblk_c;                //< Corrected energy (MeV) of largest block in cluster with highest energy
  Float_t fX;                     //< x position (m) of cluster with highest energy
  Float_t fY;                     //< y position (m) of cluster with highest energy
  Float_t fNblk;                  //< number of blocks in cluster with highest energy
  // If only storing highest energy cluster, info on other clusters may be useful
  Float_t fNclus;                 //< Number of clusters meeting threshold
  Int_t   fRowblk;                //< row of block with highest energy in highest E cluster
  Int_t   fColblk;                //< col of block with highest energy in highest E cluster

  // Blocks, where the grid is just for easy axis to the blocks by row,col,layer
  // Clusters for this event
  std::vector<SBSCalorimeterCluster*> fClusters; //[] Cluster

  // Other parameters
  Float_t    fEmin;         //< Minimum energy for a cluster center
  Int_t fMaxNclus;          //< Maximum number of clusters to store
  Bool_t fCoarseProcessed;  //< Was CourseProcessed already called?
  Bool_t fFineProcessed;    //< Was fine processed already called

  //Gain correction 
  Float_t   fConst;     // const from gain correction 
  Float_t   fSlope;     // slope for gain correction
  Float_t   fAccCharge; // accumulated charge

  // Flags for enabling and disabling various features

/*
     Int_t      GetNclust() const { return fNclust; }
     Int_t      GetNhits() const  { return fNhits; }
     Float_t    GetE(int i) const { return fE[i]; }
     Float_t    GetX(int i) const { return fX[i]; }
     Float_t    GetY(int i) const { return fY[i]; }

     Int_t      GetNBlocks() { return (fNrows * fNcols);}
     Float_t    GetBlockX( Int_t i )  { if(i < fNrows*fNcols) return fBlocks[i]->GetX(); else return 0.0;}
     Float_t    GetBlockY( Int_t i )  { if(i < fNrows*fNcols) return fBlocks[i]->GetY(); else  return 0.0;}

     Float_t    GetBlockdX()  {return fdX;}
     Float_t    GetBlockdY()  {return fdY;}
     Float_t    GetBlockdZ()  {return fdZ;}


     Float_t    GetBlockA_c( Int_t i ) const { return fA_c[i]; }

     Int_t    GetNRows() {return fNrows;}
     Int_t    GetNCols() {return fNcols;}
     Int_t    BlockColRowToNumber( Int_t col, Int_t row );

  // Blocks should have a Z!!!

  SBSCalorimeterCluster* GetClust(Int_t i) { return fClusters[i]; }

  void       AddCluster(SBSCalorimeterCluster* clus);
  void       RemoveCluster(int i);
  void       AddCluster(SBSCalorimeterCluster& clus);

  void       LoadMCHitAt( Double_t x, Double_t y, Double_t E );

  protected:

  // Per-event data
  //Float_t*   fXtarg;     // [fNClust] x position (m) of clusters in target coords
  //Float_t*   fYtarg;     // [fNClust] y position (m) of clusters in target coords
  //Float_t*   fZtarg;     // [fNClust] z position (m) of clusters in target coords
  Float_t    fTRX;       // x position of track cross point
  Float_t    fTRY;       // y position of track cross point


  Double_t   fdX;
  Double_t   fdY;
  Double_t   fdZ;

  SBSShowerBlock** fBlocks; //[fNelem] Array of blocks
  SBSShowerBlock*** fBlkGrid; //[fNrows]

  //TRotation  fDetToTarg;
  //TVector3   fDetOffset;

  // Useful derived quantities for internal usage.

  Double_t tan_angle, sin_angle, cos_angle;

  void           DeleteArrays();
  */

private:
  // Simple and quick routine to init and clear most vectors
  // (of integers, floats, doubles, etc...)
  // Reset/Init 1D vector
  template<class T>
  void InitVector(std::vector<T> &vec, T val = 0, size_t n = 0) {
    vec.resize(n);
    ResetVector(vec,val);
  }

  template<class T>
  void ResetVector(std::vector<T> &vec, T val = 0, size_t n = 0) {
    if(n > 0) {
      vec.clear();
      vec.resize(n);
    }
    for(size_t i = 0; i < vec.size(); i++) {
      vec[i] = val;
    }
  }

  // Reset 2D vector
  template<class T>
  void ResetVector(std::vector<std::vector<T> > &vec, T val = 0,
      size_t nr = 0, size_t nc = 0) {
    if(nr > 0) {
      vec.clear();
      vec.resize(nr);
    }
    for(size_t i = 0; i < vec.size(); i++) {
      ResetVector(vec[i],val,nc);
    }
  }

  ClassDef(SBSCalorimeter,0)     //Generic shower detector class
};

////////////////////////////////////////////////////////////////////////////////
// Specify some default sub-classes available to the user

#endif // SBSCalorimeter_h
