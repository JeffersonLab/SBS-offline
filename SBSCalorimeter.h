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
#include "SBSCalorimeterCluster.h"
#include "TRotation.h"
#include "TVector3.h"
#include "THaDetMap.h"

struct SBSBlockSet {
  Float_t e;
  Int_t row;
  Int_t col;
  Int_t id;
  Float_t TDCTime;
  Float_t ADCTime;
};

struct SBSCalBlocks {
  std::vector<Float_t> e;   //< []
  std::vector<Float_t> TDCTime;   //< [] 
  std::vector<Float_t> ADCTime;   //< [] 
  std::vector<Int_t> row; //< []
  std::vector<Int_t> col; //< []
  std::vector<Float_t> x; //< []
  std::vector<Float_t> y; //< []
  std::vector<Int_t>   id;      // []
  void clear() {
    e.clear();
    x.clear();
    y.clear();
    id.clear();
    row.clear();
    col.clear();
    TDCTime.clear();
    ADCTime.clear();
  }
};

struct SBSCalorimeterOutput {
  std::vector<Float_t> e;   //< []
  std::vector<Float_t> e_c;   //< []
  std::vector<Float_t> x;   //< []
  std::vector<Float_t> y;   //< []
  std::vector<Float_t> row; //< []
  std::vector<Float_t> col; //< []
  std::vector<Int_t>   n;   // Number of elements
  std::vector<Float_t> blk_e;   // block energy of max energy block
  std::vector<Float_t> blk_e_c; // block corrected energy of max energy block
  std::vector<Int_t>   id;      // block id of max energy block
};
   typedef std::vector<SBSBlockSet> SBSBlockSetList;

class SBSCalorimeter : public SBSGenericDetector {

public:
  SBSCalorimeter( const char* name, const char* description = "",
      THaApparatus* a = NULL);
  virtual ~SBSCalorimeter();

  virtual void ClearEvent();
  virtual void ClearOutputVariables();
  virtual Int_t MakeGoodBlocks();
  virtual Int_t FindClusters();

  // Get information from the main cluster
  Float_t GetE();             //< Main cluster energy
  Float_t GetECorrected();    //< Main cluster corrected energy
  Float_t GetX();             //< Main cluster energy average x
  Float_t GetY();             //< Main cluster energy average y
  Float_t GetEBlk();          //< Main cluster energy of max block in cluster
  Float_t GetEBlkCorrected(); //< Main cluster corrected energy of max block in cluster
  Int_t GetNblk();            //< Number of blocks in main cluster
  Int_t GetBlkID();           //< ID/block number of max energy block in cluster
  Int_t GetRow();             //< Main cluster row of max block
  Int_t GetCol();             //< Main cluster col of max block

  Int_t    GetNRows() {return fNrows;}
  Int_t    GetNCols(Int_t r = 0);

  // Standard apparatus re-implemented functions
  virtual Int_t      FineProcess(TClonesArray& tracks);

  Int_t GetNclust() const { return fNclus; }

  void SetDataOutputLevel(int var) { fDataOutputLevel = var; }


protected:

  virtual Int_t  ReadDatabase( const TDatime& date );
  virtual Int_t  DefineVariables( EMode mode = kDefine );
  void ClearCaloOutput(SBSCalorimeterOutput &var);

  Int_t  fNclus;        ///< Size of clusters in the output tree

  // Configuration
  Int_t  fNclublk;      ///< Max number of blocks composing a cluster
  Int_t  fNclubr;       ///< Max number of row-blocks composing a cluster
  Int_t  fNclubc;       ///< Max number of col-blocks composing a cluster

  // Mapping (see also fDetMap)

  SBSCalBlocks fGoodBlocks; // < Good block structure for tree output
  SBSBlockSetList fBlockSet; // < Vector of structure of GoodBlock info to use in FindCluster 
  typedef SBSBlockSetList::iterator fBlockSetIterator; // 
  // Output variables
  SBSCalorimeterOutput fMainclus;    //< Main cluster output
  SBSCalorimeterOutput fMainclusblk; //< Detailed info on blocks in main cluster
  SBSCalorimeterOutput fOutclus;     //< Output for all clusters

  // Clusters for this event
  std::vector<SBSCalorimeterCluster*> fClusters; // Cluster
  Float_t    fEmin;         //< Minimum energy for a cluster center
  Float_t    fXmax_dis;         //< maximum X distance from a cluster center
  Float_t    fYmax_dis;         //< maximum Y distance from a cluster center
  Float_t    fRmax_dis;         //< maximum radius from a cluster center
  Int_t fMaxNclus;          //< Maximum number of clusters to store
  Bool_t fCoarseProcessed;  //< Was CourseProcessed already called?
  Bool_t fFineProcessed;    //< Was fine processed already called

  //Gain correction 
  Float_t   fConst;     // const from gain correction 
  Float_t   fSlope;     // slope for gain correction
  Float_t   fAccCharge; // accumulated charge

  // ROOTFile output level
  Int_t fDataOutputLevel;   //  0: default only main cluster info, 1: include also blocks in main cluster, 2: all clusters,  3: all blocks regardless if they are in a cluster or not

  Float_t GetVVal(std::vector<Float_t> &v, UInt_t i = 0 );
  Int_t GetVVal(std::vector<Int_t> &v, UInt_t i = 0 );

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

inline Float_t SBSCalorimeter::GetVVal(std::vector<Float_t> &v, UInt_t i)
{
  if (v.size() > i) {
    return v[i];
  }
  return 0.0;
}

inline Int_t SBSCalorimeter::GetVVal(std::vector<Int_t> &v, UInt_t i)
{
  if (v.size() > i) {
    return v[i];
  }
  return 0.0;
}

inline Float_t SBSCalorimeter::GetE() {
  return GetVVal(fMainclus.e);
}

inline Float_t SBSCalorimeter::GetECorrected() {
  return GetVVal(fMainclus.e_c);
}

inline Float_t SBSCalorimeter::GetEBlk() {
  return GetVVal(fMainclus.blk_e);
}

inline Float_t SBSCalorimeter::GetEBlkCorrected() {
  return GetVVal(fMainclus.blk_e_c);
}

inline Float_t SBSCalorimeter::GetX() {
  return GetVVal(fMainclus.x);
}

inline Float_t SBSCalorimeter::GetY() {
  return GetVVal(fMainclus.y);
}

inline Int_t SBSCalorimeter::GetRow() {
  return GetVVal(fMainclus.row);
}

inline Int_t SBSCalorimeter::GetCol() {
  return GetVVal(fMainclus.col);
}

inline Int_t SBSCalorimeter::GetNblk() {
  return GetVVal(fMainclus.n);
}

inline Int_t SBSCalorimeter::GetBlkID() {
  return GetVVal(fMainclus.id);
}

inline Int_t SBSCalorimeter::GetNCols(Int_t r) {
  return GetVVal(fNcols,r);
}


////////////////////////////////////////////////////////////////////////////////
// Specify some default sub-classes available to the user

#endif // SBSCalorimeter_h
