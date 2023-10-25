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
#include "Helper.h"

struct SBSBlockSet {
  Double_t e;
  Double_t x;
  Double_t y;
  Int_t row;
  Int_t col;
  Int_t id;
  Double_t TDCTime;
  Double_t ADCTime;
  Bool_t InCluster;
};

struct SBSCalBlocks {
  std::vector<Double_t> e;   //< []
  std::vector<Double_t> TDCTime;   //< [] 
  std::vector<Double_t> ADCTime;   //< [] 
  std::vector<Int_t> row; //< []
  std::vector<Int_t> col; //< []
  std::vector<Double_t> x; //< []
  std::vector<Double_t> y; //< []
  std::vector<Int_t>   id;      // []
  std::vector<Int_t>   cid;
  void clear() {
    e.clear();
    x.clear();
    y.clear();
    id.clear();
    cid.clear();
    row.clear();
    col.clear();
    TDCTime.clear();
    ADCTime.clear();
  }
};

struct SBSCalorimeterOutput {
  std::vector<Double_t> e;   //< []
  std::vector<Double_t> again;   //< []
  std::vector<Double_t> atime;   //< []
  std::vector<Double_t> tdctime;   //< []
  //std::vector<Double_t> e_c;   //< []
  std::vector<Double_t> x;   //< []
  std::vector<Double_t> y;   //< []
  std::vector<Double_t> row; //< []
  std::vector<Double_t> col; //< []
  std::vector<Int_t>   n;   // Number of elements
  std::vector<Double_t> blk_e;   // block energy of max energy block
  //std::vector<Double_t> blk_e_c; // block corrected energy of max energy block
  std::vector<Int_t>   id;      // block id of max energy block
};
   
typedef std::vector<SBSBlockSet> SBSBlockSetList;

class SBSCalorimeter : public SBSGenericDetector {

public:
  explicit SBSCalorimeter( const char* name, const char* description = "",
      THaApparatus* a = nullptr);
  virtual ~SBSCalorimeter();

  virtual void Clear( Option_t* opt="" );
  virtual Int_t MakeGoodBlocks();
  virtual Int_t FindClusters();

  // Get information from the main cluster
  Double_t GetE();             //< Main cluster energy
  Double_t GetAgain();         //< Pedestal subtracted ADC integral (pC)
  Double_t GetAtime();         //< Main cluster ADC time of max block
  Double_t GetTDCtime();         //< Main cluster ADC time of max block
  /* Double_t GetECorrected();    //< Main cluster corrected energy */
  Double_t GetX();             //< Main cluster energy average x
  Double_t GetY();             //< Main cluster energy average y
  Double_t GetEBlk();          //< Main cluster energy of max block in cluster
  /* Double_t GetEBlkCorrected(); //< Main cluster corrected energy of max block in cluster */
  Int_t GetNblk();            //< Number of blocks in main cluster
  Int_t GetBlkID();           //< ID/block number of max energy block in cluster
  Int_t GetRow();             //< Main cluster row of max block
  Int_t GetCol();             //< Main cluster col of max block

  Int_t    GetNRows() {return fNrows;}
  Int_t    GetNCols(Int_t r = 0);

  // Standard apparatus re-implemented functions
  virtual Int_t      FineProcess(TClonesArray& tracks);
  std::vector<SBSBlockSet>& GetBlockSet() {return fBlockSet;}
  std::vector<SBSCalorimeterCluster*>& GetClusters() {return fClusters;}

  Int_t GetNclust() const { return fNclus; }

  void SetDataOutputLevel(int var) { fDataOutputLevel = var; }

  Double_t GetTmax() const { return fTmax; }

protected:

  virtual Int_t  ReadDatabase( const TDatime& date );
  virtual Int_t  DefineVariables( EMode mode = kDefine );
  void ClearCaloOutput(SBSCalorimeterOutput &var);

  Int_t  fNclus;        ///< Size of clusters in the output tree

  // Configuration
  Int_t  fNclublk;      ///< Max number of blocks composing a cluster
  Int_t  fNclubr;       ///< Max number of row-blocks composing a cluster
  Int_t  fNclubc;       ///< Max number of col-blocks composing a cluster

  Int_t fBestClusterIndex; //Index of best cluster in the array.
  
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
  
  Double_t    fTmax;            //< Maximum time difference for cluster block
  Double_t    fEmin;         //< Minimum energy for a cluster block
  Double_t fEmin_clusSeed; //< Minimum energy to be the seed of a cluster
  Double_t fEmin_clusTotal; //< Total energy threshold of a cluster

  Double_t    fXmax_dis;         //< maximum X distance from a cluster center
  Double_t    fYmax_dis;         //< maximum Y distance from a cluster center
  Double_t    fRmax_dis;         //< maximum radius from a cluster center
  Int_t fMaxNclus;          //< Maximum number of clusters to store
  Bool_t fCoarseProcessed;  //< Was CourseProcessed already called?
  Bool_t fFineProcessed;    //< Was fine processed already called

  //Gain correction 
  Double_t   fConst;     // const from gain correction 
  Double_t   fSlope;     // slope for gain correction
  Double_t   fAccCharge; // accumulated charge

  // ROOTFile output level
  Int_t fDataOutputLevel;   //  0: default only main cluster info, 1: include also blocks in main cluster, 2: all clusters,  3: all blocks regardless if they are in a cluster or not

  Double_t GetVVal(std::vector<Double_t> &v, UInt_t i = 0 );
  Int_t GetVVal(std::vector<Int_t> &v, UInt_t i = 0 );

private:
  void ClearOutputVariables();

  ClassDef(SBSCalorimeter,0)     //Generic shower detector class
};

inline Double_t SBSCalorimeter::GetVVal(std::vector<Double_t> &v, UInt_t i)
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

inline Double_t SBSCalorimeter::GetE() {
  return GetVVal(fMainclus.e);
}

inline Double_t SBSCalorimeter::GetAgain() {
  return GetVVal(fMainclus.again);
}

inline Double_t SBSCalorimeter::GetAtime() {
  return GetVVal(fMainclus.atime);
}

inline Double_t SBSCalorimeter::GetTDCtime() {
  return GetVVal(fMainclus.tdctime);
}

/* inline Double_t SBSCalorimeter::GetECorrected() { */
/*   return GetVVal(fMainclus.e_c); */
/* } */

inline Double_t SBSCalorimeter::GetEBlk() {
  return GetVVal(fMainclus.blk_e);
}

/* inline Double_t SBSCalorimeter::GetEBlkCorrected() { */
/*   return GetVVal(fMainclus.blk_e_c); */
/* } */

inline Double_t SBSCalorimeter::GetX() {
  return GetVVal(fMainclus.x);
}

inline Double_t SBSCalorimeter::GetY() {
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
