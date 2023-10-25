#ifndef ROOT_SBSCherenkovDetector
#define ROOT_SBSCherenkovDetector

//////////////////////////////////////////////////////////////////////////
//
// SBSCherenkovDetector
//
// The Hall A RICH
//
//////////////////////////////////////////////////////////////////////////

//#include "THaPidDetector.h"
#include "SBSGenericDetector.h"
#include "SBSCherenkov_ClusterList.h"
#include "TBits.h"
#include "TClonesArray.h"
#include <cstdint>
#include <map>

class THaTrack;
class THaBenchmark;

class SBSCherenkovDetector : public SBSGenericDetector {
  
public:

  explicit SBSCherenkovDetector( const char* name, const char* description="",
	   THaApparatus* apparatus=nullptr );
  virtual ~SBSCherenkovDetector();
  
  virtual void         Clear( Option_t* opt="" );
  virtual Int_t        Decode( const THaEvData& );
  virtual Int_t        CoarseProcess( TClonesArray& tracks );
  virtual Int_t        FineProcess( TClonesArray& tracks );

  SBSCherenkov_Hit*         GetHit(Int_t i) const  
  { return (SBSCherenkov_Hit*)fHits->At(i); } 
  SBSCherenkov_Cluster*     GetCluster(Int_t i) const 
  { return (SBSCherenkov_Cluster*)fClusters->At(i); } 
  
  SBSCherenkov_Cluster*     GetBestCluster() const 
  { return GetCluster(fBestClusterIndex); }

  
  Int_t                GetNumHits() const 
    { return fHits->GetLast()+1; }
  Int_t                GetNumClusters() const 
    { return fClusters->GetLast()+1; }

  //Int_t                GetMaxNumHits() const { return fMaxNumHits; }
  //void                 SetMaxNumHits( Int_t MaxNumHit ) 
  //{ fMaxNumHits=MaxNumHit; }
  void                 EnableClusterResolving( Bool_t flag = kTRUE )
  { fDoResolve = flag; }
  void                 EnableBenchmarks( Bool_t b = kTRUE )
  { fDoBench = b; }
  void                 PrintBenchmarks() const;
  
protected:

  TClonesArray*     fHits;          // Array of hits for each event
  TClonesArray*     fClusters;      // Clusters of hits

  SBSCherenkov_Cluster fBestCluster; //Best cluster object

  Int_t         fBestClusterIndex; //index in the array of "best" cluster

  Int_t         fNtrackMatch;

  Bool_t  fDoResolve;    // true = resolve overlapping clusters
  Bool_t  fDoTimeFilter; // true = filter the hits in each cluster with timing
  
  //Double_t fTrackX;      // x pos of Golden Track in RICH plane
  //Double_t fTrackY;      // y pos of Golden Track in RICH plane

  Bool_t         fDoBench;         //Collect detailed timing statistics
  THaBenchmark*  fBench;           //Counters for timing statistics
  
  void    DeleteClusters();
  virtual Int_t   FindClusters(){return 0;};
  virtual Int_t   MatchClustersWithTracks( TClonesArray& tracks ){return 0;};
  virtual Int_t   SelectBestCluster(Int_t nmatch = 0){return 0;};
  //Int_t   CleanClustersWithTime();

  // We will use one tmin and tmax value for all channels in the detector, the individual channel offsets will be used to align the 
  // good signal peaks at a common central value
  Double_t fHit_tmin;
  Double_t fHit_tmax;

  std::vector<Double_t> fAmpToTCoeff;
  
  virtual Int_t  ReadDatabase( const TDatime& date );
  virtual Int_t  DefineVariables( EMode mode = kDefine );

  //private:
  
  bool fMCdata;// easy way to enable/disable the use of MC data.
  
  ClassDef(SBSCherenkovDetector,0)   //The Hall A RICH
};

#endif









