#ifndef ROOT_SBSGRINCH
#define ROOT_SBSGRINCH

//////////////////////////////////////////////////////////////////////////
//
// SBSGRINCH
//
// The Hall A RICH
//
//////////////////////////////////////////////////////////////////////////

//#include "THaPidDetector.h"
#include "SBSCherenkovDetector.h"
#include "SBSCherenkov_ClusterList.h"
#include "TBits.h"
#include "TClonesArray.h"
#include <stdint.h>
#include <map>

class THaTrack;
class THaBenchmark;

//class SBSGRINCH : public THaPidDetector {
class SBSGRINCH : public SBSCherenkovDetector {
  
public:

  explicit SBSGRINCH( const char* name, const char* description="",
	   THaApparatus* apparatus=nullptr );
  virtual ~SBSGRINCH();
  
  virtual void         Clear( Option_t* opt="" );
  virtual Int_t        Decode( const THaEvData& );
  virtual Int_t        CoarseProcess( TClonesArray& tracks );
  virtual Int_t        FineProcess( TClonesArray& tracks );

protected:

  //Double_t fTrackX;      // x pos of Golden Track in RICH plane
  //Double_t fTrackY;      // y pos of Golden Track in RICH plane

  Double_t fMaxSep; // Max separation between PMT and another PMT to count as "neighbors"
  Double_t fMaxSep2; //square of fMaxSep

  Int_t fNmirror; //Number of GRINCH mirrors (define track match cuts separately for each mirror)

  Int_t fOrderTrackMatchY; // Default to 3. For now we implement GRINCH dy = pol3( track phi ); later we may get fancier:
  
  //P slope is obsolete, but for now I keep it. AJRP 10/23/23
  Double_t fTrackMatchPslope; //slope of xtrack - xGRINCH vs 1/p, default 0.1715
  //make mirror-dependent track match cuts:
  std::vector<Double_t> fTrackMatchXslope;
  std::vector<Double_t> fTrackMatchX0;
  std::vector<Double_t> fTrackMatchXsigma; 
  std::vector<Double_t> fTrackMatchYslope;
  std::vector<Double_t> fTrackMatchY0;
  std::vector<Double_t> fTrackMatchYsigma;
  std::vector<Double_t> fTrackMatchXmin; //minimum track x projection to consider for this mirror
  std::vector<Double_t> fTrackMatchXmax; //maximum track x projection to consider for this mirror

  
  
  Double_t fTrackMatchNsigmaCut; //use common cut width for each mirror

  virtual Int_t   FindClusters();
  virtual Int_t   MatchClustersWithTracks( TClonesArray& tracks );
  virtual Int_t   SelectBestCluster(Int_t nmatch=0);

  //  Int_t fBestClusterIndex; //the biggest cluster with a track match, if any matched clusters are found. Perhaps this should go with SBSCherenkov_Detector rather than SBSGRINCH

  virtual Int_t  ReadDatabase( const TDatime& date );
  virtual Int_t  DefineVariables( EMode mode = kDefine );
  
  ClassDef(SBSGRINCH,0)   //The Hall A RICH
};

#endif









