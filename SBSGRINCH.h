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

  Int_t   FindClusters();
  Int_t   MatchClustersWithTracks( TClonesArray& tracks );

  virtual Int_t  ReadDatabase( const TDatime& date );
  virtual Int_t  DefineVariables( EMode mode = kDefine );
  
  ClassDef(SBSGRINCH,0)   //The Hall A RICH
};

#endif









