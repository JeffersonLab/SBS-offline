#ifndef ROOT_SBSBBTotalShower
#define ROOT_SBSBBTotalShower

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// SBSBBTotalShower                                                          //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

//------------------------------------------------------//
//
//	Debug Definitions
//	place this section below any other head files
//
//------------------------------------------------------//
#ifdef DEBUG_LEVEL
#   undef DEBUG_LEVEL
#endif

//	DEBUG_LEVEL;	
//	=0	or not define: no debug, full speed 
//	>=1	enable debug extra warning (suggested setting)
//	>=2	above + enable debug assert
//	>=3	above + enable debug extra info
//	>=4	above + massive info (in a for or while)
#define DEBUG_LEVEL   2

#include    "DebugDef.h"
#include "THaPidDetector.h"
#include "SBSCalorimeter.h"

class SBSBBShower;
class SBSCalorimeter;

class SBSBBTotalShower : public SBSCalorimeter { //THaPidDetector {
  
 public:
  explicit SBSBBTotalShower( const char* name, const char* description = "",
		    THaApparatus* a = nullptr );
  SBSBBTotalShower( const char* name, const char* shower_name,
		    const char* preshower_name, const char* description = "",
		    THaApparatus* a = nullptr );
  virtual ~SBSBBTotalShower();

  virtual void       Clear( Option_t* opt="" );
  virtual Int_t      Decode( const THaEvData& );
  virtual Int_t      CoarseProcess( TClonesArray& tracks );
  virtual Int_t      FineProcess( TClonesArray& tracks );

  SBSBBShower* GetShower() const      { return fShower; }
  SBSBBShower* GetPreShower() const   { return fPreShower; }
  
  int PSMatchClusIdx(uint SHclusIdx){ int val = SHclusIdx<fSHclusPSclusIDmap.size() ? fSHclusPSclusIDmap[SHclusIdx] : -1; return val; }
  
  virtual EStatus    Init( const TDatime& run_time );
  virtual void       SetApparatus( THaApparatus* );
  void               LoadMCHitAt( Double_t x, Double_t y, Double_t E );
 protected:
  
  // Maximum number of clusters
  static const Int_t kMaxNClust = 16;
  
  // Subdetectors
  SBSBBShower* fShower;      // Shower subdetector
  SBSBBShower* fPreShower;   // Preshower subdetector
  

  // Parameters
  Double_t    fMaxDx;       // Maximum dx between shower and preshower centers
  Double_t    fMaxDy;       // Maximum dx between shower and preshower centers
  Double_t    fTotalSum_Threshold; //Software threshold for shower and (matched) pre-shower cluster energies

  Int_t       fPassedThreshold; //variable to indicate whether we were over threshold:
  /*
  // Per event data
  Int_t       fNclust;      // Number of clusters
  Double_t*    fE;           //[fNClust] Total shower energy
  Double_t*    fX;           //[fNClust] Total shower X
  Double_t*    fY;           //[fNClust] Total shower Y
  Int_t*      fID;          //[fNClust] ID of Presh and Shower coincidence
  */
  
  std::map<int, std::pair<int, int> > fPSSHmatchmapX;
  std::map<int, std::pair<int, int> > fPSSHmatchmapY;
  
  //key = SH cluster ID, value = PS cluster ID;
  std::vector<int> fSHclusPSclusIDmap;
  
  virtual Int_t  ReadDatabase( const TDatime& date );
  virtual Int_t  DefineVariables( EMode mode = kDefine );
  virtual Bool_t  IsPid()      { return true; }
  
 private:
  void           Setup( const char* name,  const char* desc, 
			const char* shnam, const char* psnam,
			THaApparatus* app, bool mode );
  
  ClassDef(SBSBBTotalShower,0)    //A total shower detector (shower plus preshower)
};

///////////////////////////////////////////////////////////////////////////////

#endif
