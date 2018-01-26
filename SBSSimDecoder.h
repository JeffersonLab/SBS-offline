#ifndef __SBSSimDecoder_h
#define __SBSSimDecoder_h

/////////////////////////////////////////////////////////////////////
//
//   SBSSimDecoder
//
/////////////////////////////////////////////////////////////////////

#include "SimDecoder.h"
//#include "SBSSimEvent.h"
#include "ha_compiledata.h"
#include <cassert>
#include <map>
//#include "SBSDBManager.h"
#include <stdint.h>

class THaCrateMap;

//-----------------------------------------------------------------------------
// Helper classes for making decoded event data available via global variables

class SBSSimPMTHit : public TObject {
public:
  SBSSimPMTHit() {}
  //SBSSimPMTHit(const SBSSimEvent::PMTHit& hit );

  virtual void Print( const Option_t* opt="" ) const;

  Short_t  fID;          // Hit number
  Int_t    fSource;      // MC data set source (0 = signal, >0 background)
  Char_t   fType;        // GEANT particle type (1 = primary)
  Int_t    fMCtrackID;   // GEANT track ID (if any)
  Int_t    fMCtrackPID;  // GEANT particle ID (if any)
  Short_t  fOrigVolFlag; // 
  Float_t  fXPMT;        // X coordinate of the PMT in transport coordinates
  Float_t  fYPMT;        // Y coordinate of the PMT in transport coordinates
  Float_t  fNpe;         // Number of photoelectrons
  Double_t fTime;        // Arrival time at electronics
  Double_t fTDCtime1;    // TDC rising time
  Double_t fTDCtime2;    // TDC falling time
  // Digitization results for this hit
  Short_t  fDetID;       // Detector ID
  Short_t  fChannel;     // Channel number
  Short_t  fPMTrow;      // Row number: cross reference to Channel number
  Short_t  fPMTcol;      // Column number: cross reference to Channel number
  Char_t   fVETROCID;    // VETROC ID
  uint32_t fTDC1;        // VETROC TDC "word" for rising time
  uint32_t fTDC2;        // VETROC TDC "word" for falling time
  
  ClassDef(SBSSimPMTHit,1) // A Monte Carlo hit at a GEM tracking chamber
};

class SBSSimCherCluster : public TObject {
public:
  SBSSimCherCluster() {}
  
  virtual void Print( const Option_t* opt="" ) const;
  
  //TList*     fHitList;   //List of hits belonging to this cluster
  std::vector<int> fHitList; //List of hits belonging to this cluster
  
  Int_t      fSize;      // Number of Hits in the cluster 
  
  Float_t    fXcenter;   // X mean of all hits in the list
  Float_t    fYcenter;   // Y mean of all hits in the list
  Float_t    fXcenter_w; // Weighted X mean : (Sum of x*adc)/(sum adc) of all hits in the list
  Float_t    fYcenter_w; // Weighted Y mean : (Sum of y*adc)/(sum adc) of all hits in the list
  Int_t      fNpe;       // Total number of photoelectrons 
  
  Float_t    fMeanRisingTime;  // Mean rising time of all hits in the list
  Float_t    fMeanFallingTime; // Mean falling time of all hits in the list
  Float_t    fRisingTimeRMS;   // Rising time RMS of all hits in the list
  Float_t    fFallingTimeRMS;  // Falling time RMS of all hits in the list
  
  Int_t      fMCtrackID;      // MC track ID
  Int_t      fMCtrackPID;      // MC track PID, if applicable

  ClassDef(SBSSimCherCluster,1) // A Monte Carlo hit at a GEM tracking chamber
};

//-----------------------------------------------------------------------------
// SoLID simulation decoder class
class SBSSimDecoder : public Podd::SimDecoder {
 public:
  //constructor may be inputed a data file to input some of the paramaters used by SimDecoder
  //NB: if the second file path does not select a valid file, default parameters will be used.
  SBSSimDecoder();
  virtual ~SBSSimDecoder();
  
#if ANALYZER_VERSION_CODE >= 67072 // ANALYZER_VERSION(1,6,0)
  virtual Int_t LoadEvent( const UInt_t* evbuffer );
#else
  virtual Int_t LoadEvent( const Int_t* evbuffer );
#endif
  virtual void  Clear( Option_t* opt="" );
  virtual Int_t DefineVariables( THaAnalysisObject::EMode mode =
				 THaAnalysisObject::kDefine );
  //virtual Podd::MCHitInfo GetMCHitInfo( Int_t crate, Int_t slot, Int_t chan ) const;
  
  Int_t    GetNPMThits()  const {
    return (fMCCherHits) ? fMCCherHits->GetLast()+1 : 0;
  }
  Int_t    GetNPMTclus()  const {
    return (fMCCherClus) ? fMCCherClus->GetLast()+1 : 0;
  }
  
  SBSSimPMTHit* GetPMTHit( Int_t i ) const {
    TObject* obj = fMCCherHits->UncheckedAt(i);
    assert( dynamic_cast<SBSSimPMTHit*>(obj) );
    return static_cast<SBSSimPMTHit*>(obj);
  }
  
  SBSSimCherCluster* GetPMTclus( Int_t i ) const {
    TObject* obj = fMCCherClus->UncheckedAt(i);
    assert( dynamic_cast<SBSSimCherCluster*>(obj) );
    return static_cast<SBSSimCherCluster*>(obj);
  }
  
  // Workaround for fubar THaEvData
#if ANALYZER_VERSION_CODE >= 67072  // ANALYZER_VERSION(1,6,0)
  static Int_t GetMAXSLOT() { return Decoder::MAXSLOT; }
#else
  static Int_t GetMAXSLOT() { return MAXSLOT; }
#endif

protected:
  typedef std::map<Int_t,Int_t> PMTMap_t;

  TClonesArray* fMCCherHits;
  TClonesArray* fMCCherClus;
  // Event-by-event data
  PMTMap_t      fPMTMap;   //! Map ROCKey -> index of corresponding PMT

  // Event-by-event data
#if ANALYZER_VERSION_CODE >= 67072  // ANALYZER_VERSION(1,6,0)
  Int_t DoLoadEvent( const UInt_t* evbuffer );
#else
  Int_t DoLoadEvent( const Int_t* evbuffer );
#endif

  // void  PMTtoROC( Int_t s_plane, Int_t s_sector, Int_t s_proj, Int_t s_chan,
  //		    Int_t& crate, Int_t& slot, Int_t& chan ) const;
  Int_t PMTfromROC( Int_t crate, Int_t slot, Int_t chan ) const;
  // Int_t MakeROCKey( Int_t crate, Int_t slot, Int_t chan ) const;
  
  ClassDef(SBSSimDecoder,0) // Decoder for simulated SoLID spectrometer data
};


#endif
