#ifndef __SBSSimDecoder_h
#define __SBSSimDecoder_h

/////////////////////////////////////////////////////////////////////
//
//   SBSSimDecoder
//
/////////////////////////////////////////////////////////////////////

#include "SimDecoder.h"
//#include "TSBSSimEvent.h"
#include "ha_compiledata.h"
#include "TTree.h"
#include "digsim_tree.h"

#include <cassert>
#include <map>
#include <stdint.h>

//class THaCrateMap;

//class TSBSSimMPDEncoder; // For decoding simulation GEMs
//class TDetInfo;

//-----------------------------------------------------------------------------
// SBS digitized simulation decoder class
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
  
  // a bit dumb, I know, but I don't know another way
  void SetTree(TTree *t);
  void AddDetector(std::string detname);
  
  void SetDetMapParam(const std::string detname, int cps, int spc, int fs, int fc);
  /*
  // Workaround for fubar THaEvData
#if ANALYZER_VERSION_CODE >= 67072  // ANALYZER_VERSION(1,6,0)
  static Int_t GetMAXSLOT() { return Decoder::MAXSLOT; }
#else
  static Int_t GetMAXSLOT() { return MAXSLOT; }
#endif
  */
  
protected:
  typedef std::map<Int_t,Int_t> PMTMap_t;

  // Event-by-event data
  PMTMap_t      fPMTMap;   //! Map ROCKey -> index of corresponding PMT
  
  // Event-by-event data
#if ANALYZER_VERSION_CODE >= 67072  // ANALYZER_VERSION(1,6,0)
  Int_t DoLoadEvent( const UInt_t* evbuffer );
#else
  Int_t DoLoadEvent( const Int_t* evbuffer );
#endif

  // retrive chanperslot, slotpercrate, etc...
  /*
  Int_t RetrieveDetMapParam(const char* detname, 
			    int& crateperslot, int& slotpercrate, 
			    int& firstcrate, int& firstslot);
  */
  //Int_t LoadDetector( std::map<Decoder::THaSlotData*, std::vector<UInt_t> > &map,
  //    const char *detname, TSBSSimEvent::DetectorData detdata, const int detid);
  //Int_t LoadDetector( std::map<Decoder::THaSlotData*, std::vector<UInt_t> > &map,
  //  TDetInfo& detinfo, TSBSSimEvent::DetectorData detdata);
  Int_t LoadDetector( std::map<Decoder::THaSlotData*, std::vector<UInt_t> > &map,
		      std::string detname, digsim_tree* tree); 
  //  TDetInfo& detinfo, TSBSSimEvent::DetectorData detdata);
  
  bool fCheckedForEnabledDetectors;
  std::vector<std::string> fDetectors;
  void CheckForEnabledDetectors();
  void CheckForDetector(const char *detname, short id);

  //TSBSSimMPDEncoder *fEncoderMPD;
  
  bool fTreeIsSet;
  digsim_tree* fTree;
  
  // again, probably dumb...
  std::map<std::string, int> fChansPerSlotDetMap;
  std::map<std::string, int> fSlotsPerCrateDetMap;
  std::map<std::string, int> fFirstSlotDetMap;
  std::map<std::string, int> fFirstCrateDetMap;

  Int_t ChanToROC( const std::string detname, Int_t h_chan, 
		   Int_t &crate, Int_t &slot, Int_t &chan ) const;
  //Int_t ChanFromROC( std::string detname, Int_t crate, Int_t slot, Int_t chan ) const;
  
  /*
  // void  PMTtoROC( Int_t s_plane, Int_t s_sector, Int_t s_proj, Int_t s_chan,
  //		    Int_t& crate, Int_t& slot, Int_t& chan ) const;
  // Int_t MakeROCKey( Int_t crate, Int_t slot, Int_t chan ) const;
  */
  ClassDef(SBSSimDecoder,1) // Decoder for simulated SoLID spectrometer data
};


#endif
