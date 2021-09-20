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
#include "SBSSimFile.h"//needed for SBSSimEvent
#include "TH1D.h"
//#include "digsim_tree.h"
#include "THaApparatus.h"

#include <cassert>
#include <map>
#include <stdint.h>

class THaDetMap;
class SBSSimSADCEncoder; // For decoding simulation GEMs
//class TDetInfo;

//-----------------------------------------------------------------------------
// SBS digitized simulation decoder class
class SBSSimDecoder : public Podd::SimDecoder {
 public:
  //constructor may be inputed a data file to input some of the paramaters used by SimDecoder
  //NB: if the second file path does not select a valid file, default parameters will be used.
  // MANDATORY
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
  
  // Workaround for fubar THaEvData
#if ANALYZER_VERSION_CODE >= 67072  // ANALYZER_VERSION(1,6,0)
  static Int_t GetMAXSLOT() { return Decoder::MAXSLOT; }
#else
  static Int_t GetMAXSLOT() { return MAXSLOT; }
#endif

  //Needs to be public
  struct detchaninfo{
    uint crate;
    uint slot;
    uint chan;
  detchaninfo(int cr = 0, int sl = 0, int ch = 0) ://, int apv = 0) : 
    crate(cr), slot(sl), chan(ch)//, apvnum(apv)
    {}
    virtual ~detchaninfo(){};
  };

  //strips have their own struc this way we can add as many parameters as we see fit
  struct gemstripinfo{
    uint crate;
    uint slot;
    uint apvnum;
    //uint chan_lo;
    //uint chan_hi;
    //uint axis;
  gemstripinfo(int cr = 0, int sl = 0, int apv = 0) :
	       //, int ch_lo = 0, int ch_hi = 0, int ax) : 
    crate(cr), slot(sl), apvnum(apv)
      //, chan_lo(ch_lo), chan_hi(ch_hi), axis(ax)
    {}
    virtual ~gemstripinfo(){};
  };
  
  //Utilities
  // a bit dumb, I know, but I don't know another way
  //void SetTree(TTree *t);
  //Setup all detectors for an apparatus
  //void SetDetMapParam(const std::string detname, int cps, int spc, int fs, int fc);
  
protected:
  // MANDATORY
  // Event-by-event data
#if ANALYZER_VERSION_CODE >= 67072  // ANALYZER_VERSION(1,6,0)
  Int_t DoLoadEvent( const UInt_t* evbuffer );
#else
  Int_t DoLoadEvent( const Int_t* evbuffer );
#endif

  //Utilities
  //typedef std::map<Int_t,Int_t> PMTMap_t;

  // Event-by-event data
  //PMTMap_t      fPMTMap;   //! Map ROCKey -> index of corresponding PMT
  
  // retrive chanperslot, slotpercrate, etc...
  /*
  Int_t RetrieveDetMapParam(const char* detname, 
			    int& crateperslot, int& slotpercrate, 
			    int& firstcrate, int& firstslot);
  */
  void SetDetectors();
  Int_t AddDetector(std::string detname, TDatime date);
  Int_t ReadDetectorDB(std::string detname, TDatime date);
  Int_t LoadDetector( std::map<Decoder::THaSlotData*, std::vector<UInt_t> > &map,
		      const std::string detname, 
		      const SBSSimEvent* simev); 
  
  void CheckForEnabledDetectors();
  //void CheckForDetector(const char *detname, short id);

  SBSSimSADCEncoder *fDecoderMPD;//FIXME: a bit of a kludge... 
  
  bool fCheckedForEnabledDetectors;
  std::vector<std::string> fDetectors;
  
  //bool fTreeIsSet;
  //digsim_tree* fTree;
  
  //std::map<std::string, UInt_t> fNChanDet;
  //std::map<std::string, UInt_t> fChanMapStartDet;
  //std::map<std::string, std::map<UInt_t, detchaninfo> > fInvDetMap;
  std::map<std::string, std::vector<detchaninfo> > fInvDetMap;
  std::map<std::string, std::vector<std::vector<gemstripinfo>> > fInvGEMDetMap;
  //std::map<std::string, std::vector<std::vector<gemstripinfo>> > fInvGEMDetMap;
  //std::map<std::string, std::vector<detchaninfo> > fInvDetMap_secondary;
  //std::map<std::string, std::vector< std::vector<UShort_t> > > fChanMapDet;
  
  /*
  // again, probably dumb...
  std::map<std::string, uint> fChansPerSlotDetMap;
  std::map<std::string, uint> fSlotsPerCrateDetMap;
  std::map<std::string, uint> fFirstSlotDetMap;
  std::map<std::string, uint> fFirstCrateDetMap;
  */
  
  void ChanToROC( const std::string detname, Int_t h_chan, 
		  Int_t &crate, Int_t &slot, UShort_t &chan ) const;
  
  int APVnum( const std::string detname, Int_t mod, Int_t h_chan, 
	      Int_t &crate, Int_t &slot, UShort_t &chan ) const;
  
  // TODO: function(s) that load(s) the MC track hit
  Double_t fPx;
  Double_t fPy;
  Double_t fPz;
  Double_t fVx;
  Double_t fVy;
  Double_t fVz;
  Int_t fNBBtracks;
  std::vector<Int_t> fBBtrack_Nhits;
  std::vector<Int_t> fBBtrack_TID;
  std::vector<Int_t> fBBtrack_PID;
  std::vector<Int_t> fBBtrack_MID;
  std::vector<Double_t> fBBtrack_P;
  std::vector<Double_t> fBBtrack_X;
  std::vector<Double_t> fBBtrack_Y;
  std::vector<Double_t> fBBtrack_dX;
  std::vector<Double_t> fBBtrack_dY;
  Int_t fNBBGEMhits; 
  std::vector<Int_t> fBBGEMhit_plane;
  std::vector<Int_t> fBBGEMhit_TID;
  std::vector<Int_t> fBBGEMhit_PID;
  std::vector<Int_t> fBBGEMhit_MID;
  std::vector<Double_t> fBBGEMhit_edep;
  std::vector<Double_t> fBBGEMhit_x;
  std::vector<Double_t> fBBGEMhit_y;
  
  //TH1D* h1_sizeHCal;
  //TH1D* h1_sizeGEMs;
  
  //Int_t ChanFromROC( std::string detname, Int_t crate, Int_t slot, Int_t chan ) const;
  /*
  // void  PMTtoROC( Int_t s_plane, Int_t s_sector, Int_t s_proj, Int_t s_chan,
  //		    Int_t& crate, Int_t& slot, Int_t& chan ) const;
  // Int_t MakeROCKey( Int_t crate, Int_t slot, Int_t chan ) const;
  */
  ClassDef(SBSSimDecoder,1) // Decoder for simulated SoLID spectrometer data
};


#endif
