//*-- Author :    Ole Hansen (ole@jlab.org)    9-Dec-2011

/////////////////////////////////////////////////////////////////////
//
//   SBSSimDecoder
//
//   Decoder for SoLID simulation data
//
//   Interprets event buffer from input as TSBSSimEvent objects
//   (containing digitized simulation data) and unpacks them into
//   crateslot arrays for low-level decoding by detectors.
//
/////////////////////////////////////////////////////////////////////

#include "SBSSimDecoder.h"
//#include "TSBSSimDataEncoder.h"
#include "THaCrateMap.h"
#include "THaBenchmark.h"
#include "VarDef.h"
//#include "TSBSDBManager.h"
#include "THaSlotData.h"

#include "TError.h"
#include "TSystem.h"
#include "TMath.h"
#include "TDatabasePDG.h"
#include "TRandom.h"
#include "THaVarList.h"

//#include <SBSSimFadc250Module.h>// we need not to need this

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <utility>
#include <stdexcept>

using namespace std;
//using namespace Podd;

ClassImp(SBSSimDecoder) // Implements SBSSimDecoder

//EFuchey: 2016/12/10: it is necessary to declare the TSBSDBManager as a static instance here 
// (and not make it a member) because it is used by functions whic are defined as "static inline".
//static TSBSDBManager* fManager = TSBSDBManager::GetInstance();
static const Int_t kPrimaryType = 1, kPrimarySource = 0;
// Projection types must match the definitions in TreeSearch
enum EProjType { kUPlane = 0, kVPlane =1, kXPlane = 2, kYPlane = 3};

typedef vector<int>::size_type vsiz_t;

//-----------------------------------------------------------------------------
SBSSimDecoder::SBSSimDecoder() : fCheckedForEnabledDetectors(false)
{
  // Constructor
  DefineVariables();

  gSystem->Load("libEG.so");  // for TDatabasePDG
  // Get MPD encoder for GEMs
  //fEncoderMPD = dynamic_cast<TSBSSimMPDEncoder*>(
  //  TSBSSimDataEncoder::GetEncoderByName("mpd"));
}

//-----------------------------------------------------------------------------
SBSSimDecoder::~SBSSimDecoder() {

  DefineVariables( THaAnalysisObject::kDelete );
  
}

//-----------------------------------------------------------------------------
Int_t SBSSimDecoder::DefineVariables( THaAnalysisObject::EMode mode )
{
  // Define global variables for the MC quantities. Extends the base
  // class method.
  
  const char* const here = "SBSSimDecoder::DefineVariables";
  
  if( mode == THaAnalysisObject::kDefine && fIsSetup )
    return THaAnalysisObject::kOK;
  
  SimDecoder::DefineVariables( mode );
  
  cout << "Read SBSSimDecoder variables " << endl;
  
  RVarDef vars[] = {
    { 0 }
  };

  return THaAnalysisObject::
    DefineVarsFromList( vars, THaAnalysisObject::kRVarDef,
			mode, "", this, Podd::MC_PREFIX, here );
}

//-----------------------------------------------------------------------------
void SBSSimDecoder::Clear( Option_t* opt )
{
  // Clear track and plane data

  SimDecoder::Clear(opt);   // clears fMCCherHits, fMCCherClus
  
  fPMTMap.clear(); 
}

//-----------------------------------------------------------------------------
#if ANALYZER_VERSION_CODE >= ANALYZER_VERSION(1,6,0)
int SBSSimDecoder::LoadEvent(const UInt_t* evbuffer )
#else
int SBSSimDecoder::LoadEvent(const Int_t* evbuffer )
#endif
{
  // Wrapper around DoLoadEvent so we can conveniently stop the benchmark
  // counter in case of errors

  int ret = DoLoadEvent( evbuffer );

  if( fDoBench ) fBench->Stop("physics_decode");

  return ret;
}

/*
//-----------------------------------------------------------------------------
static inline
void PMTtoROC( Int_t h_chan,
	       Int_t& crate, Int_t& slot, Int_t& chan )
{
  // Convert location parameters (row, col, chan) of the given PMT
  // to hardware channel (crate,slot,chan)
  // The (crate,slot,chan) assignment must match the detmap definition in
  // the database!  See TreeSearch/dbconvert.cxx
  // In the case of GRINCH/RICH: 
  // crate = GTP; slot = VETROC; chan = PMT. (NINOs are "transparent", in a similar way to the MPDs)
  
  //div_t d = div( h_chan, fManager->GetChanPerSlot() );
  div_t d = div( h_chan, 1 );
  slot = d.quot;
  chan = d.rem;

  d = div( slot, 1 );
  crate = d.quot;
  slot  = d.rem;
}

//-----------------------------------------------------------------------------
static inline
Int_t MakeROCKey( Int_t crate, Int_t slot, Int_t chan )
{
  return chan;// +
  //fManager->GetChanPerSlot()*( slot + fManager->GetSlotPerCrate()*crate );
}

//-----------------------------------------------------------------------------
Int_t SBSSimDecoder::PMTfromROC( Int_t crate, Int_t slot, Int_t chan ) const
{
  // Return index of digitized strip correspomding to hardware channel
  // (crate,slot,chan)

  if( fPMTMap.empty() )
    return -1;

  PMTMap_t::const_iterator found = fPMTMap.find( MakeROCKey(crate,slot,chan) );
  if( found == fPMTMap.end() )
    return -1;

  return found->second;
}
*/

//-----------------------------------------------------------------------------
static inline Int_t NumberOfSetBits( UInt_t v )
{
  // Count number of bits set in 32-bit integer. From
  // http://graphics.stanford.edu/~seander/bithacks.html#CountBitsSetParallel

  v = v - ((v >> 1) & 0x55555555);
  v = (v & 0x33333333) + ((v >> 2) & 0x33333333);
  return (((v + (v >> 4)) & 0x0F0F0F0F) * 0x01010101) >> 24;
}

//-----------------------------------------------------------------------------
#if ANALYZER_VERSION_CODE >= ANALYZER_VERSION(1,6,0)
Int_t SBSSimDecoder::DoLoadEvent(const UInt_t* evbuffer )
#else
Int_t SBSSimDecoder::DoLoadEvent(const Int_t* evbuffer )
#endif
{
  // Fill crateslot structures with Monte Carlo event data in 'evbuffer'
  
  //static const char* const here = "SBSSimDecoder::LoadEvent";

#if ANALYZER_VERSION_CODE < ANALYZER_VERSION(1,6,0)
  Bool_t fNeedInit = fgNeedInit;
#endif
  assert( fMap || fNeedInit );

  // Local copy of evbuffer pointer, used in GetMCHitInfo
  buffer = evbuffer;

  // Cast the evbuffer pointer back to exactly the event type that is present
  // in the input file (in TSBSSimFile). The pointer-to-unsigned integer is
  // needed compatibility with the standard decoder.
  //const TSBSSimEvent* simEvent = reinterpret_cast<const TSBSSimEvent*>(buffer);
  
  Int_t ret = HED_OK;
  if (first_decode || fNeedInit) {
    //fMap->print();
    if( (ret = init_cmap()) != HED_OK )
      return ret;
#if ANALYZER_VERSION_CODE < ANALYZER_VERSION(1,6,0)
    if( (ret = init_slotdata(fMap)) != HED_OK)
#else
    if( (ret = init_slotdata()) != HED_OK)
#endif
      return ret;
    first_decode = false;
  }

  if( fDoBench ) fBench->Begin("clearEvent");
  Clear();
  for( int i=0; i<fNSlotClear; i++ )
    crateslot[fSlotClear[i]]->clearEvent();
  if( fDoBench ) fBench->Stop("clearEvent");

  // FIXME: needed?
  evscaler = 0;
  event_length = 0;
  
  event_type = 1;
  //event_num = simEvent->fEvtID;
  recent_event = event_num;

  // Event weight
  //fWeight = simEvent->fWeight;

  //
  if( fDoBench ) fBench->Begin("physics_decode");
  
  //Bool_t newclus;
  //Int_t crate, slot, chan,lchan;
  
  // We must check at least once which detectors are enabled
  // before we try to load up data for that detector
  if(!fCheckedForEnabledDetectors)
    CheckForEnabledDetectors();

  std::vector<std::map<Decoder::THaSlotData*, std::vector<UInt_t> > > detmaps;
  //detmaps.resize(fDetectors.size());

  /*
  // Loop through the TSBSSimEvent vector and load the data onto
  // all declared detectors.
  for(std::vector<TSBSSimEvent::DetectorData>::const_iterator it =
      simEvent->fDetectorData.begin(); it != simEvent->fDetectorData.end();
      ++it )
  {
    for(size_t d = 0; d < fDetectors.size(); d++) {
      LoadDetector(detmaps[d], fDetectors[d], (*it));
      //LoadDetector(detmaps[d], fDetNames[d], (*it), fDetIDs[d]);
    }
    // what if we were just coding the stuff above in a function ?
    // what would this function need ? name (or CPS/SPC) +detID of det, and map ???? 
    // go for it ?
    // // OK, faisons l'exercise bete de copier en adaptant pour e.g. CDet brouillon
    // if((*it).fDetID == CDET_UNIQUE_DETID && (*it).fData.size() > 0) { // 
    //   int mod =  (*it).fChannel;
    //   //This should be *general* and work for *every* subsystem
    //   chan = mod%CPS_cdet;
    //   slot = ((mod-chan)/CPS_cdet)%SPC_cdet;//+first_slot
    //   crate = (mod-slot*CPS_cdet-chan)/SPC_cdet;//+first_crate
      
    //   Decoder::THaSlotData *sldat = crateslot[idx(crate,slot)];
    //   if(sldat) { // meaning the module is available
    // 	std::vector<UInt_t> *myev = &(map[sldat]);
    // 	myev->push_back(chan);
    // 	for(size_t k = 0; k < (*it).fData.size(); k++) {
    // 	  myev->push_back((*it).fData[k]);
    // 	}
    //   }
    //   if((*it).fData[0] == 1) {
    // 	std::cerr << "M: " << mod << ", C: " << crate << ", S: " << slot
    // 		  << ", C: " << chan << ", I: " << (*it).fData[2] << std::endl;
    //   }
    // }
  }
  */
  
  // Now call LoadSlot for the different detectors
  /*
  for(size_t d = 0; d < fDetectors.size(); d++) {
    for( std::map<Decoder::THaSlotData*, std::vector<UInt_t> >::iterator it =
	   detmaps[d].begin(); it != detmaps[d].end(); ++it) {
      //unsigned short data_type = 0, chan_mult = 0;
      //unsigned int nwords = 0;
      //TSBSSimDataEncoder::DecodeHeader(it->second.front(),data_type,chan_mult,
      //    nwords);
      if(it->first->GetModule()==0) {
	if(fDebug>0) {
	  // std::cout << "No data available for detector "
	  // 	  << fDetectors[d].DetName() << std::endl;
	}
      } else {
	it->first->GetModule()->LoadSlot(it->first,
					 it->second.data(),0,it->second.size() );
      }
      //it->first->GetModule()->LoadSlot(it->first,
      //    it->second.data(),&(it->second.back()) );
    }
  }
  */
  return HED_OK;
}

/*
Int_t SBSSimDecoder::RetrieveDetMapParam(const char* detname, 
					  int& chanperslot, int& slotpercrate, 
					  int& firstcrate, int& firstslot)
{
  // chanperslot = ((TDetInfo &)fManager->GetDetInfo("hcal")).ChanPerSlot();
  // slotpercrate = ((TDetInfo &)fManager->GetDetInfo("hcal")).SlotPerCrate();
  // firstslot = ((TDetInfo &)fManager->GetDetInfo("hcal")).FirstSlot();
  // firstcrate = ((TDetInfo &)fManager->GetDetInfo("hcal")).FirstCrate();
  TDetInfo detinfo = fManager->GetDetInfo(detname);
  chanperslot = detinfo.ChanPerSlot();
  slotpercrate = detinfo.SlotPerCrate();
  firstslot = detinfo.FirstSlot();
  firstcrate = detinfo.FirstCrate();
}
*/

Int_t SBSSimDecoder::LoadDetector( std::map<Decoder::THaSlotData*,
				   std::vector<UInt_t> > &map)
//, TDetInfo &detinfo, TSBSSimEvent::DetectorData detdata)
      //const char *detname, TSBSSimEvent::DetectorData detdata, const int detid)
{
  //TDetInfo detinfo = fManager->GetDetInfo(detname);

//Int_t SBSSimDecoder::LoadDetector( std::map<Decoder::THaSlotData*, std::vector<UInt_t> > map,
//				    TSBSSimEvent::DetectorData detdata, 
//				    const int detid, 
//				    const int chanperslot, const int slotpercrate, 
//				    const int firstcrate, const int firstslot)
//{
  //int detid = detinfo.DetUniqueId();
  Int_t crate, slot;
  unsigned int nwords = 0;
  unsigned short data_type = 0, chan = 0, chan_mult = 0;
  int lchan;
  //SimEncoder::mpd_data tmp_mpd;
  
  /*
  if(detdata.fDetID == detid && detdata.fData.size() > 1) { // Data to process
    int mod =  detdata.fChannel;
    //This should be *general* and work for *every* subsystem
    // Loop over all raw data in this event
    UInt_t j = 0;
    while(j < detdata.fData.size() ) {
      // Identify the "logical" channel number for this event
      TSBSSimDataEncoder::DecodeHeader(detdata.fData[j++],data_type,chan_mult,
          nwords);
      // The channel mapping works different for GEMs vs all other detectors
      lchan = 0;
      if(detinfo.DetType() != kGEM) {
        lchan = mod + chan_mult*detinfo.NChan();
        // Get information about this logical channel from TDetInfo
        TDigChannelInfo chinfo = detinfo.FindLogicalChannelSlot(lchan);
        crate = chinfo.crate;
        slot = chinfo.slot;
        chan = chinfo.ch; // Now this is the channel in the simulated VME module
      } else {
        fEncoderMPD->DecodeMPDHeader(&(detdata.fData[j]),tmp_mpd);
        // First, decode the header information
        TDigGEMSlot gemslot = detinfo.GetGEMSlot(chan_mult);
        crate = gemslot.GetCrate();
        slot  = gemslot.GetSlot();
        tmp_mpd.mpd_id = gemslot.GetMPDId();
        tmp_mpd.gem_id = gemslot.GetGEMId();
        tmp_mpd.adc_id = gemslot.GetADCId();
        tmp_mpd.i2c = gemslot.GetI2C();
        tmp_mpd.pos = gemslot.GetPos();
        tmp_mpd.invert = gemslot.GetInvert();
        // And now, re-encode the MPD header
        fEncoderMPD->EncodeMPDHeader(tmp_mpd,&(detdata.fData[j]),chan);
      }
      Decoder::THaSlotData *sldat = 0;
      if( crate >= 0 || slot >=  0 ) {
        sldat = crateslot[idx(crate,slot)];
      }
      // Now get the corresponding THaSlotData based on crate and slot
      // and load it with the data
      // First, check that the module is defined in the cratemap, and
      // that we have at least sufficient amount of data to match that defined
      // in the header.
      if(sldat && j+nwords-1 < detdata.fData.size()) {
        std::vector<UInt_t> *myev = &(map[sldat]);
        // First, re-encode the proper channel info into the header
        myev->push_back(TSBSSimDataEncoder::EncodeHeader(
              data_type,chan,nwords));
        for(unsigned int k = 0; k < nwords; k++) {
          myev->push_back(detdata.fData[j++]);
        }
      } else {
        std::cerr << "Yikes!! No data for " << detinfo.DetName()
          << " (mod=" << mod << ") in c: "
          << crate << " s: " << slot << " c: " << chan
          << ", lchan: " << lchan << ", mult: " << chan_mult
          << " size: " << detdata.fData.size() << ", j: " << j <<", nwords: "
          << nwords << std::endl;
      }
    }
  }
  */
  return HED_OK;
}


void SBSSimDecoder::CheckForEnabledDetectors()
{
  /*
  fDetectors = fManager->GetAllDetInfo();
  if(fDebug>0) {
    for(size_t i = 0; i < fDetectors.size(); i++) {
      std::cout << "Found detector: " << fDetectors[i].DetFullName() << ", ID: "
        << fDetectors[i].DetUniqueId() << std::endl;
    }
  }
  */
  fCheckedForEnabledDetectors = true;
}
