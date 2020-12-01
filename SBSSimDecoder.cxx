//*-- Author :    Ole Hansen (ole@jlab.org)    9-Dec-2011

/////////////////////////////////////////////////////////////////////
//
//   SBSSimDecoder
//
//   Decoder for SoLID simulation data
//
//   Interprets event buffer from input as SBSSimEvent objects
//   (containing digitized simulation data) and unpacks them into
//   crateslot arrays for low-level decoding by detectors.
//
/////////////////////////////////////////////////////////////////////

#include "SBSSimDecoder.h"
#include "SBSSimDataEncoder.h"
#include "THaCrateMap.h"
#include "THaBenchmark.h"
#include "VarDef.h"
#include "THaSlotData.h"

#include "TError.h"
#include "TSystem.h"
#include "TMath.h"
#include "TDatabasePDG.h"
#include "TRandom.h"
#include "THaVarList.h"
#include "THaDetMap.h"
#include "THaDetector.h"
#include "SBSBBShower.h"
#include "SBSBBTotalShower.h"

//#include <SBSSimFadc250Module.h>// we need not to need this
#include "TList.h"
#include "TObject.h"

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <utility>
#include <stdexcept>

using namespace std;
//using namespace Podd;

class THaAnalysisObject;

ClassImp(SBSSimDecoder) // Implements SBSSimDecoder


//static const Int_t kPrimaryType = 1, kPrimarySource = 0;
// Projection types must match the definitions in SBS-offline
//enum EProjType { kUPlane = 0, kVPlane =1, kXPlane = 2, kYPlane = 3};
//typedef vector<int>::size_type vsiz_t;

//-----------------------------------------------------------------------------
SBSSimDecoder::SBSSimDecoder()// : fCheckedForEnabledDetectors(false), fTreeIsSet(false)
{
  // Constructor
  DefineVariables();
  fDetectors.clear();
  //fTree = 0;
  // Load detectors: rely on gHaApps (please tell me it works!!!)
  cout << " Calling SBSSimDecoder! "<< endl;
  cout << " Make sure you have already declared your apparatuses and detectors, and added these to gHaApps" << endl;
  SetDetectors();
  
  gSystem->Load("libEG.so");  // for TDatabasePDG
  // Get MPD encoder for GEMs
  fEncoderMPD = dynamic_cast<SBSSimMPDEncoder*>
    (SBSSimDataEncoder::GetEncoderByName("mpd"));
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

//-----------------------------------------------------------------------------
#if ANALYZER_VERSION_CODE >= ANALYZER_VERSION(1,6,0)
Int_t SBSSimDecoder::DoLoadEvent(const UInt_t* evbuffer )
#else
Int_t SBSSimDecoder::DoLoadEvent(const Int_t* evbuffer )
#endif
{
  // Fill crateslot structures with Monte Carlo event data in 'evbuffer'
  static const char* const here = "SBSSimDecoder::LoadEvent";

#if ANALYZER_VERSION_CODE < ANALYZER_VERSION(1,6,0)
  Bool_t fNeedInit = fgNeedInit;
#endif
  assert( fMap || fNeedInit );

  // Local copy of evbuffer pointer, used in GetMCHitInfo
  buffer = evbuffer;

  // if(!fTreeIsSet){
  //   std::cerr << "SBSSimDecoder Tree not initialized correctly - exiting" << std::endl;
  //   return HED_FATAL;
  // }
  //fTree->GetEntry(GetEvNum());
  // Cast the evbuffer pointer back to exactly the event type that is present
  // in the input file (in SBSSimFile). The pointer-to-unsigned integer is
  // needed compatibility with the standard decoder.
  if(fDebug>2)std::cout << "Processing " << here << std::endl;
  
  const SBSSimEvent* simEvent = reinterpret_cast<const SBSSimEvent*>(buffer);
  
  
  
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
  event_num = simEvent->EvtID;//++;
  recent_event = event_num;

  // Event weight
  fWeight = simEvent->ev_sigma*simEvent->ev_solang;

  //
  if( fDoBench ) fBench->Begin("physics_decode");
  
  //Bool_t newclus;
  //Int_t crate, slot, chan,lchan;

  std::vector<std::map<Decoder::THaSlotData*, std::vector<UInt_t> > > detmaps;
  detmaps.resize(fDetectors.size());
  
  for(size_t d = 0; d<fDetectors.size(); d++){
    if(fDebug>2)cout << fDetectors[d] << endl;
    //SBSDigSim::UHitData_t* HitData_Det = simEvent->HitDataDet.at(fDetectors[d]);
    LoadDetector(detmaps[d], fDetectors[d], simEvent);
  }
  
  // Now call LoadSlot for the different detectors
  for(size_t d = 0; d < fDetectors.size(); d++) {
    if(fDebug>2)cout << " " << fDetectors[d] << endl;
    for( std::map<Decoder::THaSlotData*, std::vector<UInt_t> >::iterator it =
	   detmaps[d].begin(); it != detmaps[d].end(); ++it) {
      if(it->first->GetModule()==0) {
        if(fDebug>2) {
	  std::cout << "No data available for detector "
		    << fDetectors[d] << std::endl;
        }
      } else {
	if(fDebug>2){
	  std::cout << "load crate/slot: " << it->first->getCrate() << "/" << it->first->getSlot() << " it->second = {";
	  for(size_t k = 0; k<it->second.size(); k++)std::cout << it->second[k] << " ; ";
	  std::cout << " } " << std::endl;
	}
        it->first->GetModule()->LoadSlot(it->first,
					 it->second.data(),0,it->second.size() );
      }
    }
  }
  
  return HED_OK;
}


//Utilities
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
				   std::vector<UInt_t> > &map,
				   const std::string detname, 
				   const SBSSimEvent* simev)
{
  if(fDebug>1)std::cout << "SBSSimDecoder::LoadDectector(" << detname << ")" << std::endl;
  //int detid = detinfo.DetUniqueId();
  Int_t crate, slot;
  unsigned int nwords = 0;
  unsigned short data_type = 0, chan = 0, chan_mult = 0;
  int lchan;
  SimEncoder::mpd_data tmp_mpd;
  UInt_t* mpd_hdr = new UInt_t[2];
  
  Decoder::THaSlotData *sldat = 0;
  //This should be *general* and work for *every* subsystem
  // Loop over all raw data in this event
  UInt_t j = 0;
  //FIXME: we don't want that, I just set it up this way for the sake of going forward
  //Simple fix (might not be ideal): do "if(detname=="xyz")"
  //cout << detname.c_str() << endl;
  
  
  if(strcmp(detname.c_str(), "bb.ps")==0){
    //cout << " ouh " << detname.c_str() << " " << simev->Earm_BBPSTF1.nhits << " " << simev->Earm_BBPS_Dig.nchan << endl;
    for(int j = 0; j<simev->Earm_BBPS_Dig.nchan; j++){
      //cout << j << " " << simev->Earm_BBPS_Dig.chan->at(j) << " " << simev->Earm_BBPS_Dig.adc->at(j) << endl;
      lchan = simev->Earm_BBPS_Dig.chan->at(j);
      ChanToROC(detname, lchan, crate, slot, chan);
      
      if( crate >= 0 || slot >=  0 ) {
	sldat = crateslot[idx(crate,slot)];
      }
      std::vector<UInt_t> *myev = &(map[sldat]);
      
      myev->push_back(SBSSimDataEncoder::EncodeHeader(6, chan, 1));
   
      myev->push_back(simev->Earm_BBPS_Dig.adc->at(j));
      
      if(fDebug>2){
	std::cout << " j = " << j << " my ev = {";
	for(size_t k = 0; k<myev->size(); k++)std::cout << myev->at(k) << " ; ";
	std::cout << " } " << std::endl;
      }
    }
  }
  if(strcmp(detname.c_str(), "bb.sh")==0){
    //cout << " ouh " << detname.c_str() << " " << simev->Earm_BBSHTF1.nhits << " " << simev->Earm_BBSH_Dig.nchan << endl;
    for(int j = 0; j<simev->Earm_BBSH_Dig.nchan; j++){
      //cout << j << " " << simev->Earm_BBSH_Dig.chan->at(j) << " " << simev->Earm_BBSH_Dig.adc->at(j) << endl;
      lchan = simev->Earm_BBSH_Dig.chan->at(j);
      ChanToROC(detname, lchan, crate, slot, chan);
      
      if( crate >= 0 || slot >=  0 ) {
	sldat = crateslot[idx(crate,slot)];
      }
      std::vector<UInt_t> *myev = &(map[sldat]);
      
      myev->push_back(SBSSimDataEncoder::EncodeHeader(6, chan, 1));
   
      myev->push_back(simev->Earm_BBSH_Dig.adc->at(j));
      
      if(fDebug>2){
	std::cout << " j = " << j << " my ev = {";
	for(size_t k = 0; k<myev->size(); k++)std::cout << myev->at(k) << " ; ";
	std::cout << " } " << std::endl;
      }
    }
  }
  if(strcmp(detname.c_str(), "bb.hodo")==0){
    //cout << " ouh " << detname.c_str() << " " << simev->Earm_BBHodoScint.nhits << " " << simev->Earm_BBHodo_Dig.nchan << endl;
    for(int j = 0; j<simev->Earm_BBHodo_Dig.nchan; j++){
      //cout << j << " " << simev->Earm_BBHodo_Dig.chan->at(j) << " " << simev->Earm_BBHodo_Dig.adc->at(j) << endl;
      lchan = simev->Earm_BBHodo_Dig.chan->at(j)+1;
      ChanToROC(detname, lchan, crate, slot, chan);
      
      if( crate >= 0 || slot >=  0 ) {
	sldat = crateslot[idx(crate,slot)];
      }
      std::vector<UInt_t> *myev = &(map[sldat]);
      
      myev->push_back(SBSSimDataEncoder::EncodeHeader(1, chan, 2));
      
      myev->push_back(simev->Earm_BBHodo_Dig.tdc_l->at(j));
      myev->push_back(simev->Earm_BBHodo_Dig.tdc_t->at(j));
      /*
      ChanToROC(detname, lchan, crate, slot, chan);//+91 ??? that might be the trick
      if( crate >= 0 || slot >=  0 ) {
	sldat = crateslot[idx(crate,slot)];
      }
      myev = &(map[sldat]);
      
      myev->push_back(SBSSimDataEncoder::EncodeHeader(8, chan, 1));
      myev->push_back(simev->Earm_BBHodo_Dig.adc->at(j));
      */
      if(fDebug>2){
	std::cout << " j = " << j << " my ev = {";
	for(size_t k = 0; k<myev->size(); k++)std::cout << myev->at(k) << " ; ";
	std::cout << " } " << std::endl;
      }
    }
  }
  if(strcmp(detname.c_str(), "bb.grinch")==0){
    //cout << " ouh " << detname.c_str() << " " << simev->Earm_GRINCH.nhits << " " << simev->Earm_GRINCH_Dig.nchan << endl;
    for(int j = 0; j<simev->Earm_GRINCH_Dig.nchan; j++){
      //cout << j << " " << simev->Earm_GRINCH_Dig.chan->at(j) << " " << simev->Earm_GRINCH_Dig.adc->at(j) << endl;
      lchan = simev->Earm_GRINCH_Dig.chan->at(j);
      ChanToROC(detname, lchan, crate, slot, chan);
      
      if( crate >= 0 || slot >=  0 ) {
	sldat = crateslot[idx(crate,slot)];
      }
      std::vector<UInt_t> *myev = &(map[sldat]);
      
      myev->push_back(SBSSimDataEncoder::EncodeHeader(1, chan, 2));
      
      myev->push_back(simev->Earm_GRINCH_Dig.tdc_l->at(j));
      myev->push_back(simev->Earm_GRINCH_Dig.tdc_t->at(j));
      /*
      ChanToROC(detname, lchan, crate, slot, chan);//+288 ??? that might be the trick
      if( crate >= 0 || slot >=  0 ) {
	sldat = crateslot[idx(crate,slot)];
      }
      myev = &(map[sldat]);
      
      myev->push_back(SBSSimDataEncoder::EncodeHeader(8, chan, 1));
      myev->push_back(simev->Earm_GRINCH_Dig.adc->at(j));
      */
      if(fDebug>2){
	std::cout << " j = " << j << " my ev = {";
	for(size_t k = 0; k<myev->size(); k++)std::cout << myev->at(k) << " ; ";
	std::cout << " } " << std::endl;
      }
    }
  }
  
  if(strcmp(detname.c_str(), "bb.gem")==0){
    //cout << " ouh " << detname.c_str() << endl;
    for(int j = 0; j<simev->Earm_BBGEM_Dig.nstrips; j++){
      lchan = simev->Earm_BBGEM_Dig.strip->at(j);
      //???
      ChanToROC(detname, lchan, crate, slot, chan);
      
      if( crate >= 0 || slot >=  0 ) {
	sldat = crateslot[idx(crate,slot)];
      }
      std::vector<UInt_t> *myev = &(map[sldat]);
      
      //tmp_mpd = lchan/128+simev->Earm_BBGEM_Dig.module->at(j)*12;
      fEncoderMPD->EncodeMPDHeader(tmp_mpd, mpd_hdr, chan);
      myev->push_back(SBSSimDataEncoder::EncodeHeader(9, chan, 6));
      myev->push_back(simev->Earm_BBGEM_Dig.adc_0->at(j));
      myev->push_back(simev->Earm_BBGEM_Dig.adc_1->at(j));
      myev->push_back(simev->Earm_BBGEM_Dig.adc_2->at(j));
      myev->push_back(simev->Earm_BBGEM_Dig.adc_3->at(j));
      myev->push_back(simev->Earm_BBGEM_Dig.adc_4->at(j));
      myev->push_back(simev->Earm_BBGEM_Dig.adc_5->at(j));
    }
    
  }

  if(strcmp(detname.c_str(), "sbs.hcal")==0){
    //cout << " ouh " << detname.c_str() << " " << simev->Harm_HCalScint.nhits << " " << simev->Harm_HCal_Dig.nchan << endl;
    for(int j = 0; j<simev->Harm_HCal_Dig.nchan; j++){
      lchan = simev->Harm_HCal_Dig.chan->at(j);
      ChanToROC(detname, lchan, crate, slot, chan);
      //cout << lchan << " " << crate << " " << slot << " " << chan << endl;

      if( crate >= 0 || slot >=  0 ) {
	sldat = crateslot[idx(crate,slot)];
      }
      std::vector<UInt_t> *myev = &(map[sldat]);
      //cout << SBSSimDataEncoder::EncodeHeader(5, chan, 20) << endl;
      //cout << SBSSimDataEncoder::EncodeHeader(5, chan, 1) << endl;
      myev->push_back(SBSSimDataEncoder::EncodeHeader(5, chan, 20));
      myev->push_back(simev->Harm_HCal_Dig.adc_0->at(j));
      myev->push_back(simev->Harm_HCal_Dig.adc_1->at(j));
      myev->push_back(simev->Harm_HCal_Dig.adc_2->at(j));
      myev->push_back(simev->Harm_HCal_Dig.adc_3->at(j));
      myev->push_back(simev->Harm_HCal_Dig.adc_4->at(j));
      myev->push_back(simev->Harm_HCal_Dig.adc_5->at(j));
      myev->push_back(simev->Harm_HCal_Dig.adc_6->at(j));
      myev->push_back(simev->Harm_HCal_Dig.adc_7->at(j));
      myev->push_back(simev->Harm_HCal_Dig.adc_8->at(j));
      myev->push_back(simev->Harm_HCal_Dig.adc_9->at(j));
      myev->push_back(simev->Harm_HCal_Dig.adc_10->at(j));
      myev->push_back(simev->Harm_HCal_Dig.adc_11->at(j));
      myev->push_back(simev->Harm_HCal_Dig.adc_12->at(j));
      myev->push_back(simev->Harm_HCal_Dig.adc_13->at(j));
      myev->push_back(simev->Harm_HCal_Dig.adc_14->at(j));
      myev->push_back(simev->Harm_HCal_Dig.adc_15->at(j));
      myev->push_back(simev->Harm_HCal_Dig.adc_16->at(j));
      myev->push_back(simev->Harm_HCal_Dig.adc_17->at(j));
      myev->push_back(simev->Harm_HCal_Dig.adc_18->at(j));
      myev->push_back(simev->Harm_HCal_Dig.adc_19->at(j));

      ChanToROC(detname, lchan+288, crate, slot, chan);//+288 ??? that might be the trick

      //cout << lchan+288  << " " << crate << " " << slot << " " << chan << endl;
      if( crate >= 0 || slot >=  0 ) {
	sldat = crateslot[idx(crate,slot)];
      }
      myev = &(map[sldat]);
      
      myev->push_back(SBSSimDataEncoder::EncodeHeader(4, chan, 1));
      myev->push_back(simev->Harm_HCal_Dig.tdc->at(j));
    }

  }
  /*
  while(j < HitData_Det->nhits){
    //Decode header first
    lchan = 0;
    if(HitData_Det->chan->at(j)<0){
      if(fDebug>2)
	std::cout << "j = " << j << " header = " << HitData_Det->dataword->at(j) << std::endl;
      SBSSimDataEncoder::DecodeHeader(HitData_Det->dataword->at(j),
				       data_type,chan_mult,nwords);
      
      //if header if from GEM detector, also decode the MPD header
      if(detname.find("gem")!=std::string::npos){
	for(uint k = 0; k<(HitData_Det->samps_datawords->at(j)).size();k++){
	  mpd_hdr[k] = (HitData_Det->samps_datawords->at(j)).at(k);
	}
	fEncoderMPD->DecodeMPDHeader(mpd_hdr, tmp_mpd);
	//reencode header for GEMs - not sure why - to set "chan" value ?
      }

      if(nwords>0)j++;
    }
    if(fDebug>2)
      std::cout << "j = " << j << " det chan = " << HitData_Det->chan->at(j) << std::endl;
    //channel should *not* be negative (unless there's a problem with nwords...)
    assert(HitData_Det->chan->at(j)>=0);
    //determine crate/slot
    lchan = (int)HitData_Det->chan->at(j);//+chan_mult*fNChan[detname];
    ChanToROC(detname, lchan, crate, slot, chan);

    if(fDebug>2)
      std::cout << "crate " << crate  << " slot " << slot << " chan " << chan << std::endl;
    if(detname.find("gem")!=std::string::npos){
      fEncoderMPD->EncodeMPDHeader(tmp_mpd, mpd_hdr, chan);
    }

    Decoder::THaSlotData *sldat = 0;
    if( crate >= 0 || slot >=  0 ) {
      sldat = crateslot[idx(crate,slot)];
    }
    
    //save the header
    std::vector<UInt_t> *myev = &(map[sldat]);
    myev->push_back(SBSSimDataEncoder::EncodeHeader(data_type,chan,nwords));
    if(detname.find("gem")!=std::string::npos){
      for(int k = 0; k<2;k++){ myev->push_back(mpd_hdr[k]);
      }
    }
    //Then save the hits
    //nwords = n following "hits" for ECal, Cher, Scint;
    //nowrds = n following hits*n data words for HCal, GEMs
    uint i = 0;
    while(i<nwords){
      if(fDebug>2)// || detname.find("grinch")!=std::string::npos 
	std::cout << " i = " << i << " j = " << j << " dataword = " << HitData_Det->dataword->at(j) << std::endl;
      if(detname.find("gem")!=std::string::npos || 
	 detname.find("hcal")!=std::string::npos){
	//if GEM or HCal, loop on "samples datawords" 
	if(HitData_Det->adc->at(j)>-9.e5){
	  // here dataword stores the number of samples datawords
	  for(int k = 0; k<HitData_Det->dataword->at(j);k++, i++){
	    myev->push_back( (HitData_Det->samps_datawords->at(j)).at(k) );
	    if(fDebug>2)
	      std::cout << " samp " << k << " dataword = " << (HitData_Det->samps_datawords->at(j)).at(k) << std::endl;
	  }
	}else{
	  //if adc has dummy value , it is a HCal TDC
	  myev->push_back(HitData_Det->dataword->at(j));
	}
      }else{
	//straightforward for detectors other than GEMs, HCal.
	myev->push_back(HitData_Det->dataword->at(j));
      }
      i++;
      j++;
    }
    if(fDebug>2){
      std::cout << " j = " << j << " my ev = {";
      for(size_t k = 0; k<myev->size(); k++)std::cout << myev->at(k) << " ; ";
      std::cout << " } " << std::endl;
    }
  }//end loop on j
  */
  return HED_OK;
}

/*
void SBSSimDecoder::SetDetMapParam(const std::string detname, int cps, int spc, int fs, int fc)
{
  fChansPerSlotDetMap[detname] = cps;
  fSlotsPerCrateDetMap[detname] = spc;
  fFirstSlotDetMap[detname] = fs;
  fFirstCrateDetMap[detname] = fc;
}
*/

void SBSSimDecoder::CheckForEnabledDetectors()
{
  //fDetectors = fManager->GetAllDetInfo();
  if(fDebug>0) {
    for(size_t i = 0; i < fDetectors.size(); i++) {
      std::cout << "Found detector: " << fDetectors[i].c_str() << endl;
      //<< ", ID: " << fDetectors[i].DetUniqueId() << std::endl;
    }
  }
  fCheckedForEnabledDetectors = true;
}

/*
void SBSSimDecoder::SetTree(TTree *t)
{
  if(t==0)return;
  fTree = new digsim_tree(t);
  if(fTree==0)return;
  fTreeIsSet = true;
}
*/

void SBSSimDecoder::SetDetectors()
{
  TIter aiter(gHaApps);
  THaApparatus* app = 0;
  while( (app=(THaApparatus*)aiter()) ){
    TList* listdet = app->GetDetectors();
    TIter diter(listdet);
    TObject* det = 0;
    while( (det=(TObject*)diter()) ){
      cout << "Setting det " << app->GetName() << "." << det->GetName() 
	   << " into SBSSimDecoder" << endl;
      if(strcmp(app->GetDetector(det->GetName())->GetClassName(),"SBSBBTotalShower")==0){
	SBSBBTotalShower* TS = (SBSBBTotalShower*)app->GetDetector(det->GetName());
	AddDetector(Form("%s.%s",app->GetName(), TS->GetShower()->GetName()), 
		    (app->GetDetector(det->GetName()))->GetInitDate());
	AddDetector(Form("%s.%s",app->GetName(), TS->GetPreShower()->GetName()), 
		    (app->GetDetector(det->GetName()))->GetInitDate());
       }else{
	AddDetector(Form("%s.%s",app->GetName(), det->GetName()), 
		    (app->GetDetector(det->GetName()))->GetInitDate());
      }
    }
  }
}

Int_t SBSSimDecoder::AddDetector(std::string detname, TDatime date)
{
  fDetectors.push_back(detname);
  return ReadDetectorDB(detname, date);
}

Int_t SBSSimDecoder::ReadDetectorDB(std::string detname, TDatime date)
{
  //EPAF: in here the det name is the "full" det name i.e. including the spectro name
  std::string path = std::string(std::getenv("SBS"))+"/DB/";
  if(std::getenv("DB_DIR")) {
    path = std::string(std::getenv("DB_DIR"))+"/";
  }
  const string& fileName = path+"db_"+detname+".dat";
  
  const string prefix = detname+".";
  // First, open the common db file and parse info there, later, the
  // digitization specific db can be used to override any values
  FILE* file  = THaAnalysisObject::OpenFile(fileName.c_str(), date);
  
  std::vector<int> detmap,chanmap;//, detmap_adc;
  uint nchan, nlogchan = 0, chanmapstart = 0;
  
  //int cps, spc, fs, fc;
  
  bool isgem = (detname.find("gem")!=std::string::npos);
  int apv_num = -1;
  
  DBRequest request[] = {
    {"nchan", &nchan, kInt, 0, false},// 
    {"nlog_chan", &nlogchan, kInt, 0, true},// <- optional
    {"detmap", &detmap, kIntV, 0, false}, //
    {"chanmap", &chanmap, kIntV, 0, true}, // <- optional
    {"chanmap_start", &chanmapstart, kInt, 0, true}, // <- optional
    //{"detmap_adc", &detmap_adc, kIntV, 0, true}, // <- optional
    /*
    {"first_crate", &fc, kInt, 0, true},// <- optional 
    {"first_slot", &fs, kInt, 0, true},//  <- optional
    {"chan_per_slot", &cps, kInt, 0, true},//  <- optional
    {"slot_per_crate", &spc, kInt, 0, true},//  <- optional
    */
    { 0 }
  };
  Int_t err = THaAnalysisObject::LoadDB(file, date, request, prefix.c_str());
  // Could close the common file already
  fclose(file);
  
  if(nlogchan==0)nlogchan = nchan;
  
  if(err)return THaAnalysisObject::kInitError;
  
  fNChanDet[detname] = nchan;
  fChanMapStartDet[detname] = chanmapstart;
  (fInvDetMap[detname]).resize(nlogchan);
  int nparam_mod = 4;
  if(detmap[4]==-1)nparam_mod = 5;
  int crate,slot,ch_lo,ch_hi, ch_count = 0, ch_map = 0;
  for(size_t k = 0; k < detmap.size(); k+=nparam_mod) {
    crate  = detmap[k];
    slot   = detmap[k+1];
    ch_lo  = detmap[k+2];
    ch_hi  = detmap[k+3];
    /*
    if(detname.find("hodo")!=std::string::npos)
      cout << " crate " << crate << " slot " << slot 
	   << " ch_lo " << ch_lo << " ch_hi " << ch_hi << endl;
    */
    if(chanmap.empty()){
      for(int i = ch_lo; i<=ch_hi; i++, ch_count++){
	if(isgem && i%128==0){
	  apv_num++;
	  cout << crate << " " << slot << " " << i << " " << apv_num << endl;
	}
	if(ch_count>nlogchan){
	  std::cout << " number of channels defined in detmap ( >= " << ch_count << ") exceeds logical number of channels = " << nlogchan << std::endl;
	  return THaAnalysisObject::kInitError;
	}
	(fInvDetMap[detname])[ch_count]=detchaninfo(crate, slot, i, apv_num);
	/*
	if(detname.find("hodo")!=std::string::npos){
	  cout << " crate " << crate << " slot " << slot 
	       << " i " << i << " ch_count " << ch_count << endl;
	  cout << &(fInvDetMap.at(detname)).at(ch_count) << endl;
	}
	*/
      }
    }else{
      
      for(int i = ch_lo; i<=ch_hi; i++, ch_map++){
	if(ch_count>nlogchan){
	  std::cout << " number of channels defined in detmap ( >= " << ch_count << ") exceeds logical number of channels = " << nlogchan << std::endl;
	    return THaAnalysisObject::kInitError;
	}
	if(fDebug>=2)std::cout << " i = " << i << ", crate = " << crate << ", slot = " << slot <<  ", ch_count = " << ch_count << " chan = " << chanmap[ch_map]-1 << " (+" << nchan << ") " << std::endl;
	if(chanmap[ch_map]>=0){
	  if(ch_count<nchan){
	    (fInvDetMap[detname])[chanmap[ch_map]-1]=detchaninfo(crate, slot, i);
	    if(fDebug>=3)std::cout << chanmap[ch_map]-1 << " " << &(fInvDetMap.at(detname)).at(chanmap[ch_map]-1) << std::endl;
	  }else{
	    (fInvDetMap[detname])[chanmap[ch_map]+nchan-1]=detchaninfo(crate, slot, i);
	    if(fDebug>=3)std::cout <<&(fInvDetMap.at(detname)).at(chanmap[ch_map]+nchan-1) << std::endl;
	  }
	  ch_count++;
	}
      }
    }
  }
  
  /*
  fChansPerSlotDetMap[detname] = cps;
  fSlotsPerCrateDetMap[detname] = spc;
  fFirstSlotDetMap[detname] = fs;
  fFirstCrateDetMap[detname] = fc;
  */
  return(THaAnalysisObject::kOK);
}


//-----------------------------------------------------------------------------
//static inline
void SBSSimDecoder::ChanToROC(const std::string detname, Int_t h_chan,
			       Int_t& crate, Int_t& slot, UShort_t& chan )const 
{
  // Convert location parameters (row, col, chan) of the given Channel
  // to hardware channel (crate,slot,chan)
  // The (crate,slot,chan) assignment must match the detmap definition in
  // the database!  See TreeSearch/dbconvert.cxx
  // In the case of GRINCH/RICH: 
  // crate = GTP; slot = VETROC; chan = PMT. (NINOs are "transparent", in a similar way to the MPDs)
  
  /*
  int CPS = fChansPerSlotDetMap.at(detname);
  int SPC = fSlotsPerCrateDetMap.at(detname);
  int FS = fFirstSlotDetMap.at(detname);
  int FC = fFirstCrateDetMap.at(detname);
  
  //div_t d = div( h_chan, fManager->GetChanPerSlot() );
  div_t d = div( h_chan, CPS );
  slot = d.quot;
  chan = d.rem;

  d = div( slot, SPC );
  crate = d.quot+FC;
  slot  = d.rem+FS;
  */

  if(fDebug>3){
    std::cout << " " << detname << " "  << h_chan << " " << &fInvDetMap.at(detname) << " " << std::endl;
    std::cout << &(fInvDetMap.at(detname)).at(h_chan) << std::endl;
  }
  crate = ((fInvDetMap.at(detname)).at(h_chan)).crate;
  slot = ((fInvDetMap.at(detname)).at(h_chan)).slot;
  chan = ((fInvDetMap.at(detname)).at(h_chan)).chan;
  
}

int APVnum(const std::string detname, 
	   Int_t crate, Int_t slot, Int_t chan)
{
  
}

/*
//-----------------------------------------------------------------------------
static inline
Int_t MakeROCKey( Int_t crate, Int_t slot, Int_t chan )
{
  return chan;// +
  //fManager->GetChanPerSlot()*( slot + fManager->GetSlotPerCrate()*crate );
}

//-----------------------------------------------------------------------------
Int_t SBSSimDecoder::ChanFromROC( Int_t crate, Int_t slot, Int_t chan ) const
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
