/////////////////////////////////////////////////////////////////////
//
//   MPDModule
//   This is the MPD module decoder; based on SkeletonModule
//   (https://github.com/JeffersonLab/analyzer)
//
//   E. Cisbani
//   Original Version:   2015/Dec
//   
//   v5 Version based on documentation by Paulo Musico
//   Seamus Riordan
//   sriordan@anl.gov
//   Aug 31, 2018
//
//   v5 Version with online SSP zero supression, based on documentation
//   from Ben (DAQ group) which I got from Danning Di.
//   Juan Carlos Cornejo <cornejo@jlab.org> - 2018/10/23
//
/////////////////////////////////////////////////////////////////////

/*
  #define MPD_VERSION_TAG 0xE0000000
  #define MPD_EVENT_TAG   0x10000000
  #define MPD_MODULE_TAG  0x20000000
  #define MPD_ADC_TAG     0x30000000
  #define MPD_HEADER_TAG  0x40000000
  #define MPD_DATA_TAG    0x0
  #define MPD_TRAILER_TAG 0x50000000
*/

// #define SSP_DATADEF(b) ((b&0x80000000)>>31)
// #define SSP_TAG(b)     ((b&0x78000000)>>27)
// #define SSP_SAMPLE(b,c) ((b>>c)&0xFFF)|(((b>>c)&0x1000)?0xFFFFF000:0x0)

#include "MPDModule.h"
#include "THaSlotData.h"
#include "THaCrateMap.h"
#include <limits>
#include <vector>
#include <map>
#include <set>
#include <iostream>

using namespace std;

namespace Decoder {

  Module::TypeIter_t MPDModule::fgThisType =
    DoRegister( ModuleType( "Decoder::MPDModule" , 3561 ));

  MPDModule::MPDModule(Int_t crate, Int_t slot) : VmeModule(crate, slot) {
    fDebugFile=nullptr;
    Init(); //Should this be called here? not clear...
    
  }
  
  MPDModule::~MPDModule() {
    
  }
  
  void MPDModule::Init() { 
    VmeModule::Init();
    //    Config(0,25,6,16,128); // should be called by the user (but how?)
    fDebugFile=nullptr;
    Clear();
    //    fName = "MPD Module (INFN MPD for GEM and more), use Config to dynamic config";
    fName = "MPD Module";

    fOnlineZeroSuppression = false; //If this is false, then we want to calculate and subtract the common-mode from each ADC sample:
    fBlockHeader = 0;  //0000 = 0
    fBlockTrailer = 1; //0001 = 1
    fEventHeader = 2;  //0010 = 2
    fTriggerTime = 3;  //0011 = 3
    fMPDFrameHeader = 5; //0101 = 5
    fMPDEventInfo = 12; //1100 = 12
    fMPDDebugHeader = 13; //1101 = 13
    fDataNotValid = 14;  //1110 = 14
    fFillerWord = 15; //1111 = 15
    //    fAPVHeader   = 0x4;

    fSLOTID_VTP = 11;
    
    fNumSample = 6;

    fChan_CM_flags = 640; //reference channel for common-mode and zero suppression flags
    fChan_TimeStamp_low = 641;
    fChan_TimeStamp_high = 642;
    fChan_MPD_EventCount = 643;
    fChan_MPD_Debug = 644;
    
  }

  void MPDModule::Init( const char *configstr ) { //parse (optional) configuration parameters:
    Init();
    //allow users to configure the various dummy/reference channels via the crate map:
    vector<ConfigStrReq> req = { {"chan_cmflags", fChan_CM_flags},
				 {"chan_timestamp_low", fChan_TimeStamp_low},
				 {"chan_timestamp_high", fChan_TimeStamp_high},
				 {"chan_event_count", fChan_MPD_EventCount},
				 {"chan_MPD_debug", fChan_MPD_Debug} };
    ParseConfigStr(configstr, req);

    assert( fChan_CM_flags < THaCrateMap::MAXCHAN );
    assert( fChan_TimeStamp_low < THaCrateMap::MAXCHAN );
    assert( fChan_TimeStamp_high < THaCrateMap::MAXCHAN );
    assert( fChan_MPD_EventCount < THaCrateMap::MAXCHAN );
    assert( fChan_MPD_Debug < THaCrateMap::MAXCHAN );
  }

  //This version ASSUMES that there is no online zero suppression, so all 128 APV channels are present in every event!
  
  UInt_t MPDModule::LoadSlot( THaSlotData *sldat, const UInt_t *evbuffer, UInt_t pos, UInt_t len ){

    //std::cout << "Calling MPDModule::LoadSlot()... (pos, len) = " << pos << ", " << len << "..."; 
    //AJRP: LoadSlot method for MPD SSP data format to be used in Hall A during the GMN run:
    const UInt_t *datawords = &(evbuffer[pos]);

    fWordsSeen = 0;

    Int_t status;
    
    UInt_t thisword;

    UInt_t iword=0;

    //Get the slot number for this call to LoadSlot:
    UInt_t this_slot = sldat->getSlot();
    //UInt_t this_crate = sldat->getCrate();
    
    //UInt_t mask, shift;
    UInt_t slot=MAXSLOT+1, apv_id=0, apv_chan=0, fiber=0;
    //UInt_t mpd_id=0, apv_chan_40, apv_chan65;
    
    //UInt_t prev_slot=0;

    UInt_t effChan=0;
    //UInt_t effChan_old = MAXCHAN+1; //set "old" effective channel to something impossible so that the first channel will always trigger "loading" of the cm_flags:

    //The use of a map here insures that the cm_flags will only be "loaded" into the "slot" once per APV card:
    std::map<UInt_t, UInt_t> cm_flags_vs_chan; //key = "effective channel", mapped value = common-mode flags (and possibly other information)
    std::map<UInt_t, UInt_t> TimeStampL_vs_fiber; //least significant 24 bits of time stamp
    std::map<UInt_t, UInt_t> TimeStampH_vs_fiber; //most significant 24 bits of time stamp
    std::map<UInt_t, UInt_t> EventCount_vs_fiber; //20-bit "event count" variable:
    std::map<UInt_t, std::vector<UInt_t> > MPDdebugInfo_vs_chan; //TBD: key = effective channel, mapped value = vector of (3) data words encoding calculated common-mode corrections for six time samples:
    
    UInt_t TIMESTAMP_LO = 0, TIMESTAMP_HI=0, EVENT_COUNT=0;
    UInt_t eventinfo_wordcount=0;
    
    bool CM_OR=false;
    bool ENABLE_CM=false;
    bool BUILD_ALL_SAMPLES=true; 

    bool is_SSP = false; //assume VTP by default which only has the MPD data, no block header/trailer.
    
    bool found_this_slot = false;
    bool found_MPD_header = false;
    
    UInt_t mpd_strip_count = 0;
    UInt_t mpd_word_count = 0;

    //UInt_t old_type_tag=16;
    UInt_t type_tag=16; //intialize to something that is NOT one of the recognized data types

    
    //following the MPDRawParser in ROOT_GUI_multicrate, loop on all the data in the ROC bank (which corresponds to one "crate"), and populate the
    //temporary data structure above, ONLY if slot == this_slot
    // Since all the data from one slot is (or should be) in a contiguous block, we should 

    //temporary storage for information needed to store "strip hits"
    //UInt_t 
    UInt_t ADCsamples[fNumSample]; //temporary storage for ADC time samples. will this compile? Hope so... 

    UInt_t hitwords[3]; //temporary storage for the three words needed to extract the information for one "hit"

    UInt_t CMwordcount = 0;
    UInt_t CMwords[3]; //temporary storage for the three words containing the common-mode values for the six time samples of an APV:
    
    while( iword < len ){
      thisword = datawords[iword++];

      //check whether this is a data-type defining or data-type continuation word:
      UInt_t word_type = (thisword & 0x80000000)>>31;
      
      if( word_type == 1 ){ //data-type defining: extract data type from bits 30-27:
	//old_type_tag = type_tag;
	
	type_tag = (thisword & 0x78000000)>>27;
	//std::cout << "Data-type defining word, type tag, fBlockHeader = " << type_tag << ", " << fBlockHeader << std::endl;
	if( type_tag == fBlockHeader ){ //if we see a block header, this is SSP data:

	  is_SSP = true;
	  
	  found_MPD_header = false; 
	  
	  //prev_slot = slot;
	  //extract "SLOTID" from bits 26-22 and compare to this_slot
	  slot = (thisword & 0x07C00000)>>22; //7C = 0111 1100
	  if( slot == this_slot ) found_this_slot = true;

	  //std::cout << "Found block header, SLOTID = " << slot << std::endl;
	}

	if( type_tag == fBlockTrailer ){ //end of a block of SSP data. set the is_SSP flag to false unless and until we see another block header word
	  is_SSP = false; //
	  found_MPD_header = false;
	}
	
	
	//Question: ask Ben: how would/could we use the trigger time words to correct for APV trigger jitter?
	
	if( type_tag == fMPDFrameHeader ){ //extract "flags", FIBER, and MPD_ID
	  //TO-DO: figure out how to store this information in the slot data
	  //(maybe define some arbitrary "dummy" channel to hold this info, as well as trigger time words)
	  ENABLE_CM = TESTBIT( thisword, 26 );        
	  BUILD_ALL_SAMPLES = TESTBIT(thisword, 25 );
	  CM_OR = TESTBIT(thisword, 24 );

	  //NEW firmware to support up to 40 MPDs per VTP
	  //to extract bits 21-16:
	  // 0x003F0000 = 0000 0000 0011 1111 0000 0000 0000 0000
	  
	  //FIBER number was in bits 20-16:
	  //fiber = (thisword & 0x001F0000)>>16;

	  //NOW the fiber number is in bits 21-16:
	  fiber = (thisword & 0x003F0000)>>16;
	  
	  //MPD ID is in bits 0-4, but we basically ignore it:
	  //mpd_id = (thisword & 0x0000001F);

	  found_MPD_header = true;

	  //now how should we "load" the data into the "slot"?
	  //Each "hit" in this channel
	  // UInt_t flags = 2*ENABLE_CM + BUILD_ALL_SAMPLES;
	  
	  // status = sldat->loadData( fChan_CM_flags, flags, fiber );
	  // std::cout << "found MPD frame header, fiber, mpd_id, ENABLE_CM, BUILD_ALL_SAMPLES, is_SSP, SLOT, THIS_SLOT = " << fiber << ", " << mpd_id << ", "
	  // 	    << ENABLE_CM << ", " << BUILD_ALL_SAMPLES << ", "
	  // 	    << is_SSP << ", " << slot << ", " << this_slot << std::endl;
	  
	  //reset "word" and "strip" counters:
	  mpd_word_count = 0; 
	  mpd_strip_count = 0;

	  if( !is_SSP ) {
	    slot = fSLOTID_VTP; //always 11 
	    found_this_slot = ( slot == this_slot );
	  }
	}

	if( type_tag == fMPDEventInfo ){ //if we see this, next three words give coarse timestamp, fine timestamp, and event counter:
	  eventinfo_wordcount = 1;
	  //The lower 16 bits of coarse timestamp and the 8-bit fine timestamp are in
	  // bits 0-23 of thisword:
	  TIMESTAMP_LO = thisword & 0x00FFFFFF;
	  
	}

	if( type_tag == fMPDDebugHeader ){
	  CMwordcount = 0;
	  CMwords[CMwordcount] = thisword;
	}
	
      } else if( found_this_slot ){
	//data-type continuation: behavior depends on type_tag. If the most recently found "slot" 
	//doesn't match the one that we want, then do nothing. 
	//For NOW, let's focus mainly on the MPD Frame decoding and worry about anything and everything else later:

	if( type_tag == fMPDFrameHeader && found_MPD_header ){ //the data continuation words should come in bundles of three * N, where N is the total number of "hits" (i.e., fired strips) in the MPD:
	  
	  hitwords[mpd_word_count%3] = thisword; //Set hitwords before incrementing word count
	  
	  mpd_word_count++; //increment word count;

	  mpd_strip_count = mpd_word_count/3; 

	  //this should take care of the issue of missing the last "hit":
	  if( mpd_word_count%3 == 0 && mpd_strip_count > 0 ){ //extract information from the three "hit words" and "load" the data into the "slot":

	    // load up the ADC samples:
	    for( int iw=0; iw<3; iw++ ){
	      
	      // In words: take the bitwise OR of the first 12 bits of hitwords[iw] with 0xFFFFF000 or 0x0
	      // depending on the value of bit 12 (13th bit) of hitwords[iw].
	      // This is essentially implementing the two's complement representation of a 13-bit signed integer 
	      ADCsamples[2*iw] = ( hitwords[iw] & 0xFFF ) | ( ( hitwords[iw] & 0x1000 ) ? 0xFFFFF000 : 0x0 );
	      ADCsamples[2*iw+1] = ( (hitwords[iw]>>13) & 0xFFF ) | ( ( (hitwords[iw]>>13) & 0x1000 ) ? 0xFFFFF000 : 0x0 );
	    }

	    //APV_ID is in bits 26-30 of hitword[2]:
	    // 0x 0111 1100 .... = 0x7C000000
	    apv_id = (hitwords[2] & 0x7C000000)>>26;

	    //APV channel num (bits 4:0) is in bits 26:30 of hitword[0];
	    //APV channel num (bits 6:5) is in bits 26:30 of hitword[1] (but we actually only want bits 26 and 27 AFAIK:
	    UInt_t apv_chan40 = (hitwords[0] & 0x7C000000) >> 26;
	    UInt_t apv_chan65 = (hitwords[1] & 0x0C000000) >> 26;
	    apv_chan = (apv_chan65 << 5) | apv_chan40;
	    
	    effChan = (fiber << 4) | apv_id;  

	    UInt_t cm_flags = 4*CM_OR + 2*ENABLE_CM + BUILD_ALL_SAMPLES;

	    cm_flags_vs_chan[effChan] = cm_flags;
	    
	    for( int is=0; is<6; is++ ){

	      // if( ADCsamples[is] > 0xFFF && !ENABLE_CM ) {
	      // 	std::cout << "negative ADC sample encountered when CM not enabled, raw data word = " << ADCsamples[is] << ", signed int representation = " << Int_t(ADCsamples[is]) << std::endl; 
	      // 	std::cout << "This is not expected" << std::endl;
	      // }
	      // This loads each of the six ADC samples as a new "hit" into sldat, with "effChan" as the unique "channel" number,
	      // the ADC samples as the "data", and the APV channel number as the "rawData"
	      // std::cout << "decoded one strip hit: (crate, slot, fiber, apv_id, apv_chan, effChan, isample, ADCsamples[isample]) = ("
	      // 		<< sldat->getCrate() << ", " << slot << ", " << fiber << ", " << apv_id << ", " << apv_chan << ", " << effChan << ", "
	      // 		<< is <<  ", " << int(ADCsamples[is]) << ")" << std::endl;

	      //Since we need only two bits to encode the ENABLE_CM and BUILD_ALL_SAMPLES flags, we can
	      
	      
	      status = sldat->loadData( "adc", effChan, ADCsamples[is], apv_chan );
	      if( status != SD_OK ) return -1;
	    }
	  }
	  
	}

	if( type_tag == fMPDEventInfo && found_MPD_header ){
	  if( eventinfo_wordcount == 1 ){ //upper 24 bits of the coarse timestamp:
	    TIMESTAMP_HI = thisword & 0x00FFFFFF;
	    eventinfo_wordcount++;
	  } else if( eventinfo_wordcount == 2 ){
	    // event counter is in the first 20 bits:
	    EVENT_COUNT = thisword & 0x000FFFFF;

	    //fill map with most recently reported MPD fiber number as the key:
	    TimeStampL_vs_fiber[fiber] = TIMESTAMP_LO;
	    TimeStampH_vs_fiber[fiber] = TIMESTAMP_HI;
	    EventCount_vs_fiber[fiber] = EVENT_COUNT;
	    
	    eventinfo_wordcount = 0;
	  }
	}

	if( type_tag == fMPDDebugHeader && found_MPD_header ){
	  CMwordcount++;
	  CMwords[CMwordcount%3] = thisword;
	  if( CMwordcount == 2 ){ 
	    MPDdebugInfo_vs_chan[effChan].clear();
	    for( int iw=0; iw<3; iw++ ){
	      MPDdebugInfo_vs_chan[effChan].push_back( CMwords[iw] );
	    }
	  }
	}
      }
      
      fWordsSeen++;
    }

    //std::cout << "Loading common-mode flags for this event: " << std::endl;
    //now load the cm flags:
    for( auto iapv = cm_flags_vs_chan.begin(); iapv != cm_flags_vs_chan.end(); ++iapv ){

      //std::cout << "effChan = " << iapv->first << ", ENABLE_CM = " << (iapv->second / 2) << ", BUILD_ALL_SAMPLES = " << (iapv->second % 2 ) << std::endl;
      //The "flags" word is the mapped value (iapv->second)
      //The effective channel is the key (iapv->first)
      sldat->loadData( fChan_CM_flags, iapv->second, iapv->first );
    }

    //It is ASSUMED that time stamp low, time stamp high, and event count 
    for( auto ifiber = TimeStampL_vs_fiber.begin(); ifiber != TimeStampL_vs_fiber.end(); ++ifiber ){
      UInt_t fiber = ifiber->first;
      sldat->loadData( fChan_TimeStamp_low, TimeStampL_vs_fiber[fiber], fiber );
      sldat->loadData( fChan_TimeStamp_high, TimeStampH_vs_fiber[fiber], fiber );
      sldat->loadData( fChan_MPD_EventCount, EventCount_vs_fiber[fiber], fiber );
    }

    //Now loop over all APVs seen in the data and load the MPD debug headers into the appropriate dummy channel:
    for( auto iapv = MPDdebugInfo_vs_chan.begin(); iapv != MPDdebugInfo_vs_chan.end(); ++iapv ){
      UInt_t chan = iapv->first;
      for( UInt_t iw=0; iw<MPDdebugInfo_vs_chan[chan].size(); iw++ ){
	sldat->loadData( fChan_MPD_Debug, MPDdebugInfo_vs_chan[chan][iw], chan );
      }
    }
    
    //std::cout << "Finished MPDModule::LoadSlot, fWordsSeen = " << fWordsSeen << std::endl;
    //std::cout << "MPDModule::LoadSlot done." << std::endl;
    return fWordsSeen;
  
  }
  
  //
  // UInt_t MPDModule::LoadSlot( THaSlotData *sldat, const UInt_t* evbuffer, UInt_t pos, UInt_t len) {
  //   const UInt_t *p = &evbuffer[pos];
  //   //UInt_t data;
  //   fWordsSeen = 0;

  //   // From stand alone decoder
  //   // We declare an effective channel from the MPD ID 
  //   // and ADC channel
  //   Int_t ch, status;
  //   Int_t mpdID = -1;
  //   Int_t adcCh = -1;
  //   Int_t effCh = 0;

  //   UInt_t data_count = 0;

  //   //  v5 decoder (with SSP online zero suppression)
  //   UInt_t ii,jj,kk,mm; // loop indices: ii (event), jj(word), kk(mpd), mm (block)
  //   UInt_t thesewords;
  //   UInt_t hit[3] = {0};
  //   UInt_t sample_dat[6] = {0U};

  //   jj =  0;
  //   while( jj < len ){
  //     mm = jj;
  //     thesewords = p[jj++];
  //     //printf("===============================================================================\n");
  //     //printf("=    CRATE   %d   ======    SLOT   %d   =======================================\n", fCrate, fSlot);
  //     //printf("BLOCK HEADER       0x%08x\n", thesewords);
  //     //printf("Data defining? (1) %d\n", (thesewords & 0x80000000) >> 31);
  //     //printf("Type (0)           %d\n", (thesewords & 0x78000000) >> 27);
  //     //printf("BLOCK NUMBER       %d\n", (thesewords & 0x0003FF00) >> 8 );
  //     //printf("EVENT_PER_BLOCK    %d\n", (thesewords & 0x000000FF) >> 0 );
  //     //printf("\n");

  //     if( (SSP_DATADEF(thesewords) != 1) || (SSP_TAG(thesewords) != 0 )) {
  // 	fprintf(stderr, "[ERROR  MPDModule::LoadSlot, line %d] "
  // 		"BLOCK HEADER NOT FOUND\n", __LINE__);
  // 	return -1;
  //     }

  //     UInt_t nevent = (thesewords&0xFF);
  //     // Ensure we have enough data
  //     // (need at least 4 per event: 1 event header + 2 trigger words +
  //     // 1 event trailer
  //     if( nevent > 0 && jj + (nevent*4) >= len) {
  // 	fprintf(stderr, "[ERROR  MPDModule::LoadSlot, line %d] "
  // 		"NOT ENOUGH WORDS TO DECODE THIS EVENT!\n", __LINE__);
  // 	return -1;
  //     }

  //     for( ii = 0; ii < nevent; ii++ ){
  // 	thesewords = p[jj++];
  // 	//printf("EVENT HEADER       0x%08x\n", thesewords);
  // 	//printf("Data defining? (1) %d\n", (thesewords & 0x80000000) >> 31);
  // 	//printf("Type (2)           %d\n", (thesewords & 0x78000000) >> 27);
  // 	//printf("EVENT COUNT        %d\n", (thesewords & 0x3FFFFF) >> 0);
  // 	//printf("\n");
  // 	if( (SSP_DATADEF(thesewords) != 1) || (SSP_TAG(thesewords) != 2 )) {
  // 	  fprintf(stderr, "[ERROR  MPDModule::LoadSlot, line %d] EVENT HEADER NOT FOUND\n", __LINE__);
  // 	  return -1;
  // 	}

  // 	thesewords = p[jj++] & 0xFFFFFFFF;
  // 	//printf("TRIGGER TIME 1     0x%08x\n", thesewords);
  // 	//printf("Data defining? (1) %d\n", (thesewords & 0x80000000) >> 31);
  // 	//printf("Type (3)           %d\n", (thesewords & 0x78000000) >> 27);
  // 	//printf("COURSE TIME        %d\n", (thesewords & 0xFFFFFF) >> 0);
  // 	//printf("\n");
  // 	if( (SSP_DATADEF(thesewords) != 1) || (SSP_TAG(thesewords) != 3 )) {
  // 	  fprintf(stderr, "[ERROR  MPDModule::LoadSlot, line %d] TRIGGER TIME 1 WORD NOT FOUND\n", __LINE__);
  // 	  return -1;
  // 	}

  // 	thesewords = p[jj++] & 0xFFFFFFFF;
  // 	//printf("TRIGGER TIME 2     %08x\n", thesewords);
  // 	//printf("Data defining? (0) %d\n", (thesewords & 0x80000000) >> 31);
  // 	//printf("COURSE TIME        %d\n", (thesewords & 0xFFFFFF) >> 0);
  // 	//printf("\n");
  // 	if( (SSP_DATADEF(thesewords) != 0) ) {
  // 	  fprintf(stderr, "[ERROR  MPDModule::LoadSlot, line %d] TRIGGER TIME 2 WORD NOT FOUND\n", __LINE__);
  // 	  return -1;
  // 	}

  // 	// Loop through all MPD fiber headers
  // 	while( (p[jj]&0xF8000000)>>27 == 0x15 ) {
  // 	  // First word defines the tag type 5, and the MPD ID (fiber number)
  // 	  kk = 0;
  // 	  //printf("\n[Starting sample %d]\n", kk);

  // 	  thesewords = p[jj++];

  // 	  //printf("MPD HEADER        0x%08x\n", thesewords);
  // 	  //printf("Data defining? (1) %d\n", (thesewords & 0x80000000) >> 31);
  // 	  //printf("Type (5)           %d\n", (thesewords & 0x78000000) >> 27);
  // 	  //printf("MPD Fiber Number   %d\n", (thesewords & 0x0000001F) >> 0 );
  // 	  //printf("\n");

  // 	  mpdID = thesewords & 0x1F;

  // 	  // Now loop through each of the APV hits in this MPD
  // 	  while( SSP_DATADEF(p[jj]) != 1) {
  // 	    // For each one of these, we must have at least 3 more
  // 	    // words preceeding
  // 	    if(jj+2>=len) {
  // 	      fprintf(stderr, "[ERROR  MPDModule::LoadSlot, line %d] NOT ENOUGH"
  // 		      " WORDS TO DECODE APV HITS for MPD %d\n", __LINE__,
  // 		      mpdID);
  // 	      return -1;
  // 	    }
  // 	    //printf("samples: ");
  // 	    for(int h = 0; h < 3; h++) {
  // 	      hit[h] = p[jj++];
  // 	      // The samples are stored as 13-bit signed int
  // 	      // This needs to be converted back to typical 32-bit signed int
  // 	      sample_dat[h*2]   = SSP_SAMPLE(hit[h],0);
  // 	      sample_dat[h*2+1] = SSP_SAMPLE(hit[h],13);
  // 	      if(SSP_DATADEF(hit[h]) != 0) {
  // 		fprintf(stderr, "[ERROR  MPDModule::LoadSlot, line %d] MISSING"
  // 			" APV_HIT_WORD%d for APV_HIT%d of MPD=%d, word=0x%x\n", __LINE__,
  // 			h,kk,mpdID,hit[h]);
  // 		return -1;
  // 	      }
  // 	      //printf(" %d %d", sample_dat[h*3], sample_dat[h*3+1]);
  // 	    }
  // 	    //printf("APV HIT0           0x%08x\n", hit[0]);
  // 	    //printf("Data defining? (0) %d\n",   (hit[0] & 0x80000000) >> 31);
  // 	    //printf("APV_CH(b4:b0)      0x%x\n", (hit[0] & 0x7C000000) >> 26 );
  // 	    //printf("ADC_SAMP_T1        %d\n",   (hit[0] & 0x3FFE000) >> 13 );
  // 	    //printf("ADC_SAMP_T0        %d\n",   (hit[0] & 0x1FFF) >> 0 );
  // 	    //printf("APV HIT1           0x%08x\n", hit[1]);
  // 	    //printf("Data defining? (0) %d\n",   (hit[1] & 0x80000000) >> 31);
  // 	    //printf("APV_CH(b6:b5)      0x%x\n", (hit[1] & 0xC000000) >> 26 );
  // 	    //printf("ADC_SAMP_T3        %d\n",   (hit[1] & 0x3FFE000) >> 13 );
  // 	    //printf("ADC_SAMP_T2        %d\n",   (hit[1] & 0x1FFF) >> 0 );
  // 	    //printf("APV HIT2           0x%08x\n", hit[2]);
  // 	    //printf("Data defining? (0) %d\n",   (hit[2] & 0x80000000) >> 31);
  // 	    //printf("APV ID             %d\n",   (hit[2] & 0x7C000000) >> 26 );
  // 	    //printf("ADC_SAMP_T5        %d\n",   (hit[2] & 0x3FFE000) >> 13 );
  // 	    //printf("ADC_SAMP_T4        %d\n",   (hit[2] & 0x1FFF) >> 0 );
  // 	    //printf("\n");

  // 	    // Now decode the hit info

  // 	    // Strip number (APV25 channel number)
  // 	    adcCh = (hit[2]&0x7C000000)>>26;
  // 	    ch = ((hit[0]&0x7C000000)>>26) | ((hit[1]&0xC000000)>>21);
  // 	    effCh = (mpdID) << 8 | adcCh;
  // 	    for(int s = 0; s < 6; s++) {
  // 	      // the raw data will be the strip number
  // 	      status = sldat->loadData("adc",effCh, sample_dat[s], ch);
  // 	      if( status != SD_OK ) return -1;

  // 	      fWordsSeen++;
  // 	      data_count++;
  // 	    }

  // 	    kk++;
  // 	  } // apv_hit loop
  // 	} // mpd loop
  //     } //event loop

  //       // Loop over the filler words
  //     while(p[jj] == 0xF8000000) {
  // 	//printf("FILLER WORD: 0x%08x\n",p[jj]);
  // 	jj++;
  //     }

  //     // Finally, we should have the BLOCK trailer
  //     thesewords = p[jj++];
  //     //printf("BLOCK TRAILER       0x%08x\n", thesewords);
  //     //printf("Data defining? (1) %d\n", (thesewords & 0x80000000) >> 31);
  //     //printf("Type (1)           %d\n", (thesewords & 0x78000000) >> 27);
  //     //printf("NUMBER_OF_WORDS    %d\n", (thesewords & 0x003FFFFF) >> 0 );
  //     //printf("\n");
  //     if( (SSP_DATADEF(thesewords) != 1) || (SSP_TAG(thesewords) != 1 )) {
  // 	fprintf(stderr, "[ERROR  MPDModule::LoadSlot, line %d] "
  // 		"BLOCK TRAILER NOT FOUND\n", __LINE__);
  // 	return -1;
  //     }

  //     //printf("Read number of words %d expected %d\n",jj-mm,thesewords&0x3FFFFF);
  //     if((thesewords&0x3FFFFF) != (jj-mm) ) {
  // 	fprintf(stderr, "[ERROR  MPDModule::LoadSlot, line %d] NUMBER OF "
  //               "WORDS READ %d DOES NOT MATCH NUMBER EXPECTED %d\n", __LINE__,
  //               jj-mm,thesewords&0x3FFFFF);
  // 	return -1;
  //     };

  //   } // block loop


  //     /*
  //     //  v5 decoder (with no SSP zero suppression)
      
  //     int ii,jj,kk,ll;
  //     int thesewords;

  //     jj =  0;


  //     while( jj < len ){
  //     thesewords = p[jj++] & 0xFFFFFF;
  //     //printf("===============================================================================\n");
  //     //printf("=    CRATE   %d   ======    SLOT   %d   =======================================\n", fCrate, fSlot);
  //     //printf("BLOCK HEADER  %06x\n", thesewords);
  //     //printf("Good? (0)       %x\n", (thesewords & 0xe00000) >> 21);
  //     //printf("Module ID       %d\n", (thesewords & 0x1F0000) >> 16 );
  //     //printf("EVENT_PER_BLOCK %d\n", (thesewords & 0x00FF00) >> 8 );
  //     //printf("BLOCK COUNT     %d\n", (thesewords & 0x0000FF) >> 0);

  //     if( (thesewords & 0xe00000) >> 21 != 0 ){
  //     fprintf(stderr, "[ERROR  MPDModule::LoadSlot, line %d] BLOCK HEADER NOT FOUND\n", __LINE__);
  //     return -1;
  //     }

  //     mpdID = (thesewords & 0x1F0000) >> 16;

  //     int nevent = (thesewords & 0x00FF00) >> 8;

  //     for( ii = 0; ii < nevent; ii++ ){
  //     thesewords = p[jj++] & 0xFFFFFF;

  //     //printf("EVENT HEADER  %06x\n", thesewords);
  //     //printf("Good? (4)       %x\n", (thesewords & 0xF00000) >> 20);
  //     //printf("EVENT COUNT     %d\n", (thesewords & 0x0FFFFF) >> 0);
  //     if( (thesewords & 0xF00000) >> 20 != 0x4 ){
  //     fprintf(stderr, "[ERROR  MPDModule::LoadSlot, line %d] EVENT HEADER NOT FOUND\n", __LINE__);
  //     return -1;
  //     }

  //     thesewords = p[jj++] & 0xFFFFFF;
  //     //printf("TRIGGER TIME 1%06x\n", thesewords);
  //     //printf("Good? (6)       %x\n", (thesewords & 0xF00000) >> 20);
  //     //printf("COURSE TIME     %d\n", (thesewords & 0x0FFFFF) >> 0);
  //     if( (thesewords & 0xF00000) >> 20 != 0x6 ){
  //     fprintf(stderr, "[ERROR  MPDModule::LoadSlot, line %d] TRIGGER TIME 1 WORD NOT FOUND\n", __LINE__);
  //     return -1;
  //     }

  //     thesewords = p[jj++] & 0xFFFFFF;
  //     //printf("TRIGGER TIME 2%06x\n", thesewords);
  //     //printf("Good? (7)       %x\n", (thesewords & 0xF00000) >> 20);
  //     //printf("COURSE TIME     %d\n", (thesewords & 0x0FFFFF) >> 0);
  //     if( (thesewords & 0xF00000) >> 20 != 0x7 ){
  //     fprintf(stderr, "[ERROR  MPDModule::LoadSlot, line %d] TRIGGER TIME 2 WORD NOT FOUND\n", __LINE__);
  //     return -1;
  //     }

  //     kk = 0;
  //     while( ((p[jj] & 0xE00000) >> 21 ) == 0x4  ){
  //     kk++;
  //     //printf("\n[Starting sample %d]\n", kk);

  //     thesewords = p[jj++] & 0x1FFFFF;

  //     adcCh = thesewords & 0x00000F;

  //     //printf("HEADER        %06x\n", thesewords);
  //     //printf("Headergood? (0) %x\n", (thesewords & 0x1C0000) >> 18);
  //     //printf("Baselineval     %x\n", (thesewords & 0x020000) >> 17);
  //     //printf("APV HEADER      %x\n", (thesewords & 0x01FFF0) >> 4);
  //     //printf("APV ID          %x\n", (thesewords & 0x00000F) >> 0);
  //     if( (thesewords & 0x1C0000) >> 18 != 0x0 ){
  //     fprintf(stderr, "[ERROR  MPDModule::LoadSlot, line %d] DATA HEADER NOT FOUND\n", __LINE__);
  //     return -1;
  //     }

  //     // Loop while still seeing reduced data
  //     while( ((p[jj] & 0x180000) >> 19) == 0x1 ){
  //     for( ll = 0; ll < 8; ll++ ){
  //     if( ((p[jj] & 0x180000) >> 19) != 0x1 ){
  //     break;
  //     }
  //     //                      printf("%08x  ", p[jj++]);
  //     int x_data = p[jj++];

  //     data =  x_data& 0x00FFF;
  //     ch   = (x_data& 0x7F000)>>12;
  //     //printf("%3d %03x  ", ch, data);

  //     // Otherwise we have data
  //     effCh = (mpdID) << 8 | adcCh;

  //     status = sldat->loadData("adc",effCh, data, data);
  //     if( status != SD_OK ) return -1;

  //     fWordsSeen++;
  //     data_count++;
  //     }
  //     //printf("\n");
  //     }

  //     thesewords = p[jj++] & 0x1FFFFF;
  //     //printf("APV TRAILER   %06x\n", thesewords);
  //     //printf("Good? (8)       %x\n", (thesewords & 0x1E0000) >> 17);
  //     //printf("Module ID       %x\n", (thesewords & 0x01F000) >> 12);
  //     //printf("Sample Count    %x\n", (thesewords & 0x000F00) >>  8);
  //     //printf("Frame Counter   %x\n", (thesewords & 0x0000FF) >>  0);
  //     if( (thesewords & 0x1E0000) >> 17 != 0x8 ){
  //     fprintf(stderr, "[ERROR  MPDModule::LoadSlot, line %d] APV TRAILER NOT FOUND\n", __LINE__);
  //     return -1;
  //     }

  //     thesewords = p[jj++] & 0x1FFFFF;
  //     //printf("TRAILER       %06x\n", thesewords);
  //     //printf("Good? (3)       %x\n", (thesewords & 0x180000) >> 19);
  //     //printf("Baseline val    %x\n", (thesewords & 0x07FF00) >>  8);
  //     //printf("Word count      %x\n", (thesewords & 0x0000FF) >>  0);
  //     if( (thesewords & 0x180000) >> 19 != 0x3 ){
  //     fprintf(stderr, "[ERROR  MPDModule::LoadSlot, line %d] DATA TRAILER NOT FOUND\n", __LINE__);
  //     return -1;
  //     }

  //     }
  //     //printf("[ %d SAMPLES ]\n", kk);

  //     thesewords = p[jj++] & 0xFFFFFF;
  //     //printf("EVENT TRAILER %06x\n", thesewords);
  //     //printf("Good? (a)       %x\n", (thesewords & 0xF00000) >> 20);
  //     //printf("N WORDS IN EVT  %d\n", (thesewords & 0x0FFF00) >> 8);
  //     //printf("FINE TRIGGER T  %d\n", (thesewords & 0x0000FF) >> 0);
  //     if( (thesewords & 0xF00000) >> 20 != 0xa ){
  //     fprintf(stderr, "[ERROR  MPDModule::LoadSlot, line %d] EVENT TRAILER NOT FOUND\n", __LINE__);
  //     return -1;
  //     }

  //     }

  //     // Filler words to end
  //     thesewords = p[jj++] & 0xFFFFFF;
  //     while( thesewords == 0xe00000 ){
  //     thesewords = p[jj++] & 0xFFFFFF;
  //     }

  //     //printf("BLOCK TRAILER %06x\n", thesewords);
  //     //printf("Good? (2)       %x\n", (thesewords & 0xF00000) >> 20);
  //     //printf("NWORDS IN BLOCK %d\n", (thesewords & 0x0FFFFF) >> 0);
  //     if( (thesewords & 0xF00000) >> 20 != 0x2 ){
  //     fprintf(stderr, "[ERROR  MPDModule::LoadSlot, line %d] BLOCK TRAILER NOT FOUND\n", __LINE__);
  //     return -1;
  //     }
  //     }
  //     */

  //     //printf("=================  END !!! =================================\n");

  //     /*
  // 	v4 decoder

  // 	for( Int_t i = 0; i < len; i++ ){
  // 	tag = p[i] & 0xf0000000;

  // 	switch(tag) {
  // 	case MPD_MODULE_TAG:
  // 	cout <<"Module TAG"<<endl;
  // 	mpdID = p[i] & 0xffff;
  // 	break;
  // 	case MPD_ADC_TAG:
  // 	cout <<"ADC TAG"<<endl;
  // 	adcCh = p[i] & 0xff;
  // 	break;
  // 	case MPD_HEADER_TAG:
  // 	cout <<"HEADER TAG"<<endl;
  // 	header = (p[i] >> 4) & 0x1ff;
  // 	// This is following the decoder I got from Evaristo
  // 	// It doesn't seem to match the data I have
  // 	//if( (header & 0xe00) != 0xe00 ){
  // 	// APV interal memory error in header decoding
  // 	//   fprintf(stderr, "MPDModule::LoadSlot Warning: APV memory corruption 0x%03x\n", header );
  // 	//   return -1;
  // 	//}
  // 	break;
  // 	case MPD_TRAILER_TAG:
  // 	cout <<"TRAILER TAG"<<endl;
  // 	// Not sure if this is useful to save
  // 	trailer = p[i] & 0xfff;
  // 	if( (data_count % 16) != 0 ){
  // 	// Missing data
  // 	fprintf(stderr, "MPDModule::LoadSlot Warning: Missing data?\n");
  // 	return -1;
  // 	}
  // 	data_count = 0;
  // 	break;

  // 	case MPD_DATA_TAG:
  // 	cout <<"DATA TAG"<<endl;
  // 	// Not sure if this is useful to save
  // 	data = p[i] & 0xfff;
  // 	ch   = (p[i] >> 12) & 0x7f;

  // 	// Otherwise we have data
  // 	effCh = (mpdID) << 8 | adcCh;
  // 	if( fDebugFile ){
  // 	*fDebugFile << hex << "raw ev buff   "<< mpdID << "   " << adcCh <<"   "<< p[i]  <<endl;
  // 	}

  // 	status = sldat->loadData("adc",effCh, data, data);
  // 	if( status != SD_OK ) return -1;

  // 	fWordsSeen++;
  // 	data_count++;
  // 	break;

  // 	default:
  // 	// Bad tag
  // 	fprintf(stderr, "MPDModule::LoadSlot Warning: Bad Tag 0x%08x\n", tag);
  // 	return -1;

  // 	}

  // 	}
  //     */

  //   return fWordsSeen;
  // }


  // //Unclear if these are used by anything: comment for now (AJRP)
  UInt_t MPDModule::GetData( UInt_t adc, UInt_t sample, UInt_t chan) const {
    // printf("MPD GET DATA\n");
    // UInt_t idx = asc2i(adc, sample, chan);
    // if (idx >= fNumChan*fNumSample*fNumADC) { return 0; }
    // return fData[idx];
    return 0;
  }
  
  void MPDModule::Clear(const Option_t* opt) {
    VmeModule::Clear(opt);
    // fNumHits = 0;
    // for (Int_t i=0; i<fNumChan*fNumSample*fNumADC; i++) fData[i]=0;
    // for (Int_t i=0; i<fNumADC*fNumSample; i++) { 
    //   fFrameHeader[i]=0;
    //   fFrameTrailer[i]=0;
    // }
    
  }
  
  Int_t MPDModule::Decode(const UInt_t *pdat) {
    //Doesn't do anything. I suppose that's fine for now?
    
    return 0;
  }
}


ClassImp(Decoder::MPDModule)
