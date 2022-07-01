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

#include "MPDModuleVMEv4.h"
#include "THaSlotData.h"
#include <limits>
#include <vector>
#include <map>
#include <set>
#include <iostream>

using namespace std;

namespace Decoder {

  Module::TypeIter_t MPDModuleVMEv4::fgThisType =
    DoRegister( ModuleType( "Decoder::MPDModuleVMEv4" , 3560 ));

  MPDModuleVMEv4::MPDModuleVMEv4(Int_t crate, Int_t slot) : VmeModule(crate, slot) {
    fDebugFile=0;
    Init(); //Should this be called here? not clear...
    //fOnlineZeroSuppression = false; //If this is false, then we want to calculate and subtract the common-mode from each ADC sample:
    
  }
  
  MPDModuleVMEv4::~MPDModuleVMEv4() {
    
  }
  
  void MPDModuleVMEv4::Init() { 
    Module::Init();
    //    Config(0,25,6,16,128); // should be called by the user (but how?)
    fBlockHeader = 0x0;
    fAPVHeader   = 0x4;

    fNumSample = 6;
    
    fDebugFile=0;
    Clear();
    //    fName = "MPD Module (INFN MPD for GEM and more), use Config to dynamic config";
    fName = "MPD Module";
  }

  //This version ASSUMES that there is no online zero suppression, so all 128 APV channels are present in every event!
  //We may also need to code up something to handle the alternative case
  
  UInt_t MPDModuleVMEv4::LoadSlot( THaSlotData *sldat, const UInt_t *evbuffer, UInt_t pos, UInt_t len ){

    //std::cout << "Calling MPDModuleVMEv4::LoadSlot" << std::endl;
    
    //AJRP: LoadSlot method for the VME MPD4 data format used by the UVA GEM cosmic test stand ca. Jan. 2021
    const UInt_t *datawords = &(evbuffer[pos]);

    fWordsSeen = 0;

    Int_t status;
    
    UInt_t thisword;

    UInt_t iword=0;

    //Get the slot number for this call to LoadSlot:
    UInt_t this_slot = sldat->getSlot();
    
    ///bool foundslot = false;

    UInt_t thisheader;

    //UInt_t mask, shift;
    UInt_t slot=0, adc_chan=0;
    UInt_t prev_slot=0;

    bool found_adc = false;
    bool found_this_slot = false;
    
    std::map<UInt_t, std::vector<UInt_t> > RawDataByADC_Channel; //represents the raw data in ONE slot:
    
    //following the MPDRawParser in ROOT_GUI_multicrate, loop on all the data in the ROC bank (which corresponds to one "crate"), and populate the
    //temporary data structure above, ONLY if slot == this_slot
    while( iword < len ){
      thisword = datawords[iword++];

      
      
      //Extract word header from bits 22-24 of data word:
      thisheader = (thisword & 0x00E00000)>>21;

      // std::cout << "iword, thisword, thisheader, fBlockHeader, fAPVHeader, fNumSample = " << iword << ", " << thisword << ", "
      // 		<< thisheader << ", " << fBlockHeader << ", " << fAPVHeader << ", "
      // 		<< fNumSample << std::endl;
      
      //Check if new slot:
      if(thisheader == fBlockHeader){ //extract "MPDID" (slot) info from bits 17-21 of data word (NOTE that in this raw data format, "slot" and "MPD_ID" are the same thing

	prev_slot = slot;
	
	slot = (thisword & 0x001F0000)>>16;

	//	std::cout << "Found block header, slot = " << slot << std::endl;
	
	if( slot == this_slot ) found_this_slot = true; //first time we find the desired slot, set found_this_slot to true
	//new_slot = true; //Every time we encounter a block header word, we set new_slot to true.
	found_adc = false; //initialize found_adc to false when we encounter a new block header word
      }

      if( prev_slot == this_slot ) break; //we finished loading the data from the slot we actually want!
    
      if( slot == this_slot && thisheader == fAPVHeader ){ //APV data:

	UInt_t thistype = (thisword & 0x00180000)>>19; //Data type is found in bits 20-21 of data word
	
	if( thistype == 0 ){ //apv header: this word contains the ADC channel info: 
	  adc_chan = (thisword &0xf); //first four bits of data word
	  found_adc = true;

	  //  std::cout << "Found APV header, ADC channel = " << adc_chan << std::endl;
	  
	}

	if( thistype == 1 && found_adc ){ //ADC sample data:
	  RawDataByADC_Channel[adc_chan].push_back( thisword & 0x00000fff ); //Raw ADC samples are contained in bits 1-12 of data word
	  //std::cout << "Found raw data, (channel, raw ADC value) = (" << adc_chan << ", " << RawDataByADC_Channel[adc_chan].back() << ")" << std::endl;
	}
	//We ignore APV trailer and "trailer" (types "2" and "3"). We hope they aren't important
      }
    }
      

    if( found_this_slot ){ //then the current slot has data, extract it and decode it. 
      for( auto iadc = RawDataByADC_Channel.begin(); iadc != RawDataByADC_Channel.end(); ++iadc ){
	//iadc is a pointer to pair<UInt_t, std::vector<UInt_t> > (I think)
	adc_chan = iadc->first; //one APV card
	std::vector<UInt_t> ADCsamples = iadc->second; //vector of all the ADC samples:

	//Since we ignore the APV trailer word, the size of the raw data per APV card should always equal exactly 128 * number of samples:
	if( ADCsamples.size() != fNumSample*128 ) return 0;

	//In the MPD4 VME format that this method is decoding, without online zero suppression, the ADC samples are assumed to be ordered as:
	// index = ichan + 128*isamp

	//With this data version, mpdID and slot are always treated as the same thing:
	UInt_t effChan = (this_slot << 4 | adc_chan );

	for( UInt_t iAPVchan=0; iAPVchan<128; iAPVchan++ ){
	  for( UInt_t iSample=0; iSample<fNumSample; iSample++ ){

	    //SBSGEMModule::Decode() expects the raw data (whether or not it's zero-suppressed) to be ordered such that:
	    // index in raw hit array = isamp + 6*ichan (opposite of the ordering in the event buffer)
	    // The APV channel number is stored in "rawdata", while the raw ADC values are stored in "data",
	    // as the SBSGEMModule::Decode method expects!
	    status = sldat->loadData( "adc", effChan, ADCsamples[iAPVchan+128*iSample], iAPVchan );

	    if( status != SD_OK ) return -1;
	    
	    fWordsSeen++;
	  }
	}
	  
	//We'll handle common-mode subtraction in the SBSGEMModule::Decode method
     
      }
	    
	
    }
    return fWordsSeen;
  }
  
  //
  // UInt_t MPDModuleVMEv4::LoadSlot( THaSlotData *sldat, const UInt_t* evbuffer, UInt_t pos, UInt_t len) {
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
  // 	fprintf(stderr, "[ERROR  MPDModuleVMEv4::LoadSlot, line %d] "
  // 		"BLOCK HEADER NOT FOUND\n", __LINE__);
  // 	return -1;
  //     }

  //     UInt_t nevent = (thesewords&0xFF);
  //     // Ensure we have enough data
  //     // (need at least 4 per event: 1 event header + 2 trigger words +
  //     // 1 event trailer
  //     if( nevent > 0 && jj + (nevent*4) >= len) {
  // 	fprintf(stderr, "[ERROR  MPDModuleVMEv4::LoadSlot, line %d] "
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
  // 	  fprintf(stderr, "[ERROR  MPDModuleVMEv4::LoadSlot, line %d] EVENT HEADER NOT FOUND\n", __LINE__);
  // 	  return -1;
  // 	}

  // 	thesewords = p[jj++] & 0xFFFFFFFF;
  // 	//printf("TRIGGER TIME 1     0x%08x\n", thesewords);
  // 	//printf("Data defining? (1) %d\n", (thesewords & 0x80000000) >> 31);
  // 	//printf("Type (3)           %d\n", (thesewords & 0x78000000) >> 27);
  // 	//printf("COURSE TIME        %d\n", (thesewords & 0xFFFFFF) >> 0);
  // 	//printf("\n");
  // 	if( (SSP_DATADEF(thesewords) != 1) || (SSP_TAG(thesewords) != 3 )) {
  // 	  fprintf(stderr, "[ERROR  MPDModuleVMEv4::LoadSlot, line %d] TRIGGER TIME 1 WORD NOT FOUND\n", __LINE__);
  // 	  return -1;
  // 	}

  // 	thesewords = p[jj++] & 0xFFFFFFFF;
  // 	//printf("TRIGGER TIME 2     %08x\n", thesewords);
  // 	//printf("Data defining? (0) %d\n", (thesewords & 0x80000000) >> 31);
  // 	//printf("COURSE TIME        %d\n", (thesewords & 0xFFFFFF) >> 0);
  // 	//printf("\n");
  // 	if( (SSP_DATADEF(thesewords) != 0) ) {
  // 	  fprintf(stderr, "[ERROR  MPDModuleVMEv4::LoadSlot, line %d] TRIGGER TIME 2 WORD NOT FOUND\n", __LINE__);
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
  // 	      fprintf(stderr, "[ERROR  MPDModuleVMEv4::LoadSlot, line %d] NOT ENOUGH"
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
  // 		fprintf(stderr, "[ERROR  MPDModuleVMEv4::LoadSlot, line %d] MISSING"
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
  // 	fprintf(stderr, "[ERROR  MPDModuleVMEv4::LoadSlot, line %d] "
  // 		"BLOCK TRAILER NOT FOUND\n", __LINE__);
  // 	return -1;
  //     }

  //     //printf("Read number of words %d expected %d\n",jj-mm,thesewords&0x3FFFFF);
  //     if((thesewords&0x3FFFFF) != (jj-mm) ) {
  // 	fprintf(stderr, "[ERROR  MPDModuleVMEv4::LoadSlot, line %d] NUMBER OF "
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
  //     fprintf(stderr, "[ERROR  MPDModuleVMEv4::LoadSlot, line %d] BLOCK HEADER NOT FOUND\n", __LINE__);
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
  //     fprintf(stderr, "[ERROR  MPDModuleVMEv4::LoadSlot, line %d] EVENT HEADER NOT FOUND\n", __LINE__);
  //     return -1;
  //     }

  //     thesewords = p[jj++] & 0xFFFFFF;
  //     //printf("TRIGGER TIME 1%06x\n", thesewords);
  //     //printf("Good? (6)       %x\n", (thesewords & 0xF00000) >> 20);
  //     //printf("COURSE TIME     %d\n", (thesewords & 0x0FFFFF) >> 0);
  //     if( (thesewords & 0xF00000) >> 20 != 0x6 ){
  //     fprintf(stderr, "[ERROR  MPDModuleVMEv4::LoadSlot, line %d] TRIGGER TIME 1 WORD NOT FOUND\n", __LINE__);
  //     return -1;
  //     }

  //     thesewords = p[jj++] & 0xFFFFFF;
  //     //printf("TRIGGER TIME 2%06x\n", thesewords);
  //     //printf("Good? (7)       %x\n", (thesewords & 0xF00000) >> 20);
  //     //printf("COURSE TIME     %d\n", (thesewords & 0x0FFFFF) >> 0);
  //     if( (thesewords & 0xF00000) >> 20 != 0x7 ){
  //     fprintf(stderr, "[ERROR  MPDModuleVMEv4::LoadSlot, line %d] TRIGGER TIME 2 WORD NOT FOUND\n", __LINE__);
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
  //     fprintf(stderr, "[ERROR  MPDModuleVMEv4::LoadSlot, line %d] DATA HEADER NOT FOUND\n", __LINE__);
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
  //     fprintf(stderr, "[ERROR  MPDModuleVMEv4::LoadSlot, line %d] APV TRAILER NOT FOUND\n", __LINE__);
  //     return -1;
  //     }

  //     thesewords = p[jj++] & 0x1FFFFF;
  //     //printf("TRAILER       %06x\n", thesewords);
  //     //printf("Good? (3)       %x\n", (thesewords & 0x180000) >> 19);
  //     //printf("Baseline val    %x\n", (thesewords & 0x07FF00) >>  8);
  //     //printf("Word count      %x\n", (thesewords & 0x0000FF) >>  0);
  //     if( (thesewords & 0x180000) >> 19 != 0x3 ){
  //     fprintf(stderr, "[ERROR  MPDModuleVMEv4::LoadSlot, line %d] DATA TRAILER NOT FOUND\n", __LINE__);
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
  //     fprintf(stderr, "[ERROR  MPDModuleVMEv4::LoadSlot, line %d] EVENT TRAILER NOT FOUND\n", __LINE__);
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
  //     fprintf(stderr, "[ERROR  MPDModuleVMEv4::LoadSlot, line %d] BLOCK TRAILER NOT FOUND\n", __LINE__);
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
  // 	//   fprintf(stderr, "MPDModuleVMEv4::LoadSlot Warning: APV memory corruption 0x%03x\n", header );
  // 	//   return -1;
  // 	//}
  // 	break;
  // 	case MPD_TRAILER_TAG:
  // 	cout <<"TRAILER TAG"<<endl;
  // 	// Not sure if this is useful to save
  // 	trailer = p[i] & 0xfff;
  // 	if( (data_count % 16) != 0 ){
  // 	// Missing data
  // 	fprintf(stderr, "MPDModuleVMEv4::LoadSlot Warning: Missing data?\n");
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
  // 	fprintf(stderr, "MPDModuleVMEv4::LoadSlot Warning: Bad Tag 0x%08x\n", tag);
  // 	return -1;

  // 	}

  // 	}
  //     */

  //   return fWordsSeen;
  // }


  // //Unclear if these are used by anything: comment for now (AJRP)
  UInt_t MPDModuleVMEv4::GetData( UInt_t adc, UInt_t sample, UInt_t chan) const {
    // printf("MPD GET DATA\n");
    // UInt_t idx = asc2i(adc, sample, chan);
    // if (idx >= fNumChan*fNumSample*fNumADC) { return 0; }
    // return fData[idx];
    return 0;
  }
  
  void MPDModuleVMEv4::Clear(const Option_t *opt) {
    VmeModule::Clear(opt);
    // fNumHits = 0;
    // for (Int_t i=0; i<fNumChan*fNumSample*fNumADC; i++) fData[i]=0;
    // for (Int_t i=0; i<fNumADC*fNumSample; i++) { 
    //   fFrameHeader[i]=0;
    //   fFrameTrailer[i]=0;
    // }
    
  }
  
  Int_t MPDModuleVMEv4::Decode(const UInt_t *pdat) {
    //Doesn't do anything. I suppose that's fine for now?
    
    return 0;
  }


}

ClassImp(Decoder::MPDModuleVMEv4)
