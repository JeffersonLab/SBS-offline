/** \class VETROCcdet Module
    \author Stephen Wood
    \author Simona Malace
    \author Brad Sawatzky
    \author Eric Pooser

    Decoder module to retrieve Caen 1190 TDCs.  Based on CAEN 1190 decoding in
    THaCodaDecoder.C in podd 1.5.   (Written by S. Malace, modified by B. Sawatzky)
*/

#include "VETROCcdet.h"
#include "THaSlotData.h"
#include "jlabdec.h"
#include <iostream>

using namespace std;

//#define DEBUG
//#define WITH_DEBUG

namespace Decoder {

  const UInt_t NTDCCHAN = 256;
  const UInt_t MAXHIT = 100;

  Module::TypeIter_t VETROCcdetModule::fgThisType =
    DoRegister( ModuleType( "Decoder::VETROCcdetModule" , 527 ));

  VETROCcdetModule::VETROCcdetModule(Int_t crate, Int_t slot)
    : VmeModule(crate, slot), fNumHits(NTDCCHAN), fTdcData(NTDCCHAN*MAXHIT),
      fTdcOpt(NTDCCHAN*MAXHIT), slot_data(nullptr)
  {
    VETROCcdetModule::Init();
  }

  void VETROCcdetModule::Init() {
    VmeModule::Init();
    fNumHits.resize(NTDCCHAN);
    fTdcData.resize(NTDCCHAN*MAXHIT);
    fTdcOpt.resize(NTDCCHAN*MAXHIT);
    Clear();
    IsInit = true;
    fName = "VETROCcdet Module";
  }

  UInt_t VETROCcdetModule::LoadSlot( THaSlotData* sldat, const UInt_t* evbuffer,
                                   const UInt_t *pstop) {
    // This is a simple, default method for loading a slot
    const UInt_t *p = evbuffer;
    slot_data = sldat;
    fWordsSeen = 0; 		// Word count including global header  
    Int_t glbl_trl = 0;
    while(p <= pstop && glbl_trl == 0) {
      glbl_trl = Decode(p);
      fWordsSeen++;
      ++p;
    }
    return fWordsSeen;
  }

  UInt_t VETROCcdetModule::LoadSlot( THaSlotData* sldat, const UInt_t* evbuffer,
                                   UInt_t pos, UInt_t len ) {
    // Fill data structures of this class, utilizing bank structure
    // Read until out of data or until decode says that the slot is finished
    // len = ndata in event, pos = word number for block header in event
    slot_data = sldat;
    fWordsSeen = 0;
    while(fWordsSeen < len) {
      Decode(&evbuffer[pos+fWordsSeen]);
      fWordsSeen++;
    }
    return fWordsSeen;
  }

  Int_t VETROCcdetModule::Decode(const UInt_t *p) {
    //std::cout << "In this VETROCcdetModule::Decode()" << std::endl;
    Int_t glbl_trl = 0;
    UInt_t data_type_cont =(*p >> 31) & 0x1; 
    UInt_t data_type_def =(*p >> 27) & 0xF; 
    //std::cout << data_type_cont << " " << data_type_def << std::endl; 
    if (data_type_def ==0 && data_type_cont==0) data_type_def = 3; // trigger time continuatio 
    //std::cout << data_type_cont << " " << data_type_def << std::endl; 
    
    static uint32_t type_last = 15;
    static uint32_t time_last = 0;
    static int new_type = 0;
    int type_current = 0;
    generic_data_word_t gword;

    gword.raw = *p;

   if(gword.bf.data_type_defining) /* data type defining word */
     {
       new_type = 1;
       type_current = gword.bf.data_type_tag;
     }
   else
     {
       new_type = 0;
       type_current = type_last;
     }

    //std::cout << "type_current = " << type_current << std::endl;

    switch(type_current) {
    case 0: // block header
      	tdc_data.glb_hdr_slno = (*p >> 22) & 0x1F; // 
	//	tdc_data.glb_hdr_slno = (*p >> 8) & 0x3FF;       // bits 
	if (tdc_data.glb_hdr_slno == fSlot) {
#ifdef WITH_DEBUG
	if (fDebugFile)
	  *fDebugFile << "VETROCcdetModule:: Block HEADER >> data = " 
		      << hex << *p << " >> event number = " << dec 
		      <<  << " >> slot number = "  
		      << tdc_data.glb_hdr_slno << endl;
#endif
      }
       break;
    case 1: // block trailer
      glbl_trl=1;
       break;
    case 2: // event header
        tdc_data.evh_trig_num = *p & 0x7FFFFF;  // Event header trigger number
	break;
    case 3: // trigger time low 24
      if (tdc_data.glb_hdr_slno == fSlot) {
      if (data_type_cont==1) {
        tdc_data.trig_time_l = *p & 0x7FFFFF;  // Event header trigger time low 24
      } else {
        tdc_data.trig_time_h = *p & 0x7FFFFF;  // Event header trigger time high 24
	tdc_data.trig_time = (tdc_data.trig_time_h << 24) | tdc_data.trig_time_l;
      }
      }
       break;
    case 7: // TDC hit
      //cout << "Slot number = " << tdc_data.glb_hdr_slno << endl;
      
      if (tdc_data.glb_hdr_slno == fSlot) {
	//tdc_data.chan   = (*p & 0x0ff0000)>>16; // bits 23-16
	//tdc_data.raw    =  *p & 0x000ffff;      // bits 15-0
	//tdc_data.opt    = (*p & 0x04000000)>>26;      // bit 26
        Int_t group  = ((*p) & 0x07000000)>>24;
        Int_t chan   = ((*p) & 0x00F80000)>>19;
        Int_t edgeD  = ((*p) & 0x00040000)>>18;
        Int_t coarse = ((*p) & 0x0003FF00)>>8;
        Int_t two_ns = ((*p) & 0x00000080)>>7;
        Int_t fine   = ((*p) & 0x0000007F);

        tdc_data.opt = edgeD;
        tdc_data.chan = (group*32 + chan);
        tdc_data.raw = coarse*4000 + two_ns*2000 + fine*2000/109.59; // this should be the time in ps (from Tritium code)
	
	tdc_data.status = slot_data->loadData("tdc", tdc_data.chan, tdc_data.raw, tdc_data.opt);
#ifdef WITH_DEBUG
	if (fDebugFile)
	  *fDebugFile << "VETROCcdetModule:: MEASURED DATA >> data = " 
		      << hex << *p << " >> channel = " << dec
		      << tdc_data.chan << " >> edge = "
		      << tdc_data.opt  << " >> raw time = "
		      << tdc_data.raw << " >> status = "
		      << tdc_data.status << endl;
#endif
	 //std::cout << "VETROCcdetModule:: MEASURED DATA >> data = " 
	//	      << hex << *p << " >> channel = " << dec
	//	      << tdc_data.chan << " >> slot = " << tdc_data.glb_hdr_slno << " >> edge = "
	//	      << tdc_data.opt  << " >> raw time = "
	//	      << tdc_data.raw << " >> status = "
	//	      << tdc_data.status << std::endl;
        if(tdc_data.chan < NTDCCHAN &&
           fNumHits[tdc_data.chan] < MAXHIT) {
          fTdcData[tdc_data.chan * MAXHIT + fNumHits[tdc_data.chan]] = tdc_data.raw;
          fTdcOpt[tdc_data.chan * MAXHIT + fNumHits[tdc_data.chan]++] = tdc_data.opt;
        }
        if (tdc_data.status != SD_OK ) return -1;
      }
      break;
    case 14 : // invalid data
      break;
    case 15 : // buffer alignment filler word; skip
      break;
    default:	// unknown word
	cout << "unknown word for VETROCcdetModule: " << hex << (*p) << dec << endl;
	cout << "according to global header ev. nr. is: " << " " << tdc_data.glb_hdr_evno << endl;
	break;
    }
    return glbl_trl;  
  }

  UInt_t VETROCcdetModule::GetData( UInt_t chan, UInt_t hit ) const {
    if( hit >= fNumHits[chan] ) return 0;
    UInt_t idx = chan * MAXHIT + hit;
    if( idx > MAXHIT * NTDCCHAN ) return 0;
    return fTdcData[idx];
  }

  UInt_t VETROCcdetModule::GetOpt( UInt_t chan, UInt_t hit ) const {
    if( hit >= fNumHits[chan] ) return 0;
    UInt_t idx = chan * MAXHIT + hit;
    if( idx > MAXHIT * NTDCCHAN ) return 0;
    return fTdcOpt[idx];
  }
  
  void VETROCcdetModule::Clear( Option_t* ) {
    fNumHits.assign(NTDCCHAN, 0);
    fTdcData.assign(NTDCCHAN * MAXHIT, 0);
  }
}

ClassImp(Decoder::VETROCcdetModule)


