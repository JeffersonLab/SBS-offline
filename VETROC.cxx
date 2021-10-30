/** \class VETROC Module
    \author Stephen Wood
    \author Simona Malace
    \author Brad Sawatzky
    \author Eric Pooser

    Decoder module to retrieve Caen 1190 TDCs.  Based on CAEN 1190 decoding in
    THaCodaDecoder.C in podd 1.5.   (Written by S. Malace, modified by B. Sawatzky)
*/

#include "VETROC.h"
#include "THaSlotData.h"
#include <iostream>

using namespace std;

//#define DEBUG
//#define WITH_DEBUG

namespace Decoder {

  const UInt_t NTDCCHAN = 128;
  const UInt_t MAXHIT = 100;

  Module::TypeIter_t VETROCModule::fgThisType =
    DoRegister( ModuleType( "Decoder::VETROCModule" , 526 ));

  VETROCModule::VETROCModule(Int_t crate, Int_t slot)
    : VmeModule(crate, slot), fNumHits(NTDCCHAN), fTdcData(NTDCCHAN*MAXHIT),
      fTdcOpt(NTDCCHAN*MAXHIT), slot_data(nullptr)
  {
    VETROCModule::Init();
  }

  void VETROCModule::Init() {
    VmeModule::Init();
    fNumHits.resize(NTDCCHAN);
    fTdcData.resize(NTDCCHAN*MAXHIT);
    fTdcOpt.resize(NTDCCHAN*MAXHIT);
    Clear();
    IsInit = true;
    fName = "VETROC Module";
  }

  UInt_t VETROCModule::LoadSlot( THaSlotData* sldat, const UInt_t* evbuffer,
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

  UInt_t VETROCModule::LoadSlot( THaSlotData* sldat, const UInt_t* evbuffer,
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

  Int_t VETROCModule::Decode(const UInt_t *p) {
    Int_t glbl_trl = 0;
    UInt_t data_type_cont =(*p >> 31) & 0x1; 
    UInt_t data_type_def =(*p >> 27) & 0xF; 
    if (data_type_def ==0 && data_type_cont==0) data_type_def = 3; // trigger time continuatio 
    switch(data_type_def) {
    case 0: // block header
      	tdc_data.glb_hdr_slno = (*p >> 22) & 0x1F; // 
	//	tdc_data.glb_hdr_slno = (*p >> 8) & 0x3FF;       // bits 
	if (tdc_data.glb_hdr_slno == fSlot) {
#ifdef WITH_DEBUG
	if (fDebugFile)
	  *fDebugFile << "VETROCModule:: Block HEADER >> data = " 
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
    case 8: // TDC hit
      if (tdc_data.glb_hdr_slno == fSlot) {
	tdc_data.chan   = (*p & 0x0ff0000)>>16; // bits 23-16
	tdc_data.raw    =  *p & 0x000ffff;      // bits 15-0
	tdc_data.opt    = (*p & 0x04000000)>>26;      // bit 26
	tdc_data.status = slot_data->loadData("tdc", tdc_data.chan, tdc_data.raw, tdc_data.opt);
#ifdef WITH_DEBUG
	if (fDebugFile)
	  *fDebugFile << "VETROCModule:: MEASURED DATA >> data = " 
		      << hex << *p << " >> channel = " << dec
		      << tdc_data.chan << " >> edge = "
		      << tdc_data.opt  << " >> raw time = "
		      << tdc_data.raw << " >> status = "
		      << tdc_data.status << endl;
#endif
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
	cout << "unknown word for VETROC: " << hex << (*p) << dec << endl;
	cout << "according to global header ev. nr. is: " << " " << tdc_data.glb_hdr_evno << endl;
	break;
    }
    return glbl_trl;  
  }

  UInt_t VETROCModule::GetData( UInt_t chan, UInt_t hit ) const {
    if( hit >= fNumHits[chan] ) return 0;
    UInt_t idx = chan * MAXHIT + hit;
    if( idx > MAXHIT * NTDCCHAN ) return 0;
    return fTdcData[idx];
  }

  UInt_t VETROCModule::GetOpt( UInt_t chan, UInt_t hit ) const {
    if( hit >= fNumHits[chan] ) return 0;
    UInt_t idx = chan * MAXHIT + hit;
    if( idx > MAXHIT * NTDCCHAN ) return 0;
    return fTdcOpt[idx];
  }

  void VETROCModule::Clear( Option_t* ) {
    fNumHits.assign(NTDCCHAN, 0);
    fTdcData.assign(NTDCCHAN * MAXHIT, 0);
    fTdcOpt.assign(NTDCCHAN * MAXHIT, 0);
  }
}

ClassImp(Decoder::VETROCModule)


