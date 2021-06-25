//////////////////////////////////////////////////////////////////
//
//   SBSSimADC

#include "SBSSimADC.h"
#include "THaEvData.h"
#include "THaSlotData.h"
#include "TMath.h"

#include <unistd.h>
#include <stdint.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>  // for memset
#include <iostream>
#include <string>
#include <sstream>
#include <cmath>
#include <iomanip>
#include <numeric>
#include <cassert>

using namespace std;

//#define DEBUG
//#define WITH_DEBUG

namespace Decoder {

  Module::TypeIter_t SBSSimADC::fgThisType =
    //DoRegister( ModuleType( "Decoder::SBSSimADC" , 50250 ));
    DoRegister( ModuleType( "Decoder::SBSSimADC" ,  -250 ));
  Module::TypeIter_t SBSSimADC::fgType1 =
    DoRegister( ModuleType( "Decoder::SBSSimADC" , -1881 ));
  Module::TypeIter_t SBSSimADC::fgType2 =
    DoRegister( ModuleType( "Decoder::SBSSimADC" ,  -792 ));
  Module::TypeIter_t SBSSimADC::fgType3 =
    DoRegister( ModuleType( "Decoder::SBSSimADC" , -3561 ));
  //Int_t modid[4] = {-250, -1881, -792, -3561};
  
  SBSSimADC::SBSSimADC()
  {
    sadc_data.resize(NADCCHAN);
  }

  SBSSimADC::SBSSimADC(Int_t crate, Int_t slot)
    : PipeliningModule(crate, slot)//EPAF: I think we need that don't we ?
  {
    sadc_data.resize(NADCCHAN);
    IsInit = kFALSE;
    Init();
  }

  SBSSimADC::~SBSSimADC() {
#if defined DEBUG && defined WITH_DEBUG
    // delete fDebugFile; fDebugFile = 0;
#endif
  }

  /*
  Bool_t SBSSimADC::HasCapability(Decoder::EModuleType type) {
    return ( type == kSampleADC || type == kPulseIntegral || type == kPulseTime
        || type == kPulsePeak || type == kPulsePedestal || type == kCoarseTime || type == kFineTime);
  } */

  // Clear all data vectors
  void SBSSimADC::ClearDataVectors() {
    // Clear all data objects
    assert(sadc_data.size() == NADCCHAN);  // Initialization error in constructor
    for (uint32_t i = 0; i < NADCCHAN; i++) {
      sadc_data[i].integral = 0;
      sadc_data[i].samples.clear();
    }
  }

  void SBSSimADC::Clear( const Option_t* opt) {
    // Clear event-by-event data
    ClearDataVectors();
  }

  void SBSSimADC::Init() {
    Clear();
    IsInit = kTRUE;
    fName = "SBSSimADC (Simple JLab Flash ADC Simulated Module)";
  }

  void SBSSimADC::CheckDecoderStatus() const {
    std::cout << "SBSSimADC has been called" << std::endl;
  }

  UInt_t SBSSimADC::LoadSlot(THaSlotData *sldat, const UInt_t *evbuffer,
      const UInt_t *pstop) {
    Clear();
    unsigned int nwords = 0;
    unsigned short chan = 0, type;
    UInt_t raw_buff, strip;
    bool printed = false;
    bool is_first = true;
    //std::cout << "SBSSimADC load crate/slot: " << sldat->getCrate() << "/" << sldat->getSlot() << std::endl;
    while(evbuffer < pstop) {
      // First, decode the header
      SBSSimDataDecoder::DecodeHeader(*evbuffer++,type,chan,nwords);
      //std::cout << type << " " << sldat->getCrate() << " " << sldat->getSlot() << " " << chan << " " << nwords << endl;
      SBSSimDataDecoder *enc = SBSSimDataDecoder::GetEncoder(type);
      if(!enc) {
        std::cerr << "Could not find ADC decoder of type: " << type
          << ", is_first: " << is_first << std::endl;
      } else {
        if(!enc->IsADC() && !enc->IsMPD()) {
          std::cerr << "Encoder " << enc->GetName() << " of type " << type
            << " is not an ADC nor an MPD!" << std::endl;
        } else if ( nwords > 0 ) {
          //if(enc->IsFADC() || enc->IsMPD()) { // FADC with samples
	  if(enc->IsSADC()) { // FADC with samples
            SimEncoder::sadc_data tmp_sadc_data;
            //enc->DecodeFADC(tmp_sadc_data,evbuffer,nwords);
	    enc->DecodeSADC(tmp_sadc_data,evbuffer,nwords);
	    //std::cout << tmp_sadc_data.samples.size() << std::endl;
	    //std::cout << chan << " " << sadc_data[chan].samples.size() << std::endl;
	    for(size_t i = 0; i < tmp_sadc_data.samples.size(); i++) {
              raw_buff = tmp_sadc_data.samples[i];
	      //std::cout << i << " " << tmp_sadc_data.samples[i] << endl;
              sadc_data[chan].samples.push_back(tmp_sadc_data.samples[i]);
	      //std::cout << i << " " << sadc_data[chan].samples.size() << " " << raw_buff << endl;
              sldat->loadData("adc",chan,raw_buff,raw_buff);
            }
          } else if(enc->IsMPD()){
	    SimEncoder::mpd_data tmp_mpd_data;
	    enc->DecodeMPD(tmp_mpd_data,evbuffer,nwords);
	    for(size_t i = 0; i < tmp_mpd_data.samples.size(); i++) {
	      raw_buff = tmp_mpd_data.samples[i];
	      strip = tmp_mpd_data.strips[i];
	      //std::cout << i << " " << chan << " " << tmp_mpd_data.strips[i] << " " << tmp_mpd_data.samples[i] << endl;
	      sldat->loadData("adc",chan,raw_buff,strip);
	    }
	  } else if (enc->IsADC()) { // Integral of ADC
            SimEncoder::adc_data tmp_adc_data;
            enc->DecodeADC(tmp_adc_data,evbuffer,nwords);
            raw_buff = tmp_adc_data.integral;
	    //std::cout << " raw_buff " << raw_buff << endl;
            sadc_data[chan].integral = raw_buff;
            sldat->loadData("adc",chan,raw_buff,raw_buff);
          }
        }
      }
      evbuffer += nwords; // Skip ahead the number of words processed
      is_first = false;
    }
    if(printed)
      std::cerr << std::endl;
   return 0;
  }

  UInt_t SBSSimADC::LoadSlot( THaSlotData *sldat, const UInt_t *evbuffer,
                              UInt_t pos, UInt_t len) {
    return LoadSlot(sldat,evbuffer+pos,evbuffer+pos+len);
    //return SBSSimADC::LoadSlot(sldat,evbuffer,len);
  }

}

ClassImp(Decoder::SBSSimADC)
