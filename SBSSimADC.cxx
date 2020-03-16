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
    DoRegister( ModuleType( "Decoder::SBSSimADC" , 50250 ));

  SBSSimADC::SBSSimADC()
  {
    fadc_data.resize(NADCCHAN);
  }

  SBSSimADC::SBSSimADC(Int_t crate, Int_t slot)
    : PipeliningModule(crate, slot)//EPAF: I think we need that don't we ?
  {
    fadc_data.resize(NADCCHAN);
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
    assert(fadc_data.size() == NADCCHAN);  // Initialization error in constructor
    for (uint32_t i = 0; i < NADCCHAN; i++) {
      fadc_data[i].integral = 0;
      fadc_data[i].samples.clear();
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
  }

  Int_t SBSSimADC::LoadSlot(THaSlotData *sldat, const UInt_t *evbuffer,
      const UInt_t *pstop) {
    Clear();
    unsigned int nwords = 0;
    unsigned short chan = 0, type;
    UInt_t raw_buff;
    bool printed = false;
    bool is_first = true;
    //std::cout << "load crate/slot: " << sldat->getCrate() << "/" << sldat->getSlot() << std::endl;
    while(evbuffer < pstop) {
      // First, decode the header
      SBSSimDataEncoder::DecodeHeader(*evbuffer++,type,chan,nwords);
      SBSSimDataEncoder *enc = SBSSimDataEncoder::GetEncoder(type);
      if(!enc) {
        std::cerr << "Could not find ADC decoder of type: " << type
          << ", is_first: " << is_first << std::endl;
      } else {
        if(!enc->IsADC()) {
          std::cerr << "Encoder " << enc->GetName() << " of type " << type
            << " is not an ADC!" << std::endl;
        } else if ( nwords > 0 ) {
          if(enc->IsFADC()) { // FADC with samples
            SimEncoder::fadc_data tmp_fadc_data;
            enc->DecodeFADC(tmp_fadc_data,evbuffer,nwords);
            for(size_t i = 0; i < tmp_fadc_data.samples.size(); i++) {
              raw_buff = tmp_fadc_data.samples[i];
              fadc_data[chan].samples.push_back(tmp_fadc_data.samples[i]);
              sldat->loadData("adc",chan,raw_buff,raw_buff);
            }
          } else if (enc->IsADC()) { // Integral of ADC
            SimEncoder::adc_data tmp_adc_data;
            enc->DecodeADC(tmp_adc_data,evbuffer,nwords);
            raw_buff = tmp_adc_data.integral;
	    //std::cout << raw_buff << endl;
            fadc_data[chan].integral = raw_buff;
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

  Int_t SBSSimADC::LoadSlot(THaSlotData *sldat, const UInt_t *evbuffer,
      Int_t pos, Int_t len) {
    return LoadSlot(sldat,evbuffer+pos,evbuffer+pos+len);
    //return SBSSimADC::LoadSlot(sldat,evbuffer,len);
  }

}

ClassImp(Decoder::SBSSimADC)
