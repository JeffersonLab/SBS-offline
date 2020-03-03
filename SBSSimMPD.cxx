//////////////////////////////////////////////////////////////////
//
//   SBSSimMPD

#include "SBSSimMPD.h"
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

  Module::TypeIter_t SBSSimMPD::fgThisType =
    DoRegister( ModuleType( "Decoder::SBSSimMPD" , 53561 ));

  SBSSimMPD::SBSSimMPD()
  {
    //fadc_data.resize(NADCCHAN);
  }

  SBSSimMPD::SBSSimMPD(Int_t crate, Int_t slot)
    : PipeliningModule(crate, slot)//EPAF: I think we need that don't we ?
  {
    //fadc_data.resize(NADCCHAN);
    IsInit = kFALSE;
    Init();
  }

  SBSSimMPD::~SBSSimMPD() {
#if defined DEBUG && defined WITH_DEBUG
    // delete fDebugFile; fDebugFile = 0;
#endif
  }

  /*
  Bool_t SBSSimMPD::HasCapability(Decoder::EModuleType type) {
    return ( type == kSampleADC || type == kPulseIntegral || type == kPulseTime
        || type == kPulsePeak || type == kPulsePedestal || type == kCoarseTime || type == kFineTime);
  } */

  // Clear all data vectors
  void SBSSimMPD::ClearDataVectors() {
    // Clear all data objects
    //assert(fadc_data.size() == NADCCHAN);  // Initialization error in constructor
    //for (uint32_t i = 0; i < NADCCHAN; i++) {
    //  fadc_data[i].integral = 0;
    //  fadc_data[i].samples.clear();
    //}
  }

  void SBSSimMPD::Clear( const Option_t* opt) {
    // Clear event-by-event data
    ClearDataVectors();
  }

  void SBSSimMPD::Init() {
    Clear();
    IsInit = kTRUE;
    fName = "SBSSimMPD (Simple JLab Flash ADC Simulated Module)";
  }

  void SBSSimMPD::CheckDecoderStatus() const {
  }

  Int_t SBSSimMPD::LoadSlot(THaSlotData *sldat, const UInt_t *evbuffer,
      const UInt_t *pstop) {
    Clear();
    unsigned int nwords = 0;
    unsigned short chan = 0, type;
    UInt_t raw_buff;
    Int_t effCh;
    bool printed = false;
    bool is_first = true;
    SimEncoder::mpd_data mpd_data;
    std::vector<unsigned int> tmp_data;
    // Temporary variables to calculate where each strip should
    // be in in the data array
    Int_t RstripNb = 0;
    Int_t isamp;
    while(evbuffer < pstop) {
      // First, decode the header
      SBSSimDataEncoder::DecodeHeader(*evbuffer++,type,chan,nwords);
      SBSSimDataEncoder *enc = SBSSimDataEncoder::GetEncoder(type);
      if(!enc) {
        std::cerr << "Could not find ADC decoder of type: " << type
          << ", is_first: " << is_first << std::endl;
      } else {
        if(!enc->IsMPD()) {
          std::cerr << "Encoder " << enc->GetName() << " of type " << type
            << " is not an MPD!" << std::endl;
        } else if ( nwords > 0 ) {
          enc->DecodeMPD(mpd_data, evbuffer, nwords);
          effCh = (mpd_data.mpd_id<<8)|mpd_data.adc_id;
          // Well, because of this weird multiplexer thing and the
          // fact that SBS-offline expects the data this way, we are
          // going to have to re-arrange the data so that it maps
          // the way the data is read out by the multiplexer.
          // I really don't think this will ever work with zero supression....
          // But fine! It's a first try anyways
          tmp_data.clear();
          tmp_data.resize(mpd_data.samples.size());
          for(int strip = 0; strip < mpd_data.nstrips; strip++) {
            // This part comes from the APV25 multiplexer
            RstripNb= 32*(strip%4)+ 8*(int)(strip/4)- 31*(int)(strip/16);
            // Not sure which other multiplexer is responsible for this
            RstripNb=RstripNb+1+RstripNb%4-5*(((int)(RstripNb/4))%2);
            RstripNb=RstripNb+(127-2*RstripNb)*mpd_data.invert;
            //RstripPos = RstripNb + 128*mpd_data.pos;
            for(int adc_samp = 0; adc_samp < mpd_data.nsamples; adc_samp++) {
              isamp = adc_samp*mpd_data.nstrips + strip;
              tmp_data[isamp] = mpd_data.samples[mpd_data.nsamples*RstripNb
                +adc_samp];
            }
          }
          for(size_t i = 0; i < tmp_data.size(); i++) {
            raw_buff = tmp_data[i];
            if(raw_buff>4095) {
              std::cerr << "Yikes!! What's up with raw_buff? i: "
               << i << ", " << raw_buff
                << std::endl;
            }
            sldat->loadData("adc",effCh,raw_buff,raw_buff);
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

  Int_t SBSSimMPD::LoadSlot(THaSlotData *sldat, const UInt_t *evbuffer,
      Int_t pos, Int_t len) {
    return LoadSlot(sldat,evbuffer+pos,evbuffer+pos+len);
    //return SBSSimMPD::LoadSlot(sldat,evbuffer,len);
  }

}

ClassImp(Decoder::SBSSimMPD)
