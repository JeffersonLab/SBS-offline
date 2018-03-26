//////////////////////////////////////////////////////////////////
//
//   SBSSimFadc250Module
//   JLab FADC 250 Module
//
//   Author: Eric Pooser, pooser@jlab.org
//
//   For documentation pertaining to the FADC250
//   modules refer to the following link:
//   https://coda.jlab.org/drupal/system/files/pdfs/HardwareManual/fADC250/FADC250_Processing_FPGA_Firmware_ver_0x0C01_Description_Instructions.pdf
//
//   The following data type summary corresponds to firmware version 0x0C00 (06/09/2016)
//   -----------------------------------------------------------------------------------
//   0 block header
//   1 block trailer
//   2 event header
//   3 trigger time
//   4 window raw data
//   5 (reserved)
//   6 pulse raw data -> Old firmware
//   7 pulse integral -> Old firmware
//   8 pulse time     -> Old firmware
//   9 pulse parameters
//   10 pulse parameters (pedestal) -> Old firmware
//   11 (reserved)
//   12 scaler data
//   13 (reserved)
//   14 data not valid (empty module)
//   15 filler (non-data) word
//
//   Mode Summary
//   ----------------------------------------------
//   Mode 1:  data type  4          -> Old firmware
//   Mode 2:  data type  6          -> Old firmware
//   Mode 3:  data types 7 & 8      -> Old firmware
//   Mode 4:  data types 8 & 10     -> Old firmware
//   Mode 7:  data types 7 & 8 & 10 -> Old firmware
//   Mode 8:  data types 4 & 8 & 10 -> Old firmware
//   Mode 9:  data type  9
//   Mode 10: data types 4 & 9
//
/////////////////////////////////////////////////////////////////////

#include "SBSSimFadc250Module.h"
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

  Module::TypeIter_t SBSSimFadc250Module::fgThisType =
    DoRegister( ModuleType( "Decoder::SBSSimFadc250Module" , 25099 ));

  SBSSimFadc250Module::SBSSimFadc250Module()
  { fadc_data.resize(NADCCHAN);
    std::cerr << "In constructor!!!" << std::endl; }

  SBSSimFadc250Module::SBSSimFadc250Module(Int_t crate, Int_t slot)
  {
    fadc_data.resize(NADCCHAN);
    IsInit = kFALSE;
    Init();
    std::cerr << "In second constructor!!!" << std::endl;
  }

  SBSSimFadc250Module::~SBSSimFadc250Module() {
#if defined DEBUG && defined WITH_DEBUG
    // delete fDebugFile; fDebugFile = 0;
#endif
  }

  Bool_t SBSSimFadc250Module::HasCapability(Decoder::EModuleType type) {
    return ( type == kSampleADC || type == kPulseIntegral || type == kPulseTime
        || type == kPulsePeak || type == kPulsePedestal || type == kCoarseTime || type == kFineTime);
  }

  // Clear all data vectors
  void SBSSimFadc250Module::ClearDataVectors() {
    // Clear all data objects
    assert(fadc_data.size() == NADCCHAN);  // Initialization error in constructor
    for (uint32_t i = 0; i < NADCCHAN; i++) {
      fadc_data[i].clear();
    }
  }

  void SBSSimFadc250Module::Clear( const Option_t* opt) {
    // Clear event-by-event data
    ClearDataVectors();
  }

  void SBSSimFadc250Module::Init() {
    Clear();
    IsInit = kTRUE;
    fName = "SBSSimFADC250 JLab Flash ADC Module";
  }

  void SBSSimFadc250Module::CheckDecoderStatus() const {
    cout << "FADC250 Decoder has been called" << endl;
  }

  Int_t SBSSimFadc250Module::GetNumFadcEvents(Int_t chan) const {
    return 1;
  }

  Int_t SBSSimFadc250Module::GetNumFadcSamples(Int_t chan, Int_t ievent) const {
    if(chan >= int(NADCCHAN))
      return 0;
    return fadc_data[chan].samples.size();
  }

  void SBSSimFadc250Module::SetPulseSamplesVector(Int_t chan,
      std::vector<uint32_t> pulses) {
    if( chan >> int(NADCCHAN) )
      return;
    for(size_t i = 0; i < pulses.size(); i++) {
      fadc_data[chan].samples.push_back(pulses[i]);
    }
  }

  std::vector<uint32_t> SBSSimFadc250Module::GetPulseSamplesVector(
      Int_t chan) const {
    if( chan >= int(NADCCHAN) )
      return std::vector<uint32_t>();
    return fadc_data[chan].samples;
  }

  Int_t SBSSimFadc250Module::GetPulseIntegralData(Int_t chan,
      Int_t ievent) const
  {
    if( chan >= int(NADCCHAN) || ievent >=
        int(fadc_data[chan].integrals.size()))
      return 0;
    return fadc_data[chan].integrals[ievent];
  }

  Int_t SBSSimFadc250Module::LoadSlot(THaSlotData *sldat, const UInt_t *evbuffer,
      const UInt_t *pstop) {
    Clear();
    int chan = 0, type = 0, num_samples = 0;
    UInt_t raw_buff;
    bool printed = false;
    while(evbuffer < pstop) {
      // First get channel number
      chan = *evbuffer++;
      type = *evbuffer++;
      if(type == 0) { // Samples mode
        num_samples = *evbuffer++;
        for(int i = 0; i < num_samples; i++) {
          raw_buff = *evbuffer++;
          fadc_data[chan].samples.push_back(raw_buff);
          sldat->loadData("adc",chan,raw_buff,raw_buff);
        }
      } else if (type==1) { // integral of adc
        num_samples = *evbuffer++;
        for(int i = 0; i < num_samples; i++) {
          raw_buff = *evbuffer++;
          std::cerr << " [" << chan << ", " << raw_buff << "]";
          printed = true;
          fadc_data[chan].integrals.push_back(raw_buff);
          sldat->loadData("adc",chan,raw_buff,raw_buff);
        }
      }
    }
    if(printed)
      std::cerr << std::endl;
   return 0;
  }

  Int_t SBSSimFadc250Module::LoadSlot(THaSlotData *sldat, const UInt_t *evbuffer,
      Int_t pos, Int_t len) {
    return LoadSlot(sldat,evbuffer+pos,evbuffer+pos+len);
    //return SBSSimFadc250Module::LoadSlot(sldat,evbuffer,len);
/*    int chan=evbuffer[0];
    for (int i = 1; i <= len; i++) {
      fadc_data[chan].samples.push_back(evbuffer[i]);
      sldat->loadData("adc",chan,evbuffer[i],evbuffer[i]);
    }
    return 0;
    */
  }

}

ClassImp(Decoder::SBSSimFadc250Module)
