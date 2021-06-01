//////////////////////////////////////////////////////////////////
//
//   SBSSimTDC

#include "SBSSimTDC.h"
#include "SBSSimDataDecoder.h"
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

  Module::TypeIter_t SBSSimTDC::fgThisType =
    //DoRegister( ModuleType( "Decoder::SBSSimTDC" , 53204 ));
    DoRegister( ModuleType( "Decoder::SBSSimTDC" , -3204 ));
  
  Module::TypeIter_t SBSSimTDC::fgType1 =
    DoRegister( ModuleType( "Decoder::SBSSimTDC" , -1877 ));
  Module::TypeIter_t SBSSimTDC::fgType2 =
    DoRegister( ModuleType( "Decoder::SBSSimTDC" , -1190 ));
  Module::TypeIter_t SBSSimTDC::fgType3 =
    DoRegister( ModuleType( "Decoder::SBSSimTDC" , -3201 ));

  SBSSimTDC::SBSSimTDC()
  {
    tdc_data.resize(NTDCCHAN);
  }

  SBSSimTDC::SBSSimTDC(Int_t crate, Int_t slot)
    : PipeliningModule(crate, slot)//EPAF: I think we need that don't we ?
  {
    tdc_data.resize(NTDCCHAN);
    IsInit = kFALSE;
    Init();
  }

  SBSSimTDC::~SBSSimTDC() {
#if defined DEBUG && defined WITH_DEBUG
    // delete fDebugFile; fDebugFile = 0;
#endif
  }

  /*
  Bool_t SBSSimTDC::HasCapability(Decoder::EModuleType type) {
    return ( type == kSampleADC || type == kPulseIntegral || type == kPulseTime
        || type == kPulsePeak || type == kPulsePedestal || type == kCoarseTime || type == kFineTime);
  } */

  // Clear all data vectors
  void SBSSimTDC::ClearDataVectors() {
    // Clear all data objects
    assert(tdc_data.size() == NTDCCHAN);  // Initialization error in constructor
    for (uint32_t i = 0; i < NTDCCHAN; i++) {
      tdc_data[i].lead_time.clear();
      tdc_data[i].trail_time.clear();
    }
  }

  void SBSSimTDC::Clear( const Option_t* opt) {
    // Clear event-by-event data
    ClearDataVectors();
  }

  void SBSSimTDC::Init() {
    Clear();
    IsInit = kTRUE;
    fName = "SBSSimTDC Generic TDC Module";
  }

  void SBSSimTDC::CheckDecoderStatus() const {
  }

  UInt_t SBSSimTDC::LoadSlot(THaSlotData *sldat, const UInt_t *evbuffer,
      const UInt_t *pstop) {
 
    std::vector<UInt_t> evb(evbuffer, pstop);

    // Note, methods SplitBuffer, GetNextBlock  are defined in PipeliningModule

    SplitBuffer(evb);
    return LoadThisBlock(sldat, GetNextBlock());
   
  }

  UInt_t SBSSimTDC::LoadSlot(THaSlotData *sldat, const UInt_t *evbuffer,
      Int_t pos, Int_t len) {
    return LoadSlot(sldat,evbuffer+pos,evbuffer+pos+len);
    //return SBSSimTDC::LoadSlot(sldat,evbuffer,len);
  }

  UInt_t SBSSimTDC::LoadNextEvBuffer( THaSlotData *sldat) {
    // Note, GetNextBlock belongs to PipeliningModule
    return LoadThisBlock(sldat, GetNextBlock());
  }

  UInt_t SBSSimTDC::LoadThisBlock( THaSlotData *sldat, const vector<UInt_t> evb) {

    // Fill data structures of this class using the event buffer of one "event".
    // An "event" is defined in the traditional way -- a scattering from a target, etc.

    Clear();
    unsigned int nwords = 0;
    unsigned short chan = 0, type = 0;
    UInt_t raw_buff;
    SimEncoder::tdc_data tmp_tdc_data;
    //std::cout << "SBSSimTDC load crate/slot: " << sldat->getCrate() << "/" << sldat->getSlot() << std::endl;
    UInt_t index = 0;
    for( size_t i = 0; i < evb.size(); i++ ){
      // First, decode the header
      chan = type = nwords = 0;
      SBSSimDataDecoder::DecodeHeader(evb[index++],type,chan,nwords);
      //std::cout << " nwords " << nwords << " chan " << chan << " type " << type << std::endl; 
      SBSSimDataDecoder *enc = SBSSimDataDecoder::GetEncoder(type);
      if(!enc) {
	std::cerr << "Could not find TDC decoder of type: " << type
		  << std::endl;
      } else {
        if(!enc->IsTDC()) {
          std::cerr << "Encoder " << enc->GetName() << " of type " << type
		    << " is not a TDC!" << std::endl;
        } else if ( nwords > 0 ) {
	  enc->DecodeTDC(tmp_tdc_data,evb,nwords);
          //std::cerr << "Got TDC encoder for type: " << type
	  //<< ", name: " << enc->GetName() << std::endl;
	  //cout << "time array size ? " << tmp_tdc_data.time.size() << endl;
          for(size_t i = 0; i < tmp_tdc_data.time.size(); i++ ) {
            raw_buff = tmp_tdc_data.getTime(i);
            if(tmp_tdc_data.getEdge(i)) { // Trail
              tdc_data[chan].lead_time.push_back(raw_buff);
            } else { // Lead
              tdc_data[chan].trail_time.push_back(raw_buff);
            }
	    //cout << " raw_buff " << raw_buff << "; edge ? " <<  tmp_tdc_data.getEdge(i) << endl;
            // TODO: Figure out what to do with the edge information
            // I'd imagine we need to distinguish it somehow!
            sldat->loadData("tdc",chan,raw_buff,raw_buff);
          }
          tmp_tdc_data.time.clear(); // Clear it to prepare for next read
        }
      }
      index+= nwords;//evb += nwords; // Skip ahead in the buffer
    }
    return index;
    
    /*
    Clear();

    UInt_t index = 0;
    for( vsiz_t i = 0; i < evb.size(); i++ )
      DecodeOneWord(evb[index++]);

    LoadTHaSlotDataObj(sldat);

    return index;
    */
  }
  
  //UInt_t SBSSimTDC::DecodeOneWord( UInt_t pdat)
  
  
  
}

ClassImp(Decoder::SBSSimTDC)
