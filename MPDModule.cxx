/////////////////////////////////////////////////////////////////////
//
//   MPDModule
//   This is the MPD module decoder; based on SkeletonModule
//   (https://github.com/JeffersonLab/analyzer)
//
//   E. Cisbani
//   Original Version:   2015/Dec
//
/////////////////////////////////////////////////////////////////////

#define MPD_VERSION_TAG 0xE0000000
#define MPD_EVENT_TAG   0x10000000
#define MPD_MODULE_TAG  0x20000000
#define MPD_ADC_TAG     0x30000000
#define MPD_HEADER_TAG  0x40000000
#define MPD_DATA_TAG    0x0
#define MPD_TRAILER_TAG 0x50000000

#include "MPDModule.h"
#include "THaSlotData.h"

using namespace std;

namespace Decoder {

  Module::TypeIter_t MPDModule::fgThisType =
    DoRegister( ModuleType( "Decoder::MPDModule" , 3561 ));

  MPDModule::MPDModule(Int_t crate, Int_t slot) : VmeModule(crate, slot) {
    fDebugFile=0;
    Init();
  }
  
  MPDModule::~MPDModule() {
    
  }
  
  void MPDModule::Init() { 
    Module::Init();
//    Config(0,25,6,16,128); // should be called by the user
    fDebugFile=0;
    Clear("");
//    fName = "MPD Module (INFN MPD for GEM and more), use Config to dynamic config";
    fName = "MPD Module";
  }

  
  Int_t MPDModule::LoadSlot(THaSlotData *sldat, const UInt_t* evbuffer, Int_t pos, Int_t len) {
      const UInt_t *p = &evbuffer[pos];
      UInt_t data;
      fWordsSeen = 0;

      // From stand alone decoder
      // We declare an effective channel from the MPD ID 
      // and ADC channel
      Int_t tag, header, trailer, ch, status;
      Int_t mpdID = -1;
      Int_t adcCh = -1;
      Int_t effCh = 0;

      UInt_t data_count = 0;

      for( Int_t i = 0; i < len; i++ ){
          tag = p[i] & 0xf0000000;

          switch(tag) {
              case MPD_MODULE_TAG:
                  mpdID = p[i] & 0xffff;
                  break;
              case MPD_ADC_TAG:
                  adcCh = p[i] & 0xff;
                  break;
              case MPD_HEADER_TAG:
                  header = (p[i] >> 4) & 0x1ff;
                  /*
                   // This is following the decoder I got from Evaristo
                   // It doesn't seem to match the data I have
                  if( (header & 0xe00) != 0xe00 ){
                      // APV interal memory error in header decoding
                      fprintf(stderr, "MPDModule::LoadSlot Warning: APV memory corruption 0x%03x\n", header );
                      return -1;
                  }
                  */
                  break;
              case MPD_TRAILER_TAG:
                  // Not sure if this is useful to save
                  trailer = p[i] & 0xfff;
                  if( (data_count % 16) != 0 ){
                      // Missing data
                      fprintf(stderr, "MPDModule::LoadSlot Warning: Missing data?\n");
                      return -1;
                  }
                  data_count = 0;
                  break;

              case MPD_DATA_TAG:
                  data = p[i] & 0xfff;
                  ch   = (p[i] >> 12) & 0x7f;

                  // Otherwise we have data
                  effCh = (mpdID) << 8 | adcCh;
                  if( fDebugFile ){
                      *fDebugFile << hex << "raw ev buff   "<< mpdID << "   " << adcCh <<"   "<< p[i]  <<endl;
                  }
                  status = sldat->loadData("adc",effCh, data, data);
                  if( status != SD_OK ) return -1;

                  fWordsSeen++;
                  data_count++;
                  break;

              default:
                  // Bad tag
                  fprintf(stderr, "MPDModule::LoadSlot Warning: Bad Tag 0x%08x\n", tag);
                  return -1;

          }

      }

      return fWordsSeen;
  }
  
  Int_t MPDModule::GetData(Int_t adc, Int_t sample, Int_t chan) const {
    printf("MPD GET DATA\n");
    Int_t idx = asc2i(adc, sample, chan);
    if ((idx < 0 ) || (idx >= fNumChan*fNumSample*fNumADC)) { return 0; }
    return fData[idx];
  }
  
  void MPDModule::Clear(const Option_t *opt) {
    fNumHits = 0;
    for (Int_t i=0; i<fNumChan*fNumSample*fNumADC; i++) fData[i]=0;
    for (Int_t i=0; i<fNumADC*fNumSample; i++) { 
      fFrameHeader[i]=0;
      fFrameTrailer[i]=0;
    }
    
  }
  
  Int_t MPDModule::Decode(const UInt_t *pdat) {
    
    
    return 0;
  }


}

ClassImp(Decoder::MPDModule)
