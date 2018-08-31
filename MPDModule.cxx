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

      //  v5 decoder
      
      int ii,jj,kk,ll;
      int thesewords;

      jj =  0;


      while( jj < len ){
          thesewords = p[jj++] & 0xFFFFFF;
          //printf("===============================================================================\n");
          //printf("=    CRATE   %d   ======    SLOT   %d   =======================================\n", fCrate, fSlot);
          //printf("BLOCK HEADER  %06x\n", thesewords);
          //printf("Good? (0)       %x\n", (thesewords & 0xe00000) >> 21);
          //printf("Module ID       %d\n", (thesewords & 0x1F0000) >> 16 );
          //printf("EVENT_PER_BLOCK %d\n", (thesewords & 0x00FF00) >> 8 );
          //printf("BLOCK COUNT     %d\n", (thesewords & 0x0000FF) >> 0);

          if( (thesewords & 0xe00000) >> 21 != 0 ){
              fprintf(stderr, "[ERROR  MPDModule::LoadSlot, line %d] BLOCK HEADER NOT FOUND\n", __LINE__);
              return -1;
          }

          mpdID = (thesewords & 0x1F0000) >> 16;

          int nevent = (thesewords & 0x00FF00) >> 8;

          for( ii = 0; ii < nevent; ii++ ){
              thesewords = p[jj++] & 0xFFFFFF;

              //printf("EVENT HEADER  %06x\n", thesewords);
              //printf("Good? (4)       %x\n", (thesewords & 0xF00000) >> 20);
              //printf("EVENT COUNT     %d\n", (thesewords & 0x0FFFFF) >> 0);
              if( (thesewords & 0xF00000) >> 20 != 0x4 ){
                  fprintf(stderr, "[ERROR  MPDModule::LoadSlot, line %d] EVENT HEADER NOT FOUND\n", __LINE__);
                  return -1;
              }

              thesewords = p[jj++] & 0xFFFFFF;
              //printf("TRIGGER TIME 1%06x\n", thesewords);
              //printf("Good? (6)       %x\n", (thesewords & 0xF00000) >> 20);
              //printf("COURSE TIME     %d\n", (thesewords & 0x0FFFFF) >> 0);
              if( (thesewords & 0xF00000) >> 20 != 0x6 ){
                  fprintf(stderr, "[ERROR  MPDModule::LoadSlot, line %d] TRIGGER TIME 1 WORD NOT FOUND\n", __LINE__);
                  return -1;
              }

              thesewords = p[jj++] & 0xFFFFFF;
              //printf("TRIGGER TIME 2%06x\n", thesewords);
              //printf("Good? (7)       %x\n", (thesewords & 0xF00000) >> 20);
              //printf("COURSE TIME     %d\n", (thesewords & 0x0FFFFF) >> 0);
              if( (thesewords & 0xF00000) >> 20 != 0x7 ){
                  fprintf(stderr, "[ERROR  MPDModule::LoadSlot, line %d] TRIGGER TIME 2 WORD NOT FOUND\n", __LINE__);
                  return -1;
              }

              kk = 0;
              while( ((p[jj] & 0xE00000) >> 21 ) == 0x4  ){
                  kk++;
                  //printf("\n[Starting sample %d]\n", kk);

                  thesewords = p[jj++] & 0x1FFFFF;

                  adcCh = thesewords & 0x00000F;

                  //printf("HEADER        %06x\n", thesewords);
                  //printf("Headergood? (0) %x\n", (thesewords & 0x1C0000) >> 18);
                  //printf("Baselineval     %x\n", (thesewords & 0x020000) >> 17);
                  //printf("APV HEADER      %x\n", (thesewords & 0x01FFF0) >> 4);
                  //printf("APV ID          %x\n", (thesewords & 0x00000F) >> 0);
                  if( (thesewords & 0x1C0000) >> 18 != 0x0 ){
                      fprintf(stderr, "[ERROR  MPDModule::LoadSlot, line %d] DATA HEADER NOT FOUND\n", __LINE__);
                      return -1;
                  }

                  // Loop while still seeing reduced data
                  while( ((p[jj] & 0x180000) >> 19) == 0x1 ){
                      for( ll = 0; ll < 8; ll++ ){
                          if( ((p[jj] & 0x180000) >> 19) != 0x1 ){
                              break;
                          }
                          //                      printf("%08x  ", p[jj++]);
                          int x_data = p[jj++];

                          data =  x_data& 0x00FFF;
                          ch   = (x_data& 0x7F000)>>12;
                          //printf("%3d %03x  ", ch, data);

                          // Otherwise we have data
                          effCh = (mpdID) << 8 | adcCh;

                          status = sldat->loadData("adc",effCh, data, data);
                          if( status != SD_OK ) return -1;

                          fWordsSeen++;
                          data_count++;
                      }
                      //printf("\n");
                  }

                  thesewords = p[jj++] & 0x1FFFFF;
                  //printf("APV TRAILER   %06x\n", thesewords);
                  //printf("Good? (8)       %x\n", (thesewords & 0x1E0000) >> 17);
                  //printf("Module ID       %x\n", (thesewords & 0x01F000) >> 12);
                  //printf("Sample Count    %x\n", (thesewords & 0x000F00) >>  8);
                  //printf("Frame Counter   %x\n", (thesewords & 0x0000FF) >>  0);
                  if( (thesewords & 0x1E0000) >> 17 != 0x8 ){
                      fprintf(stderr, "[ERROR  MPDModule::LoadSlot, line %d] APV TRAILER NOT FOUND\n", __LINE__);
                      return -1;
                  }

                  thesewords = p[jj++] & 0x1FFFFF;
                  //printf("TRAILER       %06x\n", thesewords);
                  //printf("Good? (3)       %x\n", (thesewords & 0x180000) >> 19);
                  //printf("Baseline val    %x\n", (thesewords & 0x07FF00) >>  8);
                  //printf("Word count      %x\n", (thesewords & 0x0000FF) >>  0);
                  if( (thesewords & 0x180000) >> 19 != 0x3 ){
                      fprintf(stderr, "[ERROR  MPDModule::LoadSlot, line %d] DATA TRAILER NOT FOUND\n", __LINE__);
                      return -1;
                  }

              }
              //printf("[ %d SAMPLES ]\n", kk);

              thesewords = p[jj++] & 0xFFFFFF;
              //printf("EVENT TRAILER %06x\n", thesewords);
              //printf("Good? (a)       %x\n", (thesewords & 0xF00000) >> 20);
              //printf("N WORDS IN EVT  %d\n", (thesewords & 0x0FFF00) >> 8);
              //printf("FINE TRIGGER T  %d\n", (thesewords & 0x0000FF) >> 0);
              if( (thesewords & 0xF00000) >> 20 != 0xa ){
                  fprintf(stderr, "[ERROR  MPDModule::LoadSlot, line %d] EVENT TRAILER NOT FOUND\n", __LINE__);
                  return -1;
              }

          }

          // Filler words to end
          thesewords = p[jj++] & 0xFFFFFF;
          while( thesewords == 0xe00000 ){
              thesewords = p[jj++] & 0xFFFFFF;
          }

          //printf("BLOCK TRAILER %06x\n", thesewords);
          //printf("Good? (2)       %x\n", (thesewords & 0xF00000) >> 20);
          //printf("NWORDS IN BLOCK %d\n", (thesewords & 0x0FFFFF) >> 0);
          if( (thesewords & 0xF00000) >> 20 != 0x2 ){
              fprintf(stderr, "[ERROR  MPDModule::LoadSlot, line %d] BLOCK TRAILER NOT FOUND\n", __LINE__);
              return -1;
          }
      }

      //printf("=================  END !!! =================================\n");

      /*
         v4 decoder

      for( Int_t i = 0; i < len; i++ ){
          tag = p[i] & 0xf0000000;

          switch(tag) {
              case MPD_MODULE_TAG:
                  cout <<"Module TAG"<<endl;
                  mpdID = p[i] & 0xffff;
                  break;
              case MPD_ADC_TAG:
                  cout <<"ADC TAG"<<endl;
                  adcCh = p[i] & 0xff;
                  break;
              case MPD_HEADER_TAG:
                  cout <<"HEADER TAG"<<endl;
                  header = (p[i] >> 4) & 0x1ff;
                   // This is following the decoder I got from Evaristo
                   // It doesn't seem to match the data I have
                  //if( (header & 0xe00) != 0xe00 ){
                      // APV interal memory error in header decoding
                   //   fprintf(stderr, "MPDModule::LoadSlot Warning: APV memory corruption 0x%03x\n", header );
                   //   return -1;
                  //}
                  break;
              case MPD_TRAILER_TAG:
                  cout <<"TRAILER TAG"<<endl;
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
                  cout <<"DATA TAG"<<endl;
                  // Not sure if this is useful to save
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
  */

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
