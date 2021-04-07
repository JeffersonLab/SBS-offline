/////////////////////////////////////////////////////////////////////
//
//   HCalLED
//   Juan Carlos Cornejo <cornejo@jlab.org> - 2020/03/04
//   A very basic bank decoder to decode HCal LED data
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

#define SSP_DATADEF(b) ((b&0x80000000)>>31)
#define SSP_TAG(b)     ((b&0x78000000)>>27)
#define SSP_SAMPLE(b,c) ((b>>c)&0xFFF)|(((b>>c)&0x1000)?0xFFFFF000:0x0)

#include "SBSHCalLEDModule.h"
#include "THaSlotData.h"
#include <limits>

using namespace std;

namespace Decoder {
	
		Module::TypeIter_t HCalLED::fgThisType =
		DoRegister( ModuleType( "Decoder::HCalLED" , 4450 ));
		
		HCalLED::HCalLED(Int_t crate, Int_t slot) : VmeModule(crate, slot) {
			fDebugFile=0;
				Init();
		}
	
		HCalLED::~HCalLED() {
			
		}
	
		void HCalLED::Init() { 
			Module::Init();
				//    Config(0,25,6,16,128); // should be called by the user
				fDebugFile=0;
				Clear("");
				//    fName = "MPD Module (INFN MPD for GEM and more), use Config to dynamic config";
				fName = "MPD Module";
		}
	
		
		Int_t HCalLED::LoadSlot(THaSlotData *sldat, const UInt_t* evbuffer, Int_t pos, Int_t ) {
			const UInt_t *p = &evbuffer[pos];
			fWordsSeen=0;
			UInt_t j = 0;
			UInt_t thesewords = 0;
			// First word is index, ignore
			j++;
			fWordsSeen++;
			// Second word is LED bit
			fWordsSeen++;
			thesewords=p[j++];
                        Int_t ledbit = thesewords;
			sldat->loadData("adc",1,ledbit,ledbit);
			//fprintf(stderr,"Got LED bit: %d\n",ledbit);
			// Third word is count
			Int_t ledcount=p[j++];
			sldat->loadData("adc",2,ledcount,ledcount);
			fWordsSeen++;

return fWordsSeen;
}

Int_t HCalLED::Decode(const UInt_t* ) {
		return 0;
}


}

ClassImp(Decoder::HCalLED)
