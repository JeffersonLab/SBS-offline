/////////////////////////////////////////////////////////////////////
//
//   SBSDecodeF1TDCModule
//   JLab F1 TDC Module
//
/////////////////////////////////////////////////////////////////////
//
// - Version 1 - 2017/12/13 (Marco Carmignotto)
//    * Tested for High Resolution only
//    * Work with multiple modules in the same bank
//
//  - Version 0 - 2017/11/30 (Marco Carmignotto)
//    * Multihits
//
//  - Version -1 (Robert Michaels)
//
/////////////////////////////////////////////////////////////////////

#include "SBSDecodeF1TDCModule.h"
#include "VmeModule.h"
#include "THaEvData.h"
#include "TMath.h"
#include <iostream>
#include <string>
#include <sstream>

using namespace std;

namespace Decoder {

//  const Int_t NTDCCHAN = 32;
  const Int_t MAXHIT   = 20;

  Module::TypeIter_t SBSDecodeF1TDCModule::fgThisType =
    DoRegister( ModuleType( "Decoder::SBSDecodeF1TDCModule" , 3204 ));

SBSDecodeF1TDCModule::SBSDecodeF1TDCModule(Int_t crate, Int_t slot) : VmeModule(crate, slot) {
  fDebugFile=0;
  Init();
}

SBSDecodeF1TDCModule::~SBSDecodeF1TDCModule() {
  //if (fTdcData) delete [] fTdcData;
  //if (fNumHits) delete [] fNumHits;
  if (fTdcData.size()>0) fTdcData.clear();
  if (fNumHits.size()>0) fNumHits.clear();
  if(F1slots) delete [] F1slots;
}

void SBSDecodeF1TDCModule::Init() {
  //fTdcData = new Int_t[NTDCCHAN*MAXHIT];
  //fNumHits = new Int_t[NTDCCHAN*MAXHIT];
  fDebugFile=0;
  Clear();
  IsInit = kTRUE;
  fName = "My private F1 TDC 3204";
  nF1=0;
  F1slots = new Int_t[50];
  memset(F1slots, 0, 50*sizeof(Int_t));
  fWdcntMask=0;
  SetResolution(1);
  fNumChan = 32;
  if (fModelNum == 6401) {
    SetResolution(0);
    fNumChan = 64;
  }
}


Bool_t SBSDecodeF1TDCModule::IsSlot(UInt_t rdata)
{
  if (fDebugFile)
    *fDebugFile << "is SBSDecodeF1TDC slot ? "<<hex<<fHeader
		<<"  "<<fHeaderMask<<"  "<<rdata<<dec<<endl;
  return ((rdata != 0xffffffff) & ((rdata & fHeaderMask)==fHeader));
}

Int_t SBSDecodeF1TDCModule::GetData(Int_t chan, Int_t hit) const
{
  Int_t idx = chan*MAXHIT + hit;
  if (idx < 0 || idx > MAXHIT*fNumChan*nF1) return 0;
  return fTdcData[idx];
}

Int_t SBSDecodeF1TDCModule::GetNumHits(Int_t chan) {
  if (chan < 0 || chan > nF1*fNumChan) return -1;
  return fNumHits[chan];
}

void SBSDecodeF1TDCModule::Clear(const Option_t* opt) {
  VmeModule::Clear(opt);
  //memset(fTdcData, 0, fNumChan*MAXHIT*sizeof(Int_t));
  //memset(fNumHits, 0, fNumChan*sizeof(Int_t));
//  std::fill(fTdcData.begin(), fTdcData.end(), 0);
//  std::fill(fNumHits.begin(), fNumHits.end(), 0);
  for(Int_t i=0; i<fTdcData.size(); i++) fTdcData[i]=0;
  for(Int_t i=0; i<fNumHits.size(); i++) fNumHits[i]=0;
}

Int_t SBSDecodeF1TDCModule::LoadSlot(THaSlotData *sldat, const UInt_t *evbuffer, Int_t pos, Int_t len) {
// the 4-arg version of LoadSlot.  Let it call the 3-arg version.
// I'm not sure we need both (historical)

    return LoadSlot(sldat, evbuffer+pos, evbuffer+pos+len);

  }

Int_t SBSDecodeF1TDCModule::LoadSlot(THaSlotData *sldat, const UInt_t *evbuffer, const UInt_t *pstop) {
// This is the 3-arg version of LoadSlot
// Note, this increments evbuffer
  if (fDebugFile) *fDebugFile << "SBSDecodeF1TDCModule:: loadslot "<<endl;
  fWordsSeen = 0;

  // CAUTION: this routine re-numbers the channels
  // compared to the labelled numbering scheme on the board itself.
  // According to the labelling and internal numbering scheme,
  // the F1 module has odd numbered channels on one connector
  // and even numbered channels on the other.
  // However we usually put neighboring blocks/wires into the same
  // cable, connector etc.
  // => hana therefore uses a numbering scheme different from the module
  //
  // In normal resolution mode, the scheme is:
  //    connector 1:  ch 0 - 15
  //    conncetor 2:  ch 16 - 31
  //    connector 33: ch 32 - 47
  //    connector 34: ch 48 - 63
  //
  // In high-resolution mode, only two connectors are used since
  // two adjacent channels are internally combined and read out as the
  // internally-even numbered channels.
  // this is kind of inconvenient for the rest of the software
  // => hana therefore uses a numbering scheme different from the module
  //    connector 1:  unused
  //    connector 2:  ch 0 - 15
  //    connector 33: unused
  //    connector 34: ch 16 - 31
  //
  // In both modes:
  // it is assumed that we only get data from one single trigger
  // if the F1 is run in multiblock mode (buffered mode)
  // this might not be the case anymore - but this will be interesting anyhow
  // triggertime and eventnumber are not yet read out, they will again
  // be useful when multiblock mode (buffered mode) is used

   const UInt_t F1_HIT_OFLW = 1<<24; // bad
   const UInt_t F1_OUT_OFLW = 1<<25; // bad
   const UInt_t F1_RES_LOCK = 1<<26; // good
   const UInt_t DATA_CHK = F1_HIT_OFLW | F1_OUT_OFLW | F1_RES_LOCK;
   const UInt_t DATA_MARKER = 1<<23;
   // look at all the data
   const UInt_t *loc = evbuffer;
   Int_t fDebug=2;
   if(fDebug > 1 && fDebugFile!=0) *fDebugFile<< "Debug of (My Private) SBSDecodeF1TDC data, fResol =  "<<fResol<<"  model num  "<<fModelNum<<endl;
   // To account for multiple slots
   Int_t lastSlot = 0;
   Int_t idxSlot = -1;
   while ( loc <= pstop && IsSlot(*loc) ) {
      if ( !( (*loc) & DATA_MARKER ) ) {
	// header/trailer word, to be ignored
	 if(fDebug > 1 && fDebugFile!=0)
	    *fDebugFile<< "[" << (loc-evbuffer) << "] header/trailer  0x"
		       <<hex<<*loc<<dec<<endl;
	} else {
	    if (fDebug > 1 && fDebugFile!=0)
	       *fDebugFile<< "[" << (loc-evbuffer) << "] data            0x"
			  <<hex<<*loc<<dec<<endl;
	    Int_t chn = ((*loc)>>16) & 0x3f;  // internal channel number

	    Int_t chan;
	    if (IsHiResolution()) {
		// drop last bit for channel renumbering
		 chan=(chn >> 1);
	    } else {
		// do the reordering of the channels, for contiguous groups
		// odd numbered TDC channels from the board -> +16
		chan = (chn & 0x20) + 16*(chn & 0x01) + ((chn & 0x1e)>>1);
	    }
	     Int_t f1slot = ((*loc)&0xf8000000)>>27;
		//FIXME: cross-check slot number here
		// Make sure we have the current slot considered - important for multiple modules
		if(f1slot != lastSlot && f1slot!=30) { // we may have a new slot here
			Bool_t newslot = 1;
			for(Int_t k=0; k<nF1; k++) {
				if(f1slot == F1slots[k]) newslot=0; // we saw this slot before
			}
			if(newslot) { // increase size of vectors
				F1slots[nF1]=f1slot;
				nF1++;
				fNumHits.resize(fNumChan*nF1);
				fTdcData.resize(fNumChan*nF1*MAXHIT);
			}
			lastSlot = f1slot;
			// Updating idxSlot
			for(Int_t k=0; k<nF1; k++) {
				if(f1slot == F1slots[k]) idxSlot = k;
			}
		}
//cout << "f1slot: " << f1slot << " ";//endl;
	     if (f1slot!=30 && ((*loc) & DATA_CHK) != F1_RES_LOCK ) {
	       cout << "\tWarning: F1 TDC " << hex << (*loc) << dec;
	       cout << "\tSlot (Ch) = " << f1slot << "(" << chan << ")";
	       if ( (*loc) & F1_HIT_OFLW ) {
		  cout << "\tHit-FIFO overflow";
	       }
	       if ( (*loc) & F1_OUT_OFLW ) {
		  cout << "\tOutput FIFO overflow";
	       }
	       if ( ! ((*loc) & F1_RES_LOCK ) ) {
		  cout << "\tResolution lock failure!";
	       }
	       cout << endl;
	     }

	      Int_t raw= (*loc) & 0xffff;
	      if(fDebug > 1 && fDebugFile!=0) {
		*fDebugFile<<" int_chn chan data "<<dec<<chn<<"  "<<chan
		    <<"  0x"<<hex<<raw<<dec<<endl;
	      }
	      /*Int_t status = */sldat->loadData("tdc",chan,raw,raw);

		  if((*loc)!=0xf1daffff && nF1>0 && f1slot!=30) { // Make sure memory is allocated to save data
		    Int_t chSlot = idxSlot*fNumChan + chan;
	        Int_t idx = chSlot*MAXHIT + fNumHits[chSlot]; // multiple hits, multiple slots 
	        if (idx >= 0 && idx < MAXHIT*nF1*fNumChan  && fNumHits[chSlot]<MAXHIT) fTdcData[idx] = raw;
		    if(fNumHits[chSlot]>=MAXHIT) cout << "F1TDC warning: more than " << MAXHIT << " hits in channel " << chan << " / slot " << f1slot << ". Not taking last hit(s), although counting them." << endl;
		    (fNumHits[chSlot])++;
//cout << "chSlot=" << chSlot << " - lastSlot=" << lastSlot << " - nF1=" << nF1 << " -- ch=" << chan << " - word=0x" << hex << (*loc) << dec << " - raw=" << raw << endl;
		  }
	      fWordsSeen++;
	  }
       loc++;
   }

  return fWordsSeen;
}

}

ClassImp(Decoder::SBSDecodeF1TDCModule)
