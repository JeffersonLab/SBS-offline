#include "SBSData.h"
#include "TMath.h"
#include <iostream>
#define SBS_ADC_MODE_SINGLE 0 //< Simple ADC with only integral
#define SBS_ADC_MODE_MULTI  1 //< FADC 250 mode 7

namespace SBSData {

  /////////////////////////////////////////////////////////////////////////////
  // ADC data functions
  ADC::ADC(Float_t ped, Float_t gain, Float_t tcal, Float_t acal ) :
      fHasData(false), fMode(SBS_ADC_MODE_SINGLE)
  {
    SetPed(ped);
    SetGain(gain);
    SetTimeCal(tcal);
    SetAmpCal(acal);
  }

  void ADC::Process(Float_t val)
  {
    SingleData zero = { 0.0, 0.0 };
    SingleData integral = { val, (val-fADC.ped)*fADC.cal };
    fADC.hits.push_back({integral,zero,zero});
    fHasData = true;
    fMode = SBS_ADC_MODE_SINGLE; //< Mode gets set to simple if this function is called
  }

  void ADC::Process(Float_t integral, Float_t time, Float_t amp, Float_t ped) {
    //fADC.push_back({ped,fGlobalCal,val,val-ped,(val-ped)*fGlobalCal});
    SingleData t_integral = { integral, integral   };
    SingleData t_time     = { time, time*fADC.tcal };
    SingleData t_amp     = { amp, amp*fADC.acal };
    fADC.hits.push_back({t_integral,t_time,t_amp } );
    SetPed(ped);
    fHasData = true;
    fMode = SBS_ADC_MODE_MULTI; //< Mode gets set to multi if this function is called
  }

  void ADC::Clear()
  {
    fADC.good_hit = 0;
    fADC.hits.clear();
    fHasData = false;
  }

  /////////////////////////////////////////////////////////////////////////////
  // TDC data functions
  TDC::TDC(Float_t offset, Float_t cal, Float_t GoodTimeCut) : fHasData(false)
  {
    fEdgeIdx[0] = fEdgeIdx[1]=0;
    SetOffset(offset);
    SetCal(cal);
    SetGoodTimeCut(GoodTimeCut);
  }

  void TDC::Process(Int_t elemID, Float_t val, Float_t fedge)
  {
    Int_t edge = int(fedge);
    // std::cout << " tdc process " << val << " " << edge  << " ftdc hits size = " <<fTDC.hits.size() << " hits in edge "  << fEdgeIdx[edge]<< std::endl;
    if(edge < 0 || edge>1) {
      std::cerr << "Edge specified is not valid!" << std::endl;
      edge = 0;
    }
    size_t idx = fEdgeIdx[edge]++;
    if(idx >= fTDC.hits.size()) {
      // Must grow the hits array to accomodate the new hit
      // if ( edge ==1)   std::cout << " First edge is TE , this is not right: " << "  idx = " << idx  << " ftdc hits size = " <<fTDC.hits.size() << " hits with LE ="  << fEdgeIdx[0] << " hits with TE ="  << fEdgeIdx[1] << std::endl;
      fTDC.hits.push_back(TDCHit());
    }
    TDCHit *hit = &fTDC.hits[idx];
    hit->elemID = elemID;
    if( edge == 0 ) { // Leading edge
      hit->le.raw = val;
      hit->le.val = (val-fTDC.offset)*fTDC.cal;
    } else {
      hit->te.raw = val;
      hit->te.val = (val-fTDC.offset)*fTDC.cal;
    }
    if(fEdgeIdx[0] == fEdgeIdx[1]) { // Both leading and trailing edges now found
      hit->ToT.raw = hit->te.raw - hit->le.raw;
      hit->ToT.val = hit->te.val - hit->le.val;
   }
    if(fEdgeIdx[1] > fEdgeIdx[0]) fEdgeIdx[0] = fEdgeIdx[1]; // if TE found first force LE count to increase
    fHasData = true;
  }

  void TDC::Clear()
  {
    fEdgeIdx[0] = fEdgeIdx[1] = 0;
    fTDC.hits.clear();
    fHasData = false;
    fTDC.good_hit = 0;
  }

  /////////////////////////////////////////////////////////////////////////////
  // Waveform data functions
  Waveform::Waveform(Float_t ped, Float_t gain, Float_t ChanTomV, Float_t tcal, Float_t acal) :
    fHasData(false)
  {
    SetPed(ped);
    SetGain(gain);
    SetChanTomV(ChanTomV);
    SetTimeCal(tcal);
    SetAmpCal(acal);
  }

  void Waveform::Process(std::vector<Float_t> &vals)
  {
    //printf("vals size %d, samples raw size %d\n", vals.size(), fSamples.samples_raw.size());
    if( vals.size() != fSamples.samples_raw.size()) {
      // Resize our data vector
      fSamples.samples_raw.resize(vals.size());
      fSamples.samples.resize(vals.size());
      Clear();
    }
    //
    for(size_t i = 0; i < vals.size(); i++ ) {
      fSamples.samples_raw[i] = vals[i]*fSamples.ChanTomV;
    }
    // First determine pedestal in first four samples
    Double_t pedsum=0;
    Int_t NPedsum=GetNPedBin();
    NPedsum= TMath::Min(NPedsum,int(vals.size()));
    //       std::cout << " Npedsum = " << NPedsum << " " << GetNPedBin() << std::endl ;
    for(Int_t i = 0; i <NPedsum ; i++ ) {
      pedsum+=fSamples.samples_raw[i];
    }
    SetPed(pedsum/NPedsum);
    for(size_t i = 0; i < vals.size(); i++ ) {
      fSamples.samples[i] = (fSamples.samples_raw[i]-fSamples.ped)*fSamples.cal;
    }
    // Try and fixd sample threshold crossing above the pedestal
    Float_t ThresVal=GetThres(); // mV
    UInt_t ThresCrossBin=TMath::Max(NPedsum-1,0);
    //    std::cout << " ped = " << fSamples.ped << " thres = " << ThresVal << std::endl ;
    while ( fSamples.samples_raw[ThresCrossBin] < fSamples.ped+ThresVal && ThresCrossBin < vals.size() ) {
        ThresCrossBin++;
    }
     //
    // if threshold crossing found
    UInt_t NSB = GetNSB();
    UInt_t NSA = GetNSA();
    UInt_t FixedThresCrossBin=GetFixThresBin();
    Float_t FineTime = 0;
    Float_t max  = 0;
    Float_t sum = 0;
    Float_t sped = 0;
    
    Bool_t PeakFound= kFALSE;
      UInt_t PeakBin= 0;
      UInt_t IntMinBin= 0;
      UInt_t IntMaxBin= vals.size();
    if (ThresCrossBin < vals.size()) {
      IntMinBin= TMath::Max(ThresCrossBin-NSB,IntMinBin);
      IntMaxBin= TMath::Min(ThresCrossBin+NSA-1,IntMaxBin);
    } else {
      IntMinBin = TMath::Max(FixedThresCrossBin-NSB,IntMinBin);
      IntMaxBin = TMath::Min(FixedThresCrossBin+NSA-1,IntMaxBin);
    }
    // convert to pC, assume tcal is in ns, and 50ohm resistance
    Float_t pC_Conv = fSamples.tcal/50.;
    
    for(size_t i =IntMinBin ; i <IntMaxBin ; i++ ) {
         sped+=fSamples.ped*pC_Conv;
         sum+=fSamples.samples_raw[i]*pC_Conv;
         if ( i > ThresCrossBin && !PeakFound) {
	   if (fSamples.samples_raw[i] > max) {
	     max = fSamples.samples_raw[i];
	   } else {
             PeakFound= kTRUE;
	     PeakBin = i-1;
	   }
	 }
    }
    //
    //    std::cout << " Int = " << IntMinBin << " " << IntMaxBin<< " ThresCrossBin =   " << ThresCrossBin << " peak-found " << PeakFound << std::endl ;
    //
    Float_t VMid = (max+fSamples.ped)/2.;
    if (PeakFound) {
      for(size_t i =IntMinBin ; i <PeakBin+1 ; i++ ) {
	if (VMid >= fSamples.samples_raw[i]  && VMid < fSamples.samples_raw[i+1]) {
	FineTime = i+(VMid-fSamples.samples_raw[i])/(fSamples.samples_raw[i+1]-fSamples.samples_raw[i]);
	}
      }
    }
    
    /*
    if (ThresCrossBin==vals.size()) {
      std::cout << " Threshold = " << fThresVal << " ped = " << fSamples.ped << " element = " << ThresCrossBin << " " << " Vmid = " << VMid << "  sum = " << sum << " max = " << max << " time = " << FineTime << " tcal = " << fSamples.tcal <<  std::endl;
    }
    */
    fSamples.pulse.integral.raw = sum;
    fSamples.pulse.integral.val = (sum-sped)*fSamples.cal;
    fSamples.pulse.time.raw = FineTime;
    fSamples.pulse.time.val = (FineTime)*fSamples.tcal;
    fSamples.pulse.amplitude.raw = max;
    fSamples.pulse.amplitude.val = (max-fSamples.ped)*fSamples.cal;
    if (max==0) fSamples.pulse.amplitude.val=max;
    //
    fHasData = (vals.size() > 0);
  }

  void Waveform::Clear()
  {
    for(size_t i = 0; i < fSamples.samples.size(); i++) {
      fSamples.samples_raw[i] = fSamples.samples[i] = 0;
    }
    fHasData = false;
  }

}; // end SBSData
