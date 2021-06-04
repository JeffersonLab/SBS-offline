#include "SBSData.h"
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
  TDC::TDC(Float_t offset, Float_t cal) : fHasData(false)
  {
    fEdgeIdx[0] = fEdgeIdx[1];
    SetOffset(offset);
    SetCal(cal);
  }

  void TDC::Process(Float_t val, Int_t edge)
  {
    if(edge < 0 || edge>1) {
      std::cerr << "Edge specified is not valid!" << std::endl;
      edge = 0;
    }
    size_t idx = fEdgeIdx[edge]++;
    if(idx >= fTDC.hits.size()) {
      // Must grow the hits array to accomodate the new hit
      fTDC.hits.push_back(TDCHit());
    }
    TDCHit *hit = &fTDC.hits[idx];
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
  Waveform::Waveform(Float_t ped, Float_t gain, Float_t tcal, Float_t acal) :
    fHasData(false)
  {
    SetPed(ped);
    SetGain(gain);
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
    // TODO: Implement the FADC250 method of doing the pulse integrals
    // for now, I'll just make the integral be the whole window sum,
    // the peak is the sample with the highest number, and time is the
    // sample number of the peak * 4ns, and we'll use the user
    // provided pedestal.
    Float_t time = -404;
    Float_t max  = -404;
    Float_t sum = 0;
    Float_t sped = 0;
    
    for(size_t i = 0; i < vals.size(); i++ ) {
      fSamples.samples_raw[i] = vals[i];
      fSamples.samples[i] = (vals[i]-fSamples.ped)*fSamples.cal;
      sped+=fSamples.ped;
      sum+=vals[i];
      if(vals[i] > max) {
        max = vals[i];
        time = i;
      }
      //fSamples.data_raw_sum += fSamples.data_raw[i];
      //fSamples.data_ped_sum += fSamples.data_ped[i];
      //fSamples.data_cal_sum += fSamples.data_cal[i];
    }
    //printf("ouh!\n");
    fSamples.pulse.integral.raw = sum;
    fSamples.pulse.integral.val = (sum-sped)*fSamples.cal;
    fSamples.pulse.time.raw = time;
    fSamples.pulse.time.val = (time)*fSamples.tcal;
    fSamples.pulse.amplitude.raw = max;
    fSamples.pulse.amplitude.val = (max-fSamples.ped)*fSamples.cal;
/*
    if(vals.size()>0) {
      fSamples.data_raw_sum /= Float_t(vals.size());
      fSamples.data_ped_sum /= Float_t(vals.size());
      fSamples.data_cal_sum /= Float_t(vals.size());
    }*/
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
