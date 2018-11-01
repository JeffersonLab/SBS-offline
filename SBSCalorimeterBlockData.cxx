#include "SBSCalorimeterBlockData.h"

namespace SBSCalorimeterBlockData {

  /////////////////////////////////////////////////////////////////////////////
  // ADC data functions
  ADC::ADC(Float_t ped, Float_t gain) : fHasData(false)
  {
    SetPed(ped);
    SetGain(gain);
  }

  void ADC::Process(Float_t val)
  {
    fADC.data_raw = val;
    fADC.data_ped = val-fADC.ped;
    fADC.data_cal = fADC.data_ped*fADC.cal;
    fHasData = true;
  }

  void ADC::Clear()
  {
    fADC.data_raw = fADC.data_ped = fADC.data_cal = 0;
    fHasData = false;
  }

  /////////////////////////////////////////////////////////////////////////////
  // TDC data functions
  TDC::TDC(Float_t offset, Float_t cal) : fHasData(false)
  {
    SetOffset(offset);
    SetCal(cal);
  }

  void TDC::Process(Float_t val)
  {
    fTDC.data_raw = val;
    fTDC.data_ped = val-fTDC.ped;
    fTDC.data_cal = fTDC.data_ped*fTDC.cal;
    fHasData = true;
  }

  void TDC::Clear()
  {
    fTDC.data_raw = fTDC.data_ped = fTDC.data_cal = 0;
    fHasData = false;
  }

  /////////////////////////////////////////////////////////////////////////////
  // Samples data functions
  Samples::Samples(Float_t ped, Float_t gain) : fHasData(false)
  {
    SetPed(ped);
    SetGain(gain);
  }

  void Samples::Process(std::vector<Float_t> &vals)
  {
    if( vals.size() != fSamples.data_raw.size()) {
      // Resize our data vector
      fSamples.data_raw.resize(vals.size());
      fSamples.data_ped.resize(vals.size());
      fSamples.data_cal.resize(vals.size());
      Clear();
    }
    fSamples.data_raw_sum = fSamples.data_ped_sum = fSamples.data_cal_sum;
    for(size_t i = 0; i < vals.size(); i++ ) {
      fSamples.data_raw[i] = vals[i];
      fSamples.data_ped[i] = vals[i]-fSamples.ped;
      fSamples.data_cal[i] = fSamples.data_ped[i]*fSamples.cal;
      fSamples.data_raw_sum += fSamples.data_raw[i];
      fSamples.data_ped_sum += fSamples.data_ped[i];
      fSamples.data_cal_sum += fSamples.data_cal[i];
    }
    if(vals.size()>0) {
      fSamples.data_raw_sum /= Float_t(vals.size());
      fSamples.data_ped_sum /= Float_t(vals.size());
      fSamples.data_cal_sum /= Float_t(vals.size());
    }
    fHasData = (vals.size() > 0);
  }

  void Samples::Clear()
  {
    for(size_t i = 0; i < fSamples.data_raw.size(); i++) {
      fSamples.data_raw[i] = fSamples.data_ped[i] = fSamples.data_cal[i] = 0;
    }
    fSamples.data_raw_sum = fSamples.data_ped_sum = fSamples.data_cal_sum = 0.0;
    fHasData = false;
  }

}; // end SBSCalorimeterBlockData
