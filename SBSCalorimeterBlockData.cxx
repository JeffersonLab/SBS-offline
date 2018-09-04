#include "SBSCalorimeterBlockData.h"

namespace SBSCalorimeterBlockData {

  /////////////////////////////////////////////////////////////////////////////
  // ADC data functions
  ADC::ADC(Float_t ped, Float_t gain) : fHasADCData(false)
  {
    SetPed(ped);
    SetGain(gain);
  }

  void ADC::ProcessADC(Float_t val)
  {
    fADC.data_raw = val;
    fADC.data_ped = val-fADC.ped;
    fADC.data_cal = fADC.data_ped*fADC.cal;
    fHasADCData = true;
  }

  void ADC::ClearADC()
  {
    fADC.data_raw = fADC.data_ped = fADC.data_cal = 0;
    fHasADCData = false;
  }

  /////////////////////////////////////////////////////////////////////////////
  // TDC data functions
  TDC::TDC(Float_t offset, Float_t cal) : fHasTDCData(false)
  {
    SetTDCOffset(offset);
    SetTDCCal(cal);
  }

  void TDC::ProcessTDC(Float_t val)
  {
    fTDC.data_raw = val;
    fTDC.data_ped = val-fTDC.ped;
    fTDC.data_cal = fTDC.data_ped*fTDC.cal;
    fHasTDCData = true;
  }

  void TDC::ClearTDC()
  {
    fTDC.data_raw = fTDC.data_ped = fTDC.data_cal = 0;
    fHasTDCData = false;
  }

  /////////////////////////////////////////////////////////////////////////////
  // Samples data functions
  Samples::Samples(Float_t ped, Float_t gain) : fHasSamplesData(false)
  {
    SetSamplesPed(ped);
    SetSamplesGain(gain);
  }

  void Samples::ProcessADCSamples(std::vector<Float_t> &vals)
  {
    if( vals.size() != fSamples.data_raw.size()) {
      // Resize our data vector
      fSamples.data_raw.resize(vals.size());
      fSamples.data_ped.resize(vals.size());
      fSamples.data_cal.resize(vals.size());
      ClearSamples();
    }
    for(size_t i = 0; i < vals.size(); i++ ) {
      fSamples.data_raw[i] = vals[i];
      fSamples.data_ped[i] = vals[i]-fSamples.ped;
      fSamples.data_cal[i] = fSamples.data_ped[i]*fSamples.cal;
    }
    if(vals.size() > 0)
      fHasSamplesData = true;
  }

  void Samples::ClearSamples()
  {
    for(size_t i = 0; i < fSamples.data_raw.size(); i++) {
      fSamples.data_raw[i] = fSamples.data_ped[i] = fSamples.data_cal[i] = 0;
    }
    fHasSamplesData = false;
  }

}; // end SBSCalorimeterBlockData
