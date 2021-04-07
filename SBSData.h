#ifndef SBSDATA_H
#define SBSDATA_H

#include <vector>
#include <Rtypes.h> // Include standard ROOT types

namespace SBSData {
  ///////////////////////////////////////////////////////////////////////////////
  // Single valued data structure
  struct SingleData {
    Float_t raw; //< Raw data
    Float_t val; //< Calibrated data
  };

  ///////////////////////////////////////////////////////////////////////////////
  // Single valued Pulse ADC data
  struct PulseADCData {
    SingleData integral;
    SingleData time;
    SingleData amplitude;
  };

  ///////////////////////////////////////////////////////////////////////////////
  // Multi-function ADC data structure
  struct ADCData {
    // all public variables
    Float_t ped; //< Pedestal
    Float_t cal; //< Calibration constant for ADC integral
    Float_t tcal; //< Calibration constant for TDC value
    Float_t acal; //< Calibration constant for ADC amplitude (peak)
    std::vector<PulseADCData> hits; //< Pulse data
    UInt_t good_hit; //< Index of good hit
  };

  ///////////////////////////////////////////////////////////////////////////////
  // Waveform ADC Data
  struct WaveformData {
    // all public variables
    Float_t cal; //< Single calibration for all samples
    Float_t ped; //< Single pedestal for all samples
    Float_t tcal; //< Calibration constant for TDC value
    Float_t acal; //< Calibration constant for ADC amplitude (peak)
    std::vector<Float_t> samples_raw; //< Raw samples
    std::vector<Float_t> samples;     //< Calibrated samples
    PulseADCData         pulse;       //< Pulse information
  };

  ///////////////////////////////////////////////////////////////////////////////
  // TDC data with both leading and trailing edge info
  struct TDCHit {
    SingleData le;  //< Leading edge time
    SingleData te; //< Trailing edge time
    SingleData ToT;   //< Time-over-Threshold (if trailing edge provided)
  };

  ///////////////////////////////////////////////////////////////////////////////
  // TDC data structure
  struct TDCData {
    // all public variables
    Float_t offset; //< Time offset
    Float_t cal;    //< Conversion factor
    std::vector<TDCHit> hits;
    UInt_t good_hit; //< Index of good hit
  };

  ///////////////////////////////////////////////////////////////////////////////
  // ADC single valued
  class ADC {
    public:
      ADC(Float_t ped = 0.0, Float_t gain = 1.0, Float_t tcal = 4.0,
          Float_t acal = 1.0);
      virtual ~ADC() {};

      // Getters
      Float_t GetPed()                  const { return fADC.ped;            }
      Float_t GetGain()                 const { return fADC.cal;            }
      PulseADCData GetHit(UInt_t i)     const { return fADC.hits[i];        }
      PulseADCData GetGoodHit()         const { return fADC.hits[fADC.good_hit]; }
      SingleData GetIntegral(UInt_t i)  const { return GetHit(i).integral;  }
      SingleData GetTime(UInt_t i)      const { return GetHit(i).time;      }
      SingleData GetAmplitude(UInt_t i) const { return GetHit(i).amplitude; }
      ADCData GetADC()                  const { return fADC;                }
      Int_t GetNHits()                        { return fADC.hits.size();    }

      // Some additional helper functions for easy access to the ADC integral
      Float_t GetDataRaw(UInt_t i)      const { return GetHit(i).integral.raw; }
      Float_t GetData(UInt_t i)         const { return GetHit(i).integral.val; }
      std::vector<PulseADCData> GetAllHits()  const { return  fADC.hits; }

      // Setters
      void SetPed(Float_t var)  { fADC.ped = var; }
      void SetGain(Float_t var) { fADC.cal = var; }
      void SetTimeCal(Float_t var) { fADC.tcal = var; }
      void SetAmpCal(Float_t var) { fADC.acal = var; }
      void SetGoodHit(UInt_t i) { fADC.good_hit = i; }

      // Process data sets raw value, ped-subtracted and calibrated data
      virtual void Process(Float_t var);
      virtual void Process(Float_t integral, Float_t time,
          Float_t amp, Float_t pedestal);

      // Do we have ADC data for this event?
      Bool_t HasData() { return fHasData; }

      // Clear event
      virtual void Clear();

    protected:
      ADCData fADC; ///< ADC single-value data
      Bool_t fHasData;
      Int_t  fMode; //< ADC mode where 0 == simple, 1 == multi function
  };

  ///////////////////////////////////////////////////////////////////////////////
  // TDC single valued
  class TDC {
    public:
      TDC(Float_t offset = 0.0, Float_t cal = 1.0);
      virtual ~TDC() {};

      // Getters
      Float_t GetOffset()           const { return fTDC.offset;     }
      Float_t GetCal()              const { return fTDC.cal;        }
      TDCHit GetHit(UInt_t i)       const { return fTDC.hits[i];    }
      SingleData GetLead(UInt_t i)  const { return GetHit(i).le;    }
      SingleData GetTrail(UInt_t i) const { return GetHit(i).te;    }
      TDCHit GetGoodHit()           const { return fTDC.hits[fTDC.good_hit]; }
      Int_t GetNHits()                    { return fTDC.hits.size();    }
      std::vector<TDCHit> GetAllHits() const { return  fTDC.hits; }

      // Helper functions to get leading edge info
      Float_t GetData(UInt_t i)     const { return GetLead(i).val;  }
      Float_t GetDataRaw(UInt_t i)  const { return GetLead(i).raw;  }
      Float_t GetToT(UInt_t i)      const { return GetHit(i).ToT.val; }

      // Setters
      void SetOffset(Float_t var)  { fTDC.offset = var; }
      void SetCal(Float_t var) { fTDC.cal = var; }
      void SetGoodHit(UInt_t i) { fTDC.good_hit = i; }

      // Process data sets raw value, ped-subtracted and calibrated data
      virtual void Process(Float_t var, Int_t edge = 0);

      // Do we have TDC data for this event?
      Bool_t HasData() { return fHasData; }

      // Clear event
      virtual void Clear();

    protected:
      TDCData fTDC; ///< TDC hit container
      Bool_t fHasData;
      size_t fEdgeIdx[2]; //< Current index of the next hit data
  };

  ///////////////////////////////////////////////////////////////////////////////
  // Samples (e.g. ADC Waveform data)
  class Waveform {
    public:
      Waveform(Float_t ped = 0.0, Float_t gain = 1.0, Float_t tcal = 4.0,
          Float_t acal = 1.0);
      virtual ~Waveform() {};

      // Getters
      Float_t GetPed()  const { return fSamples.ped; }
      Float_t GetGain() const { return fSamples.cal; }
      std::vector<Float_t>& GetDataRaw() { return fSamples.samples_raw; }
      std::vector<Float_t>& GetData() { return fSamples.samples; }
      PulseADCData GetPulse() { return fSamples.pulse; }
      SingleData GetIntegral()   { return fSamples.pulse.integral; }
      SingleData GetTime()   { return fSamples.pulse.time; }
      SingleData GetAmplitude()   { return fSamples.pulse.amplitude; }

      // Setters
      void SetPed(Float_t var)  { fSamples.ped = var; }
      void SetGain(Float_t var) { fSamples.cal = var; }
      void SetTimeCal(Float_t var) { fSamples.tcal = var; }
      void SetAmpCal(Float_t var) { fSamples.acal = var; }

      // Process data sets raw value, ped-subtracted and calibrated data
      virtual void Process(std::vector<Float_t> &var);

      // Do we have samples data for this event?
      Bool_t HasData() { return fHasData; }

      // Clear event
      virtual void Clear();
    protected:
      WaveformData fSamples; ///< Samples single-value data
      Bool_t fHasData;
  };


} // end SBSData namespace

#endif // SBSDATA_H
