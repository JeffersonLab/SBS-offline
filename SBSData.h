#ifndef SBSDATA_H
#define SBSDATA_H

#include <vector>
#include <Rtypes.h> // Include standard ROOT types

namespace SBSData {
  ///////////////////////////////////////////////////////////////////////////////
  // Single valued data structure
  struct SingleData {
    Double_t raw; //< Raw data
    Double_t val; //< Calibrated data
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
    Double_t ped; //< Pedestal
    Double_t cal; //< Calibration constant for ADC integral
    Double_t tcal; //< Calibration constant for TDC value
    Double_t acal; //< Calibration constant for ADC amplitude (peak)
    Double_t trigcal; //< ratio Trig amp to FADC amp
    Double_t timeoffset; //<  time offset (ns) for ADC time
    UInt_t NSB; //< Programmed # of samples before threshold in FADC
    UInt_t NSA; //< Programmed # of samples after threshold in FADC
    Int_t NPedBin; //< Programmed # of samples used in pedestal average in FADC
    Double_t GoodTimeCut;    //< Time Cut to select good hit in multihit ADC.
    Double_t ChanTomV;    //< milliVolts/channel for FADC
    std::vector<PulseADCData> hits; //< Pulse data
    Int_t good_hit; //< Index of good hit
  };

  ///////////////////////////////////////////////////////////////////////////////
  // Waveform ADC Data
  struct WaveformData {
    // all public variables
    Double_t cal; //< Single calibration for all samples
    Double_t ped; //< Single pedestal for all samples
    Double_t thres; //< Threshold for determining the Threshold Bin.
    Double_t ChanTomV; //< Conversion of ADC channel to milliVolt.
    Double_t tcal; //< Calibration constant for TDC value
    Double_t acal; //< Calibration constant for ADC amplitude (peak)
    Double_t trigcal; //< ratio Trig amp to FADC amp
    Double_t timeoffset; //<  time offset (ns) for ADC time
    UInt_t FixThresBin; //<Fixed Threshold Bin when no peak found
    UInt_t NSB; //<Number of bins before Threshold Bin integrate when Threshold Bin found
    UInt_t NSA; //<Number of bins after Threshold Bin integrate when Threshold Bin found
    Int_t NPedBin; //<Number of bins used at beginning of smaple window
    Double_t GoodTimeCut;    //< Time Cut to select good hit in multihit ADC.
    Int_t good_hit; //< Index of good hit
    std::vector<Double_t> samples_raw; //< Raw samples
    std::vector<Double_t> samples;     //< Calibrated samples
    PulseADCData         pulse;       //< Pulse information
  };

  ///////////////////////////////////////////////////////////////////////////////
  // TDC data with both leading and trailing edge info
  struct TDCHit {
    SingleData le;  //< Leading edge time
    SingleData te; //< Trailing edge time
    SingleData ToT;   //< Time-over-Threshold (if trailing edge provided)
    Int_t elemID; //< Element ID
    UInt_t TrigTime; //< F1 trigger time
 };

  ///////////////////////////////////////////////////////////////////////////////
  // TDC data structure
  struct TDCData {
    // all public variables
    Double_t offset; //< Time offset
    Double_t cal;    //< Conversion factor
    Double_t GoodTimeCut;    //< Time Cut to select good hit in multihit TDC.
    std::vector<TDCHit> hits;
    Int_t good_hit; //< Index of good hit
  };

  ///////////////////////////////////////////////////////////////////////////////
  // ADC single valued
  class ADC {
    public:
    ADC(Double_t ped = 0.0, Double_t gain = 1.0, Double_t tcal = 4.0);
      virtual ~ADC() {};

      // Getters
      Double_t GetPed()                  const { return fADC.ped;            }
      Double_t GetGain()                 const { return fADC.cal;            }
      Double_t GetTimeCal()                 const { return fADC.tcal; }
      Double_t GetGoodTimeCut()              const { return fADC.GoodTimeCut;}
      Double_t GetAmpCal()                  const { return fADC.acal; }
      Double_t GetTrigCal()                  const { return fADC.trigcal; }
      Double_t GetTimeOffset()                  const { return fADC.timeoffset; }
      UInt_t GetNSA()                      const { return fADC.NSA; }
      UInt_t GetNSB()                      const { return fADC.NSB; }
      Int_t GetNPedBin()                      const { return fADC.NPedBin; }
      Double_t GetChanTomV() const { return fADC.ChanTomV; }
      Double_t GetGoodHitIndex()            const { return fADC.good_hit; }
      PulseADCData GetHit(UInt_t i)     const { return fADC.hits[i];        }
      PulseADCData GetGoodHit()         const { return fADC.hits[fADC.good_hit]; }
      SingleData GetIntegral(UInt_t i)  const { return GetHit(i).integral;  }
      SingleData GetTime(UInt_t i)      const { return GetHit(i).time;      }
      SingleData GetAmplitude(UInt_t i) const { return GetHit(i).amplitude; }
      ADCData GetADC()                  const { return fADC;                }
      Int_t GetNHits()                        { return fADC.hits.size();    }

      // Some additional helper functions for easy access to the ADC integral
      Double_t GetDataRaw(UInt_t i)      const { return GetHit(i).integral.raw; }
      Double_t GetTimeData(UInt_t i)         const { return GetHit(i).time.val; }
      Double_t GetData(UInt_t i)         const { return GetHit(i).integral.val; }
     std::vector<PulseADCData> GetAllHits()  const { return  fADC.hits; }

      // Setters
      void SetPed(Double_t var)  { fADC.ped = var; }
      void SetGain(Double_t var) { fADC.cal = var; }
      void SetTimeCal(Double_t var) { fADC.tcal = var; }
      void SetGoodTimeCut(Double_t var) { fADC.GoodTimeCut = var; }
      void SetAmpCal(Double_t var) { fADC.acal = var; }
      void SetTrigCal(Double_t var) { fADC.trigcal = var; }
      void SetTimeOffset(Double_t var) { fADC.timeoffset = var; }
      void SetGoodHit(Int_t i) { fADC.good_hit = i; }
      void SetChanTomV(Double_t var) { fADC.ChanTomV = var; }
      void SetADCParam(Double_t i1,Int_t i2,Int_t i3, Int_t i4,Double_t i5) { fADC.ChanTomV=i1;fADC.NSB=i2;fADC.NSA=i3;fADC.NPedBin=i4;fADC.GoodTimeCut=i5;}

      // Process data sets raw value, ped-subtracted and calibrated data
      virtual void Process(Double_t var);
      virtual void Process(Double_t integral, Double_t time,
          Double_t amp, Double_t pedestal);

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
      TDC(Double_t offset = 0.0, Double_t cal = 1.0, Double_t GoodTimeCut = 1.0);
      virtual ~TDC() {};

      // Getters
      Double_t GetOffset()           const { return fTDC.offset;     }
      Double_t GetCal()              const { return fTDC.cal;        }
      Double_t GetGoodTimeCut()              const { return fTDC.GoodTimeCut;}
      Double_t GetGoodHitIndex()            const { return fTDC.good_hit; }
      TDCHit GetHit(UInt_t i)       const { return fTDC.hits[i];    }
      SingleData GetLead(UInt_t i)  const { return GetHit(i).le;    }
      SingleData GetTrail(UInt_t i) const { return GetHit(i).te;    }
      TDCHit GetGoodHit()           const { return fTDC.hits[fTDC.good_hit]; }
      Int_t GetNHits()                    { return fTDC.hits.size();    }
      UInt_t GetTrigTime(UInt_t i)  const { return GetHit(i).TrigTime;   }
      std::vector<TDCHit> GetAllHits() const { return  fTDC.hits; }

      // Helper functions to get leading edge info
      Double_t GetData(UInt_t i)     const { return GetLead(i).val;  }
      Double_t GetDataRaw(UInt_t i)  const { return GetLead(i).raw;  }
      Double_t GetToT(UInt_t i)      const { return GetHit(i).ToT.val; }

      // Setters
      void SetOffset(Double_t var)  { fTDC.offset = var; }
      void SetCal(Double_t var) { fTDC.cal = var; }
      void SetGoodTimeCut(Double_t var) { fTDC.GoodTimeCut = var; }
      void SetGoodHit(Int_t i) { fTDC.good_hit = i; }
      
      // Process data sets raw value, ped-subtracted and calibrated data
      virtual void Process(Int_t elemID, Double_t var, Double_t edge = 0);
      virtual void ProcessSimple(Int_t elemID, Double_t var, Int_t nhit,UInt_t TrigTime);

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
      explicit Waveform( Double_t ped = 0.0, Double_t gain = 1.0, Double_t ChanTomV = 0.48828,
                       Double_t GoodTimeCut = 1.0, Double_t tcal = 4.0 );
      virtual ~Waveform() = default;;

      // Getters
      Double_t GetPed()         const { return fSamples.ped; }
      Double_t GetGain()        const { return fSamples.cal; }
      Double_t GetThres()       const { return fSamples.thres; }
      Double_t GetChanTomV()    const { return fSamples.ChanTomV; }
      Double_t GetAmpCal()      const { return fSamples.acal; }
      Double_t GetTrigCal()     const { return fSamples.trigcal; }
      Double_t GetTimeOffset()  const { return fSamples.timeoffset; }
      UInt_t GetFixThresBin()   const { return fSamples.FixThresBin; }
      UInt_t GetNSB()           const { return fSamples.NSB; }
      UInt_t GetNSA()           const { return fSamples.NSA; }
      UInt_t GetNPedBin()       const { return fSamples.NPedBin; }
      Double_t GetGoodTimeCut() const { return fSamples.GoodTimeCut;}
      Int_t GetGoodHitIndex()   const { return fSamples.good_hit; }
      std::vector<Double_t>& GetDataRaw() { return fSamples.samples_raw; }
      std::vector<Double_t>& GetData() { return fSamples.samples; }
      PulseADCData GetPulse()   const { return fSamples.pulse; }
      SingleData GetIntegral()  const { return fSamples.pulse.integral; }
      SingleData GetTime()      const { return fSamples.pulse.time; }
      SingleData GetAmplitude() const { return fSamples.pulse.amplitude; }
      Double_t GetTimeData()    const { return fSamples.pulse.time.val; }

      // Setters
      void SetValTime(Double_t var)  { fSamples.pulse.time.val = var; }
      void SetPed(Double_t var)  { fSamples.ped = var; }
      void SetGain(Double_t var) { fSamples.cal = var; }
      void SetChanTomV(Double_t var) { fSamples.ChanTomV = var; }
      void SetTimeCal(Double_t var) { fSamples.tcal = var; }
      void SetGoodTimeCut(Double_t var) { fSamples.GoodTimeCut = var; }
      void SetAmpCal(Double_t var) { fSamples.acal = var; }
      void SetTrigCal(Double_t var) { fSamples.trigcal = var; }
      void SetTimeOffset(Double_t var) { fSamples.timeoffset = var; }
      void SetGoodHit(Int_t i) { fSamples.good_hit = i; }
      void SetWaveformParam(Double_t var,Int_t i1,Int_t i2,Int_t i3, Int_t i4) { fSamples.thres = var;fSamples.FixThresBin=i1;fSamples.NSB=i2;fSamples.NSA=i3;fSamples.NPedBin=i4;}
      // Process data sets raw value, ped-subtracted and calibrated data
      virtual void Process(std::vector<Double_t> &var);

      // Do we have samples data for this event?
      Bool_t HasData() const { return fHasData; }

      // Clear event
      virtual void Clear();
    protected:
      WaveformData fSamples; ///< Samples single-value data
      Bool_t fHasData;
  };


} // end SBSData namespace

#endif // SBSDATA_H
