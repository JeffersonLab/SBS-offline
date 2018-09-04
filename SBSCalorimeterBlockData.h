#ifndef SBSCALORIMETERBLOCKDATA_H
#define SBSCALORIMETERBLOCKDATA_H

#include <vector>
#include <Rtypes.h> // Include standard ROOT types

namespace SBSCalorimeterBlockData {
  ///////////////////////////////////////////////////////////////////////////////
  // Single valued data structure
  struct SingleData {
    // all public variables
    Float_t ped;
    Float_t cal;
    Float_t data_raw;
    Float_t data_ped;
    Float_t data_cal;
  };

  ///////////////////////////////////////////////////////////////////////////////
  // multi-valued data structure
  struct MultiData {
    // all public variables
    Float_t ped;
    Float_t cal;
    std::vector<Float_t> data_raw;
    std::vector<Float_t> data_ped;
    std::vector<Float_t> data_cal;
  };


  ///////////////////////////////////////////////////////////////////////////////
  // ADC single valued
  class ADC {
    public:
      ADC(Float_t ped = 0.0, Float_t gain = 1.0);
      virtual ~ADC() {};

      // Getters
      Float_t GetPed()  const { return fADC.ped;     }
      Float_t GetGain() const { return fADC.cal;    }
      Float_t GetADCDataRaw()  const { return fADC.data_raw; }
      Float_t GetADCDataPed()  const { return fADC.data_ped; }
      Float_t GetADCDataCal()  const { return fADC.data_cal; }

      // Setters
      void SetPed(Float_t var)  { fADC.ped = var; }
      void SetGain(Float_t var) { fADC.cal = var; }

      // Process data sets raw value, ped-subtracted and calibrated data
      virtual void ProcessADC(Float_t var);

      // Do we have ADC data for this event?
      Bool_t HasADCData() { return fHasADCData; }

      // Clear event
      virtual void ClearADC();

    protected:
      SingleData fADC; ///< ADC single-value data
      Bool_t fHasADCData;
  };

  ///////////////////////////////////////////////////////////////////////////////
  // TDC single valued
  class TDC {
    public:
      TDC(Float_t offset = 0.0, Float_t cal = 1.0);
      virtual ~TDC() {};

      // Getters
      Float_t GetTDCOffset()  const { return fTDC.ped;     }
      Float_t GetTDCCal() const { return fTDC.cal;    }
      Float_t GetTDCDataRaw()  const { return fTDC.data_raw; }
      Float_t GetTDCDataPed()  const { return fTDC.data_ped; }
      Float_t GetTDCDataCal()  const { return fTDC.data_cal; }

      // Setters
      void SetTDCOffset(Float_t var)  { fTDC.ped = var; }
      void SetTDCCal(Float_t var) { fTDC.cal = var; }

      // Process data sets raw value, ped-subtracted and calibrated data
      virtual void ProcessTDC(Float_t var);

      // Do we have TDC data for this event?
      Bool_t HasTDCData() { return fHasTDCData; }

      // Clear event
      virtual void ClearTDC();

    protected:
      SingleData fTDC; ///< TDC single-value data
      Bool_t fHasTDCData;
  };

  ///////////////////////////////////////////////////////////////////////////////
  // Samples (i.e. multi-value ADC data
  class Samples {
    public:
      Samples(Float_t ped = 0.0, Float_t gain = 1.0);
      virtual ~Samples() {};

      // Getters
      Float_t GetSamplesPed()  const { return fSamples.ped; }
      Float_t GetSamplesGain() const { return fSamples.cal; }
      std::vector<Float_t>& GetSamplesDataRaw() { return fSamples.data_raw; }
      std::vector<Float_t>& GetSamplesDataPed() { return fSamples.data_ped; }
      std::vector<Float_t>& GetSamplesDataCal() { return fSamples.data_cal; }

      // Setters
      void SetSamplesPed(Float_t var)  { fSamples.ped = var; }
      void SetSamplesGain(Float_t var) { fSamples.cal = var; }

      // Process data sets raw value, ped-subtracted and calibrated data
      virtual void ProcessADCSamples(std::vector<Float_t> &var);

      // Do we have samples data for this event?
      Bool_t HasSamplesData() { return fHasSamplesData; }

      // Clear event
      virtual void ClearSamples();
    protected:
      MultiData fSamples; ///< Samples single-value data
      Bool_t fHasSamplesData;
  };


} // end SBSCalorimeterBlockData namespace

#endif // SBSCALORIMETERBLOCKDATA_H
