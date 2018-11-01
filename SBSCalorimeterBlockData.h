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
    Float_t data_raw_sum;
    Float_t data_ped_sum;
    Float_t data_cal_sum;
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
      Float_t GetDataRaw()  const { return fADC.data_raw; }
      Float_t GetDataPed()  const { return fADC.data_ped; }
      Float_t GetDataCal()  const { return fADC.data_cal; }

      // Setters
      void SetPed(Float_t var)  { fADC.ped = var; }
      void SetGain(Float_t var) { fADC.cal = var; }

      // Process data sets raw value, ped-subtracted and calibrated data
      virtual void Process(Float_t var);

      // Do we have ADC data for this event?
      Bool_t HasData() { return fHasData; }

      // Clear event
      virtual void Clear();

    protected:
      SingleData fADC; ///< ADC single-value data
      Bool_t fHasData;
  };

  ///////////////////////////////////////////////////////////////////////////////
  // TDC single valued
  class TDC {
    public:
      TDC(Float_t offset = 0.0, Float_t cal = 1.0);
      virtual ~TDC() {};

      // Getters
      Float_t GetOffset()   const { return fTDC.ped;      }
      Float_t GetCal()      const { return fTDC.cal;      }
      Float_t GetDataRaw()  const { return fTDC.data_raw; }
      Float_t GetDataPed()  const { return fTDC.data_ped; }
      Float_t GetDataCal()  const { return fTDC.data_cal; }

      // Setters
      void SetOffset(Float_t var)  { fTDC.ped = var; }
      void SetCal(Float_t var) { fTDC.cal = var; }

      // Process data sets raw value, ped-subtracted and calibrated data
      virtual void Process(Float_t var);

      // Do we have TDC data for this event?
      Bool_t HasData() { return fHasData; }

      // Clear event
      virtual void Clear();

    protected:
      SingleData fTDC; ///< TDC single-value data
      Bool_t fHasData;
  };

  ///////////////////////////////////////////////////////////////////////////////
  // Samples (i.e. multi-value ADC data
  class Samples {
    public:
      Samples(Float_t ped = 0.0, Float_t gain = 1.0);
      virtual ~Samples() {};

      // Getters
      Float_t GetPed()  const { return fSamples.ped; }
      Float_t GetGain() const { return fSamples.cal; }
      std::vector<Float_t>& GetDataRaw() { return fSamples.data_raw; }
      std::vector<Float_t>& GetDataPed() { return fSamples.data_ped; }
      std::vector<Float_t>& GetDataCal() { return fSamples.data_cal; }
      Float_t GetDataSumRaw() { return fSamples.data_raw_sum; }
      Float_t GetDataSumPed() { return fSamples.data_ped_sum; }
      Float_t GetDataSumCal() { return fSamples.data_cal_sum; }

      // Setters
      void SetPed(Float_t var)  { fSamples.ped = var; }
      void SetGain(Float_t var) { fSamples.cal = var; }

      // Process data sets raw value, ped-subtracted and calibrated data
      virtual void Process(std::vector<Float_t> &var);

      // Do we have samples data for this event?
      Bool_t HasData() { return fHasData; }

      // Clear event
      virtual void Clear();
    protected:
      MultiData fSamples; ///< Samples single-value data
      Bool_t fHasData;
  };


} // end SBSCalorimeterBlockData namespace

#endif // SBSCALORIMETERBLOCKDATA_H
