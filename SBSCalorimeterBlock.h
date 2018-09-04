#ifndef ROOT_SBSCalorimeterBlock
#define ROOT_SBSCalorimeterBlock

////////////////////////////////////////////////////////////////////////////
//                                                                        //
// SBSCalorimeterBlock                                                       //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include <TObject.h>
#include "SBSCalorimeterBlockData.h"


///////////////////////////////////////////////////////////////////////////////
// Calorimeter Blocks by default have only single-valued ADC data
class SBSCalorimeterBlock : public SBSCalorimeterBlockData::ADC,
  public TObject {

public:
  SBSCalorimeterBlock() {};
  SBSCalorimeterBlock(Float_t x, Float_t y, Float_t z,
      Int_t row, Int_t col, Int_t layer, Float_t adc_ped, Float_t adc_gain);
  virtual ~SBSCalorimeterBlock() {}

  // Getters
  Float_t GetX()     const { return fX; }
  Float_t GetY()     const { return fY; }
  Float_t GetZ()     const { return fZ; }
  Int_t   GetRow()   const { return fRow; }
  Int_t   GetCol()   const { return fCol; }
  Int_t   GetLayer() const { return fLayer; }
  Int_t   GetStat()  const { return fStat; }

  // Setters
  void SetX(Float_t var)    { fX = var; }
  void SetY(Float_t var)    { fY = var; }
  void SetZ(Float_t var)    { fZ = var; }
  void SetRow(Int_t var)    { fRow = var; }
  void SetCol(Int_t var)    { fCol = var; }
  void SetLayer(Int_t var)  { fLayer = var; }
  void SetStat(Int_t var)   { fStat = var; }

  // Sub-classes may want a more comprehensive clear
  virtual void ClearEvent();

  // Check if this block has any data
  virtual Bool_t HasData() { return HasADCData(); };

protected:
  Float_t fX;       ///< relative x position of the center
  Float_t fY;       ///< relative y position of the center
  Float_t fZ;       ///< relative z position of the center

  Int_t   fRow;     ///< Row of the block
  Int_t   fCol;     ///< Column of the block
  Int_t   fLayer;   ///< Layer of the block
  Int_t   fStat;    ///< Status: 0: not seen, 1: seen, 2: local max

  ClassDef(SBSCalorimeterBlock,1) ///< Generic shower block class single-valued
};

///////////////////////////////////////////////////////////////////////////////
// A calorimeter block with single-value and TDC data
class SBSCalorimeterBlockTDC : public SBSCalorimeterBlock,
  public SBSCalorimeterBlockData::TDC {

public:
  SBSCalorimeterBlockTDC() {};
  SBSCalorimeterBlockTDC(Float_t x, Float_t y, Float_t z,
      Int_t row, Int_t col, Int_t layer, Float_t adc_ped, Float_t adc_gain,
      Float_t tdc_offset, Float_t tdc_cal);
  virtual ~SBSCalorimeterBlockTDC() {}

  // Sub-classes may want a more comprehensive clear
  virtual void ClearEvent();


  // Check if this block has any data
  virtual Bool_t HasData() { return HasTDCData(); };

  ClassDef(SBSCalorimeterBlockTDC,1) ///< Single-valued ADC with TDC
};

///////////////////////////////////////////////////////////////////////////////
// A calorimeter block with multivalued-adc (i.e. with samples)
class SBSCalorimeterBlockSamples : public SBSCalorimeterBlock,
  public SBSCalorimeterBlockData::Samples {

public:
  SBSCalorimeterBlockSamples() {};
  SBSCalorimeterBlockSamples(Float_t x, Float_t y, Float_t z,
      Int_t row, Int_t col, Int_t layer, Float_t adc_ped, Float_t adc_gain,
      Float_t adc_ped_mult);
  virtual ~SBSCalorimeterBlockSamples() {}

  // Sub-classes may want a more comprehensive clear
  virtual void ClearEvent();

  // Check if this block has any data
  virtual Bool_t HasData() { return (SBSCalorimeterBlock::HasData()
      ||HasSamplesData()); };

  ClassDef(SBSCalorimeterBlockSamples,1) ///< Single-valued ADC with 
};

///////////////////////////////////////////////////////////////////////////////
// A calorimeter block with multivalued-adc (i.e. with samples) and TDC data
class SBSCalorimeterBlockSamplesTDC : public SBSCalorimeterBlockSamples,
  public SBSCalorimeterBlockData::TDC {

public:
  SBSCalorimeterBlockSamplesTDC() {};
  SBSCalorimeterBlockSamplesTDC(Float_t x, Float_t y, Float_t z,
      Int_t row, Int_t col, Int_t layer, Float_t adc_ped, Float_t adc_gain,
      Float_t adc_ped_mult, Float_t tdc_offset, Float_t tdc_cal);
  virtual ~SBSCalorimeterBlockSamplesTDC() {}

  // Sub-classes may want a more comprehensive clear
  virtual void ClearEvent();

  // Check if this block has any data
  virtual Bool_t HasData() { return (SBSCalorimeterBlockSamples::HasData()
      ||HasTDCData()); };
  ClassDef(SBSCalorimeterBlockSamplesTDC,1) ///< Multi-valued ADC with TDC
};

#endif
