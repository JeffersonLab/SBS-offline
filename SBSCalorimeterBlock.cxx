///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// SBSCalorimeterBlock                                                          //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include "SBSCalorimeterBlock.h"

//const double SBSCalorimeterBlock::kBig = 1.e15;

ClassImp(SBSCalorimeterBlock);
ClassImp(SBSCalorimeterBlockTDC);
ClassImp(SBSCalorimeterBlockSamples);
ClassImp(SBSCalorimeterBlockSamplesTDC);

///////////////////////////////////////////////////////////////////////////////
// Constructor for single-value ADC
SBSCalorimeterBlock::SBSCalorimeterBlock(Float_t x, Float_t y,
    Float_t z, Int_t row, Int_t col, Int_t layer,Float_t adc_ped,
    Float_t adc_gain) : SBSCalorimeterBlockData::ADC(adc_ped,adc_gain),
  fX(x), fY(y), fZ(z), fRow(row), fCol(col), fLayer(layer), fStat(0)
{
}

///////////////////////////////////////////////////////////////////////////////
// Constructor for single-value ADC + TDC
SBSCalorimeterBlockTDC::SBSCalorimeterBlockTDC(Float_t x, Float_t y, Float_t z,
    Int_t row, Int_t col, Int_t layer, Float_t adc_ped, Float_t adc_gain,
    Float_t tdc_offset, Float_t tdc_cal) :
  SBSCalorimeterBlock(x,y,z,row,col,layer,adc_ped,adc_gain),
  SBSCalorimeterBlockData::TDC(tdc_offset,tdc_cal)
{
}

///////////////////////////////////////////////////////////////////////////////
// Constructor for multi-valued ADC
SBSCalorimeterBlockSamples::SBSCalorimeterBlockSamples(Float_t x,
    Float_t y, Float_t z, Int_t row, Int_t col, Int_t layer, Float_t adc_ped,
    Float_t adc_gain, Float_t adc_ped_mult) :
  SBSCalorimeterBlock(x,y,z,row,col,layer,adc_ped,adc_gain),
  SBSCalorimeterBlockData::Samples(adc_ped*adc_ped_mult,adc_gain)
{
}

///////////////////////////////////////////////////////////////////////////////
// Constructor for multi-valued ADC + TDC
SBSCalorimeterBlockSamplesTDC::SBSCalorimeterBlockSamplesTDC(Float_t x,
    Float_t y, Float_t z, Int_t row, Int_t col, Int_t layer, Float_t adc_ped,
    Float_t adc_gain, Float_t adc_ped_mult, Float_t tdc_offset,
    Float_t tdc_cal) :
  SBSCalorimeterBlockSamples(x,y,z,row,col,layer,adc_ped,adc_gain,adc_ped_mult),
  SBSCalorimeterBlockData::TDC(tdc_offset,tdc_cal)
{
}

///////////////////////////////////////////////////////////////////////////////
// Clear event from single-valued ADC
void SBSCalorimeterBlock::ClearEvent()
{
  ClearADC();
  fStat = 0; // Reset status to 0, unseen
}

///////////////////////////////////////////////////////////////////////////////
// Clear event from single-valued ADC + TDC
void SBSCalorimeterBlockTDC::ClearEvent()
{
  SBSCalorimeterBlock::ClearEvent();
  ClearTDC();
}

///////////////////////////////////////////////////////////////////////////////
// Clear event from multi-valued ADC
void SBSCalorimeterBlockSamples::ClearEvent()
{
  SBSCalorimeterBlock::ClearEvent();
  ClearSamples();
}

///////////////////////////////////////////////////////////////////////////////
// Clear event from multi-valued ADC + TDC
void SBSCalorimeterBlockSamplesTDC::ClearEvent()
{
  SBSCalorimeterBlockSamples::ClearEvent();
  ClearTDC();
}
