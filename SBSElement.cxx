///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// SBSElement                                                                //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include <iostream>
#include "SBSElement.h"

ClassImp(SBSElement);

///////////////////////////////////////////////////////////////////////////////
// Constructor for generic Element (no-data)
SBSElement::SBSElement(Float_t x, Float_t y,
    Float_t z, Int_t row, Int_t col, Int_t layer, Int_t id) :
  fX(x), fY(y), fZ(z), fRow(row), fCol(col), fLayer(layer), fStat(0), fID(id),
  fADC(0), fTDC(0), fWaveform(0)
{
}


///////////////////////////////////////////////////////////////////////////////
// Check if this block has any ADC data
Bool_t SBSElement::HasData()
{
  return ( ( fADC && fADC->HasData() ) || ( fTDC && fTDC->HasData() ) || ( fWaveform && fWaveform->HasData() ) );
}

// Check if this block has any ADC data
Bool_t SBSElement::HasADCData()
{
  return ( ( fADC && fADC->HasData() ) || ( fWaveform && fWaveform->HasData() ) );
}


///////////////////////////////////////////////////////////////////////////////
// Clear event from generic Element (with no data)
void SBSElement::ClearEvent()
{
  fE = 0; // Reset calibrated energy for given event
  fStat = 0; // Reset status to 0, unseen
  if(fADC)
    fADC->Clear();
  if(fTDC)
    fTDC->Clear();
  if(fWaveform)
    fWaveform->Clear();
}

///////////////////////////////////////////////////////////////////////////////
// Create a single-valued ADC data structure
void SBSElement::SetADC(Float_t ped, Float_t gain)
{
  if(fADC)
    delete fADC;
  fADC = new SBSData::ADC(ped,gain);
}

///////////////////////////////////////////////////////////////////////////////
// Create a TDC data structure
void SBSElement::SetTDC(Float_t offset, Float_t cal, Float_t GoodTimeCut)
{
  if(fTDC)
    delete fTDC;
  fTDC = new SBSData::TDC(offset,cal,GoodTimeCut);
}

///////////////////////////////////////////////////////////////////////////////
// Create a multi-valued ADC data structure
void SBSElement::SetWaveform(Float_t ped, Float_t gain, Float_t ChanToMv, Float_t adc_timecut)
{
  if(fWaveform)
    delete fWaveform;
  fWaveform = new SBSData::Waveform(ped,gain,ChanToMv,adc_timecut);
}

