///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// SBSElement                                                                //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include <iostream>
#include "SBSElement.h"
#include "DataType.h"

ClassImp(SBSElement);

///////////////////////////////////////////////////////////////////////////////
// Constructor for generic Element (no-data)
SBSElement::SBSElement(Double_t x, Double_t y,
    Double_t z, Int_t row, Int_t col, Int_t layer, Int_t id) :
  fX(x), fY(y), fZ(z), fE(0), fAgain(0), fAtime(kBig), fTDCtime(kBig), fRow(row), fCol(col), fLayer(layer),
  fStat(0), fID(id), fADC(nullptr), fTDC(nullptr), fWaveform(nullptr)
{
}

///////////////////////////////////////////////////////////////////////////////
SBSElement::~SBSElement()
{
  delete fADC;
  delete fTDC;
  delete fWaveform;
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
void SBSElement::Clear( Option_t* opt )
{
  fE = 0; // Reset calibrated energy for given event
  fAgain = 0.;
  fStat = 0; // Reset status to 0, unseen
  if(fADC)
    fADC->Clear();
  if(fTDC)
    fTDC->Clear();
  if(fWaveform)
    fWaveform->Clear();
  fAtime = kBig;
  fTDCtime = kBig;
}

///////////////////////////////////////////////////////////////////////////////
// Create a single-valued ADC data structure
void SBSElement::SetADC(Double_t ped, Double_t gain)
{
  delete fADC;
  fADC = new SBSData::ADC(ped,gain);
}

///////////////////////////////////////////////////////////////////////////////
// Create a TDC data structure
void SBSElement::SetTDC(Double_t offset, Double_t cal, Double_t GoodTimeCut)
{
  delete fTDC;
  fTDC = new SBSData::TDC(offset,cal,GoodTimeCut);
}

///////////////////////////////////////////////////////////////////////////////
// Create a multi-valued ADC data structure
void SBSElement::SetWaveform(Double_t ped, Double_t gain, Double_t ChanToMv, Double_t adc_timecut)
{
  delete fWaveform;
  fWaveform = new SBSData::Waveform(ped,gain,ChanToMv,adc_timecut);
}

