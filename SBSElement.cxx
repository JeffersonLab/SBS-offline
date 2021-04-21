///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// SBSElement                                                                //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include "SBSElement.h"

ClassImp(SBSElement);

///////////////////////////////////////////////////////////////////////////////
// Constructor for generic Element (no-data)
SBSElement::SBSElement(Float_t x, Float_t y,
    Float_t z, Int_t row, Int_t col, Int_t layer, Int_t id) :
  fX(x), fY(y), fZ(z), fRow(row), fCol(col), fLayer(layer), fStat(0), fID(id),
  fADC(0), fTDC(0), fWaveform(0), fCoarseProcessed(0)
{
}


///////////////////////////////////////////////////////////////////////////////
// Check if this block has any data
Bool_t SBSElement::HasData()
{
  return ( ( fADC && fADC->HasData() ) || ( fTDC && fTDC->HasData() ) ||
      ( fWaveform && fWaveform->HasData() ) );
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
  fCoarseProcessed = false;
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
void SBSElement::SetTDC(Float_t offset, Float_t cal)
{
  if(fTDC)
    delete fTDC;
  fTDC = new SBSData::TDC(offset,cal);
}

///////////////////////////////////////////////////////////////////////////////
// Create a multi-valued ADC data structure
void SBSElement::SetWaveform(Float_t ped, Float_t gain)
{
  if(fWaveform)
    delete fWaveform;
  fWaveform = new SBSData::Waveform(ped,gain);
}

///////////////////////////////////////////////////////////////////////////////
// Coarse process this event
void SBSElement::CoarseProcess()
{

  if(fCoarseProcessed)
    return;

  // Compute the energy
  if(fWaveform) {
    // TODO: Implement the same logic as the FADC firmware logic.
    fE = fWaveform->GetIntegral().val;
  } else if ( fADC) {
    // For ADCs with multiple hits, one should mark the "good" hit like so
    // with fADC->SetGoodHit( idx );
    // TODO: Find out how to determine the "good" hit. For now, take the first one.
    fE = fADC->GetIntegral(0).val;
  } else {
    fE = 0;
  }

  if(fE < 0 ) { // Do not allowe negative energy!
    fE = 0.0;
  }

  fCoarseProcessed = true;

  // For the TDCs one should mark the "good" hit
  // with fTDC->SetGoodHit( idx );
}
