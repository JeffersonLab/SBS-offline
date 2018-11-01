///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// SBSCalorimeterBlock                                                       //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include "SBSCalorimeterBlock.h"

//const double SBSCalorimeterBlock::kBig = 1.e15;

ClassImp(SBSCalorimeterBlock);

///////////////////////////////////////////////////////////////////////////////
// Constructor for generic CalorimeterBlock (no-data)
SBSCalorimeterBlock::SBSCalorimeterBlock(Float_t x, Float_t y,
    Float_t z, Int_t row, Int_t col, Int_t layer) :
  fX(x), fY(y), fZ(z), fRow(row), fCol(col), fLayer(layer), fStat(0),
  fADC(0), fTDC(0), fSamples(0)
{
}


///////////////////////////////////////////////////////////////////////////////
// Check if this block has any data
Bool_t SBSCalorimeterBlock::HasData()
{
  return ( ( fADC && fADC->HasData() ) || ( fTDC && fTDC->HasData() ) ||
      ( fSamples && fSamples->HasData() ) );
}

///////////////////////////////////////////////////////////////////////////////
// Clear event from generic CalorimeterBlock (with no data)
void SBSCalorimeterBlock::ClearEvent()
{
  fE = 0; // Reset calibrated energy for given event
  fStat = 0; // Reset status to 0, unseen
  if(fADC)
    fADC->Clear();
  if(fTDC)
    fTDC->Clear();
  if(fSamples)
    fSamples->Clear();
}

///////////////////////////////////////////////////////////////////////////////
// Create a single-valued ADC data structure
void SBSCalorimeterBlock::SetADC(Float_t ped, Float_t gain)
{
  if(fADC)
    delete fADC;
  fADC = new SBSCalorimeterBlockData::ADC(ped,gain);
}

///////////////////////////////////////////////////////////////////////////////
// Create a TDC data structure
void SBSCalorimeterBlock::SetTDC(Float_t offset, Float_t cal)
{
  if(fTDC)
    delete fTDC;
  fTDC = new SBSCalorimeterBlockData::TDC(offset,cal);
}

///////////////////////////////////////////////////////////////////////////////
// Create a multi-valued ADC data structure
void SBSCalorimeterBlock::SetSamples(Float_t ped, Float_t gain)
{
  if(fSamples)
    delete fSamples;
  fSamples = new SBSCalorimeterBlockData::Samples(ped,gain);
}

///////////////////////////////////////////////////////////////////////////////
// Coarse process this event
void SBSCalorimeterBlock::CoarseProcess()
{
  // Compute the energy
  if(fSamples) {
    fE = fSamples->GetDataSumCal();
  } else if ( fADC) {
    fE = fADC->GetDataCal();
  } else {
    fE = 0;
  }
}
