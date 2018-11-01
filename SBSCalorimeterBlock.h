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
// Generic Calorimeter Block
class SBSCalorimeterBlock : public TObject {

public:
  SBSCalorimeterBlock() : fADC(0), fTDC(0), fSamples(0) {};
  SBSCalorimeterBlock(Float_t x, Float_t y, Float_t z,
      Int_t row, Int_t col, Int_t layer);
  virtual ~SBSCalorimeterBlock() {}

  // Getters
  Float_t GetX()     const { return fX; }
  Float_t GetY()     const { return fY; }
  Float_t GetZ()     const { return fZ; }
  Float_t GetE()     const { return fE; }
  Int_t   GetRow()   const { return fRow; }
  Int_t   GetCol()   const { return fCol; }
  Int_t   GetLayer() const { return fLayer; }
  Int_t   GetStat()  const { return fStat; }
  virtual SBSCalorimeterBlockData::ADC* ADC()         { return fADC; }
  virtual SBSCalorimeterBlockData::TDC* TDC()         { return fTDC; }
  virtual SBSCalorimeterBlockData::Samples* Samples() { return fSamples; }

  // Setters
  void SetX(Float_t var)    { fX = var; }
  void SetY(Float_t var)    { fY = var; }
  void SetZ(Float_t var)    { fZ = var; }
  void SetE(Float_t var)    { fE = var; }
  void SetRow(Int_t var)    { fRow = var; }
  void SetCol(Int_t var)    { fCol = var; }
  void SetLayer(Int_t var)  { fLayer = var; }
  void SetStat(Int_t var)   { fStat = var; }
  void SetADC(Float_t ped, Float_t gain);
  void SetTDC(Float_t offset, Float_t cal);
  void SetSamples(Float_t ped, Float_t gain);

  // Sub-classes may want a more comprehensive clear
  virtual void ClearEvent();

  // Coarse process this event for this block
  virtual void CoarseProcess();

  // Check if this block has any data
  virtual Bool_t HasData();

protected:
  Float_t fX;       ///< relative x position of the center
  Float_t fY;       ///< relative y position of the center
  Float_t fZ;       ///< relative z position of the center
  Float_t fE;       ///< calibrated energy of event in this block

  Int_t   fRow;     ///< Row of the block
  Int_t   fCol;     ///< Column of the block
  Int_t   fLayer;   ///< Layer of the block
  Int_t   fStat;    ///< Status: 0: not seen, 1: seen, 2: local max

  SBSCalorimeterBlockData::ADC *fADC;
  SBSCalorimeterBlockData::TDC *fTDC;
  SBSCalorimeterBlockData::Samples *fSamples;

  ClassDef(SBSCalorimeterBlock,1) ///< Generic shower block class (no data)
};

#endif
