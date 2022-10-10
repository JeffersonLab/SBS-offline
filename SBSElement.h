#ifndef ROOT_SBSElement
#define ROOT_SBSElement

////////////////////////////////////////////////////////////////////////////
//                                                                        //
// SBSElement                                                       //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include <TObject.h>
#include "SBSData.h"


///////////////////////////////////////////////////////////////////////////////

class SBSElement : public TObject {

public:
  SBSElement() : fADC(nullptr), fTDC(nullptr), fWaveform(nullptr) {};
  SBSElement(Double_t x, Double_t y, Double_t z,
      Int_t row, Int_t col, Int_t layer, Int_t id = 0);
  virtual ~SBSElement();

  // Getters
  Double_t GetX()     const { return fX; }
  Double_t GetY()     const { return fY; }
  Double_t GetZ()     const { return fZ; }
  Double_t GetE()     const { return fE; }
  Double_t GetAgain()     const { return fAgain; }
  Double_t GetAtime()     const { return fAtime; }
  Double_t GetTDCtime()     const { return fTDCtime; }
  Int_t   GetRow()   const { return fRow; }
  Int_t   GetCol()   const { return fCol; }
  Int_t   GetLayer() const { return fLayer; }
  Int_t   GetStat()  const { return fStat; }
  Int_t   GetID()    const { return fID; }
  virtual SBSData::ADC* ADC()         { return fADC; }
  virtual SBSData::TDC* TDC()         { return fTDC; }
  virtual SBSData::Waveform* Waveform() { return fWaveform; }

  // Setters
  void SetX(Double_t var)    { fX = var; }
  void SetY(Double_t var)    { fY = var; }
  void SetZ(Double_t var)    { fZ = var; }
  void SetE(Double_t var)    { fE = var; }
  void SetAgain(Double_t var)    { fAgain = var; }
  void SetAtime(Double_t var)    { fAtime = var; }
  void SetTDCtime(Double_t var)    { fTDCtime = var; }
  void SetRow(Int_t var)    { fRow = var; }
  void SetCol(Int_t var)    { fCol = var; }
  void SetLayer(Int_t var)  { fLayer = var; }
  void SetStat(Int_t var)   { fStat = var; }
  void SetID(Int_t var)     { fID = var; }
  void SetADC(Double_t ped, Double_t gain);
  void SetTDC(Double_t offset, Double_t cal, Double_t GoodTimeCut);
  void SetWaveform(Double_t ped, Double_t gain,Double_t ChanToMv,Double_t adc_timecut);

  // Sub-classes may want a more comprehensive clear
  virtual void Clear( Option_t* opt="" );


  // Check if this block has any data
  virtual Bool_t HasData();
  virtual Bool_t HasADCData();

protected:
  Double_t fX;       ///< relative x position of the center
  Double_t fY;       ///< relative y position of the center
  Double_t fZ;       ///< relative z position of the center
  Double_t fE;       ///< calibrated energy of event in this block
  Double_t fAgain;   ///< ADC gain coefficient (GeV/pC)
  Double_t fAtime;       ///< ADC time of event in this block
  Double_t fTDCtime;       ///< TDC time of event in this block

  Int_t   fRow;     ///< Row of the block
  Int_t   fCol;     ///< Column of the block
  Int_t   fLayer;   ///< Layer of the block
  Int_t   fStat;    ///< Status: 0: not seen, 1: seen, 2: local max
  Int_t   fID;      ///< a logical number to this element

  SBSData::ADC *fADC; //< All ADC hits
  SBSData::TDC *fTDC; //< All TDC hits
  SBSData::Waveform *fWaveform;

  ClassDef(SBSElement,1) ///< Generic shower block class (no data)
};

#endif
