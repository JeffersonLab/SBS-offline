#ifndef SBSRaster_H 
#define SBSRaster_H

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// SBSRaster                                                                 //
//  Class for a beam raster device measuring two magnet currents              //
//  which are proportional to the horizontal/vertical beam displacement.      //
//  The two planes are assumed to be decoupled.                               //
//  There is no phase shift between the current and the actual beam position  //
//  Inherited from THaRaster, make it to decode data from FADC                // 
//  only GetData modified                                                     // 
//  Original author: hanjie@jlab.org                                          // 
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include <iostream>

#include "THaEvData.h"
#include "THaDetMap.h"
#include "THaBeamDet.h"

#include "VarDef.h"
#include "VarType.h"

#include "TMath.h"
#include "TVectorT.h"

class SBSRaster : public THaBeamDet {

public:
  explicit SBSRaster( const char* name, const char* description = "",
                      THaApparatus* a = nullptr );
  virtual ~SBSRaster();

  virtual void       Clear( Option_t* ="" );
  virtual Int_t      Decode( const THaEvData& );
  virtual Int_t      Process();

  virtual TVector3 GetPosition()  const { return fPosition[2]; }
  virtual TVector3 GetDirection() const { return fDirection; }

  // As soon as someone finds a better solution, the following lines should be
  // changed. It is ridiculus to have nine methods to get the the components
  // of the beam position at various locations, but I do not know how else to
  // get them into histograms except for writing my own event loop

  Double_t GetRawPosX() { return fRawPos(0); }
  Double_t GetRawPosY() { return fRawPos(1); }
  Double_t GetRawSlopeX() { return fRawSlope(0); }
  Double_t GetRawSlopeY() { return fRawSlope(1); }

  Double_t GetPosBPMAX() { return fPosition[0](0); }
  Double_t GetPosBPMAY() { return fPosition[0](1); }
  Double_t GetPosBPMAZ() { return fPosition[0](2); }

  Double_t GetPosBPMBX() { return fPosition[1](0); }
  Double_t GetPosBPMBY() { return fPosition[1](1); }
  Double_t GetPosBPMBZ() { return fPosition[1](2); }

  Double_t GetPosTarX() { return fPosition[2](0); }
  Double_t GetPosTarY() { return fPosition[2](1); }
  Double_t GetPosTarZ() { return fPosition[2](2); }

protected:

  virtual Int_t  ReadDatabase( const TDatime& date );
  virtual Int_t  DefineVariables( EMode mode = kDefine );

  //  SBSRaster() {}
  //  SBSRaster( const SBSRaster& ) {}
  //  SBSRaster& operator=( const SBSRaster& ) { return *this; }

  typedef TVectorT<Double_t> TVectorD;

  TVectorD  fRawPos;        // current in Raster ADCs for position
  TVectorD  fRawSlope;      // current in Raster ADCs for the derivative

  TVector3  fPosition[3];   // Beam position at 1st, 2nd BPM or at the target (meters)
  TVector3  fDirection;  // Beam angle at the target (meters)

  TMatrix   fRaw2Pos[3];
  TVector3  fPosOff[3];

  TVectorD  fRasterFreq;
  TVectorD  fSlopePedestal;
  TVectorD  fRasterPedestal;

  Int_t fNfired;

  ClassDef(SBSRaster,0)   // Generic Raster class 
};

#endif
