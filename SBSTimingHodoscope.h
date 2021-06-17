//////////////////////////////////////////////////////////////////////////
//
// SBSTimingHodoscope
//
// General detector with TDC information
//
//////////////////////////////////////////////////////////////////////////

#ifndef SBSTimingHodoscope_h
#define SBSTimingHodoscope_h

#include "SBSGenericDetector.h"

class SBSTimingHodoscope : public SBSGenericDetector {
public:
  SBSTimingHodoscope( const char* name, const char* description, 
      THaApparatus* apparatus=0);

  virtual ~SBSTimingHodoscope();

  // Overwritten derived functions
  virtual Int_t   ReadDatabase( const TDatime& date );
  virtual Int_t   DefineVariables( EMode mode = kDefine );
  virtual Int_t   CoarseProcess( TClonesArray& tracks );
  virtual Int_t   FineProcess( TClonesArray& tracks );
  virtual Int_t   FindGoodHit(SBSElement *element);
  virtual void ClearEvent();

  ClassDef(SBSTimingHodoscope,5)  // Describes scintillator plane with F1TDC as a detector
};

#endif
