//////////////////////////////////////////////////////////////////////////
//
// SBSRPFarSideHodo
//
// General detector with TDC information
//
//////////////////////////////////////////////////////////////////////////

#ifndef SBSCDet_h
#define SBSCDet_h

#include "SBSGenericDetector.h"

class SBSCDet : public SBSGenericDetector {
public:
  SBSCDet( const char* name, const char* description, 
      THaApparatus* apparatus=0);

  virtual ~SBSCDet();

  // Overwritten derived functions
  virtual Int_t   ReadDatabase( const TDatime& date );
  virtual Int_t   DefineVariables( EMode mode = kDefine );
  virtual Int_t   Decode( const THaEvData& evdata );
  virtual Int_t   CoarseProcess( TClonesArray& tracks );
  virtual Int_t   FineProcess( TClonesArray& tracks );
  virtual Int_t   FindGoodHit(SBSElement *element);
  virtual void    Clear( Option_t* opt="" );

  ClassDef(SBSCDet,5)  // Describes scintillator plane with F1TDC as a detector
};

#endif
