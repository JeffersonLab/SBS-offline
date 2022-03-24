//////////////////////////////////////////////////////////////////////////
//
// SBSRPFarSideHodo
//
// General detector with TDC information
//
//////////////////////////////////////////////////////////////////////////

#ifndef SBSRPFarSideHodo_h
#define SBSRPFarSideHodo_h

#include "SBSGenericDetector.h"

class SBSRPFarSideHodo : public SBSGenericDetector {
public:
  SBSRPFarSideHodo( const char* name, const char* description, 
      THaApparatus* apparatus=0);

  virtual ~SBSRPFarSideHodo();

  // Overwritten derived functions
  virtual Int_t   ReadDatabase( const TDatime& date );
  virtual Int_t   DefineVariables( EMode mode = kDefine );
  virtual Int_t   CoarseProcess( TClonesArray& tracks );
  virtual Int_t   FineProcess( TClonesArray& tracks );
  virtual Int_t   FindGoodHit(SBSElement *element);
  virtual void    Clear( Option_t* opt="" );

  ClassDef(SBSRPFarSideHodo,5)  // Describes scintillator plane with F1TDC as a detector
};

#endif
