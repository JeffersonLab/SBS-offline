//////////////////////////////////////////////////////////////////////////
//
// SBSRPBeamSideHodo
//
// General detector with TDC information
//
//////////////////////////////////////////////////////////////////////////

#ifndef SBSRPBeamSideHodo_h
#define SBSRPBeamSideHodo_h

#include "SBSGenericDetector.h"

class SBSRPBeamSideHodo : public SBSGenericDetector {
public:
  SBSRPBeamSideHodo( const char* name, const char* description, 
      THaApparatus* apparatus=0);

  virtual ~SBSRPBeamSideHodo();

  // Overwritten derived functions
  virtual Int_t   ReadDatabase( const TDatime& date );
  virtual Int_t   DefineVariables( EMode mode = kDefine );
  virtual Int_t   CoarseProcess( TClonesArray& tracks );
  virtual Int_t   FineProcess( TClonesArray& tracks );
  virtual Int_t   FindGoodHit(SBSElement *element);
  virtual void    Clear( Option_t* opt="" );

  ClassDef(SBSRPBeamSideHodo,5)  // Describes scintillator plane with F1TDC as a detector
};

#endif
