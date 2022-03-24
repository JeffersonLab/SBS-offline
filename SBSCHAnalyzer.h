
//////////////////////////////////////////////////////////////////////////
//
// SBSCHAnalyzer
//
// General detector with TDC information
//
//////////////////////////////////////////////////////////////////////////
#ifndef SBSCHAnalyzer_h
#define SBSCHAnalyzer_h

#include "SBSGenericDetector.h"

class SBSCHAnalyzer : public SBSGenericDetector {
public:
  SBSCHAnalyzer( const char* name, const char* description, 
      THaApparatus* apparatus=0);

  virtual ~SBSCHAnalyzer();

  // Overwritten derived functions
  virtual Int_t   ReadDatabase( const TDatime& date );
  virtual Int_t   DefineVariables( EMode mode = kDefine );
  virtual Int_t   CoarseProcess( TClonesArray& tracks );
  virtual Int_t   FineProcess( TClonesArray& tracks );
  virtual Int_t   FindGoodHit(SBSElement *element);
  virtual void    Clear( Option_t* opt="" );

  ClassDef(SBSCHAnalyzer,5)  // Describes scintillator plane with F1TDC as a detector
};

#endif
