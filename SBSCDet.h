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
#include "SBSCDet_Hit.h"
#include "TBits.h"
#include "TClonesArray.h"
#include <cstdint>
#include <map>

class THaTrack;
class THaBenhmark;

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

  SBSCDet_Hit*		GetHit(Int_t i) const
  {return (SBSCDet_Hit*)fHits->At(i); }

  Int_t 		GetNumHits() const
    { return fHits->GetLast()+1; }
  
  //Int_t                GetNumClusters() const
  //  { return fClusters->GetLast()+1; }

protected:

  TClonesArray*		fHits;		// Array of hits for each event
  Double_t	fHit_tmin;
  Double_t 	fHit_tmax;
  Double_t	fHit_totmin;
  Double_t 	fHit_totmax;

  ClassDef(SBSCDet,5)  // Describes scintillator plane with F1TDC as a detector
};

#endif
