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
#include "SBSTimingHodoscopePMT.h"
#include "SBSTimingHodoscopeBar.h"

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
  virtual void    ClearEvent();
  // new functions
  Int_t   ConstructHodoscope();
  Float_t TimeWalk(Float_t time);

 protected:
  Float_t AttLength = 1.0/3.8;//380cm for ej200. units per m?
  /* speed of light in scint? what is n for ej200 */
  // n = 1.58 => n=c/v =>v=c/n
  Float_t n=1.58;
  Float_t vScint = 2.9979e8/n;

  /* std::vector<SBSTimingHodoscopePMT> fPMTMap; */
  std::vector<SBSTimingHodoscopePMT*> fPMTMapL;
  std::vector<SBSTimingHodoscopePMT*> fPMTMapR;
  std::vector<SBSTimingHodoscopeBar*> fBars;
  Int_t fTDCBarOffset;
  Int_t fADCBarOffset;
  Int_t fTDCRefLeL;
  Int_t fTDCRefLeR;

  std::vector<Int_t>   fGoodBarIDsTDC;
  std::vector<Float_t> fGoodBarTDCmean;
  std::vector<Float_t> fGoodBarTDCdiff;
  std::vector<Float_t> fGoodBarTDCvpos;
  std::vector<Float_t> fGoodBarTDCpos;
  std::vector<Float_t> fGoodBarTDCLle;
  std::vector<Float_t> fGoodBarTDCLleW;
  std::vector<Float_t> fGoodBarTDCLte;
  std::vector<Float_t> fGoodBarTDCLteW;
  std::vector<Float_t> fGoodBarTDCLtot;
  std::vector<Float_t> fGoodBarTDCLtotW;
  std::vector<Float_t> fGoodBarTDCRle;
  std::vector<Float_t> fGoodBarTDCRleW;
  std::vector<Float_t> fGoodBarTDCRte;
  std::vector<Float_t> fGoodBarTDCRteW;
  std::vector<Float_t> fGoodBarTDCRtot;
  std::vector<Float_t> fGoodBarTDCRtotW;

  std::vector<Int_t>   fGoodBarIDsADC;
  std::vector<Float_t> fGoodBarADCmean;
  std::vector<Float_t> fGoodBarADCLa;
  std::vector<Float_t> fGoodBarADCLap;
  std::vector<Float_t> fGoodBarADCLac;
  std::vector<Float_t> fGoodBarADCRa;
  std::vector<Float_t> fGoodBarADCRap;
  std::vector<Float_t> fGoodBarADCRac;

  // Mapping (see also fDetMap)
  UShort_t   fChanMapStart; ///< Starting number for block number (i.e. 0 or 1)
  // maps of time walk parameters - in row, col, lay
  std::vector<std::vector<std::vector<Float_t>>> fTimeWalkPar0;
  std::vector<std::vector<std::vector<Float_t>>> fTimeWalkPar1;

  ClassDef(SBSTimingHodoscope,5)  // Describes scintillator plane with F1TDC as a detector
};

#endif
