#ifndef _SBSScalerHelicity_h_
#define _SBSScalerHelicity_h_

//////////////////////////////////////////////////////////////////////
//
// SBSScalerHelicity
//
// Based on the Podd class HallA/THaQWEAKHelicity
// Helicity of the beam - from QWEAK electronics in delayed mode
// 
////////////////////////////////////////////////////////////////////////

#include <vector>

#include "THaHelicityDet.h"
#include "SBSScalerHelicityReader.h"

#include "TROOT.h"
#include "TTree.h"
#include "TString.h"

class TH1F;

class SBSScalerHelicity : public THaHelicityDet, public SBSScalerHelicityReader {

   public:
      SBSScalerHelicity( const char* name, const char* description,
	    THaApparatus* a = nullptr );
      SBSScalerHelicity();
      virtual ~SBSScalerHelicity();

      virtual Int_t  Begin( THaRunBase* r=nullptr );
      virtual void   Clear( Option_t* opt = "" );
      virtual Int_t  Decode( const THaEvData& evdata );
      virtual Int_t  End( THaRunBase* r=nullptr );
      virtual void   SetDebug( Int_t level );
      virtual Bool_t HelicityValid() const { return fValidHel; }

      void PrintEvent( UInt_t evtnum );

      void SetVerbosity(int v) { fVerbosity = v; } 

      void SetHelicityDelay(Int_t delay){
	if (delay<0) {
	  std::cerr << "****  Helicity Delay cannot be negative.  Force to zero." 
		    << std::endl;
	  fHelicityDelay = 0;
	}
	else if (delay>64) {
	  std::cerr << "****  Helicity Delay cannot be greater than 64.  Force to zero." 
		    << std::endl;
	} else {
	  fHelicityDelay = delay;
	}
      }


   protected:
      virtual void  FillHisto();
      void  CheckTIRvsRing( UInt_t eventnumber );
      void  LoadHelicity( UInt_t eventnumber );
      UInt_t RanBit30( UInt_t& ranseed );
      THaHelicityDet::EHelicity SetHelicity( UInt_t polarity, UInt_t phase);

      // variables that need to be read from the database
      UInt_t fOffsetTIRvsRing;
      // Offset between the ring reported value and the TIR reported value
      UInt_t fQWEAKDelay;
      // delay of helicity (# windows)
      UInt_t fMAXBIT;
      //number of bit in the pseudo random helicity generator
      std::vector<Int_t> fPatternSequence; // sequence of 0 and 1 in the pattern
      UInt_t fQWEAKNPattern; // maximum of event in the pattern
      Bool_t HWPIN;

      Int_t fHelicityDelay;

      Int_t fQrt;
      Int_t fTSettle;
      Bool_t fValidHel;

      UInt_t fHelicityLastTIR;
      UInt_t fPatternLastTIR;
      void SetErrorCode(Int_t error);
      Double_t fErrorCode;

      // Ring related data
      // UInt_t fScalerCumulative[32];
      UInt_t fScalerCumulativePlus[32];
      UInt_t fScalerCumulativeMinus[32];
      // 64-bit signed integer versions
      Long64_t fScalerCumulative[32];  // 64-bit signed integer 
      // Long64_t fScalerCumulativePlus[32];
      // Long64_t fScalerCumulativeMinus[32];

      UInt_t fRingFinalQrtHel;
      UInt_t fRingFinalEvtNum;
      UInt_t fRingFinalPatNum;
      UInt_t fRingFinalSeed;
      UInt_t fRingPattPhase;
      UInt_t fRingHelicitySum;

      Long_t fTimeStampYield;
      Long_t fTimeStampDiff;
      Long_t fScalerYield[32];
      Long_t fScalerDiff[32];


      UInt_t fRing_NSeed; //number of event collected for seed
      UInt_t fRingU3plus, fRingU3minus;
      UInt_t fRingT3plus, fRingT3minus;
      UInt_t fRingT5plus, fRingT5minus;
      UInt_t fRingT10plus, fRingT10minus;
      UInt_t fRingTimeplus, fRingTimeminus;
      UInt_t fRingSeed_reported;
      UInt_t fRingSeed_actual;
      UInt_t fRingPhase_reported;
      UInt_t fRing_reported_polarity;
      UInt_t fRing_actual_polarity;

      UInt_t fEvtype; // Current CODA event type

      Int_t fVerbosity;

      // tree to write data to 
      TTree *fHelScalerTree;
      // branch variables
      UInt_t fBranch_seed;
      Double_t fBranch_errCode;

      static const Int_t NHIST = 2;
      std::vector<TH1F*> fHisto;

      virtual Int_t DefineVariables( EMode mode = kDefine );
      virtual Int_t ReadDatabase( const TDatime& date );

      ClassDef(SBSScalerHelicity,0)   // Beam helicity from QWEAK electronics in delayed mode

};

#endif
