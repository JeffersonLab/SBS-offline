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

      UInt_t RanBit30( UInt_t& ranseed );

      // variables that need to be read from the database

      Int_t fMAXBIT; // Number of bits in the pseudorandom helicity generator
      std::vector<Int_t> fPatternSequence; // Sequence of +1 and -1 in the pattern
      Int_t fHelicityDelay; // Helicity delay in # of patterns


      Int_t fTriggerCounter; //  Reports the physics event number of the scaler event
  
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

      Int_t  fSeedbits;
      UInt_t fRingSeed_reported;
      UInt_t fRingSeed_actual;
      UInt_t fRingPhase_reported;
      UInt_t fRing_reported_polarity;
      UInt_t fRing_actual_polarity;

      UInt_t fEvtype; // Current CODA event type

      Int_t fVerbosity;

      // tree to write data to 
      TTree *fHelScalerTree;
      TTree *fAsymmetryTree;
      // branch variables
      UInt_t fBranch_seed;

      static const Int_t NHIST = 2;
      std::vector<TH1F*> fHisto;

      virtual Int_t DefineVariables( EMode mode = kDefine );
      virtual Int_t ReadDatabase( const TDatime& date );

      ClassDef(SBSScalerHelicity,0)   // Beam helicity from QWEAK electronics in delayed mode

};

#endif
