#ifndef _SBSScalerHelicityReader_h_
#define _SBSScalerHelicityReader_h_

////////////////////////////////////////////////////////////////////////
//
// SBSScalerHelicityReader
//
// Based on the Podd class HallA/THaQWEAKHelicityReader
//
//////////////////////////////////////////////////////////////////////////

#include "Rtypes.h"

class THaEvData;
class TDatime;
class TH1F;

class SBSScalerHelicityReader {

   public:
      SBSScalerHelicityReader();
      virtual ~SBSScalerHelicityReader();
      // virtual ~SBSScalerHelicityReader() = default;

      // when an event trigger the acquisition the reported helicity, T_settle and Pattern sync
      //signals are recorded by the an Input Register.
      UInt_t GetPatternTir() const { return fPatternTir; };
      UInt_t GetHelicityTir() const { return fHelicityTir; };
      UInt_t GetTSettleTir() const { return fTSettleTir; };
      UInt_t GetTimeStampTir() const { return fTimeStampTir; };
      // in parallel, these signals are recorded in the ring buffer of the helicity gated scaler
      // the maximum depth of this ring is for now 500.

      UInt_t GetRingDepth() const { return fIRing; }; // how many events are in the ring
      // I should have the depth of the ring defined in the data base
      // let call it kHelRingDepth for now
      void Print();

      void SetVerbosity(int v) { fVerbosity = v; } 

      class ROCinfo {
	 public:
	    ROCinfo() : roc(kMaxUInt), header(0), index(0) {}
	    Bool_t  valid() const;
	    UInt_t  roc;              // ROC to read out
	    UInt_t  header;           // Headers to search for (0 = ignore)
	    UInt_t  index;            // Index into buffer
      };

      class BANKinfo {
	 public:
	    BANKinfo() : roc(kMaxUInt), bankid(0), index(0), length(0) {};
	    Bool_t valid() const;
	    UInt_t roc;
	    UInt_t bankid;
	    UInt_t index;
	    UInt_t length;
      };

   protected:
      // Used by ReadDatabase
      enum EROC { kHel = 0, kTime, kRing, kROC3 };
      Int_t SetROCinfo( EROC which, UInt_t roc, UInt_t header, UInt_t index );

      const UInt_t* GetBankBuffer(const THaEvData&, SBSScalerHelicityReader::BANKinfo&);

      virtual void  Clear( Option_t* opt="" );
      virtual Int_t ReadData( const THaEvData& evdata );
      Int_t         ReadDatabase( const char* dbfilename, const char* prefix,
	    const TDatime& date, int debug_flag = 0 );
      virtual void FillHisto();
      void Begin();
      void End();

      //  These variables are the reported integrals from the FADCs
      UInt_t fFADCPatSync;
      UInt_t fFADCHelicity;
      UInt_t fFADCTSettle;
      UInt_t fFADCSpare;

      //  These variables are the reported integrals from the FADCs
      UInt_t fFADCPatSync_sbs;
      UInt_t fFADCHelicity_sbs;
      UInt_t fFADCTSettle_sbs;
      UInt_t fFADCSpare_sbs;

      UInt_t fFADCLow_min;  ///  The minimum value for FADC "low" helicity bit
      UInt_t fFADCLow_max;  ///  The maximum value for FADC "low" helicity bit
      UInt_t fFADCHigh_min; ///  The minimum value for FADC "high" helicity bit
      UInt_t fFADCHigh_max; ///  The maximum value for FADC "high" helicity bit

      UInt_t fFADCQrtHel;


      //  These variables should be removed.
      UInt_t fPatternTir;
      UInt_t fHelicityTir;
      UInt_t fTSettleTir;
      UInt_t fTimeStampTir;
      UInt_t fOldTimeStampTir;

      UInt_t fIRing;
      // the ring map is defined in 
      // http://www.jlab.org/~adaq/halog/html/1011_archive/101101100943.html
      // this map is valid after November 1rst
      enum { kHelRingDepth = 500 };
      UInt_t fScalerRing[kHelRingDepth][32]; // New for 2021
      UInt_t fHelicityRing[kHelRingDepth];
      UInt_t fPatternRing[kHelRingDepth];
      UInt_t fTimeStampRing[kHelRingDepth];

      UInt_t fT3Ring[kHelRingDepth];
      UInt_t fU3Ring[kHelRingDepth];
      UInt_t fT5Ring[kHelRingDepth];
      UInt_t fT10Ring[kHelRingDepth];

      //  New summary variables from the helicity CRL
      UInt_t fHelErrorCond;   //  Helicity ROC error condition
      UInt_t fNumEvents;      //  Number of helicity windows
      UInt_t fNumPatterns;    //  Number of patterns
      UInt_t fPatternPhase;   //  Current pattern phase
      UInt_t fSeedValue;      //  Seed value for current pattern
      UInt_t fPatternHel;     //  Current reported helicity at start of pattern
      UInt_t fEventPolarity;  //  Polarity of the current event with respect to start of pattern:  0x0==same as pattern polarity, 0x1==reversed polarity
      UInt_t fReportedQrtHel; //  Current "qrthel" value:  0x10==start of pattern, 0x1==Hel+

      // ROC information
      // If a header is zero the index is taken to be from the start of
      // the ROC (0 = first word of ROC), otherwise it's from the header
      // (0 = first word after header).
      ROCinfo  fROCinfo[kROC3+1]; // for now only work with one roc (roc 11),
      //should be extended to read helicity info from adc 

      Int_t    fQWEAKDebug;          // Debug level
      Bool_t   fHaveROCs;         // Required ROCs are defined
      Bool_t   fNegGate;          // Invert polarity of gate, TO DO implement this functionality
      Int_t    fVerbosity; 

      static const Int_t NHISTR = 12;
      TH1F*    fHistoR[12];  // Histograms

   private:
      static UInt_t FindWord( const THaEvData& evdata, const ROCinfo& info );

      ClassDef(SBSScalerHelicityReader,0) // Helper class for reading QWEAK helicity data
};

#endif

