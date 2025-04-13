#ifndef Podd_SBSHelicityDecoderBoard_h_
#define Podd_SBSHelicityDecoderBoard_h_

//////////////////////////////////////////////////////////////////////////
//
// SBSHelicityDecoderBoard
//
// In-time helicity detector - both from ADC and G0 electronics 
// (with delay=0).
// Provides redundancy for cross-checks.
//
//////////////////////////////////////////////////////////////////////////



class SBSHelicityDecoderBoard : public THaHelicity {
  
public:
  SBSHelicityDecoderBoard( const char* name, const char* description, 
	       THaApparatus* a = nullptr );
  virtual ~SBSHelicityDecoderBoard();

  virtual void   Clear( Option_t* opt = "" );
  virtual Int_t  Decode( const THaEvData& evdata );

  SBSHelicityDecoderBoard();  // For ROOT RTTI
  
protected:

  //  Configuration
  Int_t ROC_ID;
  Int_t Bank_ID;
  Int_t fHelicityDelay;   //Helicity delay in # of patterns
    
  // Event-by-event data from the HD module
  UInt_t  fSeed_Reported;
  Int_t   fNum_TStable_Fall;
  Int_t   fNum_TStable_Rise;
  Int_t   fNum_Pattern_Sync;
  Int_t   fNum_Pair_Sync;
  Int_t   fTime_since_TStable;
  Int_t   fTime_since_TSettle;
  Int_t   fLast_Duration_TStable;
  Int_t   fLast_Duration_TSettle;
  //  Values from the status at trigger
  Int_t   fPatternPhase;
  Int_t   fEventPolarity;
  Int_t   fReportedPatternHel;
  Int_t   fBit_Helicity;
  Int_t   fBit_PairSync;
  Int_t   fBit_PatSync;
  Int_t   fBit_TStable;
  //  Event history bit patterns
  UInt_t  fEvtHistory_PatSync;
  UInt_t  fEvtHistory_PairSync;
  UInt_t  fEvtHistory_ReportedHelicity;
  //  Pattern history bit pattern
  UInt_t  fPatHistory_ReportedHelicity;

  //  Calulated values
  UInt_t  fSeed_True;
  UInt_t  fHelicity;

  virtual Int_t DefineVariables( EMode mode = kDefine );
  virtual Int_t ReadDatabase( const TDatime& date );

  ClassDef(SBSHelicityDecoderBoard,2)   // Helicity information from the Helicity Decoder module

};

#endif

