#ifndef Podd_SBSHelicityDecoderBoard_h_
#define Podd_SBSHelicityDecoderBoard_h_

//////////////////////////////////////////////////////////////////////////
//
// SBSHelicityDecoderBoard
//
// Decoding and processing of the helicity informaiton from the
// Helicity Decoder board.
//
//////////////////////////////////////////////////////////////////////////

#include "THaHelicityDet.h"

class SBSHelicityDecoderBoard : public THaHelicityDet {
  
public:
  SBSHelicityDecoderBoard( const char* name, const char* description, 
	       THaApparatus* a = nullptr );
  SBSHelicityDecoderBoard();  // For ROOT RTTI
  virtual ~SBSHelicityDecoderBoard();

  virtual Int_t  Begin( THaRunBase* r=nullptr );
  virtual void   Clear( Option_t* opt = "" );
  virtual Int_t  Decode( const THaEvData& evdata );
  virtual Int_t  End( THaRunBase* r=nullptr );

  void PrintEvent( UInt_t evtnum ){};

  const char* GetDBFileName() const {return "heldecode.";};
  
  
protected:
  virtual void  FillHisto(){};
    
  //  Configuration
  Int_t ROC_ID;
  Int_t Bank_ID;
  size_t fHelicityDelay;   //Helicity delay in # of patterns
    
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
  // fHelicity and fSign belong to the base class
  //    EHelicity fHelicity;  // Beam helicity. fHelicity = fSign * decoded_helicity
  //    Int_t     fSign;      // Overall sign of beam helicity, i.e. IHWP. Default 1.


 

  virtual Int_t DefineVariables( EMode mode = kDefine );
  virtual Int_t ReadDatabase( const TDatime& date );


  Int_t  CheckPredictor(UInt_t& ranseed);
  UInt_t GetRandbit30(UInt_t& ranseed);

  
  //  The following methods 
  class BANKinfo {
  public:
    BANKinfo() : roc(kMaxUInt), bankid(0), index(0), length(0) {};
    Bool_t valid() const;
    UInt_t roc;
    UInt_t bankid;
    UInt_t index;
    UInt_t length;
  };

  const UInt_t*  GetBankBuffer( const THaEvData& evdata, BANKinfo& info );
  Bool_t DecodeBank(const UInt_t *lbuff, const BANKinfo& info);

  static const Int_t  fNumDecoderWords;
  void   FillHDVariables(uint32_t data, uint32_t index);

  ClassDef(SBSHelicityDecoderBoard,0)   // Helicity information from the Helicity Decoder module

};

#endif

