//*-- Author :    Paul King, August 2021
////////////////////////////////////////////////////////////////////////
//
// SBSScalerHelicity
//
// Based on the Podd class HallA/THaQWEAKHelicity
//
// Helicity of the beam from QWEAK electronics in delayed mode
//    +1 = plus,  -1 = minus,  0 = unknown
//
// Also supports in-time mode with delay = 0
// 
////////////////////////////////////////////////////////////////////////

#include "SBSScalerHelicity.h"
#include "THaEvData.h"
#include "TH1F.h"
#include "TMath.h"
#include <iostream>

using namespace std;

//_____________________________________________________________________________
SBSScalerHelicity::SBSScalerHelicity( const char* name, const char* description,
      THaApparatus* app ):
   THaHelicityDet( name, description, app ),
   fOffsetTIRvsRing(3), fQWEAKDelay(8), fMAXBIT(30),
   fQWEAKNPattern(4), HWPIN(true), 
   fHelicityDelay(2),
   fQrt(1), fTSettle(0),fValidHel(false),
   fHelicityLastTIR(0),fPatternLastTIR(0), fErrorCode(0), fRing_NSeed(0),
   fRingFinalEvtNum(1),fRingFinalPatNum(0),fRingFinalSeed(0),
   fRingU3plus(0),fRingU3minus(0),
   fRingT3plus(0),fRingT3minus(0),
   fRingT5plus(0),fRingT5minus(0),
   fRingT10plus(0),fRingT10minus(0),
   fRingTimeplus(0), fRingTimeminus(0),
   fRingSeed_reported(0),fRingSeed_actual(0),
   fRingPhase_reported(0),fRing_reported_polarity(0),
   fRing_actual_polarity(0), fEvtype(-1), fVerbosity(0), 
   fHelScalerTree(nullptr),fBranch_seed(0),fBranch_errCode(0),
   fHisto(NHIST, nullptr)
{
   for (UInt_t j=0; j<32; j++){
      fScalerCumulative[j] = 0;
   }
}
//_____________________________________________________________________________
   SBSScalerHelicity::SBSScalerHelicity()
: fOffsetTIRvsRing(3), fQWEAKDelay(8), fMAXBIT(30),
   fQWEAKNPattern(4), HWPIN(true),
   fHelicityDelay(2),
   fQrt(1), fTSettle(0),fValidHel(false),
   fHelicityLastTIR(0),fPatternLastTIR(0), fErrorCode(0), fRing_NSeed(0),
   fRingFinalEvtNum(1),fRingFinalPatNum(0),fRingFinalSeed(0),
   fRingU3plus(0),fRingU3minus(0),
   fRingT3plus(0),fRingT3minus(0),
   fRingT5plus(0),fRingT5minus(0),
   fRingT10plus(0),fRingT10minus(0),
   fRingTimeplus(0), fRingTimeminus(0),
   fRingSeed_reported(0),fRingSeed_actual(0),
   fRingPhase_reported(0),fRing_reported_polarity(0),
   fRing_actual_polarity(0), fEvtype(-1), fHisto(NHIST, nullptr)
{
   // Default constructor for ROOT I/O
   for (UInt_t j=0; j<32; j++){
      fScalerCumulative[j] = 0;
   }
   fFADCQrtHel = 0;

   for (UInt_t j=0; j<32; j++){
     fScalerYield[j] = 0;
     fScalerDiff[j]  = 0;
   }
   fRingPattPhase = 0;
}
//_____________________________________________________________________________
SBSScalerHelicity::~SBSScalerHelicity()
{
   RemoveVariables();

   // The tree object is owned by ROOT since it gets associated wth the output
   // file, so DO NOT delete it here. 
   if (!TROOT::Initialized()) {
      delete fHelScalerTree;
   }

   // for( Int_t i = 0; i < NHIST; ++i ) {
   //   delete fHisto[i];
   // }
}
//_____________________________________________________________________________
Int_t SBSScalerHelicity::DefineVariables( EMode mode )
{
   // Initialize global variables

   //  cout << "Called SBSScalerHelicity::DefineVariables with mode == "
   //       << mode << endl;

   // Define standard variables from base class
   Int_t ret = THaHelicityDet::DefineVariables( mode );
   if( ret )
      return ret;

   const RVarDef var[] = {
      { "hel", "True helicity for event",              "fHelicity" },
      { "lhrs.fadc.hel", "Helicity bit in LHRS FADC",  "fFADCHelicity"},
      { "lhrs.fadc.pat", "PatternSync in LHRS FADC",   "fFADCPatSync"},
      { "lhrs.fadc.tsettle", "Tsettle in LHRS FADC",   "fFADCTSettle"},
      { "errcode", "Helicity prediction error code",   "fHelErrorCond"},
      { "evtcount", "Number of helicity events",       "fNumEvents"},
      { "patcount", "Number of helicity patterns",     "fNumPatterns"},
      { "patphase", "Event phase within pattern",      "fPatternPhase"},
      { "seed", "Helicity seed value",                 "fSeedValue"},
      { nullptr }
   };
   cout << "now actually defining stuff, prefix = " << fPrefix << endl;
   return DefineVarsFromList( var, mode );
}
//_____________________________________________________________________________
void SBSScalerHelicity::PrintEvent( UInt_t evtnum )
{

   cout<<" ++++++ SBSScalerHelicity::Print ++++++\n";
   cout << dec << "--> Data for spectrometer " << GetName() << endl;
   cout << "  evtype " << fEvtype<<endl;
   cout << " event number ="<<evtnum<<endl;
   cout << " == Input register data =="<<endl;

   cout<<" +++++++++++++++++++++++++++++++++++++\n";
}
//_____________________________________________________________________________
Int_t SBSScalerHelicity::ReadDatabase( const TDatime& date )
{
   // Read general HelicityDet database values (e.g. fSign)
   Int_t st = THaHelicityDet::ReadDatabase( date );
   if( st != kOK )
      return st;

   // Read QWEAK readout parameters (ROC addresses etc.)
   st = SBSScalerHelicityReader::ReadDatabase( GetDBFileName(), GetPrefix(),
	 date, fQWEAKDebug );
   if( st != kOK )
      return st;


   // for now bypass reading the inputs from the database;
   fMAXBIT=30;
   fOffsetTIRvsRing=3;
   fQWEAKDelay=8;
   // maximum of event in the pattern, for now we are working with quartets
   // careful, the first value here should always +1
   fPatternSequence = {1,-1,-1,1};
   HWPIN=true;

   return kOK;
}
//_____________________________________________________________________________
Int_t SBSScalerHelicity::Begin( THaRunBase* )
{
   SBSScalerHelicityReader::Begin();

   fHisto[0] = new TH1F("hel.seed","hel.seed",32,-1.5,30.5);
   fHisto[1] = new TH1F("hel.error.code","hel.error.code",35,-1.5,33.5);

   // D Flay 12/9/21: I think this is where we want to set up 
   // all branches of the tree. It should be of the form "TShel.branchName"
   // or "TSsbsHel.branchName" eventually.  Let's try "hel" for now since 
   // we're only attaching to the LHRS arm at the moment in this test

   TString treeName = Form("TShel");
   TString armName  = Form("Lhel");   // for LHRS; eventually make this a user-changeable value
   TString treeInfo = Form("Helicity data plugged into LHRS");
   
   TString branchInfo;

   int j=0;
   const int NB = 130;
   TString branchName[NB];
   branchName[0]  = Form("%s.error.code"           ,armName.Data());
   branchName[1]  = Form("%s.ring.seed"            ,armName.Data());
   branchName[2]  = Form("%s.ring.seedReported"    ,armName.Data());
   branchName[3]  = Form("%s.ring.seedActual"      ,armName.Data());
   branchName[4]  = Form("%s.ring.phaseReported"   ,armName.Data());
   //   branchName[5]  = Form("%s.ring.polarityReported",armName.Data());
   //   branchName[6]  = Form("%s.ring.polarityActual"  ,armName.Data());
   branchName[5]  = Form("%s.ring.UnewPlus"        ,armName.Data());
   branchName[6] = Form("%s.ring.UnewMinus"        ,armName.Data());
   branchName[7]  = Form("%s.ring.DnewPlus"        ,armName.Data());
   branchName[8]  = Form("%s.ring.DnewMinus"       ,armName.Data());
   branchName[9] = Form("%s.ring.U1Plus"           ,armName.Data());
   branchName[10] = Form("%s.ring.U1Minus"         ,armName.Data());
   branchName[11] = Form("%s.ring.D1Plus"          ,armName.Data());
   branchName[12] = Form("%s.ring.D1Minus"         ,armName.Data());
   branchName[13] = Form("%s.ring.D3Plus"          ,armName.Data());
   branchName[14] = Form("%s.ring.D3Minus"         ,armName.Data());
   branchName[15] = Form("%s.ring.D10Plus"         ,armName.Data());
   branchName[16] = Form("%s.ring.D10Minus"        ,armName.Data());

   for(int i=0;i<32;i++){
      j = 17 + i;
      branchName[j] = Form("%s.cumulative.Ch%d",armName.Data(),i);
   }

   branchName[49] = Form("%s.hel.ErrorCode"        ,armName.Data());
   branchName[50] = Form("%s.hel.EvtNum"           ,armName.Data());
   branchName[51] = Form("%s.hel.PattNum"          ,armName.Data());
   branchName[52] = Form("%s.hel.PattPhase"        ,armName.Data());
   branchName[53] = Form("%s.hel.PatternSeed"      ,armName.Data());
   branchName[54] = Form("%s.hel.PatternPolarity"  ,armName.Data());
   branchName[55] = Form("%s.hel.EvtPolarity"      ,armName.Data());
   branchName[56] = Form("%s.hel.ReportedQrtHel"   ,armName.Data());

   branchName[57] = Form("%s.ring.FinalQrtHel"     ,armName.Data());
   branchName[58] = Form("%s.ring.FinalEvtNum"     ,armName.Data());
   branchName[59] = Form("%s.ring.FinalPatNum"     ,armName.Data());
   branchName[60] = Form("%s.ring.FinalSeed"       ,armName.Data());

   branchName[61] = Form("%s.fadc.ReportedHelicity",armName.Data());
   branchName[62] = Form("%s.fadc.PatternSync"     ,armName.Data());
   branchName[63] = Form("%s.fadc.TSettle"         ,armName.Data());
   branchName[64] = Form("%s.fadc.ReportedQrtHel"  ,armName.Data());

   branchName[65] = Form("%s.ring.PattPhase"       ,armName.Data());

   for(int i=0;i<32;i++){
      j = 66 + 2*i;
      branchName[j]   = Form("%s.Yield.Ch%d",armName.Data(),i);
      branchName[j+1] = Form("%s.Diff.Ch%d",armName.Data(),i);
   }

   if(!fHelScalerTree){
      // if the tree isn't created yet, create it
      fHelScalerTree = new TTree(treeName,treeInfo);
      fHelScalerTree->SetAutoSave(200000000);
      branchInfo = Form("%s/D",branchName[0].Data()); 
      fHelScalerTree->Branch(branchName[0].Data(),&fBranch_errCode,branchInfo.Data()); 
      branchInfo = Form("%s/i",branchName[1].Data()); 
      fHelScalerTree->Branch(branchName[1].Data(),&fRing_NSeed,branchInfo.Data());  
      branchInfo = Form("%s/i",branchName[2].Data()); 
      fHelScalerTree->Branch(branchName[2].Data(),&fRingSeed_reported,branchInfo.Data());  
      branchInfo = Form("%s/i",branchName[3].Data()); 
      fHelScalerTree->Branch(branchName[3].Data(),&fRingSeed_actual,branchInfo.Data());  
      branchInfo = Form("%s/i",branchName[4].Data()); 
      fHelScalerTree->Branch(branchName[4].Data(),&fRingPhase_reported,branchInfo.Data());  
      branchInfo = Form("%s/i",branchName[5].Data()); 
      fHelScalerTree->Branch(branchName[5].Data(),&fRing_reported_polarity,branchInfo.Data());  
      branchInfo = Form("%s/i",branchName[6].Data()); 
      fHelScalerTree->Branch(branchName[6].Data(),&fRing_actual_polarity,branchInfo.Data());  
      branchInfo = Form("%s/i",branchName[7].Data()); 
      fHelScalerTree->Branch(branchName[7].Data(),&fRingU3plus,branchInfo.Data());  
      branchInfo = Form("%s/i",branchName[8].Data()); 
      fHelScalerTree->Branch(branchName[8].Data(),&fRingU3minus,branchInfo.Data());  
      branchInfo = Form("%s/i",branchName[9].Data()); 
      fHelScalerTree->Branch(branchName[9].Data(),&fRingT3plus,branchInfo.Data());  
      branchInfo = Form("%s/i",branchName[10].Data()); 
      fHelScalerTree->Branch(branchName[10].Data(),&fRingT3minus,branchInfo.Data());  
      branchInfo = Form("%s/i",branchName[11].Data()); 
      fHelScalerTree->Branch(branchName[11].Data(),&fRingT5plus,branchInfo.Data());  
      branchInfo = Form("%s/i",branchName[12].Data()); 
      fHelScalerTree->Branch(branchName[12].Data(),&fRingT5minus,branchInfo.Data()); 
      branchInfo = Form("%s/i",branchName[13].Data()); 
      fHelScalerTree->Branch(branchName[13].Data(),&fRingT10plus,branchInfo.Data());  
      branchInfo = Form("%s/i",branchName[14].Data()); 
      fHelScalerTree->Branch(branchName[14].Data(),&fRingT10minus,branchInfo.Data());  
      branchInfo = Form("%s/i",branchName[15].Data()); 
      fHelScalerTree->Branch(branchName[15].Data(),&fRingTimeplus,branchInfo.Data());  
      branchInfo = Form("%s/i",branchName[16].Data()); 
      fHelScalerTree->Branch(branchName[16].Data(),&fRingTimeminus,branchInfo.Data());  
      for(int i=0;i<32;i++){
	 j = 17 + i; 
	 branchInfo = Form("%s/L",branchName[j].Data()); 
	 fHelScalerTree->Branch(branchName[j].Data(),&fScalerCumulative[i],branchInfo.Data()); 
      }
      branchInfo = Form("%s/i",branchName[49].Data());
      fHelScalerTree->Branch(branchName[49].Data(),&fHelErrorCond,branchInfo.Data());
      branchInfo = Form("%s/i",branchName[50].Data());
      fHelScalerTree->Branch(branchName[50].Data(),&fNumEvents,branchInfo.Data());
      branchInfo = Form("%s/i",branchName[51].Data());
      fHelScalerTree->Branch(branchName[51].Data(),&fNumPatterns,branchInfo.Data());
      branchInfo = Form("%s/i",branchName[52].Data());
      fHelScalerTree->Branch(branchName[52].Data(),&fPatternPhase,branchInfo.Data());
      branchInfo = Form("%s/i",branchName[53].Data());
      fHelScalerTree->Branch(branchName[53].Data(),&fSeedValue,branchInfo.Data());
      branchInfo = Form("%s/i",branchName[54].Data());
      fHelScalerTree->Branch(branchName[54].Data(),&fPatternHel,branchInfo.Data());
      branchInfo = Form("%s/i",branchName[55].Data());
      fHelScalerTree->Branch(branchName[55].Data(),&fEventPolarity,branchInfo.Data());
      branchInfo = Form("%s/i",branchName[56].Data());
      fHelScalerTree->Branch(branchName[56].Data(),&fReportedQrtHel,branchInfo.Data());

      branchInfo = Form("%s/i",branchName[57].Data());
      fHelScalerTree->Branch(branchName[57].Data(),&fRingFinalQrtHel,branchInfo.Data());
      branchInfo = Form("%s/i",branchName[58].Data());
      fHelScalerTree->Branch(branchName[58].Data(),&fRingFinalEvtNum,branchInfo.Data());
      branchInfo = Form("%s/i",branchName[59].Data());
      fHelScalerTree->Branch(branchName[59].Data(),&fRingFinalPatNum,branchInfo.Data());
      branchInfo = Form("%s/i",branchName[60].Data());
      fHelScalerTree->Branch(branchName[60].Data(),&fRingFinalSeed,branchInfo.Data());

      branchInfo = Form("%s/i",branchName[61].Data());
      fHelScalerTree->Branch(branchName[61].Data(),&fFADCHelicity,branchInfo.Data());
      branchInfo = Form("%s/i",branchName[62].Data());
      fHelScalerTree->Branch(branchName[62].Data(),&fFADCPatSync,branchInfo.Data());
      branchInfo = Form("%s/i",branchName[63].Data());
      fHelScalerTree->Branch(branchName[63].Data(),&fFADCTSettle,branchInfo.Data());
      branchInfo = Form("%s/i",branchName[64].Data());
      fHelScalerTree->Branch(branchName[64].Data(),&fFADCQrtHel,branchInfo.Data());

      branchInfo = Form("%s/i",branchName[65].Data());
      fHelScalerTree->Branch(branchName[65].Data(),&fRingPattPhase,branchInfo.Data());

      for(int i=0;i<32;i++){
	j = 66 + 2*i;
	branchInfo = Form("%s/L",branchName[j].Data()); 
	fHelScalerTree->Branch(branchName[j].Data(),&fScalerYield[i],branchInfo.Data());
	j = 66 + 2*i + 1;
	branchInfo = Form("%s/L",branchName[j].Data()); 
	fHelScalerTree->Branch(branchName[j].Data(),&fScalerDiff[i],branchInfo.Data());
      }


   }

   return 0;
}
//_____________________________________________________________________________
void SBSScalerHelicity::FillHisto()
{
   fHisto[0]->Fill(fRing_NSeed);
   fHisto[1]->Fill(fErrorCode);
}
//_____________________________________________________________________________
void SBSScalerHelicity::SetErrorCode(Int_t error)
{
   // used as a control for the helciity computation
   // 2^0: if the reported number of events in a  pattern is larger than fQWEAKNPattern
   // 2^1: if the offset between the ring reported value and TIR value is not fOffsetTIRvsRing
   // 2^2: if the reported time in the ring is 0
   // 2^3: if the predicted reported helicity doesn't match the reported helicity in the ring
   // 2^4: if the helicity cannot be computed using the SetHelicity routine
   // 2^5: if seed is being gathered

   if(fErrorCode==0)
      fErrorCode=(1<<error);
   // only one reported error at the time
}
//_____________________________________________________________________________
void SBSScalerHelicity::Clear( Option_t* opt ) {
   // Clear event-by-event data

   THaHelicityDet::Clear(opt);
   SBSScalerHelicityReader::Clear(opt);
   fEvtype = 0;
   fHelicity = kUnknown;

   fQrt = 0;
   fTSettle = 0;
   fRingU3plus = 0;
   fRingU3minus = 0;
   fRingT3plus = 0;
   fRingT3minus = 0;
   fRingT5plus = 0;
   fRingT5minus = 0;
   fRingT10plus = 0;
   fRingT10minus = 0;
   fRingTimeplus = 0;
   fRingTimeminus = 0;
   fErrorCode = 0;
   
   //   for (UInt_t j=0; j<32; j++){
   //     fScalerYield[j] = 0;
   //     fScalerDiff[j]  = 0;
   //   }
   //   fRingPattPhase = 0;

}
//_____________________________________________________________________________
Int_t SBSScalerHelicity::Decode( const THaEvData& evdata )
{
  static  Long_t helsign=0, patsign=0;
   // Decode Helicity data.
   // Return 1 if helicity was assigned, 0 if not, <0 if error.

   /*
    *  std::cout << "\n\nCumulative counts:   chan0, chan9, chan15:  " << std::dec
    *        << fScalerCumulative[0] << " " << fScalerCumulative[9] << " "
    *        << fScalerCumulative[15]
    *        << std::endl;
    */

   Int_t err = ReadData( evdata ); // from SBSScalerHelicityReader class
   if( err ) {
      Error( Here("SBSScalerHelicity::Decode"), "Error decoding helicity data." );
      return err;
   }
   //Here should accumulate the fScalerCumulative values, if we've read new entries.
   if (fIRing>0){
      for (UInt_t i=0; i<fIRing; i++){
	 for (UInt_t j=0; j<32; j++){
	    fScalerCumulative[j] += (Long64_t)fScalerRing[i][j];
	 }
      }
      if(fVerbosity>1){
	 std::cout << "[SBSScalerHelicity::Decode]: Cumulative counts:   chan0, chan9, chan15:  " << std::dec
	    << fScalerCumulative[0] << " " << fScalerCumulative[9] << " "
	    << fScalerCumulative[15]
	    << std::endl;
      }
   }

   if(fVerbosity>0) std::cout << "[SBSScalerHelicity::Decode]: Filling histograms... " << std::endl; 

   fEvtype = evdata.GetEvType();
   int trig_num = evdata.GetEvNum();
   fRingU3plus = trig_num;   
   /*
    *    *  TODO:  the follow funcitons were for the "old" data structure, and should
    *       *         be reconsidered if they should be adapted or removed.
    *         LoadHelicity(evdata.GetEvNum()); 
    *           if(fQWEAKDebug>1)
    *               PrintEvent(evdata.GetEvNum());
    *                 CheckTIRvsRing(evdata.GetEvNum());
    *                   if(fErrorCode==0)  
    *                       fValidHel=true;
    *                         else
    *                             fValidHel=false;
    *                               */

   // assign variables that will get to the tree 
   fBranch_seed         = fRing_NSeed; 
   fBranch_errCode      = fErrorCode;

   if (fIRing>0){
     for (UInt_t i=0; i<fIRing; i++){
       fRingFinalQrtHel = fPatternRing[i] + fHelicityRing[i];
       //  Get the sign from the reported helicity
       if (fHelicityRing[i]==0) helsign = -1;
       else                     helsign = +1;
       //  Increment event number and pattern number/phase counters.  Maintain the reported seed value.
       fRingFinalEvtNum++;
       if (fPatternRing[i]==0x10){
	 fRingPattPhase = 0;
	 fRingHelicitySum = 0;
	 fRingFinalPatNum++;
	 fRingFinalSeed = ((fRingFinalSeed<<1)&0x3ffffffe)|fHelicityRing[i];
	 //  Run the algorithm to get the pattern sign for delayed reporting
	 UInt_t tmpnewbit = fRingFinalSeed & 0x1;
	 UInt_t tmpseed = fRingFinalSeed;
	 for (size_t idelay=0; idelay<fHelicityDelay; idelay++){
	   tmpnewbit = RanBit30(tmpseed);
	 }
	 if (tmpnewbit == fHelicityRing[i]) patsign = +1;
	 else                               patsign = -1;
       } else {
	 fRingPattPhase++;
       }
       helsign *= patsign;
       fRingHelicitySum += helsign;
       //  Bulid the helicity-independent and -dependent sums. 
       if (fRingPattPhase==0){
	 fTimeStampYield = 0;
	 fTimeStampDiff  = 0;
       }
       fTimeStampYield += fTimeStampRing[i];
       fTimeStampDiff  += helsign * fTimeStampRing[i];
       for (UInt_t j=0; j<32; j++){
	 if (fRingPattPhase==0){
	   fScalerYield[j] = 0;
	   fScalerDiff[j]  = 0;
	 }
	 fScalerYield[j] += +1 * fScalerRing[i][j];
	 fScalerDiff[j]  += helsign * fScalerRing[i][j];
       }
       //  Fill histograms and tree values for each scaler read
       FillHisto();
       if(fHelScalerTree) fHelScalerTree->Fill();
     }
   }

   // UInt_t tmpseed = fRingFinalSeed;
   // if (tmpseed != fSeedValue) {
   //   std::cout << "fRingFinalSeed != fSeedValue: " << std::hex
   // 	       << tmpseed << " " << fSeedValue <<std::dec <<std::endl;
   // } else {
   //   std::cout << "ok" << std::endl;
   // }

   //  Calculate the true helicity
   if (fHelErrorCond==0){
     UInt_t tmpseed = fSeedValue;
     UInt_t tmpnewbit = fSeedValue & 0x1;
     for (size_t idelay=0; idelay<fHelicityDelay; idelay++){
       tmpnewbit = RanBit30(tmpseed);
     }
     // if (tmpseed != fSeedValue) {
     //   std::cout << "tmpseed != fSeedValue: " << std::hex
     // 		 << tmpseed << " " << fSeedValue <<std::dec <<std::endl;
     // }
     if (tmpnewbit==0) fHelicity = kMinus;
     else              fHelicity = kPlus;
     if (fPatternPhase==1 || fPatternPhase==2){
       if (tmpnewbit==0) fHelicity = kPlus;
       else              fHelicity = kMinus;
     }
   } else {
     fHelicity = kUnknown;
   }

   
   if(fVerbosity>0) std::cout << "[SBSScalerHelicity::Decode]: --> Done. " << std::endl; 

   return 0;
}
//_____________________________________________________________________________
void SBSScalerHelicity::CheckTIRvsRing( UInt_t eventnumber )
{
   // here one checks that the offset between the TIR helicity reports 
   // and the Ring report is as expected (fQWEAKOffset)
   // This is a simplified comparison: ie not every  TIR readings is 
   // compared to the Ring readings for example, and offset=3
   // for simplicity the comparison only if the  ring buffer contains at least 
   // fQWEAKOffset readings

   static const char* const here = "SBSScalerHelicity::CheckTIRvsRing";

   if(fIRing>=fOffsetTIRvsRing)
   {
      // compare the two values (last TIR) and reading in the current ring
      if(fHelicityRing[fOffsetTIRvsRing-1]!=fHelicityLastTIR
	    || fPatternRing[fOffsetTIRvsRing-1]!=fPatternLastTIR)
      {
	 fHelicity=kUnknown;
	 fRing_NSeed=0;
	 if(fQWEAKDebug>0)
	 {
	    cout<<here<<" BAD !! the offset between the helicity ring and the input register ";
	    cout << "is not what is expected: reset the seed !! event#" << eventnumber << "\n";
	 }
	 SetErrorCode(1);
	 if(fQWEAKDebug>1)
	 {
	    cout<<"=====================================\n";
	    cout<<here<<endl;
	    cout << " Event number =" << eventnumber << endl;
	    cout<<" fOffsetTIRvsRing ="<<fOffsetTIRvsRing;
	    cout<<" fHelicityLastTIR ="<<fHelicityLastTIR;
	    cout<<" fPatternLastTIR ="<<fPatternLastTIR;
	    cout<<" fIRing="<<fIRing<<endl;
	    cout<<"RING data: helicity="<<fHelicityRing[fOffsetTIRvsRing-1]
	       <<" pattern="<<fPatternRing[fOffsetTIRvsRing-1]<<endl;
	 }     
      }   
   }   
   fHelicityLastTIR= fHelicityTir;
   fPatternLastTIR=fPatternTir;
}
//_____________________________________________________________________________
Int_t SBSScalerHelicity::End( THaRunBase* )
{
   // End of run processing. Write histograms.
   SBSScalerHelicityReader::End();

   for( Int_t i = 0; i < NHIST; ++i )
      fHisto[i]->Write();

   // D Flay 12/9/21: We should be calling fHelScalerTree->Write()
   if(fHelScalerTree) fHelScalerTree->Write();  

   return 0;
}
//_____________________________________________________________________________
void SBSScalerHelicity::SetDebug( Int_t level )
{
   // Set debug level of this detector as well as the SBSScalerHelicityReader 
   // helper class.

   THaHelicityDet::SetDebug( level );
   fQWEAKDebug = level;
}
//_____________________________________________________________________________
void SBSScalerHelicity::LoadHelicity( UInt_t eventnumber )
{
   static const char* const here = "SBSScalerHelicity::LoadHelicity";

   for (UInt_t i=0;i<fIRing;i++)
   {
      //Check for the number of events between two Pattern signals
      if(fPatternRing[i]==1)
      {
	 // could  add a check that the number of pattern is what is expected. EG 4 for a quartet
	 fRingPhase_reported=0;
      }
      else
      {
	 fRingPhase_reported+=1;
      }
      if(fRingPhase_reported>fQWEAKNPattern)
      {
	 if(fQWEAKDebug>0)
	    cout << here << " Reset seed: The pattern has too many elements !! "
	       << "Should only go up to " << fQWEAKNPattern
	       << " but is now " << fRingPhase_reported
	       << " at event #=" << eventnumber
	       << endl;
	 fRing_NSeed=0;
	 SetErrorCode(0);
      }

      //Check that events in the ring is not empty, using the timestamp
      if(fTimeStampRing[i]==0)
      {
	 fRing_NSeed=0;
	 SetErrorCode(2);
      }

      // if seed has been gathered: 
      // check event by event that the seed make sense:
      // fRing_polarity!=fHelcityRing[i]
      // this should come before the section for the seed gathering::
      // do not change this order
      if(fRing_NSeed==fMAXBIT && fPatternRing[i]==1)
      {
	 fRing_reported_polarity=RanBit30(fRingSeed_reported);
	 fRing_actual_polarity=RanBit30(fRingSeed_actual);
	 if(fRing_reported_polarity!=fHelicityRing[i])
	 {
	    if(fQWEAKDebug>0)
	       cout<<here<<" Catastrophe  !!"
		  <<" predicted helicity doesn't match reported helicity !!!"
		  <<" event #="<<eventnumber
		  <<endl;
	    if(fQWEAKDebug>1)
	       cout<<" iring="<<Form("%02d", i)
		  <<" predicted helicity="<<fRing_reported_polarity
		  <<" fHelicityRing["<<i<<"]="<<fHelicityRing[i]<<endl;
	    fRing_NSeed=0;
	    SetErrorCode(3);
	 }
      }
      // Here is the seed gathering if necessary
      if(fRing_NSeed<fMAXBIT && fPatternRing[i]==1)
      {
	 SetErrorCode(5);
	 fRingSeed_reported=
	    ((fRingSeed_reported<<1)&0x3FFFFFFF)|fHelicityRing[i];
	 fRing_NSeed+=1;
	 if (fRing_NSeed==fMAXBIT)
	 {
	    fRingSeed_actual=fRingSeed_reported;
	    UInt_t advance=0;
	    //take the delay into account
	    for(UInt_t j=0;j<fQWEAKDelay;j++)
	    {
	       advance+=1;
	       if( advance == fQWEAKNPattern)
	       {
		  fRing_actual_polarity=RanBit30(fRingSeed_actual);
		  advance=0;
	       }
	    }
	 }
      }
      //now compute helicity related quantities
      EHelicity localhelicity = kUnknown;
      if(fRing_NSeed==fMAXBIT)
      {
	 localhelicity=SetHelicity(fRing_actual_polarity,fRingPhase_reported);
	 if(localhelicity==kPlus)
	 {
	    fRingU3plus+=fU3Ring[i];
	    fRingT3plus+=fT3Ring[i];
	    fRingT5plus+=fT5Ring[i];
	    fRingT10plus+=fT10Ring[i];
	    fRingTimeplus+=fTimeStampRing[i];
	 }
	 else if(localhelicity==kMinus)
	 {
	    fRingU3minus+=fU3Ring[i];
	    fRingT3minus+=fT3Ring[i];
	    fRingT5minus+=fT5Ring[i];
	    fRingT10minus+=fT10Ring[i];
	    fRingTimeminus+=fTimeStampRing[i];
	 }
	 else
	 {
	    if(fQWEAKDebug>0)
	       cout<<here<<" TROUBLE !! Local helicity doesn't make sense"
		  <<" event #="<<eventnumber<<endl;
	    fRing_NSeed=0;
	    SetErrorCode(4);
	 }

      }
      if(fQWEAKDebug>1)
      {
	 cout<<" iring="<<Form("%02d", i)
	    <<" hel,pat="
	    <<fHelicityRing[i]<<" ,"
	    <<fPatternRing[i]
	    <<" phase="<<fRingPhase_reported
	    <<" NSeed="<<fRing_NSeed
	    <<" Ring(pol reported="<<fRing_reported_polarity
	    <<", actual="<<fRing_actual_polarity
	    <<", actual hel="<<localhelicity
	    <<")"
	    <<endl;

      }
   }
   // now go and get the true helicity for the event in the TIR 
   if(fRing_NSeed==fMAXBIT)
   {
      if(fTSettleTir==1)
      {
	 fHelicity=kUnknown;
	 fTSettle=1;
      }
      else
      {
	 fTSettle=0;
	 UInt_t localfPhase=fRingPhase_reported;
	 UInt_t localfSeed=fRingSeed_actual;
	 UInt_t localfPolarity=fRing_actual_polarity;

	 for(UInt_t i=0; i<fOffsetTIRvsRing;i++)
	 {
	    localfPhase+=1;
	    if( localfPhase == fQWEAKNPattern)
	    {
	       localfPhase=0;
	       localfPolarity=RanBit30(localfSeed);
	    }
	 }
	 fHelicity=SetHelicity(localfPolarity,localfPhase);
	 if(fPatternTir==1)
	    fQrt=1;
	 else
	    fQrt=0;
	 if(fHelicity==kUnknown)
	 {
	    if(fQWEAKDebug>0)
	       std::cout<<"TROUBLE !!! when predicting the actual helicity for the CODA "
		  <<" event #="<<eventnumber<<endl;
	    fRing_NSeed=0;
	 }
      }
   }
   else
   {
      fHelicity=kUnknown;
   }
}
//_____________________________________________________________________________
   THaHelicityDet::EHelicity
SBSScalerHelicity::SetHelicity( UInt_t polarity, UInt_t phase )
{
   // here predicted_reported_helicity can have a value of 0 or 1
   // fPatternSequence[fRingPhase_reported] can have a value of 1 or -1

   THaHelicityDet::EHelicity localhel = kUnknown;

   Int_t select = static_cast<Int_t>(polarity) + fPatternSequence[phase];
   if( select == -1 || select == 2 ) {
      if( HWPIN )
	 localhel = kPlus;
      else
	 localhel = kMinus;
   } else if( select == 1 || select == 0 ) {
      if( HWPIN )
	 localhel = kMinus;
      else
	 localhel = kPlus;
   }

   // std::cout<<" ++++ SBSScalerHelicity::SetHelicity \n";
   // std::cout<<" actual ring  polarity="<<polarity
   //       <<" fRingPhase_reported="<<fPatternSequence[phase]
   //       <<" HWP="<<HWPIN
   //       <<" select="<<select
   //       <<" Helicity="<<localhel<<endl;

   return localhel;
}
//_____________________________________________________________________________
UInt_t SBSScalerHelicity::RanBit30( UInt_t& ranseed )
{

   bool bit7    = (ranseed & 0x00000040) != 0;
   bool bit28   = (ranseed & 0x08000000) != 0;
   bool bit29   = (ranseed & 0x10000000) != 0;
   bool bit30   = (ranseed & 0x20000000) != 0;

   UInt_t newbit = (bit30 ^ bit29 ^ bit28 ^ bit7) & 0x1;

   if(ranseed<=0) {
      if(fQWEAKDebug>1)
	 std::cerr<<"ranseed must be greater than zero!"<<"\n";

      newbit = 0;
   }
   ranseed =  ( (ranseed<<1) | newbit ) & 0x3FFFFFFF;
   //here ranseed is changed
   if( fQWEAKDebug > 1 ) {
      cout << "SBSScalerHelicity::RanBit30, newbit=" << newbit << "\n";
   }

   // Returns 0 or 1
   return newbit;
}
//_____________________________________________________________________________
ClassImp(SBSScalerHelicity) 
