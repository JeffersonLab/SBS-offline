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
   fMAXBIT(30),
   fHelicityDelay(2),
   fRingFinalEvtNum(1),fRingFinalPatNum(0),fRingFinalSeed(0),
   fSeedbits(0),
   fRingSeed_reported(0),fRingSeed_actual(0),
   fRingPhase_reported(0),fRing_reported_polarity(0),
   fRing_actual_polarity(0), fEvtype(-1), fVerbosity(0), 
   fHelScalerTree(nullptr),fBranch_seed(0),
   fHisto(NHIST, nullptr)
{
}
//_____________________________________________________________________________
   SBSScalerHelicity::SBSScalerHelicity()
:  fMAXBIT(30),
   fHelicityDelay(2),
   fRingFinalEvtNum(1),fRingFinalPatNum(0),fRingFinalSeed(0),
   fSeedbits(0),
   fRingSeed_reported(0),fRingSeed_actual(0),
   fRingPhase_reported(0),fRing_reported_polarity(0),
   fRing_actual_polarity(0), fEvtype(-1), fHisto(NHIST, nullptr)
{
   // Default constructor for ROOT I/O
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
   // maximum of event in the pattern, for now we are working with quartets
   // The first 8 signs are the same for quartets or octets
   // careful, the first value here should always +1
   fPatternSequence = {+1,-1,-1,+1,-1, +1, +1, -1};

   return kOK;
}
//_____________________________________________________________________________
Int_t SBSScalerHelicity::Begin( THaRunBase* )
{
   SBSScalerHelicityReader::Begin();

   fHisto[0] = new TH1F("hel.seed","hel.seed",32,-1.5,30.5);
   fHisto[1] = new TH1F("hel.error.code","hel.error.code",35,-1.5,33.5);


   TString treeName = Form("TShel");
   TString armName  = Form("Lhel");   // for LHRS; eventually make this a user-changeable value
   TString treeInfo = Form("Helicity data plugged into LHRS");
   

   if(!fHelScalerTree){
     // if the tree isn't created yet, create it
      fHelScalerTree = new TTree(treeName,treeInfo);
      fHelScalerTree->SetAutoSave(200000000);

     TString branchName;
     TString branchInfo;

      branchName = Form("%s.physics_trignum", armName.Data());
      branchInfo = Form("%s/i", branchName.Data());
      fHelScalerTree->Branch(branchName, &fTriggerCounter, branchInfo);

      branchName = Form("%s.ring.seedReported", armName.Data());
      branchInfo = Form("%s/i", branchName.Data());
      fHelScalerTree->Branch(branchName, &fRingSeed_reported, branchInfo);
      branchName = Form("%s.ring.seedActual", armName.Data());
      branchInfo = Form("%s/i",branchName.Data());
      fHelScalerTree->Branch(branchName, &fRingSeed_actual, branchInfo);
      branchName = Form("%s.ring.phaseReported", armName.Data());
      branchInfo = Form("%s/i",branchName.Data());
      fHelScalerTree->Branch(branchName, &fRingPhase_reported, branchInfo);
      branchName = Form("%s.ring.polarityReported", armName.Data());
      branchInfo = Form("%s/i",branchName.Data());
      fHelScalerTree->Branch(branchName, &fRing_reported_polarity, branchInfo);
      branchName = Form("%s.ring.polarityActual", armName.Data());
      branchInfo = Form("%s/i",branchName.Data());
      fHelScalerTree->Branch(branchName, &fRing_actual_polarity, branchInfo);
      
      branchName = Form("%s.hel.ErrorCode", armName.Data());
      branchInfo = Form("%s/i",branchName.Data());
      fHelScalerTree->Branch(branchName, &fHelErrorCond, branchInfo);
      branchName = Form("%s.hel.EvtNum", armName.Data());
      branchInfo = Form("%s/i",branchName.Data());
      fHelScalerTree->Branch(branchName, &fNumEvents, branchInfo);
      branchName = Form("%s.hel.PattNum", armName.Data());
      branchInfo = Form("%s/i",branchName.Data());
      fHelScalerTree->Branch(branchName, &fNumPatterns, branchInfo);
      branchName = Form("%s.hel.PattPhase", armName.Data());
      branchInfo = Form("%s/i",branchName.Data());
      fHelScalerTree->Branch(branchName, &fPatternPhase, branchInfo);
      branchName = Form("%s.hel.PatternSeed", armName.Data());
      branchInfo = Form("%s/i",branchName.Data());
      fHelScalerTree->Branch(branchName, &fSeedValue, branchInfo);
      branchName = Form("%s.hel.PatternPolarity", armName.Data());
      branchInfo = Form("%s/i",branchName.Data());
      fHelScalerTree->Branch(branchName, &fPatternHel, branchInfo);
      branchName = Form("%s.hel.EvtPolarity", armName.Data());
      branchInfo = Form("%s/i",branchName.Data());
      fHelScalerTree->Branch(branchName, &fEventPolarity, branchInfo);
      branchName = Form("%s.hel.ReportedQrtHel", armName.Data());
      branchInfo = Form("%s/i",branchName.Data());
      fHelScalerTree->Branch(branchName, &fReportedQrtHel, branchInfo);

      branchName = Form("%s.ring.FinalQrtHel", armName.Data());
      branchInfo = Form("%s/i",branchName.Data());
      fHelScalerTree->Branch(branchName, &fRingFinalQrtHel, branchInfo);
      branchName = Form("%s.ring.FinalEvtNum", armName.Data());
      branchInfo = Form("%s/i",branchName.Data());
      fHelScalerTree->Branch(branchName, &fRingFinalEvtNum, branchInfo);
      branchName = Form("%s.ring.FinalPatNum", armName.Data());
      branchInfo = Form("%s/i",branchName.Data());
      fHelScalerTree->Branch(branchName, &fRingFinalPatNum, branchInfo);
      branchName = Form("%s.ring.FinalSeed", armName.Data());
      branchInfo = Form("%s/i",branchName.Data());
      fHelScalerTree->Branch(branchName, &fRingFinalSeed, branchInfo);

      branchName = Form("%s.fadc.ReportedHelicity", armName.Data());
      branchInfo = Form("%s/i",branchName.Data());
      fHelScalerTree->Branch(branchName, &fFADCHelicity, branchInfo);
      branchName = Form("%s.fadc.PatternSync", armName.Data());
      branchInfo = Form("%s/i",branchName.Data());
      fHelScalerTree->Branch(branchName, &fFADCPatSync, branchInfo);
      branchName = Form("%s.fadc.TSettle", armName.Data());
      branchInfo = Form("%s/i",branchName.Data());
      fHelScalerTree->Branch(branchName, &fFADCTSettle, branchInfo);
      branchName = Form("%s.fadc.ReportedQrtHel", armName.Data());
      branchInfo = Form("%s/i",branchName.Data());
      fHelScalerTree->Branch(branchName, &fFADCQrtHel, branchInfo);

      branchName = Form("%s.ring.PattPhase", armName.Data());
      branchInfo = Form("%s/i",branchName.Data());
      fHelScalerTree->Branch(branchName, &fRingPattPhase, branchInfo);

      for(int i=0;i<32;i++){
	branchName = Form("%s.Yield.Ch%d", armName.Data(),i);
	branchInfo = Form("%s/L",branchName.Data());
	fHelScalerTree->Branch(branchName, &fScalerYield[i], branchInfo);

	branchName = Form("%s.Diff.Ch%d", armName.Data(),i);
	branchInfo = Form("%s/L",branchName.Data());
	fHelScalerTree->Branch(branchName, &fScalerDiff[i], branchInfo);
      }


   }

   return 0;
}
//_____________________________________________________________________________
void SBSScalerHelicity::FillHisto()
{
  //   fHisto[0]->Fill(fRing_NSeed);
  //   fHisto[1]->Fill(fErrorCode);
}

//_____________________________________________________________________________
void SBSScalerHelicity::Clear( Option_t* opt ) {
   // Clear event-by-event data

   THaHelicityDet::Clear(opt);
   SBSScalerHelicityReader::Clear(opt);
   fEvtype = 0;
   fHelicity = kUnknown;

   //   fQrt = 0;
   //   fTSettle = 0;

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

   if(fVerbosity>0) std::cout << "[SBSScalerHelicity::Decode]: Filling histograms... " << std::endl;

   fEvtype = evdata.GetEvType();
   fTriggerCounter = evdata.GetEvNum();

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
	 fSeedbits++;
	 fRingFinalSeed = ((fRingFinalSeed<<1)&0x3ffffffe)|fHelicityRing[i];
	 UInt_t tmpnewbit = fRingFinalSeed & 0x1;
	 UInt_t tmpseed = fRingFinalSeed;
	 if (fSeedbits>fMAXBIT && RanBit30(fRingSeed_reported)==tmpnewbit) {
	   //  The previous reported seed would predict the current heliicty;
	   //  so the predictor is good!
	   
	   //  Run the algorithm to get the pattern sign for delayed reporting
	   for (Int_t idelay=0; idelay<fHelicityDelay; idelay++){
	     tmpnewbit = RanBit30(tmpseed);
	   }
	   if (tmpnewbit == fHelicityRing[i]) patsign = +1;
	   else                               patsign = -1;
	   fRingSeed_actual = tmpseed;
	 } else {
	   patsign = 0;
	   fRingSeed_actual   = 0;
	 }
	 fRingSeed_reported = fRingFinalSeed;

	 if(fVerbosity>0){
	   std::cout << std::hex
		     << "fRingFinalSeed=="<<fRingFinalSeed
		     << "; fRingSeed_reported=="<<fRingSeed_reported
		     << "; fRingSeed_actual=="<<fRingSeed_actual
		     << std::dec << std::endl;
	 }
	 
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
     for (Int_t idelay=0; idelay<fHelicityDelay; idelay++){
       tmpnewbit = RanBit30(tmpseed);
     }
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
