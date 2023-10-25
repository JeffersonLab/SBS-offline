//*-- Author :    Paul King    August 2021
////////////////////////////////////////////////////////////////////////
//
// SBSScalerHelicityReader
//
// Based on the Podd class HallA/THaQWEAKHelicityReader
//
////////////////////////////////////////////////////////////////////////

#include "SBSScalerHelicityReader.h"
#include "THaEvData.h"
#include "TError.h"
#include "THaAnalysisObject.h"   // For LoadDB
#include <iostream>
#include "TH1F.h"

using namespace std;

//____________________________________________________________________
Bool_t SBSScalerHelicityReader::ROCinfo::valid() const
{
   return roc < Decoder::MAXROC;
}
//____________________________________________________________________
Bool_t SBSScalerHelicityReader::BANKinfo::valid() const
{
   return roc < Decoder::MAXROC;
}
//____________________________________________________________________
SBSScalerHelicityReader::SBSScalerHelicityReader()
   : fPatternTir(0), fHelicityTir(0), fTSettleTir(0), fTimeStampTir(0),
   fOldTimeStampTir(0), fIRing(0),
   fHelicityRing{}, fPatternRing{},
   fTimeStampRing{}, fT3Ring{}, fU3Ring{}, fT5Ring{}, fT10Ring{},
   fQWEAKDebug(0),      // Debug level
   fHaveROCs(false),   // Required ROCs are defined
   fNegGate(false),    // Invert polarity of gate, so that 0=active
   fFADCLow_min(30000),  fFADCLow_max(60000),
   fFADCHigh_min(75000), fFADCHigh_max(95000),
   fVerbosity(0),
   fHistoR{}
{
   // Default constructor
}
//____________________________________________________________________
SBSScalerHelicityReader::~SBSScalerHelicityReader()
{
   // Destructor

   // Histograms will be deleted by ROOT
   // for( Int_t i = 0; i < NHISTR; ++i ) {
   //   delete fHistoR[i];
   // }
}
//____________________________________________________________________
void SBSScalerHelicityReader::Print()
{
   cout<<"================================================\n";
   cout<<" SBSScalerHelicityReader::Print() \n";
   cout<<endl;
   cout<<"fPatternTir, fHelicityTir, fTSettleTir="<< fPatternTir
      <<" , "<<fHelicityTir<<" , "<<fTSettleTir<<endl;
   cout<<"fTimeStampTir ="<<fTimeStampTir<<endl;
   cout<<"fIRing="<<fIRing<<endl<<endl;
   for (UInt_t j=0;j<fIRing;j++)
   {
      cout<<j<<"Pattern, helicity, time, T3, U3, T5, T10="
	 <<fHelicityRing[j]<<" , "
	 <<fPatternRing[j]<<" , "
	 <<fTimeStampRing[j]<<" , "
	 <<fT3Ring[j]<<" , "
	 <<fU3Ring[j]<<" , "
	 <<fT5Ring[j]<<" , "
	 <<fT10Ring[j]<<endl;
   }
   if(fIRing==0)
   {
      cout<<0<<"Pattern, helicity, time, T3, U3, T5, T10="
	 <<fHelicityRing[0]<<" , "
	 <<fPatternRing[0]<<" , "
	 <<fTimeStampRing[0]<<" , "
	 <<fT3Ring[0]<<" , "
	 <<fU3Ring[0]<<" , "
	 <<fT5Ring[0]<<" , "
	 <<fT10Ring[0]<<endl;
   }
   cout<<"================================================\n";
}
//____________________________________________________________________
void SBSScalerHelicityReader::Clear( Option_t* )
{
   fPatternTir=3;
   fHelicityTir=3;
   fTSettleTir=3;
   fTimeStampTir=0;
   fIRing=0;
   for(int i=0;i<kHelRingDepth;i++)
   {
      fHelicityRing[i]=3;
      fPatternRing[i]=3;
      fTimeStampRing[i]=0;
      fT3Ring[i]=0;
      fU3Ring[i]=0;
      fT5Ring[i]=0;
      fT10Ring[i]=0;
   }
}
//____________________________________________________________________
UInt_t SBSScalerHelicityReader::FindWord( const THaEvData& evdata,
      const ROCinfo& info )
   // find the index of the word we are looking for given a header
   // or simply return the index already stored in ROC info (in that case
   // the header stored in ROC info needs to be 0 
{
   UInt_t len = evdata.GetRocLength(info.roc);
   if (len <= 4)
      return -1;

   UInt_t i = 0;
   if( info.header == 0 )
      i = info.index;
   else {
      for( ; i<len &&
	    (evdata.GetRawData(info.roc, i) & 0xffff000) != info.header;
	    ++i) {}
      i += info.index;
   }
   return (i < len) ? i : kMaxUInt;
}
//____________________________________________________________________
const UInt_t* SBSScalerHelicityReader::GetBankBuffer( const THaEvData& evdata,
      BANKinfo& info )
{
   UInt_t loff; // Local offset
   UInt_t banklen, bankhead;
   UInt_t len = evdata.GetRocLength(info.roc);
   if (len <= 4) {
      info.length = 0;
      return nullptr;
   }
   const UInt_t* buff = evdata.GetRawDataBuffer(info.roc);
   UInt_t rochead = buff[1];
   if ((rochead&0xff00)!=0x1000){
      std::cout <<  "ROC is not made of banks"  <<std::endl;
      info.length = 0;
      return nullptr;
   }
   loff = 2;
   while (loff < len-2){
      banklen  = buff[loff];
      bankhead = buff[loff+1];
      if (((bankhead&0xffff0000)>>16)==info.bankid) break;
      loff += 1+banklen;
   }
   if ( (((bankhead&0xffff0000)>>16)==info.bankid) &&
	 ((loff+1+banklen) <= len) ){
      buff+=loff;
      info.length = banklen;
   } else {
      info.length = 0;
      return nullptr;
   }
   return buff;
}
//_____________________________________________________________________________
Int_t SBSScalerHelicityReader::ReadDatabase( const char* /*dbfilename*/,
      const char* /*prefix*/,
      const TDatime& /*date*/,
      int /*debug_flag*/ )
{
   // Read parameters from database:  ROC info (detector map), QWEAK delay value
   // TODO: for now I will bypass this call and just hardcode the data I need
   //  static const char* const here = "SBSScalerHelicityReader::ReadDatabase";

   // for now bypass the reading from the data base

   SetROCinfo(kHel,11,0,3);
   SetROCinfo(kTime,11,0,4);
   SetROCinfo(kRing,11,0,0);

   return THaAnalysisObject::kOK;
}

//_____________________________________________________________________________
void SBSScalerHelicityReader::Begin()
{
   // static const char* const here = "SBSScalerHelicityReader::Begin";
   // cout<<here<<endl;

   fHistoR[0]=new TH1F("hel.Pattern.TIR","hel.Pattern.TIR",5,-0.75, 1.75);
   fHistoR[1]=new TH1F("hel.TSettle.TIR","hel.TSettle.TIR",5,-0.75, 1.75);
   fHistoR[2]=new TH1F("hel.Reported.Helicity.TIR","hel.Reported.Helicity.TIR"
	 ,5,-0.75, 1.75);
   fHistoR[3]=new TH1F("hel.dTimestamp.TIR","hel.dTimestamp.TIR",1000,-50,49950);
   fHistoR[4]=new TH1F("hel.Pattern.Ring","hel.Pattern.Ring",5,-0.75, 1.75);
   fHistoR[5]=new TH1F("hel.Reported.Helicity.Ring","hel.Reported.Helicity.Ring"
	 ,5,-0.75, 1.75);
   fHistoR[6]=new TH1F("hel.Timestamp.Ring","hel.Timestamp.Ring",100,-10,490);
   fHistoR[7]=new TH1F("hel.T3.Ring","hel.T3.Ring",53,-1.5, 50.5);
   fHistoR[8]=new TH1F("hel.U3.Ring","hel.U3.Ring",100,-0.5, 99.5);
   fHistoR[9]=new TH1F("hel.T5.Ring","hel.T5.Ring",53,-1.5, 50.5);
   fHistoR[10]=new TH1F("hel.T10.Ring","hel.T10.Ring",53,-1.5, 50.5);
   fHistoR[11]=new TH1F("hel.NRing","hel.NRing",503,-1.5, 501.5);
}
//____________________________________________________________________
void SBSScalerHelicityReader::End()
{
   // static const char* const here = "SBSScalerHelicityReader::End";
   // cout<<here<<endl;

   for(auto & i : fHistoR)
      i->Write();
}
//____________________________________________________________________
Int_t SBSScalerHelicityReader::ReadData( const THaEvData& evdata )
{
   // Obtain the present data from the event for QWEAK helicity mode.
   const UInt_t *lbuff;
   static UInt_t evt_in_pattern;
   static Int_t evtnum_offset=0;
   static Int_t patnum_offset=0;

   static const char* here = "SBSScalerHelicityReader::ReadData";

   fIRing=0;

   //  std::cout<<" kHel, kTime, kRing="<< kHel<<" "<<kTime<<" "<<kRing<<endl;
   //     for (int jk=0; jk<3; jk++)
   //     {
   //       std::cout<<" which="<<jk
   //           <<" roc="<<fROCinfo[jk].roc
   //           <<" header="<<fROCinfo[jk].header
   //           <<" index="<<fROCinfo[jk].index 
   //           <<endl;
   //     }

   //  std::cout<<" fHaveROCs="<<fHaveROCs<<endl;
   if( !fHaveROCs ) {
      ::Error( here, "ROC data (detector map) not properly set up." );
      return -1;
   }

   //  LHRS FADC helicity bits
   fFADCSpare    = evdata.GetData( Decoder::kPulseIntegral,10,14,12,0 );
   fFADCHelicity = evdata.GetData( Decoder::kPulseIntegral,10,14,15,0 );
   fFADCPatSync  = evdata.GetData( Decoder::kPulseIntegral,10,14,14,0 );
   fFADCTSettle  = evdata.GetData( Decoder::kPulseIntegral,10,14,13,0 );

   UInt_t hel, patsync, tsettle;
   UInt_t hel_sbs, patsync_sbs, tsettle_sbs;
   if (fFADCHelicity>=fFADCLow_min && fFADCHelicity<=fFADCLow_max){
     hel = 0;
   }
   else if (fFADCHelicity>=fFADCHigh_min && fFADCHelicity<=fFADCHigh_max){
     hel = 1;
   }
   else {
     hel = 0x1000;
   }
   if (fFADCPatSync>=fFADCLow_min && fFADCPatSync<=fFADCLow_max){
     patsync = 0;
   }
   else if (fFADCPatSync>=fFADCHigh_min && fFADCPatSync<=fFADCHigh_max){
     patsync = 0x10;
   }
   else {
     patsync = 0x2000;
   }
   if (fFADCTSettle>=fFADCLow_min && fFADCTSettle<=fFADCLow_max){
     tsettle = 0;
   }
   else if (fFADCTSettle>=fFADCHigh_min && fFADCTSettle<=fFADCHigh_max){
     tsettle = 0x100;
   }
   else {
     tsettle = 0x4000;
   }
   fFADCQrtHel = hel + patsync + tsettle;


   //  SBS FADC helicity bits
   fFADCSpare_sbs    = evdata.GetData( Decoder::kPulseIntegral,1,20,15,0 );
   fFADCHelicity_sbs = evdata.GetData( Decoder::kPulseIntegral,1,20,14,0 );
   fFADCPatSync_sbs  = evdata.GetData( Decoder::kPulseIntegral,1,20,13,0 );
   fFADCTSettle_sbs  = evdata.GetData( Decoder::kPulseIntegral,1,20,12,0 );
   //  Compare to LHRS
   if (fFADCHelicity_sbs>=fFADCLow_min && fFADCHelicity_sbs<=fFADCLow_max){
     hel_sbs = 0;
   }
   else if (fFADCHelicity_sbs>=fFADCHigh_min && fFADCHelicity_sbs<=fFADCHigh_max){
     hel_sbs = 1;
   }
   else {
     hel_sbs = 0x1000;
   }
   if (fFADCPatSync_sbs>=fFADCLow_min && fFADCPatSync_sbs<=fFADCLow_max){
     patsync_sbs = 0;
   }
   else if (fFADCPatSync_sbs>=fFADCHigh_min && fFADCPatSync_sbs<=fFADCHigh_max){
     patsync_sbs = 0x10;
   }
   else {
     patsync_sbs = 0x2000;
   }
   if (fFADCTSettle_sbs>=fFADCLow_min && fFADCTSettle_sbs<=fFADCLow_max){
     tsettle_sbs = 0;
   }
   else if (fFADCTSettle_sbs>=fFADCHigh_min && fFADCTSettle_sbs<=fFADCHigh_max){
     tsettle_sbs = 0x100;
   }
   else {
     tsettle_sbs = 0x4000;
   }
   







   UInt_t hroc = fROCinfo[kHel].roc;
   UInt_t len = evdata.GetRocLength(hroc);
   //  if (len <= 4) 
   //    return -1;


   BANKinfo bankinfo;
   //  Set the bank info for the ring buffer entries
   bankinfo.roc = 10;
   bankinfo.bankid = 0xab01;

   lbuff = GetBankBuffer(evdata, bankinfo);
   if (bankinfo.length>2){
      UInt_t nchan    = (lbuff[2]&0x00ff0000)>>16;
      if (nchan==0)  nchan=32;
      UInt_t numread  = lbuff[2]&0xffff;
      if (bankinfo.length != (numread*(nchan+2))+2){
	 std::cout <<"Mismatch between bank length ("<<std::dec<<bankinfo.length
	    <<") and expected from #chan and #reads (totaled=="
	    <<(numread*(nchan+2))+2<<"; #chan+2=="<<nchan+2
	    <<"; #reads==" <<numread <<std::endl;
      } else{
	 UInt_t index = 3;
	 if(fVerbosity>0) std::cout << std::dec << "[SBSScalerHelicityReader::ReadDatabase]: Numread = " << numread << std::endl;
	 fIRing = numread;
	 for (UInt_t i=0; i<numread; i++){
	   UInt_t qrthel     = lbuff[index+1];
	   fTimeStampRing[i] = lbuff[index+0];
	   fHelicityRing[i]  = (qrthel & 0x1);
	   fPatternRing[i]   = (qrthel & 0x10);
	   for (UInt_t j=0; j<nchan; j++){
	     fScalerRing[i][j]=lbuff[index+2+j];
	   }
	   /*
	     std::cout << "\n\nindex : Clock, qrthel, chan0, chan9, chan15:  "
	     << std::dec << index <<" : "<< std::hex 
	     << lbuff[index] << " " << lbuff[index+1] << " "
	     << lbuff[index+2] << " " << lbuff[index+11] 
	     << " " << lbuff[index+17] << std::endl;
	   */
	   index += nchan+2;
	 }
      }
   }
   //  Set the bank info for the helicity summary
   bankinfo.roc = 10;  //LHRS helicity scaler crate
   bankinfo.bankid = 0xab11;

   Bool_t found_lhrs_bank = kFALSE;
   lbuff = GetBankBuffer(evdata, bankinfo);
   if (bankinfo.length==9){
      /* std::cout <<std::hex << lbuff[0] << " " <<lbuff[1]
       *           <<std::endl;
       */
      found_lhrs_bank = kTRUE;
      fHelErrorCond   = lbuff[2];
      fNumEvents      = lbuff[3];
      fNumPatterns    = lbuff[4];
      fPatternPhase   = lbuff[5];
      fSeedValue      = lbuff[6];
      fPatternHel    = (fSeedValue & 0x1);
      fEventPolarity  = lbuff[7];
      fReportedQrtHel = lbuff[8];
      //  lbuff[9] is currently empty.
      /*
	 std::cout << std::dec << "errorcode=="<<fHelErrorCond
	 <<"; #event, #pattern=="<<fNumEvents<<", "<<fNumPatterns
	 <<"; pat_phase=="<<fPatternPhase
	 <<"; helbit=="<<fPatternHel
	 <<"; qrthel=="<<std::hex<<fReportedQrtHel<<std::dec
	 <<"; polarity=="<<fEventPolarity
	 <<std::endl;
	 std::cout <<std::dec << "LHRS FADC data: " 
	 << fFADCSpare    << ", HEL="
	 << fFADCHelicity << ", QRT="
	 << fFADCPatSync  << ", MPS="
	 << fFADCTSettle  << std::endl;
	 std::cout <<std::dec << "SBS FADC data: " 
	 << fFADCSpare_sbs    << ", HEL="
	 << fFADCHelicity_sbs << ", QRT="
	 << fFADCPatSync_sbs  << ", MPS="
	 << fFADCTSettle_sbs  << std::endl;
      */
   }
   // Check the SBS crate copy too.
   bankinfo.roc = 1;  //SBS helciity scaler crate
   bankinfo.bankid = 0xab11;

   lbuff = GetBankBuffer(evdata, bankinfo);
   if (bankinfo.length==9){
      /* std::cout <<std::hex << lbuff[0] << " " <<lbuff[1]
       *           <<std::endl;
       */
      UInt_t sbs_HelErrorCond   = lbuff[2];
      UInt_t sbs_NumEvents      = lbuff[3];
      UInt_t sbs_NumPatterns    = lbuff[4];
      UInt_t sbs_PatternPhase   = lbuff[5];
      UInt_t sbs_SeedValue      = lbuff[6];
      UInt_t sbs_PatternHel     = (sbs_SeedValue & 0x1);
      UInt_t sbs_EventPolarity  = lbuff[7];
      UInt_t sbs_ReportedQrtHel = lbuff[8];

      if (fHelErrorCond!=sbs_HelErrorCond) {
	if(fVerbosity>0) std::cerr << here << " LHRS/SBS Mismatch for fHelErrorCond: "
		  << fHelErrorCond << " " << sbs_HelErrorCond << std::endl;
	fHelErrorCond |= 0x01010000;
      }
      if (fNumEvents!=(sbs_NumEvents-evtnum_offset)) {
	if(fVerbosity>0) std::cerr << here << " LHRS/SBS Mismatch for fNumEvents: "
		  << fNumEvents << " " << sbs_NumEvents << std::endl;
	evtnum_offset = sbs_NumEvents - fNumEvents;
	fHelErrorCond |= 0x01020000;
      }
      if (fNumPatterns!=(sbs_NumPatterns-patnum_offset)) {
	if(fVerbosity>0) std::cerr << here << " LHRS/SBS Mismatch for fNumPatterns: "
		  << fNumPatterns << " " << sbs_NumPatterns << std::endl;
	patnum_offset = sbs_NumPatterns - fNumPatterns;
	fHelErrorCond |= 0x01040000;
      }
      if (fPatternPhase!=sbs_PatternPhase) {
	if(fVerbosity>0) std::cerr << here << " LHRS/SBS Mismatch for fPatternPhase: "
		  << fPatternPhase << " " << sbs_PatternPhase << std::endl;
	fHelErrorCond |= 0x01080000;
      }
      if (fSeedValue!=sbs_SeedValue) {
	if(fVerbosity>0) std::cerr << here << " LHRS/SBS Mismatch for fSeedValue: "
		  << fSeedValue << " " << sbs_SeedValue << std::endl;
	fHelErrorCond |= 0x01100000;
      }
      if (fPatternHel!=sbs_PatternHel) {
	if(fVerbosity>0) std::cerr << here << " LHRS/SBS Mismatch for fPatternHel: "
		  << fPatternHel << " " << sbs_PatternHel << std::endl;
	fHelErrorCond |= 0x01200000;
      }
      if (fEventPolarity!=sbs_EventPolarity) {
	if(fVerbosity>0) std::cerr << here << " LHRS/SBS Mismatch for fEventPolarity: "
		  << fEventPolarity << " " << sbs_EventPolarity << std::endl;
	fHelErrorCond |= 0x01400000;
      }
      if (fReportedQrtHel!=sbs_ReportedQrtHel) {
	if(fVerbosity>0) std::cerr << here << " LHRS/SBS Mismatch for fReportedQrtHel: "
		  << fReportedQrtHel << " " << sbs_ReportedQrtHel << std::endl;
	fHelErrorCond |= 0x01800000;
      }
   }

   //  Compare the FADC helicity bits & set error code if disagree.
   if (hel!=hel_sbs || patsync!=patsync_sbs || tsettle!=tsettle_sbs){
     if(fVerbosity>0) std::cerr << here << " Mismatch between LHRS and SBS FADC helicity bits: "
	       << "(hel, patsync,tsettle) " << hel<<":"<<hel_sbs <<" "
	       << patsync<<":"<<patsync_sbs <<" "
	       << tsettle<<":"<<tsettle_sbs << std::endl;
     fHelErrorCond |= 0x10000000;
   }
   if (fHelErrorCond!=16 && tsettle==0 &&  hel!=(fPatternHel^fEventPolarity)){
     if(fVerbosity>0) std::cerr << here << " Mismatch between FADC helicity bit and value expected from predictor: "
	       << hel << " " << "fPatternHel==" << fPatternHel
	       << " fEventPolarity=="<<fEventPolarity << " "
	       << (fPatternHel^fEventPolarity) << std::endl;
     fHelErrorCond |= 0x20000000;
   }
   if (fHelErrorCond!=16 && tsettle==0 &&  (patsync==16)!=(fPatternPhase==0)){
     if(fVerbosity>0) std::cerr << here << " Mismatch between FADC pattern sync bit and value expected from predictor: "
	       << (patsync==16) << " " << (fPatternPhase==0) << std::endl;
     fHelErrorCond |= 0x40000000;
   }
   if (tsettle!=0){
     //  Set an error bit if we're in the tsettle.  Maybe don't need this...
     fHelErrorCond |= 0x80000000;
   }



   // if(fIRing>kHelRingDepth)
   //   {
   //     ::Error( here, "Ring depth to large ");
   //     return -1;
   //   }
   // for(UInt_t i=0;i<fIRing; i++)
   //   {
   //     fTimeStampRing[i]=evdata.GetRawData(hroc,index++);
   //     data=evdata.GetRawData(hroc,index++);
   //     fHelicityRing[i]=(data & 0x1);
   //     fPatternRing[i]= (data & 0x10)>>4;
   //     fT3Ring[i]=evdata.GetRawData(hroc,index++);
   //     fU3Ring[i]=evdata.GetRawData(hroc,index++);
   //     fT5Ring[i]=evdata.GetRawData(hroc,index++);
   //     fT10Ring[i]=evdata.GetRawData(hroc,index++);
   //   }

   //    Print();

   FillHisto();

   return 0;
}
//____________________________________________________________________
void SBSScalerHelicityReader::FillHisto()
{
   // static const char* here = "SBSScalerHelicityReader::FillHisto";
   // cout<<here<<endl;

   fHistoR[0]->Fill(fPatternTir);
   fHistoR[1]->Fill(fTSettleTir);
   fHistoR[2]->Fill(fHelicityTir);
   fHistoR[3]->Fill(fTimeStampTir-fOldTimeStampTir);
   fOldTimeStampTir=fTimeStampTir;
   for(UInt_t i=0;i<fIRing;i++)
   {
      fHistoR[4]->Fill(fPatternRing[i]);
      fHistoR[5]->Fill(fHelicityRing[i]);
      fHistoR[6]->Fill(fTimeStampRing[i]);
      fHistoR[7]->Fill(fT3Ring[i]);
      fHistoR[8]->Fill(fU3Ring[i]);
      fHistoR[9]->Fill(fT5Ring[i]);
      fHistoR[10]->Fill(fT10Ring[i]);
   }
   fHistoR[11]->Fill(fIRing);
}
//TODO: this should not be needed once LoadDB can fill fROCinfo directly
//____________________________________________________________________
Int_t SBSScalerHelicityReader::SetROCinfo( EROC which, UInt_t roc,
      UInt_t header, UInt_t index )
{

   // Define source and offset of data.  Normally called by ReadDatabase
   // of the detector that is a SBSScalerHelicityReader.
   //
   // "which" is one of  { kHel, kTime, kROC2, kROC3 }.
   // You must define at least the kHel and kTime ROCs.
   // Returns <0 if parameter error, 0 if success

   if( which<kHel || which>kROC3 )
      return -1;
   if( roc <= 0 || roc > 255 )
      return -2;

   fROCinfo[which].roc    = roc;
   fROCinfo[which].header = header;
   fROCinfo[which].index  = index;

   fHaveROCs = ( fROCinfo[kHel].roc > 0 && fROCinfo[kTime].roc > 0 );

   return 0;
}
//____________________________________________________________________
ClassImp(SBSScalerHelicityReader)
