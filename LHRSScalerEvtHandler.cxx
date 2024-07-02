//////////////////////////////////////////////////////////////////
//
//   ScalerScalerEvtHandler
//
//   Event handler for Hall A scalers (THaScalerEvtHandler)
//   R. Michaels,  Sept, 2014
//   
//   Lineage: 
//   - Based on TriScalerEvtHandler for tritium experiments (Hanjie Liu, UMass) 
//   - Adapted for SBS experiments (David Flay, JLab) 
//
//   This class does the following
//      For a particular set of event types (here, event type 140)
//      decode the scalers and put some variables into global variables.
//      The global variables can then appear in the Podd output tree T.
//      In addition, a tree "TS" is created by this class; it contains
//      just the scaler data by itself.  Note, the "fName" is concatenated
//      with "TS" to ensure the tree is unqiue; further, "fName" is
//      concatenated with the name of the global variables, for uniqueness.
//      The list of global variables and how they are tied to the
//      scaler module and channels is in the scaler.map file, or could
//      be hardcoded here.   
//      NOTE: if you don't have the scaler map file (e.g. db_LeftScalevt.dat)
//      there will be no variable output to the Trees.
//
//   To use in the analyzer, your setup script needs something like this
//       gHaEvtHandlers->Add (new LHRSScalerEvtHandler("Left","HA scaler event type 140"));
//
//   To enable debugging you may try this in the setup script
// 
//     LHRSScalerEvtHandler *lscaler = new LHRSScalerEvtHandler("Left","HA scaler event type 140");
//     lscaler->SetDebugFile("LeftScaler.txt");
//     gHaEvtHandlers->Add (lscaler);
//
/////////////////////////////////////////////////////////////////////

#include "LHRSScalerEvtHandler.h"

#include "THaAnalysisObject.h"
#include "THaEvtTypeHandler.h"
#include "THaCodaData.h"
#include "THaEvData.h"
#include "THaVarList.h"
#include "THaString.h"
#include "THaAnalyzer.h"

#include "GenScaler.h"
#include "Scaler3800.h"
#include "Scaler3801.h"
#include "Scaler1151.h"
#include "Scaler560.h"

#include "VarDef.h"
#include "Textvars.h" 

#include "TNamed.h"
#include "TMath.h"
#include "TString.h"

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <sstream>

using namespace std;
using namespace Decoder;
using namespace THaString;

static const UInt_t ICOUNT    = 1;
static const UInt_t IRATE     = 2;
static const UInt_t ICURRENT  = 3;
static const UInt_t ICHARGE   = 4;
static const UInt_t ITIME     = 5;
static const UInt_t ICUT      = 6;
static const UInt_t MAXCHAN   = 32;
static const UInt_t MAXTEVT   = 5000;
static const UInt_t defaultDT = 4;

LHRSScalerEvtHandler::LHRSScalerEvtHandler(const char *name, const char* description)
  : THaEvtTypeHandler(name,description),evcount(0),fNormIdx(-1),fNormSlot(-1),
    dvars(0),fScalerTree(0),fUseFirstEvent(kTRUE),
    fClockChan(-1),fClockFreq(-1),fLastClock(0),fClockOverflows(0),
    fTotalTime(0),fPrevTotalTime(0),fDeltaTime(-1),
    dvarsFirst(0),dvars_prev_read(0),fPhysicsEventNumber(-1),
    fNumBCMs(0),fbcm_Current_Threshold_Index(0),fbcm_Current_Threshold(0),
    fBCM_Gain(0),fBCM_Offset(0),fBCM_SatOffset(0),fBCM_SatQuadratic(0),fBCM_delta_charge(0)
{
  rdata = new UInt_t[MAXTEVT];
  fDebugFile = nullptr; // initialize the pointer to null
  // for by hand calculation of rates 
  scal_prev_read.clear();
  scal_present_read.clear();
  scal_overflows.clear(); 
}
//______________________________________________________________________________
LHRSScalerEvtHandler::~LHRSScalerEvtHandler()
{
  delete [] rdata;
  if (fScalerTree) {
    delete fScalerTree;
  }
  // added by D Flay 11/9/21
  delete [] dvars_prev_read;
  delete [] dvars;
  delete [] dvarsFirst;
  delete [] fBCM_Gain;
  delete [] fBCM_Offset;
  delete [] fBCM_SatOffset;
  delete [] fBCM_SatQuadratic;
  delete [] fBCM_delta_charge;
}
//______________________________________________________________________________
Int_t LHRSScalerEvtHandler::End( THaRunBase* r)
{
  if (fScalerTree) fScalerTree->Write();
  //Insert here the addition of summary filling
  THaAnalyzer* analyzer = THaAnalyzer::GetInstance();
  if(analyzer!=nullptr){// check that the analyzer actually exists... otherwise, skip
    const char* summaryfilename = analyzer->GetSummaryFileName();
    cout << "LHRSScalerEvtHandler Summary in " << summaryfilename << endl;
    if( strcmp(summaryfilename,"")!=0  ) {
      ofstream ostr(summaryfilename, std::ofstream::app);
      if( ostr ) {
	// Write to file via cout
	//streambuf* cout_buf = cout.rdbuf();
	//cout.rdbuf(ostr.rdbuf());
	TDatime now;
	ostr << "LHRS scalers Summary " //<< fRun->GetNumber()
	     << " completed " << now.AsString()
	     << endl << " count " << evcount << endl
	     << endl;

	for (UInt_t i = 0; i < scalerloc.size(); i++) {
	  TString name = scalerloc[i]->name; 
	  //tinfo = name + "/D";
	  //fScalerTree->Branch(name.Data(), &dvars[i], tinfo.Data(), 4000);
	  ostr << " Scaler " << name.Data() <<  " value: " << dvars[i] << endl;
	}	
	//std::vector<Decoder::GenScaler*> scalers;
	//std::vector<ScalerVar*> scalerloc;
	ostr << endl;
	
	//cout.rdbuf(cout_buf);
	ostr.close();
	
      }
      
    }
  }
  return 0;
}
//______________________________________________________________________________
Int_t LHRSScalerEvtHandler::Analyze(THaEvData *evdata)
{
  Int_t lfirst=1;

  if ( !IsMyEvent(evdata->GetEvType()) ) return -1;

  if (fDebugFile) {
    *fDebugFile << endl << "---------------------------------- "<<endl<<endl;
    *fDebugFile << "\nEnter LHRSScalerEvtHandler  for fName = "<<fName<<endl;
    EvDump(evdata);
  }

  if (lfirst && !fScalerTree) {

    lfirst = 0; // Can't do this in Init for some reason

    TString sname1 = "TS";
    TString sname2 = sname1 + fName;
    TString sname3 = fName + "  Scaler Data";

    if (fDebugFile) {
      *fDebugFile << "\nAnalyze 1st time for fName = "<<fName<<endl;
      *fDebugFile << sname2 << "      " <<sname3<<endl;
    }

    fScalerTree = new TTree(sname2.Data(),sname3.Data());
    fScalerTree->SetAutoSave(200000000);

    TString name, tinfo;

    name = "evcount";
    tinfo = name + "/D";
    fScalerTree->Branch(name.Data(), &evcount, tinfo.Data(), 4000);

    // create the physics event number branch 
    fScalerTree->Branch("evnum",&fPhysicsEventNumber,"evnum/L");

    for (UInt_t i = 0; i < scalerloc.size(); i++) {
      name = scalerloc[i]->name; 
      tinfo = name + "/D";
      fScalerTree->Branch(name.Data(), &dvars[i], tinfo.Data(), 4000);
    }

  }  // if (lfirst && !fScalerTree)

  // Parse the data, load local data arrays.

  Int_t ndata = evdata->GetEvLength();
  if(fDebugFile) *fDebugFile << "NDATA = " << dec << ndata << std::endl;
  if (ndata >= static_cast<Int_t>(MAXTEVT)) {
	  cout << "LHRSScalerEvtHandler:: ERROR: Event length crazy "<<endl;
	  ndata = MAXTEVT-1;
  }

  if (fDebugFile) *fDebugFile<<"\n\nLHRSScalerEvtHandler :: Debugging event type "<<dec<<evdata->GetEvType()<<endl<<endl;

  // get the physics event number 
  fPhysicsEventNumber = evdata->GetEvNum(); 

  // local copy of data
  // NOTE: event is ASCII, not 32-bit binary! We need to convert ASCII to 32-bit binary  
  for (Int_t i=0; i<ndata; i++) rdata[i] = evdata->GetRawData(i);

  Int_t nskip=0;
  UInt_t *P = rdata;
  // UInt_t *Pstop = rdata+ndata;
  int k=0;

  Int_t ifound=0;
  Int_t itimeout=0;
  UInt_t NScalers=0;

  // Added by D Flay (10/27/21) for parsing data
  std::string word[MAXTEVT]; 
  UInt_t A[MAXTEVT]; 
  int NWORDS=0;

  // skip the first 4 words because it looks like the first word associated with 
  // the scalers starts there... 
  P = P + 4;

  // do the conversion 
  char *pc      = (char *)P;
  NWORDS        = ParseData(pc,word,A);
  UInt_t *p     = A;  
  UInt_t *pstop = p + ndata - 4;  
  if (fDebugFile) *fDebugFile << "number of words: " << NWORDS << endl;
  // char msg[200];  

  AnalyzeBuffer(ndata,rdata); 

  while (p < pstop && k < ndata) {
    if (fDebugFile) {
      // *fDebugFile << "p  and  pstop  "<< k++ << "   " << p << "   " << pstop << "   data = " << hex << *p << "   " << dec << endl;
      *fDebugFile << "k = " << k << " p  = " << p << " pstop = " << pstop << " data = " << hex << *p << " (hex), " << dec << endl;
    }
    nskip = 1;
    itimeout=0;
    NScalers = scalers.size();
    if (fDebugFile)*fDebugFile << "**** NUM SCALERS = " << NScalers << std::endl;
    for (UInt_t j=0; j<NScalers; j++) {
       // bump pointer until scaler found, and don't decode if already found for this event.
       if (scalerloc[j]->found) continue;
       if (fDebugFile) *fDebugFile << "Slot " << scalers[j]->GetSlot() << endl;
       while (p < pstop) {
          if (scalers[j]->IsSlot(*p) == kTRUE) {
             scalerloc[j]->found=kTRUE;
             ifound = 1;
             goto found1;
          }
          p++;
          if (itimeout++ > 5000) { // avoid infinite loop
             std::cout << "LHRSScalerEvtHandler:: cannot find a slot "<< std::endl;
             goto giveup1;
          }
       }
       found1:
	    if(p==pstop && ifound==0) break;
             if (fDebugFile)*fDebugFile << "\n[LHRSScalerEvtHandler::Analyze]: FOUND EVENT 140!" << std::endl;
	    nskip = scalers[j]->Decode(p);
	    if (fDebugFile && nskip > 1) {
		    *fDebugFile << "\n===== Scaler # "<<j<<"     fName = "<<fName<<"   nskip = "<<nskip<<endl;
		    scalers[j]->DebugPrint(fDebugFile);
	    }
	    if (nskip > 1) goto continue1;
    }
    continue1:
       p = p + nskip;
    k++; 
  }
       
  giveup1:
    if (fDebugFile) {
      *fDebugFile << "Finished with decoding.  "<<endl;
      *fDebugFile << "   Found flag   =  "<<ifound<<endl;
    }

  // L-HRS has headers which are different from R-HRS, but both are
  // event type 140 and come here.  If you found no headers, it was
  // the other arms event type.  (The arm is fName).
  if (!ifound) return 0;

  // // The correspondance between dvars and the scaler and the channel
  // // will be driven by a scaler.map file, or could be hard-coded.
  // for (size_t i = 0; i < scalerloc.size(); i++) {
  // 	size_t ivar  = scalerloc[i]->ivar;
  // 	size_t idx   = scalerloc[i]->index;
  // 	size_t ichan = scalerloc[i]->ichan;
  // 	if (fDebugFile) *fDebugFile << "Debug dvars i = "<<i<<", var = "<<ivar<<", index = "<<idx<<", ch = "<<ichan<<endl;
  // 	if( (ivar<scalerloc.size()) && (idx<scalers.size()) && (ichan<MAXCHAN) ){
  // 		if (scalerloc[ivar]->ikind == ICOUNT) dvars[ivar] = scalers[idx]->GetData(ichan);
  // 		if (scalerloc[ivar]->ikind == IRATE)  dvars[ivar] = scalers[idx]->GetRate(ichan);
  // 		if (fDebugFile) *fDebugFile << "   dvars kind = "<<scalerloc[ivar]->ikind<<", value = "<<dvars[ivar]<<endl;
  // 	}else{
  // 		cout << "LHRSScalerEvtHandler:: ERROR:: incorrect index "<<ivar<<"  "<<idx<<"  "<<ichan<<endl;
  // 	}
  // }
  
  // By-hand calculation of rates (added by D Flay 11/9/21) 
  UInt_t thisClock = scalers[fNormIdx]->GetData(fClockChan);

  // FIXME: empirically found maxima (is this right?) 
  UInt_t CLOCK_MAX  = 106680229;
  // UInt_t SCALER_MAX = 816862727;  
  // UInt_t delta=0; // for overflow accumulation

  if(thisClock<fLastClock){  // Count clock scaler wrap arounds
    fClockOverflows++;
    if(fLastClock>CLOCK_MAX) CLOCK_MAX = fLastClock; // update CLOCK_MAX if necessary 
    // delta = kMaxUInt - fLastClock;  
    if(fDebugFile){
       *fDebugFile << "*** CLOCK OVERFLOW! ***" << std::endl;
       *fDebugFile << "cntr       = " << fClockOverflows << std::endl;
       *fDebugFile << "this clock = " << thisClock       << std::endl;
       *fDebugFile << "last clock = " << fLastClock      << std::endl;
       // *fDebugFile << "CLOCK_MAX  = " << CLOCK_MAX       << std::endl;
       // *fDebugFile << "delta      = " << delta           << std::endl;
       *fDebugFile << "kMaxUInt   = " << kMaxUInt        << std::endl;
       *fDebugFile << "***********************" << std::endl;
    }
    // fTotalTime = (thisClock + fLastClock)/fClockFreq; // ( thisClock + ( ( 1. + (Double_t)kMaxUInt )*fClockOverflows - fLastClock ) )/fClockFreq;
  }

  // FIXME 
  fTotalTime = ( thisClock + ( ( (Double_t)fClockOverflows )*kMaxUInt + fClockOverflows ) )/fClockFreq; // this is definitely not right.
  // fTotalTime = ( thisClock + ( ( (Double_t)fClockOverflows )*CLOCK_MAX + fClockOverflows) )/fClockFreq; // this is closer 
  // fTotalTime = ( thisClock + ( ( (Double_t)fClockOverflows )*CLOCK_MAX + (Double_t)fClockOverflows)*delta )/fClockFreq;  // not right
  fDeltaTime = fTotalTime - fPrevTotalTime;

  if(fDebugFile){
     *fDebugFile << "======== Time Check ========" << std::endl;
     *fDebugFile << "Clock frequency = " << fClockFreq     << std::endl;
     *fDebugFile << "Current clock   = " << thisClock      << std::endl;
     *fDebugFile << "Previous clock  = " << fLastClock     << std::endl;
     *fDebugFile << "Current time    = " << fTotalTime     << std::endl;
     *fDebugFile << "Previous time   = " << fPrevTotalTime << std::endl;
     *fDebugFile << "delta time      = " << fDeltaTime     << std::endl;
     *fDebugFile << "============================" << std::endl;
  }

  if(fDeltaTime==0){
     cout << " *******************  Severe Warning   ****************************" << endl;
     cout << " [LHRSScalerEvtHandler]: Found fDeltaTime is zero!!   " << endl;
     cout << " ******************* Alert DAQ experts ****************************" << endl;
  }

  // set up for next event 
  fLastClock     = thisClock;
  fPrevTotalTime = fTotalTime;

  UInt_t scalerData=0;

  Double_t rate=0;
  Double_t scal_current=0;

  Int_t nscal=0;

  for(size_t i=0;i<scalerloc.size();i++){
     size_t ivar  = scalerloc[i]->ivar;
     size_t idx   = scalerloc[i]->index;
     size_t ichan = scalerloc[i]->ichan;
     if (fDebugFile) *fDebugFile << "event " << evcount << " Debug dvars i = "<<i<<", var = "<<ivar<<", index = "<<idx<<", ch = "<<ichan<<endl;
     if(evcount==0){
        if( (ivar<scalerloc.size()) && (idx<scalers.size()) && (ichan<MAXCHAN) ){
           if(fUseFirstEvent){
             scalerData = scalers[idx]->GetData(ichan); 
             if(scalerloc[ivar]->ikind==ICOUNT){
               dvars[ivar] = scalerData;
               scal_present_read.push_back(scalerData);
               scal_prev_read.push_back(0);
               scal_overflows.push_back(0);
               dvarsFirst[ivar] = 0.0;
             }
             if(scalerloc[ivar]->ikind==IRATE){
		scalerData = scalers[idx]->GetData(ichan); 
		rate       = scalerData/fDeltaTime;
                dvars[ivar]      = rate;
                dvarsFirst[ivar] = rate;
                if(fDebugFile){
		   *fDebugFile << "  RATE CALC ivar " << ivar << " diff = " << scalerData 
                               << " dtime = " << fDeltaTime << " rate = " << rate << std::endl;
                }
             }
	     if(scalerloc[ivar]->ikind == ICURRENT || scalerloc[ivar]->ikind == ICHARGE){
                Int_t bcm_ind=-1;
                for(Int_t itemp =0; itemp<fNumBCMs;itemp++){
                   size_t match = string(scalerloc[ivar]->name.Data()).find(string(fBCM_Name[itemp]));
                   if (match!=string::npos){
                        bcm_ind=itemp;
                   }
                }
                if(scalerloc[ivar]->ikind == ICURRENT){
                   dvars[ivar]=0.;
                   if (bcm_ind != -1) {
                        dvars[ivar]=((scalers[idx]->GetData(ichan))/fDeltaTime-fBCM_Offset[bcm_ind])/fBCM_Gain[bcm_ind];
                        dvars[ivar]=dvars[ivar]+fBCM_SatQuadratic[bcm_ind]*TMath::Power(TMath::Max(dvars[ivar]-fBCM_SatOffset[bcm_ind],0.0),2.0);
                   
                   }
                   if (bcm_ind==fbcm_Current_Threshold_Index) scal_current= dvars[ivar];
                }
                if(scalerloc[ivar]->ikind == ICHARGE){
                   if(bcm_ind != -1){
                      Double_t cur_temp=((scalers[idx]->GetData(ichan))/fDeltaTime-fBCM_Offset[bcm_ind])/fBCM_Gain[bcm_ind];
                      cur_temp=cur_temp+fBCM_SatQuadratic[bcm_ind]*TMath::Power(TMath::Max(cur_temp-fBCM_SatOffset[bcm_ind],0.0),2.0);
                      fBCM_delta_charge[bcm_ind]=fDeltaTime*cur_temp;
                      dvars[ivar]+=fBCM_delta_charge[bcm_ind];
                   }
                }
                //      printf("1st event %i index %i fBCMname %s scalerloc %s offset %f gain %f computed %f\n",evcount, bcm_ind, fBCM_Name[bcm_ind],scalerloc[ivar]->name.Data(),fBCM_Offset[bcm_ind],fBCM_Gain[bcm_ind],dvars[ivar]);
                //                
	     }
           }else{
              // not using first event
              if(scalerloc[ivar]->ikind==ICOUNT){
		 scalerData = scalers[idx]->GetData(ichan); 
                 dvarsFirst[ivar] = scalerData;
                 scal_present_read.push_back(dvarsFirst[ivar]);
                 scal_prev_read.push_back(0);
              }
              if(scalerloc[ivar]->ikind==IRATE){
		 scalerData = scalers[idx]->GetData(ichan); 
		 rate       = scalerData/fDeltaTime;
                 dvarsFirst[ivar] = rate;
                 if(fDebugFile) *fDebugFile << "  RATE CALC ivar " << ivar << " diff = " << scalerData << " dtime = " << fDeltaTime << " rate = " << rate << std::endl;
              }
	      if(scalerloc[ivar]->ikind==ICURRENT || scalerloc[ivar]->ikind==ICHARGE){
                 Int_t bcm_ind=-1;
                 for(Int_t itemp =0; itemp<fNumBCMs;itemp++){
                    size_t match = string(scalerloc[ivar]->name.Data()).find(string(fBCM_Name[itemp]));
                    if(match!=string::npos){
                         bcm_ind=itemp;
                    }
                 }
                 if(scalerloc[ivar]->ikind == ICURRENT){
                    dvarsFirst[ivar]=0.0;
                    if(bcm_ind != -1){
                       dvarsFirst[ivar]=((scalers[idx]->GetData(ichan))/fDeltaTime-fBCM_Offset[bcm_ind])/fBCM_Gain[bcm_ind];
                       dvarsFirst[ivar]=dvarsFirst[ivar]+fBCM_SatQuadratic[bcm_ind]*TMath::Power(TMath::Max(dvars[ivar]-fBCM_SatOffset[bcm_ind],0.0),2.);
                    }
                    if(bcm_ind==fbcm_Current_Threshold_Index) scal_current= dvarsFirst[ivar];
                 }
                 if(scalerloc[ivar]->ikind == ICHARGE){
                    if(bcm_ind != -1){
                       Double_t cur_temp=((scalers[idx]->GetData(ichan))/fDeltaTime-fBCM_Offset[bcm_ind])/fBCM_Gain[bcm_ind];
                       cur_temp=cur_temp+fBCM_SatQuadratic[bcm_ind]*TMath::Power(TMath::Max(cur_temp-fBCM_SatOffset[bcm_ind],0.0),2.);
                       fBCM_delta_charge[bcm_ind]=fDeltaTime*cur_temp;
                       dvarsFirst[ivar]+=fBCM_delta_charge[bcm_ind];
                    }
                 }
	      }
           }
        }
     }else{
	// evcount != 0
        if( (ivar<scalerloc.size()) && (idx<scalers.size()) && (ichan<MAXCHAN) ){
           if(scalerloc[ivar]->ikind==ICOUNT) {
	      scalerData = scalers[idx]->GetData(ichan);
	      rate       = 0; 
              if(scalerData<scal_prev_read[nscal]){
                scal_overflows[nscal]++;
                if(fDebugFile){
                   *fDebugFile << "*** OVERFLOW ENCOUNTERED! ***" << std::endl;
                   *fDebugFile << "scal_overflows[" << nscal << "] = " << scal_overflows[nscal] << std::endl;
                   *fDebugFile << "scal_prev_read[" << nscal << "] = " << scal_prev_read[nscal] << std::endl;
                   *fDebugFile << "scalerData = " << scalerData << std::endl;
                   *fDebugFile << "kMaxUInt = " << kMaxUInt << std::endl; 
                   // *fDebugFile << "SCALER_MAX = " << SCALER_MAX << std::endl; 
                   *fDebugFile << "*****************************" << std::endl; 
                }
                dvars[ivar] = scalerData + (1+((Double_t)kMaxUInt))*scal_overflows[nscal] - dvarsFirst[ivar];
                // SCALER_MAX = scal_prev_read[nscal]; // ok...
                // if(scal_prev_read[nscal]>SCALER_MAX) SCALER_MAX = scal_prev_read[nscal];  
                // if(scal_prev_read[nscal]>kMaxUInt){
                //    dvars[ivar] = scalerData + (1+((Double_t)kMaxUInt))*scal_overflows[nscal] - dvarsFirst[ivar];
                // }else{
		//    // dvars[ivar] = scalerData + scal_prev_read[nscal];
                //    dvars[ivar] = scalerData + (1+((Double_t)SCALER_MAX))*scal_overflows[nscal] - dvarsFirst[ivar];
                // }
              }else{
                dvars[ivar] = scalerData;
              }
              scal_present_read[nscal] = dvars[ivar]; // scalerData;
              nscal++;
           }
           if(scalerloc[ivar]->ikind==IRATE){
	      scalerData  = scalers[idx]->GetData(ichan);
	      rate        = 0; 
              UInt_t diff = 0;
              if(scalerData<scal_prev_read[nscal-1]){
                if(scal_prev_read[nscal-1]>kMaxUInt){
                   diff = (kMaxUInt-(scal_prev_read[nscal-1] - 1)) + scalerData;
                }else{
		   diff = (scal_prev_read[nscal-1] - 1) + scalerData; 
                   // diff = (SCALER_MAX-(scal_prev_read[nscal-1] - 1)) + scalerData;  
                }
                if(fDebugFile){
                   *fDebugFile << "*** OVERFLOW ENCOUNTERED! ***" << std::endl; 
                   *fDebugFile << "scal_prev_read[" << nscal-1 << "] = " << scal_prev_read[nscal-1] << std::endl; 
                   *fDebugFile << "scalerData = " << scalerData << std::endl; 
                   *fDebugFile << "diff = " << diff << std::endl; 
                   *fDebugFile << "kMaxUInt = " << kMaxUInt << std::endl; 
                   // *fDebugFile << "SCALER_MAX = " << SCALER_MAX << std::endl; 
                   *fDebugFile << "*****************************" << std::endl; 
                }
              }else{
                diff = scalerData - scal_prev_read[nscal-1];
              }
              rate        = diff/fDeltaTime;
              dvars[ivar] = rate;
              if(fDebugFile){
                 *fDebugFile << "  RATE CALC ivar " << ivar << " scalerData = " << scalerData 
                             << " scal_prev_read = " << scal_prev_read[nscal-1] 
                             << " diff = " << diff 
                             << " dtime = " << fDeltaTime << " rate = " << rate << std::endl;
              }
           }
	   if(scalerloc[ivar]->ikind == ICURRENT || scalerloc[ivar]->ikind == ICHARGE)
	   {
	      Int_t bcm_ind=-1;
	      for(Int_t itemp =0; itemp<fNumBCMs;itemp++)
	      {
		 size_t match = string(scalerloc[ivar]->name.Data()).find(string(fBCM_Name[itemp]));
		 if (match!=string::npos)
		 {
		    bcm_ind=itemp;
		 }
	      }
	      if (scalerloc[ivar]->ikind == ICURRENT) {
		 dvars[ivar]=0;
		 if (bcm_ind != -1) {
		    UInt_t scaldata = scalers[idx]->GetData(ichan);
		    UInt_t diff;
		    if(scaldata < scal_prev_read[nscal-1]) {
		       diff = (kMaxUInt-(scal_prev_read[nscal-1] - 1)) + scaldata;
		    } else {
		       diff = scaldata - scal_prev_read[nscal-1];
		    }
		    dvars[ivar]=0.;
		    if (fDeltaTime>0) {
		       Double_t cur_temp=(diff/fDeltaTime-fBCM_Offset[bcm_ind])/fBCM_Gain[bcm_ind];
		       cur_temp=cur_temp+fBCM_SatQuadratic[bcm_ind]*TMath::Power(TMath::Max(cur_temp-fBCM_SatOffset[bcm_ind],0.0),2.);

		       dvars[ivar]=cur_temp;
		    }
		 }
		 if (bcm_ind == fbcm_Current_Threshold_Index) scal_current= dvars[ivar];
	      }
	      if (scalerloc[ivar]->ikind == ICHARGE) {
		 if (bcm_ind != -1) {
		    UInt_t scaldata = scalers[idx]->GetData(ichan);
		    UInt_t diff;
		    if(scaldata < scal_prev_read[nscal-1]) {
		       diff = (kMaxUInt-(scal_prev_read[nscal-1] - 1)) + scaldata;
		    } else {
		       diff = scaldata - scal_prev_read[nscal-1];
		    }
		    fBCM_delta_charge[bcm_ind]=0;
		    if (fDeltaTime>0)  {
		       Double_t cur_temp=(diff/fDeltaTime-fBCM_Offset[bcm_ind])/fBCM_Gain[bcm_ind];
		       cur_temp=cur_temp+fBCM_SatQuadratic[bcm_ind]*TMath::Power(TMath::Max(cur_temp-fBCM_SatOffset[bcm_ind],0.0),2.);
		       fBCM_delta_charge[bcm_ind]=fDeltaTime*cur_temp;
		    }
		    dvars[ivar]+=fBCM_delta_charge[bcm_ind];
		 }
	      }
	   }
	   if (fDebugFile) *fDebugFile << "   dvars  "<<scalerloc[ivar]->ikind<<"  "<<dvars[ivar]<<endl;
	}
     } // end of evcount if-else  
     // if(fDebugFile) *fDebugFile << "ivar " << ivar << " counts = " << scalerData << " rate = " << rate << std::endl;
  } // end of for loop 

  evcount = evcount + 1.0;

  // set up for next read 
  for(size_t j=0;j<scal_prev_read.size();j++) scal_prev_read[j]=scal_present_read[j];
  
  for(size_t j=0;j<scalers.size();j++){
     scalers[j]->Clear("");
     scalerloc[j]->found=kFALSE;
  }
  
  if (fDebugFile) *fDebugFile << "scaler tree ptr  "<<fScalerTree<<endl;
  
  if (fScalerTree) fScalerTree->Fill();

  
  return 1;
}
//______________________________________________________________________________
Int_t LHRSScalerEvtHandler::AnalyzeBuffer(Int_t ndata,UInt_t *rdata){
   // added by D. Flay to analyze data
   UInt_t *P = rdata;

   // rdata is actually ASCII, have to convert 
   
   std::string word[MAXTEVT]; 
   UInt_t A[MAXTEVT]; 

   // skip the first 4 words because it looks like the first word associated with 
   // the scalers starts there... 
   P = P + 4;

   // do the conversion 
   char *pc      = (char *)P;
   int NWORDS    = ParseData(pc,word,A);
   UInt_t *p     = A;  
   // UInt_t *pstop = p + *p - 4;
    
   char msg[200];  

   if(fDebugFile){
      *fDebugFile << "========== D FLAY TEST FUNCTION ==========" << std::endl;
      *fDebugFile << "**** parsed int array = " << p << ", NWORDS = " << dec << NWORDS << endl; 
      for(int ii=0;ii<NWORDS;ii++){
         sprintf(msg,"   word index i = %03d, word = %s, int = %u, hex = %02x",ii,word[ii].c_str(),p[ii],p[ii]); 
         *fDebugFile << msg << endl; 
      }
   }

   // if(fDebugFile){
   //    *fDebugFile << "--- Increment through pointer ---" << std::endl;
   //    *fDebugFile << "start addr = " << p     << std::endl;
   //    *fDebugFile << "end addr = "   << pstop << std::endl;
   // }

   // int NS = scalers.size(); 
   // int k=0;
   // while(p<pstop && k < ndata){
   //    if(fDebugFile){
   //       *fDebugFile << "   ptr addr = " << p << std::endl;  
   //    }
   //    // loop over scalers
   //    for(int j=0;j<NS;j++){
   //       if(scalerloc[j]->found) continue;
   //       
   //    } 
   //    p++; // increment pointer 
   //    k++;  
   // }

   if(fDebugFile) *fDebugFile << "========== END D FLAY TEST FUNCTION ==========" << std::endl;

   return 0;
}
//______________________________________________________________________________
THaAnalysisObject::EStatus LHRSScalerEvtHandler::Init(const TDatime& date)
{

  ReadDatabase(date);

  const int LEN = 200;
  char cbuf[LEN];

  fStatus = kOK;
  fNormIdx = -1;

  // std::cout << "[LHRSScalerEventHandler::Init]: Initializing " << fName << "..." << std::endl;

  eventtypes.push_back(140);  // what events to look for

  TString dfile;
  dfile = fName + "scaler.txt";

  // Parse the map file which defines what scalers exist and the global variables.
  TString sname0 = "Scalevt";
  TString sname;
  sname = fName+sname0;

  FILE *fi = Podd::OpenDBFile(sname.Data(), date);
  if ( !fi ) {
     cout << "Cannot find db file for "<<fName<<" (file = " << sname << ") scaler event handler"<<endl;
     return kFileError;
  }

  string::size_type minus1 = string::npos;
  string::size_type pos1;
  const string scomment = "#";
  const string svariable = "variable";
  const string smap = "map";
  vector<string> dbline;

  while( fgets(cbuf, LEN, fi) != NULL) {
     std::string sinput(cbuf);
     if (fDebugFile) *fDebugFile << "string input "<<sinput<<endl;
     dbline = Podd::vsplit(sinput);
     if (dbline.size() > 0) {
             pos1 = FindNoCase(dbline[0],scomment);
             if (pos1 != minus1) continue;
             pos1 = FindNoCase(dbline[0],svariable);
             if (pos1 != minus1 && dbline.size()>4) {
                string sdesc = "";
                for (UInt_t j=5; j<dbline.size(); j++) sdesc = sdesc+" "+dbline[j];
                Int_t islot = atoi(dbline[1].c_str());
                Int_t ichan = atoi(dbline[2].c_str());
                Int_t ikind = atoi(dbline[3].c_str());
                if (fDebugFile)
                        *fDebugFile << "add var "<<dbline[1]<<"   desc = "<<sdesc<<"    islot= "<<islot<<"  "<<ichan<<"  "<<ikind<<endl;
                TString tsname(dbline[4].c_str());
                TString tsdesc(sdesc.c_str());
                AddVars(tsname,tsdesc,islot,ichan,ikind);
             }
             pos1 = FindNoCase(dbline[0],smap);
             if (pos1 != minus1 && dbline.size()>6) {
           	  Int_t imodel, icrate, islot, inorm;
           	  UInt_t header, mask;
           	  char cdum[20];
           	  sscanf(sinput.c_str(),"%s %d %d %d %x %x %d \n",cdum,&imodel,&icrate,&islot, &header, &mask, &inorm);
           	  if (fNormSlot >= 0 && fNormSlot != inorm) cout << "LHRSScalerEvtHandler:: WARN: contradictory norm slot "<<inorm<<endl;
           	  fNormSlot = inorm;  // slot number used for normalization.  This variable is not used but is checked.
           	  Int_t clkchan = -1;
           	  Double_t clkfreq = 1;
           	  if (dbline.size()>8) {
           		  clkchan = atoi(dbline[7].c_str());
           		  clkfreq = 1.0*atoi(dbline[8].c_str());
			  // save to the class's private variables 
			  fClockChan = clkchan; 
			  fClockFreq = clkfreq; 
           	  }
           	  if (fDebugFile) {
           		  *fDebugFile << "map line "<<dec<<imodel<<"  "<<icrate<<"  "<<islot<<endl;
           		  *fDebugFile <<"   header  0x"<<hex<<header<<"  0x"<<mask<<dec<<"  "<<inorm<<"  "<<clkchan<<"  "<<clkfreq<<endl;
           	  }
           	  switch (imodel) {
           		  case 560:
           			  scalers.push_back(new Scaler560(icrate, islot));
           			  break;
           		  case 1151:
           			  scalers.push_back(new Scaler1151(icrate, islot));
           			  break;
           		  case 3800:
           			  scalers.push_back(new Scaler3800(icrate, islot));
           			  break;
           		  case 3801:
           			  scalers.push_back(new Scaler3801(icrate, islot));
           			  break;
           		  default: 
           			  std::cout << "LHRSScalerEvtHandler:: ERROR: Invalid model " << imodel << std::endl;
           	  }
           	  if (scalers.size() > 0) {
           		  UInt_t idx = scalers.size()-1;
           		  scalers[idx]->SetHeader(header, mask);
           		  // The normalization slot has the clock in it, so we automatically recognize it.
           		  // fNormIdx is the index in scaler[] and 
           		  // fNormSlot is the slot#, checked for consistency
           		  if (clkchan >= 0) {  
           			  int clk_rc = scalers[idx]->SetClock(defaultDT, clkchan, clkfreq);
           			  fNormIdx = idx;
           			  if (islot != fNormSlot) cout << "LHRSScalerEvtHandler:: WARN: contradictory norm slot ! "<<islot<<endl;  
           			  if (fDebugFile) *fDebugFile <<"Setting scaler clock: dt = " << defaultDT <<", channel = "<<clkchan<<", freq = "<<clkfreq<<", fNormIdx = "<<fNormIdx<<", fNormSlot = "<<fNormSlot<<", slot = "<<islot<<", SetClock return value = "<<clk_rc<<endl;
           		  }
           	  }
             }
     }
  }

  // need to do LoadNormScaler after scalers created and if fNormIdx found.
  Int_t nscalers = static_cast<Int_t>(scalers.size());
  if ( fNormIdx >= 0 && fNormIdx < nscalers ) {
     for (Int_t i = 0; i < nscalers; i++) {
        if (i==fNormIdx) continue;
        scalers[i]->LoadNormScaler(scalers[fNormIdx]);
        if(fDebugFile) *fDebugFile << "==> Scaler " << i << ": Loaded normalization scaler ptr = " << scalers[fNormIdx] << std::endl;
     }
  }

#ifdef HARDCODED
  // This code is superseded by the parsing of a map file above.  It's another way ...
  if (fName == "Left") {
	  AddVars("TSbcmu1", "BCM x1 counts", 1, 4, ICOUNT);
	  AddVars("TSbcmu1r","BCM x1 rate",  1, 4, IRATE);
	  AddVars("TSbcmu3", "BCM u3 counts", 1, 5, ICOUNT);
	  AddVars("TSbcmu3r", "BCM u3 rate",  1, 5, IRATE);
  } else {
	  AddVars("TSbcmu1", "BCM x1 counts", 0, 4, ICOUNT);
	  AddVars("TSbcmu1r","BCM x1 rate",  0, 4, IRATE);
	  AddVars("TSbcmu3", "BCM u3 counts", 0, 5, ICOUNT);
	  AddVars("TSbcmu3r", "BCM u3 rate",  0, 5, IRATE);
  }
#endif

DefVars(); 

#ifdef HARDCODED
   // This code is superseded by the parsing of a map file above.  It's another way ...
   if (fName == "Left") {
   	scalers.push_back(new Scaler1151(1,0));
   	scalers.push_back(new Scaler3800(1,1));
   	scalers.push_back(new Scaler3800(1,2));
   	scalers.push_back(new Scaler3800(1,3));
   	scalers[0]->SetHeader(0xabc00000, 0xffff0000);
   	scalers[1]->SetHeader(0xabc10000, 0xffff0000);
   	scalers[2]->SetHeader(0xabc20000, 0xffff0000);
   	scalers[3]->SetHeader(0xabc30000, 0xffff0000);
   	scalers[0]->LoadNormScaler(scalers[1]);
   	scalers[1]->SetClock(4, 7, 1024);
   	scalers[2]->LoadNormScaler(scalers[1]);
   	scalers[3]->LoadNormScaler(scalers[1]);
   } else {
   	scalers.push_back(new Scaler3800(2,0));
   	scalers.push_back(new Scaler3800(2,0));
   	scalers.push_back(new Scaler1151(2,1));
   	scalers.push_back(new Scaler1151(2,2));
   	scalers[0]->SetHeader(0xceb00000, 0xffff0000);
   	scalers[1]->SetHeader(0xceb10000, 0xffff0000);
   	scalers[2]->SetHeader(0xceb20000, 0xffff0000);
   	scalers[3]->SetHeader(0xceb30000, 0xffff0000);
   	scalers[0]->SetClock(4, 7, 1024);
   	scalers[1]->LoadNormScaler(scalers[0]);
   	scalers[2]->LoadNormScaler(scalers[0]);
   	scalers[3]->LoadNormScaler(scalers[0]);
   }
#endif

   // Verify that the slots are not defined twice
   for (UInt_t i1=0; i1 < scalers.size()-1; i1++) {
      for (UInt_t i2=i1+1; i2 < scalers.size(); i2++) {
         if (scalers[i1]->GetSlot()==scalers[i2]->GetSlot())
              cout << "LHRSScalerEvtHandler:: WARN:  same slot defined twice"<<endl;
      }
   }
   // Identify indices of scalers[] vector to variables.
   for (UInt_t i=0; i < scalers.size(); i++) {
      for (UInt_t j = 0; j < scalerloc.size(); j++) {
         if (scalerloc[j]->islot==static_cast<UInt_t>(scalers[i]->GetSlot()))
              scalerloc[j]->index = i;
      }
   }

   if(fDebugFile) {
      *fDebugFile << "LHRSScalerEvtHandler:: Name of scaler bank "<<fName<<endl;
      for (UInt_t i=0; i<scalers.size(); i++) {
         *fDebugFile << "Scaler  #  "<<i<<endl;
         scalers[i]->SetDebugFile(fDebugFile);
         scalers[i]->DebugPrint(fDebugFile);
      }
   }
   for (size_t j=0; j<scalers.size(); j++) {
      scalers[j]->Clear("");
      scalerloc[j]->found=kFALSE;
   }

   return kOK;
}
//______________________________________________________________________________
void LHRSScalerEvtHandler::AddVars(TString name, TString desc, Int_t islot,
				  Int_t ichan, Int_t ikind)
{
  // need to add fName here to make it a unique variable.  (Left vs Right HRS, for example)
  // TString name1 = fName + name;
  TString name1 = Form("%s.%s",fName.Data(),name.Data());
  TString desc1 = fName + desc;
  // We don't yet know the correspondence between index of scalers[] and slots.
  // Will put that in later.
  ScalerVar *loc = new ScalerVar(name1, desc1, 0, islot, ichan, ikind);
  loc->ivar = scalerloc.size();  // ivar will be the pointer to the dvars array.
  scalerloc.push_back(loc);
}
//______________________________________________________________________________
void LHRSScalerEvtHandler::DefVars()
{
  // called after AddVars has finished being called.
  Int_t Nvars = scalerloc.size();
  if (Nvars == 0) return;
  dvars           = new Double_t[Nvars];  // dvars is a member of this class
  dvarsFirst      = new Double_t[Nvars];  // dvarsFirst is a member of this class
  dvars_prev_read = new UInt_t[Nvars];    // dvars_prev_read is a member of this class
  memset(dvars, 0, Nvars*sizeof(Double_t));
  memset(dvarsFirst, 0, Nvars*sizeof(Double_t));
  memset(dvars_prev_read, 0, Nvars*sizeof(UInt_t));
  if (gHaVars) {
  	if(fDebugFile) *fDebugFile << "LHRSScalerEvtHandler:: Have gHaVars "<<gHaVars<<endl;
  } else {
  	cout << "No gHaVars ?!  Well, that's a problem !!"<<endl;
  	return;
  }
  if(fDebugFile) *fDebugFile << "LHRSScalerEvtHandler:: scalerloc size "<<scalerloc.size()<<endl;
  const Int_t* count = 0;
  for (UInt_t i = 0; i < scalerloc.size(); i++) {
  	gHaVars->DefineByType(scalerloc[i]->name.Data(), scalerloc[i]->description.Data(),
  			&dvars[i], kDouble, count);
  }
}
//______________________________________________________________________________
Int_t LHRSScalerEvtHandler::ParseData(char *msg,std::string *word,UInt_t *word_int){
   // loop through the message (msg) and convert into data words 
   // - input:  a char array to parse (i.e., scaler data)  
   // - output: std::string array (word) and int array (word_int)     
   char data[200],subword[200];
   strcpy(data,"");
   strcpy(subword,"");

   // std::cout << "Message to decode: " << std::endl;
   // std::cout << msg << std::endl;

   char *pEnd;
   std::string myStr;

   int j=0;
   int length = strlen(msg);
   for(int i=0;i<length;i++){
      if(msg[i]=='\n'){
	 // now have a full word
	 word[j]     = data;
	 // determine if this is the header
	 for(int k=0;k<3;k++){
	    sprintf(subword,"%s%c",subword,data[k]);
	 }
	 myStr = subword;
	 if(myStr.compare("abc")==0){
	    // this is the header
	    word_int[j] = std::strtoul(data,&pEnd,16);
	 }else{
	    // this is the scaler counts
	    word_int[j] = std::strtoul(data,&pEnd,10);
	 }
	 // increment the index on the word array
	 j++;
	 // empty the constructed word 
	 strcpy(data,"");
	 strcpy(subword,"");
      }else{
	 // not a new line, build the word  
	 sprintf(data,"%s%c",data,msg[i]);
      }
   }

   return j; // return the number of words 
}
//______________________________________________________________________________
Int_t LHRSScalerEvtHandler::ReadDatabase(const TDatime& date){
   char prefix[2];
   prefix[0]='g';
   prefix[1]='\0';
   fNumBCMs = 0;

// #ifdef HALLCPARM

   DBRequest list [] = {
      {"NumBCMs",&fNumBCMs,kInt,0,1},
      {0}
   };

  TString sname = "db_LeftBCM.dat";
  std::cout << "Trying to load database file " << sname << std::endl;

  FILE *file = Podd::OpenDBFile(sname.Data(), date);
  // FILE* file = OpenFile( date );
   if( !file )
      return kInitError;

   Int_t err = kOK;

   if(!err){
      err = LoadDB( file, date,list,fPrefix);
   }
   // DBRequest list[]={
   //    {"NumBCMs",&fNumBCMs, kInt, 0, 1},
   //    {0}
   // };
   // gHcParms->LoadParmValues((DBRequest*)&list, prefix);
   std::cout << "[LHRSScalerEvtHandler::ReadDatabase]: Number of BCMs = " << fNumBCMs << std::endl;

   if(fNumBCMs>0) {
     fBCM_Gain         = new Double_t[fNumBCMs];
     fBCM_Offset       = new Double_t[fNumBCMs];
     fBCM_SatOffset    = new Double_t[fNumBCMs];
     fBCM_SatQuadratic = new Double_t[fNumBCMs];
     fBCM_delta_charge = new Double_t[fNumBCMs];
     std::string bcm_namelist;
     DBRequest list2[]={ 
       {"BCM_Gain"                   , fBCM_Gain                    , kDouble, (UInt_t) fNumBCMs    },
       {"BCM_Offset"                 , fBCM_Offset                  , kDouble, (UInt_t) fNumBCMs    },
       {"BCM_SatQuadratic"           , fBCM_SatQuadratic            , kDouble, (UInt_t) fNumBCMs, 1 },
       {"BCM_SatOffset"              , fBCM_SatOffset               , kDouble, (UInt_t) fNumBCMs, 1 },
       {"BCM_Names"                  , &bcm_namelist                , kString                       },
       {"BCM_Current_threshold"      , &fbcm_Current_Threshold      , kDouble, 0                , 1 },
       {"BCM_Current_threshold_index", &fbcm_Current_Threshold_Index, kInt   , 0                , 1 },
       {0}
     };
     fbcm_Current_Threshold = 0.0;
     fbcm_Current_Threshold_Index = 0;
     for(Int_t i=0;i<fNumBCMs;i++) {
       fBCM_SatOffset[i]=0.;
       fBCM_SatQuadratic[i]=0.;
     }
     err = LoadDB(file,date,list2,fPrefix);
     // gHcParms->LoadParmValues((DBRequest*)&list2, prefix);
     string myStr;
     std::vector<string> bcm_names = Podd::vsplit(bcm_namelist);
     for(Int_t i=0;i<fNumBCMs;i++) {
	myStr = "Left.bcm." + bcm_names[i] + ".current";
	fBCM_Name.push_back(myStr);
	fBCM_delta_charge[i]=0.;
     }
     // print what we have
     std::cout << "LOADED FROM DATABASE: " << std::endl;
     for(Int_t i=0;i<fNumBCMs;i++){
	std::cout << Form("%s: offset = %.3lf Hz, gain = %.3lf Hz/uA",fBCM_Name[i].c_str(),fBCM_Offset[i],fBCM_Gain[i]) << std::endl;
     }
   }
// #endif

   fTotalTime=0.;
   fPrevTotalTime=0.;
   fDeltaTime=-1.;

   return kOK; 
}
//______________________________________________________________________________
ClassImp(LHRSScalerEvtHandler)
