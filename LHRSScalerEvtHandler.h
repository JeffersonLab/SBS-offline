#ifndef LHRSScalerEvtHandler_H
#define LHRSScalerEvtHandler_H

///////////////////////////////////////////////////////////////////
//
//   LHRSScalerEvtHandler
//   Class to handle Hall A scaler events (type 140)
//   Adapted by David Flay (flay@jlab.org), based on code from:  
//   - THaEvtHandler (author Robert Michaels, rom@jlab.org )
//   - TriScalerEvtHandler (author Hanjie Liu, hanjieliu@umass.edu )
//
/////////////////////////////////////////////////////////////////////

#include "THaEvtTypeHandler.h"
#include "Decoder.h"

#include <string>
#include <vector>

#include "TTree.h"
#include "TString.h"  

class ScalerVar { // Utility class used by LHRSScalerEvtHandler
public:
	ScalerVar(TString nm, TString desc, Int_t idx, Int_t sl, Int_t ich, Int_t iki) :
		name(nm), description(desc), index(idx), islot(sl), ichan(ich), ikind(iki) { };
	~ScalerVar();
	TString name, description;
	UInt_t index, islot, ichan, ivar, ikind;
	Bool_t found;
};

class LHRSScalerEvtHandler : public THaEvtTypeHandler {

public:

   LHRSScalerEvtHandler(const char* name, const char* description);
   virtual ~LHRSScalerEvtHandler();

   virtual Int_t Analyze(THaEvData *evdata);
   virtual EStatus Init( const TDatime& run_time);
   virtual Int_t End( THaRunBase* r=0 );


private:

   void AddVars(TString name, TString desc, Int_t iscal, Int_t ichan, Int_t ikind);
   void DefVars();

   Int_t ParseData(char *msg,std::string *word,UInt_t *word_int);
   Int_t AnalyzeBuffer(Int_t ndata,UInt_t *rdata);
   Int_t ReadDatabase(const TDatime& date); 

   std::vector<Decoder::GenScaler*> scalers;
   std::vector<ScalerVar*> scalerloc;
   Double_t evcount;
   UInt_t *rdata;
   Int_t fNormIdx, fNormSlot;
   Double_t *dvars;
   TTree *fScalerTree;

   // added by D Flay 
   Bool_t fUseFirstEvent; 
   Int_t fClockChan;
   Double_t fClockFreq;
   UInt_t fLastClock;
   Int_t fClockOverflows;
   Double_t fTotalTime;
   Double_t fPrevTotalTime;
   Double_t fDeltaTime;
   Double_t *dvarsFirst;
   UInt_t *dvars_prev_read;
   Long64_t fPhysicsEventNumber;
   Int_t fNumBCMs;
   Int_t fbcm_Current_Threshold_Index;
   Double_t fbcm_Current_Threshold;
   Double_t *fBCM_Gain;
   Double_t *fBCM_Offset;
   Double_t *fBCM_SatOffset;
   Double_t *fBCM_SatQuadratic;
   Double_t *fBCM_delta_charge;
   std::vector<std::string> fBCM_Name;
   std::vector<UInt_t> scal_prev_read;
   std::vector<UInt_t> scal_present_read;
   std::vector<UInt_t> scal_overflows;

   LHRSScalerEvtHandler(const LHRSScalerEvtHandler& fh);
   LHRSScalerEvtHandler& operator=(const LHRSScalerEvtHandler& fh);

   ClassDef(LHRSScalerEvtHandler,0)  // Scaler Event handler

};

#endif
