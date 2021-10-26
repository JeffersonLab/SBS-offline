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

   std::vector<Decoder::GenScaler*> scalers;
   std::vector<ScalerVar*> scalerloc;
   Double_t evcount;
   UInt_t *rdata;
   Int_t fNormIdx, fNormSlot;
   Double_t *dvars;
   TTree *fScalerTree;

   LHRSScalerEvtHandler(const LHRSScalerEvtHandler& fh);
   LHRSScalerEvtHandler& operator=(const LHRSScalerEvtHandler& fh);

   ClassDef(LHRSScalerEvtHandler,0)  // Scaler Event handler

};

#endif
