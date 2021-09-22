/* #ifndef ROOT_SBSTimingHodoscopeBar */
/* #define ROOT_SBSTimingHodoscopeBar */
#ifndef SBSTimingHodoscopeBar_h
#define SBSTimingHodoscopeBar_h


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// SBSTimingHodoscopeBar                                                     //
//                                                                           //
// Class to represent a timing hodoscope bar                                 //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

/* class SBSTimingHodoscopeBar; */
//class SBSTimingHodoscopePMT;

#include "TObject.h"
#include "TRef.h"
#include "SBSTimingHodoscopePMT.h"

//380cm light attenuation length from data sheet for EJ200 (but doesn't say for what energy)

class SBSTimingHodoscopeBar : public TObject {
public: 
	SBSTimingHodoscopeBar(
			      Int_t barnum,
			      SBSTimingHodoscopePMT* leftpmt,
			      SBSTimingHodoscopePMT* rightpmt,
			      Int_t baroff);

	virtual ~SBSTimingHodoscopeBar();

	Int_t                  GetBarNum() const {return fBarNum;}
	SBSTimingHodoscopePMT* GetLPMT() { return fLPMT; }
	SBSTimingHodoscopePMT* GetRPMT() { return fRPMT; }
	Int_t                  GerBarOff() const {return fBarOff;}

	void SetBarNum( Int_t barnum ) {fBarNum=barnum;}
	void SetBarOff( Int_t baroff ) {fBarOff=baroff;}
	/* void SetLPMT( SBSTimingHodoscopePMT* leftpmt) {fLPMT=leftpmt;} */
	/* void SetRPMT( SBSTimingHodoscopePMT* rightpmt) {fRPMT=rightpmt;} */
	void ClearEvent();


protected: 
	Int_t fBarNum;
	Int_t fBarOff;
	SBSTimingHodoscopePMT* fLPMT;
	SBSTimingHodoscopePMT* fRPMT;

public:
	ClassDef(SBSTimingHodoscopeBar,5) // Scintillator bar

};
////////////////////////////////////////////////////////////////////////////////

#endif
