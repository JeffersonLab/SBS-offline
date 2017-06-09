//////////////////////////////////////////////////////////////////////////
//                                                                          
// SBSTimingHodoscope                                                            
//                                                                          
// Class for trigger plane  consisting of multiple                         
// paddles with phototubes on both ends.                                    
//                                                                          
//////////////////////////////////////////////////////////////////////////
//	
//	Author : Copy from AGen Lib
//	Modify History:
//		Jin Huang <mailto:jinhuang@jlab.org>    July 2007	
//			make database file supporting comment
//			disable Multi Hit which is for neutron detection
//
//////////////////////////////////////////////////////////////////////////

#ifndef ROOT_SBSTimingHodoscope
#define ROOT_SBSTimingHodoscope


#include "THaNonTrackingDetector.h"
#include "THaSubDetector.h"
#include "THaApparatus.h"
#include "TClonesArray.h"

//#include "THaNeutronApp.h"

class THaScintPMT;
class THaScintHit;
class THaScintBar;
class THaAdcHit;
class THaTdcHit;
class THaSubDetector;
class THaPartialHit;
//class THaMultiHit;
//class THaNeutronDetector;


//------------------------------------------------------//
//
//	Debug Defines
//	place this section below any other head files
//
//------------------------------------------------------//
//
#ifdef DEBUG_LEVEL
#undef DEBUG_LEVEL
#endif
//	DEBUG_LEVEL;	
//	=0	or not define: no debug, full speed 
//	>=1	enable debug extra warning (suggested setting)
//	>=2	above + enable debug assert
//	>=3	above + enable debug extra info
//	>=4	massive info
//  >=5 +decode check
//
#define DEBUG_LEVEL   2

#include    "DebugDef.h"
//------------------------------------------------------//

//whether build hit which do not have all tdc and adc hit info
#define BUILD_PARTIAL_HIT 0

//save flat info
#define SAVEFLAT 1

//whether delete hit, whose ypos is 2 times away from the edge of scint bar
#define CUT_ON_YPOS 1

//------------------------------------------------------//

class SBSTimingHodoscope : public THaSubDetector {

public:

	SBSTimingHodoscope( const char* name, const char* description, 
		THaDetectorBase* parent);

	SBSTimingHodoscope( );

	virtual ~SBSTimingHodoscope();
	virtual Int_t        InitOutput( THaOutput* output );

	virtual Int_t      Decode( const THaEvData& );
	virtual EStatus    Init( const TDatime& run_time); 
	virtual Int_t      CoarseProcess( TClonesArray& tracks );
	virtual Int_t      FineProcess( TClonesArray& tracks );
	//Int_t              ConstructTracks( TClonesArray * tracks = NULL ) {return 0;} 
	Bool_t             IsTracking() { return kFALSE; }
	virtual Bool_t     IsPid()      { return kFALSE; }

	//function for jump through lines starting with #
	char* ReadNumberSignStartComment( FILE* fp, char *buf, const int len );

	// alternative reconstructions
	Int_t              BuildAllBars( TClonesArray& trks ); // original, partial hits etc.
	// newer:
	Int_t              BuildCompleteBars( TClonesArray& trks ); // only complete hits
	//Int_t              CombineHits( TClonesArray& trks );       // build plane clusters

	// Public Get and Set functions for Hit Information
	Int_t GetNBars() const { return fBars->GetLast()+1; }

	// Per Event
	Int_t GetNHits()     const { return fHits->GetLast()+1; }
	Int_t GetNRefHits()  const { return fRefHits->GetLast()+1;}
	Int_t GetNLtHits()   const { return fLtHits->GetLast()+1; }
	Int_t GetNRtHits()   const { return fRtHits->GetLast()+1; } 
	Int_t GetNLaHits()   const { return fLaHits->GetLast()+1; }
	Int_t GetNRaHits()   const { return fRaHits->GetLast()+1; }
	Int_t GetNPartHits() const { return fPartHits->GetLast()+1;}
	//Int_t GetNCombinedHits() const { return fCombHits->GetLast()+1; }

	TClonesArray* GetBars() const {return fBars;}
	THaScintBar* GetBar(Int_t i) const
	{return (THaScintBar*)fBars->At(i);}

	TClonesArray* GetHits() const {return fHits;}
	THaScintHit* GetHit(Int_t i) const
	{return (THaScintHit*)fHits->At(i);}

	//TClonesArray* GetCombinedHits() const {return fCombHits;}
	//THaMultiHit* GetCombHit(Int_t i) const
	//  {return (THaMultiHit*)fCombHits->At(i);}

	TClonesArray* GetRefHits() const {return fRefHits;}
	const THaTdcHit* GetRefHit(Int_t i) const
	{return (THaTdcHit*)fRefHits->At(i);}

	TClonesArray* GetLtHits() const {return fLtHits;}
	const THaTdcHit* GetLtHit(Int_t i) const
	{return (THaTdcHit*)fLtHits->At(i);}

	TClonesArray* GetRtHits() const {return fRtHits;}
	const THaTdcHit* GetRtHit(Int_t i) const
	{return (THaTdcHit*)fRtHits->At(i);}

	TClonesArray* GetLaHits() const {return fLaHits;}
	const THaAdcHit* GetLaHit(Int_t i) const
	{return (THaAdcHit*)fLaHits->At(i);}

	TClonesArray* GetRaHits() const {return fRaHits;}
	const THaAdcHit* GetRaHit(Int_t i) const
	{return (THaAdcHit*)fRaHits->At(i);}  

	// return matching Tdc from bar ptr on side, n'th hit
	const THaTdcHit* GetBarHitT(const char side, const THaScintBar *const ptr,
		const int n=0) const;

	// return matching Adc from bar ptr on side, n'th hit
	const THaAdcHit* GetBarHitA(const char side, const THaScintBar *const ptr,
		const int n=0) const;

	Int_t          GetNRefCh()    const { return fRefCh->GetLast()+1; }
	TClonesArray*  GetRefCh()     const { return fRefCh; }
	const THaScintPMT*   GetRefCh(Int_t i) const
	{ return (THaScintPMT*)fRefCh->At(i);}

	Int_t AreRefChOkay() const { return ( fRefOkay ? 1 : 0 ); } 

	TClonesArray* GetPartHits() const {return fPartHits;}
	const THaPartialHit* GetPartHit(Int_t i) const
	{return (THaPartialHit*)fPartHits->At(i);}  

protected:
	
	Int_t GetParameter( FILE* file, const TString tag, Double_t* value  );
	Int_t GetTable( FILE* file, const TString tag, Double_t* value,
		const Int_t maxval, int* first, int* last );

	TClonesArray* fBars;  
	Int_t fNBars;

	TClonesArray* fHits;     // Hits in paddles
	//TClonesArray* fCombHits; //   combined hit information for plane

	TClonesArray* fRefHits;  

	TClonesArray* fLaHits;
	TClonesArray* fRaHits;
	TClonesArray* fLtHits;
	TClonesArray* fRtHits;

	TClonesArray* fPartHits;

	// same for reference channels 
	TClonesArray*  fRefCh;        // reference channels for pipeline tdc
	Int_t fNRefHits;  // Number of refch that were hit
	Bool_t fRefOkay; 

	//Double_t fThreshold; //   energy threshhold for a hit to be used

	// Per-event data
	Int_t      fLTNhit;    // Number of Left paddles TDC times */
	Int_t      fRTNhit;    // Number of Right paddles TDC times  */
	Int_t      fLANhit;    // Number of Left paddles ADC amplitudes */
	Int_t      fRANhit;    // Number of Right paddles ADC amplitudes */
	Int_t      fNPaddlesHit;
	// per-event: index in clonesarray of first hit on bar
	Int_t     *fLtIndex;  //![fNBars]   //   for left-tdc
	Int_t     *fRtIndex;  //![fNBars]   //       right-tdc
	Int_t     *fLaIndex;  //![fNBars]   //       left-adc
	Int_t     *fRaIndex;  //![fNBars]   //       right-adc

	Double_t* fLE;        // [fNBars]    
	Double_t* fRE;        // [fNBars]
	Double_t* fLrawA;     // [fNBars]    
	Double_t* fRrawA;     // [fNBars]
	Double_t* fLpedcA;    // [fNBars]    
	Double_t* fRpedcA;    // [fNBars]

	Double_t* fLT;        // [fNBars]
	Double_t* fRT;        // [fNBars]
	Int_t* fLTcounter;    // [fNBars]
	Int_t* fRTcounter;    // [fNBars]

	Int_t* hitcounter;    // [fNBars] 
	Double_t* Energy;     // [fNBars] 		
	Double_t* TDIFF;      // [fNBars] 
	Double_t* TOF;        // [fNBars] 
	Double_t* T_tot;      // [fNBars] 
	Double_t* Yt_pos;     // [fNBars] 
	Double_t* Ya_pos;     // [fNBars] 	
	Double_t* Y_pred;     // [fNBars]
	Double_t* Y_dev;      // [fNBars]
	Double_t* fAngle;     // [fNBars]

	TVector3  fXax;

    //event statistic
    Int_t   fEventCount;            //how many event processed

    Int_t   fErrorReferenceChCount; //how many event got no reference channel
    Double_t fErrorReferenceChRateWarningThreshold; 
    //Threshold of ave. error reference per event that will pop up warnings
    Bool_t  fTooManyErrRefCh;       //flag whether there are too much error reference


	void           ClearEvent();
	void           DeleteArrays();
	virtual Int_t  ReadDatabase( const TDatime& date );
	virtual Int_t  DefineVariables( EMode mode = kDefine );
	virtual  Double_t TimeWalkCorrection(
        THaScintPMT* pmt,
        Double_t ADC,
        Double_t time);
	enum ESide { kLeft = 0, kRight = 1 };

	SBSTimingHodoscope& operator=( const SBSTimingHodoscope& ) {return *this; }

public:

	ClassDef(SBSTimingHodoscope,4) // Describes scintillator plane with F1TDC as a subdetector
};

#endif
