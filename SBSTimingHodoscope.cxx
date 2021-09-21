//////////////////////////////////////////////////////////////////////////
//
// SBSTimingHodoscope class implementation
// rmontgom@jlab.org july 2021 drafting
//////////////////////////////////////////////////////////////////////////

#include "SBSTimingHodoscope.h"

ClassImp(SBSTimingHodoscope);

/*
 * SBSTimingHodoscope constructor.
 *
 * Use a TDC with trailing edge info, default is no ADC, but available for
 * commissioning only
 */
SBSTimingHodoscope::SBSTimingHodoscope( const char* name, const char* description,
    THaApparatus* apparatus ) : SBSGenericDetector(name,description,apparatus)
{
  SetModeTDC(SBSModeTDC::kTDC); //  A TDC with leading & trailing edge info
  SetModeADC(SBSModeADC::kNone); // Default is No ADC, but can be re-enabled later
  // SBSGenericDetector::SetStoreRawHits(true); //raw doesn't atm include all tdc hits - should request this
}

///////////////////////////////////////////////////////////////////////////////
/// Read SBSTimingHodoscope Database
Int_t SBSTimingHodoscope::ReadDatabase( const TDatime& date )
{
  // We can use this name here for logs
  //static const char* const here = "ReadDatabase()";

  // Call parent class's ReadDatabase first
  Int_t err = SBSGenericDetector::ReadDatabase(date);
  if(err)
    return err;

  // If we want to add any new variables, uncomment the following and add
  // the new variables we want to read from the database
  FILE* file = OpenFile( date );
  if( !file ) return kFileError;

  std::cout << "******** Detector " << GetName() << " ReadDatabase ********" << std::endl;

  // get time walk or other parameters from database file
  // Read mapping/geometry/configuration parameters
  fChanMapStart = 0;
  fTDCBarOffset = 0;
  Int_t tdcbaroff = 0;
  fADCBarOffset = 0;
  Int_t adcbaroff = 0;
  std::vector<Float_t> timewalkpar0;
  std::vector<Float_t> timewalkpar1;
  DBRequest config_request[] = {
    { "timewalk0map",    &timewalkpar0,      kFloatV, 0, true }, //parameter for time walk correction
    { "timewalk1map",    &timewalkpar1,      kFloatV, 0, true }, //parameter for time walk correction
    { "tdcbaroffset",    &tdcbaroff,         kInt,    0, true }, //to allow for cycling through sections
    { "adcbaroffset",    &adcbaroff,         kInt,    0, true }, //to allow for cycling through sections
    { 0 } ///< Request must end in a NULL
  };
  err = LoadDB( file, date, config_request, fPrefix );
  if(err)
    return err;

  // std::cout << "fNelem " << fNelem << std::endl;
  // std::cout << "timewalkpar0.size() " << timewalkpar0.size() << std::endl;
  // std::cout << "timewalkpar1.size() " << timewalkpar1.size() << std::endl;

  // assign the bar offsets
  fTDCBarOffset = tdcbaroff;
  fADCBarOffset = adcbaroff;

  if( WithTDC() ){
    if(timewalkpar0.size()!=fNelem || timewalkpar1.size()!=fNelem){
      Error( Here("ReadDatabase"),
	     "N elements for hodoscope time walk maps != fNelem");
      return kInitError;
    }//error on time walk map size
    else{
      // put the time walk parameters into maps
      fTimeWalkPar0.clear();
      fTimeWalkPar1.clear();
      // std::cout << "fNelem " << fNelem << std::endl;
      // std::cout << "fNrows " << fNrows << std::endl;
      // std::cout << "fNRefElem " << fNRefElem << std::endl;
      fTimeWalkPar0.resize(fNrows);
      fTimeWalkPar1.resize(fNrows);
      int rr=0;
      int cc=0;
      int ll=0;
      int k=0;
      for(int r = 0; r < fNrows; r++) {
	rr = r+fChanMapStart;
	// std::cout << "On row " << r << " fNcols " << fNcols[r] << std::endl;
	fTimeWalkPar0[r].resize(fNcols[r]);
	fTimeWalkPar1[r].resize(fNcols[r]);
	for(int c = 0; c < fNcols[r]; c++) {
	  cc = c+fChanMapStart;
	  for(int l = 0; l < fNlayers; l++, k++) {
	    // std::cout << "On col " << c << " fNlayers " << fNlayers << std::endl;
	    // std::cout << "k " << k << std::endl;
	    fTimeWalkPar0[r][c].resize(fNlayers);
	    fTimeWalkPar1[r][c].resize(fNlayers);
	    ll = l+fChanMapStart;
	    fTimeWalkPar0[rr][cc][ll] = timewalkpar0[k];
	    fTimeWalkPar1[rr][cc][ll] = timewalkpar1[k];
	    // std::cout << "timewalkpar0[k] " << timewalkpar0[k] << " timewalkpar1[k] " << timewalkpar1[k] << std::endl;
	  }//lay
	}//col
      }//row
    }// if time walk map size>0
  }// if tdc then get time walk into a grid if needed

  // Make sure to call parent class so that the generic variables can be read
  // return SBSGenericDetector::ReadDatabase(date);
  if(err)
    return err;
  
  // call the function to build the bars
  SBSTimingHodoscope::ConstructHodoscope();
  
  // All is well that ends well
  return kOK;
}

//_____________________________________________________________________________
Int_t SBSTimingHodoscope::DefineVariables( EMode mode )
{
  // Initialize global variables
  Int_t err = SBSGenericDetector::DefineVariables(mode);
  if(err) {
    return err;
  }

  if(WithTDC()){
    RVarDef vars[] = {
      { "tdcbaroff",    "Starting bar for TDC readout",       "fTDCBarOffset" },
      { "tdcbarid",     "TDC Hit Bar ID",                     "fGoodBarIDsTDC"},
      { "barmeantime",  "Bar Mean Time [ns]",                 "fGoodBarTDCmean"},
      { "bartimediff",  "Bar Time Diff [ns]",                 "fGoodBarTDCdiff"},
      { "bartimehitpos","Bar Time Hit pos from L [m]",        "fGoodBarTDCpos"},
      { "barvpos",      "Bar vertical position [m]",          "fGoodBarTDCvpos"},
      { "L.le",         "Left pmt time LE [ns]",              "fGoodBarTDCLle"},
      { "L.leW",        "Left pmt time LE walk corr [ns]",    "fGoodBarTDCLleW"},
      { "L.te",         "Left pmt time TE [ns]",              "fGoodBarTDCLte"},
      { "L.teW",        "Left pmt time TE walk corr [ns]",    "fGoodBarTDCLteW"},
      { "L.tot",        "Left pmt tot [ns]",                  "fGoodBarTDCLtot"},
      { "L.totW",       "Left pmt tot walk corr [ns]",        "fGoodBarTDCLtotW"},
      { "R.le",         "Right pmt time LE [ns]",             "fGoodBarTDCRle"},
      { "R.leW",        "Right pmt time LE walk corr [ns]",   "fGoodBarTDCRleW"},
      { "R.te",         "Right pmt time TE [ns]",             "fGoodBarTDCRte"},
      { "R.teW",        "Right pmt time TE walk corr [ns[",   "fGoodBarTDCRteW"},
      { "R.tot",        "Right pmt tot [ns[",                 "fGoodBarTDCRtot"},
      { "R.totW",       "Right pmt tot walk corr [ns]",       "fGoodBarTDCRtotW"},
      { 0 }
    };
    err = DefineVarsFromList( vars, mode );
    if(err)
      return err;
  }// if we're in a tdc event
  if(WithADC()){
    RVarDef vars[] = {
      { "adcbaroff",    "Starting bar for ADC readout",       "fADCBarOffset" },
      { "adcbarid",     "ADC Hit Bar ID",                     "fGoodBarIDsADC"},
      { "adcmean",      "ADC Hit Bar Mean [bins]",            "fGoodBarADCmean"},
      { "L.a",          "Left ADC [bins]",                    "fGoodBarADCLa"},
      { "L.ap",         "Left ADC ped corr [bins]",           "fGoodBarADCLap"},
      { "L.ac",         "Left ADC ped corr [GeV]",            "fGoodBarADCLac"},
      { "R.a",          "Right ADC [bins]",                   "fGoodBarADCRa"},
      { "R.ap",         "Right ADC ped corr [bins]",          "fGoodBarADCRap"},
      { "R.ac",         "Right ADC ped corr [GeV]",           "fGoodBarADCRac"},
      { 0 }
    };
    err = DefineVarsFromList( vars, mode );
    if(err)
      return err;
  }// adc mode

  // Finally go back
  return err;
}

/*
 * FindGoodHit()
 */
Int_t SBSTimingHodoscope::FindGoodHit(SBSElement *blk)
{
  Int_t GoodHit=0;  

  // if (blk->TDC()&& blk->HasData()) {
  //   blk->TDC()->SetGoodHit(0);
  //   GoodHit=1;
  // }
  return GoodHit;
}

Int_t SBSTimingHodoscope::CoarseProcess( TClonesArray& tracks )
{
  if(fCoarseProcessed)
    return 0;

  // Call the parent class so that it can prepare the data structure on the
  // event it just read from file
  // SBSGenericDetector::CoarseProcess(tracks);
  // Call parent class, and exit if an error is encountered
  Int_t err = SBSGenericDetector::CoarseProcess(tracks);
  if(err)
    return err;

  // the parent class ie gen det finds the good hits we need already
  // note it does this by finding one good tdc hit based on a good timing window cut
  // (as well as demanding LE, TE in certain order)
  // it finds the tdc hit closest to the good timing window cut
  // so we need to be careful about setting this number in the db file
  // and figuring out the value in the expt
  
  // now loop through bars to find good hits in bars
  // should we move this code into findgoodhit?

  fGoodBarIDsTDC.clear();
  fGoodBarTDCmean.clear();
  fGoodBarTDCdiff.clear();
  fGoodBarTDCpos.clear();
  fGoodBarTDCLle.clear();
  fGoodBarTDCLleW.clear();
  fGoodBarTDCLte.clear();
  fGoodBarTDCLteW.clear();
  fGoodBarTDCLtot.clear();
  fGoodBarTDCLtotW.clear();
  fGoodBarTDCRle.clear();
  fGoodBarTDCRleW.clear();
  fGoodBarTDCRte.clear();
  fGoodBarTDCRteW.clear();
  fGoodBarTDCRtot.clear();
  fGoodBarTDCRtotW.clear();
  fGoodBarIDsADC.clear();
  fGoodBarADCmean.clear();
  fGoodBarADCLa.clear();
  fGoodBarADCLap.clear();
  fGoodBarADCLac.clear();
  fGoodBarADCRa.clear();
  fGoodBarADCRap.clear();
  fGoodBarADCRac.clear();

  Int_t NBars = (Int_t)fBars.size();
  // std::cout << "fTDCBarOffset " << fTDCBarOffset << std::endl;
  // std::cout << "fADCBarOffset " << fADCBarOffset << std::endl;
  // std::cout << "NBars " << NBars << " fElements.size()/2 " << fElements.size()/2 << std::endl;
  if(NBars!=(fElements.size()/2)){
    Error( Here("CoarseProcess"),
	   "hodoscope #bars length not of correct size ie !=#elements/2");
    return kInitError;
  }
  for(Int_t BarInc=0; BarInc<NBars; BarInc++){
    SBSTimingHodoscopeBar *bar = fBars[BarInc];
    SBSTimingHodoscopePMT *pmtL = bar->GetLPMT();
    SBSTimingHodoscopePMT *pmtR = bar->GetRPMT();
    SBSElement *elL = pmtL->GetPMTElement();
    SBSElement *elR = pmtR->GetPMTElement();
    
    if(WithTDC()){
      if(elL->TDC()->HasData() && elR->TDC()->HasData()){
	Int_t bar = BarInc;
	// don't need to add offset to tdc since all readout simultaneously
	fGoodBarIDsTDC.push_back(bar);
	
	// left hit
	const SBSData::TDCHit &hitL = elL->TDC()->GetGoodHit();
	fGoodBarTDCLle.push_back(hitL.le.val);//.raw is tdc bin, val is corrected using offset and ns/bin
	// call the time walk function here once implemented
	fGoodBarTDCLleW.push_back(hitL.le.val);//will replace this with time walk corrected
	fGoodBarTDCLte.push_back(hitL.te.val);
	// call the time walk function here once implemented
	fGoodBarTDCLteW.push_back(hitL.te.val);//will replace this with time walk corrected
	fGoodBarTDCLtot.push_back(hitL.te.val-hitL.le.val);
	fGoodBarTDCLtotW.push_back(hitL.te.val-hitL.le.val);//replace this with walk corrected times
	
	// right hit
	const SBSData::TDCHit &hitR = elR->TDC()->GetGoodHit();
	fGoodBarTDCRle.push_back(hitR.le.val);//.raw is tdc bin, val is corrected using offset and ns/bin
	// call the time walk function here once implemented
	fGoodBarTDCRleW.push_back(hitR.le.val);//will replace this with time walk corrected
	fGoodBarTDCRte.push_back(hitR.te.val);
	// call the time walk function here once implemented
	fGoodBarTDCRteW.push_back(hitR.te.val);//will replace this with time walk corrected
	fGoodBarTDCRtot.push_back(hitR.te.val-hitR.le.val);
	fGoodBarTDCRtotW.push_back(hitR.te.val-hitR.le.val);//replace this with walk corrected times
	
	// bar properties
	Float_t barmeantime = (hitL.le.val + hitR.le.val)/2.0;
	fGoodBarTDCmean.push_back(barmeantime);
	Float_t bartimediff = (hitL.le.val - hitR.le.val);
	fGoodBarTDCdiff.push_back(bartimediff);
	// convert to position? effective velocity times time? should we divide by 2? yes
	Float_t HorizPos = 0.5 * (bartimediff*1.0e-9) * vScint; // position from L based on timediff and in m
	fGoodBarTDCpos.push_back(HorizPos);
	fGoodBarTDCpos.push_back(elR->GetY());
	
      }// tdc hit on both pmts
    }// with tdc
    // adc events
    else if(WithADC()){
      if(elL->ADC()->HasData() && elR->ADC()->HasData()){
	Int_t bar = fADCBarOffset + BarInc;
	if(bar>89){
	  Error( Here("CoarseProcess"),
		 "hodoscope bar id after adding adc offset id from db >89");
	  return kInitError;
	}
	fGoodBarIDsADC.push_back(bar);

	// left hit
	Float_t pedL=elL->ADC()->GetPed()*1.0;
	const SBSData::PulseADCData &hitL = elL->ADC()->GetGoodHit(); // this is assuming simple adc type should add a check
	fGoodBarADCLa.push_back(hitL.integral.raw);
	fGoodBarADCLap.push_back(hitL.integral.raw-pedL);
	fGoodBarADCLac.push_back(hitL.integral.val);
	
	// right hit
	Float_t pedR=elR->ADC()->GetPed()*1.0;
	const SBSData::PulseADCData &hitR = elR->ADC()->GetGoodHit();
	fGoodBarADCRa.push_back(hitR.integral.raw);
	fGoodBarADCRap.push_back(hitR.integral.raw-pedR);
	fGoodBarADCRac.push_back(hitR.integral.val);
	
	// bar properties
	Float_t barmeanadc = (hitL.integral.raw + hitR.integral.raw) / 2.0;
	// at moment dont use ped corrected for debug
	fGoodBarADCmean.push_back(barmeanadc);
      }// adc hit on both pmts
    }// with adc
  }// bar loop


  fCoarseProcessed = 1;
  return 0;
}

Int_t SBSTimingHodoscope::FineProcess( TClonesArray& tracks )
{

  if(fFineProcessed)
    return 0;

  // Do more detailed processing here.  Parent class does nothing, so no need
  // to call it.
  // We can prepare more detailed output if we want.

  fFineProcessed = 1;
  return 0;
}
/*
 * ConstructHodoscope()
 * called in read database 
 */
Int_t SBSTimingHodoscope::ConstructHodoscope()
{
  Int_t nElements = fElements.size();
  // std::cout << "n elements " << nElements << std::endl;
  if( nElements%2!=0 ) {
    Error( Here("ConstructHodoscope"),
	   "N elements for hodoscope is not even, need an even number for a 2 sided detector analysis.");
    return kInitError;
  }
  // make the pmt objects
  fPMTMapL.clear();
  fPMTMapR.clear();
  if( fNrows!=2 ) {
    Error( Here("ConstructHodoscope"),
	   "fNrows for hodoscope is not 2, which we need for left and right.");
    return kInitError;
  }
  int p=0;
  for(int r = 0; r < fNrows; r++) {
    for(int c = 0; c < fNcols[r]; c++) {
      for(int l = 0; l < fNlayers; l++, p++) {
	// get the element
	SBSElement *blk2 = fElementGrid[r][c][l];
	Int_t column = blk2->GetCol();
	Int_t row = blk2->GetRow();
	// std::cout << "row " << r << " col " << c << " l " << l  << " p " << p << std::endl;
	// std::cout << "element column " << column << " and row " << row << std::endl;
	if(row==0){//left// could also use r index
	  if( WithTDC() ){
	    SBSTimingHodoscopePMT *pmtL = new SBSTimingHodoscopePMT(fElements.at(p),
								    fTimeWalkPar0[r][c][l],
								    fTimeWalkPar0[r][c][l],
								    column, row, p);
	    // column = bar and row = side and p = ID in felements in case needed for debug
	    fPMTMapL.push_back(pmtL);
	  }// if tdc
	  else if( WithADC() ){ // if we have adc we use adc db file and therefore need a temp val for timewalk
	    SBSTimingHodoscopePMT *pmtL = new SBSTimingHodoscopePMT(fElements.at(p),
								    -9999,
								    -9999,
								    column, row, p);
	    fPMTMapL.push_back(pmtL);
	  }// if adc only
	}//left
	if(row==1){//right
	  if( WithTDC() ){
	    SBSTimingHodoscopePMT *pmtR = new SBSTimingHodoscopePMT(fElements.at(p),
								    fTimeWalkPar0[r][c][l],
								    fTimeWalkPar0[r][c][l],
								    column, row, p);
	    // column = bar and row = side and p = ID in felements in case needed
	    fPMTMapR.push_back(pmtR);
	  }// if tdc
	  else if( WithADC() ){ // if we have adc we use adc db file and therefore need a temp val for timewalk
	    SBSTimingHodoscopePMT *pmtR = new SBSTimingHodoscopePMT(fElements.at(p),
								    -9999,
								    -9999,
								    column, row, p);
	    fPMTMapR.push_back(pmtR);
	  }// if adc only
	}//right
      }//lay
    }//col
  }//row ie side of hodo
  // std::cout << "n elements in left pmt array " << fPMTMapL.size() << std::endl;
  // std::cout << "n elements in right pmt array " << fPMTMapR.size() << std::endl;
  // std::cout << "n elements " << nElements << std::endl;
  
  // now we have the arrays of left and right pmts, we need to make the bars
  if( fPMTMapL.size()!=fPMTMapR.size() || fPMTMapL.size()!=(nElements/2) || fPMTMapR.size()!=(nElements/2) ) {
    Error( Here("ConstructHodoscope"),
	   "PMT arrays for constructing hodoscope bars not of correct length");
    return kInitError;
  }
  fBars.clear();
  for(Int_t BarInc=0; BarInc<(Int_t)(nElements/2); BarInc++){
    // bar constructor is barid, pmt left, pmt right, bar offset mostly for adc sections
    if( WithTDC() ){
      // std::cout << fPMTMapL.at(BarInc)->GetTdcFlag() << std::endl;
      SBSTimingHodoscopeBar *bar = new SBSTimingHodoscopeBar(BarInc, fPMTMapL.at(BarInc),
							     fPMTMapR.at(BarInc), fTDCBarOffset);
      fBars.push_back(bar);
    }//tdc
    else if( WithADC() ){
      SBSTimingHodoscopeBar *bar = new SBSTimingHodoscopeBar(BarInc, fPMTMapL.at(BarInc),
							     fPMTMapR.at(BarInc), fADCBarOffset);
      fBars.push_back(bar);
    }//adc
  }//bar loop
  // std::cout << "We have filled " << fBars.size() << " bars" << std::endl;

}// construct hodo
/*
 * TimeWalk()
 */
 Float_t TimeWalk(Float_t time){
   return 0.0;
 }
/*
 * ClearEvent()
 * called at the end of every event
 */
void SBSTimingHodoscope::ClearEvent()
{
  // If we defined any new variables that we need to clear prior to the next event
  // clear them here:
  fGoodBarIDsTDC.clear();
  fGoodBarTDCmean.clear();
  fGoodBarTDCdiff.clear();
  fGoodBarTDCpos.clear();
  fGoodBarTDCvpos.clear();
  fGoodBarTDCLteW.clear();
  fGoodBarTDCLtot.clear();
  fGoodBarTDCLtotW.clear();
  fGoodBarTDCRle.clear();
  fGoodBarTDCRleW.clear();
  fGoodBarTDCRte.clear();
  fGoodBarTDCRteW.clear();
  fGoodBarTDCRtot.clear();
  fGoodBarTDCRtotW.clear();
  fGoodBarIDsADC.clear();
  fGoodBarADCmean.clear();
  fGoodBarADCLa.clear();
  fGoodBarADCLap.clear();
  fGoodBarADCLac.clear();
  fGoodBarADCRa.clear();
  fGoodBarADCRap.clear();
  fGoodBarADCRac.clear();

  // Make sure to call parent class's ClearEvent() also!
  SBSGenericDetector::ClearEvent();
}


/*
 * Generic SBSTimingHodoscope destructor
 */
SBSTimingHodoscope::~SBSTimingHodoscope()
{
  // Delete any new objects/instances created here
  ClearEvent();
}
