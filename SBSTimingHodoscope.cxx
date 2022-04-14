//////////////////////////////////////////////////////////////////////////
//
// SBSTimingHodoscope class implementation
//
//////////////////////////////////////////////////////////////////////////

#include "SBSTimingHodoscope.h"
#include "Helper.h"

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
  // SBSGenericDetector::SetStoreRawHits(true);
  fDataOutputLevel = 0;//default
  
  //default values for bar quality
  fHorizPosBarCut = 0.3;//m most sensible default value
  fTimeRef = 80.0;// ?
  fTimeBarCut = 10.0;//
  
  //default values for clustering parameters:
  fClusMaxSize = 5;
  fMaxYposDiffCluster = 0.15;//m
  fMaxTimeDiffCluster = 10.0;//ns
  
  fTrackMatchCutX = 0.05;
  fTrackMatchCutY = 0.15;

  fTDCBarOffset = 0;
  fADCBarOffset = 32;

  fTDCWinMin = -20.;
  fTDCWinMax = 20.;
  fTotMin = 7.;
  fTotMax= 30.;

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
  fIsInit = false;

  // If we want to add any new variables, uncomment the following and add
  // the new variables we want to read from the database
  FILE* file = OpenFile( date );
  if( !file ) return kFileError;

  std::cout << "******** Detector " << GetName() << " ReadDatabase ********" << std::endl;

  // get time walk or other parameters from database file
  // Read mapping/geometry/configuration parameters
  std::vector<Double_t> timewalkpar0;
  std::vector<Double_t> timewalkpar1;
  DBRequest config_request[] = {
    { "timewalk0map",    &timewalkpar0,      kDoubleV, 0, false }, //parameter for time walk correction
    { "timewalk1map",    &timewalkpar1,      kDoubleV, 0, false }, //parameter for time walk correction
    { "tdcbaroffset",    &fTDCBarOffset,         kInt,    0, true }, //to allow for cycling through sections
    { "adcbaroffset",    &fADCBarOffset,         kInt,    0, true }, //to allow for cycling through sections
    { "tdcwindowmin", &fTDCWinMin, kDouble, 0, true }, //parameter for tdc window min
    { "tdcwindowmax", &fTDCWinMax, kDouble, 0, true }, //parameter for tdc window max
    { "tdctotmin", &fTotMin, kDouble, 0, true }, //parameter for tdc window max
    { "tdctotmax", &fTotMax, kDouble, 0, true }, //parameter for tdc window max
    { 0 } ///< Request must end in a NULL
  };
  err = LoadDB( file, date, config_request, fPrefix );
  if(err) {
    fclose(file);
    return err;
  }
  
  DBRequest barquality_params[] = {
    { "horizposbarcut", &fHorizPosBarCut, kDouble, 0, true }, //parameter for bar horizontal position selection
    // what is below shall basically be redundant with fTDCWinMin and TDCWinMax
    { "timeref",      &fTimeRef,      kDouble, 0, true }, //parameter for time reference
    { "timebarcut", &fTimeBarCut, kDouble, 0, true }, //parameter for time bar selection
    { 0 } ///< Request must end in a NULL
  };

  err = LoadDB( file, date, barquality_params, fPrefix );
  if(err) {
    fclose(file);
    return err;
  }
  
  DBRequest clustering_params[] = {
    { "maxclussize",      &fClusMaxSize,        kInt, 0, true }, //parameter for max cluster size
    { "maxyposdiff_clus", &fMaxYposDiffCluster, kDouble, 0, true }, //parameter for max position difference for bar to be incorporated into a cluster
    { "maxtimediff_clus", &fMaxTimeDiffCluster, kDouble, 0, true }, //parameter for max time difference for bar to be incorporated into a cluster
    { 0 } ///< Request must end in a NULL
  };

  err = LoadDB( file, date, clustering_params, fPrefix );
  if(err) {
    fclose(file);
    return err;
  }
  
  DBRequest trackmatch_params[] = {
    { "vscint",         &fvScint, kDoubleV, 0, 1, 1},
    { "tdiffoffset",    &ftDiff0, kDoubleV, 0, 1, 1},
    { "trackmatchcutX",  &fTrackMatchCutX, kDouble, 0, 1, 1},
    { "trackmatchcutY",  &fTrackMatchCutY, kDouble, 0, 1, 1},
    { 0 }
  };
  err = LoadDB( file, date, trackmatch_params, fPrefix );
  if(err) {
    fclose(file);
    return err;
  }

  // If problem with vscint or tdiff0 db entries, set defaults
  if( fvScint.size() != (UInt_t)fNelem/2 ) {
    std::cout << " vScint vector too small " << fvScint.size() << " # of bars =" << fNelem/2 << std::endl;
    fvScint.resize((UInt_t)fNelem/2);
    std::fill(fvScint.begin(), fvScint.end(),  0.178 ); // m/ns
  }
  if( ftDiff0.size() != (UInt_t)fNelem/2 ) {
    std::cout << " tDiff0 vector too small " << ftDiff0.size() << " # of bars =" << fNelem/2 << std::endl;
    ftDiff0.resize((UInt_t)fNelem/2);
    std::fill(ftDiff0.begin(), ftDiff0.end(),  -1.35 ); // ns
  }

  std::vector<Double_t> ypos;//position of element
  std::vector<DBRequest> vr;
  vr.push_back({ "ypos", &ypos,    kDoubleV, 0, 1 });
  vr.push_back({0});
  err = LoadDB( file, date, vr.data(), fPrefix );
  if(err) {
    fclose(file);
    return err;
  }
  if (ypos.size()>0) {
    if (ypos.size() == (UInt_t)fNelem) {
      for (Int_t ne=0;ne<fNelem;ne++) {
	SBSElement* blk= fElements[ne];
	blk->SetY(ypos[ne]);
      }
    } else {
      std::cout << " ypos vector too small " << ypos.size() << " # of elements =" << fNelem << std::endl;
    }
  }
  
  fclose(file);

  if( WithTDC() || WithADC()){
    if(timewalkpar0.size()!=(UInt_t)fNelem || timewalkpar1.size()!=(UInt_t)fNelem){
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
      int k=0;
      for(int r = 0; r < fNrows; r++) {
	// std::cout << "On row " << r << " fNcols " << fNcols[r] << std::endl;
	fTimeWalkPar0[r].resize(fNcols[r]);
	fTimeWalkPar1[r].resize(fNcols[r]);
	for(int c = 0; c < fNcols[r]; c++) {
          fTimeWalkPar0[r][c].resize(fNlayers);
          fTimeWalkPar1[r][c].resize(fNlayers);
	  for(int l = 0; l < fNlayers; l++, k++) {
	    // std::cout << "On col " << c << " fNlayers " << fNlayers << std::endl;
	    // std::cout << "k " << k << std::endl;
	    fTimeWalkPar0[r][c][l] = timewalkpar0[k];
	    fTimeWalkPar1[r][c][l] = timewalkpar1[k];
	    // std::cout << "timewalkpar0[k] " << timewalkpar0[k] << " timewalkpar1[k] " << timewalkpar1[k] << std::endl;
	  }//lay
	}//col
      }//row
    }// if time walk map size>0
  }// if tdc then get time walk into a grid if needed

  // call the function to build the bars
  err = SBSTimingHodoscope::ConstructHodoscope();
  if(err)
    return err;

  // Make sure to call parent class so that the generic variables can be read
  // return SBSGenericDetector::ReadDatabase(date);
  
  // All is well that ends well
  fIsInit = true;
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

  if(WithTDC() ){
    RVarDef vars[] = {
      // { "tdcbaroff",    "Starting bar for TDC readout",       "fTDCBarOffset" },
      { "bar.ngoodbars",    "Number of good bars",             "GetGoodBarsSize()"},
      { "bar.tdc.id",    "TDC Hit Bar ID",                     "fGoodBarIDsTDC"},
      { "bar.tdc.meantime",  "Bar Mean Time [ns]",                 "fGoodBarTDCmean"},
      { "bar.tdc.timediff",  "Bar Time Diff [ns]",                 "fGoodBarTDCdiff"},
      { "bar.tdc.timehitpos","Bar Time Hit pos from L [m]",        "fGoodBarTDCpos"},
      { "bar.tdc.vpos",      "Bar vertical position [m]",          "fGoodBarTDCvpos"},
      { "bar.tdc.L.le",         "Left pmt time LE [ns]",              "fGoodBarTDCLle"},
      { "bar.tdc.L.leW",        "Left pmt time LE walk corr [ns]",    "fGoodBarTDCLleW"},
      { "bar.tdc.L.te",         "Left pmt time TE [ns]",              "fGoodBarTDCLte"},
      { "bar.tdc.L.teW",        "Left pmt time TE walk corr [ns]",    "fGoodBarTDCLteW"},
      { "bar.tdc.L.tot",        "Left pmt tot [ns]",                  "fGoodBarTDCLtot"},
      { "bar.tdc.L.totW",       "Left pmt tot walk corr [ns]",        "fGoodBarTDCLtotW"},
      { "bar.tdc.R.le",         "Right pmt time LE [ns]",             "fGoodBarTDCRle"},
      { "bar.tdc.R.leW",        "Right pmt time LE walk corr [ns]",   "fGoodBarTDCRleW"},
      { "bar.tdc.R.te",         "Right pmt time TE [ns]",             "fGoodBarTDCRte"},
      { "bar.tdc.R.teW",        "Right pmt time TE walk corr [ns[",   "fGoodBarTDCRteW"},
      { "bar.tdc.R.tot",        "Right pmt tot [ns[",                 "fGoodBarTDCRtot"},
      { "bar.tdc.R.totW",       "Right pmt tot walk corr [ns]",       "fGoodBarTDCRtotW"},
      { 0 }
    };
    err = DefineVarsFromList( vars, mode );
    if(err)
      return err;
  }// if we're in a tdc event
  if(WithADC() ){
    RVarDef vars[] = {
      // { "adcbaroff",    "Starting bar for ADC readout",       "fADCBarOffset" },
      { "bar.adc.id",     "ADC Hit Bar ID",                     "fGoodBarIDsADC"},
      { "bar.adc.mean",    "ADC Hit Bar Mean [bins]",            "fGoodBarADCmean"},
      { "bar.adc.L.a",          "Left ADC [bins]",                    "fGoodBarADCLa"},
      { "bar.adc.L.ap",         "Left ADC ped corr [bins]",           "fGoodBarADCLap"},
      { "bar.adc.L.ac",         "Left ADC ped corr [GeV]",            "fGoodBarADCLac"},
      { "bar.adc.R.a",          "Right ADC [bins]",                   "fGoodBarADCRa"},
      { "bar.adc.R.ap",         "Right ADC ped corr [bins]",          "fGoodBarADCRap"},
      { "bar.adc.R.ac",         "Right ADC ped corr [GeV]",           "fGoodBarADCRac"},
      { 0 }
    };
    err = DefineVarsFromList( vars, mode );
    if(err)
      return err;
  }// adc mode
  
  //if(fDataOutputLevel>1){
  RVarDef vars_clus[] = {
    { "allclus.size",  "cluster size",          "fOutClus.n"},
    { "allclus.id", "cluster max bar id",     "fOutClus.id" },
    { "allclus.xmean", "cluster mean X",        "fOutClus.x"},
    { "allclus.ymean", "cluster mean Y",        "fOutClus.y"},
    { "allclus.tmean", "cluster mean T",        "fOutClus.t"},
    { "allclus.totmean", "cluster mean ToT",    "fOutClus.tot"},
    { "allclus.tdiff", "cluster max bar tdiff", "fOutClus.tdiff"},
    { "allclus.itrack", "track index", "fOutClus.trackindex"},
    { 0 }
  };
  err = DefineVarsFromList( vars_clus, mode );
    //}
  
  //if(fDataOutputLevel>0){
  RVarDef vars_clushits[] = {
    { "clus.bar.tdc.id",       "main clus TDC Hit Bar ID",     "fMainClusBars.id"},
    { "clus.bar.tdc.meantime", "main clus Bar Mean Time [ns]", "fMainClusBars.t"},
    { "clus.bar.tdc.meantot",  "main clus Bar Mean ToT [ns]",  "fMainClusBars.tot"},
    { "clus.bar.tdc.timediff", "main clus Bar Time Diff [ns]", "fMainClusBars.tdiff"},
    { "clus.bar.tdc.timehitpos", "main clus Bar Time Hit pos from L [m]",  "fMainClusBars.y"},
    { "clus.bar.tdc.vpos",       "main clus Bar vertical position [m]",    "fMainClusBars.x"},
    { "clus.bar.tdc.itrack",  "main clus Bar track index", "fMainClusBars.trackindex" },
    { 0 }
  };

 
  err = DefineVarsFromList( vars_clushits, mode );
    //}
  
  RVarDef vars_mainclus[] = {
    { "nclus",   "number of clusters", "GetNClusters()"},
    { "clus.id",  "cluster max bar id",       "fMainClus.id"},
    { "clus.size", "cluster size",       "fMainClus.n"},
    { "clus.xmean", "cluster mean X",     "fMainClus.x"},
    { "clus.ymean",  "cluster mean Y",     "fMainClus.y"},
    { "clus.tmean",   "cluster mean T",     "fMainClus.t"},
    { "clus.totmean",  "cluster mean ToT",   "fMainClus.tot"},
    { "clus.tdiff", "cluster max bar tdiff", "fMainClus.tdiff"},
    { "clus.trackindex", "cluster track index", "fMainClus.trackindex"},
    { 0 }
  };
  err = DefineVarsFromList( vars_mainclus, mode );
  
  if(err)
    return err;
  
  // Finally go back
  return err;
}

/*
 * FindGoodHit()
 */
// Int_t SBSTimingHodoscope::FindGoodHit(SBSElement *blk)
// {
//   Int_t GoodHit=0;  

//   // if (blk->TDC()&& blk->HasData()) {
//   //   blk->TDC()->SetGoodHit(0);
//   //   GoodHit=1;
//   // }
//   return GoodHit;
// }

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
  fGoodBarTDCvpos.clear();
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
  //std::cout << "fTDCWinMin: " << fTDCWinMin << std::endl;
  //std::cout << "fTDCWinMax: " << fTDCWinMax << std::endl;
  if((UInt_t)NBars!=(fElements.size()/2)){
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
    Double_t Ltwalk0 = pmtL->GetTimeWalkPar0();
    Double_t Ltwalk1 = pmtL->GetTimeWalkPar1();
    Double_t Rtwalk0 = pmtR->GetTimeWalkPar0();
    Double_t Rtwalk1 = pmtR->GetTimeWalkPar1();

    if(WithTDC()){
      if(elL->TDC()->HasData() && elR->TDC()->HasData()){
	
	//get left and right LE times 
	const SBSData::TDCHit &hitL = elL->TDC()->GetGoodHit();
	const SBSData::TDCHit &hitR = elR->TDC()->GetGoodHit();
	Double_t LEl = hitL.le.val;
	Double_t LEr = hitR.le.val;

	if(fTDCWinMin < LEl && LEl < fTDCWinMax && fTDCWinMin < LEr && LEr < fTDCWinMax
	   && (hitL.te.val-hitL.le.val) > fTotMin && (hitR.te.val-hitR.le.val) > fTotMin 
	   && (hitL.te.val-hitL.le.val) < fTotMax && (hitR.te.val-hitR.le.val) < fTotMax ) {

	//Int_t bar = BarInc;// why redeclare the index?
	// don't need to add offset to tdc since all readout simultaneously
	
	// left hit
	// fGoodBarTDCLle.push_back(hitL.le.val);//commenting this here: this is done later after another quality check
	// fGoodBarTDCLle.push_back(hitL.le.raw);
	//.raw is tdc bin, val is corrected using offset and ns/bin
	Double_t LleW = SBSTimingHodoscope::TimeWalk(hitL.le.val,
				(hitL.te.val-hitL.le.val),
				Ltwalk0, Ltwalk1);
	Double_t LteW = SBSTimingHodoscope::TimeWalk(hitL.te.val,
				(hitL.te.val-hitL.le.val),
				Ltwalk0, Ltwalk1);
	
	// right hit
	// fGoodBarTDCRle.push_back(hitR.le.val);//commenting this here: this is done later after another quality check
	Double_t RleW = SBSTimingHodoscope::TimeWalk(hitR.le.val,
				(hitR.te.val-hitR.le.val),
				Rtwalk0, Rtwalk1);
	Double_t RteW = SBSTimingHodoscope::TimeWalk(hitR.te.val,
				(hitR.te.val-hitR.le.val),
				Rtwalk0, Rtwalk1);
	
	// bar properties
	Double_t barmeantime = (LleW + RleW)/2.0;
	Double_t bartimediff = (LleW - RleW);
	// convert to position? effective velocity times time? should we divide by 2? yes
	// Double_t HorizPos = 0.5 * (bartimediff*1.0e-9) * vScint; // position from L based on timediff and in m
	//Assuming bartimediff is in ns, then horizontal position is

	//AJRP: the -sign is added because the time difference as defined here has
	// a negative correlation with the y of transport coordinates
	// The offset aligns this quantity with the GEM track projection to the hodoscope
	// fvScint ~= 0.454c is the average effective propagation speed as
	// measured by the GEM-TH correlation.
	Double_t HorizPos = -0.5 * (bartimediff-ftDiff0.at(BarInc)) * fvScint.at(BarInc); // position from L based on timediff and in m. 
	
	// check basic quality before pushing
	//if(fabs(HorizPos)>fHorizPosBarCut)continue;
	//if(fabs(barmeantime-fTimeRef)>fTimeBarCut)continue;
	
	// Grouping all of this here makes it easier to apply additional cuts, if needed.
	fGoodBarIDsTDC.push_back(BarInc);
	
	fGoodBarTDCLle.push_back(hitL.le.val);
	fGoodBarTDCLleW.push_back(LleW);
	fGoodBarTDCLte.push_back(hitL.te.val);
	fGoodBarTDCLteW.push_back(LteW);
	fGoodBarTDCLtot.push_back(hitL.te.val-hitL.le.val);
	fGoodBarTDCLtotW.push_back(LteW-LleW);

	fGoodBarTDCRle.push_back(hitR.le.val);//.raw is tdc bin, val is corrected using offset and ns/bin
	fGoodBarTDCRleW.push_back(RleW);
	fGoodBarTDCRte.push_back(hitR.te.val);
	fGoodBarTDCRteW.push_back(RteW);
	fGoodBarTDCRtot.push_back(hitR.te.val-hitR.le.val);
	fGoodBarTDCRtotW.push_back(RteW-RleW);
	
	fGoodBarTDCmean.push_back(barmeantime);
	fGoodBarTDCdiff.push_back(bartimediff);
	
	fGoodBarTDCpos.push_back(HorizPos);
	fGoodBarTDCvpos.push_back(elR->GetY());
	
	//std::cout << ((hitL.te.val-hitL.le.val)+(hitR.te.val-hitR.le.val))/2. << " " 
	//	  << ((LteW-LleW)+(RteW-RleW))/2. << std::endl;
	
	bar->SetMeanTime(barmeantime);
	bar->SetMeanToT( ((hitL.te.val-hitL.le.val)+(hitR.te.val-hitR.le.val))/2. );
	bar->SetTimeDiff(bartimediff);
	bar->SetHitPos(HorizPos);
	bar->SetElementPos(elR->GetY());
	bar->SetLeftHit(hitL);
	bar->SetRightHit(hitR);
	}
      }// tdc hit on both pmts
    }// with tdc
    // adc events
    else if(WithADC()){
      if(elL->ADC()->HasData() && elR->ADC()->HasData()){
	Int_t bar = fADCBarOffset + BarInc;
	fGoodBarIDsADC.push_back(bar);

	// left hit
	Double_t pedL=elL->ADC()->GetPed()*1.0;
	const SBSData::PulseADCData &hitL = elL->ADC()->GetGoodHit(); // this is assuming simple adc type should add a check
	fGoodBarADCLa.push_back(hitL.integral.raw);
	fGoodBarADCLap.push_back(hitL.integral.raw-pedL);
	fGoodBarADCLac.push_back(hitL.integral.val);
	
	// right hit
	Double_t pedR=elR->ADC()->GetPed()*1.0;
	const SBSData::PulseADCData &hitR = elR->ADC()->GetGoodHit();
	fGoodBarADCRa.push_back(hitR.integral.raw);
	fGoodBarADCRap.push_back(hitR.integral.raw-pedR);
	fGoodBarADCRac.push_back(hitR.integral.val);
	
	// bar properties
	Double_t barmeanadc = (hitL.integral.raw + hitR.integral.raw) / 2.0;
	// at moment dont use ped corrected for debug
	fGoodBarADCmean.push_back(barmeanadc);
      }// adc hit on both pmts
    }// with adc
  }// bar loop

  DoClustering();
  
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
  // Clustering here? 
  // Wait, if I understand the code, 
  // the way the information is stored in the vectors is by increasing index always.
  /*int nclusters = */

  //Moved "DoClustering() call to CoarseProcess, so we can use it in track search constraint calculation
  //DoClustering();

  //fill output here:
  //if(fDataOutputLevel>1){
  for(int i = 0; i<GetNClusters(); i++){
    fOutClus.id.push_back(GetCluster(i)->GetMaxBarID());
    fOutClus.n.push_back(GetCluster(i)->GetSize());
    fOutClus.x.push_back(GetCluster(i)->GetXmean());
    fOutClus.y.push_back(GetCluster(i)->GetYmean());
    fOutClus.t.push_back(GetCluster(i)->GetTmean());
    fOutClus.tot.push_back(GetCluster(i)->GetToTmean());
    fOutClus.tdiff.push_back(GetCluster(i)->GetTdiff());
    fOutClus.trackindex.push_back(-1);
  }
    //}
  
  int clus_idx = -1;
  Int_t n_trk = tracks.GetLast()+1;
  for( Int_t t = 0; t < n_trk; t++ ) {
    auto* theTrack = static_cast<THaTrack*>( tracks.At(t) ); 
    clus_idx = MatchTrack(theTrack);
    //anything to select the best track???
    if( clus_idx >= 0 ){
      fOutClus.trackindex[clus_idx] = t;
    }
  }

  //clus_idx is the index in the cluster array of 

  //fill main cluster variables with clusters on good tracks:
  for( int i=0; i<GetNClusters(); i++ ){
    if( fOutClus.trackindex[i] >= 0 ){
      fMainClus.id.push_back(GetCluster(i)->GetMaxBarID());
      fMainClus.n.push_back(GetCluster(i)->GetSize());
      fMainClus.x.push_back(GetCluster(i)->GetXmean());
      fMainClus.y.push_back(GetCluster(i)->GetYmean());
      fMainClus.t.push_back(GetCluster(i)->GetTmean());
      fMainClus.tot.push_back(GetCluster(i)->GetToTmean());
      fMainClus.tdiff.push_back(GetCluster(i)->GetTdiff());
      fMainClus.trackindex.push_back( fOutClus.trackindex[i] );
      //AJRP: we should fill these variables regardless;
      // output is contolled by the odef file:
      //if(fDataOutputLevel>0){
      for(int j = 0; j<GetCluster(i)->GetSize(); j++){
	fMainClusBars.id.push_back(GetCluster(i)->GetElement(j)->GetBarNum());
	fMainClusBars.n.push_back(1);
	fMainClusBars.t.push_back(GetCluster(i)->GetElement(j)->GetMeanTime());
	fMainClusBars.tot.push_back(GetCluster(i)->GetElement(j)->GetMeanToT());
	fMainClusBars.tdiff.push_back(GetCluster(i)->GetElement(j)->GetTimeDiff());
	fMainClusBars.x.push_back(GetCluster(i)->GetElement(j)->GetElementPos());
	fMainClusBars.y.push_back(GetCluster(i)->GetElement(j)->GetHitPos());
	fMainClusBars.trackindex.push_back( fOutClus.trackindex[i] );
      }
    }
  }
  
  

  fFineProcessed = 1;
  return 0;
}

Int_t SBSTimingHodoscope::DoClustering()
{
  //int prev_baridx = -10;
  int halfclussize = fClusMaxSize/2;
  
  //vector for local maxima?
  std::vector<int> localmax_idx;
  localmax_idx.clear();
  
  for(int i = 0; i<(int)fGoodBarIDsTDC.size(); i++){
    int baridx = fGoodBarIDsTDC[i];
    SBSTimingHodoscopeBar* Bar = fBars[baridx];
    //std::cout << "id: " << baridx << " mean ToT: " 
    //<< fGoodBarTDCRtotW[i]+fGoodBarTDCLtotW[i] << std::endl;
    //find local maximum first:
    // check for the "middle" elements
    if(0<i && i+1<(int)fGoodBarIDsTDC.size()){
      // if the considered element has larger Time over Threshold than both its direct neighbors, it is considered a local maximum
      if( (fGoodBarTDCRtotW[i]+fGoodBarTDCLtotW[i])>(fGoodBarTDCRtotW[i-1]+fGoodBarTDCLtotW[i-1]) && (fGoodBarTDCRtotW[i]+fGoodBarTDCLtotW[i])>(fGoodBarTDCRtotW[i+1]+fGoodBarTDCLtotW[i+1]) ){
	//std::cout << " case 1 " << std::endl;
	localmax_idx.push_back(i);
	SBSTimingHodoscopeCluster* clus = new SBSTimingHodoscopeCluster(fClusMaxSize, Bar);
	fClusters.push_back(clus);
      }else if( fGoodBarIDsTDC[i+1]-baridx>1 && (fGoodBarTDCRtotW[i]+fGoodBarTDCLtotW[i])>(fGoodBarTDCRtotW[i-1]+fGoodBarTDCLtotW[i-1]) ){
	// if the considered element has no direct neighbor on one side, it only has to have a larger ToT than it's other neighbor to be considered a local maximum
	//std::cout << " case 2 " << std::endl;
	localmax_idx.push_back(i);
	SBSTimingHodoscopeCluster* clus = new SBSTimingHodoscopeCluster(fClusMaxSize, Bar);
	fClusters.push_back(clus);
      }else if( baridx-fGoodBarIDsTDC[i-1]>1 && (fGoodBarTDCRtotW[i]+fGoodBarTDCLtotW[i])>(fGoodBarTDCRtotW[i+1]+fGoodBarTDCLtotW[i+1]) ){
	//std::cout << " case 3 " << std::endl;
	localmax_idx.push_back(i);
	SBSTimingHodoscopeCluster* clus = new SBSTimingHodoscopeCluster(fClusMaxSize, Bar);
	fClusters.push_back(clus);
      }else if(baridx-fGoodBarIDsTDC[i-1]>1 &&  fGoodBarIDsTDC[i+1]-baridx>1){
	// if the element had no direct neighbors, it is a local maximum
	//std::cout << " case 4 " << std::endl;
	localmax_idx.push_back(i);
	SBSTimingHodoscopeCluster* clus = new SBSTimingHodoscopeCluster(fClusMaxSize, Bar);
	fClusters.push_back(clus);
      }
    }else if(i==0 && i+1<(int)fGoodBarIDsTDC.size()){
      if( (fGoodBarTDCRtotW[i]+fGoodBarTDCLtotW[i])>(fGoodBarTDCRtotW[i+1]+fGoodBarTDCLtotW[i+1]) || fGoodBarIDsTDC[i+1]-baridx>1 ){
	//check for first element
	//std::cout << " case 5 " << std::endl;
	localmax_idx.push_back(i);
	SBSTimingHodoscopeCluster* clus = new SBSTimingHodoscopeCluster(fClusMaxSize, Bar);
	fClusters.push_back(clus);
      }
    }else if(i==(int)fGoodBarIDsTDC.size()-1){
      if( (fGoodBarTDCRtotW[i]+fGoodBarTDCLtotW[i])>(fGoodBarTDCRtotW[i-1]+fGoodBarTDCLtotW[i-1]) || baridx-fGoodBarIDsTDC[i-1]>1 ){
	//check for last element
	//std::cout << " case 6 " << std::endl;
	localmax_idx.push_back(i);
	SBSTimingHodoscopeCluster* clus = new SBSTimingHodoscopeCluster(fClusMaxSize, Bar);
	fClusters.push_back(clus);
      }
    }
  }
  
  //localmax_idx and cluster idx should be the same...
  for(size_t i = 0; i<localmax_idx.size(); i++){
    // setting min and max for the subindex that we will loop on 
    // to aggregate elements around the localmax 
    int jmin = std::max(localmax_idx[i]-halfclussize, 0);
    int jmax = std::min(localmax_idx[i]+halfclussize, int(fGoodBarIDsTDC.size()-1));
    //std::cout << localmax_idx[i] << " " << fGoodBarIDsTDC[localmax_idx[i]]
    //<< " " << jmin << " " << jmax << std::endl; 
    for(int j = jmin; j<jmax; j++){
      if(j==localmax_idx[i])continue;//local maximum: bar is already included
      // if the following condition is not met, 
      // it means that there is a non-"good" bar between the local max
      // and the considered bar.
      //std::cout << fGoodBarIDsTDC[j]-fGoodBarIDsTDC[localmax_idx[i]] 
      //<< " " << j-localmax_idx[i] << std::endl;
      if(fGoodBarIDsTDC[j]-fGoodBarIDsTDC[localmax_idx[i]]!=j-localmax_idx[i])continue;
      // check that the new element to be added is compatible with the local maximum in terms of y position and time 
      if(fGoodBarTDCpos[j]-fGoodBarTDCpos[localmax_idx[i]]<fMaxYposDiffCluster &&
	 fGoodBarTDCmean[j]-fGoodBarTDCmean[localmax_idx[i]]<fMaxTimeDiffCluster){
	
	SBSTimingHodoscopeBar* Bar = fBars[fGoodBarIDsTDC[j]];
	fClusters[i]->AddElement(Bar);
      }

    }
  }
  
  /*
  for(size_t i = 0; i<fClusters.size(); i++){
    std::cout << " cluster " << i 
	      << ", size " << fClusters[i]->GetSize() 
	      << " time " << fClusters[i]->GetTmean() 
	      << " position " << fClusters[i]->GetYmean() 
	      << std::endl;
    for(int j = 0; j<fClusters[i]->GetSize(); j++){
      std::cout << " bar id " << fClusters[i]->GetElement(j)->GetBarNum() 
		<< " time " << fClusters[i]->GetElement(j)->GetMeanTime() 
		<< " position " << fClusters[i]->GetElement(j)->GetHitPos()
		<< std::endl;
    }
    std::cout << std::endl;
  }
  */
  //   // this below takes advantage of the fact that the bars are already sorted by ID/geometry i.e. two "good" adjacent bars are guaranteed to be stored back-to-back.
  //   if(baridx-prev_baridx==1 && fClusters.size()>0){
  //     if(!fClusters[fClusters.size()-1]->AddElement(Bar)){
  // 	SBSTimingHodoscopeCluster* clus = new SBSTimingHodoscopeCluster(fClusMaxSize, Bar);
  // 	fClusters.push_back(clus);
  //     }
  //   }else{
  //   }
  //   prev_baridx = baridx;
  // }
  return fClusters.size();
}

Int_t SBSTimingHodoscope::MatchTrack(THaTrack* the_track)
{
  double x_track = the_track->GetX(GetOrigin().Z());
  double y_track = the_track->GetY(GetOrigin().Z());

  int bestmatch=-1;

  double minxdiff = 1.e36;

  //If time resolution is ~100 ps, then position resolution should be
  //something like 100 ps * vScint ~= 3 cm
  
  for(int i=0; i<GetNClusters(); i++){
    SBSTimingHodoscopeCluster* clus = GetCluster(i);
    // if(clus->GetXmean()-clus->GetSize()*SizeCol()/2.<x_track && 
    //    x_track<clus->GetXmean()+clus->GetSize()*SizeCol()/2. &&
    //    fabs( clus->GetYmean() - y_track ) <= SizeRow()/2. ){
    // this is more idiot-proof:
    if( fabs( clus->GetXmean()-x_track ) <= fTrackMatchCutX &&
	fabs( clus->GetYmean()-y_track ) <= fTrackMatchCutY ){
    
      //later when things are calibrated we can do something like: 
      //the_track->SetTime(clus->GetTmean());
      if( bestmatch < 0 || fabs( clus->GetXmean() - x_track ) < minxdiff ){
	minxdiff = fabs(clus->GetXmean() - x_track );
	bestmatch = i;
      }
      //choose the cluster with the smallest difference between the cluster and track
      // X (vertical) positions
    }  
  }
  return bestmatch;
}

/*
 * ConstructHodoscope()
 * called in read database 
 */
Int_t SBSTimingHodoscope::ConstructHodoscope()
{
  Int_t nElements = fElements.size();
  std::cout << "n elements " << nElements << std::endl;
  if( nElements%2!=0 ) {
    Error( Here("ConstructHodoscope"),
	   "N elements for hodoscope is not even, need an even number for a 2 sided detector analysis.");
    return kInitError;
  }
  // make the pmt objects
  if( fNrows!=2 ) {
    Error( Here("ConstructHodoscope"),
	   "fNrows for hodoscope is not 2, which we need for left and right.");
    return kInitError;
  }
  const int nbars = nElements/2;
  DeleteContainer(fPMTMapL);
  DeleteContainer(fPMTMapR);
  fPMTMapL.reserve(nbars);
  fPMTMapR.reserve(nbars);
  int p=0;
  for(int r = 0; r < fNrows; r++) {
    for(int c = 0; c < fNcols[r]; c++) {
      for(int l = 0; l < fNlayers; l++, p++) {
	// get the element
	SBSElement *blk2 = fElementGrid[r][c][l];
	Int_t column = blk2->GetCol();
	Int_t row = blk2->GetRow();
	if(p>=(int)fElements.size()){
	  std::cout << "row " << r << " col " << c << " l " << l  << " p " << p << "/" << fElements.size() << std::endl;
	  std::cout << "element column " << column << " and row " << row << std::endl;
	}
	Double_t tw0 = fTimeWalkPar0[r][c][l];
	Double_t tw1 = fTimeWalkPar1[r][c][l];
	if(row==0){//left// could also use r index
	  if( WithTDC() ){
	    SBSTimingHodoscopePMT *pmtL = new SBSTimingHodoscopePMT(fElements.at(p),
								    tw0,
								    tw1,
								    column, row, p);
	    // column = bar and row = side and p = ID in felements in case needed for debug
	    fPMTMapL.push_back(pmtL);
	  }// if tdc
	  else if( WithADC() ){ // if we have adc we use adc db file and therefore need a temp val for timewalk
	    SBSTimingHodoscopePMT *pmtL = new SBSTimingHodoscopePMT(fElements.at(p),
								    tw0,
								    tw1,
								    column, row, p);
	    fPMTMapL.push_back(pmtL);
	  }// if adc only
	}//left
	if(row==1){//right
	  if( WithTDC() ){
	    SBSTimingHodoscopePMT *pmtR = new SBSTimingHodoscopePMT(fElements.at(p),
								    tw0,
								    tw1,
								    column, row, p);
	    // column = bar and row = side and p = ID in felements in case needed
	    fPMTMapR.push_back(pmtR);
	  }// if tdc
	  else if( WithADC() ){ // if we have adc we use adc db file and therefore need a temp val for timewalk
	    SBSTimingHodoscopePMT *pmtR = new SBSTimingHodoscopePMT(fElements.at(p),
								    tw0,
								    tw1,
								    column, row, p);
	    fPMTMapR.push_back(pmtR);
	  }// if adc only
	}//right
      }//lay
    }//col
  }//row ie side of hodo

  std::cout << "n elements in left pmt array " << fPMTMapL.size() << std::endl;
  std::cout << "n elements in right pmt array " << fPMTMapR.size() << std::endl;
  std::cout << "n elements " << nElements << ", nbars " << nbars << std::endl;
  
  // now we have the arrays of left and right pmts, we need to make the bars
  if( fPMTMapL.size()!=fPMTMapR.size() || fPMTMapL.size()!=(UInt_t)nbars || fPMTMapR.size()!=(UInt_t)nbars) {
    Error( Here("ConstructHodoscope"),
	   "PMT arrays for constructing hodoscope bars not of correct length");
    return kInitError;
  }
  DeleteContainer(fBars);
  fBars.reserve(nbars);

  for(Int_t BarInc=0; BarInc<nbars; BarInc++){
    //std::cout << BarInc << " " << nbars << std::endl;
    // bar constructor is barid, pmt left, pmt right, bar offset mostly for adc sections
    if( WithTDC() ){
      //std::cout << fPMTMapL.at(BarInc)->GetTdcFlag() << std::endl;
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

  //std::cout << "We have filled " << fBars.size() << " bars" << std::endl;
  
  return kOK;
}// construct hodo
/*
 * TimeWalk()
 */
Double_t SBSTimingHodoscope::TimeWalk(Double_t time, Double_t tot, Double_t timewalk0, Double_t timewalk1){
  // at the moment LE versus tot is fit with straight line in calibration
  // tc = LE + tcor
  // where tcor = [0]*TOT + [1], where [0] and [1] are from database
  Double_t tcorr = time - (timewalk0*tot + timewalk1);
  return tcorr;
 }
/*
 * Clear()
 * called at the end of every event
 */
void SBSTimingHodoscope::Clear( Option_t* opt )
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
  
  DeleteContainer(fClusters);  
  /*
  fClusterMult.clear();
  fClusterXmean.clear();
  fClusterYmean.clear();
  fClusterTmean.clear();
  fClusterToTmean.clear();
  */
  
  ClearHodoOutput(fMainClus);
  ClearHodoOutput(fMainClusBars);
  ClearHodoOutput(fOutClus);
  
  // Make sure to call parent class's Clear() also!
  SBSGenericDetector::Clear(opt);
}

void SBSTimingHodoscope::ClearHodoOutput(SBSTimingHodoscopeOutput &out)
{
  out.n.clear();
  out.id.clear();
  out.x.clear();
  out.y.clear();
  out.t.clear();
  out.tot.clear();
  out.tdiff.clear();
  out.trackindex.clear();
}

/*
 * Generic SBSTimingHodoscope destructor
 */
SBSTimingHodoscope::~SBSTimingHodoscope()
{
  // Delete any new objects/instances created here
  DeleteContainer(fBars);
  DeleteContainer(fPMTMapL);
  DeleteContainer(fPMTMapR);
  DeleteContainer(fClusters);
}

SBSTimingHodoscopeCluster* SBSTimingHodoscope::GetCluster(int i)
{
  if(i<GetNClusters())return fClusters[i];
  
  return nullptr;
}
  
