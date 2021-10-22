//////////////////////////////////////////////////////////////////////////
//
// SBSTimingHodoscope class implementation
//
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
  // SBSGenericDetector::SetStoreRawHits(true);
  fDataOutputLevel = 0;//default
  fClusMaxSize = 5;
  fvScint = 0.454*0.299792458; //m/ns
  //Defaults for track match cuts:
  ftDiff0 = -1.35;  //ns
  
  fTrackMatchCutX = 0.05;
  fTrackMatchCutY = 0.15;
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
  fTDCBarOffset = 0;
  Int_t tdcbaroff = 0;
  fADCBarOffset = 0;
  Int_t adcbaroff = 0;
  std::vector<Double_t> timewalkpar0;
  std::vector<Double_t> timewalkpar1;
  DBRequest config_request[] = {
    { "timewalk0map",    &timewalkpar0,      kDoubleV, 0, false }, //parameter for time walk correction
    { "timewalk1map",    &timewalkpar1,      kDoubleV, 0, false }, //parameter for time walk correction
    { "tdcbaroffset",    &tdcbaroff,         kInt,    0, true }, //to allow for cycling through sections
    { "adcbaroffset",    &adcbaroff,         kInt,    0, true }, //to allow for cycling through sections
    { 0 } ///< Request must end in a NULL
  };
  err = LoadDB( file, date, config_request, fPrefix );
  if(err) {
    fclose(file);
    return err;
  }

  DBRequest trackmatch_params[] = {
    { "vscint",          &fvScint, kDouble, 0, 1, 1},
    { "tdiffoffset",     &ftDiff0, kDouble, 0, 1, 1},
    { "trackmatchcutX",  &fTrackMatchCutX, kDouble, 0, 1, 1},
    { "trackmatchcutY",  &fTrackMatchCutY, kDouble, 0, 1, 1},
    { 0 }
  };

  err = LoadDB( file, date, trackmatch_params, fPrefix );
  if(err) {
    fclose(file);
    return err;
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
    if (ypos.size() == fNelem) {
      for (Int_t ne=0;ne<fNelem;ne++) {
	SBSElement* blk= fElements[ne];
	fElements[ne]->SetY(ypos[ne]);
      }
    } else {
      std::cout << " ypos vector too small " << ypos.size() << " # of elements =" << fNelem << std::endl;
    }
  }
  
  DBRequest misc_request[] = {
    { "maxclussize",    &fClusMaxSize,      kInt, 0, true }, //parameter for time walk correction
    { 0 } ///< Request must end in a NULL
  };
  err = LoadDB( file, date, misc_request, fPrefix );
  if(err) {
    fclose(file);
    return err;
  }
  fclose(file);

  // std::cout << "fNelem " << fNelem << std::endl;
  // std::cout << "timewalkpar0.size() " << timewalkpar0.size() << std::endl;
  // std::cout << "timewalkpar1.size() " << timewalkpar1.size() << std::endl;

  // assign the bar offsets
  fTDCBarOffset = tdcbaroff;
  fADCBarOffset = adcbaroff;

  if( WithTDC() || WithADC()){
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
    Double_t Ltwalk0 = pmtL->GetTimeWalkPar0();
    Double_t Ltwalk1 = pmtL->GetTimeWalkPar1();
    Double_t Rtwalk0 = pmtR->GetTimeWalkPar0();
    Double_t Rtwalk1 = pmtR->GetTimeWalkPar1();

    if(WithTDC()){
      if(elL->TDC()->HasData() && elR->TDC()->HasData()){
	//Int_t bar = BarInc;// why redeclare the index?
	// don't need to add offset to tdc since all readout simultaneously
	fGoodBarIDsTDC.push_back(BarInc);
	
	// left hit
	const SBSData::TDCHit &hitL = elL->TDC()->GetGoodHit();
	fGoodBarTDCLle.push_back(hitL.le.val);
	// fGoodBarTDCLle.push_back(hitL.le.raw);
	//.raw is tdc bin, val is corrected using offset and ns/bin
	Double_t LleW = SBSTimingHodoscope::TimeWalk(hitL.le.val,
				(hitL.te.val-hitL.le.val),
				Ltwalk0, Ltwalk1);
	fGoodBarTDCLleW.push_back(LleW);
	fGoodBarTDCLte.push_back(hitL.te.val);
	Double_t LteW = SBSTimingHodoscope::TimeWalk(hitL.te.val,
				(hitL.te.val-hitL.le.val),
				Ltwalk0, Ltwalk1);
	fGoodBarTDCLteW.push_back(LteW);
	fGoodBarTDCLtot.push_back(hitL.te.val-hitL.le.val);
	fGoodBarTDCLtotW.push_back(LteW-LleW);
	
	// right hit
	const SBSData::TDCHit &hitR = elR->TDC()->GetGoodHit();
	fGoodBarTDCRle.push_back(hitR.le.val);//.raw is tdc bin, val is corrected using offset and ns/bin
	Double_t RleW = SBSTimingHodoscope::TimeWalk(hitR.le.val,
				(hitR.te.val-hitR.le.val),
				Rtwalk0, Rtwalk1);
	fGoodBarTDCRleW.push_back(RleW);
	fGoodBarTDCRte.push_back(hitR.te.val);
	Double_t RteW = SBSTimingHodoscope::TimeWalk(hitR.te.val,
				(hitR.te.val-hitR.le.val),
				Rtwalk0, Rtwalk1);
	fGoodBarTDCRteW.push_back(RteW);
	fGoodBarTDCRtot.push_back(hitR.te.val-hitR.le.val);
	fGoodBarTDCRtotW.push_back(RteW-RleW);
	
	// bar properties
	Double_t barmeantime = (LleW + RleW)/2.0;
	fGoodBarTDCmean.push_back(barmeantime);
	Double_t bartimediff = (LleW - RleW);
	fGoodBarTDCdiff.push_back(bartimediff);
	// convert to position? effective velocity times time? should we divide by 2? yes
	// Double_t HorizPos = 0.5 * (bartimediff*1.0e-9) * vScint; // position from L based on timediff and in m
	//Assuming bartimediff is in ns, then horizontal position is

	//AJRP: the -sign is added because the time difference as defined here has
	// a negative correlation with the y of transport coordinates
	// The offset aligns this quantity with the GEM track projection to the hodoscope
	// fvScint ~= 0.454c is the average effective propagation speed as
	// measured by the GEM-TH correlation.
	Double_t HorizPos = -0.5 * (bartimediff-ftDiff0) * fvScint; // position from L based on timediff and in m. 
	fGoodBarTDCpos.push_back(HorizPos);
	fGoodBarTDCvpos.push_back(elR->GetY());
	
	bar->SetMeanTime(barmeantime);
	bar->SetMeanToT( ((hitL.te.val-hitL.le.val)+(hitR.te.val-hitR.le.val))/2. );
	bar->SetTimeDiff(bartimediff);
	bar->SetHitPos(HorizPos);
	bar->SetElementPos(elR->GetY());
	bar->SetLeftHit(hitL);
	bar->SetRightHit(hitR);
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
  int nclusters = DoClustering();
  
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
  int prev_baridx = -10;
  for(int i = 0; i<fGoodBarIDsTDC.size(); i++){
    int baridx = fGoodBarIDsTDC[i];
    SBSTimingHodoscopeBar* Bar = fBars[baridx];
    // this below takes advantage of the fact that the bars are already sorted by ID/geometry i.e. two "good" adjacent bars are guaranteed to be stored back-to-back.
    if(baridx-prev_baridx==1 && fClusters.size()>0){
      if(!fClusters[fClusters.size()-1]->AddElement(Bar)){
	SBSTimingHodoscopeCluster* clus = new SBSTimingHodoscopeCluster(fClusMaxSize, Bar);
	fClusters.push_back(clus);
      }
    }else{
      SBSTimingHodoscopeCluster* clus = new SBSTimingHodoscopeCluster(fClusMaxSize, Bar);
      fClusters.push_back(clus);
    }
    prev_baridx = baridx;
  }
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
	if(p>=fElements.size()){
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

  const int nbars = nElements/2;
  std::cout << "n elements in left pmt array " << fPMTMapL.size() << std::endl;
  std::cout << "n elements in right pmt array " << fPMTMapR.size() << std::endl;
  std::cout << "n elements " << nElements << ", nbars " << nbars << std::endl;
  
  // now we have the arrays of left and right pmts, we need to make the bars
  if( fPMTMapL.size()!=fPMTMapR.size() || fPMTMapL.size()!=nbars || fPMTMapR.size()!=nbars) {
    Error( Here("ConstructHodoscope"),
	   "PMT arrays for constructing hodoscope bars not of correct length");
    return kInitError;
  }
  fBars.clear();
  //I think until here everything is fine... right?
  
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
  
  fClusters.clear();
  
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
  
  // Make sure to call parent class's ClearEvent() also!
  SBSGenericDetector::ClearEvent();
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
  ClearEvent();
}

SBSTimingHodoscopeCluster* SBSTimingHodoscope::GetCluster(int i)
{
  if(i<GetNClusters())return fClusters[i];
  
  return nullptr;
}
  
