// *-- Author :    Guido Maria Urciuoli   12 March 2001

//////////////////////////////////////////////////////////////////////////
//
// SBSRICH
//
// The RICH detector
// Written by Guido Maria Urciuoli, INFN
// Adapted for Hall A Analyzer by Ole Hansen, JLab
// Adapted to SBS by Seamus Riordan, ANL and Eric Fuchey, UConn
//
//////////////////////////////////////////////////////////////////////////

#include "TMath.h"
#include "SBSGRINCH.h"
#include "THaTrack.h"
#include "THaEvData.h"
#include "THaGlobals.h"
#include "THaDetMap.h"
#include "THaSpectrometer.h"
#include "TError.h"
#include "VarDef.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include "THaBenchmark.h"

using namespace std;
// cout<<"stay and stop here!!!"<<endl;
//_____________________________________________________________________________
SBSGRINCH::SBSGRINCH( const char* name, const char* description, 
		  THaApparatus* apparatus ) :
  THaPidDetector(name,description,apparatus), 
  // fMIPs(0),
  // fMaxxMIP(100000), fMinxMIP(-100000), fMaxyMIP(100000), fMinyMIP(-100000),
  fDoResolve(false), fDoTimeFilter(false), //fNseg(0), fXseg(0),
  fTrackX(kBig), fTrackY(kBig)
{
  //keep this line first
  fBench = new THaBenchmark;

  // Normal constructor with name and description

  fHits             = new TClonesArray("SBSGRINCH_Hit",1000);
  fClusters         = new TClonesArray("SBSGRINCH_Cluster",100);
  // fResolvedHits     = new TClonesArray("SBSGRINCH_Hit",1000);
  // fResolvedClusters = new TClonesArray("SBSGRINCH_Cluster",100);
  
  Clear();
  //fDebug=1;
}

//_____________________________________________________________________________
SBSGRINCH::~SBSGRINCH()
{
  // Destructor. Remove variables from global list and free up the memory
  // allocated by us.
  Clear();// so the prgram doesn't complain when deleting clusters
  RemoveVariables();
  delete fHits;
  // delete fResolvedHits;
  delete fClusters;
  // delete fResolvedClusters;
  // delete [] fMIPs;
  // delete [] fXseg;
  delete fBench;
}

//_____________________________________________________________________________
void SBSGRINCH::Clear( Option_t* opt )
{
  // Reset event-by-event data
  if(fDebug)cout << "SBSGRINCH::Clear() " << endl; 
  if( fDoBench ) fBench->Begin("Clear");
  THaPidDetector::Clear(opt);
  if(fDebug)cout << "Clear hits() " << endl; 
  fHits->Clear();
  //fResolvedHits->Clear();
  DeleteClusters();
  if( fDoBench ) fBench->Stop("Clear");
  
  /*
  // reinitialize the map_chan_tdcs
  map_chan_tdcs.clear();
  int a = pow(2, 31), b = a-1;
  pair<int, int> dummy_pair(a, b);
  for(int i = 0; i<fNPMTs; i++){
    map_chan_tdcs.insert(pair<int, pair<int, int > >(i, dummy_pair) );
  }
  */
}

//_____________________________________________________________________________
Int_t SBSGRINCH::ReadDatabase( const TDatime& date )
{
  // 
  // Read the database for this detector.
  // This function is called once at the beginning of the analysis.
  // 'date' contains the date/time of the run being analyzed.
  //
  static const char* const here = "ReadDatabase";
  cout << here << endl;
  // Open the database file
  FILE* fi = OpenFile( date );
  if( !fi ) return kFileError;
 
  // to be inserted in the data base 
  /*
  fiducial_zone_range = 0.05;
  cluster_distribution_sigma = 0.025;
  minimum_chi2_degree_of_freedom = 4;
  acceptable_chi2_prob = 0.01;
  clear_noise_trial_maximum_number = 5;
  epsilon = 1.;
  */
  
  // Storage and default values for non-Double_t and non-member 
  // data from database. Note: These must be Double_t
  // Double_t debug = 0, do_resolve = false,
  //   maxdist, hit_max_number=0, MIP_through_interception = 0;
  
  vector<Int_t>* detmap = 0;
  try {
    // Oddly, putting this container on the stack may cause a stack overflow
    detmap = new vector<Int_t>;
    // Set up a table of tags to read and locations to store values.
    const DBRequest tags[] = {
      {"detmap",       detmap,         kIntV,   0, false },
      {"zpos",     &fZCkovIn,      kDouble, 0, false },
      //{"n_radiator",   &fNradiator,    kDouble, 0, 1},
      //{"l_radiator",   &fLradiator,    kDouble, 0, 1},
      {"nchan",        &fNPMTs,        kInt,    0, false },
      {"nrows",        &fNPMTrows,     kInt,    0, false },
      {"ncols",        &fNPMTcolsMax,  kInt,    0, false },
      {"pmtdistx",     &fPMTdistX,     kDouble, 0, 1},
      {"pmtdisty",     &fPMTdistY,     kDouble, 0, 1},
      {"x_tcpmt",      &fX_TCPMT,      kDouble, 0, 1},
      {"y_tcpmt",      &fY_TCPMT,      kDouble, 0, 1},
	
      // { "l_emission", &l_emission, kDouble, 1, false},
      // { "z_ckovin", &Z_CkovIn, kDouble, 1, false},
      // { "maxdist",    &maxdist, kDouble, 1, true},
      // { "hit_max_number", &hit_max_number, kDouble, 1, true },
      // //{ "MIP_through_interception", &MIP_through_interception, kDouble, 1, true },
      // { "fiducial_zone_range", &fiducial_zone_range, kDouble, 1, true },
      // { "cluster_distribution_sigma", &cluster_distribution_sigma, kDouble, 1, true  },
      // { "acceptable_chi2_prob", &acceptable_chi2_prob, kDouble, 1, true  },
      // { "minimum_chi2_degree_of_freedom", &minimum_chi2_degree_of_freedom, kDouble, 1, true  },
      // { "clear_noise_trial_maximum_number", &clear_noise_trial_maximum_number, kDouble, 1, true },
      // { "epsilon", &epsilon, kDouble, 1, true  },
      // { "do_resolve", &do_resolve, kDouble, 1, true  },
      { "do_timefilter", &fDoTimeFilter, kInt, 1, true  },
      { "debug",         &fDebug,        kInt, 1, true  },
      { "mc_data",       &fMCdata,       kInt, 0, true  },// flag for MC data
      { 0 }
    };
      
    // Read all the tag/value pairs. 
    Int_t err = LoadDB( fi, date, tags );
      
    // Complain and quit if any required values missing.
    if( err ) {    
      if( err>0 )
	Error( Here(here), "Required tag %s%s missing in the "
	       "database.\nDetector initialization failed.",
	       fPrefix, tags[err-1].name );
      fclose(fi);
      return kInitError;
    }

    fPMTmatrixHext = (fNPMTcolsMax-1)*fPMTdistY;
    fPMTmatrixVext = (fNPMTrows-1)*fPMTdistX;
    
    if(fDebug){
      cout << "GRINCH params: " << endl;

      cout << "zpos_front " << fZCkovIn  << endl;
      /*	   << " n_radiator " << fNradiator 
	   << " l_radiator " << fLradiator << endl;
      */
      cout << "npmts " << fNPMTs 
	   << " npmtrows " << fNPMTrows
	   << " npmtcolsmax " << fNPMTcolsMax << endl;
      cout << "pmtdistx " << fPMTdistX
	   << " pmtdisty " << fPMTdistY
	   << " x_tcpmt " << fX_TCPMT
	   << " y_tcpmt " << fY_TCPMT << endl;
      cout << "fPMTmatrixHext " << fPMTmatrixHext
	   << " fPMTmatrixVext " << fPMTmatrixVext << endl;
	
      if(detmap){
	cout << "detmap size = " << detmap->size() << endl;
	for(int i = 0; i<detmap->size(); i++)
	  cout << "detmap[i] = " << detmap->at(i) << endl;
      }
    }
    
    UInt_t flags = 0;//THaDetMap::kFillPlane;//TestBit(kHaveRefChans) ? THaDetMap::kFillRefChan : 0;
    //cout<<"asd################asda"<<flags<<endl;
    // Parse the detector map of the data channels
    if( FillDetMap( *detmap, flags, here ) <= 0 )
      return kInitError;
    
    delete detmap;
  }
  // Catch exceptions here so that we can close the file and clean up
  catch(...) {
    delete detmap;
    fclose(fi);
    throw;
  }
  
  fMaxNumHits = fNPMTs*10;
  
  cout << "Max num hits  = " << fMaxNumHits << ", nPMTs  = " << fNPMTs << endl;
  
  /*
  //initialize channel map here:
  int a = pow(2, 31), b = a-1;
  pair<int, int> dummy_pair(a, b);
  //dummy_pair = std::make_pair<a, b>;
  for(int i = 0; i<fNPMTs; i++){
    map_chan_tdcs.insert(pair<int, pair<int, int > >(i, dummy_pair) );
  }
  */
  
  // TString line, tag;
  // Int_t retval = kOK;
  // int n = 1;

  // fDetMap->Clear();
  /*
  int status;
  while(1) {
    Int_t crate, slot, lo, hi;
    //tag = TString(Form("detmap_%02d", n));
    tag = "detmap";
    if( status = LoadDBvalue( fi, date, tag, line ) != 0 )
      break;
    ISTR inp(line.Data());
    inp >> crate >> slot >> lo >> hi;
    if( !inp )
      goto bad_data;
    if( fDetMap->AddModule( crate, slot, lo, hi ) < 0 )
      goto bad_data;
    n++;
  }
  if(n>1)
    goto done;
  
  bad_data:    
   Error( Here(here), "Problem reading database entry \"%s\". "
  	 "Detector initialization failed.", tag.Data() );
   retval = kInitError;
  
  done:
  fclose(fi);
  return retval;  
  */

  fIsInit = true;
  
  fclose(fi);
  return kOK;
}

//_____________________________________________________________________________
Int_t SBSGRINCH::DefineVariables( EMode mode )
{
  // Define (or delete) global variables of the detector
  
  if( mode == kDefine && fIsSetup ) return kOK;
  fIsSetup = ( mode == kDefine );

  //Hits hits
  RVarDef var1[] = {
    { "nhits",      " number of PMT hits", "GetNumHits()"                    },
    { "hit.pmtnum", " Hit PMT num",        "fHits.SBSGRINCH_Hit.GetPMTNum()" },
    { "hit.xhit",   " PMT hit X",          "fHits.SBSGRINCH_Hit.GetX()"      },
    { "hit.yhit",   " PMT hit y",          "fHits.SBSGRINCH_Hit.GetY()"      },
    { "hit.row",    " PMT hit row",        "fHits.SBSGRINCH_Hit.GetRow()"    },
    { "hit.col",    " PMT hit column",     "fHits.SBSGRINCH_Hit.GetCol()"    },
    { "hit.adc",    " PMT hit ADC",        "fHits.SBSGRINCH_Hit.GetADC()"    },
    { "hit.tdc_r",  " PMT hit TDC rise",   "fHits.SBSGRINCH_Hit.GetTDC_r()"  },
    { "hit.tdc_f",  " PMT hit TDC fall",   "fHits.SBSGRINCH_Hit.GetTDC_f()"  },
    { 0 }
  };
  DefineVarsFromList( var1, mode, "" );// (re)define path here...
  
  RVarDef var2[] = {
    { "nclus",        " number of GRINCH PMT clusters",  "GetNumClusters()"                                 },
    { "clus.size",    " GRINCH cluster size",            "fClusters.SBSGRINCH_Cluster.GetNHits()"           },
    { "clus.x_mean",  " GRINCH cluster X center",        "fClusters.SBSGRINCH_Cluster.GetXcenter()"         },
    { "clus.y_mean",  " GRINCH cluster Y center",        "fClusters.SBSGRINCH_Cluster.GetYcenter()"         },
    { "clus.tr_mean", " GRINCH cluster mean lead time",  "fClusters.SBSGRINCH_Cluster.GetMeanRisingTime()"  },
    { "clus.tf_mean", " GRINCH cluster mean trail time", "fClusters.SBSGRINCH_Cluster.GetMeanFallingTime()" },
    { "clus.adc",     " GRINCH cluster total charge",    "fClusters.SBSGRINCH_Cluster.GetCharge()"          },
    { 0 }
  };
  DefineVarsFromList( var2, mode, "" );// (re)define path here...
    
  return kOK;
}

/*
//_____________________________________________________________________________
//
Double_t SBSGRINCH::Cherenkov_Angle(double mass, double momentum)  const
{
//
//   To calculate the Cherenkov angle for a particle of a certain momentum    
//   and mass 

  return acos(sqrt(1+(mass*mass)/(momentum*momentum))/n_radiator);

}
*/
//_____________________________________________________________________________
Int_t SBSGRINCH::Decode( const THaEvData& evdata )
{
  //Decode RICH data and fill hit array
  if(fDebug){
    cout << "SBSGRINCH::Decode " << endl;
    cout << " Is Init ? " << fIsInit << " Is Physics Trigger ? " << evdata.IsPhysicsTrigger() << endl;
  }
  if( !fIsInit ) return -255;
  if( !evdata.IsPhysicsTrigger() ) return -1;

  Clear();

  if( fDoBench ) fBench->Begin("Decode");

  Int_t nHit = 0;
  SBSGRINCH_Hit* theHit;
  
  bool edge;
  Int_t channel;
  ushort tdctime_raw;
  bool col0_ismaxsize = false;
  
  if(fabs((-fPMTmatrixHext/2.0-fY_TCPMT)/fY_TCPMT)<1.0e-4)col0_ismaxsize = true;
  bool col_ismaxsize;

  //int gchannel;
  int row, col;
  double X, Y;
  double TDC_r, TDC_f;
  double ADC;
  
  if(fDebug)cout << fDetMap->GetSize() << endl;
  for( UInt_t i = 0; i < fDetMap->GetSize(); i++ ) {
    THaDetMap::Module* d = fDetMap->GetModule( i );
    
    if(fDebug)
      cout << "crate " << d->crate << " slot " << d->slot 
	   << " num chans " << evdata.GetNumChan(d->crate, d->slot) << endl;
    
    for( Int_t j = 0; j < evdata.GetNumChan( d->crate, d->slot ); j++) {
      if(col0_ismaxsize){
	col_ismaxsize = true;
      }else{
	col_ismaxsize = false;
      }
      
      UInt_t chan = evdata.GetNextChan( d->crate, d->slot, j );
      
      if( chan > d->hi || chan < d->lo ) continue; // Not one of my channels
      
      // Get the data.
      UInt_t nhit = evdata.GetNumHits( d->crate, d->slot, chan );
      
      if( GetNumHits()+nhit > fMaxNumHits ) {
	Warning("Decode", "Too many hits! Should never ever happen! "
		"Event skipped.");
	cout << "(fMaxNumHits = " << fMaxNumHits << ")" << endl;
	fHits->Clear();
	if( fDoBench ) fBench->Stop("Decode");
	return -2;
      }
      
      if(fDebug)
	cout << "chan " << chan << " nhits = " << nhit << endl;
      
      for (UInt_t hit = 0; hit < nhit; hit++) {
	
	// Fill hit array
	// UInt_t data = evdata.GetData(d->crate,d->slot,chan,hit);
	// UInt_t rawdata = evdata.GetRawData(d->crate,d->slot,chan,hit);
	//assuming "data" is the vetroc word...     
	uint32_t data = evdata.GetData(d->crate,d->slot,chan,hit);
	// NB: the VETROC words are re-added 2^31 to convert them back from int to Uint... 
	// This *should be temporary* !!!
	//data+=pow(2, 31);
	/*
	uint32_t rawdata = evdata.GetRawData(d->crate,d->slot,chan,hit)+pow(2, 31);
	// FIX ME next line is for testing and can be removed later on
	if ( (rawdata & 0xfff) != data ) {
	  cout<< "Something strange happened, "
	    "raw data and data are not consistent" 
	      << endl;
       	}
	*/
	
	//cout << data << endl;
	
	//No we don't want to do that here
	//DecipherVetrocWord(data, edge, channel, tdctime_raw);
	//cout << edge << " " << channel << " " << tdctime_raw << endl;
	
	//if(channel<0)continue;
	// a priori, the channel in the vetroc word and the channel must be equal at this point
	//assert(chan==channel);
	
	// then I add the slot (VETROC) number * 128 
	// => unique channel (NB: 4 VETROC for GRINCH, 16 VETROC for RICH ?)
	channel = chan+d->slot*128;
	
	edge = hit%2;//DUMB!!! but I do it for the time being
	
	if(fDebug)cout << " Channel " << channel 
		       << " GRINCH data " << data 
		       << " edge " << edge << endl;
	
	bool neednewhit = true; 
	
	TIter next( fHits );
	while( SBSGRINCH_Hit* pHit = static_cast<SBSGRINCH_Hit*>( next() )) {
	  if(pHit->GetPMTNum()==channel){
	    if(edge==0 && !pHit->TDC_risSet()){
	      pHit->SetTDC_r(data);
	      neednewhit = false;
	      if(pHit->TDC_fisSet()){
		ADC = (pHit->GetTDC_f()-pHit->GetTDC_r());
		theHit->SetADC(ADC);
	      }
	      break;
	    }
	    if(edge==1 && !pHit->TDC_fisSet()){
	      pHit->SetTDC_f(data);
	      neednewhit = false;
	      if(pHit->TDC_risSet()){
		ADC = (pHit->GetTDC_f()-pHit->GetTDC_r());
		theHit->SetADC(ADC);
	      }
	      break;
	    }
	    
	  }
	}
	if(fDebug)cout << "need to creat new hit ??? " << neednewhit << endl;
	
	// check hit list first
	//for(int i = 0; i<)
	if(neednewhit){
	  theHit = new( (*fHits)[nHit++] ) SBSGRINCH_Hit();
	  theHit->SetPMTNum(channel);
	  
	  div_t d = div(channel, (2*fNPMTcolsMax-1));
	  
	  row = 2*d.quot;
	  col = d.rem;
	  if(d.rem>fNPMTcolsMax-1 && col0_ismaxsize){
	    row+=1;
	    col-=fNPMTcolsMax;
	    col_ismaxsize = false;
	  }
	  
	  if(d.rem>fNPMTcolsMax-2 && !col0_ismaxsize){
	    row+=1;
	    col-=fNPMTcolsMax-1;
	    col_ismaxsize = true;
	  }
	  
	  X = row*fPMTdistX-fPMTmatrixVext/2.0;
	  if(col_ismaxsize){
	    Y = col*fPMTdistX-fPMTmatrixHext/2.0;
	  }else{
	    Y = col*fPMTdistX-(fPMTmatrixHext-fPMTdistX)/2.0;
	  }
	  
	  theHit->SetRow(row);
	  theHit->SetCol(col);
	  theHit->SetX(X);
	  theHit->SetY(Y);
	  if(edge==0){
	    theHit->SetTDC_r(TDC_r);
	  }else{
	    theHit->SetTDC_f(TDC_f);
	  }
	  theHit->SetFlag(0);
	  theHit->SetVeto(0);
	}
	
	/*
	if(fDebug)
	  cout << "hit " << hit << ": channel " << channel << " edge " << edge 
	       << ":  => " << map_chan_tdcs[channel].first << " <=? time " << data 
	       << " <=? " << map_chan_tdcs[channel].second << endl;
	
	//map_chan_tdcs[channel].resize(map_chan_tdcs[channel].size()+1);
	if(edge==0 && data<map_chan_tdcs[channel].second){
	  map_chan_tdcs[channel].first = data;
	}
	if(edge==1 && data>map_chan_tdcs[channel].first){
	  map_chan_tdcs[channel].second = data;
	}
	*/
      }
    }
  }
  if(fDebug)cout << "Finished decoding" << endl;
  /*
  for(int channel = 0; channel<fNPMTs; channel++){
    //fill fHits !!!
    //only if both egdes have been initialized and the falling edge is above the rising edge
    //
    
    if(map_chan_tdcs[channel].first>=0 && map_chan_tdcs[channel].second>=map_chan_tdcs[channel].first){
      if(fDebug)cout << "channel: " << map_chan_tdcs[channel].first << " " << map_chan_tdcs[channel].second << endl;
      if(col0_ismaxsize){
	col_ismaxsize = true;
      }else{
      	col_ismaxsize = false;
      }
      
      theHit = new( (*fHits)[nHit++] ) SBSGRINCH_Hit();
      theHit->SetPMTNum(channel);
      
      div_t d = div(channel, (2*fNPMTcolsMax-1));
      
      //cout << d.quot << " " << d.rem << endl;
      row = 2*d.quot;
      col = d.rem;
      if(d.rem>fNPMTcolsMax-1 && col0_ismaxsize){
	row+=1;
	col-=fNPMTcolsMax;
	col_ismaxsize = false;
      }
      if(d.rem>fNPMTcolsMax-2 && !col0_ismaxsize){
	row+=1;
	col-=fNPMTcolsMax-1;
	col_ismaxsize = true;
      }
      
      X = row*fPMTdistX-fPMTmatrixVext/2.0;
      if(col_ismaxsize){
	Y = col*fPMTdistX-fPMTmatrixHext/2.0;
      }else{
	Y = col*fPMTdistX-(fPMTmatrixHext-fPMTdistX)/2.0;
      }
      TDC_r = map_chan_tdcs[channel].first;
      TDC_f = map_chan_tdcs[channel].second;
      ADC = (TDC_f-TDC_r);
      
      if(fDebug) {
	cout << "Event number " << evdata.GetEvNum() << endl;
	cout << " building PMT hits: channel " << channel 
	     << " tdc 0 " <<  map_chan_tdcs[channel].first
	     << " tdc 1 " <<  map_chan_tdcs[channel].second
	     << endl;
	
	cout << "PMT row " << row << " (X = " << X 
	     << " m); col " << col << " (Y = " << Y << " m)" << endl;
	cout << "TDC_r " << TDC_r << ", TDC_f " << TDC_f << endl;
      }
      
      theHit->SetRow(row);
      theHit->SetCol(col);
      theHit->SetX(X);
      theHit->SetY(Y);
      theHit->SetADC(ADC);
      theHit->SetTDC_r(TDC_r);
      theHit->SetTDC_f(TDC_f);
      theHit->SetFlag(0);
      theHit->SetVeto(0);
    }
  }
  */
  if( fDoBench ) fBench->Stop("Decode");
  return GetNumHits();
}

//_____________________________________________________________________________
Int_t SBSGRINCH::CoarseProcess( TClonesArray& tracks )
{
  // Coarse processing of the RICH. We do nothing here since we
  // need precise tracking to be able to do anything useful.

  return 0;
}


//_____________________________________________________________________________
Int_t SBSGRINCH::FineProcess( TClonesArray& tracks )
{
  if(fDebug)cout << "Begin Fine Process" << endl;
  // The main RICH processing method. Here we
  //  
  // - Identify clusters
  // - Attempt to match a traceto these
  // - Calculate particle probabilities

  if( !fIsInit ) return -255;

  if (fMaxNumHits<GetNumHits()) {  
    return -25;
  }  
  
  // Clusters reconstructed here
  if( FindClusters() == 0 ) { 
    return -1;
  }
  
  if(fDebug)cout << "DoTimeFitter ? " << fDoTimeFilter << " Tracks::GetLast ? " << tracks.GetLast() << endl;

  // if( fDoResolve )
  //   ResolveClusters();
  //cout<<"stay and stop here!!!"<<endl;
  if( fDoBench ) fBench->Begin("FineProcess");

  if(fDoTimeFilter) CleanClustersWithTime();
  
  // Clusters matched with tracks here (obviously if there are any tracks to match)
  if(tracks.GetLast()>0){
    MatchClustersWithTracks(tracks);
  }
  
  if( fDoBench ) fBench->Stop("FineProcess");
  if(fDebug)cout << "Done fine process " << endl;
  return 0;
}

/*
//__________________________________________________________________________
void SBSGRINCH::DecipherVetrocWord(uint32_t VetrocWord, Bool_t& edge, Short_t& chan, UShort_t& time)
{
  //"strip down" the vetroc word to obtain the information

  // bool header[8]; // = {0, 0, 0, 0, 0, 0, 1, 1};
  // bool tdc[16];
  uint16_t header = 0;
  uint16_t channel = 0;
  uint16_t tdc_time = 0;

  bool ibit;
  
  for(int i = 0; i<16; i++){
    if(i<8){
      ibit = (VetrocWord >> (i+24) ) &1UL;
      header ^= (-ibit ^ header) & (1U << i);
      // cout << "vetroc word bit " << i+24 << " = " << ((VetrocWord >> (i+24) ) &1UL) 
      // 	   << ", header bit " << i << " = " << ((header >> (i) ) & 1)<< endl;
      ibit = (VetrocWord >> (i+16) ) &1UL;
      channel ^= (-ibit ^ channel) & (1U << i);
      // cout << "vetroc word bit " << i+16 << " = " << ((VetrocWord >> (i+16) ) &1UL) 
      // 	   << ", channel bit " << i << " = " << ((channel >> (i) ) & 1)<< endl;
    }
    ibit = (VetrocWord >> (i) ) &1UL;
    tdc_time ^= (-ibit ^ tdc_time) & (1U << i);
    // cout << "vetroc word bit " << i << " = " << ((VetrocWord >> (i) ) & 1UL) 
    // 	 << ", time bit " << i << " = " << ((tdc_time >> (i) ) & 1U)<< endl;
  }
  
  //cout << header << " " << channel << " " << tdc_time << endl;
  
  chan = channel;
  time = tdc_time;
  
  if(header == 192){
    edge = 0;
  }else if(header == 196){
    edge = 1;
  }else{
    chan = -kBig;
    time = -kBig;
  }
  
    uint32_t TDCvetrocWord0, TDCvetrocWord1;
    TDCvetrocWord0 = TDCvetrocWord1 = 0;
    
    
    //cout << "Vetroc words: init: " << TDCvetrocWord0 << " " << TDCvetrocWord1 << endl;
    for(int i = 0; i<16; i++){
    if(i<8){
    TDCvetrocWord0 ^= (-header[i] ^ TDCvetrocWord0) & (1 << (i+24));
    TDCvetrocWord1 ^= (-header[i] ^ TDCvetrocWord1) & (1 << (i+24));
    channel[i] = (ipmt >> i) & 1;
    TDCvetrocWord0 ^= (-channel[i] ^ TDCvetrocWord0) & (1 << (i+16));
    TDCvetrocWord1 ^= (-channel[i] ^ TDCvetrocWord1) & (1 << (i+16));
    }
    tdc[i] = (fTDCArrays.at(idet).first[ipmt] >> i) & 1;
    TDCvetrocWord0 ^= (-tdc[i] ^ TDCvetrocWord0) & (1 << i);
    tdc[i] = (fTDCArrays.at(idet).second[ipmt] >> i) & 1;
    TDCvetrocWord1 ^= (-tdc[i] ^ TDCvetrocWord1) & (1 << i);
    }
    
    TDCvetrocWord1 ^= (-1 ^ TDCvetrocWord1) & (1 << 26);
    //cout << "Vetroc words: test: " << TDCvetrocWord0 << " " << TDCvetrocWord1 << endl;
    
    hit.fTDC[0] = TDCvetrocWord0;
    hit.fTDC[1] = TDCvetrocWord1;
  
  
  return;
}
*/

//__________________________________________________________________________
void SBSGRINCH::DeleteClusters()
{
  //Delete all clusters
  if(fDebug)cout << "Clear Clusters" << endl;
  fClusters->Clear("C");
  //fResolvedClusters->Clear("C");
  return;
}

//__________________________________________________________________________
Int_t SBSGRINCH::FindClusters()
{
  // Group the hits that are currently in the array into clusters.
  // Return number of clusters found.
 
  // // minimum distance between two pads.
  const double par = sqrt(fPMTdistX*fPMTdistX+fPMTdistY*fPMTdistY);//2.0; 

  // // maximum distance in X between two fired pads to be in the same cluster.
  // const double par2 = 0.1;//PAD_SIZE_X+0.1;  

  // // maximum distance in Y between two fired pads to be in the same cluster.
  // const double par3 = 0.1;//PAD_SIZE_Y+0.1;

  DeleteClusters();

  if( fDoBench ) fBench->Begin("FindClusters");

  SBSGRINCH_Hit* theHit;
  SBSGRINCH_Cluster* theCluster;
  Int_t nHits  = GetNumHits();
  Int_t nClust = GetNumClusters();
  Int_t hit_1stclusmatch;
  
  //loop on hits
  for( int k=0; k<nHits; k++ ) {
    if( !(theHit = GetHit(k))) continue;
    hit_1stclusmatch = -1;
    // HitFlag not equal 0: The Hit was alredy processed
    // HitFlag equal 0: insert the Hit as first element of the cluster
    
    if( theHit->GetFlag() == 0 ) {
      // we loop on *all* existing clusters to check which of them the hit can already be associated to.
      for(int i_cl = 0; i_cl<nClust; i_cl++){
	theCluster = GetCluster(i_cl);
	if(fDebug)cout << " i_cl " << i_cl << " <? nClust " << nClust << " => theCluster ? " << theCluster << endl;
	if(theCluster->IsNeighbor(theHit, par)){
	  if(theHit->GetFlag() == 0){
	    // if it is the first cluster the hit may be associated to,
	    // the hit is inserted to the cluster, and the cluster number is stored
	    theCluster->Insert(theHit);
	    theHit->SetFlag(1);
	    hit_1stclusmatch = i_cl;
	  }else{
	    SBSGRINCH_Cluster* refCluster = GetCluster(hit_1stclusmatch);
	    // if it is not the first cluster the hit is associated to, 
	    // the cluster is added to the first cluster associated to the hit
	    // and removed after
	    refCluster->MergeCluster(*theCluster);
	    // theCluster->AddCluster(GetCluster(hit_1stclusmatch));
	    theCluster->Clear("F");
	    //cout << "i_cl " << i_cl << " fClusters->At(i_cl) ? " << fClusters->At(i_cl) << endl;
	    fClusters->RemoveAt(i_cl);
	    fClusters->Compress();
	    //cout << "i_cl " << i_cl << " fClusters->At(i_cl) ? " << fClusters->At(i_cl) << endl;
	    //i_cl--;//??
	    nClust--;
	  }
	}
      }
      
      // if the hit cannot be asscoiated to any existing cluster, we create a new one
      if(theHit->GetFlag() == 0){
	theCluster = new( (*fClusters)[nClust++] ) SBSGRINCH_Cluster();
	theHit->SetFlag(1);
	theCluster->Insert( theHit );
      }
      
    }
  }
  // now find if there are clusters genetrated by overlapping of clusters
  // (number of overlapping clusters = number of local maximums in a cluster).
  // if( fDoResolve ) {
  //   for( int k = 0; k < nClust; k++ ) {
  //     if( !(theCluster = GetCluster(k))) continue;
  //     //theCluster->FindLocalMaximumNumber();
  //   }
  // }
  if(fDebug)cout << "Finished clustering " << endl;
  if( fDoBench ) fBench->Stop("FindClusters");
  return nClust;
}

//_____________________________________________________________________________
Int_t SBSGRINCH::MatchClustersWithTracks( TClonesArray& tracks )
{
  // std::vector<bool> TrackIsAssociated;
  // TrackIsAssociated.resize(tracks.GetSize());
  SBSGRINCH_Cluster* theCluster;
  THaTrack* theTrack;
  Int_t nClust = GetNumClusters();

  // for(int l = 0; l<tracks.GetSize(); l++){
  //   TrackIsAssociated[l] = false;
  // }
  Int_t Nassociated = 0;
  for( int k = 0; k < nClust; k++ ) {
    if( !(theCluster = GetCluster(k))) continue;
    for(int l = 0; l<tracks.GetSize(); l++){
      theTrack = (THaTrack*)tracks.At(l);
      //Assuming the TClonesArray is indeed a THaTrack TClonesArray
      
      double Xclus = theCluster->GetXcenter()*100.0;// in cm...
      double Yclus = theCluster->GetYcenter()*100.0;// in cm...
      
      double Xtrk = theTrack->GetX(fZCkovIn)*100.0;// in cm...
      double dXtrk = theTrack->GetTheta()*100.0;// in cm...
      double Ytrk = theTrack->GetY(fZCkovIn)*100.0;// in cm...
      
      /* 
	 // fx0,1,2_min,max (6); fx0,1,2_p0,1_corr (6); fsigma_x (1);
	 // fdx0,1,2_min,max (6); fdx0,1,2_p0,1_corr (6); fsigma_dx (1);
	 // fy_p0,1,2_corr (3); f(5?)sigma_y (1);
	 // fnsigmas
	 // => 31 parameters... but "only" 9 lines in the DB (if we pass vectors...)
	 // => still a lot but might be worth it if we use something similar for the RICH
	 
	 bool xsel_0 = ( (fx_[0]<=Xtrk && Xtrk<=fx_[1]) && fabs( Xclus-(fx_corr[0]+fx_corr[1]*Xtrk) )<fsigma_x*fnsigmas );
	 bool xsel_1 = ( (fx_[2]<=Xtrk && Xtrk<=fx_[3]) && fabs( Xclus-(fx_corr[2]+fx_corr[3]*Xtrk) )<fsigma_x*fnsigmas );
	 bool xsel_2 = ( (fx_[4]<=Xtrk && Xtrk<=fx_[5]) && fabs( Xclus-(fx_corr[4]+fx_corr[5]*Xtrk) )<fsigma_x*fnsigmas );
	 bool xsel = xsel_0 || xsel_1 || xsel_2;
	 
	 bool dxsel_0 = ( (fx_[0]<=dXtrk && dXtrk<=fx_[1]) && fabs( (Xclus-Xtrk)-(fdx_corr[0]+fdx_corr[1]*Xtrk) )<fsigma_dx*fnsigmas );
	 bool dxsel_1 = ( (fx_[2]<=dXtrk && dXtrk<=fx_[3]) && fabs( (Xclus-Xtrk)-(fdx_corr[2]+fdx_corr[3]*Xtrk) )<fsigma_dx*fnsigmas );
	 bool dxsel_2 = ( (fx_[4]<=dXtrk && dXtrk<=fx_[5]) && fabs( (Xclus-Xtrk)-(fdx_corr[4]+fdx_corr[5]*Xtrk) )<fsigma_dx*fnsigmas );
	 bool dxsel = dxsel_0 || dxsel_1 || dxsel_2;
	 
	 bool ysel = ( fabs( Yclus-(fy_corr[0]+fy_corr[1]*Ytrk+fy_corr[2]*Ytrk*Ytrk) )<fsigma_y*fnsigmas );
      */
      
      bool xsel_0 = ( (-75.<=Xtrk && Xtrk<=-50.) && fabs( Xclus-(36.91610+1.621343*Xtrk) )<2.541292*5 );
      bool xsel_1 = ( (-55.<=Xtrk && Xtrk<= 50.) && fabs( Xclus-( 7.05105+1.365355*Xtrk) )<2.541292*5 );
      bool xsel_2 = ( ( 45.<=Xtrk && Xtrk<= 80.) && fabs( Xclus-(-7.44171+1.277790*Xtrk) )<2.541292*5 );
      bool xsel = xsel_0 || xsel_1 || xsel_2;
 
      bool dxsel_0 = ( (-0.260<=dXtrk && dXtrk<=-0.085) && fabs( (Xclus-Xtrk)-( 15.4112   +109.019*dXtrk) )<1.63893*5 );
      bool dxsel_1 = ( (-0.180<=dXtrk && dXtrk<= 0.220) && fabs( (Xclus-Xtrk)-(  0.0295942+123.143*dXtrk) )<1.63893*5 );
      bool dxsel_2 = ( ( 0.080<=dXtrk && dXtrk<= 0.300) && fabs( (Xclus-Xtrk)-(-13.8429   +100.121*dXtrk) )<1.63893*5 );
      bool dxsel = dxsel_0 || dxsel_1 || dxsel_2;
      
      bool ysel = ( fabs( Yclus-(-2.20127e-01-0.091446*Ytrk+ 0.00846439 *Ytrk*Ytrk) )<1.04718*5 );
      
      // compare
      if( ysel && xsel && dxsel ){
	theCluster->SetTrack(theTrack);
       	Nassociated++;
      }
    }
    
  }
  return(Nassociated);
}

//__________________________________________________________________________
Int_t SBSGRINCH::CleanClustersWithTime()// pass as an option...
{
  SBSGRINCH_Cluster* theCluster;
  Int_t nClust = GetNumClusters();

  for( int k = 0; k < nClust; k++ ) {
    if( !(theCluster = GetCluster(k))) continue;
    
    
  }
  
  return(1);
}

/*
//__________________________________________________________________________
Int_t SBSGRINCH::ResolveClusters()
{
  // Resolve the clusters. To be calle only after FindClusters.
  // Return the new number of clusters.

  if( fDoBench ) fBench->Begin("ResolveClusters");

  Int_t nClust = GetNumClusters();

  fResolvedHits->Clear("C");
  fResolvedClusters->Clear("C");
  Int_t nResolvedClusters = 0;

  for( int k = 0; k < nClust; k++ ) {
    SBSGRINCH_Cluster* theCluster = GetCluster(k);
    if( !theCluster ) continue;
    Int_t Number = 1;//theCluster->GetLocalMaximumNumber();
    if( (Number > 1) ) {// && !(theCluster->IsMIP())
      TIter next(theCluster->GetHitList());
      // Make a resolved cluster for each local maximum of the cluster
      while( const SBSGRINCH_Hit* pHit = static_cast<SBSGRINCH_Hit*>(next())) {
	if( pHit->GetVeto() != 0 ) continue;
	// found a local maximum
	SBSGRINCH_Cluster* the_Resolved_Cluster = 
	  new( (*fResolvedClusters)[nResolvedClusters++] ) SBSGRINCH_Cluster();

	// the real work is done here -
	theCluster->FindResolvedClusterElements
	  ( pHit, the_Resolved_Cluster, fResolvedHits );
      }
    }
    else {
      new( (*fResolvedClusters)[nResolvedClusters++] ) 
	SBSGRINCH_Cluster( *theCluster );
    }
  }
  if( fDoBench ) fBench->Stop("ResolveClusters");
  return nResolvedClusters;

}
*/

//_____________________________________________________________________________
void SBSGRINCH::PrintBenchmarks() const
{
  // Print benchmark results
  
  if( !fDoBench )
    return;

  fBench->Print("Clear");
  fBench->Print("Decode");
  fBench->Print("FineProcess");
  fBench->Print("FindClusters");
  fBench->Print("FindMIP");
  //fBench->Print("ResolveClusters");

  cout << endl << "Breakdown of time spent in FineProcess:" << endl;
  fBench->Print("RecoAng");
  fBench->Print("Get_phi_photon");
  fBench->Print("Get_theta_photon");
  
  return;
}

ClassImp(SBSGRINCH)

