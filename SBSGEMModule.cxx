#include <iostream>

#include "SBSGEMModule.h"
#include "TDatime.h"
#include "THaEvData.h"

//This should not be hard-coded, I think, but read in from the database (or perhaps not, if it never changes?)
const int APVMAP[128] = {1, 33, 65, 97, 9, 41, 73, 105, 17, 49, 81, 113, 25, 57, 89, 121, 3, 35, 67, 99, 11, 43, 75, 107, 19, 51, 83, 115, 27, 59, 91, 123, 5, 37, 69, 101, 13, 45, 77, 109, 21, 53, 85, 117, 29, 61, 93, 125, 7, 39, 71, 103, 15, 47, 79, 111, 23, 55, 87, 119, 31, 63, 95, 127, 0, 32, 64, 96, 8, 40, 72, 104, 16, 48, 80, 112, 24, 56, 88, 120, 2, 34, 66, 98, 10, 42, 74, 106, 18, 50, 82, 114, 26, 58, 90, 122, 4, 36, 68, 100, 12, 44, 76, 108, 20, 52, 84, 116, 28, 60, 92, 124, 6, 38, 70, 102, 14, 46, 78, 110, 22, 54, 86, 118, 30, 62, 94, 126};


SBSGEMModule::SBSGEMModule( const char *name, const char *description,
			    THaDetectorBase* parent ):
  THaSubDetector(name,description,parent)
{
  // FIXME:  To database
  //Set Default values for fZeroSuppress and fZeroSuppressRMS:
  fZeroSuppress    = kFALSE;
  fZeroSuppressRMS = 5.0;

  //Set default values for decode map parameters:
  fN_APV25_CHAN = 128;
  fN_MPD_TIME_SAMP = 6;
  fMPDMAP_ROW_SIZE = 9;
    
  // for( Int_t i = 0; i < N_MPD_TIME_SAMP; i++ ){
  //     fadc[i] = NULL;
  // }

  return;
}

SBSGEMModule::~SBSGEMModule() {
  // if( fStrip ){
  //     fadc0 = NULL;
  //     fadc1 = NULL;
  //     fadc2 = NULL;
  //     fadc3 = NULL;
  //     fadc4 = NULL;
  //     fadc5 = NULL;

  //     for( Int_t i = 0; i < N_MPD_TIME_SAMP; i++ ){
  //         delete fadc[i];
  //         fadc[i] = NULL;
  //     }
  //     delete fPedestal;
  //     fPedestal = NULL;
  //     delete fStrip;
  //     fStrip = NULL;
  // }

  return;
}

Int_t SBSGEMModule::ReadDatabase( const TDatime& date ){
  std::cout << "[SBSGEMModule::ReadDatabase]" << std::endl;

  Int_t status;

  FILE* file = OpenFile( date );
  if( !file ) return kFileError;

  std::vector<Double_t> rawped;
  std::vector<Double_t> rawrms;

  const DBRequest request[] = {
    { "chanmap",        &fChanMapData,        kIntV, 0, 0},
    { "ped",            &rawped,        kDoubleV, 0, 1},
    { "rms",            &rawrms,        kDoubleV, 0, 1},
    {}
  };
  status = LoadDB( file, date, request, fPrefix );
  fclose(file);

  Int_t nentry = fChanMapData.size()/fMPDMAP_ROW_SIZE;
  for( Int_t mapline = 0; mapline < nentry; mapline++ ){
    mpdmap_t thisdata;
    thisdata.crate  = fChanMapData[0+mapline*MPDMAP_ROW_SIZE];
    thisdata.slot   = fChanMapData[1+mapline*MPDMAP_ROW_SIZE];
    thisdata.mpd_id = fChanMapData[2+mapline*MPDMAP_ROW_SIZE];
    thisdata.gem_id = fChanMapData[3+mapline*MPDMAP_ROW_SIZE];
    thisdata.adc_id = fChanMapData[4+mapline*MPDMAP_ROW_SIZE];
    thisdata.i2c    = fChanMapData[5+mapline*MPDMAP_ROW_SIZE];
    thisdata.pos    = fChanMapData[6+mapline*MPDMAP_ROW_SIZE];
    thisdata.invert = fChanMapData[7+mapline*MPDMAP_ROW_SIZE];
    thisdata.axis   = fChanMapData[8+mapline*MPDMAP_ROW_SIZE];
    fMPDmap.push_back(thisdata);
  }

  std::cout << fName << " mapped to " << nentry << " APV25 chips" << std::endl;

  // FIXME:  make sure to delete if already initialized
  fStrip    = new Int_t [N_APV25_CHAN*nentry];


  for( Int_t i = 0; i < N_MPD_TIME_SAMP; i++ ){
    fadc[i] = new Int_t [N_APV25_CHAN*nentry];
    for( Int_t j = 0; j < N_MPD_TIME_SAMP; j++ ){
      fadc[i][j] = 0.0;
    }
  }
  fadc0 = fadc[0];
  fadc1 = fadc[1];
  fadc2 = fadc[2];
  fadc3 = fadc[3];
  fadc4 = fadc[4];
  fadc5 = fadc[5];

  fPedestal = new Double_t [N_APV25_CHAN*nentry];
  fRMS      = new Double_t [N_APV25_CHAN*nentry];
  for( Int_t i = 0; i < N_APV25_CHAN*nentry; i++ ){
    // FIXME needs to read in pedestal map
    fPedestal[i] = 0.0;
    fRMS[i] = 0.0;
  }

  for( UInt_t i = 0; i < rawped.size(); i++ ){
    if( (i % 2) == 1 ) continue;
    int idx = (int) rawped[i];
	
    if( idx < N_APV25_CHAN*nentry ){
      fPedestal[idx] = rawped[i+1];
    } else {
		
      std::cout << "[SBSGEMModule::ReadDatabase]  WARNING: " << " strip " << idx  << " listed but not enough strips in cratemap" << std::endl;
    }
  }

  for( UInt_t i = 0; i < rawrms.size(); i++ ){
    if( (i % 2) == 1 ) continue;
    int idx = (int) rawrms[i];
    if( idx < N_APV25_CHAN*nentry ){
      fRMS[idx] = rawrms[i+1];
    } else {
      std::cout << "[SBSGEMModule::ReadDatabase]  WARNING: " << " strip " << idx  << " listed but not enough strips in cratemap" << std::endl;
    }
  }


  return 0;
}

Int_t SBSGEMModule::DefineVariables( EMode mode ) {
  if( mode == kDefine and fIsSetup ) return kOK;
  fIsSetup = ( mode == kDefine );

  RVarDef vars[] = {
    { "nch",   "Number of channels",   "fNch" },
    { "strip", "Strip number mapping", "fStrip" },
    { "adc0", "Strip number mapping", "fadc0" },
    { "adc1", "Strip number mapping", "fadc1" },
    { "adc2", "Strip number mapping", "fadc2" },
    { "adc3", "Strip number mapping", "fadc3" },
    { "adc4", "Strip number mapping", "fadc4" },
    { "adc5", "Strip number mapping", "fadc5" },
    { 0 },
  };


  Int_t ret = DefineVarsFromList( vars, mode );

  if( ret != kOK )
    return ret;
    
  return kOK;
    
}

void    SBSGEMModule::Clear( Option_t* opt){
  fNch = 0;
  fIsDecoded = false;
  
  return;
}

Int_t   SBSGEMModule::Decode( const THaEvData& evdata ){
  //    std::cout << "[SBSGEMModule::Decode " << fName << "]" << std::endl;

  int i;

  fNch = 0;
  for (std::vector<mpdmap_t>::iterator it = fMPDmap.begin() ; it != fMPDmap.end(); ++it){
    Int_t effChan = it->mpd_id << 8 | it->adc_id;
    // Find channel for this crate/slot

    Int_t nchan = evdata.GetNumChan( it->crate, it->slot );

    //        printf("nchan = %d\n", nchan );

    for( Int_t ichan = 0; ichan < nchan; ++ichan ) {
      Int_t chan = evdata.GetNextChan( it->crate, it->slot, ichan );
      if( chan != effChan ) continue; // not part of this detector


      Int_t nsamp = evdata.GetNumHits( it->crate, it->slot, chan );
      assert(nsamp%N_MPD_TIME_SAMP==0);
      Int_t nstrips = nsamp/N_MPD_TIME_SAMP;

      // Loop over all the strips and samples in the data
      Int_t isamp = 0;
      for( Int_t istrip = 0; istrip < nstrips; ++istrip ) {
	assert(isamp<nsamp);
	Int_t strip = evdata.GetRawData(it->crate, it->slot, chan, isamp);
	assert(strip>=0&&strip<128);
	// Madness....   copy pasted from stand alone decoder
	// I bet there's a more elegant way to express this
	//Int_t RstripNb= 32*(strip%4)+ 8*(int)(strip/4)- 31*(int)(strip/16);
	//RstripNb=RstripNb+1+RstripNb%4-5*(((int)(RstripNb/4))%2);
	// New: Use a pre-computed array from Danning to skip the above
	// two steps.
	Int_t RstripNb = APVMAP[strip];
	RstripNb=RstripNb+(127-2*RstripNb)*it->invert;
	Int_t RstripPos = RstripNb + 128*it->pos;
	strip = RstripPos;

	fStrip[fNch] = strip;
	for(Int_t adc_samp = 0; adc_samp < N_MPD_TIME_SAMP; adc_samp++) {

	  fadc[adc_samp][fNch] =  evdata.GetData(it->crate, it->slot,
						 chan, isamp++) - fPedestal[strip];

	  assert( ((UInt_t) fNch) < fMPDmap.size()*N_APV25_CHAN );
	}
	assert(strip>=0); // Make sure we don't end up with negative strip numbers!

	// Zero suppression
	if(!fZeroSuppress ||
	   ( fRMS[strip] > 0.0 && fabs(fadc[2][fNch])/
	     fRMS[strip] > fZeroSuppressRMS ) ){
	  fNch++;
	}
      }
    }

  }

  //    std::cout << fName << " channels found  " << fNch << std::endl;

  fIsDecoded = true;
  
  return 0;
}

void SBSGEMModule::find_2Dhits(){
  find_clusters(false); //u strips
  find_clusters(true); //v strips
}

void SBSGEMModule::find_clusters_1D(bool axis){
  
  if( !fIsDecoded ){
    cout << "find_clusters invoked before decoding for GEM Module " << GetName() << ", doing nothing" << endl;
    return;
  }

  UShort_t maxsep = axis ? fMaxNeighborsV_totalcharge : fMaxNeighborsU_totalcharge;
  UShort_t maxsepcoord = axis ? fMaxNeighborsV_hitpos : fMaxNeighborsU_hitpos; 

  UInt_t Nstrips = axis ? fNstripsV : fNstripsU;
  
  Double_t pitch = axis ? fVStripPitch : fUStripPitch;
  
  std::set<UShort_t> striplist;  //sorted list of strips for 1D clustering
  std::map<UShort_t, UInt_t> hitindex; //key = strip ID, mapped value = index in decoded hit array, needed to access the other information efficiently:
  
  for( int ihit=0; ihit<fNstrips_hit; ihit++ ){
    if( fAxis[ihit] == axis && fKeepStrip[ihit] ){
      bool newstrip = (striplist.insert( fStrip[ihit] ) ).second;

      if( newstrip ){ //should always be true:
	hitindex[fStrip[ihit]] = ihit;
      }
    }
  }

  std::set<UShort_t> localmaxima;
  std::map<UShort_t,bool> islocalmax;

  for( std::set<UShort_t>::iterator i=striplist.begin(); i != striplist.end(); ++i ){
    int strip = *i;
    //int hitidx = hitindex[strip];
    islocalmax[strip] = false;

    double sumstrip = fADCsums[hitindex[strip]];
    double sumleft = 0.0;
    double sumright = 0.0;

    if( striplist.find( strip - 1 ) != striplist.end() ){
      sumleft = fADCsums[hitindex[strip-1]]; //if strip - 1 is found in strip list, hitindex is guaranteed to have been initialized above
    }
    if( striplist.find( strip + 1 ) != striplist.end() ){
      sumright = fADCsums[hitindex[strip+1]];
    }

    if( sumstrip >= sumleft && sumstrip >= sumright ){ //new local max:
      islocalmax[strip] = true;
      localmaxima.insert( strip );
    } 
  } // end loop over list of strips along this axis:

  vector<int> peakstoerase; 

  //now calculate "prominence" for all peaks and erase "insignificant" peaks:

  for( std::set<int>::iterator i=localmaxima.begin(); i != localmaxima.end(); ++i ){
    int stripmax = *i;

    double ADCmax = fADCsums[hitindex[stripmax]];
    double prominence = ADCmax;

    int striplo = stripmax, striphi = stripmax;
    double ADCminright=ADCmax, ADCminleft=ADCmax;

    bool higherpeakright=false,higherpeakleft=false;
    int peakright = -1, peakleft = -1;

    while( striplist.find( striphi+1 ) != striplist.end() ){
      striphi++;

      Double_t ADCtest = fADCsums[hitindex[striphi]];
      
      if( ADCtest < ADCminright && !higherpeakright ){ //as long as we haven't yet found a higher peak to the right, this is the lowest point between the current maximum and the next higher peak to the right:
	ADCminright = ADCtest;
      }

      if( islocalmax[striphi] && ADCtest > ADCmax ){ //then this peak is in a contiguous group with another higher peak to the right:
	higherpeakright = true;
	peakright = striphi;
      }
    }

    while( striplist.find(striplo-1) != striplist.end() ){
      striplo--;
      Double_t ADCtest = fADCsums[hitindex[striplo]];
      if( ADCtest < ADCminleft && !higherpeakleft ){ //as long as we haven't yet found a higher peak to the left, this is the lowest point between the current maximum and the next higher peak to the left:
	ADCminleft = ADCtest;
      }

      if( islocalmax[striplo] && ADCtest > ADCmax ){ //then this peak is in a contiguous group with another higher peak to the left:
	higherpeakleft = true;
	peakleft = striplo;
      }
    }

    double sigma_sum = sqrt( double(fN_MPD_TIME_SAMP) )*fZeroSuppressRMS; //~25-50 ADC

    bool peak_close = false;
    if( !higherpeakleft ) ADCminleft = 0.0;
    if( !higherpeakright ) ADCminright = 0.0;

    if( higherpeakright || higherpeakleft ){ //this peak is contiguous with higher peaks on either the left or right or both:
      prominence = ADCmax - std::max( ADCminleft, ADCminright );

      if( higherpeakleft && std::abs( peakleft - stripmax ) <= 2*maxsep ) peak_close = true;
      if( higherpeakright && std::abs( peakright - stripmax ) <= 2*maxsep ) peak_close = true;

      if( peak_close && (prominence < fThreshold_2ndMax_nsigma * sigma_sum ||
			 prominance/ADCmax < fThresh_2ndMax_fraction ) ){
	peakstoerase.push_back( stripmax );
      }
    }
  }

  //Erase "insignificant" peaks (those in contiguous grouping with higher peak with prominence below thresholds):
  for( int ipeak=0; ipeak<peakstoerase.size(); ipeak++ ){
    localmaxima.erase( peakstoerase[ipeak] );
    islocalmax[peakstoerase[ipeak]] = false;
  }

  //Cluster formation and cluster splitting from remaining local maxima:
  for( std::set<int>::iterator i = localmax.begin(); i != localmax.end(); ++i ){
    int stripmax = *i;
    int striplo = stripmax;
    int striphi = stripmax;
    double ADCmax = fADCsums[hitindex[stripmax]];

    while( striplist.find( striplo-1 ) != striplist.end() &&
	   stripmax - striplo < maxsep ){
      striplo--;
    }

    while( striplist.find( striphi+1 ) != striplist.end() &&
	   striphi - stripmax < maxsep ){
      striphi++;
    }

    int nstrips = striphi-striplo+1;

    double sumx = 0.0, sumx2 = 0.0, sumADC = 0.0, sumt = 0.0, sumt2 = 0.0;
    double sumwx = 0.0;

    map<int,double> splitfraction;
    vector<double> stripADCsum(nstrips);

    for( int istrip=striplo; istrip<=striphi; istrip++ ){
      int nmax_strip = 1;
      double sumweight = ADCmax/(1.0 + pow( (stripmax-istrip)*pitch/fSigma_hitshape, 2 ) );
      double maxweight = sumweight;
      for( int jstrip=istrip-maxsep; jstrip<=istrip+maxsep; jstrip++ ){
	if( localmaxima.find( jstrip ) != localmaxima.end() && jstrip != stripmax ){
	  sumweight += fADCsums[hitindex[jstrip]]/( 1.0 + pow( (jstrip-istrip)*pitch/fSigma_hitshape, 2 ) );
	}
      }
   
      splitfraction[istrip] = maxweight/sumweight;

      double hitpos = (istrip + 0.5 - 0.5*Nstrips) * pitch; //local hit position along direction measured by these strips
      double ADCstrip = fADCsums[hitindex[istrip]] * splitfraction[istrip];
      double tstrip = fTmean[hitindex[istrip]];

      
      
    
  }
  
}

void    SBSGEMModule::Print( Option_t* opt) const{
  return;
}

Int_t   SBSGEMModule::Begin( THaRunBase* r){
  return 0;
}

Int_t   SBSGEMModule::End( THaRunBase* r){
  return 0;
}

