#include <iostream>

#include "SBSGEMPlane.h"
#include "TDatime.h"
#include "THaDetMap.h"
#include "THaEvData.h"

const int APVMAP[128] = {1, 33, 65, 97, 9, 41, 73, 105, 17, 49, 81, 113, 25, 57, 89, 121, 3, 35, 67, 99, 11, 43, 75, 107, 19, 51, 83, 115, 27, 59, 91, 123, 5, 37, 69, 101, 13, 45, 77, 109, 21, 53, 85, 117, 29, 61, 93, 125, 7, 39, 71, 103, 15, 47, 79, 111, 23, 55, 87, 119, 31, 63, 95, 127, 0, 32, 64, 96, 8, 40, 72, 104, 16, 48, 80, 112, 24, 56, 88, 120, 2, 34, 66, 98, 10, 42, 74, 106, 18, 50, 82, 114, 26, 58, 90, 122, 4, 36, 68, 100, 12, 44, 76, 108, 20, 52, 84, 116, 28, 60, 92, 124, 6, 38, 70, 102, 14, 46, 78, 110, 22, 54, 86, 118, 30, 62, 94, 126};

SBSGEMPlane::SBSGEMPlane( const char *name, const char *description,
			  THaDetectorBase* parent, bool ismc ):
    THaSubDetector(name,description,parent),
    fNch(0),fStrip(NULL),fPedestal(NULL),fIsMC(ismc)
{
    // FIXME:  To database
    fZeroSuppress    = kFALSE;
    fZeroSuppressRMS = 5.0;

        for( Int_t i = 0; i < N_MPD_TIME_SAMP; i++ ){
            fadc[i] = NULL;
        }

        return;
}

SBSGEMPlane::~SBSGEMPlane() {
    if( fStrip ){
        fadc0 = NULL;
        fadc1 = NULL;
        fadc2 = NULL;
        fadc3 = NULL;
        fadc4 = NULL;
        fadc5 = NULL;

        for( Int_t i = 0; i < N_MPD_TIME_SAMP; i++ ){
            delete fadc[i];
            fadc[i] = NULL;
        }
        delete fPedestal;
        fPedestal = NULL;
        delete fStrip;
        fStrip = NULL;
    }

    return;
}

Int_t SBSGEMPlane::ReadDatabase( const TDatime& date ){
  std::cout << "[SBSGEMPlane::ReadDatabase " << fName << "]" << std::endl;

    Int_t status;

    static const char* const here = "ReadDatabase";
    FILE* file = OpenFile( date );
    if( !file ) return kFileError;

    std::vector<Double_t> rawped;
    std::vector<Double_t> rawrms;

    //#ifdef MCDATA
    if(fIsMC){
      const DBRequest request[] = {
	{ "nchan",            &fNch,        kInt, 0, 1},
	{ "detmap",        &fChanMapData,        kIntV, 0, 0},
	{ "ped",            &rawped,        kDoubleV, 0, 1},
	{ "rms",            &rawrms,        kDoubleV, 0, 1},
	{}
      };
      status = LoadDB( file, date, request, fPrefix );
      fclose(file);
    
      //std::cout << fPrefix << " " << fNch << " " << fChanMapData.size() << " " << rawped.size() << " " << rawrms.size() << std::endl;
    
      UInt_t flags = 0;
    
      if( FillDetMap( fChanMapData, flags, here ) <= 0 )
	return kInitError;
    
      fStrip    = new Int_t [fNch];
    
      for( Int_t i = 0; i < N_MPD_TIME_SAMP; i++ ){
        fadc[i] = new Int_t [fNch];
	//std::cout << i << " " << fNch << std::endl;
        for( Int_t j = 0; j < fNch; j++ ){
	  fadc[i][j] = 0.0;
        }
      }
      fadc0 = fadc[0];
      fadc1 = fadc[1];
      fadc2 = fadc[2];
      fadc3 = fadc[3];
      fadc4 = fadc[4];
      fadc5 = fadc[5];
    
      fPedestal = new Double_t [fNch];
      fRMS      = new Double_t [fNch];

      for( Int_t i = 0; i < fNch; i++ ){
	fPedestal[i] = 0.0;
	fRMS[i] = 0.0;
      }
    
      if(rawped.size()>=fNch/128){
	for( Int_t i = 0; i < fNch; i++ )fPedestal[i] = rawped[i/128];
      }else{
	for( Int_t i = 0; i < fNch; i++ )fPedestal[i] = rawped[0];
      }

      if(rawrms.size()>=fNch/128){
	for( Int_t i = 0; i < fNch; i++ )fRMS[i] = rawrms[i/128];
      }else{
	for( Int_t i = 0; i < fNch; i++ )fRMS[i] = rawrms[0];
      }
    }else{
      //#endif
      const DBRequest request[] = {
        { "chanmap",        &fChanMapData,        kIntV, 0, 0},
        { "ped",            &rawped,        kDoubleV, 0, 1},
        { "rms",            &rawrms,        kDoubleV, 0, 1},
        {}
      };
      status = LoadDB( file, date, request, fPrefix );
      fclose(file);

      Int_t nentry = fChanMapData.size()/MPDMAP_ROW_SIZE;
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
        fMPDmap.push_back(thisdata);
      }

      std::cout << fName << " mapped to " << nentry << " APV25 chips" << std::endl;

      // FIXME:  make sure to delete if already initialized
      fStrip    = new Int_t [fNch];


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
		
	  std::cout << "[SBSGEMPlane::ReadDatabase]  WARNING: " << " strip " << idx  << " listed but not enough strips in cratemap" << std::endl;
	}
      }

      for( UInt_t i = 0; i < rawrms.size(); i++ ){
        if( (i % 2) == 1 ) continue;
        int idx = (int) rawrms[i];
	if( idx < N_APV25_CHAN*nentry ){
	  fRMS[idx] = rawrms[i+1];
	} else {
	  std::cout << "[SBSGEMPlane::ReadDatabase]  WARNING: " << " strip " << idx  << " listed but not enough strips in cratemap" << std::endl;
	}
      }
      //#ifdef MCDATA
    }
    //#endif
    return 0;
}

Int_t SBSGEMPlane::DefineVariables( EMode mode ) {
    if( mode == kDefine and fIsSetup ) return kOK;
      fIsSetup = ( mode == kDefine );

      RVarDef vars[] = {
          { "nch",   "Number of channels",   "fNchEff" },
          { "strip", "Strip number mapping", "fStrip_" },
          { "adc0", "ADC sample 0", "fadc_0" },
          { "adc1", "ADC sample 1", "fadc_1" },
          { "adc2", "ADC sample 2", "fadc_2" },
          { "adc3", "ADC sample 3", "fadc_3" },
          { "adc4", "ADC sample 4", "fadc_4" },
          { "adc5", "ADC sample 5", "fadc_5" },
	  /*
          { "strip", "Strip number mapping", "fStrip" },
          { "adc0", "ADC sample 0", "fadc0" },
          { "adc1", "ADC sample 1", "fadc1" },
          { "adc2", "ADC sample 2", "fadc2" },
          { "adc3", "ADC sample 3", "fadc3" },
          { "adc4", "ADC sample 4", "fadc4" },
          { "adc5", "ADC sample 5", "fadc5" },
	  */
          { 0 },
      };


      Int_t ret = DefineVarsFromList( vars, mode );

      if( ret != kOK )
          return ret;

      return kOK;

}

void    SBSGEMPlane::Clear( Option_t* opt){
  fNchEff = 0;
  memset(fStrip, 0, fNch*sizeof(Int_t));
  for(int i = 0; i<N_MPD_TIME_SAMP; i++){
    memset(fadc[i], 0, fNch*sizeof(Int_t));
  }
  memset(fadc0, 0, fNch*sizeof(Int_t));
  memset(fadc1, 0, fNch*sizeof(Int_t));
  memset(fadc2, 0, fNch*sizeof(Int_t));
  memset(fadc3, 0, fNch*sizeof(Int_t));
  memset(fadc4, 0, fNch*sizeof(Int_t));
  memset(fadc5, 0, fNch*sizeof(Int_t));
  
  fStrip_.clear();
  fadc_0.clear();
  fadc_1.clear();
  fadc_2.clear();
  fadc_3.clear();
  fadc_4.clear();
  fadc_5.clear();
  
  THaSubDetector::Clear(opt);
}

Int_t   SBSGEMPlane::Decode( const THaEvData& evdata ){
  //std::cout << "[SBSGEMPlane::Decode " << fName << "]" << std::endl;
    //#ifdef MCDATA
    if(fIsMC){
      UInt_t nmodules = fDetMap->GetSize();
      //std::cout << "nmodules " << nmodules << std::endl;
      int strip0 = 0;
      for( Int_t i = 0; i < nmodules; i++ ) {
	THaDetMap::Module* d = fDetMap->GetModule( i );
	//if(evdata.GetNumChan( d->crate, d->slot )>0)std::cout << fName << " " << d->crate << " " << d->slot << " " << evdata.GetNumChan( d->crate, d->slot ) << std::endl;
	for( UInt_t j = 0; j < evdata.GetNumChan( d->crate, d->slot ); j++) {
	
    UInt_t chan = evdata.GetNextChan( d->crate, d->slot, j );
    //std::cout << j << " chan " << chan << " first " << d->first << " lo " << d->lo << " hi " << d->hi << std::endl;
	  if( chan > d->hi || chan < d->lo ) continue;    // Not one of my channels.
	  int strip = strip0 + chan-d->lo;
	  assert(strip<fNch);
	
	  UInt_t nsamps = evdata.GetNumHits(d->crate, d->slot, chan);
	  if(nsamps!=N_MPD_TIME_SAMP)continue;
	  //std::cout << nsamps << " strip: " << strip << std::endl;
	  
	  fStrip_.push_back(strip);
	  fadc_0.push_back( evdata.GetData(d->crate, d->slot, chan, 0) - fPedestal[strip] );
	  fadc_1.push_back( evdata.GetData(d->crate, d->slot, chan, 1) - fPedestal[strip] );
	  fadc_2.push_back( evdata.GetData(d->crate, d->slot, chan, 2) - fPedestal[strip] );
	  fadc_3.push_back( evdata.GetData(d->crate, d->slot, chan, 3) - fPedestal[strip] );
	  fadc_4.push_back( evdata.GetData(d->crate, d->slot, chan, 4) - fPedestal[strip] );
	  fadc_5.push_back( evdata.GetData(d->crate, d->slot, chan, 5) - fPedestal[strip] );

	  fNchEff++;
	  fadc0[strip] = evdata.GetData(d->crate, d->slot, chan, 0) - fPedestal[strip];
	  fadc1[strip] = evdata.GetData(d->crate, d->slot, chan, 1) - fPedestal[strip];
	  fadc2[strip] = evdata.GetData(d->crate, d->slot, chan, 2) - fPedestal[strip];
	  fadc3[strip] = evdata.GetData(d->crate, d->slot, chan, 3) - fPedestal[strip];
	  fadc4[strip] = evdata.GetData(d->crate, d->slot, chan, 4) - fPedestal[strip];
	  fadc5[strip] = evdata.GetData(d->crate, d->slot, chan, 5) - fPedestal[strip];

	  // std::cout << fadc0[strip] << " " << fadc1[strip] << " " 
	  // 	    << fadc2[strip] << " " << fadc3[strip] << " " 
	  // 	    << fadc4[strip] << " " << fadc5[strip] << std::endl; 
	}// end loop on j
	strip0+= d->hi-d->lo;
      }//end loop on modules
    }else{
      //#endif
      
      fNch = 0;
      for (std::vector<mpdmap_t>::iterator it = fMPDmap.begin() ; it != fMPDmap.end(); ++it){
        Int_t effChan = it->mpd_id << 8 | it->adc_id;
        // Find channel for this crate/slot

        UInt_t nchan = evdata.GetNumChan( it->crate, it->slot );

	//        printf("nchan = %d\n", nchan );
	//if(nchan)std::cout << GetName() << " " << it->crate << " " << it->slot << " " << nchan << std::endl;
	  
	
        for( UInt_t ichan = 0; ichan < nchan; ++ichan ) {
          UInt_t chan = evdata.GetNextChan( it->crate, it->slot, ichan );
          if( chan != effChan ) continue; // not part of this detector


          UInt_t nsamp = evdata.GetNumHits( it->crate, it->slot, chan );
          assert(nsamp%N_MPD_TIME_SAMP==0);
          UInt_t nstrips = nsamp/N_MPD_TIME_SAMP;

	  // Loop over all the strips and samples in the data
	  Int_t isamp = 0;
	  for( UInt_t istrip = 0; istrip < nstrips; ++istrip ) {
	    assert(isamp<nsamp);
      UInt_t strip = evdata.GetRawData(it->crate, it->slot, chan, isamp);
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

	    fStrip[strip] = strip;
	    fNchEff++;
	    fStrip_.push_back(strip);
	    for(UInt_t adc_samp = 0; adc_samp < N_MPD_TIME_SAMP; adc_samp++) {
	      fadc[adc_samp][strip] =  evdata.GetData(it->crate, it->slot,
						     chan, isamp++) - fPedestal[strip];

	      assert( ((UInt_t) fNch) < fMPDmap.size()*N_APV25_CHAN );
	    }
	    assert(strip>=0); // Make sure we don't end up with negative strip numbers!
	    fadc_0.push_back(fadc[0][strip]);
	    fadc_1.push_back(fadc[1][strip]);
	    fadc_2.push_back(fadc[2][strip]);
	    fadc_3.push_back(fadc[3][strip]);
	    fadc_4.push_back(fadc[4][strip]);
	    fadc_5.push_back(fadc[5][strip]);
	    
	    // Zero suppression
	    if(!fZeroSuppress ||
	       ( fRMS[strip] > 0.0 && fabs(fadc[2][fNch])/
		 fRMS[strip] > fZeroSuppressRMS ) ){
	      fNch++;
	    }
	  }
        }

      }
      //#ifdef MCDATA
    }
    //#endif
//    std::cout << fName << " channels found  " << fNch << std::endl;

    return 0;
}

void    SBSGEMPlane::Print( Option_t* opt) const{
    return;
}

Int_t   SBSGEMPlane::Begin( THaRunBase* r){
    return 0;
}

Int_t   SBSGEMPlane::End( THaRunBase* r){
    return 0;
}

