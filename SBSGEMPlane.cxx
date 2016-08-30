#include <iostream>

#include "SBSGEMPlane.h"
#include "TDatime.h"
#include "THaEvData.h"


SBSGEMPlane::SBSGEMPlane( const char *name, const char *description,
    THaDetectorBase* parent ):
    THaSubDetector(name,description,parent),
    fNch(0),fStrip(NULL),fPedestal(NULL)
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
    std::cout << "[SBSGEMPlane::ReadDatabase]" << std::endl;

    Int_t status;

    FILE* file = OpenFile( date );
    if( !file ) return kFileError;

    std::vector<Double_t> rawped;
    std::vector<Double_t> rawrms;

    const DBRequest request[] = {
        { "chanmap",        &fChanMapData,        kIntV},
        { "ped",            &rawped,        kDoubleV},
        { "rms",            &rawrms,        kDoubleV},
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

    // Find max strip number
    int maxstrip = -1e9;
    for( UInt_t i = 0; i < rawped.size(); i++ ){
        // Odd values are pedestals themselves
        if( (i % 2) == 1 ) continue;
        if( rawped[i] > maxstrip ){
            maxstrip = rawped[i];
        }
    }

    fPedestal = new Double_t [maxstrip+1];
    fRMS      = new Double_t [maxstrip+1];
    for( Int_t i = 0; i < N_APV25_CHAN*nentry; i++ ){
        // FIXME needs to read in pedestal map
        fPedestal[i] = 0.0;
        fRMS[i] = 0.0;
    }

    for( UInt_t i = 0; i < rawped.size(); i++ ){
        if( (i % 2) == 1 ) continue;
        int idx = (int) rawped[i];
        fPedestal[idx] = rawped[i+1];
    }

    for( UInt_t i = 0; i < rawrms.size(); i++ ){
        if( (i % 2) == 1 ) continue;
        int idx = (int) rawrms[i];
        fRMS[idx] = rawrms[i+1];
    }


    return 0;
}

Int_t SBSGEMPlane::DefineVariables( EMode mode ) {
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

void    SBSGEMPlane::Clear( Option_t* opt){
    fNch = 0;
    return;
}

Int_t   SBSGEMPlane::Decode( const THaEvData& evdata ){
//    std::cout << "[SBSGEMPlane::Decode " << fName << "]" << std::endl;

    fNch = 0;
    for (std::vector<mpdmap_t>::iterator it = fMPDmap.begin() ; it != fMPDmap.end(); ++it){
        Int_t effChan = it->mpd_id << 8 | it->adc_id;
        // Find channel for this crate/slot

        Int_t nchan = evdata.GetNumChan( it->crate, it->slot );
        for( Int_t ichan = 0; ichan < nchan; ++ichan ) {
            Int_t chan = evdata.GetNextChan( it->crate, it->slot, ichan );
            if( chan != effChan ) continue; // not part of this detector


            Int_t nsamp = evdata.GetNumHits( it->crate, it->slot, chan );

//            std::cout << fName << " MPD " << it->mpd_id << " ADC " << it->adc_id << " found " << nsamp << std::endl;

//            std::cout << nsamp << " samples detected" << std::endl;

            assert( (nsamp/N_APV25_CHAN) == N_MPD_TIME_SAMP );
            for( Int_t strip = 0; strip < N_APV25_CHAN; ++strip ) {
                // data is packed like this
                // [ts1 of 128 chan] [ts2 of 128chan] ... [ts6 of 128chan]
                
                // Madness....   copy pasted from stand alone decoder
                // I bet there's a more elegant way to express this
                Int_t RstripNb= 32*(strip%4)+ 8*(int)(strip/4)- 31*(int)(strip/16);
                RstripNb=RstripNb+1+RstripNb%4-5*(((int)(RstripNb/4))%2);
                RstripNb=RstripNb+(127-2*RstripNb)*it->invert;

                Int_t RstripPos = RstripNb + 128*it->pos;

                /*
                if( it->adc_id == 10 ){
                std::cout << "ADC " << it->adc_id << " final strip pos: " << RstripPos << std::endl;
                }
                */

                fStrip[fNch] = RstripPos;



                for( Int_t adc_samp = 0; adc_samp < N_MPD_TIME_SAMP; adc_samp++ ){
                    int isamp = adc_samp*N_APV25_CHAN + strip;

                    assert(isamp < nsamp);

                    fadc[adc_samp][fNch] =  evdata.GetData(it->crate, it->slot, chan, isamp) -
                                            fPedestal[RstripPos];

                    assert( fNch < fMPDmap.size()*N_APV25_CHAN );
                }

                // Zero suppression
                if( !fZeroSuppress ||  
                      ( fRMS[RstripPos] > 0.0 && fabs(fadc[2][fNch])/fRMS[RstripPos] > fZeroSuppressRMS ) ){
                    fNch++;
                }
            }
        }

    }

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

