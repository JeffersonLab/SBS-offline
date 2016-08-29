#include <iostream>

#include "SBSGEMPlane.h"
#include "TDatime.h"
#include "THaEvData.h"

#define MPDMAP_ROW_SIZE 8


SBSGEMPlane::SBSGEMPlane( const char *name, const char *description,
        THaDetectorBase* parent ):
    THaSubDetector(name,description,parent) {
        return;
}

SBSGEMPlane::~SBSGEMPlane() {
    return;
}

Int_t SBSGEMPlane::ReadDatabase( const TDatime& date ){
    std::cout << "[SBSGEMPlane::ReadDatabase]" << std::endl;

    Int_t status;

    FILE* file = OpenFile( date );
    if( !file ) return kFileError;

    const DBRequest request[] = {
        { "chanmap",        &fChanMapData,        kIntV},
        {}
    };
    status = LoadDB( file, date, request, fPrefix );

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

    return 0;
}


void    SBSGEMPlane::Clear( Option_t* opt){
    return;
}

Int_t   SBSGEMPlane::Decode( const THaEvData& evdata ){
    std::cout << "[SBSGEMPlane::Decode " << fName << "]" << std::endl;
    std::cout << evdata.GetNumChan( 11, 0 ) << std::endl;

    for (std::vector<mpdmap_t>::iterator it = fMPDmap.begin() ; it != fMPDmap.end(); ++it){
        Int_t effChan = it->mpd_id << 8 | it->adc_id;
        // Find channel for this crate/slot

        Int_t nchan = evdata.GetNumChan( it->crate, it->slot );
        for( Int_t ichan = 0; ichan < nchan; ++ichan ) {
            Int_t chan = evdata.GetNextChan( it->crate, it->slot, ichan );
            if( chan != effChan ) continue; // not part of this detector

            Int_t nsamp = evdata.GetNumHits( it->crate, it->slot, chan );

            for( Int_t isamp = 0; isamp < nsamp; ++isamp ) {
                std::cout << it->mpd_id << " " << it->adc_id << " "  << std::hex << evdata.GetData(it->crate, it->slot, chan, isamp) << std::endl;
                std::cout << std::dec;
            }
        }

    }

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

