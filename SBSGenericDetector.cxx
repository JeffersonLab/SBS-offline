////////////////////////////////////////////////////////////////////////////////
//
// SBSGenericDetector
//
//
////////////////////////////////////////////////////////////////////////////////

#include "SBSGenericDetector.h"

#include "THaEvData.h"
#include "THaDetMap.h"
#include "VarDef.h"
#include "VarType.h"
#include "THaTrack.h"
#include "TClonesArray.h"
#include "TDatime.h"
#include "TMath.h"
#include "SBSManager.h"
#include "THaCrateMap.h"

#include <cstring>
#include <iostream>
#include <iomanip>

ClassImp(SBSGenericDetector);

///////////////////////////////////////////////////////////////////////////////
/// SBSGenericDetector constructor
///
/// The default is to have single-valued ADC with no TDC information
/// Sub-classes can change this accordingly.
SBSGenericDetector::SBSGenericDetector( const char* name, const char* description,
    THaApparatus* apparatus ) :
  THaNonTrackingDetector(name,description,apparatus), fNrows(0),fNcolsMax(0),
  fNlayers(0), fModeADC(SBSModeADC::kADCSimple), fModeTDC(SBSModeTDC::kNone),
  fDisableRefADC(true),fDisableRefTDC(true),
  fConst(1.0), fSlope(0.0), fAccCharge(0.0), fStoreRawHits(false),
  fStoreEmptyElements(true), fIsMC(true)
{
  // Constructor.
  fCoarseProcessed = 0;
  fFineProcessed = 0;
}

///////////////////////////////////////////////////////////////////////////////
/// Default Destructor
SBSGenericDetector::~SBSGenericDetector()
{
  // Destructor. Removes internal arrays and global variables.

  if( fIsSetup )
    RemoveVariables();
  if( fIsInit ) {
    // What should be cleaned?
    for(Int_t i = 0; i < fNelem; i++) {
      delete fElements[i];
    }
  }

  ClearEvent();
}

///////////////////////////////////////////////////////////////////////////////
/// SetModeADC
void SBSGenericDetector::SetModeADC(SBSModeADC::Mode mode)
{
  fModeADC = mode;
  // Only the multi-function ADC is expected to have reference ADC
  SetDisableRefADC(fModeADC != SBSModeADC::kADC); 
}

///////////////////////////////////////////////////////////////////////////////
/// Read SBSGenericDetector Database
Int_t SBSGenericDetector::ReadDatabase( const TDatime& date )
{

  // Read this detector's parameters from the database file 'fi'.
  // This function is called by THaDetectorBase::Init() once at the
  // beginning of the analysis.
  // 'date' contains the date/time of the run being analyzed.

  // We can use this name here for logs
  static const char* const here = "ReadDatabase()";

  FILE* file = OpenFile( date );
  if( !file ) return kFileError;

  // Read in required geometry variables, which include fOrigin and fSize
  Int_t err = ReadGeometry( file, date, true );
  if( err ) {
    fclose(file);
    return err;
  }

  // Some temporary variables which we'll use to read in the database
  std::vector<Int_t> detmap, chanmap;
  std::vector<Float_t> xyz, dxyz;
  Float_t angle = 0.0;
  Int_t model_in_detmap = 0;

  Int_t nrows = 1, nlayers = 1;
  // In the following two optional parameters a pattern can be specified, where
  // if the number of entries corresponds to one or more row but less than
  // the total nnrows, the pattern will repeat until all nrows are described.
  // [Optional] specify variable number or columns per row, 1 entry per row
  // [Optional] specify starting offset for rows. Expects 3*n entries, where
  // 0 <= n <= nrows
  std::vector<Float_t> row_offset_pattern;
  // Specify number of columns per row.  If less than nrows entries provided
  // the pattern will repeat to fill up nrows entries
  std::vector<Int_t> ncols;
  
  bool is_mc;
  
  // Read mapping/geometry/configuration parameters
  fChanMapStart = 0;
  DBRequest config_request[] = {
    { "detmap",       &detmap,  kIntV }, ///< Detector map
   { "model_in_detmap", &model_in_detmap,  kInt, 0, true }, ///< Optional Does detector map have module numbers?
    { "chanmap",      &chanmap, kIntV,    0, true }, ///< Optional channel map
    { "start_chanmap",&fChanMapStart, kInt, 0, true}, ///< Optional start of channel numbering
    { "nrows",        &nrows,   kInt, 1, true }, ///< [Optional] Number of rows in detector
    { "ncols",        &ncols,   kIntV, 0, false }, ///< Number of columns in detector
    { "nlayers",       &nlayers,  kInt, 1, true }, ///< [Optional] Number of layers/divisions in each element of the detector
    { "angle",        &angle,   kFloat,  0, true },
    { "xyz",           &xyz,      kFloatV, 3 },  ///< If only 3 values specified, then assume as stating point for fist block and distribute according to dxyz
    { "dxdydz",         &dxyz,     kFloatV, 3},  ///< element spacing (dx,dy,dz)
    { "row_offset_pattern",        &row_offset_pattern,   kFloatV, 0, true }, ///< [Optional] conflicts with ncols
    { "is_mc",      &is_mc, kInt,    0, true }, ///< Optional channel map
    { 0 } ///< Request must end in a NULL
  };
  err = LoadDB( file, date, config_request, fPrefix );

  if(is_mc){// if this is simulated data, we do not care about the reference channel
    fIsMC = true;
    fDisableRefADC = true;
    fDisableRefTDC = true;
  }
  
  // Sanity checks (make sure there were no inconsistent values entered.
  if( !err && (nrows <= 0 || ncols.size() <= 0 || int(ncols.size()) > nrows 
        || nlayers <= 0) ) {
    Error( Here(here), "Illegal number of rows, columns and/or layers: %d %d %d"
        ". Must be > 0. Please fix the database.", nrows, int(ncols.size()), nlayers);
    return kInitError;
  }
  // Padd the ncols vector with a repeating pattern if not enough entries were
  // specified.
  Int_t ntemp = ncols.size();
  for(Int_t r = ntemp, i = 0; r < nrows; r++,i++) {
    if(ncols[i%ntemp]<=0) {
      Error( Here(here), "ncols cannot have negative entries!");
      return kInitError;
    }
    ncols.push_back(ncols[i%ntemp]);
  }

  // Padd the row_offset_pattern if not enough rows were specified.
  if(row_offset_pattern.size() > 0) {
    if(int(row_offset_pattern.size()) > 3*nrows || row_offset_pattern.size()%3 != 0) {
      Error( Here(here), "Inconsistent number of entries in row_offset_pattern "
          " specified.  Expected 3*nrows = %d but got %d",3*nrows,
          int(row_offset_pattern.size()));
      return kInitError;
    }

    ntemp = row_offset_pattern.size()/3;
    for(int r = ntemp, i = 0; r < nrows; r++,i++) {
      row_offset_pattern.push_back(row_offset_pattern[(i*3  )%ntemp]);
      row_offset_pattern.push_back(row_offset_pattern[(i*3+1)%ntemp]);
      row_offset_pattern.push_back(row_offset_pattern[(i*3+2)%ntemp]);
    }
  }

  std::vector< std::vector<int> > vnelems;
  vnelems.resize(nrows);
  Int_t nelem = 0;
  for(int r = 0; r < nrows; r++) {
    vnelems[r].resize(ncols[r]);
    for(int c = 0; c < ncols[r]; c++) {
      vnelems[r][c] = nlayers;
      nelem+=nlayers;
    }
  }
  // Safety check, make sure we didn't somehow change number of entries
  assert(int(ncols.size()) == nrows);

    // Reinitialization only possible for same basic configuration
  if( !err ) {
    if( fIsInit) {
      if(nelem != fNelem || nrows != fNrows ||
          ncols.size() != fNcols.size() || nlayers != fNlayers ) {
        Error( Here(here), "Cannot re-initalize with different number of rows, "
            "cols or layers. nelem(%d vs %d), nrows( %d vs %d), ncols(%d vs %d),"
            " nlayers(%d vs %d). Detector not re-initialized.",
            fNelem, nelem, fNrows, nrows, int(fNcols.size()), int(ncols.size()),
            fNlayers, nlayers);
        err = kInitError;
      } else {
        for(int r = 0; r < nrows; r++) {
          if(fNcols[r] != ncols[r]) {
            Error( Here(here), "Cannot re-initalize with different number of "
                " columns ( %d != %d ) for row %d.",fNcols[r],ncols[r], r);
            return kInitError;
          }
        }
      }
    }
    fNelem   = nelem;
    fNrows   = nrows;
    fNcols.clear();
    for(int r = 0; r < nrows; r++) {
      fNcols.push_back(ncols[r]);
      if(ncols[r]>fNcolsMax)
        fNcolsMax = ncols[r];
    }
    fNlayers = nlayers;
  }

  if(err)
    return err;

  // Find out how many channels got skipped:
  int nskipped = 0;
    int nrefchans = 0;
  if( !chanmap.empty() ) {
    for(auto i : chanmap) {
      if (i == -1) nskipped++;
      if (i == -1000) nrefchans++;
    }
  }
  // Clear out the old channel map before reading a new one
  fChanMap.clear();
  Int_t detmap_flags = THaDetMap::kFillRefIndex; // Specify reference index/channel
  if(model_in_detmap) {
    detmap_flags |= THaDetMap::kFillModel;
  }
  if( FillDetMap(detmap, detmap_flags, here) <= 0 ) {
    err = kInitError;  // Error already printed by FillDetMap
  } else {
    nelem = fDetMap->GetTotNumChan() - nskipped - nrefchans ; // Exclude skipped channels in count

    if( WithTDC() && WithADC() ) {
      if(nelem != 2*fNelem ) {
        Error( Here(here), "Number of crate module channels (%d) "
            "inconsistent with 2 channels per block (%d, expected)", nelem,
            fNelem );
        err = kInitError;
      }
    } else if ( nelem != fNelem) {
      Error( Here(here), "Number of crate module channels (%d) "
	     "inconsistent with number of blocks (%d) nskipped (%d) nrefchans (%d) ", nelem, fNelem , nskipped,nrefchans);
      err = kInitError;
    }
  }

  if(err)
    return err;

  if( !chanmap.empty() ) {
    // If a map is found in the database, ensure it has the correct size
    Int_t cmapsize = chanmap.size();
    if( WithTDC() && WithADC() ) {
      if(cmapsize - nskipped- nrefchans != 2*fNelem ) {
        Error( Here(here), "Number of logical channel to detector block map (%d) "
            "inconsistent with 2 channels per block (%d, expected)", cmapsize,
            2*fNelem );
        err = kInitError;
      }
    } else if ( cmapsize - nskipped- nrefchans != fNelem) {
      Error( Here(here), "Number of logical channel to detector block map (%d) "
          "inconsistent with number of blocks (%d)", cmapsize, fNelem );
      err = kInitError;
    }
  }
  Int_t NRefTDCElem=0;
  Int_t NRefADCElem=0;
  std::vector<Int_t> RefMode; //< Reftime MODE tdc =0 , adc = 1
  if( !err ) {
    // Here, we will build our "local" channel map and check to make sure
    // we have the right number of adc channels, and when using TDCs, the
    // right number of TDC channels.
    // The map we are interested in is module channel to block number, where
    // the numbering of the blocks starts on the top left corner when standing
    // behind the detector and facing the target. We turn this into a row
    // and column, layer as appropriate.
    //
    //
    Int_t nmodules = GetDetMap()->GetSize();
    assert( nmodules > 0 );
    fChanMap.resize(nmodules);
    fModuleRefTimeFlag.resize(nmodules);
    fRefChanMap.resize(nmodules);
    fRefChanLo.resize(nmodules);
    fRefChanHi.resize(nmodules);
    Decoder::THaCrateMap *cratemap = SBSManager::GetInstance()->GetCrateMap();
    Int_t kr = 0,ka = 0, kt = 0, km = 0, k = 0;
    for( Int_t i = 0; i < nmodules && !err; i++) {
      Int_t krmod = 0;
      fModuleRefTimeFlag[i] = kFALSE;
      fRefChanLo[i]=0;
      fRefChanHi[i]=0;
      THaDetMap::Module *d = GetDetMap()->GetModule(i);
      // If the model number was not listed in the detmap section, fill it
      // in from the crate map
      if(!model_in_detmap) {
        d->SetModel(cratemap->getModel(d->crate,d->slot));
      }
      if( model_in_detmap) {
	if (d->GetModel() == 526) {
	  d->MakeTDC(); 
	} else {
          Error( Here(here), "Need to modify SBSGenericDetector to specify whether TDC or ADC for module %d.", i);   
	}
      }
      if(!d->IsADC() && !d->IsTDC()) {
        // An unknown module was specified, complain and exit
        Error( Here(here), "Unknown module specified for module %d.", i);
        err = kInitError;
        continue;
      }
      Int_t nchan = d->GetNchan();
      // If this module has the reference channels, just create instances
      // of reference elements that will contain their data.
      if ( nchan > 0 ) {
        fChanMap[i].resize(nchan);
        fRefChanMap[i].resize(nchan);
        for(Int_t chan = 0; chan < nchan && !err; chan++,k++) {
          if(!chanmap.empty() ) {
            // Skip disabled channels
   	    fRefChanMap[i][chan]=-1;
            if(chanmap[k]==-1) {
              fChanMap[i][chan] = -1;
              continue;
            }
 	    if(chanmap[k]==-1000) {
              fModuleRefTimeFlag[i] = kTRUE;
              if(d->IsADC()) {
	        RefMode.push_back(1);
		fDisableRefADC = kFALSE;
		NRefADCElem++;
              } else {
	        RefMode.push_back(0);
		fDisableRefTDC = kFALSE;
		NRefTDCElem++;
              }
              fChanMap[i][chan] = chanmap[k];
    	      fRefChanMap[i][chan]=kr;
  	      if ( krmod==0) fRefChanLo[i]=chan;
  	      if ( krmod>=0) fRefChanHi[i]=chan;
              kr++;
	      krmod++;
	    } else {
            km = chanmap[k] - fChanMapStart;
	    }
          } else {
            km = d->IsADC() ? ka : kt;
          }
          // Count ADC and TDC channels separately
          if(d->IsADC()) {
            if(chanmap[k]!=-1000) ka++;
          } else {
            if(chanmap[k]!=-1000) kt++;
          }
          assert( km < fNelem );
          assert( km >= 0);
          if(chanmap[k]!=-1000) fChanMap[i][chan] = km;
        }
      } else {
        Error( Here(here), "No channels defined for module %d.", i);
        fChanMap.clear();
        err = kInitError;
      }
    }
    if(WithADC() && ka != fNelem) {
        Error( Here(here), "Inconsistent ADC channels, found %d, expected %d.", ka,fNelem);
        return kInitError;
    }
    if(WithTDC() && kt != fNelem) {
        Error( Here(here), "Inconsistent TDC channels, found %d, expected %d.", kt,fNelem);
        return kInitError;
    }
  }

  // At this point, if an error has been encountered, don't bother continuing,
  // complain and return the error now.
  if(err)
    return err;
  //
  // ADC and TDC reference time parameters
  std::vector<Float_t> reftdc_offset, reftdc_cal, reftdc_GoodTimeCut;
  std::vector<Float_t> refadc_ped, refadc_gain, refadc_conv, refadc_thres;
  std::vector<Int_t> refadc_FixThresBin,refadc_NSB,refadc_NSA,refadc_NPedBin;

  // Read calibration parameters
  // Read adc pedestal and gains, and tdc offset and calibration
  // (should be organized by logical channel number, according to channel map)
  std::vector<Float_t>  tdc_offset, tdc_cal, tdc_GoodTimeCut;
  std::vector<Float_t> adc_ped, adc_gain, adc_conv, adc_thres;
  std::vector<Int_t> adc_FixThresBin,adc_NSB,adc_NSA,adc_NPedBin;
  std::vector<DBRequest> vr;
  if(WithADC()) {
    vr.push_back({ "adc.pedestal", &adc_ped,    kFloatV, 0, 1 });
    vr.push_back({ "adc.gain",     &adc_gain,   kFloatV, 0, 1 });
    vr.push_back({ "adc.conv",     &adc_conv,   kFloatV, 0, 1 });
    vr.push_back({ "adc.thres",     &adc_thres,   kFloatV, 0, 1 });
    vr.push_back({ "adc.FixThresBin",     &adc_FixThresBin,   kIntV, 0, 1 });
    vr.push_back({ "adc.NSB",     &adc_NSB,   kIntV, 0, 1 });
    vr.push_back({ "adc.NSA",     &adc_NSA,   kIntV, 0, 1 });
    vr.push_back({ "adc.NPedBin",     &adc_NPedBin,   kIntV,0, 1 });
    if (!fDisableRefADC) {
    vr.push_back({ "refadc.pedestal", &refadc_ped,    kFloatV, 0, 1 });
    vr.push_back({ "refadc.gain",     &refadc_gain,   kFloatV, 0, 1 });
    vr.push_back({ "refadc.conv",     &refadc_conv,   kFloatV, 0, 1 });
    vr.push_back({ "refadc.thres",     &refadc_thres,   kFloatV, 0, 1 });
    vr.push_back({ "refadc.FixThresBin",     &refadc_FixThresBin,   kIntV, 0, 1 });
    vr.push_back({ "refadc.NSB",     &refadc_NSB,   kIntV, 0, 1 });
    vr.push_back({ "refadc.NSA",     &refadc_NSA,   kIntV, 0, 1 });
    vr.push_back({ "refadc.NPedBin",     &refadc_NPedBin,   kIntV,0, 1 });
    }
  }
  if(WithTDC()) {
    vr.push_back({ "tdc.offset",   &tdc_offset, kFloatV, 0, 1 });
    vr.push_back({ "tdc.calib",    &tdc_cal,    kFloatV, 0, 1 });
    vr.push_back({ "tdc.GoodTimeCut",    &tdc_GoodTimeCut,    kFloatV, 0, 1 });
    if (!fDisableRefTDC) {
    vr.push_back({ "reftdc.offset",   &reftdc_offset, kFloatV, 0, 1 });
    vr.push_back({ "reftdc.calib",    &reftdc_cal,    kFloatV, 0, 1 });
    vr.push_back({ "reftdc.GoodTimeCut",    &reftdc_GoodTimeCut,    kFloatV, 0, 1 });
    }
  };
  vr.push_back({0});
  err = LoadDB( file, date, vr.data(), fPrefix );


  // We are done reading from the file, so we can safely close it now
  fclose(file);

  // Again, no need to continue on errors
  if( err )
    return err;

  // What does this do again?!?!
  DefineAxes( angle*TMath::DegToRad() );

  // Check that there were either only 1 calibratoin value specified per key
  // or fNelements
    if (!fDisableRefTDC) {
       if(reftdc_offset.size() == 0) { // set all offset to zero
         ResetVector(reftdc_offset,Float_t(0.0),NRefTDCElem);
  } else if(reftdc_offset.size() == 1) { // expand vector to specify calibration for all elements
    Float_t temp=reftdc_offset[0];
    ResetVector(reftdc_offset,temp,NRefTDCElem);    
  } else if ( reftdc_offset.size() != NRefTDCElem ) {
    Error( Here(here), "Inconsistent number of reftdc.offset  specified. Expected "
	   "%d but got %d",NRefTDCElem,int(reftdc_offset.size()));
    return kInitError;
  }
       //
       if(reftdc_GoodTimeCut.size() == 0) { // set all GoodTimeCut to zero
         ResetVector(reftdc_GoodTimeCut,Float_t(0.0),NRefTDCElem);
  } else if(reftdc_GoodTimeCut.size() == 1) { // expand vector to specify calibration for all elements
    Float_t temp=reftdc_GoodTimeCut[0];
    ResetVector(reftdc_GoodTimeCut,temp,NRefTDCElem);    
  } else if ( reftdc_GoodTimeCut.size() != NRefTDCElem ) {
    Error( Here(here), "Inconsistent number of reftdc.GoodTimeCut  specified. Expected "
	   "%d but got %d",NRefTDCElem,int(reftdc_GoodTimeCut.size()));
    return kInitError;
  }
       //
       if(reftdc_cal.size() == 0) { // set all cal to 0.1
         ResetVector(reftdc_cal,Float_t(0.1),NRefTDCElem);
  } else if(reftdc_cal.size() == 1) { // expand vector to specify calibration for all elements
    Float_t temp=reftdc_cal[0];
    ResetVector(reftdc_cal,temp,NRefTDCElem);    
  } else if ( reftdc_cal.size() != NRefTDCElem ) {
    Error( Here(here), "Inconsistent number of reftdc.cal specified. Expected "
	   "%d but got %d",NRefTDCElem,int(reftdc_cal.size()));
    return kInitError;
  }
    }
    //
    if (!fDisableRefADC) {
  if(refadc_ped.size() == 0) { // set all ped to zero
    ResetVector(refadc_ped,Float_t(0.0),NRefADCElem);
  } else if(refadc_ped.size() == 1) { // expand vector to specify calibration for all elements
    Float_t temp=refadc_ped[0];
    ResetVector(refadc_ped,temp,NRefADCElem);    
  } else if ( refadc_ped.size() != NRefADCElem ) {
    Error( Here(here), "Inconsistent number of adc.ped specified. Expected "
	   "%d but got %d",NRefADCElem,int(refadc_ped.size()));
    return kInitError;
  }

  if(refadc_gain.size() == 0) { // set all gain to 1
     ResetVector(refadc_gain,Float_t(1.0),NRefADCElem);
  } else if(refadc_gain.size() == 1) { // expand vector to specify calibration for all elements
    Float_t temp=refadc_gain[0];
    ResetVector(refadc_gain,temp,NRefADCElem);    
  } else if ( refadc_gain.size() != NRefADCElem ) {
    Error( Here(here), "Inconsistent number of adc.gain specified. Expected "
        "%d but got %d",int(refadc_gain.size()),NRefADCElem);
    return kInitError;
  }

  if(refadc_thres.size() == 0) { // expand vector to specify calibration for all elements
    ResetVector(refadc_thres,Float_t(1.0),NRefADCElem);    
  } else if(refadc_thres.size() == 1) { // expand vector to specify calibration for all elements
    Float_t temp=refadc_thres[0];
    ResetVector(refadc_thres,temp,NRefADCElem);    
    std::cout << "set all elements  thres = " << refadc_thres[0] << std::endl;
  } else if ( refadc_thres.size() != NRefADCElem ) {
    Error( Here(here), "Inconsistent number of refadc.thres specified. Expected "
        "%d but got %d",int(refadc_thres.size()),NRefADCElem);
    return kInitError;
  }

  if(refadc_conv.size() == 0) { // expand vector to specify calibration for all elements
    ResetVector(refadc_conv,Float_t(1.0),NRefADCElem);    
  } else if(refadc_conv.size() == 1) { // expand vector to specify calibration for all elements
    Float_t temp=refadc_conv[0];
    ResetVector(refadc_conv,temp,NRefADCElem);    
    std::cout << "set all elements  conv = " << refadc_conv[0] << std::endl;
  } else if ( refadc_conv.size() != NRefADCElem ) {
    Error( Here(here), "Inconsistent number of refadc.conv specified. Expected "
        "%d but got %d",int(refadc_conv.size()),NRefADCElem);
    return kInitError;
  }

  if(refadc_NSB.size() == 0) { // expand vector to specify calibration for all elements
    ResetVector(refadc_NSB,3,NRefADCElem);    
  } else if(refadc_NSB.size() == 1) { // expand vector to specify calibration for all elements
    Int_t temp=refadc_NSB[0];
    ResetVector(refadc_NSB,temp,NRefADCElem);    
    std::cout << "set all elements  NSB = " << refadc_NSB[0] << std::endl;
  } else if ( refadc_NSB.size() != NRefADCElem ) {
    Error( Here(here), "Inconsistent number of refadc.NSB specified. Expected "
        "%d but got %d",int(refadc_NSB.size()),NRefADCElem);
    return kInitError;
  }

  if(refadc_NSA.size() == 0) { // expand vector to specify calibration for all elements
    ResetVector(refadc_NSA,10,NRefADCElem);    
  } else if(refadc_NSA.size() == 1) { // expand vector to specify calibration for all elements
    Int_t temp=refadc_NSA[0];
    ResetVector(refadc_NSA,temp,NRefADCElem);    
    std::cout << "set all elements  NSA = " << refadc_NSA[0] << std::endl;
  } else if ( refadc_NSA.size() != NRefADCElem ) {
    Error( Here(here), "Inconsistent number of refadc.NSA specified. Expected "
        "%d but got %d",int(refadc_NSA.size()),NRefADCElem);
    return kInitError;
  }

  if(refadc_NPedBin.size() == 0) { // expand vector to specify calibration for all elements
    ResetVector(refadc_NPedBin,4,NRefADCElem);    
  } else if(refadc_NPedBin.size() == 1) { // expand vector to specify calibration for all elements
    Int_t temp=refadc_NPedBin[0];
    ResetVector(refadc_NPedBin,temp,NRefADCElem);    
    std::cout << "set all elements  NPedBin = " << refadc_NPedBin[0] << std::endl;
  } else if ( refadc_NPedBin.size() != NRefADCElem ) {
    Error( Here(here), "Inconsistent number of adc.NPedBin specified. Expected "
        "%d but got %d",int(refadc_NPedBin.size()),NRefADCElem);
    return kInitError;
  }

  if(refadc_FixThresBin.size() == 0) { // expand vector to specify calibration for all elements
    ResetVector(refadc_FixThresBin,10,NRefADCElem);    
  } else if(refadc_FixThresBin.size() == 1) { // expand vector to specify calibration for all elements
    Int_t temp=refadc_FixThresBin[0];
    ResetVector(refadc_FixThresBin,temp,NRefADCElem);    
    std::cout << "set all elements  FixThresBin = " << refadc_FixThresBin[0] << std::endl;
  } else if ( refadc_FixThresBin.size() != NRefADCElem ) {
    Error( Here(here), "Inconsistent number of adc.FixThresBin specified. Expected "
        "%d but got %d",int(refadc_FixThresBin.size()),NRefADCElem);
    return kInitError;
  }
    }
    //
    fNRefElem = NRefADCElem + NRefTDCElem;
    // Set the reference time elements
    if (RefMode.size() > 0) {
    fRefElements.clear();
    fRefElements.resize(fNRefElem);
    for (Int_t nr=0;nr<fNRefElem;nr++) {
      SBSElement *el = MakeElement(0,0,0,nr,0,0,nr);
      if (RefMode[nr] ==0) {
	el->SetTDC(reftdc_offset[nr],reftdc_cal[nr],reftdc_GoodTimeCut[nr]);
      } else {
            if( fModeADC == SBSModeADC::kWaveform ) {
              el->SetWaveform(refadc_ped[nr],refadc_gain[nr],refadc_conv[nr]);
	      SBSData::Waveform *wave = el->Waveform();
              wave->SetWaveformParam(refadc_thres[nr],refadc_FixThresBin[nr],refadc_NSB[nr],refadc_NSA[nr],refadc_NPedBin[nr]);
            } else {
              el->SetADC(refadc_ped[nr],refadc_gain[nr]);
            }
      }
          fRefElements[nr] = el;
    }
    }
    //
  if(tdc_offset.size() == 0) { // set all ped to zero
    ResetVector(tdc_offset,Float_t(0.0),fNelem);
  } else if(tdc_offset.size() == 1) { // expand vector to specify calibration for all elements
    Float_t temp=tdc_offset[0];
    ResetVector(tdc_offset,temp,fNelem);    
  } else if ( tdc_offset.size() != fNelem ) {
    Error( Here(here), "Inconsistent number of adc.ped specified. Expected "
	   "%d but got %d",fNelem,int(tdc_offset.size()));
    return kInitError;
  }

  if(tdc_GoodTimeCut.size() == 0) { // set all ped to zero
    ResetVector(tdc_GoodTimeCut,Float_t(0.0),fNelem);
  } else if(tdc_GoodTimeCut.size() == 1) { // expand vector to specify calibration for all elements
    Float_t temp=tdc_GoodTimeCut[0];
    ResetVector(tdc_GoodTimeCut,temp,fNelem);    
  } else if ( tdc_GoodTimeCut.size() != fNelem ) {
    Error( Here(here), "Inconsistent number of adc.ped specified. Expected "
	   "%d but got %d",fNelem,int(tdc_GoodTimeCut.size()));
    return kInitError;
  }

  if(tdc_cal.size() == 0) { // set all ped to zero
    ResetVector(tdc_cal,Float_t(0.1),fNelem);
  } else if(tdc_cal.size() == 1) { // expand vector to specify calibration for all elements
    Float_t temp=tdc_cal[0];
    ResetVector(tdc_cal,temp,fNelem);    
  } else if ( tdc_cal.size() != fNelem ) {
    Error( Here(here), "Inconsistent number of adc.ped specified. Expected "
	   "%d but got %d",fNelem,int(tdc_cal.size()));
    return kInitError;
  }

  if(adc_ped.size() == 0) { // set all ped to zero
    ResetVector(adc_ped,Float_t(0.0),fNelem);
  } else if(adc_ped.size() == 1) { // expand vector to specify calibration for all elements
    Float_t temp=adc_ped[0];
    ResetVector(adc_ped,temp,fNelem);    
  } else if ( adc_ped.size() != fNelem ) {
    Error( Here(here), "Inconsistent number of adc.ped specified. Expected "
	   "%d but got %d",fNelem,int(adc_ped.size()));
    return kInitError;
  }

  if(adc_gain.size() == 0) { // set all gain to 1
     ResetVector(adc_gain,Float_t(1.0),fNelem);
  } else if(adc_gain.size() == 1) { // expand vector to specify calibration for all elements
    Float_t temp=adc_gain[0];
    ResetVector(adc_gain,temp,fNelem);    
  } else if ( adc_gain.size() != fNelem ) {
    Error( Here(here), "Inconsistent number of adc.gain specified. Expected "
        "%d but got %d",int(adc_gain.size()),fNelem);
    return kInitError;
  }

  if(adc_thres.size() == 0) { // expand vector to specify calibration for all elements
    ResetVector(adc_thres,Float_t(1.0),fNelem);    
  } else if(adc_thres.size() == 1) { // expand vector to specify calibration for all elements
    Float_t temp=adc_thres[0];
    ResetVector(adc_thres,temp,fNelem);    
    std::cout << "set all elements  thres = " << adc_thres[0] << std::endl;
  } else if ( adc_thres.size() != fNelem ) {
    Error( Here(here), "Inconsistent number of adc.thres specified. Expected "
        "%d but got %d",int(adc_thres.size()),fNelem);
    return kInitError;
  }

  if(adc_conv.size() == 0) { // expand vector to specify calibration for all elements
    ResetVector(adc_conv,Float_t(1.0),fNelem);    
  } else if(adc_conv.size() == 1) { // expand vector to specify calibration for all elements
    Float_t temp=adc_conv[0];
    ResetVector(adc_conv,temp,fNelem);    
    std::cout << "set all elements  conv = " << adc_conv[0] << std::endl;
  } else if ( adc_conv.size() != fNelem ) {
    Error( Here(here), "Inconsistent number of adc.conv specified. Expected "
        "%d but got %d",int(adc_conv.size()),fNelem);
    return kInitError;
  }

  if(adc_NSB.size() == 0) { // expand vector to specify calibration for all elements
    ResetVector(adc_NSB,3,fNelem);    
  } else if(adc_NSB.size() == 1) { // expand vector to specify calibration for all elements
    Int_t temp=adc_NSB[0];
    ResetVector(adc_NSB,temp,fNelem);    
    std::cout << "set all elements  NSB = " << adc_NSB[0] << std::endl;
  } else if ( adc_NSB.size() != fNelem ) {
    Error( Here(here), "Inconsistent number of adc.NSB specified. Expected "
        "%d but got %d",int(adc_NSB.size()),fNelem);
    return kInitError;
  }

  if(adc_NSA.size() == 0) { // expand vector to specify calibration for all elements
    ResetVector(adc_NSA,10,fNelem);    
  } else if(adc_NSA.size() == 1) { // expand vector to specify calibration for all elements
    Int_t temp=adc_NSA[0];
    ResetVector(adc_NSA,temp,fNelem);    
    std::cout << "set all elements  NSA = " << adc_NSA[0] << std::endl;
  } else if ( adc_NSA.size() != fNelem ) {
    Error( Here(here), "Inconsistent number of adc.NSA specified. Expected "
        "%d but got %d",int(adc_NSA.size()),fNelem);
    return kInitError;
  }

  if(adc_NPedBin.size() == 0) { // expand vector to specify calibration for all elements
    ResetVector(adc_NPedBin,4,fNelem);    
  } else if(adc_NPedBin.size() == 1) { // expand vector to specify calibration for all elements
    Int_t temp=adc_NPedBin[0];
    ResetVector(adc_NPedBin,temp,fNelem);    
    std::cout << "set all elements  NPedBin = " << adc_NPedBin[0] << std::endl;
  } else if ( adc_NPedBin.size() != fNelem ) {
    Error( Here(here), "Inconsistent number of adc.NPedBin specified. Expected "
        "%d but got %d",int(adc_NPedBin.size()),fNelem);
    return kInitError;
  }

  if(adc_FixThresBin.size() == 0) { // expand vector to specify calibration for all elements
    ResetVector(adc_FixThresBin,10,fNelem);    
  } else if(adc_FixThresBin.size() == 1) { // expand vector to specify calibration for all elements
    Int_t temp=adc_FixThresBin[0];
    ResetVector(adc_FixThresBin,temp,fNelem);    
    std::cout << "set all elements  FixThresBin = " << adc_FixThresBin[0] << std::endl;
  } else if ( adc_FixThresBin.size() != fNelem ) {
    Error( Here(here), "Inconsistent number of adc.FixThresBin specified. Expected "
        "%d but got %d",int(adc_FixThresBin.size()),fNelem);
    return kInitError;
  }


  // Before finishing, prepare vectors that will hold variable output data
  if( !fIsInit ) {
    fElements.clear();
    fElements.resize(fNelem);
    Float_t x = 0;
    Float_t y = 0;
    Float_t z = 0;
    Int_t k = 0;
    // the next three variables are the row,col,layer number starting
    // at fChanMapStart
    int rr = 0;
    int cc = 0;
    int ll = 0;
    fElementGrid.resize(fNrows);
    for(int r = 0; r < fNrows; r++) {
      rr = r;
      fElementGrid[r].resize(fNcols[r]);
      for(int c = 0; c < fNcols[r]; c++) {
        cc = c;
        for(int l = 0; l < fNlayers; l++, k++) {
          fElementGrid[r][c].resize(fNlayers);
          ll = l;
          //k = blkidx(r,c,l);
          x = xyz[0] - c*dxyz[0];
          y = xyz[1] - r*dxyz[1];
          z = xyz[2] - l*dxyz[2];
          SBSElement *e = MakeElement(x,y,z,rr,cc,ll,k+fChanMapStart);
          if( WithADC() ) {
            if( fModeADC == SBSModeADC::kWaveform ) {
              e->SetWaveform(adc_ped[k],adc_gain[k],adc_conv[k]);
	      SBSData::Waveform *wave = e->Waveform();
              wave->SetWaveformParam(adc_thres[k],adc_FixThresBin[k],adc_NSB[k],adc_NSA[k],adc_NPedBin[k]);
            } else {
              e->SetADC(adc_ped[k],adc_gain[k]);
            }
          }
          if( WithTDC() ) { // TDC info
            e->SetTDC(tdc_offset[k],tdc_cal[k],tdc_GoodTimeCut[k]);
          }
          fElements[k] = e;
          fElementGrid[r][c][l] = e;
        }
      }
    }
  }

  // All is well that ends well
  return kOK;
}

//_____________________________________________________________________________
Int_t SBSGenericDetector::DefineVariables( EMode mode )
{
  // Initialize global variables

  if( mode == kDefine && fIsSetup ) return kOK;
  if( !( WithADC() || WithTDC() ) ) {
    Error( Here("DefineVariables"),
        "GenericDetector %s defined with no data payload.",GetName());
    return kInitError;
  }
  fIsSetup = ( mode == kDefine );

  // Most of these variables were used previously, and so I'll leave
  // them here to preserve old analysis macros people may have written.
  // This includes things like fE, fNblk, fE_c, etc...
  RVarDef vars[] = {
    { "nhits", "Nhits",  "fNhits" },
    { "nrefhits", "Number of reference time hits",  "fNRefhits" },
    { "ngoodhits", "NGoodhits",  "fNGoodhits" },
    { "row", "Row for block in data vectors",  "fGood.row" },
    { "col", "Col for block in data vectors",  "fGood.col" },
    { "elemID", "Element ID for block in data vectors",  "fGood.elemID" },
    { "layer", "Layer for block in data vectors",  "fGood.layer" },
    { 0 }
  };

  Int_t err = DefineVarsFromList( vars, mode );
  if( err != kOK)
    return err;

 

  std::vector<RVarDef> ve;

  // TDC Reference Time variables 
  if(WithTDC() && !fDisableRefTDC) {
    ve.push_back({ "Ref.tdc", "Ref Time Calibrated TDC value", "fRefGood.t" });
    ve.push_back({ "Ref.tdc_mult", "Ref Time # hits in channel", "fRefGood.t_mult" });
    if(fModeTDC != SBSModeTDC::kTDCSimple) {
      // We have trailing edge and Time-Over-Threshold info to store
      ve.push_back({"Ref.tdc_te","Ref Time Calibrated TDC trailing info","fRefGood.t_te"});
      ve.push_back({"Ref.tdc_tot","Ref Time  Time Over Threshold","fRefGood.t_ToT"});
    }
    if(fStoreRawHits) {
      ve.push_back({ "Ref.hits.elemID",   "Ref Time ALL index",  "fRefRaw.elemID" });
      ve.push_back({ "Ref.hits.t",   "Ref Time All TDC leading edge times",  "fRefRaw.t" });
      if(fModeTDC != SBSModeTDC::kTDCSimple) {
        ve.push_back({ "Ref.hits.t_te",   "Ref Time All TDC trailing edge times",  "fRefRaw.t_te" });
        ve.push_back({ "Ref.hits.t_tot",  "Ref Time All TDC Time-over-threshold",  "fRefRaw.t_ToT" });
      }
    }
  }
//
  // Do we have an ADC? Then define ADC variables
  if(WithADC()) {
    // Register variables in global list
    ve.push_back({ "ped", "Pedestal for block in data vectors",  "fGood.ped" }),
     ve.push_back( {"a","ADC integral", "fGood.a"} );
     ve.push_back( {"a_mult","ADC # hits in channel", "fGood.a_mult"} );
    ve.push_back( {"a_p","ADC integral - ped", "fGood.a_p"} );
    ve.push_back( {"a_c","(ADC integral - ped)*gain", "fGood.a_c"} );
    if(fModeADC != SBSModeADC::kADCSimple) {
      ve.push_back( {"a_amp","ADC pulse amplitude", "fGood.a_amp"} );
      ve.push_back( {"a_amp_p","ADC pulse amplitude -ped", "fGood.a_amp_p"} );
      ve.push_back( {"a_time","ADC pulse time", "fGood.a_time"} );
    }
    if(fStoreRawHits) {
      ve.push_back({ "hits.a",   "All ADC inntegrals",  "fRaw.a" });
      ve.push_back({ "hits.a_amp",   "All ADC amplitudes",  "fRaw.a_amp" });
      ve.push_back({ "hits.a_time",   "All ADC pulse times",  "fRaw.a_time" });
    }
  }

  // Are we using TDCs? If so, define variables for TDCs
  if(WithTDC()) {
    ve.push_back({ "tdc", "Calibrated TDC value", "fGood.t" });
    ve.push_back({ "tdc_mult", "TDC # of hits per channel", "fGood.t_mult" });
    if(fModeTDC != SBSModeTDC::kTDCSimple) {
      // We have trailing edge and Time-Over-Threshold info to store
      ve.push_back({"tdc_te","Calibrated TDC trailing info","fGood.t_te"});
      ve.push_back({"tdc_tot","Time Over Threshold","fGood.t_ToT"});
    }
    if(fStoreRawHits) {
      ve.push_back({ "hits.elemID",   "All TDC Element ID",  "fRaw.elemID" });
      ve.push_back({ "hits.t",   "All TDC leading edge times",  "fRaw.t" });
      if(fModeTDC != SBSModeTDC::kTDCSimple) {
        ve.push_back({ "hits.t_te",   "All TDC trailing edge times",  "fRaw.t_te" });
        ve.push_back({ "hits.t_tot",  "All TDC Time-over-threshold",  "fRaw.t_ToT" });
      }
    }
  }

  // Are we using multi-valued ADCs? Then define the samples variables
  if(fModeADC == SBSModeADC::kWaveform) {
    ve.push_back({ "samps_idx", "Index in samples vector for given row-col module",
        "fGood.sidx" });
    ve.push_back({ "nsamps" , "Number of samples for given row-col",
        "fGood.nsamps"});
    ve.push_back({ "samps", "Calibrated ADC samples",  "fGood.samps" });
    ve.push_back({ "samps_elemID", "Calibrated ADC samples",  "fGood.samps_elemID" });
  }

  ve.push_back({0}); // Needed to specify the end of list
  return DefineVarsFromList( ve.data(), mode );
}

//_____________________________________________________________________________
Int_t SBSGenericDetector::Decode( const THaEvData& evdata )
{
  // Decode data

  // Clear last event
  ClearEvent();
  //static const char* const here = "Decode()";
  // Loop over modules for the reference time
  SBSElement *blk = 0;
  if (!fDisableRefADC || !fDisableRefTDC) {
  for( UShort_t imod = 0; imod < fDetMap->GetSize(); imod++ ) {
    if (!fModuleRefTimeFlag[imod]) continue;
    THaDetMap::Module *d = fDetMap->GetModule( imod );
    for(Int_t ihit = 0; ihit < evdata.GetNumChan( d->crate, d->slot ); ihit++) {
      Int_t chan = evdata.GetNextChan( d->crate, d->slot, ihit );
      if( chan > fRefChanHi[imod] || chan < fRefChanLo[imod]||  fChanMap[imod][chan - d->lo] >= -1) continue;
	blk = fRefElements[ fRefChanMap[imod][chan]];
	fNRefhits++;
         if(d->IsADC()) {
	   DecodeADC(evdata,blk,d,chan,kTRUE);
         } else if ( d->IsTDC()) {
	    DecodeTDC(evdata,blk,d,chan,kTRUE);
         }
    }
  }
  }

  // Loop over all modules decode accordingly
  blk = 0;
  for( UShort_t imod = 0; imod < fDetMap->GetSize(); imod++ ) {
    THaDetMap::Module *d = fDetMap->GetModule( imod );
	for(Int_t ihit = 0; ihit < evdata.GetNumChan( d->crate, d->slot ); ihit++) {
      // Get the next available channel, skipping the ones that do not belong
      // to our detector
      Int_t chan = evdata.GetNextChan( d->crate, d->slot, ihit );
      if( chan > d->hi || chan < d->lo || fChanMap[imod][chan - d->lo] == -1 || fChanMap[imod][chan - d->lo] == -1000)
        continue;
       fNhits++;

      // Get the block index for this crate,slot,channel combo
      blk = fElements[ fChanMap[imod][chan - d->lo] ];
      if(d->IsADC()) {
        DecodeADC(evdata,blk,d,chan,kFALSE);
      } else if ( d->IsTDC()) {
        DecodeTDC(evdata,blk,d,chan,kFALSE);
      }
    }
  }

  return fNhits;
  }
////
Int_t SBSGenericDetector::DecodeADC( const THaEvData& evdata,
    SBSElement *blk, THaDetMap::Module *d, Int_t chan,Bool_t IsRef)
{
  Int_t nhit = evdata.GetNumHits(d->crate, d->slot, chan);
  if(nhit <= 0  || !WithADC() || !blk)    return 0;
  // If not a reference element then determine the reference time to use
  Int_t reftime=0; 
  if (!IsRef && !fDisableRefADC && d->refindex>=0) {
     SBSElement *refblk = fRefElements[d->refindex];
     if(fModeADC == SBSModeADC::kWaveform && refblk->ADC()->HasData()) {
        SBSData::Waveform *wave = refblk->Waveform();
	//Now only one pulse found per sample 
       reftime = wave->GetTime().val;
      refblk->ADC()->SetGoodHit(0);
    }
  }
  if(fModeADC != SBSModeADC::kWaveform) {
    // Process all hits in this channel
    if(fModeADC == SBSModeADC::kADCSimple) { // Single ADC value 
        blk->ADC()->Process( evdata.GetData(d->crate, d->slot, chan, 0));
    } else if (fModeADC == SBSModeADC::kADC) { // mode==7 in FADC250
      // here integral, time, peak, and pedestal are provided
      Float_t integral,time,peak,pedestal;
      Int_t lnhit = nhit/4; // Real number of hits
      for(Int_t ihit = 0; ihit < lnhit; ihit++) {
        integral = evdata.GetData(d->crate, d->slot, chan,           ihit);
        time     = 0.0625*evdata.GetData(d->crate, d->slot, chan,   lnhit + ihit);
        peak     = evdata.GetData(d->crate, d->slot, chan, 2*lnhit + ihit);
        pedestal = evdata.GetData(d->crate, d->slot, chan, 3*lnhit + ihit);
        blk->ADC()->Process(integral,time,peak,pedestal);
	blk->ADC()->SetPed(pedestal);
 	
      }
    }
  } else {
    std::vector<Float_t> samples;
    samples.resize(nhit);
    for(Int_t i = 0; i < nhit; i++) {
      samples[i] = evdata.GetData(d->crate, d->slot, chan, i);
    }
    //std::cout << std::endl;
    //std::cout << blk << std::endl;
    //std::cout << blk->Waveform() << std::endl;
    blk->Waveform()->Process(samples);
    //std::cout << "ouh" << std::endl;
    samples.clear();
    SBSData::Waveform *wave = blk->Waveform();
    
    wave->SetValTime(wave->GetTime().val- reftime);
  }
  return nhit;
}


Int_t SBSGenericDetector::DecodeTDC( const THaEvData& evdata,
    SBSElement *blk, THaDetMap::Module *d, Int_t chan,Bool_t IsRef)
{
  if(!WithTDC() || !blk->TDC())
    return 0;
  //
  Int_t nhit = evdata.GetNumHits(d->crate, d->slot, chan);
  Float_t reftime  = 0;
  //
  //
  if(!IsRef && !fDisableRefTDC && d->refindex>=0) {
     SBSElement *refblk = fRefElements[d->refindex];
    if(!refblk->TDC()->HasData()) {
            std::cout << "Error reference TDC channel has no hits! refindex = " << d->refindex << " num ref tot = " << fNRefhits << " size = " << fRefElements.size() << std::endl;
    } else {
       Int_t nhits = refblk->TDC()->GetNHits(); 
       Float_t MinDiff = 10000.;
       UInt_t HitIndex = 0;
       Float_t RefCent = refblk->TDC()->GetGoodTimeCut();
       for (UInt_t ih=0;ih<nhits;ih++) {
	 if (abs(refblk->TDC()->GetData(ih)-RefCent) < MinDiff) {
           HitIndex = ih;
	   MinDiff = abs(refblk->TDC()->GetData(ih)-RefCent);
	 }
       }      
       reftime = refblk->TDC()->GetDataRaw(HitIndex);
      refblk->TDC()->SetGoodHit(HitIndex);
    }
  }
  
  if(fIsMC)reftime = 1000;
  
  Int_t edge = 0;
  Int_t elemID=blk->GetID();
  for(Int_t ihit = 0; ihit < nhit; ihit++) {
    edge = 0; // Default is to not have any trailing info
    if(fModeTDC != SBSModeTDC::kTDCSimple) { // trailing edge info stored on raw data variable
      edge = evdata.GetRawData(d->crate, d->slot, chan, ihit);
      //std::cout << ihit << " " << evdata.GetData(d->crate, d->slot, chan, ihit) - refval << " " << edge << std::endl;
    }
    if (edge ==1 && ihit ==0) continue; // skip first hit if trailing edge
    if (fModeTDC != SBSModeTDC::kTDCSimple && edge ==0 && ihit == nhit-1)  continue; // skip last hit if leading edge
    blk->TDC()->Process(elemID,
        evdata.GetData(d->crate, d->slot, chan, ihit) - reftime, edge);
  }
  if (!blk->TDC()->HasData()) {
          Float_t val= evdata.GetData(d->crate, d->slot, chan, 0);
	  blk->TDC()->Process(elemID,val - reftime , edge);
	  /*
            if (nhit==1)  {	 
	  std::cout << "Only one hit in time but not LE ref index = "  << d->refindex << " nhits = " << nhit  << " val = " << val << " event num = " << evdata.GetEvNum() << std::endl;
       } else if (nhit==2) {
	  std::cout << "Only two hits in time but not LE ref index = "  << d->refindex << " nhits = " << nhit  << " val = " << val << " event num = " << evdata.GetEvNum()<< std::endl;	 
       } else {
	  /std::cout << "More than two hits in time but not LE ref index = "  << d->refindex << " nhits = " << nhit  << " val = " << val << " event num = " << evdata.GetEvNum()<< std::endl;
       }	 
	  */
       }

  return nhit;
}

//_____________________________________________________________________________
void SBSGenericDetector::ClearEvent()
{
  // Call our version in case sub-classes have re-implemented it
  SBSGenericDetector::ClearOutputVariables();
  fNhits = 0;
  fNRefhits = 0;
  fNGoodhits = 0;
  fCoarseProcessed = 0;
  fFineProcessed = 0;
  for(size_t k = 0; k < fElements.size(); k++) {
    fElements[k]->ClearEvent();
  }
  for(size_t k = 0; k < fRefElements.size(); k++) {
    fRefElements[k]->ClearEvent();
  }
}

//_____________________________________________________________________________
void SBSGenericDetector::Clear(Option_t* opt)
{
  // Call our version in case sub-classes have re-implemented it
  SBSGenericDetector::ClearOutputVariables();
  fNhits = 0;
  fNRefhits = 0;
  fNGoodhits = 0;
  fCoarseProcessed = 0;
  fFineProcessed = 0;
  for(size_t k = 0; k < fElements.size(); k++) {
    fElements[k]->ClearEvent();
  }
  for(size_t k = 0; k < fRefElements.size(); k++) {
    fRefElements[k]->ClearEvent();
  }
}

Int_t SBSGenericDetector::CoarseProcess(TClonesArray& )// tracks)
{
  // Make sure we haven't already been called in this event
  if(fCoarseProcessed) return 0;

  // Pack simple data for output to the tree, and call CoarseProcess on all elements
  SBSElement *blk = 0;
  size_t nsamples;
  size_t idx;
  // Reference time
  for(Int_t k = 0; k < fNRefElem; k++) {
    blk = fRefElements[k];
    if(!blk) continue;
    fRefGood.row.push_back(blk->GetRow());
    fRefGood.col.push_back(blk->GetCol());
    fRefGood.layer.push_back(blk->GetLayer());
    fRefGood.elemID.push_back(blk->GetID());
    if(WithTDC() && blk->TDC()) {
      if ( blk->TDC()->HasData()) {
        const SBSData::TDCHit &hit = blk->TDC()->GetGoodHit();
        fRefGood.t.push_back(hit.le.val);
        fRefGood.t_mult.push_back(blk->TDC()->GetNHits());
        if(fModeTDC == SBSModeTDC::kTDC) { // has trailing info
          fRefGood.t_te.push_back(hit.te.val);
          fRefGood.t_ToT.push_back(hit.ToT.val);
        }
      } else if ( fStoreEmptyElements ) {
        fRefGood.t_mult.push_back(0);
        fRefGood.t.push_back(kBig);
        if(fModeTDC == SBSModeTDC::kTDC) {
          fRefGood.t_te.push_back(kBig);
          fRefGood.t_ToT.push_back(kBig);
        }
      }
        if(fStoreRawHits) {
            const std::vector<SBSData::TDCHit> &hits = blk->TDC()->GetAllHits();
            for( const auto &hit : hits) {
              fRefRaw.elemID.push_back(hit.elemID);
              fRefRaw.t.push_back(hit.le.val);
               if(fModeTDC == SBSModeTDC::kTDC) { // has trailing info
              fRefRaw.t_te.push_back(hit.te.val);
              fRefRaw.t_ToT.push_back(hit.ToT.val);
	       }
              }
	}
    }
  }

  
  
  for(Int_t k = 0; k < fNelem; k++) {
    blk = fElements[k];
    if(!blk)
      continue;
 
    // If the above did not define the good hit, the sub-class is expected
    // to use re-implement the following function to find the good hit.
    fNGoodhits=fNGoodhits+FindGoodHit(blk);

    // Skip blocks that have no new data (unless allowed by the user)
     if(!blk->HasData() && !fStoreEmptyElements)
      continue;

    fGood.row.push_back(blk->GetRow());
    fGood.col.push_back(blk->GetCol());
    fGood.layer.push_back(blk->GetLayer());
    fGood.elemID.push_back(blk->GetID());
    if(WithTDC() && blk->TDC()) {
      if(blk->TDC()->HasData() ) {
        const SBSData::TDCHit &hit = blk->TDC()->GetGoodHit();
        fGood.t.push_back(hit.le.val);
        fGood.t_mult.push_back(blk->TDC()->GetNHits());
        if(fModeTDC == SBSModeTDC::kTDC) { // has trailing info
          fGood.t_te.push_back(hit.te.val);
          fGood.t_ToT.push_back(hit.ToT.val);
        }
      } else if ( fStoreEmptyElements ) {
        fGood.t.push_back(kBig);
        fGood.t_mult.push_back(0);
        if(fModeTDC == SBSModeTDC::kTDC) {
          fGood.t_te.push_back(kBig);
          fGood.t_ToT.push_back(kBig);
        }
      }
      if(fStoreRawHits) {
            const std::vector<SBSData::TDCHit> &hits = blk->TDC()->GetAllHits();
            for( const auto &hit : hits) {
              fRaw.elemID.push_back(hit.elemID);
              fRaw.t.push_back(hit.le.val);
               if(fModeTDC == SBSModeTDC::kTDC) { // has trailing info
              fRaw.t_te.push_back(hit.te.val);
              fRaw.t_ToT.push_back(hit.ToT.val);
	       }
              }
      }
    }

    if(WithADC()) {
      if(fModeADC != SBSModeADC::kWaveform) {
          if(blk->ADC()->HasData() ){
	  if (blk->ADC()->GetGoodHitIndex() >=0) {
          Float_t ped=blk->ADC()->GetPed();
          fGood.ped.push_back(ped);
          const SBSData::PulseADCData &hit = blk->ADC()->GetGoodHit();
          fGood.a.push_back(hit.integral.raw);
          fGood.a_mult.push_back(blk->ADC()->GetNHits());
          fGood.a_p.push_back(hit.integral.raw-ped);
          fGood.a_c.push_back(hit.integral.val);
          if(fModeADC == SBSModeADC::kADC) { // Amplitude and time are also available
            fGood.a_amp.push_back(hit.amplitude.raw*2000./4096.);
            fGood.a_amp_p.push_back((hit.amplitude.val-ped)*2000./4096.);
            fGood.a_time.push_back(hit.time.val);
          }
	  }
          // Now store all the hits if specified the by user
          if(fStoreRawHits) {
            const std::vector<SBSData::PulseADCData> &hits = blk->ADC()->GetAllHits();
            for( const auto &hit : hits) {
              fRaw.a.push_back(hit.integral.val);
              fRaw.a_amp.push_back(hit.amplitude.val);
              fRaw.a_time.push_back(hit.time.val);
             }
          }
	  } else if (fStoreEmptyElements) {
          fGood.a.push_back(0.0);
          fGood.a_mult.push_back(0);
          fGood.a_p.push_back(0.0);
          fGood.a_c.push_back(0.0);
          if(fModeADC == SBSModeADC::kADC) {
            fGood.a_amp.push_back(0.0);
            fGood.a_amp_p.push_back(0.0);
            fGood.a_time.push_back(0.0);
          }
         }
      } else { // Waveform mode
        SBSData::Waveform *wave = blk->Waveform();
	if(wave->HasData()) {		
        std::vector<Float_t> &s_r =wave->GetDataRaw();
        std::vector<Float_t> &s_c = wave->GetData();
        nsamples = s_r.size();
        idx = fGood.samps.size();
        fGood.sidx.push_back(idx);
        fGood.samps_elemID.push_back(k);
        fGood.nsamps.push_back(nsamples);
        fGood.samps.resize(idx+nsamples);
        for(size_t s = 0; s < nsamples; s++) {
          fGood.samps[idx+s]   = s_c[s];
        }
        fGood.ped.push_back(wave->GetPed());
        fGood.a.push_back(wave->GetIntegral().raw);
        fGood.a_p.push_back(wave->GetIntegral().val);
        fGood.a_c.push_back(wave->GetIntegral().val);
        fGood.a_amp.push_back(wave->GetAmplitude().raw);
        fGood.a_amp_p.push_back(wave->GetAmplitude().val);
        fGood.a_time.push_back(wave->GetTime().val);
       } else if (fStoreEmptyElements) {
          fGood.a.push_back(0.0);
          fGood.a_mult.push_back(0);
          fGood.a_p.push_back(0.0);
          fGood.a_c.push_back(0.0);
            fGood.a_amp.push_back(0.0);
            fGood.a_amp_p.push_back(0.0);
            fGood.a_time.push_back(0.0);
	}
      }
    }
  }

  fCoarseProcessed = 1;
  return 0;
}

//
Int_t SBSGenericDetector::FindGoodHit(SBSElement *blk)
{
  Int_t GoodHit=0;  
  if (blk->TDC()) blk->TDC()->SetGoodHit(-1);
  if (blk->TDC()&& blk->HasData()) {
       Int_t nhits = blk->TDC()->GetNHits(); 
       Float_t MinDiff = 10000.;
       UInt_t HitIndex = 0;
       Float_t GoodTimeCut = blk->TDC()->GetGoodTimeCut();
       for (UInt_t ih=0;ih<nhits;ih++) {
	 if (abs(blk->TDC()->GetData(ih)-GoodTimeCut) < MinDiff) {
           HitIndex = ih;
	   MinDiff = abs(blk->TDC()->GetData(ih)-GoodTimeCut);
	 }
       }      
      blk->TDC()->SetGoodHit(HitIndex);
    GoodHit=1;
  }
  if (blk->ADC()&& blk->HasData()) {		
    blk->ADC()->SetGoodHit(0);
    GoodHit=1;
    
  }
  return GoodHit;
}

//_____________________________________________________________________________
Int_t SBSGenericDetector::FineProcess(TClonesArray&)//tracks)
{
  fFineProcessed = 1;
  return 0;
}

void SBSGenericDetector::ClearOutputVariables()
{
  fGood.clear();
  fRaw.clear();
  fRefGood.clear();
  fRefRaw.clear();
}


///////////////////////////////////////////////////////////////////////////////
/// SBSGenericDetector constructor
SBSElement* SBSGenericDetector::MakeElement(Float_t x, Float_t y, Float_t z,
    Int_t row, Int_t col, Int_t layer, Int_t id)
{
  return new SBSElement(x,y,z,row,col,layer, id);
}
