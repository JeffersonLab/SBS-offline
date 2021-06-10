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
  fDisableRefADC(false),fDisableRefTDC(false),
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
    { "model_in_detmap", &model_in_detmap,  kInt, 0, true }, ///< Does detector map have module numbers?
    { "chanmap",      &chanmap, kIntV,    0, true }, ///< Optional channel map
    { "start_chanmap",&fChanMapStart, kInt, 0, true}, ///< Optional start of channel numbering
    { "nrows",        &nrows,   kInt, 1, true }, ///< Number of rows in detector
    { "ncols",        &ncols,   kIntV, 0, false }, ///< [Optional] Number of columns in detector
    { "nlayers",       &nlayers,  kInt, 1, true }, ///< [Optional] Number of layers/divisions in each element of the detector
    { "angle",        &angle,   kFloat,  0, true },
    { "xyz",           &xyz,      kFloatV, 3 },  ///< If only 3 values specified, then assume as stating point for fist block and distribute according to dxyz
    { "dxdydz",         &dxyz,     kFloatV, 3, true },  ///< element spacing (dx,dy,dz)
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
  if( !chanmap.empty() ) {
    for(auto i : chanmap) {
      if (i < 0)
        nskipped++;
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
    int nrefchans = 0;
    // Find out how many channels in the detmap are reference channels
    for( Int_t i = 0; i < GetDetMap()->GetSize(); i++) {
      THaDetMap::Module *d = GetDetMap()->GetModule(i);
      Int_t nchan = d->GetNchan();
      if(d->refindex == -1) {
        nrefchans += nchan;
      }
    }
    nelem = fDetMap->GetTotNumChan() - nskipped - nrefchans; // Exclude skipped channels in count

    if( WithTDC() && WithADC() ) {
      if(nelem != 2*fNelem ) {
        Error( Here(here), "Number of crate module channels (%d) "
            "inconsistent with 2 channels per block (%d, expected)", nelem,
            fNelem );
        err = kInitError;
      }
    } else if ( nelem != fNelem) {
      Error( Here(here), "Number of crate module channels (%d) "
          "inconsistent with number of blocks (%d)", nelem, fNelem );
      err = kInitError;
    }
  }

  if(err)
    return err;

  if( !chanmap.empty() ) {
    // If a map is found in the database, ensure it has the correct size
    Int_t cmapsize = chanmap.size();
    if( WithTDC() && WithADC() ) {
      if(cmapsize - nskipped != 2*fNelem ) {
        Error( Here(here), "Number of logical channel to detector block map (%d) "
            "inconsistent with 2 channels per block (%d, expected)", cmapsize,
            2*fNelem );
        err = kInitError;
      }
    } else if ( cmapsize - nskipped != fNelem) {
      Error( Here(here), "Number of logical channel to detector block map (%d) "
          "inconsistent with number of blocks (%d)", cmapsize, fNelem );
      err = kInitError;
    }
  }
  if( !err ) {
    // Here, we will build our "local" channel map and check to make sure
    // we have the right number of adc channels, and when using TDCs, the
    // right number of TDC channels.
    // The map we are interested in is module channel to block number, where
    // the numbering of the blocks starts on the top left corner when standing
    // behind the detector and facing the target. We turn this into a row
    // and column, layer as appropriate.
    Int_t nmodules = GetDetMap()->GetSize();
    assert( nmodules > 0 );
    fChanMap.resize(nmodules);
    Decoder::THaCrateMap *cratemap = SBSManager::GetInstance()->GetCrateMap();
    Int_t ka = 0, kt = 0, km = 0, k = 0;
    for( Int_t i = 0; i < nmodules && !err; i++) {
      THaDetMap::Module *d = GetDetMap()->GetModule(i);
      // If the model number was not listed in the detmap section, fill it
      // in from the crate map
      if(!model_in_detmap) {
        d->SetModel(cratemap->getModel(d->crate,d->slot));
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
      if( d->refindex == -1) {
        for(Int_t chan = 0; chan < nchan; chan++) {
          SBSElement *el = MakeElement(0,0,0,0,0,0,-1);
          if(d->IsADC()) {
            el->SetADC(0.,1.);
          } else {
            el->SetTDC(0.,1.);
          }
          fRefElements.push_back(el);
        }
      } else if ( nchan > 0 ) {
        fChanMap[i].resize(nchan);
        for(Int_t chan = 0; chan < nchan && !err; chan++,k++) {
          if(!chanmap.empty() ) {
            // Skip disabled channels
            if(chanmap[k]<0) {
              fChanMap[i][chan] = -1;
              continue;
            }
            km = chanmap[k] - fChanMapStart;
          } else {
            km = d->IsADC() ? ka : kt;
          }
          // Count ADC and TDC channels separately
          if(d->IsADC()) {
            // Check the reference channel (if any) is valid
            if(!fDisableRefADC && d->refindex >=0 ) {
              // Make sure that if a reference channel/index was specified, there
              // have already been sufficient number of reference blocks created
              if(d->refindex >= int(fRefElements.size())) {
                Error( Here(here), "Cannot find reference module with index %d.",
                    d->refindex);
                err = kInitError;
                continue;
              } else if ( !fRefElements[d->refindex]->ADC() ) {
                Error( Here(here), "Error reference index %d is not an ADC!",
                    d->refindex);
                err = kInitError;
                continue;
              }
            }
            ka++;
          } else {
            // Check the reference channel (if any) is valid
            if(!fDisableRefTDC && d->refindex >=0 ) {
              // Make sure that if a reference channel/index was specified, there
              // have already been sufficient number of reference blocks created
              if(d->refindex >= int(fRefElements.size())) {
                Error( Here(here), "Cannot find reference module with index %d.",
                    d->refindex);
                err = kInitError;
                continue;
              } else if ( !fRefElements[d->refindex]->TDC() ) {
                Error( Here(here), "Error reference index %d is not an TDC!",
                    d->refindex);
                err = kInitError;
                continue;
              }
            }
            kt++;
          }
          assert( km < fNelem );
          assert( km >= 0);
          fChanMap[i][chan] = km;
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
  if(err)
    return err;

  // At this point, if an error has been encountered, don't bother continuing,
  // complain and return the error now.
  if(err)
    return err;

  // Read calibration parameters
  std::vector<Float_t> adc_ped, adc_gain, tdc_offset, tdc_cal;
  ResetVector(adc_ped,Float_t(0.0),fNelem);
  ResetVector(adc_gain,Float_t(1.0),fNelem);
  ResetVector(tdc_offset,Float_t(0.0),fNelem);
  ResetVector(tdc_cal,Float_t(1.0),fNelem);

  // Read adc pedestal and gains, and tdc offset and calibration
  // (should be organized by logical channel number, according to channel map)
  std::vector<DBRequest> vr;
  if(WithADC()) {
    vr.push_back({ "adc.pedestal", &adc_ped,    kFloatV, UInt_t(fNelem), 1 });
    vr.push_back({ "adc.gain",     &adc_gain,   kFloatV, UInt_t(fNelem), 1 });
  }
  if(WithTDC()) {
    vr.push_back({ "tdc.offset",   &tdc_offset, kFloatV, UInt_t(fNelem), 1 });
    vr.push_back({ "tdc.calib",    &tdc_cal,    kFloatV, UInt_t(fNelem), 1 });
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
  if(adc_ped.size() == 1) { // expand vector to specify calibration for all elements
    InitVector(adc_ped,adc_ped[0]);
  } else if ( adc_ped.size() != fNelem ) {
    Error( Here(here), "Inconsistent number of adc.ped specified. Expected "
        "%d but got %d",int(adc_ped.size()),fNelem);
    return kInitError;
  }

  // The user could have specified a single value for all calibrations,
  // in which case, we should propogate that

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
      rr = r+fChanMapStart;
      fElementGrid[r].resize(fNcols[r]);
      for(int c = 0; c < fNcols[r]; c++) {
        cc = c+fChanMapStart;
        for(int l = 0; l < fNlayers; l++, k++) {
          fElementGrid[r][c].resize(fNlayers);
          ll = l+fChanMapStart;
          //k = blkidx(r,c,l);
          x = xyz[0] - c*dxyz[0];
          y = xyz[1] - r*dxyz[1];
          z = xyz[2] - l*dxyz[2];
          SBSElement *e = MakeElement(x,y,z,rr,cc,ll,k);
          if( WithADC() ) {
            if( fModeADC == SBSModeADC::kWaveform ) {
              e->SetWaveform(adc_ped[k],adc_gain[k]);
            } else {
              e->SetADC(adc_ped[k],adc_gain[k]);
            }
          }
          if( WithTDC() ) { // TDC info
            e->SetTDC(tdc_offset[k],tdc_cal[k]);
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
    { "row", "Row for block in data vectors",  "fGood.row" },
    { "col", "Col for block in data vectors",  "fGood.col" },
    { "layer", "Layer for block in data vectors",  "fGood.layer" },
    { 0 }
  };

  Int_t err = DefineVarsFromList( vars, mode );
  if( err != kOK)
    return err;

  std::vector<RVarDef> ve;

  // Do we have an ADC? Then define ADC variables
  if(WithADC()) {
    // Register variables in global list
    ve.push_back( {"a","ADC integral", "fGood.a"} );
    if(fModeADC != SBSModeADC::kADCSimple) {
      ve.push_back( {"a_amp","ADC pulse amplitude", "fGood.a_amp"} );
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
    if(fModeTDC != SBSModeTDC::kTDCSimple) {
      // We have trailing edge and Time-Over-Threshold info to store
      ve.push_back({"tdc_te","Calibrated TDC trailing info","fGood.t_te"});
      ve.push_back({"tdc_tot","Time Over Threshold","fGood.t_ToT"});
    }

    if(fStoreRawHits) {
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

  SBSElement *blk = 0;
  // Loop over all modules in the calorimeter and decode accordingly
  for( UShort_t imod = 0; imod < fDetMap->GetSize(); imod++ ) {
    THaDetMap::Module *d = fDetMap->GetModule( imod );

    for(Int_t ihit = 0; ihit < evdata.GetNumChan( d->crate, d->slot ); ihit++) {
      fNhits++;

      // Get the next available channel, skipping the ones that do not belong
      // to our detector
      Int_t chan = evdata.GetNextChan( d->crate, d->slot, ihit );
      if( chan > d->hi || chan < d->lo || fChanMap[imod][chan - d->lo] == -1)
        continue;

      // Get the block index for this crate,slot,channel combo
      blk = fElements[ fChanMap[imod][chan - d->lo] ];
      if(d->IsADC()) {
        DecodeADC(evdata,blk,d,chan);
      } else if ( d->IsTDC()) {
        DecodeTDC(evdata,blk,d,chan);
      }
    }
  }

  return fNhits;
}

Int_t SBSGenericDetector::DecodeADC( const THaEvData& evdata,
    SBSElement *blk, THaDetMap::Module *d, Int_t chan)
{
  Int_t nhit = evdata.GetNumHits(d->crate, d->slot, chan);
  //std::cout << d->crate << " " << d->slot << " " << chan << std::endl;
  if(nhit <= 0  || !WithADC() || !blk)
    return 0;

  // TODO: Get the mode from the data stream (whenever that is implemented).
  // For now, we'll just use the user-specified fModeADC
  //Int_t mode = evdata.GetModule(d->crate, d->slot)->GetMode();
  
  if(fModeADC != SBSModeADC::kWaveform) {
    // Process all hits in this channel
    if(fModeADC == SBSModeADC::kADCSimple) { // Single ADC value (FADC250 mode 1)
      for(Int_t ihit = 0; ihit < nhit; ihit++) {
        blk->ADC()->Process( evdata.GetData(d->crate, d->slot, chan, ihit));
      }
    } else if (fModeADC == SBSModeADC::kADC) { // mode==7 in FADC250
      // here integral, time, peak, and pedestal are provided
      Float_t integral,time,peak,pedestal;
      Int_t lnhit = nhit/4; // Real number of hits
      for(Int_t ihit = 0; ihit < lnhit; ihit++) {
        integral = evdata.GetData(d->crate, d->slot, chan,           ihit);
        time     = evdata.GetData(d->crate, d->slot, chan,   lnhit + ihit);
        peak     = evdata.GetData(d->crate, d->slot, chan, 2*lnhit + ihit);
        pedestal = evdata.GetData(d->crate, d->slot, chan, 3*lnhit + ihit);
        blk->ADC()->Process(integral,time,peak,pedestal);
      }
    }
  } else {
    // As a first assumption, when using multi-valued ADCs then each "hit"
    // will correspond to the number of samples taken.
    std::vector<Float_t> samples;
    samples.resize(nhit);
    //std::cout << " nhit = " << nhit << ": ";  
    for(Int_t i = 0; i < nhit; i++) {
      samples[i] = evdata.GetData(d->crate, d->slot, chan, i);
      //std::cout << "  " << samples[i];
    }
    //std::cout << std::endl;
    //std::cout << blk << std::endl;
    //std::cout << blk->Waveform() << std::endl;
    blk->Waveform()->Process(samples);
    //std::cout << "ouh" << std::endl;
    samples.clear();
  }

  return nhit;
}

Int_t SBSGenericDetector::DecodeTDC( const THaEvData& evdata,
    SBSElement *blk, THaDetMap::Module *d, Int_t chan)
{
  if(!WithTDC() || !blk->TDC())
    return 0;

  Int_t nhit = evdata.GetNumHits(d->crate, d->slot, chan);
  Float_t refval  = 0;
  //std::cout << " crate " << d->crate << " slot " << d->slot << " chan " << chan << " nhit " << nhit << std::endl;
  if(!fDisableRefTDC && d->refindex>=0) {
    if(!fRefElements[d->refindex]->TDC()->HasData()) {
      std::cout << "Error reference TDC channel has no hits!" << std::endl;
    } else {
      // TODO: What should we do if reference channel has multiple hits?
      // For now, just take the first one
      refval = fRefElements[d->refindex]->TDC()->GetDataRaw(0);
    }
  }
  
  if(fIsMC)refval = 1000;
  
  Int_t edge = 0;
  for(Int_t ihit = 0; ihit < nhit; ihit++) {
    edge = 0; // Default is to not have any trailing info
    if(fModeTDC != SBSModeTDC::kTDCSimple) { // trailing edge info stored on raw data variable
      edge = evdata.GetRawData(d->crate, d->slot, chan, ihit);
      //std::cout << ihit << " " << evdata.GetData(d->crate, d->slot, chan, ihit) - refval << " " << edge << std::endl;
    }
    blk->TDC()->Process(
        evdata.GetData(d->crate, d->slot, chan, ihit) - refval, edge);
  }

  return nhit;
}

//_____________________________________________________________________________
void SBSGenericDetector::ClearEvent()
{
  // Call our version in case sub-classes have re-implemented it
  SBSGenericDetector::ClearOutputVariables();
  fNhits = 0;
  fCoarseProcessed = 0;
  fFineProcessed = 0;
  for(size_t k = 0; k < fElements.size(); k++) {
    fElements[k]->ClearEvent();
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
  for(Int_t k = 0; k < fNelem; k++) {
    blk = fElements[k];
    if(!blk)
      continue;

    blk->CoarseProcess();
    // If the above did not define the good hit, the sub-class is expected
    // to use re-implement the following function to find the good hit.
    FindGoodHit(blk);

    // Skip blocks that have no new data (unless allowed by the user)
    if(!blk->HasData() && !fStoreEmptyElements)
      continue;

    fGood.row.push_back(blk->GetRow());
    fGood.col.push_back(blk->GetCol());
    fGood.layer.push_back(blk->GetLayer());
    if(WithTDC() && blk->TDC()) {
      if(blk->TDC()->HasData()) {
        const SBSData::TDCHit &hit = blk->TDC()->GetGoodHit();
        fGood.t.push_back(hit.le.val);
        if(fModeTDC == SBSModeTDC::kTDC) { // has trailing info
          fGood.t_te.push_back(hit.te.val);
          fGood.t_ToT.push_back(hit.ToT.val);
        }
      } else if ( fStoreEmptyElements ) {
        fGood.t.push_back(0.0);
        if(fModeTDC == SBSModeTDC::kTDC) {
          fGood.t_te.push_back(0.0);
          fGood.t_ToT.push_back(0.0);
        }
      }
    }

    if(WithADC()) {
      if(fModeADC != SBSModeADC::kWaveform) {
        if(blk->ADC()->HasData()){
          const SBSData::PulseADCData &hit = blk->ADC()->GetGoodHit();
          fGood.a.push_back(hit.integral.val);
          if(fModeADC == SBSModeADC::kADC) { // Amplitude and time are also available
            fGood.a_amp.push_back(hit.amplitude.val);
            fGood.a_time.push_back(hit.time.val);
          }

          // Now store all the hits if specified the by user
          if(fStoreRawHits) {
            const std::vector<SBSData::PulseADCData> &hits = blk->ADC()->GetAllHits();
            for( const auto &hit : hits) {
              fRaw.a.push_back(hit.integral.val);
              fRaw.a_amp.push_back(hit.amplitude.val);
              fRaw.a_time.push_back(hit.time.val);
              // Do the same for the raw data
              //fRaw_raw.a.push_back(hit.integral.raw);
              //fRaw_raw.a_amp.push_back(hit.amplitude.raw);
              //fRaw_raw.a_time.push_back(hit.time.raw);
            }
          }
        } else if (fStoreEmptyElements) {
          fGood.a.push_back(0.0);
          if(fModeADC == SBSModeADC::kADC) {
            fGood.a_amp.push_back(0.0);
            fGood.a_time.push_back(0.0);
          }
        }
      } else { // Waveform mode
        SBSData::Waveform *wave = blk->Waveform();
        std::vector<Float_t> &s_r =wave->GetDataRaw();
        std::vector<Float_t> &s_c = wave->GetData();
        nsamples = s_r.size();
        idx = fGood.samps.size();
        fGood.sidx.push_back(idx);
        fGood.nsamps.push_back(nsamples);
        fGood.samps.resize(idx+nsamples);
        for(size_t s = 0; s < nsamples; s++) {
          fGood.samps[idx+s]   = s_c[s];
        }
        fGood.a.push_back(wave->GetIntegral().val);
        fGood.a_amp.push_back(wave->GetAmplitude().val);
        fGood.a_time.push_back(wave->GetTime().val);
      }
    }
  }

  fCoarseProcessed = 1;
  return 0;
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
}


///////////////////////////////////////////////////////////////////////////////
/// SBSGenericDetector constructor
SBSElement* SBSGenericDetector::MakeElement(Float_t x, Float_t y, Float_t z,
    Int_t row, Int_t col, Int_t layer, Int_t id)
{
  return new SBSElement(x,y,z,row,col,layer, id);
}
