///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// SBSCalorimeter                                                               //
//                                                                           //
// Shower counter class, describing a generic segmented shower detector      //
// (preshower or shower).                                                    //
// Currently, only the "main" cluster, i.e. cluster with the largest energy  //
// deposition is considered. Units of measurements are MeV for energy of     //
// shower and centimeters for coordinates.                                   //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "SBSCalorimeter.h"

//#include "THaBBe.h"
//#include "THaGlobals.h"

#include "THaEvData.h"
#include "THaDetMap.h"
#include "VarDef.h"
#include "VarType.h"
#include "THaTrack.h"
#include "TClonesArray.h"
#include "TDatime.h"
#include "TMath.h"

#include <cstring>
#include <iostream>
#include <iomanip>
#define CLUSTER_BLOCK_RADIUS 1

ClassImp(SBSCalorimeter);

///////////////////////////////////////////////////////////////////////////////
/// SBSCalorimeter constructor
///
/// The default is to have single-valued ADC with no TDC information
/// Sub-classes can change this accordingly.
SBSCalorimeter::SBSCalorimeter( const char* name, const char* description,
    THaApparatus* apparatus ) :
  THaNonTrackingDetector(name,description,apparatus), fNrows(0), fNcols(0),
  fNlayers(0), fWithTDC(false), fWithADCSamples(false)
  //, fNChan(NULL), fChanMap(NULL)
{
  /*
  // Constructor.
  fCoarseProcessed = 0;
  fFineProcessed = 0;
  */
}

///////////////////////////////////////////////////////////////////////////////
/// Default Destructor
SBSCalorimeter::~SBSCalorimeter()
{
  // Destructor. Removes internal arrays and global variables.

  if( fIsSetup )
    RemoveVariables();
  if( fIsInit ) {
    // What should be cleaned?
    for(Int_t i = 0; i < fNelem; i++) {
      delete fBlocks[i];
    }
  }

  fDataOut.ClearEvent();
}

///////////////////////////////////////////////////////////////////////////////
/// Read SBSCalorimeter Database
Int_t SBSCalorimeter::ReadDatabase( const TDatime& date )
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
  std::vector<Double_t> xyz, dxyz;
  Int_t ncols = 1, nrows = 1, nlayers = 1;
  Double_t angle = 0.0;
  Double_t adc_ped_mult = 1;

  // Read mapping/geometry/configuration parameters
  DBRequest config_request[] = {
    { "detmap",       &detmap,  kIntV }, ///< Detector map
    { "chanmap",      &chanmap, kIntV,    0, true }, ///< Optional channel map
    { "chanmap_start",&fChanMapStart, kInt, 0, true}, ///< Optional start of channel numbering
    { "ncols",        &ncols,   kInt, 1, true }, ///< Number of columns in detector
    { "nrows",        &nrows,   kInt, 1, true }, ///< Number of rows in detector
    { "nlayers",       &nlayers,  kInt,1,true }, ///< [Optional] Number of layers/divisions in each module/block of the detector
    { "angle",        &angle,   kDouble,  0, true },
    { "xyz",           &xyz,      kDoubleV, 3 },  ///< center pos of block 1
    { "dxdydz",         &dxyz,     kDoubleV, 3 },  ///< block spacing (dx,dy,dz)
    //{ "emin",         &fEmin,   kDouble },
    { "adc_ped_mult", &adc_ped_mult,  kDouble, 0, true },
    //{ "sum_gain_mult",&sum_gain_mult,  kDouble, 0, true },
    { 0 } ///< Request must end in a NULL
  };
  err = LoadDB( file, date, config_request, fPrefix );

  // Sanity checks (make sure there were no inconsistent values entered.
  if( !err && (nrows <= 0 || ncols <= 0 || nlayers <= 0) ) {
    Error( Here(here), "Illegal number of rows, columns and/or layers: %d %d %d"
        ". Must be > 0. Please fix the database.", nrows, ncols, nlayers);
    err = kInitError;
  }

  Int_t nelem = ncols * nrows * nlayers;
  // Also compute the max possible cluster size (which at most should be 3x3)
  // TODO: Don't hard code a max of 3x3 cluster!
  Int_t nclbl = TMath::Min( 3, nrows ) * TMath::Min( 3, ncols );

  // Reinitialization only possible for same basic configuration
  if( !err ) {
    if( fIsInit && (nelem != fNelem || nrows != fNrows || ncols != fNcols ||
         nlayers != fNlayers ) ) {
      Error( Here(here), "Cannot re-initalize with different number of blocks or "
          "blocks per cluster (was: %d, now: %d). Detector not re-initialized.",
          fNelem, nelem );
      err = kInitError;
    } else {
      fNelem   = nelem;
      fNrows   = nrows;
      fNcols   = ncols;
      fNlayers = nlayers;
      fNclublk = nclbl;
    }
  }

  if(err)
    return err;

  // Clear out the old channel map before reading a new one
  fChanMap.clear();
  if( FillDetMap(detmap, 0, here) <= 0 ) {
    err = kInitError;  // Error already printed by FillDetMap
  } else {
    nelem = fDetMap->GetTotNumChan();
    if( fWithTDC ) {
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
    if( fWithTDC  ) {
      if(cmapsize != 2*fNelem ) {
        Error( Here(here), "Number of logical channel to detector block map (%d) "
            "inconsistent with 2 channels per block (%d, expected)", cmapsize,
            2*fNelem );
        err = kInitError;
      }
    } else if ( cmapsize != fNelem) {
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
    Bool_t makeADC = true;
    for( Int_t i = 0, k = 0; i < nmodules && !err; i++) {
      THaDetMap::Module *d = GetDetMap()->GetModule(i);
      Int_t nchan = d->GetNchan();
      if( nchan > 0 ) {
        fChanMap[i].resize(nchan);
        // To simplify finding out which channels are ADCs and which are TDCs
        // we'll just require that the user make the first fNelem channels
        // correspond to the TDCs, and the other fNelem to the TDCs (if in use).
        if(makeADC) {
          d->MakeADC();
        } else {
          d->MakeTDC();
        }
        for(Int_t chan = 0; chan < nchan; chan++) {
          assert( k < fNelem );
          fChanMap[i][chan] = chanmap.empty() ? k : chanmap[k] - fChanMapStart;
          k++;
        }
      } else {
        Error( Here(here), "No channels defined for module %d.", i);
        fChanMap.clear();
        err = kInitError;
      }
      // The check for TDC was done above already. So if we reach fNelem,
      // then it means that we are now labeling the TDC channels for remaining
      // modules.
      if(k>=fNelem) {
        k = 0;
        makeADC = false;
      }
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
  DBRequest calib_request[] = {
    { "adc.pedestal", &adc_ped,    kFloatV, UInt_t(fNelem), 1 },
    { "adc.gain",     &adc_gain,   kFloatV, UInt_t(fNelem), 1 },
    { "tdc.offset",   &tdc_offset, kFloatV, UInt_t(fNelem), 1 },
    { "tdc.calib",    &tdc_cal,    kFloatV, UInt_t(fNelem), 1 },
    { 0 }
  };
  err = LoadDB( file, date, calib_request, fPrefix );
  fclose(file);

  // Again, no need to continue on errors
  if( err )
    return err;

  // What does this do again?!?!
  DefineAxes( angle*TMath::DegToRad() );

  // Before finishing, prepare vectors that will hold variable output data
  if( !fIsInit ) {
    fBlocks.clear();
    Float_t x = 0;
    Float_t y = 0;
    Float_t z = 0;
    Int_t k = 0;
    // the next three variables are the row,col,layer number starting
    // at fChanMapStart
    int rr = 0;
    int cc = 0;
    int ll = 0;
    for(int r = 0; r < fNrows; r++) {
      rr = r+fChanMapStart;
      for(int c = 0; c < fNcols; c++) {
        cc = c+fChanMapStart;
        for(int l = 0; l < fNlayers; l++) {
          ll = l+fChanMapStart;
          k = blkidx(r,c,l);
          x = xyz[0] + c*dxyz[0];
          y = xyz[1] + r*dxyz[1];
          z = xyz[2] + l*dxyz[2];
          if(!fWithTDC) { // No TDC information
            if(!fWithADCSamples) { // Single-Valued ADC
              fBlocks.push_back( new SBSCalorimeterBlock(x,y,z,rr,cc,ll,
                    adc_ped[k],adc_gain[k]) );
            } else { // Multi-valued ADC
              fBlocks.push_back( new SBSCalorimeterBlockSamples(x,y,z,rr,cc,ll,
                    adc_ped[k],adc_gain[k],adc_ped_mult) );
            }
          } else { // With TDC information
            if(!fWithADCSamples) { // Single-valued ADC + TDC
              fBlocks.push_back( new SBSCalorimeterBlockTDC(x,y,z,rr,cc,ll,
                    adc_ped[k],adc_gain[k],tdc_offset[k],tdc_cal[k]) );
            } else { // Multi-valued ADC + TDC
              fBlocks.push_back( new SBSCalorimeterBlockSamplesTDC(x,y,z,rr,cc,
                    ll,adc_ped[k],adc_gain[k],adc_ped_mult,tdc_offset[k],
                    tdc_cal[k]) );
            }
          }
        }
      }
    }
  }

  // All is well that ends well
  return kOK;
}

//_____________________________________________________________________________
Int_t SBSCalorimeter::DefineVariables( EMode mode )
{
  // Initialize global variables

  if( mode == kDefine && fIsSetup ) return kOK;
  fIsSetup = ( mode == kDefine );

  // Register variables in global list
  RVarDef vars[] = {
    { "row", "Row for block in data vectors",  "fDataOut.fRow" },
    { "col", "Col for block in data vectors",  "fDataOut.fCol" },
    { "layer", "Layer for block in data vectors",  "fDataOut.fLayer" },

    { "a",   "Raw ADC amplitudes",  "fDataOut.fA" },
    { "a_p", "Ped-subtracted ADC amplitudes",  "fDataOut.fA_p" },
    { "a_c", "Calibrated ADC amplitudes",  "fDataOut.fA_c" },
    /*
    //{ "nhit",   "Number of hits",                     "fNhits" },
    { "nclust", "Number of clusters",                 "fNclust" },
    { "e",      "Energy (MeV) of largest cluster",    "fE" },
    { "e_c",    "Corrected Energy (MeV) of largest cluster",    "fE_c" },
    { "x",      "x-position (m) of largest cluster", "fX" },
    { "y",      "y-position (m) of largest cluster", "fY" },
    //{ "targ.x", "x-position (m) of largest cluster in target coords", "fXtarg" },
    //{ "targ.y", "y-position (m) of largest cluster in target coords", "fYtarg" },
    //{ "targ.z", "z-position (m) of largest cluster in target coords", "fZtarg" },
    { "mult",   "Multiplicity of largest cluster",    "fMult" },
    { "nblk",   "Numbers of blocks in main cluster",  "fNblk" },
    { "eblk",   "Energies of blocks in main cluster", "fEblk" },
    //     { "trx",    "track x-position in det plane",      "fTRX" },
    //     { "try",    "track y-position in det plane",      "fTRY" },
    */
    { 0 }
  };
  Int_t err = DefineVarsFromList( vars, mode );
  if( err != kOK)
    return err;

  // Are we using TDCs? If so, define variables for TDCs
  if(fWithTDC) {
    RVarDef vars_tdc[] = {
      { "tdc", "Raw TDC value", "fDataOut.fTDC" },
      { "tdc_c", "Calibrated TDC value", "fDataOut.fTDC_c" },
      { 0 }
    };
    err = DefineVarsFromList( vars_tdc, mode );
    if( err != kOK)
      return err;
  }

  // Are we using multi-valued ADCs? Then define the samples variables
  if(fWithADCSamples) {
    RVarDef vars_samps[] = {
      { "samps_idx", "Index in samples vector for given row-col module",
        "fDataOut.fSampsIdx" },
      { "nsamps" , "Number of samples for given row-col",
        "fDataOut.fNSamps"},
      { "samps",   "RAW ADC samples",  "fDataOut.fSamps" },
      { "samps_p", "Pedestal-subtracted ADC samples",  "fDataOut.fSamps_p" },
      { "samps_c", "Calibrated ADC samples",  "fDataOut.fSamps_c" },
      { 0 }
    };
    err =  DefineVarsFromList(vars_samps, mode);
  }

  return err;
}

//_____________________________________________________________________________
Int_t SBSCalorimeter::Decode( const THaEvData& evdata )
{
  // Decode data and store into the following local data structure:
  //
  // fNhits           - Number of hits on HCal
  // fASamples[][]    - 2D Array of ADC samples/values for each block
  // fASamplesPed[][] - 2D Array of ped subtracted fASamples[][]
  // fASamplesCal[][] - 2D array of ped subtracted and calibrated fASamples[][]
  //
  // (The following are presently now being used)
  // fAsum_p          -  Sum of shower blocks ADC minus pedestal values;
  // fAsum_c          -  Sum of shower blocks corrected ADC values;

  // Clear last event
  ClearEvent();
  //static const char* const here = "Decode()";

  SBSCalorimeterBlock *blk = 0;
  // Loop over all modules in the calorimeter and decode accordingly
  for( UShort_t imod = 0; imod < fDetMap->GetSize(); imod++ ) {
    THaDetMap::Module *d = fDetMap->GetModule( imod );

    for(Int_t ihit = 0; ihit < evdata.GetNumChan( d->crate, d->slot ); ihit++) {
      fNhits++;

      // Get the next available channel, skipping the ones that do not belong
      // to our detector
      Int_t chan = evdata.GetNextChan( d->crate, d->slot, ihit );
      if( chan > d->hi || chan < d->lo ) continue;

      // Get the block index for this crate,slot,channel combo
      blk = fBlocks[fChanMap[imod][chan - d->lo]];
      if(d->IsADC()) {
        DecodeADC(evdata,blk,d,chan);
      } else if ( d->IsTDC()) {
        DecodeTDC(evdata,blk,d,chan);
      }
    }
  }

  return fNhits;
}

Int_t SBSCalorimeter::DecodeADC( const THaEvData& evdata,
    SBSCalorimeterBlock *blk, THaDetMap::Module *d, Int_t chan)
{
  Int_t nhit = evdata.GetNumHits(d->crate, d->slot, chan);
  if(nhit <= 0 )
    return 0;
  if(!fWithADCSamples) {
    // TODO: Again, no clue what to do for multiple hits
    // For now, just take the first one
    blk->ProcessADC( evdata.GetData(d->crate, d->slot, chan, 0));
  } else {
    // As a first assumption, when using multi-valued ADCs then each "hit"
    // will correspond to the number of samples taken.
    std::vector<Float_t> samples;
    samples.resize(nhit);
    for(Int_t i = 0; i < nhit; i++) {
      samples[i] = evdata.GetData(d->crate, d->slot, chan, i);
    }
    dynamic_cast<SBSCalorimeterBlockSamples*>(blk)->ProcessADCSamples( samples );
    samples.clear();
  }

  return 1; // Ignore multiple hits for now
}

Int_t SBSCalorimeter::DecodeTDC( const THaEvData& evdata,
    SBSCalorimeterBlock *blk, THaDetMap::Module *d, Int_t chan)
{
  // TODO: Again, no clue what to do for multiple hits
  // For now, just take the first one
  Int_t nhit = evdata.GetNumHits(d->crate, d->slot, chan);
  if(nhit > 0 ) {
    if(!fWithADCSamples) {
      SBSCalorimeterBlockTDC *blk2 = dynamic_cast<SBSCalorimeterBlockTDC*>(blk);
      blk2->ProcessTDC( evdata.GetData(d->crate, d->slot, chan, 0) );
    } else {
      SBSCalorimeterBlockSamplesTDC *blk2 =
        dynamic_cast<SBSCalorimeterBlockSamplesTDC*>(blk);
      blk2->ProcessTDC( evdata.GetData(d->crate, d->slot, chan, 0) );
    }
  }

  return 1; // Ignore multiple hits for now
}

//_____________________________________________________________________________
void SBSCalorimeter::ClearEvent()
{
  fDataOut.ClearEvent();
  fNhits = 0;
  for(size_t k = 0; k < fBlocks.size(); k++) {
    fBlocks[k]->ClearEvent();
  }
}

Int_t SBSCalorimeter::CoarseProcess(TClonesArray& tracks)
{
  // Make sure to clear old data out
  fDataOut.ClearEvent();
  // Pack data for output to the tree
  SBSCalorimeterBlock *blk = 0;
  SBSCalorimeterBlockTDC *blk_tdc = 0;
  SBSCalorimeterBlockSamples *blk_samps = 0;
  SBSCalorimeterBlockSamplesTDC *blk_tdc_samps = 0;
  size_t nsamples;
  size_t idx;
  for(Int_t k = 0; k < fNelem; k++) {
    blk = fBlocks[k];

    // Skip blocks that have no new data
    if(!blk->HasData())
      continue;

    fDataOut.fRow.push_back(blk->GetRow());
    fDataOut.fCol.push_back(blk->GetCol());
    fDataOut.fLayer.push_back(blk->GetLayer());
    fDataOut.fA.push_back(blk->GetADCDataRaw());
    fDataOut.fA_p.push_back(blk->GetADCDataPed());
    fDataOut.fA_c.push_back(blk->GetADCDataCal());
    if(fWithTDC) {
      if(!fWithADCSamples) {
        blk_tdc = dynamic_cast<SBSCalorimeterBlockTDC*>(blk);
        fDataOut.fTDC.push_back(blk_tdc->GetTDCDataRaw());
        fDataOut.fTDC_c.push_back(blk_tdc->GetTDCDataCal());
      } else {
        blk_tdc_samps = dynamic_cast<SBSCalorimeterBlockSamplesTDC*>(blk);
        fDataOut.fTDC.push_back(blk_tdc_samps->GetTDCDataRaw());
        fDataOut.fTDC_c.push_back(blk_tdc_samps->GetTDCDataCal());
        std::vector<Float_t> &s_r = blk_tdc_samps->GetSamplesDataRaw();
        std::vector<Float_t> &s_p = blk_tdc_samps->GetSamplesDataPed();
        std::vector<Float_t> &s_c = blk_tdc_samps->GetSamplesDataCal();
        nsamples = s_r.size();
        idx = fDataOut.fSamps.size();
        fDataOut.fSampsIdx.push_back(idx);
        fDataOut.fNSamps.push_back(nsamples);
        fDataOut.fSamps.resize(idx+nsamples);
        fDataOut.fSamps_p.resize(idx+nsamples);
        fDataOut.fSamps_c.resize(idx+nsamples);
        for(size_t s = 0; s < nsamples; s++) {
          fDataOut.fSamps[idx+s]   = s_r[s];
          fDataOut.fSamps_p[idx+s] = s_p[s];
          fDataOut.fSamps_c[idx+s] = s_c[s];
        }
      }
    } else if (fWithADCSamples) {
      blk_samps = dynamic_cast<SBSCalorimeterBlockSamples*>(blk);
      std::vector<Float_t> &s_r = blk_samps->GetSamplesDataRaw();
      std::vector<Float_t> &s_p = blk_samps->GetSamplesDataPed();
      std::vector<Float_t> &s_c = blk_samps->GetSamplesDataCal();
      nsamples = s_r.size();
      idx = fDataOut.fSamps.size();
      fDataOut.fSampsIdx.push_back(idx);
      fDataOut.fNSamps.push_back(nsamples);
      fDataOut.fSamps.resize(idx+nsamples);
      fDataOut.fSamps_p.resize(idx+nsamples);
      fDataOut.fSamps_c.resize(idx+nsamples);
      for(size_t s = 0; s < nsamples; s++) {
        fDataOut.fSamps[idx+s]   = s_r[s];
        fDataOut.fSamps_p[idx+s] = s_p[s];
        fDataOut.fSamps_c[idx+s] = s_c[s];
      }
    }
  }
  return 0;
}

Int_t SBSCalorimeter::FineProcess(TClonesArray& tracks)
{
  return 0;
}


/*
//_____________________________________________________________________________
void SBSCalorimeter::DeleteArrays()
{
  // Delete member arrays. Internal function used by destructor.

  delete [] fNChan; fNChan = 0;
  UShort_t mapsize = fDetMap->GetSize();
  if( fChanMap ) {
    for( UShort_t i = 0; i<mapsize; i++ )
      delete [] fChanMap[i];
  }
  delete [] fChanMap; fChanMap = 0;
  delete [] fBlockX;  fBlockX  = 0;
  delete [] fBlockY;  fBlockY  = 0;
  delete [] fPed;     fPed     = 0;
  delete [] fGain;    fGain    = 0;
  delete [] fA;       fA       = 0;
  delete [] fA_p;     fA_p     = 0;
  delete [] fA_c;     fA_c     = 0;
  delete [] fNblk;    fNblk    = 0;
  delete [] fEblk;    fEblk    = 0;
  delete [] fBlocks;  fBlocks  = 0;
  for (int i=0;i<fNrows;i++) {
    delete [] fBlkGrid[i]; fBlkGrid[i] = 0;
  }
  delete [] fBlkGrid; fBlkGrid = 0;
  delete [] fClusters; fClusters = 0;
  delete [] fX; fX = 0;
  delete [] fY; fY = 0;
  //delete [] fXtarg; fXtarg = 0;
  //delete [] fYtarg; fYtarg = 0;
  //delete [] fZtarg; fZtarg = 0;
  delete [] fE; fE = 0;
  delete [] fE_c; fE_c = 0;
  delete [] fMult; fMult = 0;
}

//_____________________________________________________________________________
Int_t SBSCalorimeter::CoarseProcess(TClonesArray& tracks)
{
  // Reconstruct Clusters in shower detector and copy the data 
  // into the following local data structure:
  //
  // fNclust        -  Number of clusters in shower;
  // fE             -  Energy (in MeV) of the "main" cluster;
  // fX             -  X-coordinate (in cm) of the cluster;
  // fY             -  Y-coordinate (in cm) of the cluster;
  // fMult          -  Number of blocks in the cluster;
  // fNblk[0]...[5] -  Numbers of blocks composing the cluster;
  // fEblk[0]...[5] -  Energies in blocks composing the cluster;
  // fTRX;          -  X-coordinate of track cross point with shower plane
  // fTRY;          -  Y-coordinate of track cross point with shower plane
  //

  if( fCoarseProcessed ) return 0;

  Int_t col, row;
  Int_t colmax=0, rowmax=0;
  Double_t  energy_prev = 0.0;

# if not defined(_WIN32)//Win32 compiler do not support variable as array size
  Double_t energyDep[fNcols][fNrows];
# else
  Double_t energyDep[100][100];
# endif

  Double_t energyTotal = 0.0;
  SBSCalorimeterCluster cluster(9);

  //  for( col = 0; col < fNcols; col++ )
  //     {
  //       for( row = 0; row < fNrows; row++ )
  // 	{
  // 	  energyDep[col][row] = 0.0;
  // 	}
  //     }

  //  cout << "Energy Deposition:" << endl <<"___________________________________________________" << endl;
  for( row = 0; row < fNrows; row++ )
  {
    for( col = 0; col < fNcols; col++ )
    {
      energyDep[col][row] = fA_c[BlockColRowToNumber(col,row)]; 

      //	  cout << energyDep[col][row] << " ";
      if( energyDep[col][row] < 0.0 ) 
        energyDep[col][row] = 0.0;
      energyTotal += energyDep[col][row];
    }
    //      cout << endl;
  }

  for( row = 0; row < fNrows; row++ )
  {
    for( col = 0; col < fNcols; col++ )
    {
      if(energyDep[col][row]>energy_prev)
      {
        energy_prev=energyDep[col][row];
        colmax = col;
        rowmax = row;
      }
    }
  }


  //  cout <<"___________________________________________________" << endl;

  Int_t i, j, k=0;
  Double_t energyClusterTotal = 0.0;
  //  Double_t energyClusterGreatest = 0.0;

  Int_t clusterRow = 0;
  Int_t clusterCol = 0;

  //  for( row = 0; row < fNrows; row++ )
  //     {
  //       for( col = 0; col < fNcols; col++ )
  // 	{
  // 	  energyClusterTotal = 0.0;
  // 	  for( i = row-CLUSTER_BLOCK_RADIUS; i <= row+CLUSTER_BLOCK_RADIUS; i++ )
  // 	    {
  // 	      for( j = col-CLUSTER_BLOCK_RADIUS; j <= col+CLUSTER_BLOCK_RADIUS; j++)
  // 		{
  // 		  if( (i >= 0 && i < fNrows ) && ( j >=0 && j < fNcols ) ){   
  // 		    energyClusterTotal += energyDep[j][i];
  // 		  }
  // 		}
  // 	    }

  // 	  if( energyClusterTotal > energyClusterGreatest )
  // 	    {
  // 	      energyClusterGreatest = energyClusterTotal;
  // 	      clusterRow = row;
  // 	      clusterCol = col;
  // 	    }
  // 	}
  //     }
  energyClusterTotal = 0.0; 
  Int_t
    mnrow=TMath::Max(rowmax-CLUSTER_BLOCK_RADIUS,0),
    mxrow=TMath::Min(rowmax+CLUSTER_BLOCK_RADIUS,fNrows-1),
    mncol=TMath::Max(colmax-CLUSTER_BLOCK_RADIUS,0),
    mxcol=TMath::Min(colmax+CLUSTER_BLOCK_RADIUS,fNcols-1);

  for( i = mnrow; i <= mxrow; i++ )
  {
    for( j = mncol; j <= mxcol; j++)
    {
      energyClusterTotal += energyDep[j][i];
      fEblk[k] = energyDep[j][i];
      k++;
    }
  }

  //  cout <<"___________________________________________________" << endl;

  Double_t energyCluster = energyClusterTotal;
  Double_t X, Y;

  if( energyCluster < 0.0 ) return 0;

  //  cout << "Got a cluster!" << endl;
  X = fBlockX[BlockColRowToNumber(colmax, rowmax)];
  Y = fBlockY[BlockColRowToNumber(colmax, rowmax)];

  Double_t energyX = 0.0;
  Double_t energyY = 0.0;

  Int_t  blockcounter = 0;
  for( i = clusterRow-CLUSTER_BLOCK_RADIUS; i <= clusterRow + CLUSTER_BLOCK_RADIUS; i++ )
  {
    for( j = clusterCol-CLUSTER_BLOCK_RADIUS; j <= clusterCol + CLUSTER_BLOCK_RADIUS; j++ )
    {
      if( (i >= 0 && i < fNrows ) && ( j >=0 && j < fNcols ) )
      {
        energyX += energyDep[j][i]*fBlockX[BlockColRowToNumber(j,i)];
        energyY += energyDep[j][i]*fBlockY[BlockColRowToNumber(j,i)];

        cluster.AddBlock( fBlocks[BlockColRowToNumber(j,i)] );
        blockcounter++;
      }
    }
  }

  //  cout << energyCluster << " " << energyX/energyCluster << " " << energyY/ energyCluster << " " << cluster.GetMult() << endl;

  cluster.SetE( energyCluster );
  //cluster.SetX( energyX/energyCluster );
  cluster.SetX( X+fOrigin.X() );
  cluster.SetY( Y+fOrigin.Y() );
  //cluster.SetY( energyY/energyCluster );

  AddCluster(cluster);  

  //  cout << "Added - we now have " << fNclust << endl;

  fCoarseProcessed = 1;
  return 0;

}

//_____________________________________________________________________________
Int_t SBSCalorimeter::FineProcess(TClonesArray& tracks)
{

  if( fFineProcessed ) return 0;

  // Fine Shower processing.

  //   cout << endl << fNclust << " clusters " << GetName()  << endl;
  //   for (int i=0;i<fNclust;i++) {
  //     cout << setw(2) << i << setw(7) << setprecision(1) 
  // 	 << fClusters[i]->GetE() << setw(8) << setprecision(3) 
  // 	 << fClusters[i]->GetX() << setw(8) << fClusters[i]->GetY() 
  // 	 << setw(4) << fClusters[i]->GetMult() << endl;
  //   }

  TVector3 clusterpoint;

  for (int i=0;i<fNclust;i++) {
    //    cout << fClusters[i]->GetE() << " " << fClusters[i]->GetX() << " " << fClusters[i]->GetY() <<fClusters[i]->GetMult()  << endl; 
    fE[i] = fClusters[i]->GetE();
    fE_c[i] = fClusters[i]->GetE()*(gconst + gslope*acc_charge);
    fX[i] = fClusters[i]->GetX();
    fY[i] = fClusters[i]->GetY();
    fMult[i] = fClusters[i]->GetMult();

    //clusterpoint.SetXYZ( fX[i], fY[i], fOrigin.Z() );
    //clusterpoint.Transform(fDetToTarg);
    //clusterpoint = clusterpoint + fDetOffset;

    //fXtarg[i] = clusterpoint.X();
    //fYtarg[i] = clusterpoint.Y();
    //fZtarg[i] = clusterpoint.Z();
    //// We want the shower coordinates in target coordinates

  }

  fFineProcessed = 1;
  return 0;
}

///////////////////////////////////////////////////////////////////////////////

void SBSCalorimeter::AddCluster(SBSCalorimeterCluster* clus) {
  fClusters[fNclust++]=clus;
}

void SBSCalorimeter::AddCluster(SBSCalorimeterCluster& clus) {

  fClusters[fNclust] = new SBSCalorimeterCluster(clus.GetNMaxBlocks());
  fClusters[fNclust]->SetE(clus.GetE());
  fClusters[fNclust]->SetX(clus.GetX());
  fClusters[fNclust]->SetY(clus.GetY());
  fClusters[fNclust++]->SetMult(clus.GetMult());
}

void SBSCalorimeter::RemoveCluster(int i) {
  fNclust--;
  for (int j=i;j<fNclust;j++) fClusters[j]=fClusters[j+1];
}

Int_t SBSCalorimeter::BlockColRowToNumber( Int_t col, Int_t row )
{
  return col*fNrows + row;
}

*/

void SBSCalorimeter::OutputData::ClearEvent()
{
  fRow.clear();
  fCol.clear();
  fA.clear();
  fA_p.clear();
  fA_c.clear();
  fTDC.clear();
  fTDC_c.clear();
  fSampsIdx.clear();
  fNSamps.clear();
  fSamps.clear();
  fSamps_p.clear();
  fSamps_c.clear();
}
