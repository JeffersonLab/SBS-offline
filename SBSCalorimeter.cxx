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
  //THaShower(name,description,apparatus), fNrows(0), fNcols(0),
  fNlayers(0), fWithTDC(false), fWithADCSamples(false), fWithADC(true),
  fMaxNclus(10), fConst(1.0), fSlope(0.0), fAccCharge(0.0)
{
  // Constructor.
  fCoarseProcessed = 0;
  fFineProcessed = 0;
  fEmin = 1.0; // 1 MeV minimum
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

  ClearEvent();
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
  std::vector<Float_t> xyz, dxyz;
  Int_t ncols = 1, nrows = 1, nlayers = 1;
  Float_t angle = 0.0;
  std::vector<Int_t> cluster_dim; // Default is 3x3 if none specified

  // Read mapping/geometry/configuration parameters
  fChanMapStart = 0;
  DBRequest config_request[] = {
    { "detmap",       &detmap,  kIntV }, ///< Detector map
    { "chanmap",      &chanmap, kIntV,    0, true }, ///< Optional channel map
    { "start_chanmap",&fChanMapStart, kInt, 0, true}, ///< Optional start of channel numbering
    { "ncols",        &ncols,   kInt, 1, true }, ///< Number of columns in detector
    { "nrows",        &nrows,   kInt, 1, true }, ///< Number of rows in detector
    { "nlayers",       &nlayers,  kInt,1,true }, ///< [Optional] Number of layers/divisions in each module/block of the detector
    { "angle",        &angle,   kFloat,  0, true },
    { "xyz",           &xyz,      kFloatV, 3 },  ///< center pos of block 1
    { "dxdydz",         &dxyz,     kFloatV, 3 },  ///< block spacing (dx,dy,dz)
    { "emin",         &fEmin,   kFloat, 0, false }, ///< minimum energy threshold
    { "cluster_dim",   &cluster_dim,   kIntV, 0, true }, ///< cluster dimensions (2D)
    { "nmax_cluster",   &fMaxNclus,   kInt, 0, true }, ///< maximum number of clusters to store
    { "const", &fConst, kFloat, 0, true }, ///< const from gain correction 
    { "slope", &fConst, kFloat, 0, true }, ///< slope for gain correction 
    { "acc_charge", &fAccCharge, kFloat, 0, true }, ///< accumulated charge
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

      // Also compute the max possible cluster size (which at most should be
      // cluster_dim x cluster_dim)
      if(cluster_dim.size() == 0) {
        cluster_dim.push_back(3);
        cluster_dim.push_back(3);
      } else if (cluster_dim.size() < 2) {
        cluster_dim.push_back(cluster_dim[0]);
      }
      if(cluster_dim[0] < 1)
        cluster_dim[0] = 3;
      if(cluster_dim[1] < 1)
        cluster_dim[1] = 3;
      fNclubr = TMath::Min( cluster_dim[0], nrows );
      fNclubc = TMath::Min( cluster_dim[1], ncols );
      fNclublk = fNclubr*fNclubc;
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
    UInt_t nmodules = GetDetMap()->GetSize();
    assert( nmodules > 0 );
    fChanMap.resize(nmodules);
    Bool_t makeADC = true;
    for( UInt_t i = 0, k = 0; i < nmodules && !err; i++) {
      THaDetMap::Module *d = GetDetMap()->GetModule(i);
      UInt_t nchan = d->GetNchan();
      if( nchan > 0 ) {
        fChanMap[i].resize(nchan);
        // To simplify finding out which channels are ADCs and which are TDCs
        // we'll just require that the user make the first fNelem channels
        // correspond to the ADCs, and the other fNelem to the TDCs (if in use).
        if(makeADC) {
          d->MakeADC();
        } else {
          d->MakeTDC();
        }
        for(UInt_t chan = 0; chan < nchan; chan++) {
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
    fBlocks.resize(fNelem);
    Float_t x = 0;
    Float_t y = 0;
    Float_t z = 0;
    Int_t k = 0;
    // the next three variables are the row,col,layer number starting
    // at fChanMapStart
    int rr = 0;
    int cc = 0;
    int ll = 0;
    fBlocksGrid.resize(fNrows);
    for(int r = 0; r < fNrows; r++) {
      rr = r+fChanMapStart;
      fBlocksGrid[r].resize(fNcols);
      for(int c = 0; c < fNcols; c++) {
        cc = c+fChanMapStart;
        for(int l = 0; l < fNlayers; l++) {
          fBlocksGrid[r][c].resize(fNlayers);
          ll = l+fChanMapStart;
          k = blkidx(r,c,l);
          x = xyz[0] - c*dxyz[0];
          y = xyz[1] - r*dxyz[1];
          z = xyz[2] - l*dxyz[2];
          SBSCalorimeterBlock *blk = new SBSCalorimeterBlock(x,y,z,rr,cc,ll);
          if(fWithADC) { // Single-valued ADC version
            blk->SetADC(adc_ped[k],adc_gain[k]);
          }
          if(fWithTDC) { // TDC info
            blk->SetTDC(tdc_offset[k],tdc_cal[k]);
          }
          if(fWithADCSamples) {
            blk->SetSamples(adc_ped[k],adc_gain[k]);
          }
          fBlocks[k] = blk;
          fBlocksGrid[r][c][l] = blk;
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
  if( !( fWithADC || fWithTDC || fWithADCSamples) ) {
    Error( Here("DefineVariables"),
        "Calorimeter %s defined with no data payload.",GetName());
    return kInitError;
  }
  fIsSetup = ( mode == kDefine );

  // Most of these variables were used previously, and so I'll leave
  // them here to preserve old analysis macros people may have written.
  // This includes things like fE, fNblk, fE_c, etc...
  RVarDef vars[] = {
    { "row", "Row for block in data vectors",  "fRow" },
    { "col", "Col for block in data vectors",  "fCol" },
    { "layer", "Layer for block in data vectors",  "fLayer" },
    { "nhit",   "Number of hits",                     "fNhits" },
    { "nclust", "Number of clusters meeting threshold", "fNclus" },
    { "e",      "Energy (MeV) of largest cluster",    "fE" },
    { "e_c",    "Corrected Energy (MeV) of largest cluster",    "fE_c" },
    { "eblk",   "Energy (MeV) of highest energy block in the largest cluster",    "fEblk" },
    { "eblk_c", "Corrected Energy (MeV) of highest energy block in the largest cluster",    "fEblk_c" },
    { "rowblk", "Row of block with highest energy in the largest cluster",    "fRowblk" },
    { "colblk", "Col of block with highest energy in the largest cluster",    "fColblk" },
    { "x",      "x-position (mm) of largest cluster", "fX" },
    { "y",      "y-position (mm) of largest cluster", "fY" },
    { "nblk",   "Number of blocks in the largest cluster",    "fNblk" },
    { "clus.e", "Energy (MeV) of clusters", "fEclus"},
    { "clus.e_c","Energy (MeV) of clusters", "fEclus_c"},
    { "clus.eblk", "Energy (MeV) block with highest E of clusters", "fEclusBlk"},
    { "clus.eblk_c","Energy (MeV) block with highest E of clusters", "fEclusBlk_c"},
    { "clus.x", "x-position (m) of clusters", "fXclus"},
    { "clus.y", "x-position (m) of clusters", "fYclus"},
    { "clus.nblk","Number of blocks in the clusters",    "fNblkclus" },
    { "clus.row", "Row of block with highest E in the cluster",    "fRowblkclus" },
    { "clus.col", "Col of block with highest E in the cluster",    "fColblkclus" },
    //{ "targ.x", "x-position (m) of largest cluster in target coords", "fXtarg" },
    //{ "targ.y", "y-position (m) of largest cluster in target coords", "fYtarg" },
    //{ "targ.z", "z-position (m) of largest cluster in target coords", "fZtarg" },
    //     { "trx",    "track x-position in det plane",      "fTRX" },
    //     { "try",    "track y-position in det plane",      "fTRY" },
    { 0 }
  };

  Int_t err = DefineVarsFromList( vars, mode );
  if( err != kOK)
    return err;

  // Do we have a single valued ADC? Then define single-valued adc variables
  if(fWithADC||fWithADCSamples) {
    // Register variables in global list
    RVarDef vars_adc[] = {
      { "a",   "Raw ADC amplitudes",  "fA" },
      { "a_p", "Ped-subtracted ADC amplitudes",  "fA_p" },
      { "a_c", "Calibrated ADC amplitudes",  "fA_c" },
      { 0 }
    };
    err = DefineVarsFromList( vars_adc, mode );
    if( err != kOK)
      return err;
  }

  // Are we using TDCs? If so, define variables for TDCs
  if(fWithTDC) {
    RVarDef vars_tdc[] = {
      { "tdc", "Raw TDC value", "fTDC" },
      { "tdc_c", "Calibrated TDC value", "fTDC_c" },
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
        "fSampsIdx" },
      { "nsamps" , "Number of samples for given row-col",
        "fNsamps"},
      { "samps",   "RAW ADC samples",  "fSamps" },
      { "samps_p", "Pedestal-subtracted ADC samples",  "fSamps_p" },
      { "samps_c", "Calibrated ADC samples",  "fSamps_c" },
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
  for( UInt_t imod = 0; imod < fDetMap->GetSize(); imod++ ) {
    THaDetMap::Module *d = fDetMap->GetModule( imod );

    for(UInt_t ihit = 0; ihit < evdata.GetNumChan( d->crate, d->slot ); ihit++) {
      fNhits++;

      // Get the next available channel, skipping the ones that do not belong
      // to our detector
      UInt_t chan = evdata.GetNextChan( d->crate, d->slot, ihit );
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
  UInt_t nhit = evdata.GetNumHits(d->crate, d->slot, chan);
  //std::cout << d->crate << " " << d->slot << " " << chan << std::endl;
  if(nhit <= 0 )
    return 0;
  if(fWithADC && !fWithADCSamples && blk->ADC()) {
    // TODO: Again, no clue what to do for multiple hits
    // For now, just take the first one
    blk->ADC()->Process( evdata.GetData(d->crate, d->slot, chan, 0));
  }
  if(fWithADCSamples && !fWithADC && blk->Samples()) {
    // As a first assumption, when using multi-valued ADCs then each "hit"
    // will correspond to the number of samples taken.
    std::vector<Float_t> samples;
    samples.resize(nhit);
    for(Int_t i = 0; i < nhit; i++) {
      samples[i] = evdata.GetData(d->crate, d->slot, chan, i);
    }
    blk->Samples()->Process(samples);
    samples.clear();
  }

  return 1; // Ignore multiple hits for now
}

Int_t SBSCalorimeter::DecodeTDC( const THaEvData& evdata,
    SBSCalorimeterBlock *blk, THaDetMap::Module *d, Int_t chan)
{
  if(!fWithTDC || !blk->TDC())
    return 0;

  // TODO: Again, no clue what to do for multiple hits
  // For now, just take the first one
  UInt_t nhit = evdata.GetNumHits(d->crate, d->slot, chan);
  if(nhit > 0 ) {
    blk->TDC()->Process( evdata.GetData(d->crate, d->slot, chan, 0) );
  }

  return 1; // Ignore multiple hits for now
}

//_____________________________________________________________________________
void SBSCalorimeter::ClearEvent()
{
  ClearOutputVariables();
  fNhits = 0;
  fCoarseProcessed = 0;
  fFineProcessed = 0;
  for(size_t k = 0; k < fBlocks.size(); k++) {
    fBlocks[k]->ClearEvent();
  }
  for(size_t k = 0; k < fClusters.size(); k++) {
    if(fClusters[k])
      delete fClusters[k];
  }
  fClusters.clear();
}

Int_t SBSCalorimeter::CoarseProcess(TClonesArray& )// tracks)
{
  // Make sure we haven't already been called in this event
  if(fCoarseProcessed) return 0;

  // Pack simple data for output to the tree, and call CoarseProcess on blocks
  SBSCalorimeterBlock *blk = 0;
  size_t nsamples;
  size_t idx;
  for(Int_t k = 0; k < fNelem; k++) {
    blk = fBlocks[k];
    blk->CoarseProcess();

    // Skip blocks that have no new data
    if(!blk->HasData())
      continue;

    fRow.push_back(blk->GetRow());
    fCol.push_back(blk->GetCol());
    fLayer.push_back(blk->GetLayer());
    if(fWithADC && blk->ADC() && !fWithADCSamples) {
      fA.push_back(blk->ADC()->GetDataRaw());
      fA_p.push_back(blk->ADC()->GetDataPed());
      fA_c.push_back(blk->ADC()->GetDataCal());
    }

    if(fWithTDC && blk->TDC()) {
      fTDC.push_back(blk->TDC()->GetDataRaw());
      fTDC_c.push_back(blk->TDC()->GetDataCal());
    }

    if (fWithADCSamples && blk->Samples()) {
      std::vector<Float_t> &s_r = blk->Samples()->GetDataRaw();
      std::vector<Float_t> &s_p = blk->Samples()->GetDataPed();
      std::vector<Float_t> &s_c = blk->Samples()->GetDataCal();
      nsamples = s_r.size();
      idx = fSamps.size();
      fSampsIdx.push_back(idx);
      fNsamps.push_back(nsamples);
      fSamps.resize(idx+nsamples);
      fSamps_p.resize(idx+nsamples);
      fSamps_c.resize(idx+nsamples);
      for(size_t s = 0; s < nsamples; s++) {
        fSamps[idx+s]   = s_r[s];
        fSamps_p[idx+s] = s_p[s];
        fSamps_c[idx+s] = s_c[s];
      }
      fA.push_back(blk->Samples()->GetDataSumRaw());
      fA_p.push_back(blk->Samples()->GetDataSumPed());
      fA_c.push_back(blk->Samples()->GetDataSumCal());
    }
  }

  // Now, find as many clusters as meet the minimum energy
  for(Int_t r = 0; r <= fNrows-fNclubr; r++) {
    for(Int_t c = 0; c <= fNcols-fNclubc; c++) {
      for(Int_t l = 0; l < fNlayers; l++) {

        // Now perform the sum
        SBSCalorimeterCluster *clus = new SBSCalorimeterCluster(fNclublk);
        for(Int_t rr = 0; rr < fNclubr; rr++) {
          for(Int_t cc = 0; cc < fNclubc; cc++) {
            blk = fBlocksGrid[r+rr][c+cc][l];
            if(blk->GetE()>0)
              clus->AddBlock(blk);
          }
        }
        if(clus->GetE() >= fEmin) {
          fClusters.push_back(clus);
          fNclus++;
        } else {
          // Cluster did not meet the minimum threshold, so delete now to
          // avoid memory leak later.
          delete clus;
        }
      }
    }
  }
  fCoarseProcessed = 1;
  return 0;
}

//_____________________________________________________________________________
Int_t SBSCalorimeter::FineProcess(TClonesArray&)//tracks)
{

  if( fFineProcessed ) return 0;
  // Now, sort the clusters by energy
  std::sort(fClusters.rbegin(),fClusters.rend(),SBSCalorimeterClusterCompare());
  // Get information on the cluster with highest energy (useful even if
  // fMaxNclus is zero, i.e., storing no vector of clusters)
  if(fClusters.size()>0) {
    fE    = fClusters[0]->GetE();
    fE_c  = fClusters[0]->GetE()*(fConst + fSlope*fAccCharge);
    fX    = fClusters[0]->GetX();
    fY    = fClusters[0]->GetY();
    fNblk = fClusters[0]->GetMult();
    fEblk    = fClusters[0]->GetEblk();
    fEblk_c  = fClusters[0]->GetEblk()*(fConst + fSlope*fAccCharge);
    fRowblk = fClusters[0]->GetRow();
    fColblk = fClusters[0]->GetCol();
  }
  // Now store the remaining clusters (up to fMaxNclus)
  Int_t nres = TMath::Min(Int_t(fMaxNclus),Int_t(fClusters.size()));
  fEclus.reserve(nres);
  fEclus_c.reserve(nres);
  fXclus.reserve(nres);
  fYclus.reserve(nres);
  fNblkclus.reserve(nres);
  fEclusBlk.reserve(nres);
  fEclusBlk_c.reserve(nres);
  fRowblkclus.reserve(nres);
  fColblkclus.reserve(nres);
  for(size_t i = 0; i < fClusters.size(); i++) {
    if(fNclus <= fMaxNclus) { // Keep adding them until we reach fMaxNclus
      fEclus.push_back(fClusters[i]->GetE());
      fXclus.push_back(fClusters[i]->GetX());
      fYclus.push_back(fClusters[i]->GetY());
      fEclus_c.push_back(fClusters[i]->GetE()*(fConst + fSlope*fAccCharge));
      fNblkclus.push_back(fClusters[i]->GetMult());
      fEclusBlk.push_back(fClusters[i]->GetEblk());
      fEclusBlk_c.push_back(fClusters[i]->GetEblk()*(fConst + fSlope*fAccCharge));
      fRowblkclus.push_back(fClusters[i]->GetRow());
      fColblkclus.push_back(fClusters[i]->GetCol());
    }
  }


  /*
  TVector3 clusterpoint;

  for (int i=0;i<fNclust;i++) {
    //clusterpoint.SetXYZ( fX[i], fY[i], fOrigin.Z() );
    //clusterpoint.Transform(fDetToTarg);
    //clusterpoint = clusterpoint + fDetOffset;

    //fXtarg[i] = clusterpoint.X();
    //fYtarg[i] = clusterpoint.Y();
    //fZtarg[i] = clusterpoint.Z();
    //// We want the shower coordinates in target coordinates

  }
  */

  fFineProcessed = 1;
  return 0;
}

void SBSCalorimeter::ClearOutputVariables()
{
  fRow.clear();
  fCol.clear();
  fLayer.clear();
  fA.clear();
  fA_p.clear();
  fA_c.clear();
  fTDC.clear();
  fTDC_c.clear();
  fSampsIdx.clear();
  fNsamps.clear();
  fSamps.clear();
  fSamps_p.clear();
  fSamps_c.clear();
  fEclus.clear();
  fEclus_c.clear();
  fXclus.clear();
  fYclus.clear();
  fNblkclus.clear();
  fEclusBlk.clear();
  fEclusBlk_c.clear();
  fRowblkclus.clear();
  fColblkclus.clear();
  fE = fE_c = fX = fY = fNblk = fNclus = 0;
}

void SBSCalorimeter::SetWithADC(Bool_t var)
{
  fWithADC = var;
  if(var)
    fWithADCSamples = false;
}

void SBSCalorimeter::SetWithADCSamples(Bool_t var)
{
  fWithADCSamples = var;
  if(var)
    fWithADC = false;
}
