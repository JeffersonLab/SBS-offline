///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// SBSBBShower                                                               //
//                                                                           //
// Shower counter class, describing a generic segmented shower detector      //
// (preshower or shower).                                                    //
// Currently, only the "main" cluster, i.e. cluster with the largest energy  //
// deposition is considered. Units of measurements are MeV for energy of     //
// shower and centimeters for coordinates.                                   //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "SBSBBShower.h"

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
#include <cassert>
#ifdef HAS_SSTREAM
 #include <sstream>
 #define OSSTREAM ostringstream
#else
 #include <strstream>
 #define OSSTREAM ostrstream
#endif

using namespace std;

ClassImp(SBSBBShower)

//_____________________________________________________________________________
SBSBBShower::SBSBBShower( const char* name, const char* description,
                         THaApparatus* apparatus ) :
//THaPidDetector(name,description,apparatus), fNChan(NULL), fChanMap(NULL)
THaShower(name,description,apparatus), //fNChan(NULL), fChanMap(NULL)
  fE_cl(0), fX_cl(0), fY_cl(0), fMult_cl(0), fNblk_cl(0), fEblk_cl(0), fE_cl_corr(0), 
  fBlocks(0), fClusters(0), fMCdata(0)//, fBlkGrid(0)
{
  // Constructor.
  fCoarseProcessed = 0;
  fFineProcessed = 0;
  //Clear();
}


//_____________________________________________________________________________
Int_t SBSBBShower::ReadDatabase( const TDatime& date )
{
    // Read this detector's parameters from the database file 'fi'.
    // This function is called by THaDetectorBase::Init() once at the
    // beginning of the analysis.
    // 'date' contains the date/time of the run being analyzed.
  
  //cout << "******readdatabase *******" << endl;

  static const char* const here = "ReadDatabase()";
  
  FILE* file = OpenFile( date );
  if( !file ) return kFileError;

  // Read fOrigin and fSize (required!)
  Int_t err = ReadGeometry( file, date, true );
  if( err ) {
//cout << " readgeo err" << endl;
    fclose(file);
    return err;
  }

  vector<Int_t> detmap, chanmap;
  vector<Double_t> xy, dxy;
  Int_t ncols, nrows;
  fMaxNClust=10;
  // Read mapping/geometry/configuration parameters
  DBRequest config_request[] = {
    { "detmap",       &detmap,      kIntV     },
    { "chanmap",      &chanmap,     kIntV,    0, 1 },
    { "ncols",        &ncols,       kInt      },
    { "nrows",        &nrows,       kInt      },
    { "xy",           &xy,          kDoubleV, 2 },  // center pos of block 1
    { "dxdy",         &dxy,         kDoubleV, 2 },  // dx and dy block spacings
    { "thr_adc",      &fThrADC,     kDouble,     0, 1 },
    { "emin",         &fEmin,       kFloat,   0, 1 },
    { "clus_rad",     &fClusRadius, kFloat,   0, 1 },
    { "maxclust",     &fMaxNClust,  kInt,     0, 1 },
    { "mc_data",      &fMCdata,     kInt,     0, 1 },// flag for MC data
    { 0 }
  };
  //cout << " call Loaddb" << endl;
  err = LoadDB( file, date, config_request, fPrefix );

  fClusBlockRadX = Int_t(fClusRadius/dxy[0]);
  fClusBlockRadY = Int_t(fClusRadius/dxy[1]);
  
  //cout << " nrows  " << nrows << " " <<  ncols << endl;
  // Sanity checks
  if( !err && (nrows <= 0 || ncols <= 0) ) {
    Error( Here(here), "Illegal number of rows or columns: %d %d. Must be 0. "
	   "Fix database.", nrows, ncols );
    err = kInitError;
  }

  Int_t nelem = ncols * nrows; 
  Int_t nclbl = TMath::Min( 3, nrows ) * TMath::Min( 3, ncols );
  
  // Reinitialization only possible for same basic configuration
  if( !err ) {
    if( fIsInit && nelem != fNelem ) {
      Error( Here(here), "Cannot re-initalize with different number of blocks or "
	     "blocks per cluster (was: %d, now: %d). Detector not re-initialized.",
	     fNelem, nelem );
      err = kInitError;
    } else {
      fNelem = nelem;
      fNrows = nrows;
      fNcols = ncols;
      fNclublk = nclbl;
    }
  }

  if( !err ) {
    // Clear out the old channel map before reading a new one
//cout << " clea map " << endl;
    fChanMap.clear();
    if( FillDetMap(detmap, 0, here) <= 0 ) {
      err = kInitError;  // Error already printed by FillDetMap
      //    } else if( (nelem = fDetMap->GetTotNumChan()) != fNelem ) {
      //Error( Here(here), "Number of detector map channels (%d) "
      //	     "inconsistent with number of blocks (%d)", nelem, fNelem );
      //err = kInitError;
    }
  }
//cout << " filled map" << endl;
  if( !err ) {
    if( !chanmap.empty() ) {
      // If a map is found in the database, ensure it has the correct size
      //Int_t cmapsize = chanmap.size();
      //if( cmapsize != fNelem ) {
      //	Error( Here(here), "Channel map size (%d) and number of detector "
      //	       "channels (%d) must be equal. Fix database.", cmapsize, fNelem );
      //	err = kInitError;
      //}
    }
    if( !err ) {
      // Set up the new channel map
      Int_t nmodules = fDetMap->GetSize();
      if(fDebug>=2)cout << "Set up the new channel map " << nmodules << endl;
      assert( nmodules > 0 );
      fChanMap.resize(nmodules);
      for( Int_t i=0, k=0; i < nmodules && !err; i++ ) {
	THaDetMap::Module* module = fDetMap->GetModule(i);
	Int_t nchan = module->hi - module->lo + 1;
        if(fDebug>=2)cout << " nchan = " << nchan << endl;
	if( nchan > 0 ) {
	  fChanMap.at(i).resize(nchan);
	  for( Int_t j=0; j<nchan; ++j ) {
	    if(fDebug>=2)cout << " k = " << k << " " << nchan*nmodules << endl;
	    assert( k < nmodules*nchan );
	    fChanMap.at(i).at(j) = chanmap.empty() ? k : chanmap[k]-1;
	    if(fDebug>=2)cout << " k = " << k << " " << nchan*nmodules << " " << chanmap[k] << " " << fChanMap.at(i).at(j)<< endl;
	    ++k;
	  }
	} else {
	  Error( Here(here), "No channels defined for module %d.", i);
	  fChanMap.clear();
	  err = kInitError;
	}
      }
    }
  }

//cout << " close file" << file << endl;
  if( err ) {
    fclose(file);
    return err;
  }

  // Dimension arrays
  //FIXME: use a structure!
  UInt_t nval = fNelem;
  if( !fIsInit ) {
    // Geometry
    fBlockX = new Float_t[ nval ];
    fBlockY = new Float_t[ nval ];

    // Calibrations
    fPed    = new Float_t[ nval ];
    fGain   = new Float_t[ nval ];

    // Per-event data
    fA    = new Float_t[ nval ];
    fA_p  = new Float_t[ nval ];
    fA_c  = new Float_t[ nval ];
    fNblk = new Int_t[ fNclublk ];
    fEblk = new Float_t[ fNclublk ];

    fE_cl = new Float_t[ fMaxNClust ];
    fX_cl = new Float_t[ fMaxNClust ];
    fY_cl = new Float_t[ fMaxNClust ];
    fMult_cl = new Int_t[ fMaxNClust ];
    fE_cl_corr = new Float_t[ fMaxNClust ];
    
    fNblk_cl = new Int_t*[ fMaxNClust ];
    fEblk_cl = new Float_t*[ fMaxNClust ];
    for(Int_t k = 0; k<fMaxNClust; k++){
      fNblk_cl[k] = new Int_t[ fNclublk ];
      fEblk_cl[k] = new Float_t[ fNclublk ];
    }
    if(fMCdata){
      fE_cl_res = new Float_t[ fMaxNClust ];
      fX_cl_res = new Float_t[ fMaxNClust ];
      fY_cl_res = new Float_t[ fMaxNClust ];
    }
    
    fIsInit = true;
  }

//cout << " define geometry" << endl;
  // Compute block positions
  for( int c=0; c<ncols; c++ ) {
    for( int r=0; r<nrows; r++ ) {
      int k = nrows*c + r;
      // Units are meters
      fBlockX[k] = xy[0] + r*dxy[0];
      fBlockY[k] = xy[1] - c*dxy[1];
      if(fDebug){
	cout << " k " << k << " r " << r << " c " << c 
	     << " x " << xy[0] << " dx " << dxy[0] << " => " <<  fBlockX[k]
	     << " y " << xy[1] << " dy " << dxy[1] << " => " <<  fBlockY[k] 
	     << endl;
      }
    }
  }

  // Read calibration parameters

  // Set DEFAULT values here
  // Default ADC pedestals (0) and ADC gains (1)
//cout << " set ped memory " << endl;
  memset( fPed, 0, nval*sizeof(fPed[0]) );
//cout << " set gain" << endl;
  for( UInt_t i=0; i<nval; ++i ) { fGain[i] = 1.0; }

  // Read ADC pedestals and gains (in order of logical channel number)
  DBRequest calib_request[] = {
    { "pedestal",    fPed,   kFloat, nval, 1 },
    { "adc_calib",        fGain,  kFloat, nval, 1 },
    { 0 }
  };
  err = LoadDB( file, date, calib_request, fPrefix );
  fclose(file);
  if( err )
    return err;

#ifdef WITH_DEBUG
  // Debug printout
  if ( fDebug 2 ) {
    const UInt_t N = static_cast<UInt_t>(fNelem);
    Double_t pos[3]; fOrigin.GetXYZ(pos);
    DBRequest list[] = {
      { "Number of blocks",       &fNelem,     kInt        },
      { "Detector center",        pos,         kDouble, 3  },
      { "Detector size",          fSize,       kDouble, 3  },
      { "Channel map",            &chanmap,    kIntV       },
      { "Position of block 1",    &xy,         kDoubleV    },
      { "Block x/y spacings",     &dxy,        kDoubleV    },
      { "Minimum cluster energy", &fEmin,      kFloat,  1  },
      { "ADC pedestals",          fPed,        kFloat,  N  },
      { "ADC pedestals",          fPed,        kFloat,  N  },
      { "ADC gains",              fGain,       kFloat,  N  },
      { 0 }
    };
    DebugPrint( list );
  }
#endif
  fClusters = new SBSBBShowerCluster*[fMaxNClust];
  fBlocks = new SBSShowerBlock*[fNelem];
  
  cout << " fEmin " << fEmin << " fClusterRadius " << fClusRadius 
       << " fClusBlockRadX " << fClusBlockRadX 
       << " fClusBlockRadY " << fClusBlockRadY << endl;
  
  for(int k = 0; k<fNelem;k++){
    int row = k%nrows;
    int col = (k-row)/nrows;
    if(fDebug){
      cout << " k " << k << " row " << row << " col " << col 
	   << " gain " << fGain[k] << " ped " << fPed[k] << endl;
    }
    fBlocks[k] = new SBSShowerBlock(fBlockX[k], fBlockY[k], fPed[k], fGain[k], row, col);
  }
  //cout << " retruning" << endl;
  return kOK;
}


//_____________________________________________________________________________
Int_t SBSBBShower::DefineVariables( EMode mode )
{
  // Initialize global variables

  if( mode == kDefine && fIsSetup ) return kOK;
  fIsSetup = ( mode == kDefine );

  // Register variables in global list

  RVarDef vars[] = {
    { "nhit",   "Number of hits",                     "fNhits" },
    { "a",      "Raw ADC amplitudes",                 "fA" },
    { "a_p",    "Ped-subtracted ADC amplitudes",      "fA_p" },
    { "a_c",    "Calibrated ADC amplitudes",          "fA_c" },
    { "asum_p", "Sum of ped-subtracted ADCs",         "fAsum_p" },
    { "asum_c", "Sum of calibrated ADCs",             "fAsum_c" },
    { "nclust", "Number of clusters",                 "fNclust" },
    { "e_m",      "Energy (MeV) of main (high E) cluster",    "fE" },
    { "e_c_m",    "Corrected Energy (MeV) of main (high E) cluster",    "fE_corr" },
    { "x_m",      "x-position (m) of main (high E) cluster", "fX" },
    { "y_m",      "y-position (m) of main (high E) cluster", "fY" },
    { "mult_m",   "Multiplicity of main (high E) cluster",    "fMult" },
    { "nblk_m",   "Numbers of blocks in main (high E) cluster",  "fNblk" },
    { "eblk_m",   "Energies of blocks in main (high E) cluster", "fEblk" },
    //     { "trx",    "track x-position in det plane",      "fTRX" },
    //     { "try",    "track y-position in det plane",      "fTRY" },
    { "e",      "Energy (MeV) of all clusters", "fE_cl" },
    { "e_c",    "Corrected Energy (MeV) of all clusters", "fE_cl_corr" },
    { "x",      "x-position (m) of all clusters", "fX_cl" },
    { "y",      "y-position (m) of all clusters", "fY_cl" },
    { "mult",   "Multiplicity of all clusters",    "fMult_cl" },
    { "nblk",   "Numbers of blocks in all clusters",  "fNblk_cl" },
    { "eblk",   "Energies of blocks in all clusters", "fEblk_cl" },
    { 0 }
  };
  Int_t ret = DefineVarsFromList( vars, mode );
  if( ret != kOK )
    return ret;
  
  if(fMCdata){
    RVarDef varsmc[] = {
      { "e_m_res", "Energy resolution of main cluster",    "fEres" },
      { "x_m_res", "x-position resolution (m) of main cluster", "fXres" },
      { "y_m_res", "y-position resolution  (m) of main cluster", "fYres" },
      { "e_res", "Energy resolution of all clusters",    "fE_cl_res" },
      { "x_res", "x-position resolution (m) of all clusters", "fX_cl_res" },
      { "y_res", "y-position resolution  (m) of all clusters", "fY_cl_res" },
      { 0 }
    };
    
    ret = DefineVarsFromList( varsmc, mode );
    if( ret != kOK )
      return ret;
  }
  
  return ret;
}

//_____________________________________________________________________________
SBSBBShower::~SBSBBShower()
{
    // Destructor. Removes internal arrays and global variables.

    if( fIsSetup )
        RemoveVariables();
    if( fIsInit )
        DeleteArrays();
}

//_____________________________________________________________________________
void SBSBBShower::DeleteArrays()
{
  // Delete member arrays. Internal function used by destructor.
  cout << "SBSBBShower::DeleteArrays() " << endl;
  // ? THaShower::DeleteArrays()
  //delete [] fNChan; fNChan = 0;
  //UShort_t mapsize = fDetMap->GetSize();
  // if( fChanMap ) {
  //     for( UShort_t i = 0; i<mapsize; i++ )
  //         delete [] fChanMap[i];
  // }
  //delete [] fChanMap; fChanMap = 0;
  fChanMap.clear();
  delete [] fBlockX;  fBlockX  = 0;
  delete [] fBlockY;  fBlockY  = 0;
  delete [] fPed;     fPed     = 0;
  delete [] fGain;    fGain    = 0;
  delete [] fA;       fA       = 0;
  delete [] fA_p;     fA_p     = 0;
  delete [] fA_c;     fA_c     = 0;
  delete [] fBlocks;  fBlocks  = 0;
  //for (int i=0;i<fNrows;i++) {
  //    delete [] fBlkGrid[i]; fBlkGrid[i] = 0;
  //}
  //delete [] fBlkGrid; fBlkGrid = 0;
  delete [] fClusters; fClusters = 0;
  //delete [] fXtarg; fXtarg = 0;
  //delete [] fYtarg; fYtarg = 0;
  //delete [] fZtarg; fZtarg = 0;
  delete [] fE_cl; fE_cl = 0;
  delete [] fX_cl; fX_cl = 0;
  delete [] fY_cl; fY_cl = 0;
  delete [] fMult_cl; fMult_cl = 0;
  /*
  for(int i = 0; i<fNclust; i++){
    delete [] fNblk_cl[i]; fNblk_cl = 0;
    delete [] fEblk_cl[i]; fEblk_cl = 0;
  }
  */
  delete [] fNblk_cl; fNblk_cl = 0;
  delete [] fEblk_cl; fEblk_cl = 0; 
  delete [] fNblk; fNblk = 0;
  delete [] fEblk; fEblk = 0;
  delete [] fE_cl_corr; fE_cl_corr = 0;
  if(fMCdata){
    delete [] fE_cl_res; fE_cl_res = 0;
    delete [] fX_cl_res; fX_cl_res = 0;
    delete [] fY_cl_res; fY_cl_res = 0;
  }
}

//_____________________________________________________________________________
inline
void SBSBBShower::Clear( Option_t* opt )//Event()//
{
  // Reset all local data to prepare for next event.
  
  if(fDebug)cout << " resetting variables for new event" << endl;
  
  fCoarseProcessed = 0;
  fFineProcessed = 0;
  
  const int lsh = fNelem*sizeof(Float_t);
  const int lshh = fMaxNClust*sizeof(Float_t);
  const int lsc = fNclublk*sizeof(Float_t);
  const int lsi = fNclublk*sizeof(Int_t);
  //const int lsch = fNclublk*fMaxNClust*sizeof(Float_t);
  //const int lsih = fNclublk*fMaxNClust*sizeof(Int_t);
  const int lsj = fMaxNClust*sizeof(Int_t);
  
  if(fDebug){
    cout << " lsh " << lsh << " lshh " << lshh << " lsc " << lsc   << " lsi " << lsi  << " lsj " << lsj << endl;
    cout << " lsh " << lsh << " sizeof(fA) " << sizeof(*fA)*fNelem << " " << sizeof(*fA_p)*fNelem << " " << sizeof(*fA_c)*fNelem << endl;
    cout << " lshh " << lshh << " sizeof(fE_cl) " << sizeof(*fE_cl)*fMaxNClust << " " << sizeof(*fE_cl_corr)*fMaxNClust << " " << sizeof(*fX_cl)*fMaxNClust << " " << sizeof(*fY_cl)*fMaxNClust << endl;
    cout << " lsc " << lsc << " sizeof(fEblk) " << sizeof(*fEblk)*fNclublk << " * " << fMaxNClust << " = " << sizeof(**fEblk_cl)*fNclublk*fMaxNClust << endl;
    cout << " lsi " << lsi << " sizeof(fNblk) " << sizeof(*fNblk)*fNclublk << " * " << fMaxNClust << " = " << sizeof(**fNblk_cl)*fNclublk*fMaxNClust << endl;
    cout << " lsj " << lsj << " sizeof(fMult_cl) " << sizeof(*fMult_cl)*fMaxNClust << endl;
  }
  
  fNhits = 0;
  memset( fA, 0, lsh );
  memset( fA_p, 0, lsh );
  memset( fA_c, 0, lsh );
  memset( fE_cl, 0, lshh );
  memset( fE_cl_corr, 0, lshh );
  memset( fX_cl, 0, lshh );
  memset( fY_cl, 0, lshh );
  //memset( fXtarg, 0, lshh );
  //memset( fYtarg, 0, lshh );
  //memset( fZtarg, 0, lshh );
  fAsum_p = 0.0;
  fAsum_c = 0.0;
  fNclust = 0;
  memset( fNblk, 0, lsi );
  memset( fEblk, 0, lsc );
  for(int i = 0; i<fMaxNClust; i++){
    memset( fNblk_cl[i], 0, lsi );
    memset( fEblk_cl[i], 0, lsc );
  }
  memset( fMult_cl, 0, lsj );
  fE = 0.0;
  fX = 0.0;
  fY = 0.0;
  fMult = 0.0;
  fTRX = 0.0;
  fTRY = 0.0;
  
  if(fMCdata){
    fEres = 0.0;
    fXres = 0.0;
    fYres = 0.0;
    memset( fE_cl_res, 0, lshh );
    memset( fX_cl_res, 0, lshh );
    memset( fY_cl_res, 0, lshh );
  }
  
  for (int i=0;i<fNelem;i++) 
    if(fBlocks[i])fBlocks[i]->ClearEvent();
  for (int i=0;i<fNclust;i++) 
    if(fClusters[i])fClusters[i]->ClearEvent();
  
  
  /*
    for (int i=0;i<fNrows;i++)
    for (int j=0;j<fNcols;j++)  
    fBlkGrid[i][j]->ClearEvent();
  */
  if(fDebug)cout << "Done Clearing" << endl;
}

//_____________________________________________________________________________
Int_t SBSBBShower::Decode( const THaEvData& evdata )
{
  // Decode shower data, scale the data to energy deposition
  // ( in MeV ), and copy the data into the following local data structure:
  //
  // fNhits           -  Number of hits on shower;
  // fA[]             -  Array of ADC values of shower blocks;
  // fA_p[]           -  Array of ADC minus ped values of shower blocks;
  // fA_c[]           -  Array of corrected ADC values of shower blocks;
  // fAsum_p          -  Sum of shower blocks ADC minus pedestal values;
  // fAsum_c          -  Sum of shower blocks corrected ADC values;

  const char* const here = "Decode";
  // Loop over all modules defined for shower detector
  bool has_warning = false;
  Int_t nmodules = fDetMap->GetSize();
  for( Int_t i = 0; i < nmodules; i++ ) {
    THaDetMap::Module* d = fDetMap->GetModule( i );

    // Loop over all channels that have a hit.
    for( Int_t j = 0; j < evdata.GetNumChan( d->crate, d->slot ); j++) {
      Int_t chan = evdata.GetNextChan( d->crate, d->slot, j );
       if( chan > d->hi || chan < d->lo ) continue;    // Not one of my channels.
      Int_t nhit = evdata.GetNumHits(d->crate, d->slot, chan);
      if( nhit > 1 || nhit == 0 ) {
	OSSTREAM msg;
	msg << nhit << " hits on " << "ADC channel "
	    << d->crate << "/" << d->slot << "/" << chan;
	++fMessages[msg.str()];
	has_warning = true;
	if( nhit == 0 ) {
	  msg << ". Should never happen. Decoder bug. Call expert.";
	  Warning( Here(here), "Event %d: %s", evdata.GetEvNum(),
		   msg.str().c_str() );
	  continue;
	}
#ifdef WITH_DEBUG
	if( fDebug>0 ) {
	  Warning( Here(here), "Event %d: %s", evdata.GetEvNum(),
		   msg.str().c_str() );
	}
#endif
      }
      // Get the data. If multiple hits on a channel, take the first (ADC)
      Int_t data = evdata.GetData( d->crate, d->slot, chan, 0 );
      
      Int_t jchan = (d->reverse) ? d->hi - chan : chan-d->lo;
      //cout << " jchan = " <<  jchan << " " << chan << " " << d->hi << " chanmap = " << fChanMap[i][jchan]<< endl;
       if( jchan<0 || jchan>d->hi ) {
	Error( Here(here), "Illegal detector channel: %d", jchan );
	continue;
      }
#ifdef NDEBUG
      Int_t k = fChanMap[i][jchan];
#else
      Int_t k = fChanMap.at(i).at(jchan);
#endif
      //      cout << " k = " << k << " " << i << " " << jchan << endl;
      if (k<0) continue; // < 0 means channel is not used
      if( k>fNelem ) {
	Error( Here(here), "Bad array index: %d. Your channel map is "
	       "invalid. Data skipped.", k );
	continue;
      }
      
      // Copy the data and apply calibrations
      //cout << "k ? " << k << " data ? " << data << " fA[k] ???" << fA[k] << endl;
      if( (Float_t)data - fPed[k] > fThrADC ){
	fA[k]   = (Float_t)data;                   // ADC value
	fA_p[k] = (Float_t)data - fPed[k];         // ADC minus ped
	fA_c[k] = fA_p[k] * fGain[k];//FIXME: should be calibration...     // ADC corrected for simu: SH 6.64734e-01 PS 1.36180e+00
	fAsum_p += fA_p[k];             // Sum of ADC minus ped
	//if( fA_c[k] > 0.0 )
	fAsum_c += fA_c[k];             // Sum of ADC corrected
	fNhits++;
	fBlocks[k]->SetE(fA_c[k]);
      }else{
	fA[k] = fA_p[k] = fA_c[k] = 0.0;
      }
      //cout << " channel " << k << " data " << data << endl;
    }
  }
  if( has_warning )
    ++fNEventsWithWarnings;

#ifdef WITH_DEBUG
  if ( fDebug 3 ) {
    cout << endl << "Shower Detector " << GetPrefix() << ":" << endl;
    int ncol=3;
    for (int i=0; i<ncol; i++) {
      cout << "  Block  ADC  ADC_p  ";
    }
    cout << endl;

    for (int i=0; i<(fNelem+ncol-1)/ncol; i++ ) {
      for (int c=0; c<ncol; c++) {
	int ind = c*fNelem/ncol+i;
	if (ind < fNelem) {
	  cout << "  " << setw(3) << ind+1;
	  cout << "  "; WriteValue(fA[ind]);
	  cout << "  "; WriteValue(fA_p[ind]);
	  cout << "  ";
	} else {
	  //	  cout << endl;
	  break;
	}
      }
      cout << endl;
    }
  }
#endif
  return fNhits;
}

//_____________________________________________________________________________
Int_t SBSBBShower::CoarseProcess(TClonesArray& tracks)
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
  //return 0;
  if( fCoarseProcessed ) return 0;

  Int_t col, row;
  Int_t colmax=0, rowmax=0;
  Double_t  energy_max = 0.0;
  
  std::vector<Int_t> locrowmax, loccolmax;
  std::vector<Double_t> locenergy_max;
  
# if not defined(_WIN32)//Win32 compiler do not support variable as array size
  Double_t energyDep[fNcols][fNrows];
# else
  Double_t energyDep[100][100];
# endif
  
  Double_t energyTotal = 0.0;
  SBSBBShowerCluster cluster(fNclublk);
  
  
  //  for( col = 0; col < fNcols; col++ )
  //     {
  //       for( row = 0; row < fNrows; row++ )
  // 	{
  // 	  energyDep[col][row] = 0.0;
  // 	}
  //     }
  
  //  cout << "Energy Deposition:" << endl <<"___________________________________________________" << endl;
  for( row = 0; row < fNrows; row++ ){
    for( col = 0; col < fNcols; col++ ){
      energyDep[col][row] = fA_c[BlockColRowToNumber(col,row)]; 
      

      if( energyDep[col][row] < 0.0 ) 
	energyDep[col][row] = 0.0;
      energyTotal += energyDep[col][row];
      
      //cluster seeds
      if(energyDep[col][row]>fEmin){
	if(energyDep[col][row]>energy_max){
	  energy_max=energyDep[col][row];
	  colmax = col;
	  rowmax = row;
	}
	
	if(fDebug)cout << " col " << col << " row " << row << " block ? " << BlockColRowToNumber(col,row) << " Edep ? " << energyDep[col][row] << endl; 
	
	locrowmax.push_back(row);
	loccolmax.push_back(col);
	locenergy_max.push_back(energyDep[col][row]);
      }//end      
    }
    //      cout << endl;
  }
  
  /*
  for( row = 0; row < fNrows; row++ ){
    for( col = 0; col < fNcols; col++ ){
      // Main cluster (highest energy)
      if(energyDep[col][row]>fEmin){
	if(energyDep[col][row]>energy_max){
	  energy_max=energyDep[col][row];
	  colmax = col;
	  rowmax = row;
	}
	
	locrowmax.push_back(row);
	loccolmax.push_back(col);
	locenergy_max.push_back(energyDep[col][row]);
      }//end 
    }
  }
  */
  
  if(energy_max < fEmin){
    fCoarseProcessed = 1;
    return 0;
  }
    
  //clean the extrac cluster seeds
  for(size_t i = 0; i<locenergy_max.size(); i++){
    /*
    if(locrowmax[i]==rowmax && loccolmax[i]==colmax){
      locrowmax.erase(locrowmax.begin()+i);
      loccolmax.erase(loccolmax.begin()+i);
      locenergy_max.erase(locenergy_max.begin()+i);
    }
    */
    for(size_t j = 0; j<i; j++){
      /*
      if(locrowmax[j]==rowmax && loccolmax[j]==colmax){
	locrowmax.erase(locrowmax.begin()+j);
	loccolmax.erase(loccolmax.begin()+j);
	locenergy_max.erase(locenergy_max.begin()+j);
      }
      */
      if(abs(locrowmax[i]-locrowmax[j])<=fClusBlockRadX &&
	 abs(loccolmax[i]-loccolmax[j])<=fClusBlockRadY){
	if(locenergy_max[i]<locenergy_max[j]){
	  locrowmax.erase(locrowmax.begin()+i);
	  loccolmax.erase(loccolmax.begin()+i);
	  locenergy_max.erase(locenergy_max.begin()+i);
	}else{
	  locrowmax.erase(locrowmax.begin()+j);
	  loccolmax.erase(loccolmax.begin()+j);
	  locenergy_max.erase(locenergy_max.begin()+j);
	}
      }
    }
  }
  
  if(fDebug){
    cout << " col max ? " << colmax << " row max ? " << rowmax << " Emax ? " << energy_max << " Emin = " << fEmin << endl;
    cout << " secondary clusters seeds size " << locrowmax.size() << " " << loccolmax.size() << " " << locenergy_max.size() << endl;
  }
  
  Int_t i, j;//, k=0;
  Double_t energyClusterTotal = 0.0;
  //  Double_t energyClusterGreatest = 0.0;
  
  Int_t mnrow, mxrow, mncol, mxcol;
  
  for(size_t cls = 0; cls<locenergy_max.size(); cls++){
    energyClusterTotal = 0.0;
    
    mnrow=TMath::Max(locrowmax[cls]-fClusBlockRadX,0);
    mxrow=TMath::Min(locrowmax[cls]+fClusBlockRadX,fNrows-1);
    mncol=TMath::Max(loccolmax[cls]-fClusBlockRadY,0);
    mxcol=TMath::Min(loccolmax[cls]+fClusBlockRadY,fNcols-1);
    
    if(fDebug){
      cout << " Cluster: mnrow " << mnrow << " mxrow " << mxrow 
	   << " mncol " << mncol << " mxcol " << mxcol << endl;
    }
    for( i = mnrow; i <= mxrow; i++ ){
      for( j = mncol; j <= mxcol; j++){
	energyClusterTotal += energyDep[j][i];
	//fEblk[k] = energyDep[j][i];
	//k++;
      }
    }
    
    //Double_t energyCluster = energyClusterTotal;
    Double_t X, Y;
    
    //if( energyCluster < 0.0 ) return 0;
    
    X = fBlockX[BlockColRowToNumber(colmax, rowmax)];
    Y = fBlockY[BlockColRowToNumber(colmax, rowmax)];
    
    if(fDebug)cout << "Got a cluster! E = " << energyClusterTotal << " X = " << X << " Y = " << Y << endl;
    
    Double_t energyX = 0.0;
    Double_t energyY = 0.0;
    
    Int_t  blockcounter = 0;
    for( i = mnrow; i <= mxrow; i++ ){
      for( j = mncol; j <= mxcol; j++ ){
	if( (i >= 0 && i < fNrows ) && ( j >=0 && j < fNcols ) ){
	  energyX += energyDep[j][i]*fBlockX[BlockColRowToNumber(j,i)];
	  energyY += energyDep[j][i]*fBlockY[BlockColRowToNumber(j,i)];
	  
	  /*
	    if(i!=fBlocks[BlockColRowToNumber(j,i)]->GetRow()){
	    cout << "row " << i << " " << fBlocks[BlockColRowToNumber(j,i)]->GetRow() << endl;
	    }
	    if(j!=fBlocks[BlockColRowToNumber(j,i)]->GetCol()){
	    cout << "col " << j << " " << fBlocks[BlockColRowToNumber(j,i)]->GetCol() << endl;
	    }
	  */
	  cluster.AddBlock( fBlocks[BlockColRowToNumber(j,i)] );
	  if(fDebug){
	    cout << "Cluster " << &cluster << " Adding block row " 
		 << fBlocks[BlockColRowToNumber(j,i)]->GetRow() 
		 << " col " << fBlocks[BlockColRowToNumber(j,i)]->GetCol() 
		 << " E " << fBlocks[BlockColRowToNumber(j,i)]->GetE() 
		 << " new size " << cluster.GetSize() << " block " << blockcounter 
		 << " address " << cluster.GetBlock(blockcounter) << endl;
	  } 
	  blockcounter++;
	}
      }
    }
    
    X = energyX/energyClusterTotal;
    Y = energyY/energyClusterTotal;
    
    if(fDebug){
      cout << energyClusterTotal << " " << X << " " << Y 
	   << " " << fOrigin.X() << " " << fOrigin.Y() << " " << cluster.GetMult() << endl;
    }
    
    cluster.SetE( energyClusterTotal );
    cluster.SetX( X+fOrigin.X() );
    cluster.SetY( Y+fOrigin.Y() );
    //cluster.SetX( X );
    //cluster.SetY( Y );
    
    if(fMCdata){
      fEres = cluster.GetE();//-
      fXres = cluster.GetX();//-
      fYres = cluster.GetY();//-
    }
    
    AddCluster(cluster);
    if(fDebug)
      cout << "Added - we now have " << fNclust << endl;
  }
  
  locrowmax.clear();  
  loccolmax.clear();  
  locenergy_max.clear();  
  fCoarseProcessed = 1;
  return 0;
}

//_____________________________________________________________________________
Int_t SBSBBShower::FineProcess(TClonesArray& tracks)
{
  //return 0;
  if( fFineProcessed ) return 0;
  
  // Fine Shower processing.
  
  if(fDebug){
    cout << endl << fNclust << " clusters " << GetName()  << endl;
    for (int i=0;i<fNclust;i++) {
      cout << setw(2) << i << setw(7) << setprecision(1) 
	   << fClusters[i]->GetE() << setw(8) << setprecision(3) 
	   << fClusters[i]->GetX() << setw(8) << fClusters[i]->GetY() 
	   << setw(4) << fClusters[i]->GetMult() << endl;
    }
  }
  
  TVector3 clusterpoint;
  Double_t Emax = -10.0;
  for (int i=0;i<fNclust;i++) {
    if(fDebug){
      cout << fClusters[i] << " " << fClusters[i]->GetE() << " " 
	   << fClusters[i]->GetX() << " " << fClusters[i]->GetY() 
	   << " " << fClusters[i]->GetMult()  << endl; 
    }
    fE_cl[i] = fClusters[i]->GetE();
    fE_cl_corr[i] = fClusters[i]->GetE()*(gconst + gslope*acc_charge);
    fX_cl[i] = fClusters[i]->GetX();
    fY_cl[i] = fClusters[i]->GetY();
    fMult_cl[i] = fClusters[i]->GetMult();
 
    if(fMCdata){
      fE_cl_res[i] = 1.0 - fClusters[i]->GetE();// /
      fX_cl_res[i] = fClusters[i]->GetX();//-
      fY_cl_res[i] = fClusters[i]->GetY();//-
    }

    for(int j = 0; j<fMult_cl[i]; j++){
      if(fDebug){
	cout << " block " << j << " address " << fClusters[i]->GetBlock(j) << endl;
      }
      if(fClusters[i]->GetBlock(j)){
	if(fDebug){
	  cout << " E ? "<< fClusters[i]->GetBlock(j)->GetE() << endl;
	  cout << " col ? "<< fClusters[i]->GetBlock(j)->GetCol() << endl;
	  cout << " row ? " << fClusters[i]->GetBlock(j)->GetRow() << endl;
	}
	fNblk_cl[i][j] = BlockColRowToNumber(fClusters[i]->GetBlock(j)->GetCol(), fClusters[i]->GetBlock(j)->GetRow());
	fEblk_cl[i][j] = fClusters[i]->GetBlock(j)->GetE();
      }
    }
    
    if(fClusters[i]->GetE()>Emax){
      fE = fE_cl[i];
      fX = fX_cl[i];
      fY = fY_cl[i];
      fMult = fMult_cl[i];
      for(int j = 0; j<fMult; j++){
	if(fClusters[i]->GetBlock(j)){
	  fNblk[j] = fNblk_cl[i][j];
	  fEblk[j] = fEblk_cl[i][j];
	}
      }
      Emax = fE;
      
      if(fMCdata){
	fEres = 1.0-fE;// /
	fXres = fX;//-
	fYres = fY;//-
      }

    }
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

void SBSBBShower::AddCluster(SBSBBShowerCluster* clus) {
    fClusters[fNclust++]=clus;
}

void SBSBBShower::AddCluster(SBSBBShowerCluster& clus) {

    fClusters[fNclust] = new SBSBBShowerCluster(clus.GetNMaxBlocks());
    for(int k = 0; k<clus.GetSize(); k++){
      fClusters[fNclust]->AddBlock(clus.GetBlock(k));
    }
    fClusters[fNclust]->SetE(clus.GetE());
    fClusters[fNclust]->SetX(clus.GetX());
    fClusters[fNclust]->SetY(clus.GetY());
    if(fDebug){
      cout << " E " << fClusters[fNclust]->GetE() 
	   << " X " << fClusters[fNclust]->GetX()
	   << " Y " << fClusters[fNclust]->GetY()
	   << endl;
    }
    fClusters[fNclust++]->SetMult(clus.GetMult());
}

void SBSBBShower::RemoveCluster(int i) {
    fNclust--;
    for (int j=i;j<fNclust;j++) fClusters[j]=fClusters[j+1];
}

Int_t SBSBBShower::BlockColRowToNumber( Int_t col, Int_t row )
{
    return col*fNrows + row;
}

void SBSBBShower::LoadMCHitAt( Double_t x, Double_t y, Double_t E )
{  
  //ClearEvent();
  Clear();
  fNclust = 0;
  fClusters[fNclust] = new SBSBBShowerCluster(0);
  fClusters[fNclust]->SetE(E);
  fClusters[fNclust]->SetX(x);
  fClusters[fNclust]->SetY(y);
  fClusters[fNclust]->SetMult(0);   
  
  fE_cl[fNclust] = fClusters[fNclust]->GetE();
  fX_cl[fNclust] = fClusters[fNclust]->GetX();
  fY_cl[fNclust] = fClusters[fNclust]->GetY();
  fMult_cl[fNclust] = fClusters[fNclust]->GetMult();
  fNclust++;
}

