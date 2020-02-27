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
  fBlocks(0), fClusters(0) //, fBlkGrid(0)
{
    // Constructor.
    fCoarseProcessed = 0;
    fFineProcessed = 0;
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

  vector<Int_tdetmap, chanmap;
  vector<Double_txy, dxy;
  Int_t ncols, nrows;
  // Read mapping/geometry/configuration parameters
  DBRequest config_request[] = {
    { "detmap",       &detmap,  kIntV },
    { "chanmap",      &chanmap, kIntV,    0, 1 },
    { "ncols",        &ncols,   kInt },
    { "nrows",        &nrows,   kInt },
    { "xy",           &xy,      kDoubleV, 2 },  // center pos of block 1
    { "dxdy",         &dxy,     kDoubleV, 2 },  // dx and dy block spacings
    { "emin",         &fEmin,   kDouble },
    { 0 }
  };
//cout << " call Loaddb" << endl;
  err = LoadDB( file, date, config_request, fPrefix );


//cout << " (nrows ' " << nrows << " " <<  ncols << endl;
  // Sanity checks
  if( !err && (nrows <= 0 || ncols <= 0) ) {
    Error( Here(here), "Illegal number of rows or columns: %d %d. Must be 0. "
	   "Fix database.", nrows, ncols );
    err = kInitError;
  }

  Int_t nelem = ncols * nrows; 
  Int_t nclbl = TMath::Min( 3, nrows ) * TMath::Min( 3, ncols );


//cout << " nelem= " << nelem << " " << nclbl << " " << fNelem << endl;

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
//cout << "Set up the new channel map" << endl;
      Int_t nmodules = fDetMap->GetSize();
cout << "Set up the new channel map" << nmodules << endl;
      assert( nmodules 0 );
fChanMap.resize(nmodules);
      for( Int_t i=0, k=0; i < nmodules && !err; i++ ) {
	THaDetMap::Module* module = fDetMap->GetModule(i);
	Int_t nchan = module->hi - module->lo + 1;
        cout << " nchan = " << nchan << endl;
	if( nchan 0 ) {
	  fChanMap.at(i).resize(nchan);
	  for( Int_t j=0; j<nchan; ++j ) {
	    cout << " k = " << k << " " << nchan*nmodules << endl;
	    assert( k < nmodules*nchan );
	    fChanMap.at(i).at(j) = chanmap.empty() ? k : chanmap[k]-1;
	  	    cout << " k = " << k << " " << nchan*nmodules << " " << chanmap[k] << " " << fChanMap.at(i).at(j)<< endl;
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

    fIsInit = true;
  }

//cout << " define geometry" << endl;
  // Compute block positions
  for( int c=0; c<ncols; c++ ) {
    for( int r=0; r<nrows; r++ ) {
      int k = nrows*c + r;
      // Units are meters
      fBlockX[k] = xy[0] + r*dxy[0];
      fBlockY[k] = xy[1] + c*dxy[1];
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
    { "pedestals",    fPed,   kFloat, nval, 1 },
    { "gains",        fGain,  kFloat, nval, 1 },
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
fMaxNClust=10;
  fClusters = new SBSBBShowerCluster*[fMaxNClust];
  fBlocks = new SBSShowerBlock*[fNelem];
//cout << " retruning" << endl;
  return kOK;
  /*
    const int LEN = 100;
    char buf[LEN];
    Int_t nelem, ncols, nrows, nclbl;

    // clean out the old
    RemoveVariables();

    // Read data from database
    FILE* fi = OpenFile( date );
    if( !fi ) return kFileError;

    // Blocks, rows, max blocks per cluster
    fgets ( buf, LEN, fi ); fgets ( buf, LEN, fi );          
    fscanf ( fi, "%d%d", &ncols, &nrows );  

    nelem = ncols * nrows;
    nclbl = TMath::Min( 3, nrows ) * TMath::Min( 3, ncols );
    // Reinitialization only possible for same basic configuration 
    if( fIsInit && (nelem != fNelem || nclbl != fNclublk) ) {
        Error( Here(here), "Cannot re-initalize with different number of blocks or "
            "blocks per cluster. Detector not re-initialized." );
        fclose(fi);
        return kInitError;
    }

    if( nrows <= 0 || ncols <= 0 || nclbl <= 0 ) {
        Error( Here(here), "Illegal number of rows or columns: "
            "%d %d", nrows, ncols );
        fclose(fi);
        return kInitError;
    }
    fNelem = nelem;
    fNrows = nrows;
    fNcols = ncols;
    fNclublk = nclbl;

    // Clear out the old detector map before reading a new one
    UShort_t mapsize = fDetMap->GetSize();
    delete [] fNChan;
    if( fChanMap ) {
        for( UShort_t i = 0; i<mapsize; i++ )
            delete [] fChanMap[i];
    }
    delete [] fChanMap;
    fDetMap->Clear();

    // Read detector map

    fgets ( buf, LEN, fi ); fgets ( buf, LEN, fi );
    while (1) {
        Int_t crate, slot, first, last;
        fscanf ( fi,"%d%d%d%d", &crate, &slot, &first, &last );
        fgets ( buf, LEN, fi );
        if( crate < 0 ) break;
        if( fDetMap->AddModule( crate, slot, first, last ) < 0 ) {
            Error( Here(here), "Too many DetMap modules (maximum allowed - %d).", 
                THaDetMap::kDetMapSize);
            fclose(fi);
            return kInitError;
        }
    }

    // Set up the new channel map
    mapsize = fDetMap->GetSize();
    if( mapsize == 0 ) {
        Error( Here(here), "No modules defined in detector map.");
        fclose(fi);
        return kInitError;
    }

    fNChan = new UShort_t[ mapsize ];
    fChanMap = new UShort_t*[ mapsize ];
    for( UShort_t i=0; i < mapsize; i++ ) {
        THaDetMap::Module* module = fDetMap->GetModule(i);
        fNChan[i] = module->hi - module->lo + 1;
        if( fNChan[i] > 0 )
            fChanMap[i] = new UShort_t[ fNChan[i] ];
        else {
            Error( Here(here), "No channels defined for module %d.", i);
            delete [] fNChan; fNChan = NULL;
            for( UShort_t j=0; j<i; j++ )
                delete [] fChanMap[j];
            delete [] fChanMap; fChanMap = NULL;
            fclose(fi);
            return kInitError;
        }
    }
    // Read channel map
    fgets ( buf, LEN, fi );
    for ( UShort_t i = 0; i < mapsize; i++ ) {
        for ( UShort_t j = 0; j < fNChan[i]; j++ ) 
            fscanf (fi, "%hu", *(fChanMap+i)+j ); 
        fgets ( buf, LEN, fi );
    }
    fgets ( buf, LEN, fi );

    Float_t x,y,z;
    fscanf ( fi, "%f%f%f", &x, &y, &z );               // Detector's X,Y,Z coord
    fOrigin.SetXYZ( x, y, z );
    fgets ( buf, LEN, fi ); fgets ( buf, LEN, fi );
    fscanf ( fi, "%lf%lf%lf", fSize, fSize+1, fSize+2 );  // Sizes of det in X,Y,Z
    fdZ = TMath::Abs(fSize[2]);
    fgets ( buf, LEN, fi ); fgets ( buf, LEN, fi );

    Float_t angle;
    fscanf ( fi, "%f", &angle );                       // Rotation angle of det
    fgets ( buf, LEN, fi ); fgets ( buf, LEN, fi );
    const Double_t degrad = TMath::Pi()/180.0;
    tan_angle = TMath::Tan(angle*degrad);
    sin_angle = TMath::Sin(angle*degrad);
    cos_angle = TMath::Cos(angle*degrad);

    DefineAxes(angle*degrad);

    // Dimension arrays
    if( !fIsInit ) {
        fBlockX = new Float_t[ fNelem ];
        fBlockY = new Float_t[ fNelem ];
        fPed    = new Float_t[ fNelem ];
        fGain   = new Float_t[ fNelem ];

        // Per-event data
        fA    = new Float_t[ fNelem ];
        fA_p  = new Float_t[ fNelem ];
        fA_c  = new Float_t[ fNelem ];
        fNblk = new Int_t[ fNclublk ];
        fEblk = new Float_t[ fNclublk ];

        fIsInit = true;
    }

    fscanf ( fi, "%f%f", &x, &y );                     // Block 1 center position
    fgets ( buf, LEN, fi ); fgets ( buf, LEN, fi );
    Float_t dx, dy;
    fscanf ( fi, "%f%f", &dx, &dy );                   // Block spacings in x and y
    fdX = TMath::Abs(dx); 
    fdY = TMath::Abs(dy);
    fgets ( buf, LEN, fi ); fgets ( buf, LEN, fi );
    fscanf ( fi, "%f", &fEmin );                       // Emin thresh for center
    fgets ( buf, LEN, fi ); fgets ( buf, LEN, fi );
    fscanf ( fi, "%i", &fMaxNClust );                   // Max number of clusters
    fgets ( buf, LEN, fi ); fgets ( buf, LEN, fi );

    // Read calibrations.
    // Before doing this, search for any date tags that follow, and start reading from
    // the best matching tag if any are found. If none found, but we have a configuration
    // string, search for it.
    if( SeekDBdate( fi, date ) == 0 && fConfig.Length() > 0 && 
        SeekDBconfig( fi, fConfig.Data() )) {}

    fgets ( buf, LEN, fi );  
    // Crude protection against a missed date/config tag
    if( buf[0] == '[' ) fgets ( buf, LEN, fi );

    // Read ADC pedestals and gains (in order of logical channel number)
    for (int j=0; j<fNelem; j++)
        fscanf (fi,"%f",fPed+j);
    fgets ( buf, LEN, fi ); fgets ( buf, LEN, fi );
    for (int j=0; j<fNelem; j++) 
        fscanf (fi, "%f",fGain+j);
    fgets ( buf, LEN, fi ); fgets ( buf, LEN, fi );  //new line added
    
    // Compute block positions and creates blocks array
    fBlkGrid = new SBSShowerBlock**[fNrows];
    for (int i=0;i<fNrows;i++) fBlkGrid[i] = new SBSShowerBlock*[fNcols];
    fClusters = new SBSBBShowerCluster*[fMaxNClust];
    fBlocks = new SBSShowerBlock*[fNelem];
    for( int c=0; c<ncols; c++ ) {
        for( int r=0; r<nrows; r++ ) {
            int k = nrows*c + r;
            fBlockX[k] = x + r*dx;                         // Units are meters
            fBlockY[k] = y + c*dy;
            SBSShowerBlock* block = 
                new SBSShowerBlock(fBlockX[k],fBlockY[k],fPed[k],fGain[k],r,c);
            fBlocks[k]=block;
            fBlkGrid[r][c]=fBlocks[k];
        }
    }


    fE = new Float_t[fMaxNClust];
    fE_cl_corr = new Float_t[fMaxNClust];
    fX = new Float_t[fMaxNClust];
    fY = new Float_t[fMaxNClust]; 
    //fXtarg = new Float_t[fMaxNClust];
    //fYtarg = new Float_t[fMaxNClust];
    //fZtarg = new Float_t[fMaxNClust];
    fMult = new Int_t[fMaxNClust];
    
    //Read parameters for correcting gain drop.
    fscanf ( fi, "%f%f%f", &gconst, &gslope, &acc_charge );  
    //    cout << "gconst : " << gconst << endl;
    //    cout << "gslope : " << gslope <<endl;
    //    cout << "acc_charge : " << acc_charge <<endl;

    fclose(fi);


    //THaBBe *bb = dynamic_cast<THaBBe *> (GetApparatus());
    //if( bb )
    //{
    //    fDetToTarg = bb->GetDetectorToTarg();
    //    fDetOffset = bb->GetDetectorOffset();
    //}
    //else
    //{
    //    fDetToTarg.SetToIdentity();
    //    fDetOffset.SetXYZ( 0.0, 0.0, 0.0 );
    //}

    return kOK;
    */
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
        { "e",      "Energy (MeV) of largest cluster",    "fE" },
        { "e_c",    "Corrected Energy (MeV) of largest cluster",    "fE_cl_corr" },
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
        { 0 }
    };
    return DefineVarsFromList( vars, mode );
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
  cout << 0 << endl;
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
    for(int i = 0; i<fNclusters; i++){
      delete [] fNblk_cl[i]; fNblk_cl = 0;
      delete [] fEblk_cl[i]; fEblk_cl = 0;
    }
    delete [] fNblk_cl;
    delete [] fEblk_cl;
    delete [] fE_cl_corr; fE_cl_corr = 0;
    
}

//_____________________________________________________________________________
inline
void SBSBBShower::ClearEvent()
{
    // Reset all local data to prepare for next event.

    fCoarseProcessed = 0;
    fFineProcessed = 0;

    const int lsh = fNelem*sizeof(Float_t);
    const int lshh = fMaxNClust*sizeof(Float_t);
    const int lsc = fNclublk*sizeof(Float_t);
    const int lsi = fNclublk*sizeof(Int_t);
    const int lsj = fMaxNClust*sizeof(Int_t);

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
    memset( fNblk_cl, 0, lsi*fMaxNClust );
    memset( fEblk_cl, 0, lsc*fMaxNClust );
    memset( fMult_cl, 0, lsj );
    fE = 0.0;
    fX = 0.0;
    fY = 0.0;
    fMult = 0.0;
    fTRX = 0.0;
    fTRY = 0.0;
  
    for (int i=0;i<fNelem;i++) 
      fBlocks[i]->ClearEvent();
    for (int i=0;i<fMaxNClust;i++) 
      fBlocks[i]->ClearEvent();
     
    for (int i=0;i<fNrows;i++)
      for (int j=0;j<fNcols;j++)  
	fBlkGrid[i][j]->ClearEvent();


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
       if( chan d->hi || chan < d->lo ) continue;    // Not one of my channels.
      Int_t nhit = evdata.GetNumHits(d->crate, d->slot, chan);
      if( nhit 1 || nhit == 0 ) {
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
      if (k+1==999) continue; // 999 means channel is not used
      if( k<0 || k>fNelem ) {
	Error( Here(here), "Bad array index: %d. Your channel map is "
	       "invalid. Data skipped.", k );
	continue;
      }
      
      // Copy the data and apply calibrations
      fA[k]   = data;                   // ADC value
      fA_p[k] = data - fPed[k];         // ADC minus ped
      fA_c[k] = fA_p[k] * fGain[k];     // ADC corrected
      if( fA_p[k] 0.0 )
	fAsum_p += fA_p[k];             // Sum of ADC minus ped
      if( fA_c[k] 0.0 )
	fAsum_c += fA_c[k];             // Sum of ADC corrected
      fNhits++;
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
    return 0;
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
    SBSBBShowerCluster cluster(9);

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
Int_t SBSBBShower::FineProcess(TClonesArray& tracks)
{
     return 0;
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
    Double_t Emax = -10.0;
    for (int i=0;i<fNclust;i++) {
        //    cout << fClusters[i]->GetE() << " " << fClusters[i]->GetX() << " " << fClusters[i]->GetY() <<fClusters[i]->GetMult()  << endl; 
      fE_cl[i] = fClusters[i]->GetE();
      fE_cl_corr[i] = fClusters[i]->GetE()*(gconst + gslope*acc_charge);
      fX_cl[i] = fClusters[i]->GetX();
      fY_cl[i] = fClusters[i]->GetY();
      fMult_cl[i] = fClusters[i]->GetMult();
      for(int j = 0; j<fMult_cl[i]; j++){
	if(fClusters[i]->GetBlock(j)){
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
    fClusters[fNclust]->SetE(clus.GetE());
    fClusters[fNclust]->SetX(clus.GetX());
    fClusters[fNclust]->SetY(clus.GetY());
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
    ClearEvent();
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
