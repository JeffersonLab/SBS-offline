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
#include "SBSManager.h"
#include "THaCrateMap.h"

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
  SBSGenericDetector(name,description,apparatus),
  fMaxNclus(10), fConst(1.0), fSlope(0.0), fAccCharge(0.0)
{
  // Constructor.
  fEmin = 1.0; // 1 MeV minimum
}

///////////////////////////////////////////////////////////////////////////////
/// Default Destructor
SBSCalorimeter::~SBSCalorimeter()
{
  // Destructor. Removes internal arrays and global variables.

  if( fIsSetup )
    RemoveVariables();
//  if( fIsInit ) {
//  }

  ClearEvent();
}

///////////////////////////////////////////////////////////////////////////////
/// Read SBSCalorimeter Database
Int_t SBSCalorimeter::ReadDatabase( const TDatime& date )
{
  // Call parent class's ReadDatabase first
  Int_t err = SBSGenericDetector::ReadDatabase(date);
  if(err)
    return err;

  // Read this detector's parameters from the database file 'fi'.
  // This function is called by THaDetectorBase::Init() once at the
  // beginning of the analysis.
  // 'date' contains the date/time of the run being analyzed.

  // We can use this name here for logs
  static const char* const here = "ReadDatabase()";

  FILE* file = OpenFile( date );
  if( !file ) return kFileError;

  // Some temporary variables which we'll use to read in the database
  std::vector<Float_t> xyz, dxyz;
  Float_t angle = 0.0;
  Int_t nchan_per_element = 1, model_in_detmap = 0;
  std::vector<Int_t> cluster_dim; // Default is 3x3 if none specified

  // Read mapping/geometry/configuration parameters
  fChanMapStart = 0;
  DBRequest config_request[] = {
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

  // Reinitialization only possible for same basic configuration
  if( !err ) {
    // Compute the max possible cluster size (which at most should be
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
    // TODO: Make this smarter, now that rows could be variable
    fNclubr = TMath::Min( cluster_dim[0], fNrows);
    fNclubc = TMath::Min( cluster_dim[1], fNcols[0] );
    fNclublk = fNclubr*fNclubc;
  }

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

  // All is well that ends well
  return kOK;
}

//_____________________________________________________________________________
Int_t SBSCalorimeter::DefineVariables( EMode mode )
{
  if( mode == kDefine && fIsSetup ) return kOK;
  // Initialize parent variables first
  Int_t err = SBSGenericDetector::DefineVariables(mode);
  if(err)
    return err;

  // Most of these variables were used previously, and so I'll leave
  // them here to preserve old analysis macros people may have written.
  // This includes things like fE, fNblk, fE_c, etc...
  RVarDef vars[] = {
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
    { 0 }
  };

  return DefineVarsFromList( vars, mode );
}

//_____________________________________________________________________________
void SBSCalorimeter::ClearEvent()
{
  SBSGenericDetector::ClearEvent();
  ClearOutputVariables();
  fClusters.clear();
}

Int_t SBSCalorimeter::CoarseProcess(TClonesArray& array)// tracks)
{
  // Call parent class, and exit if an error is encountered
  Int_t err = SBSGenericDetector::CoarseProcess(array);
  if(err)
    return err;

  // Pack simple data for output to the tree, and call CoarseProcess on blocks
  SBSElement *blk = 0;
  size_t nsamples;
  size_t idx;

  // Now, find as many clusters meet the minimum energy
  for(Int_t r = 0; r <= fNrows-fNclubr; r++) {
    for(Int_t c = 0; c <= fNcols[0]-fNclubc; c++) {
      for(Int_t l = 0; l < fNlayers; l++) {

        // Now perform the sum
        SBSCalorimeterCluster *clus = new SBSCalorimeterCluster(fNclublk);
        for(Int_t rr = 0; rr < fNclubr; rr++) {
          for(Int_t cc = 0; cc < fNclubc; cc++) {
            blk = fElementGrid[r+rr][c+cc][l];
            if(blk->GetE()>0)
              clus->AddElement(blk);
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
Int_t SBSCalorimeter::FineProcess(TClonesArray& array)//tracks)
{
  Int_t err = SBSGenericDetector::FineProcess(array);
  if(err)
    return err;


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
