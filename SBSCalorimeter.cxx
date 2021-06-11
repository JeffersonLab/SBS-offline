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
  fMaxNclus(10), fConst(1.0), fSlope(0.0), fAccCharge(0.0), fDataOutputLevel(1000)
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
  //static const char* const here = "ReadDatabase()";

  FILE* file = OpenFile( date );
  if( !file ) return kFileError;

  // Some temporary variables which we'll use to read in the database
  //Float_t angle = 0.0;
  std::vector<Int_t> cluster_dim; // Default is 3x3 if none specified

  // Read mapping/geometry/configuration parameters
  fChanMapStart = 0;
  DBRequest config_request[] = {
    { "emin",         &fEmin,   kFloat, 0, false }, ///< minimum energy threshold
    { "cluster_dim",   &cluster_dim,   kIntV, 0, true }, ///< cluster dimensions (2D)
    { "nmax_cluster",   &fMaxNclus,   kInt, 0, true }, ///< maximum number of clusters to store
    { "const", &fConst, kFloat, 0, true }, ///< const from gain correction 
    { "slope", &fSlope, kFloat, 0, true }, ///< slope for gain correction 
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
    { "nclus", "Number of clusters meeting threshold", "fNclus" },
    { "e",      "Energy (MeV) of largest cluster",    "GetE()" },
    { "e_c",    "Corrected Energy (MeV) of largest cluster",    "GetECorrected()" },
    { "eblk",   "Energy (MeV) of highest energy block in the largest cluster",    "GetEBlk()" },
    { "eblk_c", "Corrected Energy (MeV) of highest energy block in the largest cluster",    "GetEBlkCorrected()" },
    { "rowblk", "Row of block with highest energy in the largest cluster",    "GetRow()" },
    { "colblk", "Col of block with highest energy in the largest cluster",    "GetCol()" },
    { "x",      "x-position (mm) of largest cluster", "GetX()" },
    { "y",      "y-position (mm) of largest cluster", "GetY()" },
    { "nblk",   "Number of blocks in the largest cluster",    "GetNblk()" },
    { "idblk",  "Logic number of block with highest energy in cluster",    "GetBlkID()" },
    {0}
  };
  err = DefineVarsFromList( vars, mode );
  if(err)
    return err;

  if(fDataOutputLevel>0) {
    // Store all blocks in main cluster
    RVarDef vars_raw[] = {
      { "clus_blk.e", "Energy of block in main cluster", "fMainclusblk.e"},
      { "clus_blk.e_c","Energy calibrated of block in main cluster", "fMainclusblk.e_c"},
      { "clus_blk.x", "x-position of block in main cluster", "fMainclusblk.x"},
      { "clus_blk.y", "y-position of block in main cluster", "fMainclusblk.y"},
      { "clus_blk.row","block row in main cluster",    "fMainclusblk.row" },
      { "clus_blk.col","block col in main cluster",    "fMainclusblk.col" },
      { "clus_blk.id","block number in main cluster",    "fMainclusblk.id" },
      { 0 }
    };
    err = DefineVarsFromList( vars_raw, mode );
    if(err)
      return err;
  }

  if(fDataOutputLevel>1) {
    // Store every cluster
    RVarDef vars_raw[] = {
      { "clus.e", "Energy of cluster", "fOutclus.e"},
      { "clus.e_c","Energy calibrated of cluster", "fOutclus.e_c"},
      { "clus.x", "x-position of cluster", "fOutclus.x"},
      { "clus.y", "y-position of cluster", "fOutclus.y"},
      { "clus.row","block row in cluster with highest energy",    "fOutclus.row" },
      { "clus.col","block col in cluster with highest energy",    "fOutclus.col" },
      { "clus.id","block number in cluster",    "fOutclus.id" },
      { "clus.nblk","number of blocks in cluster",    "fOutclus.n" },
      { "clus.eblk", "Energy of block with highest energy in cluster", "fOutclus.blk_e"},
      { "clus.eblk_c","Energy calibrated of block with highest energy in cluster", "fOutclus.blk_e_c"},
      { 0 }
    };
    err = DefineVarsFromList( vars_raw, mode );
    if(err)
      return err;
  }

  if(fDataOutputLevel>2) {
    // Store every single block
    RVarDef vars_raw[] = {
      { "blk.e", "Energy of block", "fOutblk.e"},
      { "blk.e_c","Energy calibrated of block", "fOutblk.e_c"},
      { "blk.x", "x-position of block", "fOutblk.x"},
      { "blk.y", "y-position of block", "fOutblk.y"},
      { "blk.row","block row",    "fOutblk.row" },
      { "blk.col","block col",    "fOutblk.col" },
      { "blk.id","block number",    "fOutblk.id" },
      { 0 }
    };
    err = DefineVarsFromList( vars_raw, mode );
    if(err)
      return err;
  }
  return err;
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

  // Now, find as many clusters that meet the minimum energy
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
    SBSCalorimeterCluster *clus = fClusters[0];
    fMainclus.e.push_back(clus->GetE());
    fMainclus.e_c.push_back(clus->GetE()*(fConst + fSlope*fAccCharge));
    fMainclus.x.push_back(clus->GetX());
    fMainclus.y.push_back(clus->GetY());
    fMainclus.n.push_back(clus->GetMult());
    fMainclus.blk_e.push_back(clus->GetEblk());
    fMainclus.blk_e_c.push_back(clus->GetEblk()*(fConst + fSlope*fAccCharge));
    fMainclus.id.push_back(clus->GetMaxElement()->GetID());
    fMainclus.row.push_back(clus->GetRow());
    fMainclus.col.push_back(clus->GetCol());

    if(fDataOutputLevel > 0 ) {
      for( SBSElement *blk : clus->GetElements() ) {
        fMainclusblk.e.push_back(blk->GetE());
        fMainclusblk.e.push_back(blk->GetE()*(fConst + fSlope*fAccCharge));
        fMainclusblk.x.push_back(blk->GetX());
        fMainclusblk.y.push_back(blk->GetY());
        fMainclusblk.row.push_back(blk->GetRow());
        fMainclusblk.col.push_back(blk->GetCol());
        fMainclusblk.id.push_back(blk->GetID());
      }
    }
  }

  // Should we store all the cluster info?
  if(fDataOutputLevel>1) {
    // Now store the remaining clusters (up to fMaxNclus)
    Int_t nres = TMath::Min(Int_t(fMaxNclus),Int_t(fClusters.size()));
    fOutclus.e.reserve(nres);
    fOutclus.e_c.reserve(nres);
    fOutclus.x.reserve(nres);
    fOutclus.y.reserve(nres);
    fOutclus.n.reserve(nres);
    fOutclus.blk_e.reserve(nres);
    fOutclus.blk_e_c.reserve(nres);
    fOutclus.row.reserve(nres);
    fOutclus.col.reserve(nres);
    fOutclus.id.reserve(nres);

    int nclus = 0;
    for(SBSCalorimeterCluster *cluster : fClusters) {
      if(nclus < fMaxNclus) { // Keep adding them until we reach fMaxNclus
        fOutclus.e.push_back(cluster->GetE());
        fOutclus.e_c.push_back(cluster->GetE()*(fConst + fSlope*fAccCharge));
        fOutclus.x.push_back(cluster->GetX());
        fOutclus.y.push_back(cluster->GetY());
        fOutclus.n.push_back(cluster->GetMult());
        fOutclus.blk_e.push_back(cluster->GetEblk());
        fOutclus.blk_e_c.push_back(cluster->GetEblk()*(fConst + fSlope*fAccCharge));
        fOutclus.row.push_back(cluster->GetRow());
        fOutclus.col.push_back(cluster->GetCol());
        fOutclus.id.push_back(cluster->GetMaxElement()->GetID());
      }
      nclus++;
    }
  }

  // Should we store the individual block info?
  if(fDataOutputLevel > 2 ) {
    for( SBSElement *blk : fElements ) {
      fMainclusblk.e.push_back(blk->GetE());
      fOutblk.e.push_back(blk->GetE()*(fConst + fSlope*fAccCharge));
      fOutblk.x.push_back(blk->GetX());
      fOutblk.y.push_back(blk->GetY());
      fOutblk.row.push_back(blk->GetRow());
      fOutblk.col.push_back(blk->GetCol());
      fOutblk.id.push_back(blk->GetID());
    }
  }

  fFineProcessed = 1;
  return 0;
}

void SBSCalorimeter::ClearOutputVariables()
{
  ClearCaloOutput(fMainclus);
  ClearCaloOutput(fMainclusblk);
  ClearCaloOutput(fOutclus);
  ClearCaloOutput(fOutblk);
  fNclus = 0;
}


void SBSCalorimeter::ClearCaloOutput(SBSCalorimeterOutput &out)
{
  out.e.clear();
  out.e_c.clear();
  out.x.clear();
  out.y.clear();
  out.row.clear();
  out.col.clear();
  out.n.clear();
  out.blk_e.clear();
  out.blk_e_c.clear();
  out.id.clear();
}
