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
  fMaxNclus(20), fConst(1.0), fSlope(0.0), fAccCharge(0.0), fDataOutputLevel(1000)
{
  // Constructor.
  fTmax = 1000.0; // 1000 ns maximum arrival time difference with seed to be in cluster
  fEmin = 0.001; // 1 MeV minimum energy to be in cluster (Hit threshold)  
  fEmin_clusSeed = 0.001; // 1 MeV minimum energy to be the seed of a cluster
  fEmin_clusTotal = 0.001; // Minimum total cluster energy is 1 MeV
  fXmax_dis = .30; // Maximum X (m) distance from cluster center to be included in cluster
  fYmax_dis = .30; // Maximum Y (m) distance from cluster center to be included in cluster
  fRmax_dis = .30; // Maximum Radius (m) from cluster center to be included in cluster
  fBestClusterIndex = -1;
  fClusters.reserve(20);
  ftdctw.reserve(2);
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
}

///////////////////////////////////////////////////////////////////////////////
/// Read SBSCalorimeter Database
Int_t SBSCalorimeter::ReadDatabase( const TDatime& date )
{
  // Call parent class's ReadDatabase first
  Int_t err = SBSGenericDetector::ReadDatabase(date);
  if(err)
    return err;
  fIsInit = false;

  // Read this detector's parameters from the database file 'fi'.
  // This function is called by THaDetectorBase::Init() once at the
  // beginning of the analysis.
  // 'date' contains the date/time of the run being analyzed.

  // We can use this name here for logs
  //static const char* const here = "ReadDatabase()";

  FILE* file = OpenFile( date );
  if( !file ) return kFileError;

  // Some temporary variables which we'll use to read in the database
  //Double_t angle = 0.0;
  std::vector<Int_t> cluster_dim; // Default is 3x3 if none specified

  // Read mapping/geometry/configuration parameters
  DBRequest config_request[] = {
    { "emin",         &fEmin,   kDouble, 0, true }, ///< minimum energy threshold
    { "tmax",         &fTmax,   kDouble, 0, true }, ///< maximum time difference for block
    { "tdc.tw",       &ftdctw,   kDoubleV, 0, true }, ///< tdc dt = P0/(E^P1) timewalk fit parameters
    { "emin_clSeed", &fEmin_clusSeed, kDouble, 0, true }, ///< minimum cluster seed energy
    { "emin_clTotal", &fEmin_clusTotal, kDouble, 0, true }, ///< minimum total cluster energy
    { "cluster_dim",   &cluster_dim,   kIntV, 0, true }, ///< cluster dimensions (2D)
    { "nmax_cluster",   &fMaxNclus,   kInt, 0, true }, ///< maximum number of clusters to store
    { "const", &fConst, kDouble, 0, true }, ///< const from gain correction 
    { "slope", &fSlope, kDouble, 0, true }, ///< slope for gain correction 
    { "Rmax_dis", &fRmax_dis, kDouble, 0, true }, ///< slope for gain correction 
    { "acc_charge", &fAccCharge, kDouble, 0, true }, ///< accumulated charge
    { 0 } ///< Request must end in a NULL
  };
  err = LoadDB( file, date, config_request, fPrefix );
  if(err) {
    fclose(file);
    return err;
  }

  // tdc timewalk correction: ensure that if no parameters passed in db data isn't effected
  if(ftdctw.empty()){
    ftdctw.push_back(0.);
    ftdctw.push_back(0.);
  }

  // Compute the max possible cluster size (which at most should be
  // cluster_dim x cluster_dim)
  if(cluster_dim.empty()) {
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

  //
  std::vector<Double_t> xpos,ypos;
  std::vector<Double_t> trigtoFADCratio;
  std::vector<DBRequest> vr;
  vr.push_back({ "xpos", &xpos,    kDoubleV, 0, 1 });
  vr.push_back({ "ypos", &ypos,    kDoubleV, 0, 1 });
  vr.push_back({ "trigtoFADCratio", &trigtoFADCratio,    kDoubleV, 0, 1 });
  vr.push_back({0});
  err = LoadDB( file, date, vr.data(), fPrefix );
  fclose(file);
  if(err)
    return err;

  if (!trigtoFADCratio.empty()) {
    if ((int)trigtoFADCratio.size() == fNelem) {
      for (Int_t ne=0;ne<fNelem;ne++) {
	SBSElement* blk= fElements[ne];
	if (WithADC() && fModeADC == SBSModeADC::kWaveform) {
	  SBSData::Waveform *wave = blk->Waveform();
	  Double_t gain = wave->GetGain();
	  wave->SetGain(gain*trigtoFADCratio[ne]);
	  wave->SetTrigCal(trigtoFADCratio[ne]);
	}
	if (WithADC() && fModeADC == SBSModeADC::kADC) {
	  Double_t gain=blk->ADC()->GetGain();
	  blk->ADC()->SetGain(gain*trigtoFADCratio[ne]);
	  blk->ADC()->SetTrigCal(trigtoFADCratio[ne]);
	}
      }
    } else {
      std::cout << " trigtoFADCratio vector too small " << trigtoFADCratio.size() << " # of elements =" << fNelem << std::endl;
    }
  }
  //
  if (!xpos.empty()) {
    if ((int)xpos.size() == fNelem) {
      for (Int_t ne=0;ne<fNelem;ne++) {
	fElements[ne]->SetX(xpos[ne]);
      }
    } else {
      std::cout << "  vector too small " << xpos.size() << " # of elements =" << fNelem << std::endl;
    }
  }
  //
  if (!ypos.empty()) {
    if ((int)ypos.size() == fNelem) {
      for (Int_t ne=0;ne<fNelem;ne++) {
	fElements[ne]->SetY(ypos[ne]);
      }
    } else {
      std::cout << " ypos vector too small " << ypos.size() << " # of elements =" << fNelem << std::endl;
    }
  }

  // All is well that ends well
  fIsInit = true;
  return kOK;
}

//_____________________________________________________________________________
Int_t SBSCalorimeter::DefineVariables( EMode mode )
{
  if( mode == kDefine && fIsSetup ) return kOK;
  Int_t err = SBSGenericDetector::DefineVariables(mode);
  if(err)
    return err;

  // Most of these variables were used previously, and so I'll leave
  // them here to preserve old analysis macros people may have written.
  // This includes things like fE, fNblk, fE_c, etc...
  RVarDef vars[] = {
    { "nclus", "Number of clusters meeting threshold", "fNclus" },
    { "e",      "Energy (MeV) of largest cluster",    "GetE()" },
    //{ "e_c",    "Corrected Energy (MeV) of largest cluster",    "GetECorrected()" },
    { "againblk",    "ADC gain coeff. of highest energy block in the largest cluster",  "GetAgain()" },
    { "atimeblk", "ADC time of highest energy block in the largest cluster", "GetAtime()" },
    { "tdctimeblk", "TDC time of highest energy block in the largest cluster", "GetTDCtime()" },
    { "tdctimeblk_tw", "TDC time of highest energy block in the largest cluster, timewalk corrected", "GetTDCtimeTW()" },
    { "eblk",   "Energy (MeV) of highest energy block in the largest cluster",    "GetEBlk()" },
    //{ "eblk_c", "Corrected Energy (MeV) of highest energy block in the largest cluster",    "GetEBlkCorrected()" },
    { "rowblk", "Row of block with highest energy in the largest cluster",    "GetRow()" },
    { "colblk", "Col of block with highest energy in the largest cluster",    "GetCol()" },
    { "x",      "x-position (mm) of largest cluster", "GetX()" },
    { "y",      "y-position (mm) of largest cluster", "GetY()" },
    { "nblk",   "Number of blocks in the largest cluster",    "GetNblk()" },
    { "idblk",  "Logic number of block with highest energy in cluster",    "GetBlkID()" },
    { "index", "Index of best cluster in the array of all clusters", "fBestClusterIndex" },
    {0}
  };
  err = DefineVarsFromList( vars, mode );
  if(err)
    return err;

  RVarDef vars_gb[] = {
    { "goodblock.e", "Energy of good blocks", "fGoodBlocks.e"},
    //     { "goodblock.tdc", "TDC time of good blocks", "fGoodBlocks.TDCTime"},
    { "goodblock.atime", "Energy of good blocks", "fGoodBlocks.ADCTime"},
    { "goodblock.tdctime", "TDC time of good blocks", "fGoodBlocks.TDCTime"},
    { "goodblock.tdctime_tw", "TDC time of good blocks, timewalk corrected", "fGoodBlocks.TDCTimeTW"},
    { "goodblock.row", "Row of good blocks", "fGoodBlocks.row"},
    { "goodblock.col", "Col of good blocks", "fGoodBlocks.col"},
    { "goodblock.x", "x pos (m) of good blocks", "fGoodBlocks.x"},
    { "goodblock.y", "y pos (m) of good blocks", "fGoodBlocks.y"},
    { "goodblock.id", "Element ID of good blocks", "fGoodBlocks.id"},
    { "goodblock.cid", "Cluster ID of good blocks", "fGoodBlocks.cid"},
    {0}
  };
  err = DefineVarsFromList( vars_gb, mode );
  if(err)
    return err;

  if(fDataOutputLevel>0) {
    // Store all blocks in main cluster
    RVarDef vars_raw[] = {
      { "clus_blk.e", "Energy of block in main cluster", "fMainclusblk.e"},
      //{ "clus_blk.e_c","Energy calibrated of block in main cluster", "fMainclusblk.e_c"},
      { "clus_blk.again","ADC gain coeff. of block in main cluster", "fMainclusblk.again"},
      { "clus_blk.x", "x-position of block in main cluster", "fMainclusblk.x"},
      { "clus_blk.y", "y-position of block in main cluster", "fMainclusblk.y"},
      { "clus_blk.row","block row in main cluster",    "fMainclusblk.row" },
      { "clus_blk.atime","block ADC time in main cluster",    "fMainclusblk.atime" },
      { "clus_blk.tdctime","block TDC time in main cluster",    "fMainclusblk.tdctime" },
      { "clus_blk.tdctime_tw","block TDC time in main cluster, timewalk corrected",    "fMainclusblk.tdctime_tw" },
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
      { "clus.atime", "ADC time of cluster", "fOutclus.atime"},
      { "clus.tdctime", "TDC time of cluster", "fOutclus.tdctime"},
      { "clus.tdctime_tw", "TDC time of cluster, timewalk corrected", "fOutclus.tdctime_tw"},
      //{ "clus.e_c","Energy calibrated of cluster", "fOutclus.e_c"},
      { "clus.again","ADC gain coeff. of cluster", "fOutclus.again"},
      { "clus.x", "x-position of cluster", "fOutclus.x"},
      { "clus.y", "y-position of cluster", "fOutclus.y"},
      { "clus.row","block row in cluster with highest energy",    "fOutclus.row" },
      { "clus.col","block col in cluster with highest energy",    "fOutclus.col" },
      { "clus.id","block number in cluster",    "fOutclus.id" },
      { "clus.nblk","number of blocks in cluster",    "fOutclus.n" },
      { "clus.eblk", "Energy of block with highest energy in cluster", "fOutclus.blk_e"},
      // { "clus.eblk_c","Energy calibrated of block with highest energy in cluster", "fOutclus.blk_e_c"},
      { 0 }
    };
    err = DefineVarsFromList( vars_raw, mode );
    if(err)
      return err;
  }
  return err;
}

//_____________________________________________________________________________
void SBSCalorimeter::Clear( Option_t* opt )
{
  SBSGenericDetector::Clear(opt);
  ClearOutputVariables();
  DeleteContainer(fClusters);
  fGoodBlocks.clear();
  fBlockSet.clear();
  fBestClusterIndex = -1;
}

//_____________________________________________________________________________
Int_t SBSCalorimeter::MakeGoodBlocks()
{
  // Fill the fElements which have Good hit in Block "Cluster"
  SBSElement *blk = 0;
  for(Int_t k = 0; k < fNelem; k++) {  
    blk = fElements[k];
    Bool_t ADC_HasData=kFALSE;
    Int_t ADC_GoodHitIndex=-1;
    if(fModeADC != SBSModeADC::kWaveform) {
      ADC_HasData  = blk->ADC()->HasData();
      if (ADC_HasData) ADC_GoodHitIndex = blk->ADC()->GetGoodHitIndex();
    } else {
      SBSData::Waveform *wave = blk->Waveform();
      ADC_HasData = wave->HasData();
      if (ADC_HasData) ADC_GoodHitIndex = wave->GetGoodHitIndex();
    }
    if (WithADC() && ADC_HasData) {  
      if (ADC_GoodHitIndex != -1)  {
	fGoodBlocks.row.push_back(blk->GetRow());
	fGoodBlocks.col.push_back(blk->GetCol());
	fGoodBlocks.id.push_back(blk->GetID());
	fGoodBlocks.x.push_back(blk->GetX());
	fGoodBlocks.y.push_back(blk->GetY());
	fGoodBlocks.cid.push_back(-1); //Default cluster value. Will overwrite in FindClusters
	//
	//
	Double_t blkE = 0.;
	//
	if(fModeADC != SBSModeADC::kWaveform) {
	  const SBSData::PulseADCData &ahit = blk->ADC()->GetGoodHit();
	  blk->SetE(ahit.integral.val);
	  blkE = ahit.integral.val;
	  //blk->SetAintp(ahit.integral.val / blk->ADC()->GetGain());
	  blk->SetAgain(blk->ADC()->GetGain() / blk->ADC()->GetTrigCal());
	  blk->SetAtime(ahit.time.val);
	  fGoodBlocks.e.push_back(ahit.integral.val);
	  fGoodBlocks.ADCTime.push_back(ahit.time.val);
	} else {
	  SBSData::Waveform *wave = blk->Waveform();
	  blk->SetE(wave->GetIntegral().val);
	  blkE = wave->GetIntegral().val;
	  //blk->SetAintp(wave->GetIntegral().val / wave->GetGain());
	  blk->SetAgain(wave->GetGain() / wave->GetTrigCal());
	  blk->SetAtime(wave->GetTime().val);
	  fGoodBlocks.e.push_back(wave->GetIntegral().val);
	  fGoodBlocks.ADCTime.push_back(wave->GetTime().val);
	}
	if (WithTDC() && blk->TDC()->HasData() ) { 
	  const SBSData::TDCHit &hit = blk->TDC()->GetGoodHit();
	  fGoodBlocks.TDCTime.push_back(hit.le.val);
	  
	  //timewalk correction to goodblocks
	  Double_t tdc_tw = hit.le.val;
	  if( blkE>0 )
	    tdc_tw = hit.le.val - ftdctw[0]/pow( blkE, ftdctw[1] );

	  fGoodBlocks.TDCTimeTW.push_back(tdc_tw);

	  blk->SetTDCtime(hit.le.val);
	  blk->SetTDCtimeTW(tdc_tw);
	} else {
	  fGoodBlocks.TDCTime.push_back(-1000.);
	  fGoodBlocks.TDCTimeTW.push_back(-1000);

	  blk->SetTDCtime(-1000.);
	  blk->SetTDCtimeTW(-1000.);
	}
      }
    }
  }
  // Put good blocks in fBlockSet to use in FindCluster
  //  fBlockSet.reserve(fGoodBlocks.e.size());
  fBlockSet.clear();
  for (UInt_t nb=0;nb< fGoodBlocks.e.size();nb++) {
    SBSBlockSet c1 = {fGoodBlocks.e[nb],fGoodBlocks.x[nb],fGoodBlocks.y[nb],fGoodBlocks.row[nb],fGoodBlocks.col[nb],fGoodBlocks.id[nb],fGoodBlocks.cid[nb],fGoodBlocks.TDCTime[nb],fGoodBlocks.ADCTime[nb],kFALSE};
    if (fGoodBlocks.e[nb] > fEmin) fBlockSet.push_back(c1);
  }
  std::sort(fBlockSet.begin(), fBlockSet.end(), [](const SBSBlockSet& c1, const SBSBlockSet& c2) {
      return c1.e > c2.e;});
  //
  return fGoodBlocks.e.size();
}
//_____________________________________________________________________________
Int_t SBSCalorimeter::FindClusters()
{
  // fBlockSet is initially ordered by energy in MakeGoodblocks
  fNclus = 0;
  DeleteContainer(fClusters);

  Int_t NSize = fBlockSet.size();

  Int_t NGBSize = fGoodBlocks.id.size();

  fBestClusterIndex = -1;

  double Emax = 0.0;
  //We want the "best" cluster to be the one with the largest total energy (in general)
  
  while ( NSize != 0 )  {
    //sort all remaining blocks in terms of energy
    std::sort(fBlockSet.begin(), fBlockSet.end(), [](const SBSBlockSet& c1, const SBSBlockSet& c2) { return c1.e > c2.e;});
    
    fBlockSetIterator it = fBlockSet.begin();
    
    if ( (*it).e > fEmin_clusSeed ){
      Bool_t AddingBlocksToCluster = kTRUE;
      SBSElement *blk= fElements[(*it).id-fChanMapStart] ; 
      SBSCalorimeterCluster* cluster = new SBSCalorimeterCluster(fBlockSet.size(),blk);
      fClusters.push_back(cluster);
      fBlockSet.erase(it);
      NSize--;

      while (AddingBlocksToCluster) {
	Bool_t Close=kFALSE; //Check if block within user-def radius and time to add to cluster
	fBlockSetIterator it2 = fBlockSet.begin();

	while (!Close && (it2 < fBlockSet.end())) {
	  //next block in blockset variables
	  SBSElement *blk= fElements[(*it2).id-fChanMapStart]; 
	  Int_t blk_ID = blk->GetID();
	  Double_t blkx = blk->GetX();
	  Double_t blky = blk->GetY();

	  //current cluster primary variables
	  Int_t Index = fClusters.size()-1;
	  Double_t Clus_px = fClusters[Index]->GetX();
	  Double_t Clus_py = fClusters[Index]->GetY(); 
	  Int_t Clus_pblkid = fClusters[Index]->GetElemID();
	  Double_t Clus_ptime = fClusters[Index]->GetAtime();

	  Double_t Rad = sqrt( pow((Clus_px-blkx),2) + pow((Clus_py-blky),2) );
	  Double_t tDiff = blk->GetAtime()-Clus_ptime;

	  //set up vector to add cluster id to goodblocks
	  std::vector<Int_t> goodblock_ids = {Clus_pblkid};

	  //Check each additional block in blockset to see if it can be added to current cluster or pass to next
	  Close =( Rad<fRmax_dis && fabs(tDiff)<fTmax );
	  if (Close) {

	    fClusters[Index]->AddElement(blk);
	    goodblock_ids.push_back(blk_ID);
	  } else {	
       
	    ++it2;
	  }

	  //Add all relevant cluster ids to goodblocks. With simple assignments and similar order iterations,
	  //using nested loops for clarity and negligible impact to processing time
	  for( Int_t i=0; i<NGBSize; ++i ){
	    for( size_t j=0; j<goodblock_ids.size(); ++j ){
	      if( goodblock_ids[j]==fGoodBlocks.id[i] && fGoodBlocks.cid[i]==-1){
		fGoodBlocks.cid[i]=Index;
	      }
	    }
	  }

	}
	if (it2 == fBlockSet.end()) AddingBlocksToCluster = kFALSE;
	if (Close)   {
	  fBlockSet.erase(it2);
	  NSize--;
	}	
      }
      
      // Adding total cluster energy threshold
      if ( (cluster->GetE())<fEmin_clusTotal ) fClusters.pop_back();

      //After we are done adding blocks to the cluster, we check if this is the cluster
      //with the largest total energy:
      
      if( fClusters.size() == 1 || cluster->GetE() > Emax ){
	Emax = cluster->GetE();
	fBestClusterIndex = fClusters.size()-1;
      }
    } else break; //if the block doesn't pass fEmin_clusSeed, no remaining blocks will. Break.
  }

  //Fill main cluster variables here; this may get overridden by derived classes:
  if(!fClusters.empty()) {
    SBSCalorimeterCluster *clus = fClusters[fBestClusterIndex];

    SBSElement *pblk= clus->GetElement(0);
    Double_t pblkE= pblk->GetE();

    fMainclus.e.push_back(clus->GetE());
    //fMainclus.e_c.push_back(clus->GetE()*(fConst + fSlope*fAccCharge));
    fMainclus.again.push_back(clus->GetAgain());
    fMainclus.atime.push_back(clus->GetAtime());
    fMainclus.tdctime.push_back(clus->GetTDCtime());

    //Add timewalk correction to main cluster output
    Double_t tdc_tw = clus->GetTDCtime();

    if( pblkE!=0 )
      tdc_tw = clus->GetTDCtime() - ftdctw[0]/pow( pblkE, ftdctw[1] );

    fMainclus.tdctime_tw.push_back(tdc_tw);

    fMainclus.x.push_back(clus->GetX());
    fMainclus.y.push_back(clus->GetY());
    fMainclus.n.push_back(clus->GetMult());
    fMainclus.blk_e.push_back(clus->GetEblk());
    // fMainclus.blk_e_c.push_back(clus->GetEblk()*(fConst + fSlope*fAccCharge));
    fMainclus.id.push_back(clus->GetElemID());
    fMainclus.row.push_back(clus->GetRow());
    fMainclus.col.push_back(clus->GetCol());
  }

  //
  //
  fNclus = fClusters.size();
  return fNclus;
}
//_____________________________________________________________________________
Int_t SBSCalorimeter::FineProcess(TClonesArray& array)//tracks)
{
  Int_t err = SBSGenericDetector::FineProcess(array);
  if(err)
    return err;
  // Get information on the cluster with highest energy (useful even if
  // fMaxNclus is zero, i.e., storing no vector of clusters)

  if( !(fClusters.empty()) && fBestClusterIndex >= 0 && fBestClusterIndex < (Int_t)fClusters.size() ) {
    //if( !(fClusters.empty()) ) {
    SBSCalorimeterCluster *clus = fClusters[fBestClusterIndex];
 
    if(fDataOutputLevel > 0 ) {
      for(Int_t nc=0;nc<clus->GetMult();nc++ ) {
	SBSElement *blk= clus->GetElement(nc);
        fMainclusblk.e.push_back(blk->GetE());
        //fMainclusblk.e_c.push_back(blk->GetE()*(fConst + fSlope*fAccCharge));
	fMainclusblk.again.push_back(blk->GetAgain());
        fMainclusblk.atime.push_back(blk->GetAtime());        
        fMainclusblk.tdctime.push_back(blk->GetTDCtime());   

	//Add timewalk correction to main cluster block output
	Double_t tdc_tw;
	if( blk->GetE()>0 )
	  tdc_tw = blk->GetTDCtime() - ftdctw[0]/pow( blk->GetE(), ftdctw[1] );

	fMainclusblk.tdctime_tw.push_back(tdc_tw);
     
        fMainclusblk.x.push_back(blk->GetX());
        fMainclusblk.y.push_back(blk->GetY());
        fMainclusblk.row.push_back(blk->GetRow());
        fMainclusblk.col.push_back(blk->GetCol());
        fMainclusblk.id.push_back(blk->GetID());
      }
    }
  }
 
  // store all the cluster info
  if(fDataOutputLevel>1) {
    // Now store the remaining clusters (up to fMaxNclus)
    Int_t nres = TMath::Min(Int_t(fMaxNclus),Int_t(fClusters.size()));
    fOutclus.e.reserve(nres);
    //fOutclus.e_c.reserve(nres);
    fOutclus.again.reserve(nres);
    fOutclus.atime.reserve(nres);
    fOutclus.tdctime.reserve(nres);
    fOutclus.tdctime_tw.reserve(nres);
    fOutclus.x.reserve(nres);
    fOutclus.y.reserve(nres);
    fOutclus.n.reserve(nres);
    fOutclus.blk_e.reserve(nres);
    // fOutclus.blk_e_c.reserve(nres);
    fOutclus.row.reserve(nres);
    fOutclus.col.reserve(nres);
    fOutclus.id.reserve(nres);

    int nclus = 0;
    for( const auto* cluster: fClusters ) {
      if(nclus < fMaxNclus) { // Keep adding them until we reach fMaxNclus

        fOutclus.e.push_back(cluster->GetE());
        //fOutclus.e_c.push_back(cluster->GetE()*(fConst + fSlope*fAccCharge));
        fOutclus.again.push_back(cluster->GetAgain());
        fOutclus.atime.push_back(cluster->GetAtime());
        fOutclus.tdctime.push_back(cluster->GetTDCtime());

	//Add timewalk correction to all cluster output
	Double_t tdc_tw = cluster->GetTDCtime();
	if( cluster->GetEblk()>0 )
	  tdc_tw = cluster->GetTDCtime() - ftdctw[0]/pow( cluster->GetEblk(), ftdctw[1] );

	fOutclus.tdctime_tw.push_back(tdc_tw);

	//std::cout << "Calorimeter Name: " << GetName() << ", block energy: " << cluster->GetEblk() << " outclus TDC time: " << cluster->GetTDCtime() << ", modified: " << tdc_tw << std::endl;

        fOutclus.x.push_back(cluster->GetX());
        fOutclus.y.push_back(cluster->GetY());
        fOutclus.n.push_back(cluster->GetMult());
        fOutclus.blk_e.push_back(cluster->GetEblk());
        // fOutclus.blk_e_c.push_back(cluster->GetEblk()*(fConst + fSlope*fAccCharge));
        fOutclus.row.push_back(cluster->GetRow());
        fOutclus.col.push_back(cluster->GetCol());
        fOutclus.id.push_back(cluster->GetElemID());
      }
      nclus++;
    }
  }

  fFineProcessed = 1;
  return 0;
}

void SBSCalorimeter::ClearOutputVariables()
{
  fGoodBlocks.clear();
  ClearCaloOutput(fMainclus);
  ClearCaloOutput(fMainclusblk);
  ClearCaloOutput(fOutclus);
  fBestClusterIndex = 0;
  fNclus = 0;
}



void SBSCalorimeter::ClearCaloOutput(SBSCalorimeterOutput &out)
{
  out.e.clear();
  // out.e_c.clear();
  out.again.clear();
  out.atime.clear();
  out.tdctime.clear();
  out.tdctime_tw.clear();
  out.x.clear();
  out.y.clear();
  out.row.clear();
  out.col.clear();
  out.n.clear();
  out.blk_e.clear();
  // out.blk_e_c.clear();
  out.id.clear();
}
