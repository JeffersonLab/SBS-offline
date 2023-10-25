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
  fTmax = 1000.0; // 1000 ns maximum arrival time difference with seed to be in cluster
  fEmin = 0.001; // 1 MeV minimum energy to be in cluster (Hit threshold)  
  fEmin_clusSeed = 0.001; // 1 MeV minimum energy to be the seed of a cluster
  fEmin_clusTotal = 0.001; // Minimum total cluster energy is 1 MeV
  fXmax_dis = .30; // Maximum X (m) distance from cluster center to be included in cluster
  fYmax_dis = .30; // Maximum Y (m) distance from cluster center to be included in cluster
  fRmax_dis = .30; // Maximum Radius (m) from cluster center to be included in cluster
  fBestClusterIndex = -1;
  fClusters.reserve(10);
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
    { "goodblock.tdctime", "Energy of good blocks", "fGoodBlocks.TDCTime"},
    { "goodblock.row", "Row of good blocks", "fGoodBlocks.row"},
    { "goodblock.col", "Col of good blocks", "fGoodBlocks.col"},
    { "goodblock.x", "x pos (m) of good blocks", "fGoodBlocks.x"},
    { "goodblock.y", "y pos (m) of good blocks", "fGoodBlocks.y"},
    { "goodblock.id", "Element ID of good blocks", "fGoodBlocks.id"},
    { "goodblock.cid", "Cluster index of good blocks", "fGoodBlocks.cid"},
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
	fGoodBlocks.cid.push_back(-1); //initialize good block cluster id to -1
	fGoodBlocks.x.push_back(blk->GetX());
	fGoodBlocks.y.push_back(blk->GetY());
	//
	//
	if(fModeADC != SBSModeADC::kWaveform) {
	  const SBSData::PulseADCData &ahit = blk->ADC()->GetGoodHit();
	  blk->SetE(ahit.integral.val);
	  //blk->SetAintp(ahit.integral.val / blk->ADC()->GetGain());
	  blk->SetAgain(blk->ADC()->GetGain() / blk->ADC()->GetTrigCal());
	  blk->SetAtime(ahit.time.val);
	  fGoodBlocks.e.push_back(ahit.integral.val);
	  fGoodBlocks.ADCTime.push_back(ahit.time.val);
	} else {
	  SBSData::Waveform *wave = blk->Waveform();
	  blk->SetE(wave->GetIntegral().val);
	  //blk->SetAintp(wave->GetIntegral().val / wave->GetGain());
	  blk->SetAgain(wave->GetGain() / wave->GetTrigCal());
	  blk->SetAtime(wave->GetTime().val);
	  fGoodBlocks.e.push_back(wave->GetIntegral().val);
	  fGoodBlocks.ADCTime.push_back(wave->GetTime().val);
	}
	if (WithTDC() && blk->TDC()->HasData() ) { 
	  const SBSData::TDCHit &hit = blk->TDC()->GetGoodHit();
	  fGoodBlocks.TDCTime.push_back(hit.le.val);
	  blk->SetTDCtime(hit.le.val);
	} else {
	  fGoodBlocks.TDCTime.push_back(-1000.);
	  blk->SetTDCtime(-1000.);
	}
	//	 std::cout << blk->GetID() << " set tdc time = " << blk->GetTDCtime() << " set adc time = " << blk->GetAtime() << std::endl;
      }
    }
  }
  // Put good blocks in fBlockSet to use in FindCluster
  //  fBlockSet.reserve(fGoodBlocks.e.size());
  fBlockSet.clear();
  for (UInt_t nb=0;nb< fGoodBlocks.e.size();nb++) {
    SBSBlockSet c1 = {fGoodBlocks.e[nb],fGoodBlocks.x[nb],fGoodBlocks.y[nb],fGoodBlocks.row[nb],fGoodBlocks.col[nb],fGoodBlocks.id[nb],fGoodBlocks.TDCTime[nb],fGoodBlocks.ADCTime[nb],kFALSE};
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

  fBestClusterIndex = -1;

  double Emax = 0.0;

  std::map<int,int> clusterids_by_blockid; //key = block id; mapped value = cluster ID.
  
  //Test new experimental clustering algorithm here; island algorithm: AJRP
  while( NSize != 0 ){
    std::sort(fBlockSet.begin(), fBlockSet.end(), [](const SBSBlockSet& c1, const SBSBlockSet& c2) { return c1.e > c2.e;});

    fBlockSetIterator it = fBlockSet.begin(); //This starts out as the highest-energy remaining block 
    
    if ( (*it).e > fEmin_clusSeed ){ //don't bother if the highest-energy remaining block is too low in energy

      SBSElement *blk= fElements[(*it).id-fChanMapStart] ; //here blk is the pointer to the highest energy block remaining in fBlockSet, the cluster seed
      SBSCalorimeterCluster* cluster = new SBSCalorimeterCluster(fBlockSet.size(),blk); //seed a new cluster with blk but don't add it to the cluster array yet"

      fBlockSet.erase(it); //"erase" the cluster seed from fBlockSet. 
      NSize--; //decrement the total number of remaining unused blocks in fBlockSet
      
      Int_t iblk=0;

      // This loop implements the "island" algorithm: loop on
      // all existing blocks in the cluster, and keep adding neighbors
      // of each block in the cluster until we either hit the max
      // cluster size or we run out of cluster blocks to test:
      while( iblk<cluster->GetMult() && cluster->GetMult() < cluster->GetNMaxElements() ){
	SBSElement *blk_i = cluster->GetElement( iblk ); //grab pointer to the ith block in the cluster:
	//Now loop on all the remaining unused blocks.
	Double_t xblk = blk_i->GetX();
	Double_t yblk = blk_i->GetY();
	Double_t tblk = blk_i->GetAtime();
	Double_t eblk = blk_i->GetE();

	//now loop on all the unused blocks and check if they are within Rmax of the current block. It is safer, albeit slightly less efficient, to loop on fBlockSet once and then erase any blocks that were added to the cluster in a second loop. For the erasing, it is simpler to refer to the block id than to the index in fBlockSet.
	std::set<int> blockidstoerase;
	
	for( int j=0; j<fBlockSet.size(); j++ ){
	  SBSElement *blk_j = fElements[fBlockSet[j].id-fChanMapStart];
	  Double_t xblk_j = blk_j->GetX();
	  Double_t yblk_j = blk_j->GetY();
	  Double_t eblk_j = blk_j->GetE();
	  Double_t tblk_j = blk_j->GetAtime();

	  Double_t Rad2 = pow( xblk-xblk_j, 2 ) + pow( yblk-yblk_j, 2 );
	  //Even if we later update cluster::AddElement to set cluster time to energy-weighted mean time, this will still compare the current block to the seed: 
	  Double_t tdiff = tblk_j - cluster->GetElement( 0 )->GetAtime();

	  if( Rad2 < pow(fRmax_dis,2) && fabs( tdiff ) < fTmax &&
	      cluster->GetMult() < cluster->GetNMaxElements() ){
	    cluster->AddElement( blk_j );
	    blockidstoerase.insert( fBlockSet[j].id ); //list of blocks to remove from the set
	  }
	}

	//now loop on fBlockSet again using an iterator and erase any elements that were added to the cluster.
	auto it2 = fBlockSet.begin();
	while( it2 != fBlockSet.end() ){
	  if( blockidstoerase.find( (*it2).id ) != blockidstoerase.end() ){ //this block was added to the cluster; erase from fBlockSet
	    //auto itnext = fBlockSet.erase( it2 ); //return value is an iterator to the next element
	    //it2 = itnext; //Perhaps it is safe to assign the result of fBlockSet.erase( it2 ) to it2 itself, making this division on two lines unnecessary?
	    it2 = fBlockSet.erase( it2 );
	    NSize--;
	  } else {
	    ++it2;
	  }
	}

	//We don't really need to sort fBlockSet again until we start
	//the next cluster.
	
	iblk++; 
      } //end loop over blocks in current cluster. On each iteration of this loop, all blocks in fblockset that are added to the cluster are erased from fBlockSet, ensuring that any given block can only be added to a cluster exactly once and can only be used in exactly one cluster!

      //NOW we add the cluster to the array, IF the total energy is above
      // threshold:
      if( cluster->GetE() >= fEmin_clusTotal ){
	fClusters.push_back( cluster );

	clusterids_by_blockid.clear();
	
	for( iblk=0; iblk<cluster->GetMult(); iblk++ ){

	  int blkid_i = cluster->GetElement(iblk)->GetID();
	  clusterids_by_blockid[blkid_i] = fClusters.size()-1;
	}
	
	if( fClusters.size() == 1 || cluster->GetE() > Emax ){
	  Emax = cluster->GetE();
	  fBestClusterIndex = fClusters.size()-1;
	}
      } else delete cluster; //prevent memory leak
    } else break; //end check that primary block is above cluster seed threshold. If the primary block is below cluster seed threshold, exit the loop; there are no more clusters to find; AND we would get stuck in infinite loop without this break statement;
    
  } //End loop on while( NSize != 0 ). Each iteration of this loop finds one cluster, as long as there are more clusters to find.

  //NOW clustering is finished; we need to loop on fGoodBlocks and assign the cluster IDs;
  for( int iblk=0; iblk<fGoodBlocks.id.size(); iblk++ ){
    auto foundblock = clusterids_by_blockid.find( fGoodBlocks.id[iblk] );

    if( foundblock != clusterids_by_blockid.end() ){
      fGoodBlocks.cid[iblk] = foundblock->second;
    }
  }
  //We want the "best" cluster to be the one with the largest total energy (in general)

  //For now, keeping the old clustering code here, but commented out. 
  // while ( NSize != 0 )  {
  //   std::sort(fBlockSet.begin(), fBlockSet.end(), [](const SBSBlockSet& c1, const SBSBlockSet& c2) { return c1.e > c2.e;});
    
  //   fBlockSetIterator it = fBlockSet.begin(); //This starts out as the highest-energy remaining block 
    
  //   if ( (*it).e > fEmin_clusSeed ){ //don't bother if the highest-energy remaining block is too low in energy
  //     Bool_t AddingBlocksToCluster = kTRUE;
  //     SBSElement *blk= fElements[(*it).id-fChanMapStart] ; //here blk is the pointer to the highest energy block remaining in fBlockSet, the cluster seed
  //     SBSCalorimeterCluster* cluster = new SBSCalorimeterCluster(fBlockSet.size(),blk); //seed a new cluster with blk
  //     fClusters.push_back(cluster); //add it to the array.
  //     fBlockSet.erase(it); //"erase" the cluster seed from fBlockSet. 
  //     NSize--; //decrement the total number of remaining unused blocks in fBlockSet

  //     while (AddingBlocksToCluster) {
  // 	Bool_t Close=kFALSE; //Check if block within user-def radius and time to add to cluster
  // 	fBlockSetIterator it2 = fBlockSet.begin(); //it2 is an iterator to highest-energy remaining block in fBlockSet after removing the primary block (cluster seed)
  // 	SBSElement *blk_p= fElements[(*fBlockSet.begin()).id-fChanMapStart]; //blk_p is a pointer to the element of the highest-energy remaining block, as referenced by it2. THIS COULD BE ANYWHERE IN THE CALORIMETER!
  //       Double_t pTime=blk_p->GetAtime(); //pTime is the ADC time of blk_p 

  // 	while (!Close && (it2 < fBlockSet.end())) { 
  // 	  SBSElement *blk= fElements[(*it2).id-fChanMapStart]  ;  //this is shadowing the declaration above. Unintended consquences? On first iteration this is the same as blk_p. on subsequent iterations it points somewhere else:
  // 	  Int_t Index = fClusters.size()-1;
  // 	  Double_t Rad = sqrt( pow((fClusters[Index]->GetX()-blk->GetX()),2) + pow((fClusters[Index]->GetY()-blk->GetY()),2) ); //check distance of block from cluster energy-weighted average position
  // 	  Double_t tDiff = blk->GetAtime()-pTime; //HUGE MISTAKE! we are comparing blk time with pTime, which could be some random block from elsewhere!
  // 	  //Close =( Rad<fRmax_dis && fabs(tDiff)<fTmax );
  // 	  Close = ( Rad<fRmax_dis );
  // 	  if (Close) {
  // 	    fClusters[Index]->AddElement(blk); //we found the highest-energy remaining block close enough to the cluster mean position; add it to the cluster and exit the loop without incrementing it2 or erasing it2 from fBlockSet
  // 	  } else {	       
  // 	    ++it2; //increment it2 (keep looking)
  // 	  }
  // 	}
  // 	//When we get to this point, we have either found one block that is "Close" to it2 or we have reached the end of fBlockSet. 
  // 	if (it2 == fBlockSet.end()) AddingBlocksToCluster = kFALSE; //no more blocks close enough to add to the cluster. exit
  // 	if (Close)   { //we found a block to add to the cluster and it is referenced by it2. 
  // 	  fBlockSet.erase(it2); //erase it2 from fBlockSet; it was already added to the cluster
  // 	  NSize--; //decrement NSize
  // 	}	
  //     } //At this point we have added all blocks to the cluster that are sufficiently close to the cluster mean position, which is a moving target as we add new blocks. Need to implement island algorithm.
      
  //     // Adding total cluster energy threshold
  //     if ( (cluster->GetE())<fEmin_clusTotal ) fClusters.pop_back(); //erase the found cluster from fClusters. we still have the pointer. 
  //     //could the line above interact with the lines below in a problematic way? Almost certainly yes.
      
  //     //After we are done adding blocks to the cluster, we check if this is the cluster
  //     //with the largest total energy:
      
  //     if( fClusters.size() == 1 || cluster->GetE() > Emax ){
  // 	Emax = cluster->GetE();
  // 	fBestClusterIndex = fClusters.size()-1;
  //     }
  //   } else break;
  // }
  //

  //Fill main cluster variables here; this may get overridden by derived classes:
  if(!fClusters.empty()) {
    SBSCalorimeterCluster *clus = fClusters[fBestClusterIndex];
    fMainclus.e.push_back(clus->GetE());
    //fMainclus.e_c.push_back(clus->GetE()*(fConst + fSlope*fAccCharge));
    fMainclus.again.push_back(clus->GetAgain());
    fMainclus.atime.push_back(clus->GetAtime());
    fMainclus.tdctime.push_back(clus->GetTDCtime());
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
  out.x.clear();
  out.y.clear();
  out.row.clear();
  out.col.clear();
  out.n.clear();
  out.blk_e.clear();
  // out.blk_e_c.clear();
  out.id.clear();
}
