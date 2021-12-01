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
  fEmin = 1.0; // 1 MeV minimum energy to be in cluster
  fEmin_seed = 1.0; // 1 MeV minimum energy to seed a cluster
  fTmax = 1000.0; // 1000 ns maximum arrival time difference with seed to be in cluster
  fXmax_dis = .30; // Maximum X (m) distance from cluster center to be included in cluster
  fYmax_dis = .30; // Maximum Y (m) distance from cluster center to be included in cluster
  fRmax_dis = .30; // Maximum Radius (m) from cluster center to be included in cluster
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
    { "emin",         &fEmin,   kDouble, 0, false }, ///< minimum energy threshold
    { "emin_seed",    &fEmin_seed,   kDouble, 0, false }, ///< minimum energy threshold for seed
    { "tmax",         &fTmax,   kDouble, 0, false }, ///< maximum time difference for block
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
    { "e_c",    "Corrected Energy (MeV) of largest cluster",    "GetECorrected()" },
    { "atimeblk", "ADC time of highest energy block in the largest cluster", "GetAtime()" },
    { "tdctimeblk", "TDC time of highest energy block in the largest cluster", "GetTDCtime()" },
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
    {0}
  };
  err = DefineVarsFromList( vars_gb, mode );
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
  return err;
}

//_____________________________________________________________________________
void SBSCalorimeter::ClearEvent()
{
  SBSGenericDetector::ClearEvent();
  ClearOutputVariables();
  DeleteContainer(fClusters);
  fGoodBlocks.clear();
  fBlockSet.clear();
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
	 //
	 //
	 if(fModeADC != SBSModeADC::kWaveform) {
	   const SBSData::PulseADCData &ahit = blk->ADC()->GetGoodHit();
	   blk->SetE(ahit.integral.val);
	   blk->SetAtime(ahit.time.val);
	   fGoodBlocks.e.push_back(ahit.integral.val);
            fGoodBlocks.ADCTime.push_back(ahit.time.val);
	 } else {
           SBSData::Waveform *wave = blk->Waveform();
	   blk->SetE(wave->GetIntegral().val);
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
 	//
	Int_t NSize = fBlockSet.size();
        while ( NSize != 0 )  {
             std::sort(fBlockSet.begin(), fBlockSet.end(), [](const SBSBlockSet& c1, const SBSBlockSet& c2) { return c1.e > c2.e;});

	     //SAS - Add minimum seed energy requirement here, after the blocks are sorted. Return first element of fBlockset and check it against fEmin_seed

	     const SBSBlockSet& cS = fBlockSet.begin(); 
	     Bool_t AddingBlocksToCluster = kTRUE;
	     SBSElement *blk= fElements[(*it).id-fChanMapStart] ; 
	     SBSCalorimeterCluster* cluster = new SBSCalorimeterCluster(fBlockSet.size(),blk);

	     if( cS.e < fEmin_seed ){
	       AddingBlocksToCluster = kFALSE;
	       fBlockSet.erase(it);
	       NSize--;
	       continue;
	     }	     

             fClusters.push_back(cluster);
	     fBlockSet.erase(it);
	     NSize--;
	while (AddingBlocksToCluster) {
	  Bool_t InTime=kFALSE;
	  Bool_t IsNeighbor=kFALSE;
	     fBlockSetIterator it2 = fBlockSet.begin();
	    while (!IsNeighbor && (it2 < fBlockSet.end())) {
	      SBSElement *blk= fElements[(*it2).id-fChanMapStart]  ; 
	      SBSElement *blk_p= fElements[fBlockSet.begin()-fChanMapStart];
	      Int_t Index = fClusters.size()-1;
	      Double_t Rad = sqrt( pow((fClusters[Index]->GetX()-blk->GetX()),2) + pow((fClusters[Index]->GetY()-blk->GetY()),2) );
 	      IsNeighbor =( Rad<fRmax_dis );
	      InTime=( fabs( blk->GetAtime()-blk_p->GetAtime() ) < fTmax);
     	      if (IsNeighbor&&InTime) {
		fClusters[Index]->AddElement(blk);
	      } else {	       
		  ++it2;
              }
	    }
	    if (it2 == fBlockSet.end()) AddingBlocksToCluster = kFALSE;
	    if (IsNeighbor)   {
               fBlockSet.erase(it2);
               NSize--;
	    }
	}
	}
	//
  if(!fClusters.empty()) {
    SBSCalorimeterCluster *clus = fClusters[0];
    fMainclus.e.push_back(clus->GetE());
    fMainclus.e_c.push_back(clus->GetE()*(fConst + fSlope*fAccCharge));
    fMainclus.atime.push_back(clus->GetAtime());
    fMainclus.tdctime.push_back(clus->GetTDCtime());
    fMainclus.x.push_back(clus->GetX());
    fMainclus.y.push_back(clus->GetY());
    fMainclus.n.push_back(clus->GetMult());
    fMainclus.blk_e.push_back(clus->GetEblk());
    fMainclus.blk_e_c.push_back(clus->GetEblk()*(fConst + fSlope*fAccCharge));
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
  if(!fClusters.empty()) {
    SBSCalorimeterCluster *clus = fClusters[0];
 
    if(fDataOutputLevel > 0 ) {
      for(Int_t nc=0;nc<clus->GetMult();nc++ ) {
	SBSElement *blk= clus->GetElement(nc);
        fMainclusblk.e.push_back(blk->GetE());
        fMainclusblk.e_c.push_back(blk->GetE()*(fConst + fSlope*fAccCharge));
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
    fOutclus.e_c.reserve(nres);
    fOutclus.atime.reserve(nres);
    fOutclus.tdctime.reserve(nres);
    fOutclus.x.reserve(nres);
    fOutclus.y.reserve(nres);
    fOutclus.n.reserve(nres);
    fOutclus.blk_e.reserve(nres);
    fOutclus.blk_e_c.reserve(nres);
    fOutclus.row.reserve(nres);
    fOutclus.col.reserve(nres);
    fOutclus.id.reserve(nres);

    int nclus = 0;
    for( const auto* cluster: fClusters ) {
      if(nclus < fMaxNclus) { // Keep adding them until we reach fMaxNclus
        fOutclus.e.push_back(cluster->GetE());
        fOutclus.e_c.push_back(cluster->GetE()*(fConst + fSlope*fAccCharge));
        fOutclus.atime.push_back(cluster->GetAtime());
        fOutclus.tdctime.push_back(cluster->GetTDCtime());
        fOutclus.x.push_back(cluster->GetX());
        fOutclus.y.push_back(cluster->GetY());
        fOutclus.n.push_back(cluster->GetMult());
        fOutclus.blk_e.push_back(cluster->GetEblk());
        fOutclus.blk_e_c.push_back(cluster->GetEblk()*(fConst + fSlope*fAccCharge));
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
  fNclus = 0;
}



void SBSCalorimeter::ClearCaloOutput(SBSCalorimeterOutput &out)
{
  out.e.clear();
  out.e_c.clear();
  out.atime.clear();
  out.tdctime.clear();
  out.x.clear();
  out.y.clear();
  out.row.clear();
  out.col.clear();
  out.n.clear();
  out.blk_e.clear();
  out.blk_e_c.clear();
  out.id.clear();
}
