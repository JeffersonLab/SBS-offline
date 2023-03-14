// *-- Author :    Guido Maria Urciuoli   12 March 2001

//////////////////////////////////////////////////////////////////////////
//
// SBSRICH
//
// The RICH detector
// Written by Guido Maria Urciuoli, INFN
// Adapted for Hall A Analyzer by Ole Hansen, JLab
// Adapted to SBS by Seamus Riordan, ANL and Eric Fuchey, UConn
//
//////////////////////////////////////////////////////////////////////////

#include "TMath.h"
#include "SBSGRINCH.h"
#include "THaTrack.h"
#include "THaEvData.h"
#include "THaGlobals.h"
#include "THaDetMap.h"
#include "THaSpectrometer.h"
#include "TError.h"
#include "VarDef.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include "THaBenchmark.h"

using namespace std;
// cout<<"stay and stop here!!!"<<endl;
//_____________________________________________________________________________
SBSGRINCH::SBSGRINCH( const char* name, const char* description, 
		  THaApparatus* apparatus ) :
  SBSCherenkovDetector(name,description,apparatus)
  //fTrackX(kBig), fTrackY(kBig)
{
  //keep this line first
  fBench = new THaBenchmark;

  // Normal constructor with name and description
  
  // default separation is 3.1 cm along X (rows) and Y (columns)
  // Y positions are staggered by 1.55 cm in alternating rows but they are not 
  // close-packed as if in a honeycomb layout.
  fMaxSep = 0.035; //3.5 cm is the max. center-to-center separation between "neighbor" PMTs
  // we can adjust via the DB if we want, but we'll use this as a default 

  fMaxSep2 = pow(fMaxSep,2);

  fTrackMatchXslope = 0.71; 
  fTrackMatchX0 = 0.0;
  fTrackMatchXsigma = 0.014; 
  fTrackMatchYslope = 1.8;
  fTrackMatchY0 = -0.06;
  fTrackMatchYsigma = 0.04; 
  fTrackMatchNsigmaCut = 4.5; //default to a relatively wide cut:

  Clear();
}

//_____________________________________________________________________________
SBSGRINCH::~SBSGRINCH()
{
   // Destructor. Remove variables from global list and free up the memory
  // allocated by us.
  Clear();// so the prgram doesn't complain when deleting clusters
  RemoveVariables();
  // delete fHits;
  // delete fResolvedHits;
  // delete fClusters;
  // delete fResolvedClusters;
  // delete [] fMIPs;
  // delete [] fXseg;
  // delete fBench;
}

//_____________________________________________________________________________
void SBSGRINCH::Clear( Option_t* opt )
{
  // Reset event-by-event data
  if(fDebug)cout << "SBSGRINCH::Clear() " << endl; 
  if( fDoBench ) fBench->Begin("Clear");
  SBSCherenkovDetector::Clear(opt);
  if( fDoBench ) fBench->Stop("Clear");
}

//_____________________________________________________________________________
Int_t SBSGRINCH::ReadDatabase( const TDatime& date )
{
  // 
  // Read the database for this detector.
  // This function is called once at the beginning of the analysis.
  // 'date' contains the date/time of the run being analyzed.
  //
  static const char* const here = "ReadDatabase";
  cout << here << endl;
  // Open the database file
  FILE* fi = OpenFile( date );
  if( !fi ) return kFileError;

  Int_t err = SBSCherenkovDetector::ReadDatabase(date);
  if(err) {
    return err;
  }
  fIsInit = false;
  
  //add grinch specific stuff here if needed

  const DBRequest request[] = { 
    { "maxsep", &fMaxSep, kDouble, 0, 1, 1 },
    { "trackmatch_xslope", &fTrackMatchXslope, kDouble, 0, 1, 1 },
    { "trackmatch_x0", &fTrackMatchX0, kDouble, 0, 1, 1 },
    { "trackmatch_xsigma", &fTrackMatchXsigma, kDouble, 0, 1, 1 },
    { "trackmatch_yslope", &fTrackMatchYslope, kDouble, 0, 1, 1 },
    { "trackmatch_y0", &fTrackMatchY0, kDouble, 0, 1, 1 },
    { "trackmatch_ysigma", &fTrackMatchYsigma, kDouble, 0, 1, 1 },
    { "trackmatch_nsigma_cut", &fTrackMatchNsigmaCut, kDouble, 0, 1, 1 },
    {0}
  };

  err = LoadDB( fi, date, request, fPrefix );

  fMaxSep2 = pow(fMaxSep,2);

  if( err != 0 ){
    fclose(fi);
    return err; 
  }

  fIsInit = true;
  
  fclose(fi);
  return kOK;
}

//_____________________________________________________________________________
Int_t SBSGRINCH::DefineVariables( EMode mode )
{
  // Define (or delete) global variables of the detector
  
  if( mode == kDefine && fIsSetup ) return kOK;
  //fIsSetup = ( mode == kDefine );

  Int_t err = SBSCherenkovDetector::DefineVariables(mode);
  if(err)
    return err;

  
  //Add here the cluster output variables? maybe
  RVarDef var2[] = {
    { "nclus",        " number of GRINCH PMT clusters",  "GetNumClusters()"                                  },
    { "clus.size",    " GRINCH cluster size",            "fClusters.SBSCherenkov_Cluster.GetNHits()"            },
    { "clus.x_mean",  " GRINCH cluster X center",        "fClusters.SBSCherenkov_Cluster.GetXcenter()"          },
    { "clus.y_mean",  " GRINCH cluster Y center",        "fClusters.SBSCherenkov_Cluster.GetYcenter()"          },
    { "clus.t_mean", " GRINCH cluster mean time",       "fClusters.SBSCherenkov_Cluster.GetMeanTime()"         },
    { "clus.tot_mean", " GRINCH cluster mean time over threshold", "fClusters.SBSCherenkov_Cluster.GetMeanToT()" },
    { "clus.adc",     " GRINCH cluster total charge",    "fClusters.SBSCherenkov_Cluster.GetCharge()"           },
    { "clus.trackindex", "Index of matched track",      "fClusters.SBSCherenkov_Cluster.GetTrackIndex()" },
    { "bestcluster", "Index of best cluster", "fBestClusterIndex" },
    { 0 }
  };
  DefineVarsFromList( var2, mode, "" );// (re)define path here...
  
  return kOK;
}


//_____________________________________________________________________________
Int_t SBSGRINCH::Decode( const THaEvData& evdata )
{
  //Decode RICH data and fill hit array
  if(fDebug){
    cout << "SBSGRINCH::Decode " << endl;
    cout << " Is Init ? " << fIsInit << " Is Physics Trigger ? " << evdata.IsPhysicsTrigger() << endl;
  }
  return SBSCherenkovDetector::Decode(evdata);
}

//_____________________________________________________________________________
Int_t SBSGRINCH::CoarseProcess( TClonesArray& tracks )
{
  if(fDebug)cout << "Begin Coarse Process" << endl;
  if( fDoBench ) fBench->Begin("CoarseProcess");
 
  SBSCherenkovDetector::CoarseProcess(tracks);
  
  // add clustering here: hits have been filtered by SBSCherenkovDetector::CoarseProcess()
  // clustering ought to be coded in function "FindClusters" below
  
  /*int err = */FindClusters();
  
  if( fDoBench ) fBench->Stop("CoarseProcess");
  if(fDebug)cout << "End Coarse Process" << endl;
  
  return 0;
}


//_____________________________________________________________________________
Int_t SBSGRINCH::FineProcess( TClonesArray& tracks )
{
  if(fDebug)cout << "Begin Fine Process" << endl;
  // The main GRINCH processing method. 
  // Here we attempt to match a track to each cluster

  if( !fIsInit ) return -255;

  if( fDoBench ) fBench->Begin("FineProcess");

  SBSCherenkovDetector::FineProcess(tracks);
  
  Int_t nmatch = 0;
  // Clusters matched with tracks here (only if there are any tracks to match)
  if(tracks.GetLast()>=0){
    nmatch = MatchClustersWithTracks(tracks);
  }
  fBestClusterIndex = SelectBestCluster(nmatch);
  
  if( fDoBench ) fBench->Stop("FineProcess");
  if(fDebug)cout << "Done fine process " << endl;
  
  return 0;
}


//__________________________________________________________________________
Int_t SBSGRINCH::FindClusters()
{
  // here will be implemented the clustering for the GRINCH
  // Basic idea: loop on all the good hits; group all contiguous nearest-neighbor hits into clusters based on position and/or row/col number; 
  // when this is called, the "fHits" array has already been filled by CoarseProcess
  
  Int_t ngoodhits = GetNumHits(); 

  Int_t nclust = 0;

  std::set<Int_t> list_pmts; //list of all fired PMTs
  std::map<Int_t,Int_t> hitindex_pmts; //index in hit array of fired PMTs
  //std::map<Int_t,Bool_t> hitused_pmts; //flag to indicate a PMT has already been used in a cluster

  std::map<Int_t,Int_t> clustindex_pmts; //cluster index of pmts (maybe redundant, unclear if necessary yet)

  std::vector<std::vector<Int_t> > pmtlist_clusters;

  for( int ihit=0; ihit<ngoodhits; ihit++ ){
    int pmtnum = GetHit(ihit)->GetPMTNum();
    list_pmts.insert( pmtnum );
    hitindex_pmts[pmtnum] = ihit;
    //hitused_pmts[pmtnum] = false;
    clustindex_pmts[pmtnum] = -1;
  }
  
  //ALL PMTs will be part of at least one "cluster". Later we will also consider timing here:
  if( ngoodhits > 0 ){
    for( auto ipmt : list_pmts ){
      nclust = pmtlist_clusters.size();
      //if( !hitused_pmts[ipmt] ){ //find a cluster to put this PMT into!
      // Check existing clusters to see if this is a nearest-neighbor of 
      // any PMT in an existing cluster:
      Double_t xpmt = GetHit(hitindex_pmts[ipmt])->GetX();
      Double_t ypmt = GetHit(hitindex_pmts[ipmt])->GetY();
      
      for( int iclust=0; iclust<nclust; iclust++ ){ //loop over all existing clusters
	bool isneighbor = false;
	for( int ihit=0; ihit<pmtlist_clusters[iclust].size(); ihit++ ){ //loop over pmts in the cluster:
	  //Check if this PMT is within maxsep of any hit in an existing cluster:
	  Double_t xtest = GetHit(hitindex_pmts[pmtlist_clusters[iclust][ihit]])->GetX();
	  Double_t ytest = GetHit(hitindex_pmts[pmtlist_clusters[iclust][ihit]])->GetY();
	  
	  Double_t sep2 = pow( xpmt - xtest, 2 ) + pow( ypmt - ytest, 2 );
	  if( sep2 <= fMaxSep2  ){ // this is a neighbor to at least one of the hits in an existing cluster! --> Add it!
	    clustindex_pmts[ipmt] = iclust; 
	    //hitused_pmts[ipmt] = true;
	    isneighbor = true;
	    break; //exit the loop, we don't need to consider any more pmts
	  }
	}
	if( isneighbor ){
	  pmtlist_clusters[iclust].push_back( ipmt );
	}
      }

      if( clustindex_pmts[ipmt] < 0 ){ //PMT is not a neighbor to any existing cluster; create a new one!
	std::vector<Int_t> pmtlist_temp(1,ipmt);
	pmtlist_clusters.push_back( pmtlist_temp );
      }
    }

    nclust = pmtlist_clusters.size();

    for( int iclust=0; iclust<nclust; iclust++ ){
      new( (*fClusters)[iclust] ) SBSCherenkov_Cluster();
      SBSCherenkov_Cluster *clusttemp = ( (SBSCherenkov_Cluster *) (*fClusters)[iclust] );
      
      for( int ihit=0; ihit<pmtlist_clusters[iclust].size(); ihit++ ){
	//insert the hit into the cluster. This automatically updates the cluster properties
	SBSCherenkov_Hit *hittemp = GetHit( hitindex_pmts[pmtlist_clusters[iclust][ihit]] );

	clusttemp->Insert( hittemp );

	hittemp->SetClustIndex( iclust );
      }
    }
  }

  return 0;
}

//_____________________________________________________________________________
Int_t SBSGRINCH::MatchClustersWithTracks( TClonesArray& tracks )
{
  Int_t ntracks = tracks.GetLast() + 1; 
  Int_t nclust = fClusters->GetLast()+1;

  Int_t nmatch = 0;

  for( int itrack=0; itrack<ntracks; itrack++ ){
    auto* theTrack = static_cast<THaTrack*>( tracks.At(itrack) ); 
   
    // Now we need to come up with some sort of matching criterion for GRINCH clusters and tracks. 
    // I suppose what we want for now is the cluster closest to the track along X and/or Y This could become more 
    // sophisticated over time
 
    //Project track to GRINCH entrance window: 
    double zGRINCH = GetOrigin().Z();

    double xtrack = theTrack->GetX(zGRINCH);
    double ytrack = theTrack->GetY(zGRINCH);
    
    double minxdiff = 0.0; 
    double mindiff2 = 0.0;
    int iclmin=-1;

    for( int iclust=0; iclust<nclust; iclust++ ){
      //SBSCherenkov_Cluster *clusttemp = ( (SBSCherenkov_Cluster*) (*fClusters)[iclust] );
      SBSCherenkov_Cluster *clusttemp = GetCluster( iclust );

      
      double xGRINCH = clusttemp->GetXcenter();
      double yGRINCH = clusttemp->GetYcenter();

      double xdiff =  fabs(xtrack - fTrackMatchXslope * xGRINCH - fTrackMatchX0);
      double ydiff = fabs(ytrack - fTrackMatchYslope * yGRINCH - fTrackMatchY0);

      double diff2 = pow( xdiff/fTrackMatchXsigma, 2 ) + pow( ydiff/fTrackMatchYsigma, 2 );
    
      //     if( xdiff <= fTrackMatchXcut ){
      if( diff2 <= pow( fTrackMatchNsigmaCut, 2 ) ){
	//if( iclmin < 0 || xdiff < minxdiff ){
	if( iclmin < 0 || diff2 < mindiff2 ){
	  iclmin = iclust;
	  //  minxdiff = xdiff;
	  mindiff2 = diff2;
	}
      }
    }

    if( iclmin >= 0 ){ //match:
      //SBSCherenkov_Cluster *clusttemp = ( (SBSCherenkov_Cluster*) (*fClusters)[iclmin] );
      SBSCherenkov_Cluster *clusttemp = GetCluster( iclmin );

      clusttemp->SetTrackIndex( itrack ); 
      clusttemp->SetTrack( theTrack );
      //We also need to loop on all the hits and set the track index for those hits:
      TIter next(clusttemp->GetHitList());
      SBSCherenkov_Hit *hittemp;
      while( ( hittemp = static_cast<SBSCherenkov_Hit*>( next() ) ) ){
	hittemp->SetTrackIndex( itrack );
      }

      nmatch++;
    }
    
  }

  // here will be implemented the matching of GRINCH clusters with tracks
  return nmatch;
}

Int_t SBSGRINCH::SelectBestCluster(Int_t nmatch){
  Int_t nclusters = fClusters->GetLast() + 1; 
  
  Int_t best=-1;
  Int_t maxsize=0;
  for( int iclust=0; iclust<nclusters; iclust++ ){
    SBSCherenkov_Cluster *clusttemp = GetCluster(iclust);

    if( clusttemp->GetTrackIndex() >= 0 || nmatch <= 0 ){ //either this cluster is matched with a track or there are no clusters matched with tracks and we just pick the biggest one:
      if( best < 0 || clusttemp->GetNHits() > maxsize ){
	best = iclust;
	maxsize = clusttemp->GetNHits();
      }
    }
  }
  return best;
}

ClassImp(SBSGRINCH)

