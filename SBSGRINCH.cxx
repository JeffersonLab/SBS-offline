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
  fMaxSep = 0.07; //3.5 cm is the max. center-to-center separation between "neighbor" PMTs, this is twice that number
  // we can adjust via the DB if we want, but we'll use this as a default 

  fMaxSep2 = pow(fMaxSep,2);

  fMaxTdiffClust = 30.0; //ns (default);
  fMaxTdiffRefTime = 30.0; //ns (default)

  fGRINCH_dToffset = 0.0; //ns (default);
  fGRINCH_dTsigma = 2.0; //ns (default)

  fUseRefTimeDiffCut = false; //default to false
  
  //Set up default track match cuts based on mirror number:

  fNmirror = 4; 

  fTrackMatchPslope = 0.1715;

  fOrderTrackMatchY = 3; //default values.
  
  fTrackMatchXslope.resize(fNmirror);
  fTrackMatchX0.resize(fNmirror);
  fTrackMatchXsigma.resize(fNmirror);
  fTrackMatchXmin.resize(fNmirror);
  fTrackMatchXmax.resize(fNmirror);

  fTrackMatchYslope.resize(fNmirror*fOrderTrackMatchY);
  fTrackMatchY0.resize(fNmirror);
  fTrackMatchYsigma.resize(fNmirror);

   
  //  Double_t mirrx0[] = { 0.1789, 0.0052, 0.0052, -0.2576 };

  //These parameters come from Maria's analysis:
  Double_t mirrx0[] = {0.17, 0.005, 0.0095, -0.27}; //intercept of dx versus theta
  Double_t mirrxslope[] = {1.361, 1.424, 1.3387, 1.466}; //slope of dx versus theta
  Double_t mirrxmin[] = {-0.8, -0.55, -0.05, 0.45}; //mirror boundary cuts, min
  Double_t mirrxmax[] = {-0.45, 0.05, 0.58, 0.85}; //mirror boundary cuts, max
  Double_t mirrxsigma[] = {0.013, 0.012, 0.013, 0.011}; //sigma of corrected dx distribution
  
  // Double_t mirrx0[] = { -0.1580, -0.0193, -0.0047, 0.1968};
  // Double_t mirrxslope[] = {0.71, 0.71, 0.71, 0.71};
  // Double_t mirrxmin[] = {-0.8, -0.55, -0.1, 0.45};
  // Double_t mirrxmax[] = {-0.45, 0.1, 0.6, 0.85};
  // Double_t mirrxsigma[] = {0.016, 0.014, 0.015, 0.016};

  // Double_t mirryslope[] = {0.42, 0.48, 0.173, 0.56};
  // Double_t mirry0[] = {-0.017, 0.018, 0.014, 0.046};
  // Double_t mirrysigma[] = {0.022, 0.015, 0.022, 0.023}; 
  //Let's implement Maria's 3rd-order polynomial in track phi for dy:
  //Need four parameters per mirror: 
  Double_t mirrypar[] = {  -0.002, -1.5, 10.0, -150.0,
			   0.041, -1.4, -5., -185.0,
			   0.051, -1.9, 0.0, -120.0,
			   0.081, -1.8, -25.0, 10.0 };
  Double_t mirrysigma[] = { 0.023, 0.020, 0.025, 0.020 };
  

  for( int imirr=0; imirr<fNmirror; imirr++ ){
    fTrackMatchXslope[imirr] = mirrxslope[imirr];
    fTrackMatchX0[imirr] = mirrx0[imirr];
    fTrackMatchXsigma[imirr] = mirrxsigma[imirr];
    fTrackMatchXmin[imirr] = mirrxmin[imirr];
    fTrackMatchXmax[imirr] = mirrxmax[imirr];

    //New y params; implement Maria's third-order polynomial fit versus track phi:
    fTrackMatchY0[imirr] = mirrypar[4*imirr];
    fTrackMatchYsigma[imirr] = mirrysigma[imirr];
    for( int ipar=1; ipar<=3; ipar++ ){
      fTrackMatchYslope[3*imirr + ipar-1] = mirrypar[4*imirr+ipar];
    }
  }

  fTrackMatchNsigmaCut = 4.5;

  // fTrackMatchXslope = 0.71; 
  // fTrackMatchX0 = 0.0;
  // fTrackMatchXsigma = 0.014; 
  // fTrackMatchYslope = 1.8;
  // fTrackMatchY0 = -0.06;
  // fTrackMatchYsigma = 0.04; 
  // fTrackMatchNsigmaCut = 4.5; //default to a relatively wide cut:

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

  fRefTimeIsSet = false;
  fRefTime = 0.0;
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

  int usereftimediffcut = fUseRefTimeDiffCut ? 1 : 0;
  
  const DBRequest request[] = { 
    { "maxsep", &fMaxSep, kDouble, 0, 1, 1 },
    { "maxtdiff_clust", &fMaxTdiffClust, kDouble, 0, 1, 1 },
    { "maxtdiff_ref", &fMaxTdiffRefTime, kDouble, 0, 1, 1 },
    { "usereftimediffcut", &usereftimediffcut, kInt, 0, 1, 1 },
    { "nmirror", &fNmirror, kInt, 0, 1, 1 },
    { "trackmatch_yorder", &fOrderTrackMatchY, kInt, 0, 1, 1 },
    { "trackmatch_pslope", &fTrackMatchPslope, kDouble, 0, 1, 1 },
    { "trackmatch_xslope", &fTrackMatchXslope, kDoubleV, 0, 1, 1 },
    { "trackmatch_x0", &fTrackMatchX0, kDoubleV, 0, 1, 1 },
    { "trackmatch_xsigma", &fTrackMatchXsigma, kDoubleV, 0, 1, 1 },
    { "trackmatch_yslope", &fTrackMatchYslope, kDoubleV, 0, 1, 1 },
    { "trackmatch_y0", &fTrackMatchY0, kDoubleV, 0, 1, 1 },
    { "trackmatch_ysigma", &fTrackMatchYsigma, kDoubleV, 0, 1, 1 },
    { "trackmatch_xmin", &fTrackMatchXmin, kDoubleV, 0, 1, 1 },
    { "trackmatch_xmax", &fTrackMatchXmax, kDoubleV, 0, 1, 1 },
    { "trackmatch_nsigma_cut", &fTrackMatchNsigmaCut, kDouble, 0, 1, 1 },
    { "trackmatch_dtsigma", &fGRINCH_dTsigma, kDouble, 0, 1, 1 },
    { "trackmatch_dtoffset", &fGRINCH_dToffset, kDouble, 0, 1, 1 },
    {0}
  };

  err = LoadDB( fi, date, request, fPrefix );

  fMaxSep2 = pow(fMaxSep,2);

  fUseRefTimeDiffCut = (usereftimediffcut != 0);
  
  //Check size of track match cut arrays, if they were loaded from the DB they had better all be the same size:
  bool sizecheck = ( fTrackMatchXslope.size() == fNmirror && 
		     fTrackMatchX0.size() == fNmirror && 
		     fTrackMatchXsigma.size() == fNmirror && 
		     fTrackMatchYslope.size() == fNmirror*fOrderTrackMatchY && 
		     fTrackMatchY0.size() == fNmirror && 
		     fTrackMatchYsigma.size() == fNmirror && 
		     fTrackMatchXmin.size() == fNmirror && 
		     fTrackMatchXmax.size() == fNmirror );

  if( !sizecheck ){
    Error(Here("SBSGRINCH::ReadDatabase"), TString::Format("Size of track match cut parameters does not match number of GRINCH mirrors (%d): Fix Database!",fNmirror) );
    return kInitError;
  }
		     

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
    { "allclus.size",    " GRINCH cluster size",            "fClusters.SBSCherenkov_Cluster.GetNHits()"            },
    { "allclus.trackindex", "Index of matched track",      "fClusters.SBSCherenkov_Cluster.GetTrackIndex()" },
    { "allclus.mirrorindex", "Mirror index of matched track", "fClusters.SBSCherenkov_Cluster.GetMirrorIndex()" },
    { "allclus.x_mean",  " GRINCH cluster X center",        "fClusters.SBSCherenkov_Cluster.GetXcenter()"          },
    { "allclus.y_mean",  " GRINCH cluster Y center",        "fClusters.SBSCherenkov_Cluster.GetYcenter()"          },
    { "allclus.t_mean", " GRINCH cluster mean time",       "fClusters.SBSCherenkov_Cluster.GetMeanTime()"         },
    { "allclus.t_mean_corr", " GRINCH cluster mean corrected time", "fClusters.SBSCherenkov_Cluster.GetMeanTimeCorrected()" },
    { "allclus.t_rms", "RMS of GRINCH cluster hit times", "fClusters.SBSCherenkov_Cluster.GetTimeRMS()" },
    { "allclus.t_rms_corr", "RMS of GRINCH cluster corrected times", "fClusters.SBSCherenkov_Cluster.GetTimeRMSCorrected()" },
    { "allclus.tot_mean", " GRINCH cluster mean time over threshold", "fClusters.SBSCherenkov_Cluster.GetMeanAmp()" },
    { "allclus.adc",     " GRINCH cluster total charge",    "fClusters.SBSCherenkov_Cluster.GetCharge()"           },
    { "allclus.dx", " GRINCH cluster corrected dx wrt track projection", "fClusters.SBSCherenkov_Cluster.GetDX()" },
    { "allclus.dy", " GRINCH cluster corrected dy wrt track projection", "fClusters.SBSCherenkov_Cluster.GetDY()" },
    { "bestcluster", "Index of best cluster", "fBestClusterIndex" },
    { 0 }
  };

  //Define global variables for "best cluster" properties:
  RVarDef var3[] = {
    { "clus.size",    " GRINCH best cluster size",            "fBestCluster.GetNHits()"            },
    { "clus.trackindex", "Index of best cluster matched track",      "fBestCluster.GetTrackIndex()" },
    { "clus.mirrorindex", "Mirror index of best cluster", "fBestCluster.GetMirrorIndex()" },
    { "clus.x_mean",  " GRINCH best cluster X center",        "fBestCluster.GetXcenter()"          },
    { "clus.y_mean",  " GRINCH best cluster Y center",        "fBestCluster.GetYcenter()"          },
    { "clus.t_mean", " GRINCH best cluster mean time",       "fBestCluster.GetMeanTime()"         },
    { "clus.t_mean_corr", " GRINCH best cluster mean corrected time", "fBestCluster.GetMeanTimeCorrected()" },
    { "clus.t_rms", "RMS of GRINCH best cluster hit times", "fBestCluster.GetTimeRMS()" },
    { "clus.t_rms_corr", " RMS of GRINCH best cluster corrected hit times", "fBestCluster.GetTimeRMSCorrected()" },
    { "clus.tot_mean", " GRINCH best cluster mean time over threshold", "fBestCluster.GetMeanAmp()" },
    { "clus.adc",     " GRINCH best cluster total charge",    "fBestCluster.GetCharge()"           },
    { "clus.dx", " GRINCH cluster corrected dx wrt track projection", "fBestCluster.GetDX()" },
    { "clus.dy", " GRINCH cluster corrected dy wrt track projection", "fBestCluster.GetDY()" },
    { 0 }
  };

  DefineVarsFromList( var2, mode, "" );// (re)define path here...
  DefineVarsFromList( var3, mode, "" );
  
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

  // NOTE: is there any particular reason we can't move the GRINCH clustering to FineProcess() to take advantage of
  // "reference" timing?
  // I don't see why not! 

  //GRINCH clustering moved to FineProcess so we can grab "reference" time from the
  //shower
  /*int err = *///FindClusters();
  
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

  SBSCherenkovDetector::FineProcess(tracks); //currently this does nothing

  //Moved clustering here so we can implement tighter time cuts internally to clustering:
  FindClusters();
  
  //Int_t nmatch = 0;
  fNtrackMatch = 0;
  // Clusters matched with tracks here (only if there are any tracks to match)
  if(tracks.GetLast()>=0){
    fNtrackMatch = MatchClustersWithTracks(tracks);
  }

  fBestClusterIndex = SelectBestCluster(fNtrackMatch);

  if( fBestClusterIndex >= 0 && fBestClusterIndex < fClusters->GetLast()+1 && tracks.GetLast() >= 0 ){
    fBestCluster = *(GetBestCluster());

    // std::cout << "ntrackmatch, best cluster index, best cluster track index = " 
    // 	      << fNtrackMatch << ", " << fBestClusterIndex << ", " << fBestCluster.GetTrackIndex()
    // 	      << std::endl;

  }
  
  if( fDoBench ) fBench->Stop("FineProcess");
  if(fDebug)cout << "Done fine process " << endl;
  
  return 0;
}


//__________________________________________________________________________
Int_t SBSGRINCH::FindClusters()
{
  // here will be implemented the clustering for the GRINCH
  // Basic idea: loop on all the good hits; group all contiguous nearest-neighbor hits into clusters based on position and/or row/col number (and time); 
  // when this is called, the "fHits" array has already been filled by CoarseProcess
  
  Int_t ngoodhits = GetNumHits(); 

  Int_t nclust = 0;

  //std::set<Int_t> list_pmts; //list of all fired PMTs (currently not actually used)
  std::set<Int_t> unused_pmts; //list of all PMTs that aren't part of clusters
  std::map<Int_t,Int_t> hitindex_pmts; //index in hit array of fired PMTs
  //std::map<Int_t,Bool_t> hitused_pmts; //flag to indicate a PMT has already been used in a cluster
  std::map<Int_t,Int_t> pmt_row;
  std::map<Int_t,Int_t> pmt_col; 
 
  std::map<Int_t,Double_t> xpmt;
  std::map<Int_t,Double_t> ypmt; 
  std::map<Int_t,Double_t> tpmt; 
  
  std::map<Int_t,Int_t> clustindex_pmts; //cluster index of pmts (maybe redundant, unclear if necessary yet)

  std::vector<std::set<Int_t> > pmtlist_clusters;
  std::vector<Double_t> tmean_clusters;
  
  pmtlist_clusters.clear();

  // If reference time has been set, start the first cluster
  // using the hit closest in time to the reference time. 
  
  for( int ihit=0; ihit<ngoodhits; ihit++ ){
    int pmtnum = GetHit(ihit)->GetPMTNum();
    
    pmt_row[pmtnum] = GetHit(ihit)->GetRow();
    pmt_col[pmtnum] = GetHit(ihit)->GetCol();
    xpmt[pmtnum] = GetHit(ihit)->GetX();
    ypmt[pmtnum] = GetHit(ihit)->GetY();
    tpmt[pmtnum] = GetHit(ihit)->GetTimeCorrected(); //includes trigger phase correction
    
    hitindex_pmts[pmtnum] = ihit;
    //hitused_pmts[pmtnum] = false;
    clustindex_pmts[pmtnum] = -1;

    bool keep = true;
    if( fRefTimeIsSet && fUseRefTimeDiffCut ){
      keep = fabs( tpmt[pmtnum] - fRefTime - fGRINCH_dToffset ) <= fMaxTdiffRefTime;
    }

    if( keep ){
      //list_pmts.insert( pmtnum );
      unused_pmts.insert( pmtnum );
    }
  }

  nclust = 0; //for readability

  //Start "new" improved clustering algorithm here:
  while( !( unused_pmts.empty() ) ){
    //grab the first unused PMT: 
    auto ipmt = *(unused_pmts.begin()); 
    //Each iteration of this loop is guaranteed to remove at least one pmt from the set otherwise we would get stuck in an infinite loop. The first pmt in the set is always either used to create a new cluster or added to an existing one!

    //Check if this pmt is not already in a cluster:
    //if( unused_pmts.find(ipmt) != unused_pmts.end() ){ //compare this pmt to all existing clusters. 
    if( nclust == 0 ){ //create first cluster
      std::set<Int_t> pmtlist_temp;
      pmtlist_temp.insert( ipmt );
      pmtlist_clusters.push_back( pmtlist_temp );
      clustindex_pmts[ipmt] = nclust;
      unused_pmts.erase( ipmt );

      tmean_clusters.push_back( tpmt[ipmt] );
      
      //Now loop on all OTHER unused PMTs and add any immediate neighbors of ipmt found to this cluster:
      for( auto jpmt : unused_pmts ){
	//if( jpmt != ipmt ){ //In principle this check is unnecessary since we already erased ipmt
	Double_t sep2 = pow( xpmt[ipmt] - xpmt[jpmt], 2 ) + pow( ypmt[ipmt] - ypmt[jpmt], 2 );

	Double_t dt = tpmt[jpmt] - tmean_clusters[nclust];
	
	if( sep2 <= fMaxSep2 && fabs( dt ) <= fMaxTdiffClust ){ //in-time neighbor!
	  //update candidate cluster mean time as PMTs are added:
	  double oldtmean = tmean_clusters[nclust];
	  double npmt_clust = pmtlist_clusters[nclust].size();
	  tmean_clusters[nclust] = (oldtmean * npmt_clust + tpmt[jpmt])/(npmt_clust+1.0);
	    
	  pmtlist_clusters[nclust].insert( jpmt );
	  //unused_pmts.erase( jpmt ); DON'T modify the set while iterating on it, this could have unintended consequences
	  clustindex_pmts[jpmt] = nclust;
	}
      }

      // NOW go back over all the PMTs in the cluster and erase them from the unused PMTs list. This is a safe operation
      for( auto jpmt : pmtlist_clusters[nclust] ){
	unused_pmts.erase(jpmt);
      }
      nclust++; //increment cluster count. IMPORTANT!
    } else { //nclust > 0. loop on all existing clusters. Compare this PMT to all PMTs in all existing clusters, add it to the first cluster containing an in-time neighbor of this PMT according to MaxSep!
      int iclust=0; 
      while( unused_pmts.find( ipmt ) != unused_pmts.end() && iclust < nclust ){
	//for( int iclust=0; iclust<nclust; iclust++ ){
	for( auto jpmt : pmtlist_clusters[iclust] ){ //loop on all PMTs in the cluster and check for a match:
	  Double_t sep2 = pow( xpmt[ipmt]-xpmt[jpmt],2 ) + pow( ypmt[ipmt]-ypmt[jpmt],2);
	  Double_t dt = tpmt[ipmt]-tmean_clusters[iclust];
	  
	  if( sep2 <= fMaxSep2 && fabs( dt ) <= fMaxTdiffClust ){ //iclust contains a neighbor of ipmt, add ipmt to iclust!  
	    double oldtmean = tmean_clusters[iclust];
	    double npmt_clust = pmtlist_clusters[iclust].size();
	    tmean_clusters[iclust] = (oldtmean * npmt_clust + tpmt[ipmt])/(npmt_clust+1.0);
	    
	    pmtlist_clusters[iclust].insert( ipmt );
	    clustindex_pmts[ipmt] = iclust;
	    unused_pmts.erase(ipmt);
	    break; //if we found a match, we don't need to iterate any more
	  }
	} 
	iclust++; //increment cluster index regardless of match. If no match found in existing clusters then this will eventually force us to exit the loop 
      }
     
      //check whether the current pmt was matched to any existing cluster:
      if( unused_pmts.find( ipmt ) != unused_pmts.end() ){ //this PMT was not matched to any existing cluster; start a new one!
	std::set<Int_t> pmtlist_temp;
	pmtlist_temp.insert( ipmt );
	pmtlist_clusters.push_back( pmtlist_temp );
	clustindex_pmts[ipmt] = nclust; 
	unused_pmts.erase( ipmt );

	tmean_clusters.push_back( tpmt[ipmt] );

	//after creating a new cluster, loop on all remaining unused pmts and check for immediate neighbors of this one! 
	for( auto jpmt : unused_pmts ){
	//if( jpmt != ipmt ){ //In principle this check is unnecessary since we already erased ipmt
	  Double_t sep2 = pow( xpmt[ipmt] - xpmt[jpmt], 2 ) + pow( ypmt[ipmt] - ypmt[jpmt], 2 );
	  Double_t dt = tpmt[jpmt] - tmean_clusters[nclust];
	  
	  if( sep2 <= fMaxSep2 && fabs( dt ) <= fMaxTdiffClust ){ //neighbor!
	    double oldtmean = tmean_clusters[nclust];
	    double npmt_clust = pmtlist_clusters[nclust].size();
	    tmean_clusters[nclust] = (oldtmean * npmt_clust + tpmt[jpmt])/(npmt_clust+1.0);
	    pmtlist_clusters[nclust].insert( jpmt );
	    //unused_pmts.erase( jpmt ); DON'T modify the set while iterating on it, this could have unintended consequences
	    clustindex_pmts[jpmt] = nclust;
	  }
	}
	
	//go back over all pmts added to the current cluster and erase them from the unused pmts list
	for( auto jpmt : pmtlist_clusters[nclust] ){
	  unused_pmts.erase(jpmt);
	} 
	nclust++; //increment cluster count! IMPORTANT!
      } //end check on need to start a new cluster
    } //end check on any existing clusters to start first cluster
  } //end while ( !(unused_pmts.empty()))


  //ALL PMTs will be part of at least one "cluster". Later we will also consider timing here; for now we allow clusters to grow to arbitrary size in any direction; we'll probably want to revisit that later:
  // if( ngoodhits > 0 ){
  //   for( auto ipmt : list_pmts ){
  //     nclust = pmtlist_clusters.size();
  //     //if( !hitused_pmts[ipmt] ){ //find a cluster to put this PMT into!
  //     // Check existing clusters to see if this is a nearest-neighbor of 
  //     // any PMT in an existing cluster:
  //     Double_t xpmt = GetHit(hitindex_pmts[ipmt])->GetX();
  //     Double_t ypmt = GetHit(hitindex_pmts[ipmt])->GetY();
      
  //     if( nclust == 0 ){ //first pmt starts a cluster
  // 	std::set<Int_t> pmtlist_temp;
  // 	pmtlist_temp.insert(ipmt);
  // 	pmtlist_clusters.push_back( pmtlist_temp );
  // 	clustindex_pmts[ipmt] = 0;
  //     } else {	
  // 	//if( clustindex_pmts[ipmt] < 0 ){ //if this pmt is not already part of a cluster, check it against existing clusters:
  // 	for( int iclust=0; iclust<nclust; iclust++ ){ //loop over all existing clusters
  // 	  bool isneighbor = false;
  // 	  for( int ihit=0; ihit<pmtlist_clusters[iclust].size(); ihit++ ){ //loop over pmts in the cluster:
  // 	    //Check if this PMT is within maxsep of any hit in an existing cluster:
  // 	    Double_t xtest = GetHit(hitindex_pmts[pmtlist_clusters[iclust][ihit]])->GetX();
  // 	    Double_t ytest = GetHit(hitindex_pmts[pmtlist_clusters[iclust][ihit]])->GetY();
	    
  // 	    Double_t sep2 = pow( xpmt - xtest, 2 ) + pow( ypmt - ytest, 2 );
  // 	    if( sep2 <= fMaxSep2  ){ // this is a neighbor to at least one of the hits in an existing cluster! --> Add it!
  // 	      clustindex_pmts[ipmt] = iclust; 
  // 	      //hitused_pmts[ipmt] = true;
  // 	      isneighbor = true;
  // 	      break; //exit the loop, we don't need to consider any more pmts
  // 	    }
  // 	  }
  // 	  if( isneighbor ){
  // 	    pmtlist_clusters[iclust].push_back( ipmt );
  // 	    break; //this insures that every PMT is part of one and ONLY one cluster; 
  // 	    // the first cluster containing one of its neighbors.
  // 	  }
  // 	}
  //     }
  //     // 
  //     if( clustindex_pmts[ipmt] < 0 ){ //PMT is neither the first pmt nor the neighbor to any existing cluster; create a new one!
  // 	std::vector<Int_t> pmtlist_temp(1,ipmt);
  // 	pmtlist_clusters.push_back( pmtlist_temp );
  //     }
  //   }

  nclust = pmtlist_clusters.size();
  
  for( int iclust=0; iclust<nclust; iclust++ ){
    new( (*fClusters)[iclust] ) SBSCherenkov_Cluster();
    SBSCherenkov_Cluster *clusttemp = GetCluster(iclust);
    
    //for( int ihit=0; ihit<pmtlist_clusters[iclust].size(); ihit++ ){
    for( auto ipmt : pmtlist_clusters[iclust] ){
      //insert the hit into the cluster. This automatically updates the cluster properties
      SBSCherenkov_Hit *hittemp = GetHit( hitindex_pmts[ipmt] );
      
      clusttemp->Insert( hittemp );
      
      hittemp->SetClustIndex( iclust );
    }
  }

  return 0;
}

//_____________________________________________________________________________
Int_t SBSGRINCH::MatchClustersWithTracks( TClonesArray& tracks )
{
  // what we actually want to do at this stage is to find the closest TRACK to each GRINCH CLUSTER 
  // consistent with the cuts; and assign a track to each cluster that is "close enough" to the track:
  // Later we will choose the biggest one that has a track match.

  Int_t ntracks = tracks.GetLast() + 1; 
  Int_t nclust = GetNumClusters();

  Int_t nmatch = 0;

  for( int iclust=0; iclust<nclust; iclust++ ){
    //Int_t nmatch=0; 
    Int_t itrmin=-1;
    Int_t bestmirror=-1;

    Double_t dxbest = kBig;
    Double_t dybest = kBig;
    Double_t dtbest = kBig;
    
    SBSCherenkov_Cluster *clusttemp = GetCluster( iclust );

    double xGRINCH = clusttemp->GetXcenter();
    double yGRINCH = clusttemp->GetYcenter();
    
    double tGRINCH = clusttemp->GetMeanTimeCorrected();

    double dtGRINCH = tGRINCH - fGRINCH_dToffset; //THIS assumes GRINCH time to be aligned at offset
    if( fRefTimeIsSet ) dtGRINCH -= fRefTime; //"ref time" is ordinarily from the shower
    
    for( int itrack=0; itrack<ntracks; itrack++ ){
      auto* theTrack = static_cast<THaTrack*>( tracks.At(itrack) ); 
   
      // Now we need to come up with some sort of matching criterion for GRINCH clusters and tracks. 
      // I suppose what we want for now is the cluster closest to the track along X and/or Y This could become more 
      // sophisticated over time.
      // Also consider timing!
 
      //Project track to GRINCH entrance window: 
      double zGRINCH = GetOrigin().Z();
      
      double xtrack = theTrack->GetX(zGRINCH);
      double ytrack = theTrack->GetY(zGRINCH);
      double thetatrack = theTrack->GetTheta();
      double phitrack = theTrack->GetPhi();
      double ptrack = theTrack->GetP();
      
      double minxdiff = 0.0; 
      double mindiff2 = 0.0;
      //int iclmin=-1;

      //    for( int iclust=0; iclust<nclust; iclust++ ){
      //SBSCherenkov_Cluster *clusttemp = ( (SBSCherenkov_Cluster*) (*fClusters)[iclust] );	

      if( fabs(dtGRINCH) <= fMaxTdiffRefTime || !fUseRefTimeDiffCut ){ //don't consider clusters with bad GRINCH timing:
	
	for( int imirr=0; imirr<fNmirror; imirr++ ){
	  if( xtrack >= fTrackMatchXmin[imirr] && xtrack < fTrackMatchXmax[imirr] ){
	    //double xdiff =  fabs(xtrack - fTrackMatchXslope[imirr] * xGRINCH - fTrackMatchX0[imirr]-fTrackMatchPslope/ptrack);
	    //Implement new x-matching criterion parametrized based on track theta:
	    double xdiff = xGRINCH - xtrack - ( fTrackMatchX0[imirr] + fTrackMatchXslope[imirr]*thetatrack );
	    //double ydiff = fabs(ytrack - fTrackMatchYslope[imirr] * yGRINCH - fTrackMatchY0[imirr]);
	    
	    //The lines below implement Maria's 3rd-order polynomial in track phi for the y correlation:
	    double ypredict = fTrackMatchY0[imirr];
	    for( int ipar=1; ipar<fOrderTrackMatchY; ipar++ ){
	      ypredict += fTrackMatchYslope[3*imirr + ipar-1] * pow( phitrack, ipar );
	    }
	    
	    double ydiff = yGRINCH - ytrack - ypredict; 
	    
	    double diff2 = pow( xdiff/fTrackMatchXsigma[imirr], 2 ) + pow( ydiff/fTrackMatchYsigma[imirr], 2 );

	    if( fUseRefTimeDiffCut ) diff2 += pow(dtGRINCH/fGRINCH_dTsigma,2);
	    
	    //     if( xdiff <= fTrackMatchXcut ){
	    if( diff2 <= pow( fTrackMatchNsigmaCut, 2 ) ){
	      //if( iclmin < 0 || xdiff < minxdiff ){
	      if( itrmin < 0 || diff2 < mindiff2 ){
		itrmin = itrack;
		bestmirror = imirr+1; //keep Maria's numbering convention of 1 to 4
		//  minxdiff = xdiff;
		mindiff2 = diff2;
		
		dxbest = xdiff;
		dybest = ydiff;
		dtbest = dtGRINCH;
	      }
	    }
	  }
	}
      }
    }

    if( itrmin >= 0 ){ //match:
      //SBSCherenkov_Cluster *clusttemp = ( (SBSCherenkov_Cluster*) (*fClusters)[iclmin] );
      //SBSCherenkov_Cluster *clusttemp = GetCluster( iclmin );

      auto* theTrack = static_cast<THaTrack*>( tracks.At(itrmin) ); 

      clusttemp->SetTrackIndex( itrmin ); 
      clusttemp->SetTrack( theTrack );
      clusttemp->SetMirrorIndex( bestmirror );

      clusttemp->SetDX( dxbest );
      clusttemp->SetDY( dybest );
      
      //We also need to loop on all the hits and set the track index for those hits:
      TIter next(clusttemp->GetHitList());
      SBSCherenkov_Hit *hittemp;
      while( ( hittemp = static_cast<SBSCherenkov_Hit*>( next() ) ) ){
	hittemp->SetTrackIndex( itrmin );
      }

      nmatch++;
    }
    
  }

  // here will be implemented the matching of GRINCH clusters with tracks
  return nmatch;
}

Int_t SBSGRINCH::SelectBestCluster(Int_t nmatch){
  Int_t nclusters = GetNumClusters(); 
  
  Int_t best=-1;
  Int_t maxsize=0;
  for( int iclust=0; iclust<nclusters; iclust++ ){
    SBSCherenkov_Cluster *clusttemp = GetCluster(iclust);

    if( clusttemp->GetTrackIndex() >= 0 || nmatch <= 0 ){ //either this cluster is matched with a track or there are no clusters matched with tracks and we just pick the biggest one:
      if( best < 0 || clusttemp->GetNHits() > maxsize ){
	best = iclust;
	maxsize = clusttemp->GetNHits();
	//fBestCluster = clusttemp;
      }
    }
  }

  return best;
}

void SBSGRINCH::SetRefTime( Double_t RefTime ){
  fRefTime = RefTime;
  fRefTimeIsSet = true;
}

ClassImp(SBSGRINCH)

