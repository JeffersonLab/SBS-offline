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
  
  Clear();
}

//_____________________________________________________________________________
SBSGRINCH::~SBSGRINCH()
{
   // Destructor. Remove variables from global list and free up the memory
  // allocated by us.
  Clear();// so the prgram doesn't complain when deleting clusters
  RemoveVariables();
  delete fHits;
  // delete fResolvedHits;
  delete fClusters;
  // delete fResolvedClusters;
  // delete [] fMIPs;
  // delete [] fXseg;
  delete fBench;
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

  /*
  //Add here the cluster output variables? maybe
  RVarDef var2[] = {
    { "nclus",        " number of GRINCH PMT clusters",  "GetNumClusters()"                                  },
    { "clus.size",    " GRINCH cluster size",            "fClusters.SBSCherenkov_Cluster.GetNHits()"            },
    { "clus.x_mean",  " GRINCH cluster X center",        "fClusters.SBSCherenkov_Cluster.GetXcenter()"          },
    { "clus.y_mean",  " GRINCH cluster Y center",        "fClusters.SBSCherenkov_Cluster.GetYcenter()"          },
    { "clus.tr_mean", " GRINCH cluster mean time",       "fClusters.SBSCherenkov_Cluster.GetMeanTime()"         },
    { "clus.tf_mean", " GRINCH cluster mean time over threshold", "fClusters.SBSCherenkov_Cluster.GetMeanToT()" },
    { "clus.adc",     " GRINCH cluster total charge",    "fClusters.SBSCherenkov_Cluster.GetCharge()"           },
    { 0 }
  };
  DefineVarsFromList( var2, mode, "" );// (re)define path here...
  */
  
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
  
  // Clusters matched with tracks here (obviously if there are any tracks to match)
  if(tracks.GetLast()>0){
    MatchClustersWithTracks(tracks);
  }
  
  if( fDoBench ) fBench->Stop("FineProcess");
  if(fDebug)cout << "Done fine process " << endl;
  
  return 0;
}


//__________________________________________________________________________
Int_t SBSGRINCH::FindClusters()
{
  // here will be implemented the clustering for the GRINCH
  return 0;
}

//_____________________________________________________________________________
Int_t SBSGRINCH::MatchClustersWithTracks( TClonesArray& tracks )
{
  // here will be implemented the matching of GRINCH clusters with tracks
  return 0;
}

ClassImp(SBSGRINCH)

