// *-- Author :    Guido Maria Urciuoli   12 March 2001

//////////////////////////////////////////////////////////////////////////
//
// SBSCherenkovDetector
//
// The RICH detector
// Written by Guido Maria Urciuoli, INFN
// Adapted for Hall A Analyzer by Ole Hansen, JLab
// Adapted to SBS by Seamus Riordan, ANL and Eric Fuchey, UConn
//
//////////////////////////////////////////////////////////////////////////

#include "TMath.h"
#include "SBSCherenkovDetector.h"
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
SBSCherenkovDetector::SBSCherenkovDetector( const char* name, const char* description, 
		  THaApparatus* apparatus ) :
  SBSGenericDetector(name,description,apparatus), 
  fDoResolve(false), fDoTimeFilter(false)//, 
  //fTrackX(kBig), fTrackY(kBig)
{
  //keep this line first
  fBench = new THaBenchmark;

  // Normal constructor with name and description

  fHits             = new TClonesArray("SBSCherenkov_Hit",500);
  fClusters         = new TClonesArray("SBSCherenkov_Cluster",50);
  
  Clear();
  //fDebug=1;
}

//_____________________________________________________________________________
SBSCherenkovDetector::~SBSCherenkovDetector()
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
void SBSCherenkovDetector::Clear( Option_t* opt )
{
  // Reset event-by-event data
  if(fDebug)cout << "SBSCherenkovDetector::Clear() " << endl; 
  if( fDoBench ) fBench->Begin("Clear");
  SBSGenericDetector::Clear(opt);
  if(fDebug)cout << "Clear hits() " << endl; 
  fHits->Clear();
  //fResolvedHits->Clear();
  DeleteClusters();
  if( fDoBench ) fBench->Stop("Clear");
}

//_____________________________________________________________________________
Int_t SBSCherenkovDetector::ReadDatabase( const TDatime& date )
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

  Int_t err = SBSGenericDetector::ReadDatabase(date);
  if(err) {
    return err;
  }
  fIsInit = false;

  fHit_tmin.clear();
  fHit_tmax.clear();
  fAmpToTCoeff.clear();

  std::vector<Double_t> hit_tmin;
  std::vector<Double_t> hit_tmax;
  std::vector<Double_t> amp_tot_coeffs;

  DBRequest config_request[] = {
    { "hit_mintime",           &hit_tmin,   kDoubleV, 0, false }, 
    { "hit_maxtime",           &hit_tmax,   kDoubleV, 0, false }, 
    { "amp_tot_coeffs",        &amp_tot_coeffs,   kDoubleV, 0, true }, 
    { 0 } ///< Request must end in a NULL
  };
  err = LoadDB( fi, date, config_request, fPrefix );
  
  if(hit_tmin.size()==fNelem){
    for(int i = 0; i<hit_tmin.size(); i++){
      fHit_tmin.push_back(hit_tmin[i]);
    }
  }else{
    fHit_tmin.push_back(hit_tmin[0]);
  }
  
  if(hit_tmax.size()==fNelem){
    for(int i = 0; i<hit_tmax.size(); i++){
      fHit_tmax.push_back(hit_tmax[i]);
    }
  }else{
    fHit_tmax.push_back(hit_tmax[0]);
  }

  if(amp_tot_coeffs.size()==fNelem){
    for(int i = 0; i<amp_tot_coeffs.size(); i++){
      fAmpToTCoeff.push_back(amp_tot_coeffs[i]);
    }
  }else{
    fAmpToTCoeff.push_back(1.0);
  }
  
  fIsInit = true;
  
  fclose(fi);
  return kOK;
}

//_____________________________________________________________________________
Int_t SBSCherenkovDetector::DefineVariables( EMode mode )
{
  // Define (or delete) global variables of the detector
  
  if( mode == kDefine && fIsSetup ) return kOK;
  //fIsSetup = ( mode == kDefine );

  Int_t err = SBSGenericDetector::DefineVariables(mode);
  if(err)
    return err;
  
  //Hits hits
  RVarDef var1[] = {
    { "ngoodhits",  " number of PMT hits", "GetNumHits()"                       },
    { "hit.pmtnum", " Hit PMT num",        "fHits.SBSCherenkov_Hit.GetPMTNum()" },
    { "hit.row",    " PMT hit row",        "fHits.SBSCherenkov_Hit.GetRow()"    },
    { "hit.col",    " PMT hit column",     "fHits.SBSCherenkov_Hit.GetCol()"    },
    { "hit.xhit",   " PMT hit X",          "fHits.SBSCherenkov_Hit.GetX()"      },
    { "hit.yhit",   " PMT hit y",          "fHits.SBSCherenkov_Hit.GetY()"      },
    { "hit.amp",    " PMT hit ampliutude", "fHits.SBSCherenkov_Hit.GetAmp()"    },
    { "hit.time",   " PMT hit time",       "fHits.SBSCherenkov_Hit.GetTime()"   },
    { 0 }
  };
  DefineVarsFromList( var1, mode, "" );// (re)define path here...
  
  // RVarDef var2[] = {
  //   { "nclus",        " number of GRINCH PMT clusters",  "GetNumClusters()"                                 },
  //   { "clus.size",    " GRINCH cluster size",            "fClusters.SBSCherenkov_Cluster.GetNHits()"           },
  //   { "clus.x_mean",  " GRINCH cluster X center",        "fClusters.SBSCherenkov_Cluster.GetXcenter()"         },
  //   { "clus.y_mean",  " GRINCH cluster Y center",        "fClusters.SBSCherenkov_Cluster.GetYcenter()"         },
  //   { "clus.tr_mean", " GRINCH cluster mean lead time",  "fClusters.SBSCherenkov_Cluster.GetMeanRisingTime()"  },
  //   { "clus.tf_mean", " GRINCH cluster mean trail time", "fClusters.SBSCherenkov_Cluster.GetMeanFallingTime()" },
  //   { "clus.adc",     " GRINCH cluster total charge",    "fClusters.SBSCherenkov_Cluster.GetCharge()"          },
  //   { 0 }
  // };
  // DefineVarsFromList( var2, mode, "" );// (re)define path here...
  
  return kOK;
}

//_____________________________________________________________________________
Int_t SBSCherenkovDetector::Decode( const THaEvData& evdata )
{
  //Decode RICH data and fill hit array
  if(fDebug){
    cout << "SBSCherenkovDetector::Decode " << endl;
    cout << " Is Init ? " << fIsInit << " Is Physics Trigger ? " << evdata.IsPhysicsTrigger() << endl;
  }
  if( !fIsInit ) return -255;
  if( !evdata.IsPhysicsTrigger() ) return -1;

  Clear();
  
  if( fDoBench ) fBench->Begin("Decode");
  
  SBSGenericDetector::Decode(evdata);
  
  if(fDebug)cout << "Finished decoding" << endl;
  
  if( fDoBench ) fBench->Stop("Decode");
  return GetNumHits();
}

//_____________________________________________________________________________
Int_t SBSCherenkovDetector::CoarseProcess( TClonesArray& tracks )
{
  //
  if(fDebug)cout << "Begin Coarse Process" << endl;
  if( fDoBench ) fBench->Begin("CoarseProcess");
 
  SBSGenericDetector::CoarseProcess(tracks);

  double amp, x, y;
  double tmin, tmax;
  
  Int_t nHit = 0;
  SBSCherenkov_Hit* the_hit = nullptr;
    
  for(int k = 0; k<fNGoodTDChits; k++){
    if(fHit_tmin.size()==fNelem){
      tmin = fHit_tmin[fGood.TDCelemID[k]];
    }else{
      tmin = fHit_tmin[0];
    }
    if(fHit_tmax.size()==fNelem){
      tmax = fHit_tmax[fGood.TDCelemID[k]];
    }else{
      tmax = fHit_tmax[0];
    }
    
    if(tmin<=fGood.t[k] && fGood.t[k]<=tmax){
      the_hit = new( (*fHits)[nHit++] ) SBSCherenkov_Hit();
      
      the_hit->SetPMTNum(fGood.TDCelemID[k]);
      the_hit->SetRow(fGood.TDCrow[k]);
      the_hit->SetCol(fGood.TDCcol[k]);
      the_hit->SetTime(fGood.t[k]);
      
      x = (fElements[fGood.TDCelemID[k]])->GetX();
      y = (fElements[fGood.TDCelemID[k]])->GetY();
      if(k<fAmpToTCoeff.size()){
	amp = fGood.t_ToT[k]*fAmpToTCoeff[k];
      }else{// we should be guaranteed that the array has at least one element
	amp = fGood.t_ToT[k]*fAmpToTCoeff[0];
      }
      
      the_hit->SetX(x);
      the_hit->SetY(y);
      the_hit->SetAmp(amp);
    }
  }
  
  if( fDoBench ) fBench->Stop("CoarseProcess");
  if(fDebug)cout << "End Coarse Process" << endl;
  return 0;
}


//_____________________________________________________________________________
Int_t SBSCherenkovDetector::FineProcess( TClonesArray& tracks )
{
  if(fDebug)cout << "Begin Fine Process" << endl;
  // The main RICH processing method. Here we
  //  
  // - Identify clusters
  // - Attempt to match a traceto these
  // - Calculate particle probabilities

  if( !fIsInit ) return -255;

  // Clusters reconstructed here
  if( FindClusters() == 0 ) { 
    return -1;
  }
  
  if(fDebug)cout << "DoTimeFitter ? " << fDoTimeFilter << " Tracks::GetLast ? " << tracks.GetLast() << endl;

  // if( fDoResolve )
  //   ResolveClusters();
  //cout<<"stay and stop here!!!"<<endl;
  if( fDoBench ) fBench->Begin("FineProcess");

  //if(fDoTimeFilter) CleanClustersWithTime();
  
  // Clusters matched with tracks here (obviously if there are any tracks to match)
  if(tracks.GetLast()>0){
    //MatchClustersWithTracks(tracks);
  }
  
  if( fDoBench ) fBench->Stop("FineProcess");
  if(fDebug)cout << "Done fine process " << endl;
  return 0;
}


//__________________________________________________________________________
void SBSCherenkovDetector::DeleteClusters()
{
  //Delete all clusters
  if(fDebug)cout << "Clear Clusters" << endl;
  fClusters->Clear("C");
  //fResolvedClusters->Clear("C");
  return;
}

//__________________________________________________________________________
Int_t SBSCherenkovDetector::FindClusters()
{
  DeleteClusters();

  if( fDoBench ) fBench->Begin("FindClusters");

  Int_t nClust = GetNumClusters();
  
  if(fDebug)cout << "Finished clustering " << endl;
  if( fDoBench ) fBench->Stop("FindClusters");
  return nClust;
}

//_____________________________________________________________________________
// Int_t SBSCherenkovDetector::MatchClustersWithTracks( TClonesArray& tracks )
// {
//   Int_t Nassociated = 0;
//   return(Nassociated);
// }

// //__________________________________________________________________________
// Int_t SBSCherenkovDetector::CleanClustersWithTime()// pass as an option...
// {
  
//   return(1);
// }

//_____________________________________________________________________________
void SBSCherenkovDetector::PrintBenchmarks() const
{
  // Print benchmark results
  
  if( !fDoBench )
    return;

  fBench->Print("Clear");
  fBench->Print("Decode");
  fBench->Print("FineProcess");
  fBench->Print("FindClusters");
  
  //cout << endl << "Breakdown of time spent in FineProcess:" << endl;
  
  return;
}

ClassImp(SBSCherenkovDetector)

