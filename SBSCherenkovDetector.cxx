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
  SetModeTDC(SBSModeTDC::kTDC); //  A TDC with leading & trailing edge info
  SetModeADC(SBSModeADC::kNone); // Default is No ADC, but can be re-enabled later
  //keep this line first
  fBench = new THaBenchmark;

  // Normal constructor with name and description

  fHits             = new TClonesArray("SBSCherenkov_Hit",500);
  fClusters         = new TClonesArray("SBSCherenkov_Cluster",50);

  //fEmptyCluster = new SBSCherenkov_Cluster(); //dummy instance of cluster to initialize best cluster pointer for events with no found clusters
  //  fBestCluster = nullptr;
  //fEmptyCluster->SetXcenter(kBig);
  //fEmptyCluster->SetYcenter(kBig);

  //The above lines are not necessary as we have other handles to identify events with no clusters found

  //fBestCluster = fEmptyCluster;

  //set default timing cuts wide open, these values are somewhat arbitrary:
  fHit_tmin = -1000.0;
  fHit_tmax = 4000.0; 
  
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
  fHits->Clear("C");
  //fResolvedHits->Clear();
  DeleteClusters();
  
  fNtrackMatch = 0;

  fBestCluster.Clear("F");

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

  std::vector<Double_t> xpos,ypos;
  std::vector<Double_t> amp_tot_coeffs;

  DBRequest config_request[] = {
    { "xpos", &xpos,    kDoubleV, 0, 1 },
    { "ypos", &ypos,    kDoubleV, 0, 1 },
    { "hit_mintime",           &fHit_tmin,   kDouble, 0, 1 }, 
    { "hit_maxtime",           &fHit_tmax,   kDouble, 0, 1 }, 
    { "amp_tot_coeffs",        &amp_tot_coeffs,   kDoubleV, 0, true }, 
    { 0 } ///< Request must end in a NULL
  };
  err = LoadDB( fi, date, config_request, fPrefix );
  
  fAmpToTCoeff.clear();
  if(int(amp_tot_coeffs.size())==fNelem){
    for(size_t i = 0; i<amp_tot_coeffs.size(); i++){
      fAmpToTCoeff.push_back(amp_tot_coeffs[i]);
    }
  }else{
    fAmpToTCoeff.push_back(1.0);
  }
  
  
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
    { "hit.amp",    " PMT hit amplitude", "fHits.SBSCherenkov_Hit.GetAmp()"    },
    { "hit.time",   " PMT hit time",       "fHits.SBSCherenkov_Hit.GetTime()"   },
    { "hit.clustindex", " Index of cluster to which this hit belongs", "fHits.SBSCherenkov_Hit.GetClustIndex()" },
    { "hit.trackindex", " Index of track to which this hit belongs", "fHits.SBSCherenkov_Hit.GetTrackIndex()" },
    { "ntrackmatch", "Number of track-matched clusters", "fNtrackMatch" },
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
  //double tmin, tmax;
  
  Int_t nHit = 0;
  SBSCherenkov_Hit* the_hit = nullptr;

  fNtrackMatch = 0;
  
  for(int k = 0; k<fNGoodTDChits; k++){
    //tmin = -fElements[fGood.TDCelemID[k]]->TDC()->GetGoodTimeCut();
    //tmax = +fElements[fGood.TDCelemID[k]]->TDC()->GetGoodTimeCut();
    
    //double t0 = fElements[fGood.TDCelemID[k]]->TDC()->GetGoodTimeCut();

    //    if(tmin<=fGood.t[k] && fGood.t[k]<=tmax){
    if( fHit_tmin <= fGood.t[k] && fGood.t[k] <= fHit_tmax ){
      the_hit = new( (*fHits)[nHit++] ) SBSCherenkov_Hit();
      
      the_hit->SetPMTNum(fGood.TDCelemID[k]);
      the_hit->SetRow(fGood.TDCrow[k]);
      the_hit->SetCol(fGood.TDCcol[k]);
      the_hit->SetTime(fGood.t[k]);
      
      x = (fElements[fGood.TDCelemID[k]])->GetX();
      y = (fElements[fGood.TDCelemID[k]])->GetY();
      if(k<(int)fAmpToTCoeff.size()){
	amp = fGood.t_ToT[k]*fAmpToTCoeff[k];
      }else{// we should be guaranteed that the array has at least one element
	amp = fGood.t_ToT[k]*fAmpToTCoeff[0];
      }
      
      the_hit->SetX(x);
      the_hit->SetY(y);
      the_hit->SetAmp(amp);
    }
  }
  //clustering to be done by dereived class...
  
  if( fDoBench ) fBench->Stop("CoarseProcess");
  if(fDebug)cout << "End Coarse Process" << endl;
  return 0;
}


//_____________________________________________________________________________
Int_t SBSCherenkovDetector::FineProcess( TClonesArray& tracks )
{
  // fine processing like association with tracks belong to 
  // derived classes such as SBSGRINCH
  return 0;
}


//__________________________________________________________________________
void SBSCherenkovDetector::DeleteClusters()
{
  //Delete all clusters
  if(fDebug)cout << "Clear Clusters" << endl;
  fClusters->Clear("C");
  //fResolvedClusters->Clear("C");
}

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
  
}

ClassImp(SBSCherenkovDetector)

