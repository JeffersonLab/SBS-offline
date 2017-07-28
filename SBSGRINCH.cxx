//*-- Author :    Guido Maria Urciuoli   12 March 2001

//////////////////////////////////////////////////////////////////////////
//
// SBSRICH
//
// The RICH detector
// Written by Guido Maria Urciuoli, INFN
// Adapted for Hall A Analyzer by Ole Hansen, JLab
// Adapted to SBS by Seamus Riordan, ANL
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
#if defined(HAS_SSTREAM) || (defined(__GNUC__)&&(__GNUC__ >= 3))
#include <sstream>
#define HAS_SSTREAM
#define ISTR istringstream
#else
#include <strstream>
#undef HAS_SSTREAM
#define ISTR istrstream
#endif
#include "THaBenchmark.h"

using namespace std;
// cout<<"stay and stop here!!!"<<endl;
//_____________________________________________________________________________
SBSGRINCH::SBSGRINCH( const char* name, const char* description, 
		  THaApparatus* apparatus ) :
  THaPidDetector(name,description,apparatus), 
  fMIPs(0),
  fMaxxMIP(100000), fMinxMIP(-100000), fMaxyMIP(100000), fMinyMIP(-100000),
  fDoResolve(false), fNseg(0), fXseg(0),
  fTrackX(kBig), fTrackY(kBig)
{
  //keep this line first
  fBench = new THaBenchmark;

  // Normal constructor with name and description

  fHits             = new TClonesArray("SBSGRINCH_Hit",1000);
  fClusters         = new TClonesArray("SBSGRINCH_Cluster",100);
  fResolvedHits     = new TClonesArray("SBSGRINCH_Hit",1000);
  fResolvedClusters = new TClonesArray("SBSGRINCH_Cluster",100);

  Clear();

}

//_____________________________________________________________________________
SBSGRINCH::~SBSGRINCH()
{
  // Destructor. Remove variables from global list and free up the memory
  // allocated by us.

  RemoveVariables();
  delete fHits;
  delete fResolvedHits;
  delete fClusters;
  delete fResolvedClusters;
  delete [] fMIPs;
  delete [] fXseg;
  delete fBench;
}

//_____________________________________________________________________________
void SBSGRINCH::Clear( Option_t* opt )
{
  // Reset event-by-event data

  if( fDoBench ) fBench->Begin("Clear");
  THaPidDetector::Clear(opt);
  fHits->Clear();
  fResolvedHits->Clear();
  DeleteClusters();
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

  // Open the database file
  FILE* fi = OpenFile( date );
  if( !fi ) return kFileError;
 
  // to be inserted in the data base 
  fiducial_zone_range = 0.05;
  cluster_distribution_sigma = 0.025;
  minimum_chi2_degree_of_freedom = 4;
  acceptable_chi2_prob = 0.01;
  clear_noise_trial_maximum_number = 5;
  epsilon = 1.;

  // Storage and default values for non-Double_t and non-member 
  // data from database. Note: These must be Double_t
  Double_t debug = 0, do_resolve = false,
    maxdist, hit_max_number=0, MIP_through_interception = 0;

  // Set up a table of tags to read and locations to store values.
  const DBRequest tags[] = { 
    { "l_emission", &l_emission, kDouble, 1, false},
    { "maxdist",    &maxdist, kDouble, 1, true},
    { "hit_max_number", &hit_max_number, kDouble, 1, true },
    { "MIP_through_interception", &MIP_through_interception, kDouble, 1, true },
    { "fiducial_zone_range", &fiducial_zone_range, kDouble, 1, true },
    { "cluster_distribution_sigma", &cluster_distribution_sigma, kDouble, 1, true  },
    { "acceptable_chi2_prob", &acceptable_chi2_prob, kDouble, 1, true  },
    { "minimum_chi2_degree_of_freedom", &minimum_chi2_degree_of_freedom, kDouble, 1, true  },
    { "clear_noise_trial_maximum_number", &clear_noise_trial_maximum_number, kDouble, 1, true },
    { "epsilon", &epsilon, kDouble, 1, true  },
    { "do_resolve", &do_resolve, kDouble, 1, true  },
    { "debug",      &debug, kDouble, 1, true  },
    { 0 }
  };

  // Read all the tag/value pairs. 
  Int_t err = LoadDB( fi, date, tags );

  // Complain and quit if any required values missing.
  if( err ) {    
    if( err>0 )
      Error( Here(here), "Required tag %s%s missing in the "
	     "database.\nDetector initialization failed.",
	     fPrefix, tags[err-1].name );
    fclose(fi);
    return kInitError;
  }

  // Now read the array data.
  TString line, tag;
  Int_t retval = kOK;
  int n;
  std::vector<double> ctr, xrotv, yrotv, zrotv, rad, quartz, gap;
  std::vector<int> npads;
  std::vector<double> padsize, xmip, ymip, size;

  Double_t xrot[] = {0,1,1}, yrot[] = {0,2,2}, zrot[] = {0,3,3};

  Double_t* rotdef[3] = { xrot, yrot, zrot };

  Int_t ntotal, nxpads, nypads;
  const DBRequest atags[] = {
    { "origin",       &ctr,      kDoubleV,   3, false},
    { "xrot",         &xrotv,     kDoubleV,   2, false},
    { "yrot",         &yrotv,     kDoubleV,   2, false},
    { "zrot",         &zrotv,     kDoubleV,   2, false},
    { "rad",          &rad,      kDoubleV,   4, false},
    { "quartz",       &quartz,   kDoubleV,   4, false},
    { "gap",          &gap,      kDoubleV,   3, false},
    { "npads",        &npads,    kIntV,      2, false},
    { "padsize",      &padsize,  kDoubleV,   2, false},
    { "xmip_range",   &xmip,     kDoubleV,   2, false},
    { "ymip_range",   &ymip,     kDoubleV,   2, false},
    { "size",         &size,     kDoubleV,   3, false},
    //    { "npads_adc",    npads_adc,  13, 2 },
    { 0 }
  };

  err = LoadDB( fi, date, atags, fPrefix );
  

  //--- Process the data

  if( npads.size() != 2) {
    Error( Here(here), "Illegal number of pads: nx = %d, ny = %d. "
	   "Detector initialization failed.", nxpads, nypads );
    fclose(fi);
    return kInitError;
  }

  // Check nelem for sanity
  nxpads = int(npads[0]+0.5);
  nypads = int(npads[1]+0.5);
  ntotal = nxpads*nypads;

  if( nxpads <= 0 || nypads <= 0 || ntotal <= 0 ) {
    Error( Here(here), "Illegal number of pads: nx = %d, ny = %d. "
	   "Detector initialization failed.", nxpads, nypads );
    fclose(fi);
    return kInitError;
  }
  // If the detector has already been initialized, the number of elements must not
  // change. This prevents problems with global variables that
  // point to memory that is dynamically allocated.
  if( fIsInit && ntotal != fNelem ) {
    Error( Here(here), "Cannot re-initalize with different number of pads. "
	   "(was: %d, now: %d). Detector not re-initialized.", fNelem, ntotal);
    fclose(fi);
    return kInitError;
  }
  fNelem  = ntotal;
  fNypads = nypads;
  // All ok - convert to non-Double_t types
  fDebug  = int(debug+0.5);
  fDoResolve  = bool(do_resolve);
  fMaxdist2   = maxdist*maxdist;
  fMIP_through_interception = int(MIP_through_interception);

  fMaxNumHits = fNelem;
  if (hit_max_number != 0)
  fMaxNumHits = int(hit_max_number);

  fIsInit = true;

  L_RAD = rad[0];
  n_radiator = rad[1];
  n_radiator_min = rad[2];
  n_radiator_max = rad[3];
  l_quartz = quartz[0];
  n_quartz = quartz[1];
  n_quartz_min = quartz[2];
  n_quartz_max = quartz[3];
  l_gap = gap[0];
  n_gap = gap[1];
  //  PADS_X = int(npads_adc[0]+0.5);
  //  PADS_Y = int(npads_adc[1]+0.5);
  PAD_SIZE_X = padsize[0];
  PAD_SIZE_Y = padsize[1];
  fMinxMIP = xmip[0];
  fMaxxMIP = xmip[1];
  fMinyMIP = ymip[0];
  fMaxyMIP = ymip[1];

  //FIXME: put in database

  // Geometry stuff. 
  // NB: The RICH depends heavily on correct geometry/orientation data.

  // Size of the detector. Currently only used by the event display, 
  // therefore optional
  for(int i=0; i<3; i++ )
    fSize[i] = float(size[i]);

  // Origin of RICH coordinate system. Include offsets here.
  fOrigin.SetXYZ(ctr[0],ctr[1],ctr[2]);

  ///////////////////////


  std::vector<double> rotset[3] = {xrotv, yrotv, zrotv};

  for( int j = 0; j < 3; j++ ){
      for( unsigned int i = 0; i < xrotv.size(); i++ ){
          rotdef[j][i] = rotset[j][i];
      }
  }

  // Define the x/y-plane vectors according to the rotations specified.
  // Put the rotations in the desired order
  for(int i=0; i<2; i++) {
    for(int j=i+1; j<3; j++) {
      if( (rotdef[j])[1] < (rotdef[i])[1] ) {
	Double_t* tmp = rotdef[i];
	rotdef[i] = rotdef[j];
	rotdef[j] = tmp;
      }
    }
  }
  // Carry out the rotations of the axis vectors
  fXax.SetXYZ(1,0,0);
  fYax.SetXYZ(0,1,0);
  for(int i=0; i<3; i++ ) {
    Double_t angle = (rotdef[i])[0]*TMath::Pi()/180.0;
    if( angle != 0.0) {
      int type = int((rotdef[i])[2]+0.5);
      switch(type) {
      case 1:
	fXax.RotateX(angle);
	fYax.RotateX(angle);
	break;
      case 2:
	fXax.RotateY(angle);
	fYax.RotateY(angle);
	break;
      case 3:
	fXax.RotateZ(angle);
	fYax.RotateZ(angle);
	break;
      }
    }
  }
  fZax = fXax.Cross(fYax);

  // Plane segmentation in x
  // The format of the line is: 
  //  <n_segments> <lowx1> <highx1> <offset1> <lowx2> <highx2> <offset2> ...
  // If <lowx> <= raw x_RICH < <highx>, then subtract the offset. 
  // Takes into account "the presence of not sensible materials between the 
  // planes". For positive offsets, this leads to overlapping regions of x 
  // coordinates.
  // Any raw x_RICH coordinates outside of the specified <low>-<high> ranges
  // are considered invalid.

  tag = "x_segments";
  err = LoadDBvalue( fi, date, tag, line );
  if( err == 0 && line.Length()>0 ) {
    ISTR inp(line.Data());
    Int_t nseg;
    inp >> nseg;
    if( !inp )
      goto bad_data;
    fNseg = nseg;
    delete [] fXseg; fXseg = new Double_t[3*fNseg];
    for( int i=0; i<fNseg; i++ ) {
      int k = 3*i;
      inp >> fXseg[k] >> fXseg[k+1] >> fXseg[k+2];
      if( !inp )
	goto bad_data;
    }
  }

  // Read detector map
  // Multiple lines are allowed. The format is:
  //  detmap_1  <crate> <slot> <first channel> <last channel>
  //  detmap_2   dto.
  //  etc.

  fDetMap->Clear();
  n = 1;
  int status;
  while(1) {
    Int_t crate, slot, lo, hi;
    tag = TString(Form("detmap_%02d", n));
    if( status = LoadDBvalue( fi, date, tag, line ) != 0 )
      break;
    ISTR inp(line.Data());
    inp >> crate >> slot >> lo >> hi;
    if( !inp )
      goto bad_data;
    if( fDetMap->AddModule( crate, slot, lo, hi ) < 0 )
      goto bad_data;
    n++;
  }
  if(n>1)
    goto done;
 
 bad_data:    
  Error( Here(here), "Problem reading database entry \"%s\". "
	 "Detector initialization failed.", tag.Data() );
  retval = kInitError;

 done:
  fclose(fi);
  return retval;
}

//_____________________________________________________________________________
Int_t SBSGRINCH::DefineVariables( EMode mode )
{
  // Define (or delete) global variables of the detector
  
  if( mode == kDefine && fIsSetup ) return kOK;
  fIsSetup = ( mode == kDefine );

  RVarDef var1[] = {
    { "nhit",    "Number of pads with hits",    "GetNumHits()" },
    { "nrhit",   "Number of resolved hits",     "GetNumResolvedHits()" },
    { "nclust",  "Number of clusters",          "GetNumClusters()" },
    { "nresolv", "Number of resolved clusters", "GetNumResolvedClusters()" },
    // { "ngood",   "Number of clusters within cone",   "fNClustersGood" },
    // { "nextra",  "Number of clusters outside cone",  "fNClustersOut" },
    { "track.x", "Golden Track x in RICH plane (mm)",   "fTrackX" },
    { "track.y", "Golden Track y in RICH plane (mm)",   "fTrackY" },
    { 0 }
  };
  DefineVarsFromList( var1, mode );
  
  RVarDef var2[] = {
    { "hit.i",   "Hit i index (along x_rich)",       "fI" },
    { "hit.j",   "Hit j index (along y_rich)",       "fJ" },
    { "hit.x",   "Hit x position (mm)",              "fX" },
    { "hit.y",   "Hit y position (mm)",              "fY" },
    { "hit.adc", "Hit charge (ADC value)",           "fADC" },
    { "hit.flag","Hit status flag",                  "fFlag" },
    { "hit.veto","Hit local maximum veto flag",      "fVeto" },
    { 0 }
  };
  DefineVarsFromList( var2, mode, "fHits.SBSGRINCH_Hit." );
  
  RVarDef var3[] = {
    { "rhit.i",   "Resolved hit i index (along x_rich)",  "fI" },
    { "rhit.j",   "Resolved hit j index (along y_rich)",  "fJ" },
    { "rhit.x",   "Resolved hit x position (mm)",         "fX" },
    { "rhit.y",   "Resolved hit y position (mm)",         "fY" },
    { "rhit.adc", "Resolved hit charge (ADC value)",      "fADC" },
    { "rhit.flag","Resolved hit status flag",             "fFlag" },
    { 0 }
  };
  DefineVarsFromList( var3, mode, "fResolvedHits.SBSGRINCH_Hit." );

  RVarDef var4[] = {
    { "clus.x",     "Cluster x position (mm)",        "fXcenter" },
    { "clus.y",     "Cluster y position (mm)",        "fYcenter" },
    { "clus.chrg",  "Cluster charge (ADC sum)",       "fCharge"  },
    { "clus.size",  "Cluster size (no. of hits)",     "GetNHits()"},
    { "clus.angle", "Cluster Cherenkov angle (rad)",  "fAngle"   },
    { "clus.theta", "Cluster Theta photon",           "GetTheta_photon()" },
    { "clus.phi",   "Cluster Phi photon",             "GetPhi_photon()" },
    { "clus.mip",   "Cluster is MIP flag",            "fMIPflag" },
    { "clus.pionchi2",  "cluster chi square values according to the hypothesis the MIP is a pion","Getchi2_pi()"},
    { "clus.kaonchi2",  "cluster chi square values according to the hypothesis the MIP is a kaon","Getchi2_k()"},
    { "clus.protonchi2",  "cluster chi square values according to the hypothesis the MIP is a proton","Getchi2_p()"},
    { "clus.PionChi2AnalysisFlag",  "cluster involved in chi2 test for the pion hypothesis flag","GetPionChi2AnalysisFlag()"},
    { "clus.KaonChi2AnalysisFlag",  "cluster involved in chi2 test for the kaon hypothesis flag","GetKaonChi2AnalysisFlag()"},
    { "clus.ProtonChi2AnalysisFlag",  "cluster involved in chi2 test for the proton hypothesis flag","GetProtonChi2AnalysisFlag()"},
    { 0 }
  };
  DefineVarsFromList( var4, mode, "fClusters.SBSGRINCH_Cluster." );

  RVarDef var5[] = {
    { "rclus.x",     "Resolved cluster x position (mm)",    "fXcenter" },
    { "rclus.y",     "Resolved cluster y position (mm)",    "fYcenter" },
    { "rclus.chrg",  "Resolved cluster charge (ADC sum)",   "fCharge"  },
    { "rclus.size",  "Resolved cluster size (no. of hits)", "GetNHits()"},
    { "rclus.angle", "Resolved cluster Cherenkov angle (rad)", "fAngle"   },
    { "rclus.theta", "Resolved cluster Theta photon",   "GetTheta_photon()" },
    { "rclus.phi",   "Resolved cluster Phi photon",     "GetPhi_photon()" },
    { "rclus.mip",   "Resolved cluster is MIP flag",    "fMIPflag" },
    { "rclus.pionchi2",  "Resolved cluster chi square values according to the hypothesis the MIP is a pion","Getchi2_pi()"},
    { "rclus.kaonchi2",  "Resolved Cluster chi square values according to the hypothesis the MIP is a kaon","Getchi2_k()"},
    { "rclus.protonchi2",  "Resolved cluster chi square values according to the hypothesis the MIP is a proton","Getchi2_p()"},
    { "rclus.PionChi2AnalysisFlag",  "Resolved cluster involved in chi2 test for the pion hypothesis flag","GetPionChi2AnalysisFlag()"},
    { "rclus.KaonChi2AnalysisFlag",  "Resolved cluster involved in chi2 test for the kaon hypothesis flag","GetKaonChi2AnalysisFlag()"},
    { "rclus.ProtonChi2AnalysisFlag",  "Resolved cluster involved in chi2 test for the proton hypothesis flag","GetProtonChi2AnalysisFlag()"},
    { 0 }
  };
  DefineVarsFromList( var5, mode, "fResolvedClusters.SBSGRINCH_Cluster." );

  RVarDef var6[] = {
    { "mip.size",  "Golden Track MIP size",                 "GetNHits()" },
    { "mip.x",     "Golden Track MIP x (mm)",               "fXcenter" },
    { "mip.y",     "Golden Track MIP y (mm)",               "fYcenter" },
    { "mip.chrg",  "Golden Track MIP charge (ADC sum)",     "fCharge"  },
    { "mip.theta", "Theta MIP",   "GetTheta_photon()" },
    { "mip.phi",   "Phi MIP",     "GetPhi_photon()" },
    { "mip.fictious_flag",   "greater than 0 when the MIP is the interception between the track and the pad plane (equal to the MIP_through_interception value in the data base ",     " GetFictious_Mip_Flag()" },
    { "mip.nphot",  "Clusters in pi/K/p region of this MIP","N_Photon"  },
    { "mip.angles", "Avg. Cherenkov angle of clusters in pi/K/p regions (rad)", "angle"  },
    { "mip.angles_corrected", "Avg. Cherenkov angle of clusters in pi/K/p regions (rad) - (expected Cherenkov angle) + (HRS right arm central momentum corresponding cherenkov angle)", "angle_corrected"  },
    { "mip.n_chi2_phot",  "Clusters employed in  the chi2 test for pi/K/p hypothesis for this MIP","N_chi2_Photon"  },
    { "mip.chi2", " chi square values according to the hypothesis the MIP is a pi/k/p","chi2"},
    { "mip.chi2_prob", " chi square probability according to the hypothesis the MIP is a pi/k/p","chi2_prob"},
    { "mip.n_chi2_corrected_phot",  "Clusters employed in the $Corrected$ (that is after getting rid of the niose chi2 test for pi/K/p hypothesis for this MIP","N_chi2_corrected_Photon"  },
    { "mip.chi2_corrected", " chi square probability according to the hypothesis the MIP is a pi/k/p and after noise has been cut away","chi2_corrected"},
    { "mip.chi2_corrected_prob", " chi square probability according to the hypothesis the MIP is a pi/k/p","chi2_corrected_prob"},
   { "mip.noise_cut_success", "equal to 1 when the noise was succesfully thrown away","Getunresolved_noise_cut_success()"},
   { "mip.n_Maximum_Likelihood_phot", "Cluster employed in the Maximum Likelihood calculation according to the hypothesis the MIP is a pi/k/p", "N_MaximumLikelihood_Photon"}, 
   { "mip.Maximum_Likelihood", "Makimum Likelihood values according to the hypothesis the MIP is a pi/k/p", "MaximumLikelihood"}, 
    { "mip.rnphot",  "Resolved Clusters in pi/K/p region of this MIP","ResolvedN_Photon"  },
    { "mip.rangles", "Avg. Cherenkov angle of resolved clusters in pi/K/p regions (rad)", "Resolvedangle"  },
    { "mip.rangles_corrected", "Avg. Cherenkov angle of resolved clusters in pi/K/p regions (rad) - (expected Cherenkov angle) + (HRS right arm central momentum corresponding cherenkov angle)", "Resolvedangle_corrected"  },
   { "mip.rn_chi2_phot",  "Resolved clusters employed in  the chi2 test for pi/K/p hypothesis for this MIP","ResolvedN_chi2_Photon"  },
    { "mip.rchi2", " Resolved cluster chi square values according to the hypothesis the MIP is a pi/k/p","Resolvedchi2"},
    { "mip.rchi2_prob", " Resolved cluster chi square probability according to the hypothesis the MIP is a pi/k/p","Resolvedchi2_prob"},
    { "mip.rn_chi2_corrected_phot",  "Resolved clusters employed in the $Corrected$ (that is after getting rid of the niose chi2 test for pi/K/p hypothesis for this MIP","ResolvedN_chi2_corrected_Photon"  },
    { "mip.rchi2_corrected", " Resolved cluster chi square values according to the hypothesis the MIP is a pi/k/p and after noise has been cut away","Resolvedchi2_corrected"},
    { "mip.rchi2_corrected_prob", " Resolved cluster chi square probability according to the hypothesis the MIP is a pi/k/p","Resolvedchi2_corrected_prob"},
   { "mip.rnoise_cut_success", "equal to 1 when the noise was succesfully thrown away from resolved clusters","Getresolved_noise_cut_success()"},
   { "mip.rn_Maximum_Likelihood_phot", "Resoved cluster employed in the Maximum Likelihood calculation according to the hypothesis the MIP is a pi/k/p", "N_MaximumLikelihood_Photon"}, 
   { "mip.rMaximum_Likelihood", "Resolved Cluster Makimum Likelihood values according to the hypothesis the MIP is a pi/k/p", "ResolvedMaximumLikelihood"},  
    { 0 }
  };
  DefineVarsFromList( var6, mode, "fMIP." );

  return kOK;
}

//_____________________________________________________________________________
//
Double_t SBSGRINCH::Cherenkov_Angle(double mass, double momentum)  const
{
//
//   To calculate the Cherenkov angle for a particle of a certain momentum    
//   and mass 

  return acos(sqrt(1+(mass*mass)/(momentum*momentum))/n_radiator);

}


//_____________________________________________________________________________
Int_t SBSGRINCH::Decode( const THaEvData& evdata )
{
  //Decode RICH data and fill hit array

  if( !fIsInit ) return -255;
  if( !evdata.IsPhysicsTrigger() ) return -1;

  Clear();

  if( fDoBench ) fBench->Begin("Decode");

  for( UShort_t i = 0; i < fDetMap->GetSize(); i++ ) {
    THaDetMap::Module* d = fDetMap->GetModule( i );

    // FIXME: this is really dirty, but how can I decide, which module
    //        I am actually looking at ???
    // FIXME (joh): do we still need this?
//      if (d->slot == 25) {
//        if (evdata.GetNumChan( d->crate, d->slot ) > 0) {
//  	fTIRDat=evdata.GetData(d->crate,d->slot,0,0);
//        }
//      }
//      else {
      for( Int_t j = 0; j < evdata.GetNumChan( d->crate, d->slot ); j++) {
	Int_t chan = evdata.GetNextChan( d->crate, d->slot, j );
	if( chan > d->hi || chan < d->lo ) continue; // Not one of my channels

	// Get the data.
	Int_t nhit = evdata.GetNumHits( d->crate, d->slot, chan );
	if( GetNumHits()+nhit > fMaxNumHits ) {
	  //  Warning("Decode", "Too many hits! Should never ever happen! "
	  //	  "Event skipped.");
	  fHits->Clear();
	  if( fDoBench ) fBench->Stop("Decode");
	  return -2;
	}
	
	 for (int hit = 0; hit < nhit; hit++) {

             // FIXME:  Need to fix how geometry gets mapped to x/y
             static const Int_t chx[6]={0,5,1,4,2,3};

             Int_t richADC = 0;
             Int_t richXPAD_X_ADC=6;
             Int_t richYPAD_X_ADC=80;
             Int_t xhit = richADC*richXPAD_X_ADC+chx[chan%richXPAD_X_ADC];
             Int_t yhit = (chan/richXPAD_X_ADC)%richYPAD_X_ADC;

	  
	    // Fill hit array
	    UInt_t data = evdata.GetData(d->crate,d->slot,chan,hit);
	    UInt_t rawdata = evdata.GetRawData(d->crate,d->slot,chan,hit);
	    // FIX ME next line is for testing and can be removed later on
	    if ( (rawdata & 0xfff) != data ) {
	      cout<< "Something strange happened, "
		"raw data and data are not consistent" 
		  << endl;
	    }
	    if ( (rawdata & 0xc0000000) == 0x40000000) {
	      Double_t charge = static_cast<Double_t>(data);
	      Padn2xy( xhit, yhit, charge );
	      // cout << " ADC " << ADC << " chan " << chan;
	      // cout << " Hit (X,Y) " << xhit << " " << yhit;
	      // cout << "; charge " << charge << endl; 
	    }
	  }
      }
  }
  //  }

  if( fDoBench ) fBench->Stop("Decode");
  return GetNumHits();
}


//_____________________________________________________________________________
Int_t SBSGRINCH::CoarseProcess( TClonesArray& tracks )
{
  // Coarse processing of the RICH. We do nothing here since we
  // need precise tracking to be able to do anything useful.

  return 0;
}

//_____________________________________________________________________________
Int_t SBSGRINCH::FineProcess( TClonesArray& tracks )
{
   // The main RICH processing method. Here we
  //  
  // - Identify clusters
  // - Find the MIP(s) among the clusters
  // - Calculate photon angles
  // - Calculate particle probabilities

  // To put here otherwise Mip global variables are not zeroed in case 
  // of no cluster or no Mip were found 
  fMIP.Clear("F");

  if( !fIsInit ) return -255;

  if (fMaxNumHits<GetNumHits()) {  
    return -25;
  }  

  if( FindClusters() == 0 ) { 
    return -1;
  }

  if( FindMIP( tracks ) == 0 ) {
    return -2;
  }

  if( fDoResolve )
    ResolveClusters();
  //cout<<"stay and stop here!!!"<<endl;
  if( fDoBench ) fBench->Begin("FineProcess");

    Double_t central_momentum = 0.;
    Double_t central_e_angle     = 0.;
    Double_t central_pi_angle     = 0.;
    Double_t central_kaon_angle   = 0.;
    Double_t central_proton_angle = 0;
    // FC: presently not used, it has to be debugged first
    // Double_t central_electron_angle = 0;
    //
    Double_t coefficient = log(1. + 1./(sqrt(2.*cluster_distribution_sigma*
					    epsilon*TMath::Pi())));

  // Golden Track
  THaTrack* gold = NULL;

  Int_t igold = -1;
  THaSpectrometer *ec_sp = NULL;

  if( GetApparatus() != NULL ) {
    //    gold = static_cast<THaSpectrometer*>(fApparatus)->GetGoldenTrack();
    ec_sp = (THaSpectrometer *) GetApparatus();
    
    gold = ec_sp->GetGoldenTrack();
  
    // central momentum corresponding cherenkov angles 
    
    central_momentum = ec_sp->GetPcentral();
    central_e_angle     = Cherenkov_Angle(m_el,central_momentum);
    central_pi_angle     = Cherenkov_Angle(m_pi,central_momentum);
    central_kaon_angle   = Cherenkov_Angle(m_ka,central_momentum);
    central_proton_angle = Cherenkov_Angle(m_pr,central_momentum);
    //
    //  FC: the line below has to be declared somewhere else
    //   central_electron_angle = Cherenkov_Angle(m_el,central_momentum);
  }

  // For each track, find its angles of incidence with respect to
  // the RICH plane. We only need to do this once for each event.
  // Make a unit vector opposite to the track direction. The reason
  // for the - sign is that, with the definition of the RICH coordinate
  // system, fZax is opposite to the track, and we want theta=0 for 
  // normal incidence.
  // 2008/11/01 (FC): what stated above should not be anymore, now the 
  // RICH coordinate system corresponds to the track direction.
  //
  Int_t nTracks = tracks.GetLast()+1;
  struct TrackInfo_t {
    Double_t theta, phi;
  };
  TrackInfo_t* trkifo = new TrackInfo_t[nTracks];
  for( int i = 0; i<nTracks; i++ ) {
    THaTrack* trk = static_cast<THaTrack*>( tracks.At(i) );
    if( !trk ) continue;
    if( trk == gold )  igold = i;
    TVector3 tv( trk->GetTheta(), trk->GetPhi(), 1.0 );
    tv = tv.Unit();

    // MIP angles are in spherical coordinates in RICH system, units rad
    double cos_theta_mip = -fZax.Dot(tv);
    double theta_mip = acos(cos_theta_mip);
    double phi_mip = 0.0;
    if (cos_theta_mip != 1.) {
      double sin_theta_mip = sqrt(1.-cos_theta_mip*cos_theta_mip);
      phi_mip = acos( (-fXax.Dot(tv))/sin_theta_mip );
       if(-fYax.Dot(tv) < 0.) {
	 phi_mip = 2.*TMath::Pi() - phi_mip;
      }
    }

    trkifo[i].theta = theta_mip;
    trkifo[i].phi   = phi_mip;
  }

  // For each cluster, calculate Cherenkov angles of the cluster 
  // with respect to ALL of the MIPs. Store the average results with
  // each MIP
  Int_t nClusters = GetNumClusters();
  Int_t nResolvedClusters = GetNumResolvedClusters();
  for( int k = 0; k < (nClusters+nResolvedClusters); k++ ) {
    SBSGRINCH_Cluster* theCluster = 
      ( k < nClusters ) ? GetCluster(k) : GetResolvedCluster(k-nClusters); 
    if( !theCluster || theCluster->IsMIP() ) 
      continue;  // ignore MIPs

    // Find the MIP closest to this cluster. This is not necessarily
    // the MIP that generated the cluster.
    Double_t mindist = 1e10;
    Int_t imip = -1;
    SBSGRINCH_Cluster* theMIP_mindist = NULL;
    for( int i = 0; i<nTracks; i++ ) {
      if( !fMIPs[i] ) continue;
      Double_t dist = theCluster->Dist( fMIPs[i] );
      if( dist < mindist ) {
	mindist = dist;
	theMIP_mindist = fMIPs[i];
	imip = i;
      }
    }

    Double_t  x_photon = theCluster->GetXcenter();
    Double_t  y_photon = theCluster->GetYcenter();

    Double_t chi2_pi = 0;
    Double_t chi2_kaon = 0;
    Double_t chi2_proton = 0;
    Double_t MaximumLikelihood_pi = 0;
    Double_t MaximumLikelihood_kaon = 0;
    Double_t MaximumLikelihood_proton = 0;
    // FC: presently not used, it needs debug first
    //    Double_t MaximumLikelihood_electron = 0;

    for( int i = 0; i<nTracks; i++ ) {
      //    MIPtrack = static_cast<THaTrack*>( tracks.At(imip) );
      THaTrack* MIPtrack = static_cast<THaTrack*>( tracks.At(i) );
      SBSGRINCH_Cluster* theMIP = fMIPs[i];
      if( !MIPtrack || !theMIP ) continue;

      Double_t p_mip = MIPtrack->GetP();  // in GeV!
      if( p_mip == 0.0 ) continue;  // Track not fully reconstructed
      Double_t x_mip = theMIP->GetXcenter();
      Double_t y_mip = theMIP->GetYcenter();
      Double_t theta_mip = trkifo[i].theta;
      Double_t phi_mip   = trkifo[i].phi;

      // Calculate the avarage angle and  the maximum and minimum
      // angles allowed geometrically.

      Double_t theta_photon, phi_photon;
      Double_t  anglemax = RecoAng( x_photon, y_photon,theta_photon,
				    phi_photon, x_mip, y_mip, theta_mip,
				    phi_mip, -1 );
      Double_t  anglemin = RecoAng( x_photon, y_photon,theta_photon,
				    phi_photon, x_mip, y_mip, theta_mip,
				    phi_mip, 1 );
      Double_t  angle = RecoAng( x_photon, y_photon,theta_photon,
				 phi_photon, x_mip, y_mip, theta_mip,
				 phi_mip, 0 );

      //if the Mip analyzed is the closest to the cluster suppose
      //the cluster as generated by this MIP and set the cluster 
      //parameters accordingly. In any case this is for statistical 
      //purposes only. The analysis is performed through the Mip clusters
      //whose no-zero parameters are set afterwards)
      if (theMIP_mindist == theMIP) {
	theCluster->SetTrack( MIPtrack );
	theCluster->SetMIP( theMIP ); 
	theCluster->SetTheta_photon( theta_photon );
	theCluster->SetPhi_photon( phi_photon );
	theCluster->SetAngle( angle );
      }

      // Calculate the three Cherenkov angles expected in case the MIP 
      // was generated by a pion, a kaon and a proton respectively.

      Int_t ResolvedFlag = 0;
      if(k >= nClusters) ResolvedFlag = 1;

      Double_t e_angle      = Cherenkov_Angle(m_el,p_mip);
      Double_t pi_angle     = Cherenkov_Angle(m_pi,p_mip);
      Double_t kaon_angle   = Cherenkov_Angle(m_ka,p_mip);
      Double_t proton_angle = Cherenkov_Angle(m_pr,p_mip);

      if(anglemin > (angle - fiducial_zone_range)) 
	anglemin = angle - fiducial_zone_range;

      if(anglemax < (angle + fiducial_zone_range)) 
	anglemax = angle + fiducial_zone_range;


      // if the cluster is inside one of the three fiducial regions 
      // of the MIP analyzed (one region for each kind of particle)
      // consider it in the calculations for the angle average (and hence 
      // of the momentum of the particle) of the MIP itself.

      if ( (anglemin < pi_angle) && (anglemax > pi_angle) )  {
	theMIP->Insert_Photon(0, angle, ResolvedFlag, pi_angle, 
			      central_pi_angle);
      }
      if ( (anglemin < kaon_angle) && (anglemax > kaon_angle) ) {
	theMIP->Insert_Photon(1, angle, ResolvedFlag, kaon_angle,
			      central_kaon_angle);
      }
      if ( (anglemin < proton_angle) && (anglemax > proton_angle) ){
	theMIP->Insert_Photon(2, angle, ResolvedFlag, proton_angle,
			      central_proton_angle);
      }
      // FC: presently not used, it has to be debugged first
      //
      //   if ( (anglemin < electron_angle) && (anglemax > electron_angle) ){
      //	theMIP->Insert_Photon(2, angle, ResolvedFlag, electron_angle,
      //			      central_electron_angle);
      // }
      if ( ((anglemin < pi_angle) && (anglemax > pi_angle)) 
	   || ( (anglemin < kaon_angle) && (anglemax > kaon_angle) )
	   ||( (anglemin < proton_angle) && (anglemax > proton_angle) )
	   // FC: presently not used, has to be debugged first
	   //  ||( (anglemin < electron_angle) && (anglemax > electron_angle) )
	   ) {
 
	chi2_pi = ((angle - pi_angle)*(angle - pi_angle))/ 
	  ((cluster_distribution_sigma)*(cluster_distribution_sigma));
	chi2_kaon = ((angle - kaon_angle)*(angle - kaon_angle))/ 
	  ((cluster_distribution_sigma)*(cluster_distribution_sigma));
	chi2_proton = ((angle - proton_angle)*(angle - proton_angle))/ 
	  ((cluster_distribution_sigma)*(cluster_distribution_sigma));
	MaximumLikelihood_pi = 
	  coefficient*log(exp(-((angle - pi_angle)*(angle - pi_angle)) 
			      /((cluster_distribution_sigma)*
				(cluster_distribution_sigma)*2))
			  *coefficient+1);
	MaximumLikelihood_kaon = 
	  coefficient*log(exp(-((angle - kaon_angle)*(angle - kaon_angle)) 
			      /((cluster_distribution_sigma)*
				(cluster_distribution_sigma)*2))
			  *coefficient+1);
	MaximumLikelihood_proton = 
	  coefficient*log(exp(-((angle - proton_angle)*(angle - proton_angle)) 
			      /((cluster_distribution_sigma)*
				(cluster_distribution_sigma)*2))
			  *coefficient+1);
	//if the Mip analyzed is the closest to the cluster suppose
	//the cluster as generated by this MIP and set the single cluster 
	//chi square value accordingly. In case more Mips (and tracks)
	//should be considered, all cluster parameters should be put
	//in an array
	if (theMIP_mindist == theMIP) {
	  theCluster->SetPionChi2AnalysisFlag(kTRUE);
	  theCluster->SetKaonChi2AnalysisFlag(kTRUE);
	  theCluster->SetProtonChi2AnalysisFlag(kTRUE);
	  theCluster->Insert_chi2(0, chi2_pi);
	  theCluster->Insert_chi2(1, chi2_kaon);
	  theCluster->Insert_chi2(2, chi2_proton);
	}
	theMIP->Insert_chi2(0, chi2_pi, ResolvedFlag);
	theMIP->Insert_chi2(1, chi2_kaon, ResolvedFlag);
	theMIP->Insert_chi2(2, chi2_proton, ResolvedFlag);
	theMIP->Insert_MaximumLikelihood(0, MaximumLikelihood_pi, 
					 ResolvedFlag);
	theMIP->Insert_MaximumLikelihood(1, MaximumLikelihood_kaon, 
					 ResolvedFlag);
	theMIP->Insert_MaximumLikelihood(2, MaximumLikelihood_proton, 
					 ResolvedFlag);
      }
    }
  }
  
  // Save the results for the Golden Track MIP
  if( igold >= 0 && fMIPs[igold] ) 
    {
      fMIP = *fMIPs[igold];
      fMIP.SetTheta_photon( trkifo[igold].theta );
      fMIP.SetPhi_photon( trkifo[igold].phi );
      for( Int_t ResolvedFlag = 0; ResolvedFlag < 2; ResolvedFlag++ ) 
	{
	  fMIP.Setchi2_prob(0, fMIP.Getchi2_pi(ResolvedFlag),
			    fMIP.GetN_chi2_phot_pi(ResolvedFlag),
			   ResolvedFlag);
	  fMIP.Setchi2_prob(1, fMIP.Getchi2_k(ResolvedFlag),
			    fMIP.GetN_chi2_phot_k(ResolvedFlag),
			   ResolvedFlag);
	  fMIP.Setchi2_prob(2, fMIP.Getchi2_p(ResolvedFlag),
			    fMIP.GetN_chi2_phot_p(ResolvedFlag),
			   ResolvedFlag);
	}
      for( Int_t ResolvedFlag = 0; ResolvedFlag < 2; ResolvedFlag++ ) 
	{
	  fMIP.Insert_N_chi2_corrected_Photon(0, 
					 fMIP.GetN_chi2_phot_pi(ResolvedFlag),
					      ResolvedFlag);
	  fMIP.Insert_N_chi2_corrected_Photon(1, 
					 fMIP.GetN_chi2_phot_k(ResolvedFlag),
					      ResolvedFlag);
	  fMIP.Insert_N_chi2_corrected_Photon(2, 
                                         fMIP.GetN_chi2_phot_p(ResolvedFlag),
					      ResolvedFlag);
	  fMIP.Insert_chi2_corrected(0, fMIP.Getchi2_pi(ResolvedFlag),
				     ResolvedFlag);
	  fMIP.Insert_chi2_corrected(1, fMIP.Getchi2_k(ResolvedFlag), 
				     ResolvedFlag);
	  fMIP.Insert_chi2_corrected(2, fMIP.Getchi2_p(ResolvedFlag), 
				     ResolvedFlag);
	  fMIP.Setchi2_corrected_prob(0, fMIP.Getchi2_pi(ResolvedFlag),
				      fMIP.GetN_chi2_phot_pi(ResolvedFlag),
			   ResolvedFlag);
	  fMIP.Setchi2_corrected_prob(1, fMIP.Getchi2_k(ResolvedFlag),
				      fMIP.GetN_chi2_phot_k(ResolvedFlag),
			   ResolvedFlag);
	  fMIP.Setchi2_corrected_prob(2, fMIP.Getchi2_p(ResolvedFlag),
				      fMIP.GetN_chi2_phot_p(ResolvedFlag),
			   ResolvedFlag);
	  Int_t flag = 1;
	  Int_t BadFlag = 0;
	  Int_t TrialNumber = 0;
	  while(flag) {
	    TrialNumber++;
	    if ( fMIP.GetN_chi2_corrected_phot_pi(ResolvedFlag) < 1 )
	      {
		// Number of photons in the analysis null, skip away
		// note: at the moment number of photons equal for each 
		// hypothesis  on the kind of particle crossing the Rich  
		flag = 0;
	      }
	    else  
	      {
		if( (TMath::Prob(fMIP.Getchi2_corrected_pi(ResolvedFlag),
				 fMIP.GetN_chi2_corrected_phot_pi(ResolvedFlag))
		     < acceptable_chi2_prob) 
		    && 
		    (TMath::Prob(fMIP.Getchi2_corrected_k(ResolvedFlag),
				 fMIP.GetN_chi2_corrected_phot_k(ResolvedFlag)) 
		     < acceptable_chi2_prob) 
		    && 
		    (TMath::Prob(fMIP.Getchi2_corrected_p(ResolvedFlag),
				 fMIP.GetN_chi2_corrected_phot_pi(ResolvedFlag)) 
		     < acceptable_chi2_prob))
		  {
		    // no chi square value acceptable; try to cut away noise
		    if ( (fMIP.GetN_chi2_corrected_phot_pi(ResolvedFlag) <
			  minimum_chi2_degree_of_freedom) && 
			 (fMIP.GetN_chi2_corrected_phot_k(ResolvedFlag) < 
			  minimum_chi2_degree_of_freedom) && 
			 (fMIP.GetN_chi2_corrected_phot_p(ResolvedFlag) < 
			  minimum_chi2_degree_of_freedom) ) {
		      // too few clusters remained for the analysis; 
		      // give up
		      flag = 0;
		    } else {
		      // try to get rid of the noise 
		      if (TrialNumber > clear_noise_trial_maximum_number) {
			flag = 0; 
		      } else {
			Int_t success = ClearNoise(igold, ResolvedFlag);
			if ( success <= -2) flag = 0;
			// no succes or ambigous solution
			if ( success == -1) BadFlag = 1;
			if ((success == 0) && (BadFlag == 1) ) flag = 0;
			// ambigous solution
		      }
		    }
		  } else {
		  // everithing is o.k.
		  fMIP.Setnoise_cut_success(1, ResolvedFlag);
		  flag = 0;
		}
	      }
	  }
	}
    }
  delete [] trkifo;
  
  if( fDoBench ) fBench->Stop("FineProcess");
  return 0;
}


//_____________________________________________________________________________
//***************************************************************************
// int ReadData()
//    ... read the certain number of data set of pad-number and charge
//   for ONE event.
// 
//**************************************************************************
// THIS FUNCTION IS FOR TEST PURPOSES ONLY.  FOR REAL DATA, USE Decode()
//
Int_t SBSGRINCH::ReadData( FILE* infile )
{
  //Read a single event from input file 'infile' into internal hit array.
  //Clear all previous hit information.

  char line[100];
  int count = 0;
  Clear();

  fgets(line,100,infile);
#ifdef HAS_SSTREAM
  istringstream sline(line);
#else
  istrstream sline(line,100);
#endif
  sline >> count;
  if( sline.good() ) {
    if( fDebug > 0 ) printf( "*****Hit count: %5d\n", count );
  } else {
    return 0;
  }
  
  if( count >= fNelem ) {
    Warning("ReadData", "Too many hits! Event skipped.");
    Clear();
    return 0;
  }

  for( int i=0; i<count; i++ ) {
    int i_pad, j_pad; 
    double charg;
    fgets(line,100,infile);
    sline.clear(); sline.seekg(0);
    sline >> i_pad >> j_pad >> charg;
    if( !sline.good() ) break;
    Padn2xy( i_pad, j_pad, charg );
    if( fDebug > 0 ) printf("%5d %5d %5f\n",i_pad,j_pad,charg);
  }

  return GetNumHits();
}

//__________________________________________________________________________
//**************************************************************************
// void Padn2xy( int i_pad, int j_pad, double charg )
//   Define the ith fHits (X and Y coordinate, I and J (serial number  
//   in X and in Y raw respectively, charge, serial number from the pad  
//   number (pad) and the charge collected (charg).
//      
//   i_pad and j_pad start counting at 0.
//
//**************************************************************************
void SBSGRINCH::Padn2xy( Int_t i_pad, Int_t j_pad, Double_t charg ) 
{
  //Add a hit to the internal hit array, given i and j indices.

  Int_t i = GetNumHits();
  Int_t icharg = int(charg + 0.5);
  Float_t i_X = i_pad*PAD_SIZE_X + PAD_SIZE_X/2;
  for( int i=0; i<fNseg; i++ ) 
    {
      int k = 3*i;
      if(i_X >= fXseg[k] -  fXseg[k+2] && i_X < fXseg[k+1] - fXseg[k+2]) 
	{
	  // apply the offset relative to the plane and quit
	  i_X +=  fXseg[k+2];
	  break;
	}
    }
  SBSGRINCH_Hit* h = new( (*fHits)[i] ) 
    SBSGRINCH_Hit( i, icharg, i_pad, j_pad, 
		 i_X + PAD_SIZE_X/2,
		 j_pad*PAD_SIZE_Y + PAD_SIZE_Y/2 );
  if( !h ) {
    // urgh
    Error("Padn2xy", "Failed to allocate new hit, "
	  "index,ipad,jpad = %d %d %d ", i, i_pad, j_pad );
    return;
  }
  h->SetFlag( 0 );
}  


//__________________________________________________________________________
void SBSGRINCH::DeleteClusters()
{
  //Delete all clusters

  fClusters->Clear("C");
  fResolvedClusters->Clear("C");
  return;
}

//__________________________________________________________________________
Int_t SBSGRINCH::FindClusters()
{
  // Group the hits that are currently in the array into clusters.
  // Return number of clusters found.
 
  // minimum distance between two pads.
  const double par1 = 2.0; 

  // maximum distance in X between two fired pads to be in the same cluster.
  const double par2 = PAD_SIZE_X+0.1;  

  // maximum distance in Y between two fired pads to be in the same cluster.
  const double par3 = PAD_SIZE_Y+0.1;  

  DeleteClusters();

  if( fDoBench ) fBench->Begin("FindClusters");

  SBSGRINCH_Hit* theHit;
  SBSGRINCH_Cluster* theCluster;
  Int_t nHits  = GetNumHits();
  Int_t nClust = 0;

  for( int k=0; k<nHits; k++ ) {
    if( !(theHit = GetHit(k))) continue;

    // HitFlag not equal 0: The Hit was alredy processed
    // HitFlag equal 0: insert the Hit as first element of the cluster

    if( theHit->GetFlag() == 0 ) {
  
      theCluster = new( (*fClusters)[nClust++] ) SBSGRINCH_Cluster();
      theHit->SetFlag(1);
      theCluster->Insert( theHit );
      
      //---Scanning all the hit pads

      int flag = 0;
      while( flag==0 ) {
	flag = 1;
	SBSGRINCH_Hit* hit;
	for( int i=0; i<nHits; i++ ) {
	  if( !(hit = GetHit(i))) continue;
	  if( hit->GetFlag() == 0 && 
	      theCluster->Test( hit, par1, par2, par3 )) {
	      
	    // the hit belongs to the Cluster

	    theCluster->Insert( hit );

	    //tagging to avoid double processing
	    hit->SetFlag(1);  

	    flag = 0;
	    // at least one new element in the Cluster; check again.
	    // FIXME: this loop is O(2*nHits^2) - 2e6 iterations for 1k hits!!
	    // ->extremely inefficient
	  }
	}
      }
    }
  }
  // now find if there are clusters genetrated by overlapping of clusters
  // (number of overlapping clusters = number of local maximums in a cluster).
  if( fDoResolve ) {
    for( int k = 0; k < nClust; k++ ) {
      if( !(theCluster = GetCluster(k))) continue;
      theCluster->FindLocalMaximumNumber();
    }
  }
  if( fDoBench ) fBench->Stop("FindClusters");
  return nClust;
}

//__________________________________________________________________________
Int_t SBSGRINCH::FindMIP( const TClonesArray& tracks )
{

  // Find the MIPs among the clusters. A MIP is a cluster generated 
  // by a track rather than by Cherenkov photons. Every track 
  // crossing the RICH plane is associated with a MIP. 
  
  // FIXME: if no tracks match well, assign MIP to the largest
  // charge cluster
  // FIXME: save some indicator how well track matches MIP 
  
  if( fDoBench ) fBench->Begin("FindMIP");
  
  // Clear the list of MIPs - we'll be making a new one
  delete [] fMIPs; fMIPs = 0;
  fMIP.Clear("F");
  fTrackX = fTrackY = kBig;
  
  Int_t ntracks = tracks.GetLast()+1;
  Int_t nClust  = GetNumClusters();
  // Quit if no tracks or no clusters. Nothing to do.
  if( ntracks == 0 || nClust == 0 ) return 0;
  
  fMIPs  = new SBSGRINCH_Cluster* [ ntracks ];
  for( int i=0; i<ntracks; i++ ) fMIPs[i] = 0;
  Int_t nfound = 0;
  
  // Find the closest Cluster for each track.
  
  for( Int_t k=0; k<ntracks; k++ ) {
    THaTrack* theTrack = static_cast<THaTrack*>( tracks.At(k) );
    if( theTrack == NULL ) continue;
    // Golden Track?
    THaTrack* tr = NULL;
    THaSpectrometer *ec_sp = NULL;
    if( GetApparatus() != NULL)
      ec_sp = (THaSpectrometer *) GetApparatus(); 
    
    tr = ec_sp->GetGoldenTrack();
    bool is_golden = ( theTrack == tr );
    
    // Find the coordinates of the track crossing point in the RICH plane.
    // trackX/Y: Coordinates of the crossing point in the RICH system.
    Double_t pathlength, trackX, trackY;
    CalcTrackIntercept( theTrack, pathlength, trackX, trackY );
    
    // RICH coordinates are in mm
    trackX *= -1000.;
    trackY *= -1000.;
    
    if( is_golden || ntracks == 1 ) {
      fTrackX = trackX;
      fTrackY = trackY;
    }

    SBSGRINCH_Cluster* theCluster;
    SBSGRINCH_Cluster* minDistClust   = NULL;
    SBSGRINCH_Cluster* maxChargeClust = NULL;
    SBSGRINCH_Cluster* theMIP         = NULL;
    
    // cout << " -trackX " << -trackX << " -trackY " << -trackY;
    //cout  << " fXseg[0] " << fXseg[0] << " fXseg[13] " << fXseg[13] << endl;
    //cout << " PAD_SIZE_Y*fNypads " << PAD_SIZE_Y*fNypads << endl; 
    
    if(( trackX  < fXseg[0] || trackX > fXseg[3*fNseg-2] ) ||  //EC: invece di 13 metterei comunque parametrico: fNseg-2 
       //
       // DEBUGGING
       //
       // ( trackX < 1050. ) ||
       ( trackY < 0 || trackY > PAD_SIZE_Y*fNypads) ) continue;
    // skip if the track does not hit the RICH pad plane
    
    //cout << " (fMIP_through_interception " << fMIP_through_interception;
    //cout << endl;
    
    if (fMIP_through_interception == 3)
      {
	// The MIP is forced to be the interception between the 
	// track ant the pad plane 
	// create an additional fictious cluster that is supposed to be a 
	// Mip and whose coordinates are equal to the interception of the 
	// track with the pad plane
	// NOTE: This cluster must be handled carefully. No physical hit and
	// a "virtual" charge of 10000 associated with hit. 
	theCluster = new( (*fClusters)[nClust++] ) SBSGRINCH_Cluster();
	theCluster->SetXcenter(trackX);
	theCluster->SetYcenter(trackY);
	theCluster->SetCharge(10000);
	theMIP = theCluster;
	nfound++;
	theMIP->MakeMIP();
	theMIP->SetFictious_MIP_Flag(fMIP_through_interception);
      } else {
	// try to find a Mip among the clusters.
	// First of all take in account the shift in XMIP position due to the 
	// mis-positioning of the RICH and the presence of not sensible 
	// material between planes.
	bool found = false;
	for( int i=0; i<fNseg; i++ ) {
	   int kk = 3*i;  //EC: questo potrebbe essere il problema! c'e' gia' un loop su k; questa variabile va cambiata -> metto kk
	    if( trackX >= fXseg[kk] && trackX < fXseg[kk+1] ) {
	    // If within x-region of the plane, then apply offset and quit
	    // trackX -= fXseg[k+2]; 
	    // now X coordinate of the pads are corrected for the offset 
	    // the istruction above is hence not more necessary
	    // cout << " trackX " << trackX <<  " fXseg[k] " << fXseg[k];
	    // cout << " fXseg[k] - Value_to_subtract + fXseg[k+2] ";
	    // cout << fXseg[k] - Value_to_subtract + fXseg[k+2];
	    // cout << " fXseg[k+1] " <<  fXseg[k+1] << endl;
	    found = true;
	    break;
	      }
	}
	// this track outside of any defined plane? -> ignore the track
	// (if theMip search algortihm allow to do it
	if( ( !found && ( fMIP_through_interception == 0) ) && fNseg>0 )
	  continue;
	
	Double_t         maxcharge      = 0.0;
	Double_t         mindist        = 1e10;
	Int_t            ncloseClusters = 0;
	
	for( int j=0; j<nClust; j++ ) {
	  if( !(theCluster = GetCluster(j))) continue;
	  // For any given track, find the cluster(s) within the
	  // search radius fMaxdist around the point where the
	  // track crosses the RICH plane.
	  
	  Double_t dx    = theCluster->GetXcenter() - trackX;
	  Double_t dy    = theCluster->GetYcenter() - trackY;
	  Double_t dist2 = dx*dx + dy*dy;
	  
	  if( dist2 < fMaxdist2 ) {
	    ncloseClusters++;
	    Double_t charge  = theCluster->GetCharge();
	    if( charge > maxcharge ) {
	      maxChargeClust = theCluster;
	      maxcharge      = charge;
	    }
	    if( dist2 < mindist ) {
	      minDistClust = theCluster;
	      mindist      = dist2;
	    }	  
	  }
	}//end loop over clusters
	// Simple logic for now: 
	// If any clusters are in the search area, then the MIP is the 
	// one with the maximum charge
	if( ncloseClusters > 0 )
	  theMIP = maxChargeClust;
	
	if( theMIP ) {
	  nfound++;
	  theMIP->MakeMIP();
	} else {
	  // if the MIP was supposed hit the pad plane,
	  // act according to the value of  fMIP_through_interception
	  if( (fMIP_through_interception == 2) || 
	      ( (fMIP_through_interception == 1) && !found ) ) {
	    // create an additional fictious cluster that is supposed to be a 
	    // Mip and whose coordinates are equal to the interception of the 
	    // track with the pad plane
	    // NOTE: This cluster must be handled carefully. No physical hit and
	    // a virtual charge of 10000 associated with hit. 
	    theCluster = new( (*fClusters)[nClust++] ) SBSGRINCH_Cluster();
	    theCluster->SetXcenter(trackX);
	    theCluster->SetYcenter(trackY);
	    theCluster->SetCharge(10000);
	    theMIP = theCluster;
	    nfound++;
	    theMIP->MakeMIP();
	    if(!found)
	      {
		theMIP->SetFictious_MIP_Flag(1); // the Mip has hit a not 
		// sensible part of the RICH
	      } 
	    else
	      {
		theMIP->SetFictious_MIP_Flag(2); // the Mip has hit a  
		// sensible part of the RICH
	      } 
	  }
	}
      }
    // Save the MIP found for this track, if any
    fMIPs[k] = theMIP;
    
  } //end loop over tracks  
  
  if( fDoBench ) fBench->Stop("FindMIP");
  return nfound;
}
 
//__________________________________________________________________________
Int_t SBSGRINCH::ResolveClusters()
{
  // Resolve the clusters. To be calle only after FindClusters.
  // Return the new number of clusters.

  if( fDoBench ) fBench->Begin("ResolveClusters");

  Int_t nClust = GetNumClusters();

  fResolvedHits->Clear("C");
  fResolvedClusters->Clear("C");
  Int_t nResolvedClusters = 0;

  for( int k = 0; k < nClust; k++ ) {
    SBSGRINCH_Cluster* theCluster = GetCluster(k);
    if( !theCluster ) continue;
    Int_t Number = theCluster->GetLocalMaximumNumber();
    if( (Number > 1) && !(theCluster->IsMIP()) ) {
      TIter next(theCluster->GetHitList());
      // Make a resolved cluster for each local maximum of the cluster
      while( const SBSGRINCH_Hit* pHit = static_cast<SBSGRINCH_Hit*>(next())) {
	if( pHit->GetVeto() != 0 ) continue;
	// found a local maximum
	SBSGRINCH_Cluster* the_Resolved_Cluster = 
	  new( (*fResolvedClusters)[nResolvedClusters++] ) SBSGRINCH_Cluster();

	// the real work is done here -
	theCluster->FindResolvedClusterElements
	  ( pHit, the_Resolved_Cluster, fResolvedHits );
      }
    }
    else {
      new( (*fResolvedClusters)[nResolvedClusters++] ) 
	SBSGRINCH_Cluster( *theCluster );
    }
  }

  if( fDoBench ) fBench->Stop("ResolveClusters");
  return nResolvedClusters;

}

//_____________________________________________________________________________
//**************************************************************************
// The function for calculating "phi_photon";the azimuthal angle of 
// the C-photon.
//**************************************************************************
Double_t SBSGRINCH::Get_phi_photon( Double_t x_photon, Double_t y_photon,
				 Double_t x_mip, Double_t y_mip,
				 Double_t theta_mip, Double_t phi_mip, 
				 Int_t Calculation_kind) const
{

  // Calculation kind: 
  //    0:  medium angle ("normal" calculation with the photon supposed 
  //                      generated in the middle of the radiator).
  //   -1: maximum possible angle (photon generated with minimum energy at the
  //                               end of the radiator)
  //    1: minimum possible angle (photon generated with the maximum energy 
  //                              at the beginning of the radiator)

  if( fDoBench ) fBench->Begin("Get_phi_photon");

  //calculation for "phi_photon"
  double sint = sin(theta_mip);
  double cost = cos(theta_mip);
  double tant = sint/cost;

  Double_t answer = 0.0;
  double length = (L_RAD+l_quartz+l_gap)*tant;
  double nr2 = n_radiator*n_radiator;
  double nq2 = n_quartz*n_quartz;
  double ng2 = n_gap*n_gap;
  double sq1=nq2-nr2*sint*sint;
  double sq2=ng2-nr2*sint*sint;
  double term1 = L_RAD*tant +
    l_quartz*n_radiator*sint/sqrt(sq1) +
    l_gap*n_radiator*sint/sqrt(sq2);
  double term2  = sin(phi_mip);
  double term3  = cos(phi_mip);

  double x_origin = x_mip - length*term3;  //the entrance point of the MIP
  double y_origin = y_mip - length*term2;

 double l_used = l_emission;
  if (Calculation_kind == -1) 
    {  
      l_used = L_RAD;
    }
  if (Calculation_kind == 1) 
    {  
      l_used = 0.;
    }

  term1 = l_used*tant;

  double ypart = y_photon - y_origin - term1*term2;
  double xpart = x_photon - x_origin - term1*term3;
  double hyp   = sqrt( xpart*xpart + ypart*ypart );
  if(fabs(xpart) > hyp+0.0001) {
    // it should happen only for ypart << xpart because of calculation rounding
    // double ratio = xpart/hyp;
    // cout << "WARNING xpart/hyp = " << ratio << 
    // ", (its module should be < 1)" << endl;
    answer = 0.0;
  } else {
    answer = acos(xpart/hyp);
  }
  if(ypart<0 && xpart>0) answer = - answer + 2.*M_PI;
  if(ypart<0 && xpart<0) answer = - answer + 2.*M_PI;
  //***These 2 lines above are to get around the problem with "arccos".
  //If these are missing, phi_photon calculation will return wrong result. 

  if( fDebug > 2 ) {
    printf("********************************************\n");
    printf("* CHERENKOV ANGLE RECONSTRUCTION           *\n");
    printf("********************************************\n");
    printf("***phi_photon***\n");
    printf("   phi_photon = atan((%f - %f*%f*%f) \n",
	   y_photon,l_used,term1,term2);
    printf("           / (%f - %f*%f*%f)\n",x_photon,l_used,term1,term3); 
    printf("   phi_photon = %f DEG \n",answer/M_PI*180);
  }

  if( fDoBench ) fBench->Stop("Get_phi_photon");
  return answer;
}
//___________________________________________________________________________

//***************************************************************************
//the function for computing "dist_a";the distance between the MIP impact pad
//and the origin.
//***************************************************************************
Double_t SBSGRINCH::Get_a( Double_t theta_mip, Int_t Calculation_kind) const
{

  // Calculation kind: 
  //    0:  medium angle ("normal" calculation with the photon supposed 
  //                      generated in the middle of the radiator).
  //   -1: maximum possible angle (photon generated with minimum energy 
  //                               at the end of the radiator)
  //    1: minimum possible angle (photon generated with the maximum energy 
  //                               at the beginning of the radiator)

  //  double nr2 = n_radiator*n_radiator;
  //  double nq2 = n_quartz*n_quartz;
  //  double ng2 = n_gap*n_gap;
  //  double sq1=nq2-nr2*sin(theta_mip)*sin(theta_mip);
  //  double sq2=ng2-nr2*sin(theta_mip)*sin(theta_mip);

  double l_used = l_emission;
  if (Calculation_kind == -1) 
    {  
      l_used = L_RAD;
    }
  if (Calculation_kind == 1) 
    {  
      l_used = 0.;
    }

  //  double answer1 = (L_RAD - l_used)*tan(theta_mip) +
  //    l_quartz*n_radiator*sin(theta_mip)/sqrt(sq1) +
  //    l_gap*n_radiator*sin(theta_mip)/sqrt(sq2);

  double answer = (L_RAD - l_used + l_quartz + l_gap) * tan(theta_mip);
  
  //cout << " nominal distance " << answer1 << " effective distance ";
  //cout << answer << " theta_mip " << theta_mip <<  endl;

  //double theta_effective = atan(answer/(L_RAD - l_used + l_quartz + l_gap));
  // double theta_inclination = theta_mip - theta_effective;
 
  // cout << " theta effective " << theta_effective; 
  // cout << " rich inclination angle " << theta_inclination << endl; 

  if( fDebug > 2 ) {
    printf("***dist_a : MIP-origin***\n");
    printf("   dist_a = (%f - %f + %f + %f)*tan(%f)\n"
	   ,L_RAD,l_emission,l_quartz,l_gap,theta_mip/M_PI*180);
    printf("   dist_a = %f mm\n",answer);
  }

  return answer;
}

//___________________________________________________________________________
//***************************************************************************
//The function for computing "dist_b";distance between the photon pad and 
//the origin on the pad plane that is refered as "b" in [1] and
//for computing "dist_R";the distance between the photon pad 
//and the MIP impact pad.

//This function calls : Get_a
//                      Get_phi_photon
//***************************************************************************
Double_t SBSGRINCH::Get_b( Double_t x_photon, Double_t y_photon,
			Double_t x_mip, Double_t y_mip,
			Double_t theta_mip, Double_t phi_mip, 
			Int_t Calculation_kind) const
{

  // Calculation kind: 
  //    0:  medium angle ("normal" calculation with the photon supposed 
  //                      generated in the middle of the radiator).
  //   -1: maximum possible angle (photon generated with minimum energy 
  //                               at the end of the radiator)
  //    1: minimum possible angle (photon generated with the maximum energy 
  //                               at the beginning of the radiator)

  double dist_a = Get_a(theta_mip, Calculation_kind);

  double term1  = x_photon - x_mip;
  double term2  = y_photon - y_mip;
  double dist_R = 0.0;

  double phi_photon = Get_phi_photon
    (x_photon,y_photon,x_mip,y_mip,theta_mip,phi_mip,Calculation_kind);

  if( fDebug > 0 ) {
    dist_R = sqrt( term1*term1 + term2*term2 );
    printf("dist_a = %f , dist_R = %f , phi_photon = %f\n",dist_a,
	   dist_R, phi_photon);
  }
  //Computing "dist_b(answer)"

  double term3 = 0.0;
  double term4 = 0.0; 
  //cout << "phi_photon " << phi_photon << endl;
    
  double cosp = cos(phi_photon);
  double sinp = sin(phi_photon);
  if(fabs(cosp) > 0.0003) 
    term3 = (term1 + dist_a*cos(phi_mip))/cosp;
  if(fabs(sinp) > 0.0003) {
    term4 = (term2 + dist_a*sin(phi_mip))/sinp;
    if(term3 == 0) {
      term3 = term4; // phi = 90 degree.
    }
  } else
    term4 = term3; // phi = 0 degree.



  Double_t answer = (term3 + term4)/2.0;

  // average between the two (in principle identical) values.

  if( fDebug > 0 ) {
    printf("term3 term4 ");
    printf("%4f,%4f",term3,term4);
    if( fDebug > 2 ) {
      printf("***dist_R: photon-MIP ***\n");
      printf("dist_R = sqrt(%f^2 + %f^2)\n",term1,term2);
      printf("dist_R = %f mm\n",answer);
      printf("***dist_b: photon-origin***\n");
      printf("   dist_b = %f mm\n",answer);
    }
  }

  return answer;
}

//___________________________________________________________________________
//***************************************************************************
//The function for computing the "theta_photon"; the polar angle of the photon
//with respect NOT to the MIP, but to the z-axis that is parpendicular to the
//pad plane.
//This function calls :  Get_b()
//***************************************************************************
Double_t SBSGRINCH::Get_theta_photon(Double_t x_photon, Double_t y_photon,
				  Double_t x_mip, Double_t y_mip,
				  Double_t theta_mip, Double_t phi_mip, 
				  Int_t Calculation_kind)  const
{
  // Calculation kind: 
  //    0:  medium angle ("normal" calculation with the photon supposed 
  //                      generated in the middle of the radiator).
  //   -1: maximum possible angle (photon generated with minimum energy 
  //                              at the end of the radiator)
  //    1: minimum possible angle (photon generated with the maximum energy 
  //                              at the beginning of the radiator)

  static const double twopi = 2.0*TMath::Pi();

  if( fDoBench ) fBench->Begin("Get_theta_photon");

  double TOLERANCE=0.0001;  //tolerance in the precision of calculating 
  double term1=0.0,term2=0.0,term3=0.0;  //for clarity of calculation
  double sq1,sq2;  //terms in square root
  double nr2 = n_radiator*n_radiator;
  double nq2 = n_quartz*n_quartz;
  double ng2 = n_gap*n_gap;

  double l_used = l_emission;
  if (Calculation_kind == -1) 
    {  
      l_used = L_RAD;
      nq2 = n_quartz_min*n_quartz_min;
      nr2 = n_radiator_min*n_radiator_min;
    }
  if (Calculation_kind == 1) 
    {  
      l_used = 0.;
      nq2 = n_quartz_max*n_quartz_max;
      nr2 = n_radiator_max*n_radiator_max;
    }


  //INITIALIZING
  //used in computing the "theta_photon" as a candidate for 
  //the "dist_b"
  double trial = 0.0;

  //arbitrary
  double theta_photon = 45./180.*M_PI;   
  
  //upper and lower value of BRACKETING(calculating)
  double upper = M_PI/2.;
  double lower = 0.0;

  double dist_b = Get_b(x_photon, y_photon, x_mip, y_mip, theta_mip, phi_mip, 
		       Calculation_kind);
  
  //Loop for computing "theta_photon" solving the eq. for the "dist_b"

  // due to precision problems in solving the equation for theta_photon 
  // (see below), increase TOLERANCE when dist_b is big otherwise the program 
  //will loop without an end !!

  //cout << dist_b << endl;
  if(dist_b > 300.) 
    {TOLERANCE = 0.0002;}
  if(dist_b > 400.)
    { TOLERANCE = 0.0003;}
  if(dist_b > 500.)
    {TOLERANCE = 0.0004;}
  if(dist_b > 550.)
    {TOLERANCE = 0.0005;}
  if(dist_b > 600.)
    {TOLERANCE = 0.0007;}
  int IterationNumber = 0;

  /*************************************************************************
  // GUIDO part:

  while( fabs(dist_b-trial)>TOLERANCE ) {
    IterationNumber ++;
    //pre-calculation
    double sint = sin(theta_photon);
    double sint2 = sint*sint;
    sq1=nq2-nr2*sint2;
    sq2=ng2-nr2*sint2;
    if (sq1 < 0.0001 || sq2 < 0.0001) {
      //avoid doing sqrt(negative) or sq1/sq2 equal to zero.
      while(sq1 < 0.0001 || sq2 < 0.0001) {
	theta_photon -= TOLERANCE;
	sint = sin(theta_photon);
	sint2 = sint*sint;
	sq1=nq2-nr2*sint2;
	sq2=ng2-nr2*sint2;
      }
    }
    double tant = sint/sqrt(1.-sint2); //Avoid computing tan(theta_photon)
    double q = theta_photon/twopi - floor(theta_photon/twopi);
    if( q>0.25 && q<0.75 ) tant *= -1.;
    term1 = (L_RAD-l_used)*tant;
    term2 = l_quartz*n_radiator*sint/sqrt(sq1);
    term3 = l_gap*n_radiator*sint/sqrt(sq2);
    if( fDebug > 1 ) {
      printf("term2=%f term3=%f\n",term2,term3);
    }
    trial =  term1 + term2 + term3;
 
    if(trial > dist_b) {
      upper = theta_photon;
      theta_photon = theta_photon - (upper - lower)/2;
    } else {
      lower = theta_photon;
      theta_photon = theta_photon + (upper - lower)/2;
    }
    if(IterationNumber > 1000) {
      // cout << "LOOP ! Trial = " << trial << ", dist_b = " << dist_b << endl;
      // cout <<  "upper = " << upper << ", lower = " << lower << endl;
      break;
    }
  }

  // End of Guido part
  **************************************************************************/

  /*************************************************************************/
  // Vahe part

  while( fabs(dist_b-trial)>TOLERANCE ) {
    IterationNumber ++;
    //pre-calculation
    double sint = sin(theta_photon);
    double sint2 = sint*sint;
    sq1=nq2-nr2*sint2;
    sq2=ng2-nr2*sint2;

    double tmp_sin2;
    double tmp_th;
    double tmp_sin;
    
    if (sq1 < 0.0001 || sq2 < 0.0001) {
      //avoid doing sqrt(negative) or sq1/sq2 equal to zero.
      if(sq1<0.0001&&sq2<0.0001) {
	if(sq1>sq2) 
	  {
	    tmp_sin2=((ng2-0.0001)/nr2  );
	    tmp_sin=sqrt(tmp_sin2);
	    tmp_th =asin(tmp_sin);
	    
	    sq1=nq2-nr2*tmp_sin2;
	    sq2=ng2-nr2*tmp_sin2;
	    sint=tmp_sin;
	    sint2=tmp_sin2;
	    theta_photon=tmp_th;	      
	  }
	else
	  {
	    tmp_sin2=( (nq2-0.0001)/nr2  );
	    tmp_sin=sqrt(tmp_sin2);
	    tmp_th  =asin(tmp_sin);
	    
	    sq1=nq2-nr2*tmp_sin2;
	    sq2=ng2-nr2*tmp_sin2;
	    sint=tmp_sin;
	    sint2=tmp_sin2;
	    theta_photon=tmp_th;
	  }
	if(fDebug>0 )	cout<<"vahe1 sq1 sq2 = "<<sq1<<" "<<sq2<<" theta = "<<tmp_th<<endl;

	// char tt=getchar();
      }
      else if(sq1 >0.0001 && sq2 <0.0001)
	{
	  tmp_sin2=( (ng2-0.0001)/nr2  );
	  tmp_sin = sqrt(tmp_sin2);
	  tmp_th  =asin(tmp_sin);
	  
	  sq1=nq2-nr2*tmp_sin2;
	  sq2=ng2-nr2*tmp_sin2;
	  sint =tmp_sin;
	  sint2=tmp_sin2;
	  theta_photon=tmp_th;
	  if(fDebug>0 )  cout<<"vahe2 sq1 sq2 = "<<sq1<<" "<<sq2<<" theta = "<<tmp_th<<endl;
	}
      else if(sq1<0.0001 && sq2>0.0001)
	{
	  tmp_sin2=( (nq2-0.0001)/nr2  );
	  tmp_sin =sqrt(tmp_sin2);
	  tmp_th  =asin(tmp_sin);
	  
	  sq1=nq2-nr2*tmp_sin2;
	  sq2=ng2-nr2*tmp_sin2;
	  sint =tmp_sin;
	  sint2=tmp_sin2;
	  theta_photon=tmp_th;	  
	  if(fDebug>0 ) cout<<"vahe3 sq1 sq2 = "<<sq1<<" "<<sq2<<" theta = "<<tmp_th<<endl;

	} 
    }
    double tant = sint/sqrt(1.-sint2); //Avoid computing tan(theta_photon)
    double q = theta_photon/twopi - floor(theta_photon/twopi);
    if( q>0.25 &&  q<0.75 ) tant *= -1.;
    term1 = (L_RAD-l_used)*tant;
    term2 = l_quartz*n_radiator*sint/sqrt(sq1);
    term3 = l_gap*n_radiator*sint/sqrt(sq2);
    if( fDebug > 1 ) {
      //      printf("term2=%f term3=%f\n",term2,term3);
    }
    trial =  term1 + term2 + term3;
 
    if(trial > dist_b) {
      upper = theta_photon;
      theta_photon = theta_photon - (upper - lower)/2;
    } else {
      lower = theta_photon;
      theta_photon = theta_photon + (upper - lower)/2;
    }
    if(IterationNumber > 1000) {
      // cout << "LOOP ! Trial = " << trial << ", dist_b = " << dist_b << endl;
      // cout <<  "upper = " << upper << ", lower = " << lower << endl;
      break;
    }
  }

  // END of Vahe part
  /*************************************************************************/

  if( fDebug > 0 )
    printf(" theta photon = %5f ",theta_photon);

  if( fDebug > 2 ) {
    printf("***theta_photon***\n");
    printf("   theta_photon = %f DEG\n",theta_photon/M_PI*180);
    printf("  'dist_b' calculated with theta_photon obtained above is...\n");
    printf("   %f + %f + %f = %f : should be as close as original 'dist_b:%f'\n"
	   ,term1,term2,term3,trial,dist_b);
  }

  if( fDoBench ) fBench->Stop("Get_theta_photon");
  return theta_photon;
}


//_____________________________________________________________________________
//**************************************************************************
//double RecoAng(v1,v2,&v3,&v4,v5,v6,v7,v8)
//    double  v1: x-coordinate of the C-photon hit-pad (mm)
//    double  v2: y-coordinate of the C-photon hit-pad (mm)
//    double  v3: theta angle of the C-photon hit-pad (rad) (to be calculated)
//    double  v4: phi angle of the C-photon hit-pad (rad) (to be calculated)
//    double  v3: x- of the MIP hit-pad (mm)
//    double  v4: y- of the MIP hit-pad (mm)
//    double  v5: the polar angle of the MIP (rad)
//    double  v6:the azimuthal angle of the MIP (rad)
//  NOTE: Mind the unit of each parameter is correct.
//              
// This function performe a series of calculation to return the reconstructed
// Cherenkov angle(radians). 
//**************************************************************************

Double_t SBSGRINCH::RecoAng( Double_t x_photon, Double_t y_photon,
                          Double_t &theta_photon, Double_t & phi_photon, 
			  Double_t x_mip, Double_t y_mip,
			  Double_t theta_mip, Double_t phi_mip, 
			  Int_t Calculation_Kind)  const
{

  //    Calculation kind: 
  //    0:  medium angle ("normal" calculation with the photon supposed 
  //                      generated in the middle of the radiator).
  //   -1: maximum possible angle (photon generated with minimum energy 
  //                              at the end of the radiator)
  //    1: minimum possible angle (photon generated with the maximum energy 
  //                              at the beginning of the radiator)


  if( fDoBench ) fBench->Begin("RecoAng");

  Double_t ReturnValue = -100.0;

  theta_photon = Get_theta_photon(x_photon,y_photon,x_mip,y_mip,
				  theta_mip,phi_mip, Calculation_Kind);
  phi_photon = Get_phi_photon
    (x_photon,y_photon,x_mip,y_mip,theta_mip,phi_mip, Calculation_Kind);

  //NOTE: the variable "phi_photon" required in the next is defined 
  //      during the procedure for theta_photon.

  //The last calculation to get the reconstructed angle.

  Double_t term1 = sin(theta_mip)*sin(theta_photon);
  // Double_t term1 = sin(theta_mip);
  Double_t term2 = cos(phi_photon-phi_mip);
  Double_t term3 = cos(theta_mip);
  Double_t term4 = cos(theta_photon);
  Double_t term5 = term1*term2 + term3*term4;
  
  if( fabs(term5)<=1.0 ) 
    ReturnValue = acos(term5);

  if( fDebug>0 ) {
    printf(" cherenkov angle = %5f ",ReturnValue);
    printf("  %5f ",phi_photon);
    printf("  %5f ",ReturnValue);
    //printf("\n");
  }

  if( fDoBench ) fBench->Stop("RecoAng");
  return ReturnValue;

}

//_____________________________________________________________________________
//**************************************************************************
//int ClearNoise(int igold, int ResolvedFlag)

// This function try to cut away from the chi square analysis clusters that 
// are generated by random noise.
// It should be called when all the three chi squared test (one for each 
// hypothesis on the kind of particle crossing the Rich) are not reseanable. 
// The algorithm try to see if at least one reduced chi square reached a 
// reasonable probability (more than a defined parameter 
// "acceptable_chi2_prob") when one cluster is cut away from the analysis.   
//
// returned values and operation performed by ClearNoise on the clusters
//
// success = 0                no chi square values is fixed after getting rid 
//                            of one cluster.
//                            Maybe the clusters generated by the noise are
//                            two or more. There is however one cluster whose 
//                            contribution to each of the three chi squares 
//                            considered is bigger than any contribution from 
//                            other clusters (let us call this cluster A in 
//                            the following).
//                            ClearNoise 
//                              1) sets for the mip 
//                                a) chi2_corrected[i] = chi2_corrected[i] + 
//                                   - chi2A[i]
//                                  (i = 0 for pion, i = 1 for kaon and i = 2 
//                                  for proton; chi2_corrected[] is equal to 
//                                  chi2[] the first time ClearNoise is 
//                                  called; chi2A[i] is equal to the three
//                                  contributions of A to the three chi2 
//                                  values).
//                                b) N_chi2_Photon[i] = N_chi2_Photon[i] - 1.
//                              2) sets for the cluster A 
//                                 PionChi2AnalysisFlag = 0
//                                 KaonChi2AnalysisFlag = 0
//                                 ProtonChi2AnalysisFlag = 0
//                                 (A will not be considered further in the
//                                 analysis).
// success = -1               no chi square values is fixed after getting 
//                            rid of one cluster and there is no cluster that
//                            gives the worst contribution to all the three
//                            chi square values. Let us call A, B and C 
//                            the clusters that give the worst contribution
//                            to pion chi square, kaon chi square and proton
//                            chi square respectively and chi2A, chi2B and 
//                            chi2C these contributions. 
//                            ClearNoise 
//                              1) sets for the mip
//                                a) chi2_corrected[0] = chi2_corrected[0] + 
//                                   - chi2A[0]
//                                   chi2_corrected[1] = chi2_corrected[1] + 
//                                   - chi2B[1]
//                                   chi2_corrected[2] = chi2_corrected[2] + 
//                                   - chi2C[2]
//                                  (i = 0 for pion, i = 1 for kaon and i = 2 
//                                  for proton; chi2_corrected[] is equal to 
//                                  chi2[] the first time ClearNoise is 
//                                  called; chi2A[0], chi2B[1], chi2C[2] are 
//                                  equal to the contributions of A B, and C 
//                                  to the pion, kaon and proton chi square
//                                  values respectively.
//                                b) N_chi2_Photon[i] = N_chi2_Photon[i] - 1.
//                              2) sets for the cluster A
//                                           PionChi2AnalysisFlag = 0
//                                 for the cluster B
//                                           KaonChi2AnalysisFlag = 0
//                                 for the cluster C
//                                           ProtonChi2AnalysisFlag = 0
// success = -2               there are multiple solutions; that is:  
//                            cutting away one cluster (let us say cluster A)
//                            one reduced chi square (pion chi square for 
//                            example) is o.k., cutting away another cluster
//                            (let us say cluster B) another chi square 
//                            (kaon chi square for instance) is o.k.
//                            The action taken is the same seen in the case
//                            of success = -1.
// succes = -3                the number of clusters the analysis is performed
//                            is not bigger than one or than the minimum 
                              // desired
//                            chi square degree of freedom. No other steps  
//                            are performed.
// success > 0                one cluster generated by noise has been found.
//                            chi2_corrected[] array values are set equal
//                            to chi[2] array minus the contribution of the 
//                            cluster from noise. N_chi2_corrected_photon[] 
//                            is set equal to N_chi2_corrected_photon - 1.
//                            success is set equal to the number of chi 
//                            square value fixed.  
//  
// NOTE: because unresolved and resolved cluster analysis could be completly  
// indipendent the function needs ResolvedFlag flag to know if it has to
// work with unresolved or resolved clusters.
// igold is simple the index of the MIP ClearMoise has to deal with.
// fMIP sould be already inizialed befeore calling ClearNoise
//**************************************************************************

Int_t SBSGRINCH::ClearNoise(Int_t igold, Int_t ResolvedFlag)
{
  Int_t success = 0;
  Int_t Index;
  Int_t Pion_success = 0;
  Int_t Kaon_success = 0;
  Int_t Proton_success = 0;
  // FC: presently not used, it has to be debugged first
  // Int_t Electron_success = 0;
  Int_t Pion_Index = 0;
  Int_t Kaon_Index = 0;
  Int_t Proton_Index = 0;
  //  Int_t Electron_Index = 0;
  // FC: presently not used, it has to be debugged first
  Int_t Pion_chi2_degree_of_freedom;
  Int_t Kaon_chi2_degree_of_freedom;
  Int_t Proton_chi2_degree_of_freedom;
  Double_t Pion_Chi2;
  Double_t Kaon_Chi2;
  Double_t Proton_Chi2;
  Double_t Pion_Chi2_Max = 0.0;
  Double_t Kaon_Chi2_Max = 0.0;
  Double_t Proton_Chi2_Max = 0.0;

  Pion_chi2_degree_of_freedom = fMIP.GetN_chi2_corrected_phot_pi(ResolvedFlag);
  Kaon_chi2_degree_of_freedom = fMIP.GetN_chi2_corrected_phot_k(ResolvedFlag);
  Proton_chi2_degree_of_freedom = fMIP.GetN_chi2_corrected_phot_p(ResolvedFlag);

  if ( (Pion_chi2_degree_of_freedom == 1) || (Pion_chi2_degree_of_freedom <=
					      minimum_chi2_degree_of_freedom) )
    {
      // only one cluster in the analysis. No further steps possible or
      // minimum analysis cluster number desired reached.
      // Please note: at the moment the number of clusters involved in the
      // analysis is the same for each hypothesis on the kind of particle 
      // crossing the Rich.
      success = -3;
      return success;
    }

  Pion_Chi2 = fMIP.Getchi2_corrected_pi(ResolvedFlag);
  Kaon_Chi2 = fMIP.Getchi2_corrected_k(ResolvedFlag);
  Proton_Chi2 = fMIP.Getchi2_corrected_p(ResolvedFlag);

  Int_t nClusters = GetNumClusters();
  Int_t nResolvedClusters = GetNumResolvedClusters();
  Int_t istart = 0;
  Int_t istop = nClusters;
  if (ResolvedFlag == 1)
    {
      istart = nClusters;
      istop = nClusters+nResolvedClusters;
    }
  for( int k = istart; k < istop; k++ ) 
    {
      SBSGRINCH_Cluster* theCluster = 
	( k < nClusters ) ? GetCluster(k) : GetResolvedCluster(k-nClusters); 
      if( !theCluster || theCluster->IsMIP() ) 
	continue;  // ignore MIPs
      // Int_t ResolvedFlag = 0;
      if(k >= nClusters) ResolvedFlag = 1;
      if ( (theCluster->GetPionChi2AnalysisFlag()) && 
	  (theCluster->Getchi2_pi() >  Pion_Chi2_Max) )
      {
	// the cluster is involved in the pion analysis and is a candidate to 
	// be "The noise one". 
	Pion_Chi2_Max = theCluster->Getchi2_pi();
	Pion_Index = k; 
      }
      if ( (theCluster->GetKaonChi2AnalysisFlag()) && 
	  (theCluster->Getchi2_k() >  Kaon_Chi2_Max) )
      {
	// the cluster is involved in the kaon analysis and is a candidate to 
	// be "The noise one". 
	Kaon_Chi2_Max = theCluster->Getchi2_k();
	Kaon_Index = k; 
      }
      if ( (theCluster->GetProtonChi2AnalysisFlag()) && 
	  (theCluster->Getchi2_p() >  Proton_Chi2_Max) )
      {
	// the cluster is involved in the proton analysis and is a candidate 
	// to be "The noise one". 
	Proton_Chi2_Max = theCluster->Getchi2_p();
	Proton_Index = k; 
      } 
    }
      if ( TMath::Prob((Pion_Chi2 - Pion_Chi2_Max),
		       (Pion_chi2_degree_of_freedom - 1)) 
	   > acceptable_chi2_prob)
	Pion_success = 1;
      if ( TMath::Prob((Kaon_Chi2 - Kaon_Chi2_Max),
		       (Kaon_chi2_degree_of_freedom - 1)) 
	   > acceptable_chi2_prob)
	Kaon_success = 1;
      if ( TMath::Prob((Proton_Chi2 - Proton_Chi2_Max),
		       (Proton_chi2_degree_of_freedom -1)) 
	   > acceptable_chi2_prob)
	Proton_success = 1;
      success =  Pion_success + Kaon_success + Proton_success;
      // throw away the worst contribution to chi square values
      fMIP.Insert_N_chi2_corrected_Photon(0, Pion_chi2_degree_of_freedom - 1,
					  ResolvedFlag);
      fMIP.Insert_N_chi2_corrected_Photon(1, Kaon_chi2_degree_of_freedom - 1,
					  ResolvedFlag);
      fMIP.Insert_N_chi2_corrected_Photon(2, Proton_chi2_degree_of_freedom-1,
					  ResolvedFlag);
      Index = Pion_Index;
      SBSGRINCH_Cluster* theCluster = ( Index < nClusters ) ? GetCluster(Index) 
	: GetResolvedCluster(Index-nClusters);
      theCluster->SetPionChi2AnalysisFlag(kFALSE);
      fMIP.Insert_chi2_corrected(0, Pion_Chi2 - Pion_Chi2_Max, ResolvedFlag);
      Index = Kaon_Index;
      theCluster = ( Index < nClusters ) ? GetCluster(Index) 
	: GetResolvedCluster(Index-nClusters);
      theCluster->SetKaonChi2AnalysisFlag(kFALSE);
      fMIP.Insert_chi2_corrected(1, Kaon_Chi2 - Kaon_Chi2_Max, ResolvedFlag);
      Index = Proton_Index;
      theCluster = ( Index < nClusters ) ? GetCluster(Index) 
	: GetResolvedCluster(Index-nClusters);
      theCluster->SetProtonChi2AnalysisFlag(kFALSE);
      fMIP.Insert_chi2_corrected(2, Proton_Chi2 - Proton_Chi2_Max, 
				 ResolvedFlag);
      fMIP.Setchi2_corrected_prob(0, fMIP.Getchi2_corrected_pi(ResolvedFlag),
				fMIP.GetN_chi2_corrected_phot_pi(ResolvedFlag),
				  ResolvedFlag);
      fMIP.Setchi2_corrected_prob(1, fMIP.Getchi2_corrected_k(ResolvedFlag),
				fMIP.GetN_chi2_corrected_phot_k(ResolvedFlag),
				  ResolvedFlag);
      fMIP.Setchi2_corrected_prob(2, fMIP.Getchi2_corrected_p(ResolvedFlag),
				fMIP.GetN_chi2_corrected_phot_p(ResolvedFlag),
				  ResolvedFlag);
      if(success == 0)
	{	     
	  if( (Pion_Index != Kaon_Index) || (Pion_Index != Proton_Index)
	      || (Kaon_Index != Proton_Index) )
	    success = -1;
	      // There is no cluster that gives the worst contribution 
	      // to all the chi square values
	} 
      else
	{	     
	  if(success > 1)
	    {
	      if( (Pion_Index != Kaon_Index) || (Pion_Index != Proton_Index)
		  || (Kaon_Index != Proton_Index) )
	    success = -2;
	      // There is no cluster that gives the worst contribution 
	      // to all the chi square values. There is no unambiguos 
	      // solution
	    }
	}
      return success;
}

//_____________________________________________________________________________
// Old code, left here for reference, will not be compiled
#if 0
Int_t SBSGRINCH::FineProcess_old( TClonesArray& tracks )
{
  // The main RICH processing method. Here we
  //  
  // - Identify clusters
  // - Find the MIP(s) among the clusters
  // - Calculate photon angles
  // - Calculate particle probabilities

  Float_t theta_photon = 0, phi_photon= 0;

  if (fMaxNumHits<fNHits) {
    return -25;
  }  

  if( FindClusters() == 0 ) { 
    return -1;
  }

  if( FindMIP( tracks ) == 0 ) {
    return -2;
  }

  // For each cluster, find the associated MIP.
  // This is equivalent to identifying the Cherenkov rings.
  // FIXME: Here we make a crude approximation: for each cluster,
  // find the *closest* MIP.  This is correct only in the case of 
  // non-overlapping rings.

  Int_t last_imip = -1;
  Float_t x_mip = 0.0, y_mip = 0.0, p_mip = 0.0;
  Float_t theta_mip = 0.0, phi_mip = 0.0;
  THaTrack* MIPtrack = NULL;

  SBSGRINCH_Cluster* theCluster = fClusters;
  for( int k = 0; k < fNClusters; k++, theCluster++ ) {
    if( theCluster->IsMIP() ) continue;  // ignore MIPs

    Float_t mindist         = 1e10;
    SBSGRINCH_Cluster* theMIP = NULL;
    Int_t imip              = 0;
    for( int i = 0; i < tracks.GetLast()+1; i++ ) {

      if( !fMIPs[i] ) continue;
      Float_t dist = theCluster->Dist( fMIPs[i] );
      if( dist < mindist ) {
	mindist = dist;
	theMIP = fMIPs[i];
	imip = i;
      }
    }

    theCluster->SetMIP( theMIP );
    if( !theMIP ) continue;

    if( imip != last_imip ) {
      last_imip = imip;
      x_mip = theMIP->GetXcenter();
      y_mip = theMIP->GetYcenter();

      // Get momentum and angle of the track associated with each MIP

      MIPtrack = static_cast<THaTrack*>( tracks.At(imip) );
      //FIXME: wrong, Px etc are at the target!
      TVector3 tv( MIPtrack->GetLabPx(), 
		   MIPtrack->GetLabPy(), 
		   MIPtrack->GetLabPz() );
      p_mip = MIPtrack->GetP();
      if( p_mip == 0.0 ) continue;  // Urgh
      tv *= 1.0/p_mip;   // Unit vector along track direction

      // MIP angles are in spherical coordinates in RICH system, units rad

      theta_mip = acos( fZax.Dot(tv) );
      if (theta_mip == 0) 
	{
	  phi_mip = 0;
	}
      else
	{
	  phi_mip = acos( (fXax.Dot(tv))/sin(theta_mip));
	  if(fYax.Dot(tv) < 0)
	    {
	      phi_mip = 2.*3.1415926536 - phi_mip;
	    }
	  //phi_mip   = acos( (tv-fZax.Dot(tv)*fZax).Dot(fXax) );
	} 
    }      

    theCluster->SetTrack( MIPtrack );

    Float_t  x_photon = theCluster->GetXcenter();
    Float_t  y_photon = theCluster->GetYcenter();
    Float_t  angle = RecoAng( x_photon, y_photon,theta_photon, phi_photon,
			      x_mip, y_mip, theta_mip, phi_mip, 
			      0);
    theCluster->SetTheta_photon( theta_photon);
    theCluster->SetPhi_photon( phi_photon);
    theCluster->SetAngle( angle );
  }

  // FIXME: fill a PID structure here and associate it with MIPtrack...

  return 0;
}
#endif

//_____________________________________________________________________________

//_____________________________________________________________________________
void SBSGRINCH::PrintBenchmarks() const
{
  // Print benchmark results
  
  if( !fDoBench )
    return;

  fBench->Print("Clear");
  fBench->Print("Decode");
  fBench->Print("FineProcess");
  fBench->Print("FindClusters");
  fBench->Print("FindMIP");
  fBench->Print("ResolveClusters");

  cout << endl << "Breakdown of time spent in FineProcess:" << endl;
  fBench->Print("RecoAng");
  fBench->Print("Get_phi_photon");
  fBench->Print("Get_theta_photon");
  
  return;
}

ClassImp(SBSGRINCH)

