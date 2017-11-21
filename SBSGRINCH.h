#ifndef ROOT_SBSGRINCH
#define ROOT_SBSGRINCH

//////////////////////////////////////////////////////////////////////////
//
// SBSGRINCH
//
// The Hall A RICH
//
//////////////////////////////////////////////////////////////////////////

#include "THaPidDetector.h"
#include "SBSGRINCH_ClusterList.h"
#include "TBits.h"
#include "TClonesArray.h"

class THaTrack;
class THaBenchmark;

const double m_el = 0.5110034e-3; // electron mas in GeV
 const double m_pi = 139.57018e-3; // pion mass in GeV
// const double m_pi = 0.5110034e-3; // FC: NOT TRUE, JUST FOR TESTING!
const double m_ka = 493.677e-3;   // kaon mass in GeV
const double m_pr = 938.272e-3;   // proton mass in GeV

class SBSGRINCH : public THaPidDetector {
  
public:

  SBSGRINCH( const char* name, const char* description="", 
	   THaApparatus* apparatus=NULL );
  virtual ~SBSGRINCH();
  
  virtual void         Clear( Option_t* opt="" );
  virtual Int_t        Decode( const THaEvData& );
  virtual Int_t        CoarseProcess( TClonesArray& tracks );
  virtual Int_t        FineProcess( TClonesArray& tracks );

  void                 ReadBadPads(Char_t* infilename);
  Int_t                ReadData( FILE *infile );
  SBSGRINCH_Hit*         GetHit(Int_t i) const 
    { return (SBSGRINCH_Hit*)fHits->At(i); }
  SBSGRINCH_Hit*         GetResolvedHit(Int_t i) const 
    { return (SBSGRINCH_Hit*)fResolvedHits->At(i); }
  SBSGRINCH_Cluster*     GetCluster(Int_t i) const 
    { return (SBSGRINCH_Cluster*)fClusters->At(i); }
  SBSGRINCH_Cluster*     GetResolvedCluster(Int_t i) const 
    { return (SBSGRINCH_Cluster*)fResolvedClusters->At(i); }
  Int_t                GetNumHits() const 
    { return fHits->GetLast()+1; }
  Int_t                GetNumClusters() const 
    { return fClusters->GetLast()+1; }
  Int_t                GetNumResolvedHits() const 
    { return fResolvedHits->GetLast()+1; }
  Int_t                GetNumResolvedClusters() const 
    { return fResolvedClusters->GetLast()+1; }

  // FIX ME the latter ones return non sense if decode has not been processed

  //  Int_t                GetTIRDat() const { return fTIRDat; }
  Int_t                GetMaxNumHits() const { return fMaxNumHits; }
  void                 SetMaxNumHits( Int_t MaxNumHit ) 
    { fMaxNumHits=MaxNumHit; }
  void                 EnableClusterResolving( Bool_t flag = kTRUE )
  { fDoResolve = flag; }
  void                 EnableBenchmarks( Bool_t b = kTRUE )
  { fDoBench = b; }
  void                 PrintBenchmarks() const;

protected:

  Int_t             fNypads;  // Number of pads along y (transverse)

  TClonesArray*     fHits;          // Array of hits for each event
  TClonesArray*     fClusters;      // Clusters of hits
  TClonesArray*     fResolvedHits;  // Hits of resolved clusters
  TClonesArray*     fResolvedClusters; // Resolved clusters

  //RICH parameters from database
  Double_t L_RAD,l_quartz,l_gap;    //length of radiator,quartz,proxiity gap
  Double_t l_emission;              //photon emission depth in the radiator.
  Double_t n_radiator,n_quartz,n_gap; //the refraction indices 
  Double_t n_quartz_min, n_quartz_max; 
  // Minimum and maximun refraction index value for the quartz in the range 
  // of Cherenkov photon energy the PMT are sensible at.
  Double_t n_radiator_min, n_radiator_max; 
  // Minimum and maximun refraxion index value for the radiator in the range 
  // of Cherenkov photon energy the PMT are sensible at.
  Double_t fiducial_zone_range; 
  // angular range of the fiducial zone around the expected angle for each 
  // kind of particle
  Double_t cluster_distribution_sigma;
  // sigma of single cluster angular distribution.
  //Double_t PMTinterdist;// distance between two PMTs in a row, or between 2 rows of PMTs
  
  
  Int_t   fMaxNumHits;            


  Bool_t  fDoResolve;    // true = resolve overlapping clusters
  Int_t   fNseg;         // Number of x segments
  Double_t* fXseg;       // Array of x segmentation boudaries and offsets

  Double_t fTrackX;      // x pos of Golden Track in RICH plane
  Double_t fTrackY;      // y pos of Golden Track in RICH plane

  Bool_t         fDoBench;         //Collect detailed timing statistics
  THaBenchmark*  fBench;           //Counters for timing statistics

  void    DeleteClusters();
  Int_t   FindClusters();
  Int_t   ResolveClusters();
  Double_t Cherenkov_Angle(double mass, double momentum) const;

  virtual Int_t  ReadDatabase( const TDatime& date );
  virtual Int_t  DefineVariables( EMode mode = kDefine );

private:

  // Fix me: to insert in the data base
  Double_t minimum_chi2_degree_of_freedom; // minum number of degree of freedom
                                      // (that is clusters) on desires
                                      // perform chisquare test
  Double_t clear_noise_trial_maximum_number;
  Double_t acceptable_chi2_prob; // the probability a reduced chi2 is 
                                 // "acceptable"
  Double_t epsilon; // epsilon parameter in the Maximum Likeood Algorithm 

  ClassDef(SBSGRINCH,0)   //The Hall A RICH
};

#endif









