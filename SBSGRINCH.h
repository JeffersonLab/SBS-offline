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
#include <stdint.h>
#include <map>

class THaTrack;
class THaBenchmark;

const double m_el = 0.5110034e-3; // electron mas in GeV
 const double m_pi = 139.57018e-3; // pion mass in GeV
// const double m_pi = 0.5110034e-3; // FC: NOT TRUE, JUST FOR TESTING!
const double m_ka = 493.677e-3;   // kaon mass in GeV
const double m_pr = 938.272e-3;   // proton mass in GeV
static const Float_t kBig = 1e38;

class SBSGRINCH : public THaPidDetector {
  
public:

  SBSGRINCH( const char* name, const char* description="", 
	   THaApparatus* apparatus=NULL );
  virtual ~SBSGRINCH();
  
  virtual void         Clear( Option_t* opt="" );
  virtual Int_t        Decode( const THaEvData& );
  virtual Int_t        CoarseProcess( TClonesArray& tracks );
  virtual Int_t        FineProcess( TClonesArray& tracks );

  SBSGRINCH_Hit*         GetHit(Int_t i) const 
    { return (SBSGRINCH_Hit*)fHits->At(i); }
  /* SBSGRINCH_Hit*         GetResolvedHit(Int_t i) const  */
  /*   { return (SBSGRINCH_Hit*)fResolvedHits->At(i); } */
  SBSGRINCH_Cluster*     GetCluster(Int_t i) const 
    { return (SBSGRINCH_Cluster*)fClusters->At(i); }
  /* SBSGRINCH_Cluster*     GetResolvedCluster(Int_t i) const  */
  /*   { return (SBSGRINCH_Cluster*)fResolvedClusters->At(i); } */
  Int_t                GetNumHits() const 
    { return fHits->GetLast()+1; }
  Int_t                GetNumClusters() const 
    { return fClusters->GetLast()+1; }
  /* Int_t                GetNumResolvedHits() const  */
  /*   { return fResolvedHits->GetLast()+1; } */
  /* Int_t                GetNumResolvedClusters() const  */
  /*   { return fResolvedClusters->GetLast()+1; } */
  
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

  //Int_t             fNypads;  // Number of pads along y (transverse)

  TClonesArray*     fHits;          // Array of hits for each event
  TClonesArray*     fClusters;      // Clusters of hits
  /* TClonesArray*     fResolvedHits;  // Hits of resolved clusters */
  /* TClonesArray*     fResolvedClusters; // Resolved clusters */


  //RICH parameters from database
  /*
  Double_t Z_CkovIn; // Z of the Cherenkov entrance window in the spectrometer
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

  Double_t fNradiator;     // radiator index of refraction;
  Double_t fLradiator;     // radiator length on central ray;
  */
  Double_t fZCkovIn;       // Z of the entrance window in the spectrometer central ray;
  Int_t    fNPMTs;         // number of PMTs
  Int_t    fNPMTrows;      // number of PMT rows
  Int_t    fNPMTcolsMax;   // max number of PMT columns 
  Double_t fPMTmatrixHext; // horizontal extension, in m, of the PMT matrix (from lower PMT center to higher PMT center)
  Double_t fPMTmatrixVext; // vertical extension, in m, of the PMT matrix (from left PMT center to right PMT center)
  Double_t fPMTdistX;      // X distance between the center of 2 PMT tubes in consecutive rows, in m
  Double_t fPMTdistY;      // Y distance between the center of 2 PMT tubes in consecutive columns, in m
  Double_t fX_TCPMT;       // X position of the top close PMT center in the PMT matrix (transport coord)
  Double_t fY_TCPMT;       // Y position of the top close PMT center in the PMT matrix (transport coord)
  
  Int_t   fMaxNumHits;
  
  Bool_t  fDoResolve;    // true = resolve overlapping clusters
  Bool_t  fDoTimeFilter; // true = filter the hits in each cluster with timing
  /* Int_t   fNseg;         // Number of x segments */
  /* Double_t* fXseg;       // Array of x segmentation boudaries and offsets */
  
  Double_t fTrackX;      // x pos of Golden Track in RICH plane
  Double_t fTrackY;      // y pos of Golden Track in RICH plane

  Bool_t         fDoBench;         //Collect detailed timing statistics
  THaBenchmark*  fBench;           //Counters for timing statistics
  
  Bool_t GetEdge(uint32_t &dataword);
  //void    DecipherVetrocWord(uint32_t VetrocWord, Bool_t& edge, Short_t& chan, UShort_t& time);
  void    DeleteClusters();
  Int_t   FindClusters();
  Int_t   MatchClustersWithTracks( TClonesArray& tracks );
  Int_t   CleanClustersWithTime();
  //Int_t   ResolveClusters();
  //Double_t Cherenkov_Angle(double mass, double momentum) const;

  virtual Int_t  ReadDatabase( const TDatime& date );
  virtual Int_t  DefineVariables( EMode mode = kDefine );

private:
  
  std::map< int, std::pair< int, int > > map_chan_tdcs;
  
  /*
  // Fix me: to insert in the data base
  Double_t minimum_chi2_degree_of_freedom; // minum number of degree of freedom
                                      // (that is clusters) on desires
                                      // perform chisquare test
  Double_t clear_noise_trial_maximum_number;
  Double_t acceptable_chi2_prob; // the probability a reduced chi2 is 
                                 // "acceptable"
  Double_t epsilon; // epsilon parameter in the Maximum Likeood Algorithm 
  */
  ClassDef(SBSGRINCH,0)   //The Hall A RICH
};

#endif









