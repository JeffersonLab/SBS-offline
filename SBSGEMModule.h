#ifndef SBSGEMModule_H
#define SBSGEMModule_H

#include "THaSubDetector.h"
#include <vector>
#include <set>
#include <map>

//using namespace std;

class THaDetectorBase;
class THaEvData;
class THaRunBase;
class TH1F;
class TH2F;
class TClonesArray;

namespace SBSGEM {
  enum GEMaxis_t { kUaxis=0, kVaxis };
}

struct mpdmap_t {
  UInt_t crate;
  UInt_t slot;
  UInt_t mpd_id;
  UInt_t gem_id;
  UInt_t adc_id;
  UInt_t i2c;
  UInt_t pos;
  UInt_t invert;
  UInt_t axis; //needed to add axis to the decode map
};

//Clustering results can be held in a simple C struct, as a cluster is just a collection of basic data types (not even arrays):
//Each module will have an array of clusters as a data member:
struct sbsgemhit_t { //2D reconstructed hits
  //Module and layer info: I think we don't need here because it's already associated with the module containing the cluster!
  //Track info:
  Bool_t keep;     //Should this cluster be considered for tracking? We use this variable to implement "cluster quality" cuts (thresholds, XY ADC and time correlation, etc.)
  Bool_t ontrack;  //Is this cluster on any track?
  Int_t trackidx; //Index of track containing this cluster (within the array of tracks found by the parent SBSGEMTracker
  Int_t iuclust;  //Index in (1D) U cluster array of the "U" cluster used to define this 2D hit.
  Int_t ivclust;  //Index in (1D) V cluster array of the "V" cluster used to define this 2D hit.
  
  //Strip info: should we store a list of strips? We shouldn't need to if the clustering algorithm requires all strips in a cluster to be contiguous:
  //Some of this information is probably redundant with 1D clustering results stored in sbsgemcluster_t, but it may or may not be more efficient to store a copy here:
  /* UShort_t nstripu; //Number of strips in cluster along "U" axis */
  /* UShort_t nstripv; //Number of strips in cluster along "V" axis */
  /* UShort_t ustriplo; //lowest u strip index in cluster */
  /* UShort_t ustriphi; //highest u strip index in cluster */
  /* UShort_t ustripmaxADC; //index of U strip with largest ADC in cluster */
  /* UShort_t vstriplo; //lowest v strip index in cluster */
  /* UShort_t vstriphi; //highest v strip index in cluster; */
  /* UShort_t vstripmaxADC; //index of V strip with largest ADC in cluster */
  //Coordinate info:
  Double_t umom;  //umean - u of center of strip with max ADC (could be used to refine hit coordinate reconstruction, probably unnecessary), in units of strip pitch
  Double_t vmom;  //vmean - v of center of strip with max ADC (could be used to refine hit coordinate reconstruction, probably unnecessary), in units of strip pitch
  Double_t uhit;  //Reconstructed U position of cluster in local module/strip coordinates
  Double_t vhit;  //Reconstructed V position of cluster in local module/strip coordinates
  Double_t xhit;  //Reconstructed "X" position of cluster in local module/strip coordinates
  Double_t yhit;  //Reconstructed "Y" position of cluster in local module/strip coordinates
  Double_t xghit; //"Global" X coordinate of hit (in coordinate system of parent SBSGEMTracker: usually spectrometer TRANSPORT coordinate)
  Double_t yghit; //"Global" Y coordinate of hit (in coordinates of parent SBSGEMTracker)
  Double_t zghit; //"Global" Z coordinate of hit (in coordinates of parent SBSGEMTracker)
  //
  Double_t Ehit;  //Sum of all ADC values on all strips in the cluster: actually 1/2*( ADCX + ADCY ); i.e., average of cluster sums in X and Y
  /* Double_t Euhit; //Sum of all ADC values on U strips in the cluster; */
  /* Double_t Evhit; //Sum of all ADC values on V strips in the cluster; */
  Double_t thit;  //Average of ADC-weighted mean U strip time and V strip time
  /* Double_t tuhit; //Average time of U strips in cluster; */
  /* Double_t tvhit; //Average time of V strips in cluster; */
  Double_t thitcorr; //"Corrected" hit time (with any trigger or other corrections we might want to apply)
  /* Double_t tuhitcorr; //"Corrected" U hit time (with any trigger time or other corrections) */
  /* Double_t tvhitcorr; //"Corrected" V hit time */
  Double_t ADCasym; //(ADCU-ADCV)/(ADCU+ADCV)
  Double_t tdiff;   //tu - tv
  Double_t corrcoeff_clust; //"Cluster" level XY correlation coefficient
  Double_t corrcoeff_strip; //Correlation coefficient of best XY strip pair, used to "seed" the cluster
  
};
  
struct sbsgemcluster_t {  //1D clusters;
  UInt_t nstrips;
  UInt_t istriplo;
  UInt_t istriphi;
  UInt_t istripmax;
  std::vector<double> ADCsamples; //cluster-summed ADC samples (accounting for split fraction)
  Double_t hitpos_mean;  //ADC-weighted mean coordinate along the direction measured by the strip
  Double_t hitpos_sigma; //ADC-weighted RMS coordinate deviation from the mean along the direction measured by the strip
  Double_t clusterADCsum; //Sum of ADCs over all samples on all strips
  std::vector<double> stripADCsum; //Sum of individual strip ADCs over all samples on all strips; accounting for split fraction
  Double_t t_mean; //reconstructed hit time
  Double_t t_sigma; //unclear what we might use this for
  //Do we want to store the individual strip ADC Samples with the 1D clustering results? I don't think so; as these can be accessed via the decoded strip info.
  std::vector<UInt_t> hitindex; //position in decoded hit array of each strip in the cluster:
  bool keep;
};


//Should these be hardcoded? Probably not!
/* #define N_APV25_CHAN    128 */
/* #define N_MPD_TIME_SAMP 6 */
/* #define MPDMAP_ROW_SIZE 8 */

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//      GEM module class: a GEM module consists of one part of a tracking unit;
//      Basic assumptions:
//      1. Rectangular active area
//      2. Strip readout based on two non-parallel strip orientations, denoted "U" and "V", with some constant number of time samples
//      3. A logical GEM layer or GEM plane will consist of one or more GEM modules, and a GEM "Tracker"
//         will consist of one or more GEM "layers" or "planes". We need a separate "module" class, but we probably don't need a separate "GEM layer" class,
//         since the layer will merely be an integer index that logically groups one or more modules together as a single plane for purposes of tracking.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Is the use of THaSubDetector appropriate here? Or should we just use THaDetector or similar?
  
class SBSGEMModule : public THaSubDetector {
 public:

  SBSGEMModule( const char *name, const char *description = "",
		THaDetectorBase* parent = 0 );

  virtual ~SBSGEMModule();

  virtual void    Clear( Option_t* opt="" ); //should be called once per event
  virtual Int_t   Decode( const THaEvData& );
  virtual void    Print( Option_t* opt="" ) const;

  virtual Int_t   ReadDatabase(const TDatime& );
  virtual Int_t   DefineVariables( EMode mode );
  //We are overriding THaDetectorBase's ReadGeometry method for SBSGEMModule, because we use a different definition of the angles:
  virtual Int_t   ReadGeometry( FILE* file, const TDatime& date, 
			      Bool_t required = true );
  
  virtual Int_t   Begin( THaRunBase* r=0 );
  virtual Int_t   End( THaRunBase* r=0 );

  //Don't call this method directly, it is called by find_2Dhits. Call that instead:
  void find_clusters_1D(SBSGEM::GEMaxis_t axis, Double_t constraint_center=0.0, Double_t constraint_width=1000.0); //Assuming decode has already been called; this method is fast so we probably don't need to implement constraint points and widths here, or do we?
  void find_2Dhits(); // Version with no arguments assumes no constraint points
  void find_2Dhits(TVector2 constraint_center, TVector2 constraint_width); // Version with TVector2 arguments 

  // fill the 2D hit arrays from the 1D cluster arrays:
  void fill_2D_hit_arrays(); 

  //Filter 1D hits by criteria possibly to include ADC threshold, cluster size
  void filter_1Dhits(SBSGEM::GEMaxis_t axis);
  
  //Filter 2D hits by criteria possibly to include ADC X/Y asymmetry, cluster size, time correlation, (lack of) overlap, possibly others:
  void filter_2Dhits(); 
  
  //Utility function to calculate correlation coefficient between U and V time samples:
  Double_t CorrCoeff( int nsamples, std::vector<double> Usamples, std::vector<double> Vsamples );

  //Utility functions to compute "module local" X and Y coordinates from U and V (strip coordinates) to "transport" coordinates (x,y) and vice-versa:
  TVector2 UVtoXY( TVector2 UV );
  TVector2 XYtoUV( TVector2 XY );

  //function to convert from APV channel number to strip number ordered by position:
  Int_t GetStripNumber( UInt_t rawstrip, UInt_t pos, UInt_t invert );

  void PrintPedestals( std::ofstream &dbfile, std::ofstream &daqfile_ped, std::ofstream &daqfile_cmr );
  
  bool fIsDecoded;

  //UShort_t GetLayer() const { return fLayer; }

  //std::vector<sbsgemhit_t> GetHitList() { return fHits; }
  
  //If we are going to declare all these data members private, we will need to write public getters and setters for at least the information required by the tracker classes:
  //private:

  //Let's just make all this stuff public because we're lazy. A lot of this is used by the track-finding routines
  //anyway
  
  //Decode map information: 
  std::vector<mpdmap_t>    fMPDmap; //this may need to be modified
  std::vector<Int_t>       fChanMapData;

  //some convenience maps: 
  std::map<Int_t, Int_t> fAPVch_by_Ustrip;
  std::map<Int_t, Int_t> fAPVch_by_Vstrip;
  std::map<Int_t, Int_t> fMPDID_by_Ustrip;
  std::map<Int_t, Int_t> fMPDID_by_Vstrip;
  std::map<Int_t, Int_t> fADCch_by_Ustrip;
  std::map<Int_t, Int_t> fADCch_by_Vstrip;
  
  Bool_t fPedestalMode;
  Double_t fZeroSuppressRMS;
  Bool_t fZeroSuppress;
  //Moved to the MPD module class:
  Bool_t fOnlineZeroSuppression; //this MIGHT be redundant with fZeroSuppress (or not)

  //move these to trackerbase:
  //Double_t fSigma_hitpos;   //sigma parameter controlling resolution entering track chi^2 calculation
  //Double_t fSigma_hitshape; //Sigma parameter controlling hit shape for cluster-splitting algorithm.
    
  //Note on ROOT data types:
  // Char_t/UChar_t is one byte signed unsigned = 8 bits = 0..255 for unsigned
  // Short_t/UShort_t is two bytes signed/unsigned = 16 bits = 0..65535 for unsigned
  // Int_t/UInt_t is four bytes signed/unsigned = 32 bits
  // Long_t/ULong_t is eight bytes = 64 bit
  // Float_t/Double_t is four bytes/eight bytes

  UChar_t fN_APV25_CHAN;     //number of APV25 channels, default 128
  UChar_t fN_MPD_TIME_SAMP;  //number of MPD time samples, default = 6
  UShort_t fMPDMAP_ROW_SIZE; //MPDMAP_ROW_SIZE: default = 9, let's not hardcode
  UShort_t fNumberOfChannelInFrame; //default 128, not clear if this is used: This might not be constant, in fact

  Double_t fSamplePeriod; //for timing calculations: default = 25 ns.

  //variables defining rectangular track search region constraint (NOTE: these will change event-to-event, they are NOT constant!)
  Double_t fxcmin, fxcmax;
  Double_t fycmin, fycmax;
  
  
  //BASIC DECODED STRIP HIT INFO:
  //By the time the information is populated here, the ADC values are already assumed to be pedestal/common-mode subtracted and/or zero-suppressed as appropriate:
  UInt_t fNstrips_hit; //total Number of strips fired (after common-mode subtraction and zero suppression)
  UInt_t fNstrips_hitU; //total number of U strips fired
  UInt_t fNstrips_hitV; //total number of V strips fired
  
  //Map strip indices in this array:
  //key = U or V  strip number, mapped value is position of that strip's information in the "decoded strip" arrays below:
  std::map<UInt_t, UInt_t> fUstripIndex; 
  std::map<UInt_t, UInt_t> fVstripIndex; 
  
  std::vector<UInt_t> fStrip;  //Strip index of hit (these could be "U" or "V" generalized X and Y), assumed to run from 0..N-1
  std::vector<SBSGEM::GEMaxis_t>  fAxis;  //We just made our enumerated type that has two possible values, makes the code more readable (maybe)
  std::vector<std::vector<Double_t> > fADCsamples; //2D array of ADC samples by hit: Outer index runs over hits; inner index runs over ADC samples
  std::vector<std::vector<Int_t> > fRawADCsamples; //2D array of raw (non-baseline-subtracted) ADC values.
  std::vector<Double_t> fADCsums;
  std::vector<Double_t> fStripADCavg;
  std::vector<UInt_t> fStripIsU; // is this a U strip? 0/1
  std::vector<UInt_t> fStripIsV; // is this a V strip? 0/1
  std::vector<bool> fKeepStrip; //keep this strip?
  std::vector<UInt_t> fMaxSamp; //APV25 time sample with maximum ADC;
  std::vector<Double_t> fADCmax; //largest ADC sample on the strip:
  std::vector<Double_t> fTmean; //ADC-weighted mean strip time:
  std::vector<Double_t> fTsigma; //ADC-weighted RMS deviation from the mean
  std::vector<Double_t> fTcorr; //Strip time with all applicable corrections; e.g., trigger time, etc.
  
  ////// (1D and 2D) Clustering results (see above for definition of struct sbsgemcluster_t and struct sbsgemhit_t):

  UInt_t fNclustU; // number of U clusters found
  UInt_t fNclustV; // number of V clusters found
  std::vector<sbsgemcluster_t> fUclusters; //1D clusters along "U" direction
  std::vector<sbsgemcluster_t> fVclusters; //1D clusters along "V" direction

  UInt_t fN2Dhits; // number of 2D hits found in region of interest:
  std::vector<sbsgemhit_t> fHits; //2D hit reconstruction results

  /////////////////////// Global variables that are more convenient for ROOT Tree/Histogram Output (to the extent needed): ///////////////////////
  //Raw strip info:
  //std::vector<UInt_t> fStripAxis; //redundant with fAxis, but more convenient for Podd global variable definition
  std::vector<Double_t> fADCsamples1D; //1D array to hold ADC samples; should end up with dimension fNstrips_hit*fN_MPD_TIME_SAMP
  std::vector<Int_t> fStripTrackIndex; // If this strip is included in a cluster that ends up on a good track, we want to record the index in the track array of the track that contains this strip.
  std::vector<Int_t> fRawADCsamples1D;
  /////////////////////// End of global variables needed for ROOT Tree/Histogram output ///////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////
  //  Define some global variables for tests/cuts                                             //
  //bool fDidHit; //a good track passed within the active area of this module AND a hit occurred within some small distance from the projected track
  //bool fShouldHit; //a good track passed within the active area of this module
  //moved the above to TrackerBase
  //////////////////////////////////////////////////////////////////////////////////////////////
  
  //Constant, module-specific parameters:
  UShort_t fModule; // Module index within a tracker. Should be unique! Since this is a GEM module class, this parameter should be unchanging
  UShort_t fLayer;  // Layer index of this module within a tracker. Since this is a GEM module class, this parameter should be unchanging

  UInt_t fNstripsU; // Total number of strips in this module along the generic "U" axis
  UInt_t fNstripsV; // Total number of strips in this module along the generic "V" axis

  //Pedestal means and RMS values for all channels:
  std::vector<Double_t> fPedestalU, fPedRMSU; 
  std::vector<Double_t> fPedestalV, fPedRMSV;

  //To be determined from channel map/strip count information:
  UShort_t fNAPVs_U; //Number of APV cards per module along "U" strip direction; this is typically 8, 10, or 12, but could be larger for U/V GEMs
  UShort_t fNAPVs_V; //Number of APV cards per module along "V" strip direction; 
  //UInt_t fNTimeSamplesADC; //Number of ADC time samples (this could be variable in principle, but should be the same for all strips within a module within a run) redundant with fN_MPD_TIME_SAMP

  // Should we define the parameters controlling cluster-finding in the Module class? Why not:
  std::vector<Double_t> fUgain; // "gain match" coefficients for U strips by APV card;
  std::vector<Double_t> fVgain; // "gain match" coefficients for V strips by APV card;
    
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  //     CLUSTERING PARAMETERS (to be read from database and/or given sensible default values)        //
  ///////////////////////////////////////////////////////////////////////////////////////////////////

  Double_t fThresholdSample; //Threshold on the (gain-matched and pedestal-subtracted) max. sample on a strip to keep that strip for clustering 
  Double_t fThresholdStripSum; //Threshold on the sum of (pedestal-subtracted) ADC samples on a strip
  Double_t fThresholdClusterSum; //Threshold on the sum of (pedestal-subtracted) ADC samples 

  Double_t fADCasymCut;       // Filtering criterion for ADC X/Y (or U/V) asymmetry
  Double_t fTimeCutUVdiff;    // Filtering criterion for ADC X/Y (or U/V) time difference (this is distinct from any timing cuts relative to
  //trigger or reference timing at the individual strip level
  //Parameters controlling cluster splitting and insignificant peak elimination based on "peak prominence" calculation
  Double_t fThresh_2ndMax_nsigma;   //Number of sigmas above noise level for minimum peak prominence in splitting overlapping clusters
  Double_t fThresh_2ndMax_fraction; //Peak prominence threshold as a fraction of peak height for splitting overlapping clusters

  UShort_t fMaxNeighborsU_totalcharge; //Only strips within +/- fMaxNeighborsU of the peak can be added to a cluster for total charge calculation
  UShort_t fMaxNeighborsV_totalcharge; //Only strips within +/- fMaxNeighborsV of the peak can be added to a cluster for total charge calculation

  //Only strips within these limits around the peak can be used for hit position reconstruction
  UShort_t fMaxNeighborsU_hitpos; 
  UShort_t fMaxNeighborsV_hitpos; 

  Double_t fSigma_hitshape; //Sigma parameter controlling hit shape for cluster-splitting algorithm.
  
  //GEOMETRICAL PARAMETERS:
  Double_t fUStripPitch;    //strip pitch along U, will virtually always be 0.4 mm
  Double_t fVStripPitch;    //strip pitch along V, will virtually always be 0.4 mm
  Double_t fUAngle;         //Angle between U strips and "X" axis of TRANSPORT coordinates;
  Double_t fVAngle;         //Angle between V strips and "X" axis of TRANSPORT coordinates;
  Double_t fPxU;            //U Strip X projection = cos( UAngle );
  Double_t fPyU;            //U Strip Y projection = sin( UAngle );
  Double_t fPxV;            //V Strip X projection = cos( VAngle );
  Double_t fPyV;            //V Strip Y projection = sin( VAngle );
  
  //Note: These size variables are redundant with THaDetectorBase
  //Double_t fLx;             //Full width of module active area along "X" (usually + dispersive direction)
  //Double_t fLy;             //Full width of module active area along "Y" (usually + dispersive direction);
	
  //Alignment: These alignment variables are redundant with THaDetectorBase, so re-use those:
  /* TVector3 fModPos;        //Module position in BigBite TRANSPORT coordinates (or other global coordinate system); position within the overall tracking unit: */
  /* Double_t fModAlphaX;     //Module rotational offset about X axis of BigBite TRANSPORT coordinates (or other global coordinate system)  */
  /* Double_t fModAlphaY;     //Module rotational offset about Y axis of BigBite TRANSPORT coordinates (or other global coordinate system)  */
  /* Double_t fModAlphaZ;     //Module rotational offset about Z axis of BigBite TRANSPORT coordinates (or other global coordinate system)  */
  /* TRotation fModRot;       //Module rotation matrix to go from local --> global */
  /* TRotation fModRotInv;    //Module rotation matrix to go from global --> local */

  //For geometry parameters, we will re-use the THaDetectorBase functionalities as much as possible
       
  Bool_t fIsMC;//we kinda want this guy no matter what don't we...

  //Efficiency histograms:
  TH1F *fhdidhitx;
  TH1F *fhdidhity;
  TH2F *fhdidhitxy;

  TH1F *fhshouldhitx;
  TH1F *fhshouldhity;
  TH2F *fhshouldhitxy;

  //we should let the user configure this: this is set at the "tracker level" which then propagates down to all the modules:
  bool fMakeEfficiencyPlots;

  //Pedestal plots: only generate if pedestal mode = true:
  
  TH2F *hrawADCs_by_stripU; //raw adcs by strip, no corrections, filled for each SAMPLE:
  TH2F *hrawADCs_by_stripV; //raw adcs by strip, no corrections, filled for each SAMPLE:
  TH2F *hcommonmode_subtracted_ADCs_by_stripU; //common-mode subtracted ADCS without ped subtraction
  TH2F *hcommonmode_subtracted_ADCs_by_stripV; 
  TH2F *hpedestal_subtracted_ADCs_by_stripU; //common-mode AND pedestal subtracted ADCs
  TH2F *hpedestal_subtracted_ADCs_by_stripV;

  TH2F *hpedestal_subtracted_rawADCs_by_stripU; //pedestal-subtracted ADCs w/o common-mode correction
  TH2F *hpedestal_subtracted_rawADCs_by_stripV;

  //Summed over all strips, pedestal-subtracted (but not common-mode subtracted) ADCs:
  TH1F *hpedestal_subtracted_rawADCsU;
  TH1F *hpedestal_subtracted_rawADCsV;

  //Summed over all strips, pedestal-subtracted and  common-mode subtracted ADCs:
  TH1F *hpedestal_subtracted_ADCsU;
  TH1F *hpedestal_subtracted_ADCsV;

  //in pedestal-mode analysis, we 
  TH2F *hcommonmode_mean_by_APV_U;
  TH2F *hcommonmode_mean_by_APV_V;
  
  //Comment out for now, uncomment later if we deem these interesting:
  // TClonesArray *hrawADCs_by_strip_sampleU;
  // TClonesArray *hrawADCs_by_strip_sampleV;
  // TClonesArray *hcommonmode_subtracted_ADCs_by_strip_sampleU;
  // TClonesArray *hcommonmode_subtracted_ADCs_by_strip_sampleV;

  // TClonesArray *hpedestal_subtracted_ADCs_by_strip_sampleU;
  // TClonesArray *hpedestal_subtracted_ADCs_by_strip_sampleV;
  
  ClassDef(SBSGEMModule,0);

};




#endif
