#include "THaSubDetector.h"
#include <vector>

using namespace std;

class THaDetectorBase;
class THaEvData;
class THaRunBase;


struct mpdmap_t {
    UInt_t crate;
    UInt_t slot;
    UInt_t mpd_id;
    UInt_t gem_id;
    UInt_t adc_id;
    UInt_t i2c;
    UInt_t pos;
    UInt_t invert;
};

//Clustering results can be held in a simple C struct, as a cluster is just a collection of basic data types (not even arrays):
//Each module will have an array of clusters as a data member:
struct sbsgemcluster_t {
  //Module and layer info: I think we don't need here because it's already associated with the module containing the cluster!
  //Track info:
  Bool_t keep;     //Should this cluster be considered for tracking? We use this variable to implement "cluster quality" cuts (thresholds, XY ADC and time correlation, etc.)
  Bool_t ontrack;  //Is this cluster on any track?
  UInt_t trackidx; //Index of track containing this cluster (within the array of tracks found by the parent SBSGEMTracker
  //Strip info: should we store a list of strips? We shouldn't need to if the clustering algorithm requires all strips in a cluster to be contiguous:
  UShort_t nstripu; //Number of strips in cluster along "U" axis
  UShort_t nstripv; //Number of strips in cluster along "V" axis
  UShort_t ustriplo; //lowest u strip index in cluster
  UShort_t ustriphi; //highest u strip index in cluster
  UShort_t ustripmaxADC; //index of U strip with largest ADC in cluster
  UShort_t vstriplo; //lowest v strip index in cluster
  UShort_t vstriphi; //highest v strip index in cluster;
  UShort_t vstripmaxADC; //index of V strip with largest ADC in cluster
  //Coordinate info:
  Double_t umom;  //umean - u of center of strip with max ADC (could be used to refine hit coordinate reconstruction, probably unnecessary)
  Double_t vmom;  //vmean - v of center of strip with max ADC (could be used to refine hit coordinate reconstruction, probably unnecessary)
  Double_t uhit;  //Reconstructed U position of cluster in local module/strip coordinates
  Double_t vhit;  //Reconstructed V position of cluster in local module/strip coordinates
  Double_t xhit;  //Reconstructed "X" position of cluster in local module/strip coordinates
  Double_t yhit;  //Reconstructed "Y" position of cluster in local module/strip coordinates
  Double_t xghit; //"Global" X coordinate of hit (in coordinate system of parent SBSGEMTracker: usually spectrometer TRANSPORT coordinate)
  Double_t yghit; //"Global" Y coordinate of hit (in coordinates of parent SBSGEMTracker)
  Double_t zghit; //"Global" Z coordinate of hit (in coordinates of parent SBSGEMTracker)
  //
  Double_t Ehit;  //Sum of all ADC values on all strips in the cluster: actually 1/2*( ADCX + ADCY ); i.e., average of cluster sums in X and Y
  Double_t Euhit; //Sum of all ADC values on U strips in the cluster;
  Double_t Evhit; //Sum of all ADC values on V strips in the cluster;
  Double_t thit;  //Average of ADC-weighted mean U strip time and V strip time
  Double_t tuhit; //Average time of U strips in cluster;
  Double_t tvhit; //Average time of V strips in cluster;
  Double_t thitcorr; //"Corrected" hit time (with any trigger or other corrections we might want to apply)
  Double_t tuhitcorr; //"Corrected" U hit time (with any trigger time or other corrections)
  Double_t tvhitcorr; //"Corrected" V hit time
  Double_t ADCasym; //(ADCU-ADCV)/(ADCU+ADCV)
  Double_t tdiff;   //tu - tv
  Double_t corrcoeff_clust; //"Cluster" level XY correlation coefficient
  Double_t corrcoeff_strip; //Correlation coefficient of best XY strip pair, used to "seed" the cluster
  
}

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

        virtual void    Clear( Option_t* opt="" );
        virtual Int_t   Decode( const THaEvData& );
        virtual void    Print( Option_t* opt="" ) const;

        virtual Int_t   ReadDatabase(const TDatime& );
        virtual Int_t   DefineVariables( EMode mode );

        virtual Int_t   Begin( THaRunBase* r=0 );
        virtual Int_t   End( THaRunBase* r=0 );

    private:

	//Decode map information: 
        std::vector<mpdmap_t>    fMPDmap;
        std::vector<Int_t>       fChanMapData;

        Double_t fZeroSuppressRMS;
        Bool_t fZeroSuppress;

	//Note on ROOT data types:
	// Char_t/UChar_t is one byte signed unsigned = 8 bits = 0..255 for unsigned
	// Short_t/UShort_t is two bytes signed/unsigned = 16 bits = 0..65535 for unsigned
	// Int_t/UInt_t is four bytes signed/unsigned = 32 bits
	// Long_t/ULong_t is eight bytes = 64 bit
	// Float_t/Double_t is four bytes/eight bytes

	UChar_t fN_APV25_CHAN;     //number of APV25 channels, default 128
	UChar_t fN_MPD_TIME_SAMP;  //number of MPD time samples, default = 6
	UShort_t fMPDMAP_ROW_SIZE; //MPDMAP_ROW_SIZE: default = 8, let's not hardcode
	
	
	//BASIC DECODED STRIP HIT INFO:
	UInt_t fNstrips_hit; //total Number of strips fired (after common-mode subtraction and zero suppression)
	
	std::vector<UShort_t> fStrip;  //Strip index of hit (these could be "U" or "V" generalized X and Y), assumed to run from 0..N-1
	std::vector<Bool_t>  fAxis;    //The use of Bool_t here reflects the (probably safe) assumption that there will never be more than 2 strip orientations present in one GEM module. 
	std::vector<std::vector<Int_t> > fADCsamples; //2D array of ADC samples by hit: Outer index runs over hits; inner index runs over ADC samples

	//Basic decoded strip info (analogous to "moduledata_t" structure from standalone code) for input to clustering algorithm:
	std::set<Int_t> fUstripList;  //List of unique "U" strips fired;
	std::set<Int_t> fVstripList;  //List of unique "V" strips fired;
	std::map<Int_t, Float_t> fADCsum_Ustrips; //Key = U strip index, value = sum of all ADC samples
	std::map<Int_t, Float_t> fADCsum_Vstrips; //Key = V strip index, value = sum of all ADC samples
	//Since we allow re-use of strips in multiple clusters, these arrays are no longer relevant
	/* std::map<Int_t, Bool_t> fUStripInCluster; //Key = U strip index, value = strip in cluster? */
	/* std::map<Int_t, Bool_t> fVStripInCluster; //Key = V strip index, value = strip in cluster? */
	/* std::map<Int_t, UInt_t> fUStripClusterIdx; //Key = U strip index, value = index in cluster/hit array */
	/* std::map<Int_t, UInt_t> fVStripClusterIdx; //Key = V strip index, value = inded in cluster/hit array */
	std::map<Int_t, std::vector<Float_t> > fADCsamp_Ustrips; //Key = U strip index, value = vector of ADC samples (converted to float)
	std::map<Int_t, std::vector<Float_t> > fADCsamp_Vstrips; //Key = V strip index, value = vector of ADC samples (converted to float)
	/* std::map<Int_t, Int_t> fBestMatch_Ustrips; //Key = U strip index, value = Best matching V strip index (by correlation coefficient) */
	/* std::map<Int_t, Int_t> fBestMatch_Vstrips; //Key = V strip index, value = Best matching U strip index (by correlation coefficient) */
	/* std::map<Int_t, Float_t> fBestCor_Ustrips; //Key = U strip index, value = Best correlation coefficient with best matching V strip */
	/* std::map<Int_t, Float_t> fBestCor_Vstrips; //Key = V strip index, value = Best correlation coefficient with best matching U strip */

	/* std::set<Int_t> fUstripList_filtered; //"Filtered" U strip list, based on whatever criteria we want */
	/* set::set<Int_t> fVstripList_filtered; //"Filtered" V strip list, based on whatever criteria we want */

	std::map<Int_t,Float_t> fADCmax_Ustrips; //Key = U strip index, value = maximum ADC sample on this strip
	std::map<Int_t,Float_t> fADCmax_Vstrips; //Key = V strip index, value = maximum ADC sample on this strip
	std::map<Int_t,Int_t>   fIsampmax_Ustrips; //Key = U strip index, value = time sample where maximum ADC occurs on this strip;
	std::map<Int_t,Int_t>   fIsampmax_Vstrips; //Key = V strip index, value = time sample where max. ADC occurs on this strip;

	std::map<Int_t,Float_t> fTmean_Ustrips;    //Key = U strip idx, value = ADC-weighted mean strip time (may add trigger jitter correction later on)
	std::map<Int_t,Float_t> fTmean_Vstrips;    //Key = V strip idx, value = ADC-weighted mean strip time (may add trigger jitter correction later on)

	std::map<Int_t,Float_t> fTsigma_Ustrips;    //Key = U strip idx, value = RMS of ADC-weighted mean strip time (may add trigger jitter correction later on)
	std::map<Int_t,Float_t> fTsigma_Vstrips;    //Key = V strip idx, value = RMS of ADC-weighted mean strip time (may add trigger jitter correction later on)

	//Strip timing with any relevant "corrections" applied (such as trigger time):
	std::map<Int_t,Float_t> fTmeanCorr_Ustrips;    //Key = U strip idx, value = ADC-weighted mean strip time (may add trigger jitter correction later on)
	std::map<Int_t,Float_t> fTmeanCorr_Vstrips;    //Key = V strip idx, value = ADC-weighted mean strip time (may add trigger jitter correction later on)

	std::map<Int_t,Float_t> fTsigmaCorr_Ustrips;    //Key = U strip idx, value = RMS of ADC-weighted mean strip time (may add trigger jitter correction later on)
	std::map<Int_t,Float_t> fTsigmaCorr_Vstrips;    //Key = V strip idx, value = RMS of ADC-weighted mean strip time (may add trigger jitter correction later on)
	////// (2D) Clustering results (see above for definition of struct sbsgemcluster_t:
	std::vector<sbsgemcluster_t> fClusters; 
	
	//Constant, module-specific parameters:
	UShort_t fModule; // Module index within a tracker. Since this is a GEM module class, this parameter should be unchanging
	UShort_t fLayer;  // Layer index of this module within a tracker. Since this is a GEM module class, this parameter should be unchanging

	UInt_t fNstripsU; // Total number of strips in this module along the generic "U" axis
	UInt_t fNstripsV; // Total number of strips in this module along the generic "V" axis
	
	//UInt_t fNTimeSamplesADC; //Number of ADC time samples (this could be variable in principle, but should be the same for all strips within a module within a run) redundant with fN_MPD_TIME_SAMP

	
	
	//GEOMETRICAL PARAMETERS:
	Double_t fUStripPitch;    //strip pitch along U, will virtually always be 0.4 mm
	Double_t fVStripPitch;    //strip pitch along V, will virtually always be 0.4 mm
	Double_t fUAngle;         //Angle between U strips and "X" axis of TRANSPORT coordinates;
	Double_t fVangle;         //Angle between V strips and "X" axis of TRANSPORT coordinates;
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

	//We will use 
       
	
	
        ClassDef(SBSGEMModule,0)

};




