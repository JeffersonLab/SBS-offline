#ifndef SBSGEMTRACKERBASE_H
#define SBSGEMTRACKERBASE_H

#include <vector>
#include <map>
#include <set>
#include "SBSGEMModule.h"
#include "TVector3.h"
#include "TVector2.h"
//#include <THaTrackingDetector.h>

//class THaRunBase;
//class THaApparatus;
//class THaEvData;
class SBSGEMModule;
//class THaCrateMap;

//This class is not going to inherit from THaAnything or from TObject.
//Instead, this class is only going to contain the common data members and methods needed by SBSGEMSpectrometerTracker and SBSGEMPolarimeterTracker, largely following the stand-alone clustering and track finding codes. The database reading and initialization will be taken care of by the derived classes: 
//Base class for GEM tracking assembly (of either the "tracking" or "non-tracking" flavor)


class SBSGEMTrackerBase {
public:
  void Clear(); //clear out all the event-specific data structures
  
protected:
  SBSGEMTrackerBase(); //only derived classes can construct me.
  virtual ~SBSGEMTrackerBase(); 

  bool fclustering_done;
  bool ftracking_done;
  
  //1D and 2D clustering: 
  void hit_reconstruction();
  //track-finding: 
  void find_tracks();
  
  //Utility methods: initialization:
  void InitLayerCombos();
  void InitGridBins(); //initialize 

  void InitHitList(); //Initialize (unchanging) "hit list" arrays used by track-finding: this only happens at the beginning of tracking
  void InitFreeHitList(); //Initialize "free hit list" arrays used on each track-finding iteration

  //Retrieve the global position of a hit by module and hit index:
  TVector3 GetHitPosGlobal( int modidx, int clustidx );

  //Calculate the intersection of a track with the plane of a module's active area:
  TVector3 TrackIntersect( int module, TVector3 track_origin, TVector3 track_direction );

  //Calculates the position of the track's intersection with a module in "module" coordinates U/V (the ones measured by the strips)
  TVector2 GetUVTrack( int module, TVector3 track_origin, TVector3 track_direction );
  
  //Utility method to iterate over combinations of hits in layers, used by find_tracks()
  bool GetNextCombo( const std::set<int> &layers, const std::map<int,std::vector<int> > &hitlist, std::map<int,int> &hitcounter, std::map<int,int> &hitcombo, bool &firstcombo=false );

  // Utility method to take a list of hits mapped by layer as input, and give track parameters and chi2 as output.
  // This relies on the "hit list" and "free hit list" information also being sensibly populated
  //FitTrack calculates chi2 and residuals as well as best fit parameters
  void FitTrack( const std::map<int,int> &hitcombo, double &xtrack, double &ytrack, double &xptrack, double &yptrack, double &chi2ndf, vector<double> &uresid, vector<double> &vresid );
  //CalcLineOfBestFit only calculates the track parameters, does not calculate chi2 or residuals:
  void CalcLineOfBestFit( const std::map<int,int> &hitcombo, double &xtrack, double &ytrack, double &xptrack, double &yptrack );

  // routine to fit the best track to a set of hits, without the overhead of chi2 calculation, useful for "exclusive residuals" calculation:
  //void FitTrackNoChisquaredCalc( const std::map<int,int> &hitcombo, double &xtrack, double &ytrack, double &xptrack, double &yptrack );
  
  // Method to add a new Track to the track arrays: this takes the best hit combination and the parameters of the line of best fit to those hits
  // and the (already calculated) chi2 and fills the tracking results arrays: best fit parameters, inclusive and exclusive tracking residuals, and hit lists by track:
  void AddTrack( const std::map<int,int> &hitcombo, double *BestTrack, double chi2ndf, vector<double> &uresid, vector<double> &vresid );
  
  //Data members:
  std::vector <SBSGEMModule *> fModules; //array of SBSGEMModules:

  bool fOnlinePedestalSubtraction; //Flag specifying whether pedestal subtraction has been done "online" (maybe this should be module-specific? probably not)
  
  bool fIsMC;

  int fNmodules; //Total number of modules
  int fNlayers;  //total number of tracking layers

  int fTrackingAlgorithmFlag; //Choose track algorithm

  int fMinHitsOnTrack; //default = 3; cannot be less than 3, cannot be more than total number of layers
  
  long fMaxHitCombinations; //default = 100000; skip 
  
  // The use of maps here instead of vectors may be slightly algorithmically inefficient, but it DOES guarantee that the maps are
  //  (a) sorted by increasing layer index, which, generally speaking, for a sensibly constructed database, will also be in ascending order of Z.
  //  (b) each unique logical tracking layer index can only occur exactly once
  std::map<int,int> fNumModulesByLayer; //key = unique layer ID (logical tracking layer), mapped value = number of modules per layer
  std::map<int, std::set<int> > fModuleListByLayer;  //key = unique layer ID, mapped value = list of unique modules associated with this layer

  std::map<int, std::vector<std::vector<int> > > fLayerCombinations; //key = minimum hit requirement to form a track. Mapped value = 2D array of layer combinations at a given minimum hit requirement.

  std::map<int, double> fZavgLayer; //Average z position by layer. This is MAINLY used by the efficiency determination and plotting machinery, not so much by the reconstruction algorithms.

  //"Grid bins" for fast track-finding algorithm(s): define limits of layer active area:
  std::map<int, double> fXmin_layer, fXmax_layer, fYmin_layer, fYmax_layer;
  //Grid bin size: smaller values should give faster track finding, but bins should be large compared to GEM spatial resolution.
  double fGridBinWidthX, fGridBinWidthY; //Default bin size = 10 mm for both. we are using same grid bin width at all layers:
  double fGridEdgeToleranceX, fGridEdgeToleranceY; 
  std::map<int, int> fGridNbinsX_layer, fGridNbinsY_layer;          //In the standalone code, these are typically derived from the grid bin size and the layer active area dimensions.
  //These variables are arguably redundant with the ones above, but as defined, these include a bit of extra "slop" to account for resolution, misalignments, etc.
  std::map<int, double> fGridXmin_layer, fGridYmin_layer, fGridXmax_layer, fGrid_Ymax_layer;
  
  double fTrackChi2Cut; //chi2/NDF cut for track validity

  bool fUseConstraint;
  // "Constraint points" to restrict the search region for track-finding:
  TVector3 fConstraintPoint_Front;
  TVector3 fConstraintPoint_Back;

  TVector2 fConstraintWidth_Front;
  TVector2 fConstraintWidth_Back;

  Double_t fSigma_hitpos;   //sigma parameter controlling resolution entering track chi^2 calculation
  Double_t fSigma_hitshape; //Sigma parameter controlling hit shape for cluster-splitting algorithm.
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //            DATA members to hold the track information (at least temporarily, will eventually                 //
  //            be passed to the THaSpectrometer tracks TClonesArray for SBSGEMSpectrometerTracker                //
  //            We'll need to figure out how to handle things for SBSGEMPolarimeterTracker                        //
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  //////////////////// "Hit list" arrays used by tracking algorithm /////////////////////////

  /////////// unchanging arrays that are filled at the beginning of track-finding: ///////////
  //Define and initialize the "hit lists" and other arrays that we will need to do the tracking:
  //These "hit list" arrays will NOT be modified throughout the track-finding procedure, except at the beginning:
  std::set<int> layers_with_2Dhits; //list of tracking layers with at least one 2D hit reconstructed (within the track search region)
  std::map<int,int> N2Dhits_layer; // key = layer, mapped value = total number of reconstructed 2D hits
  std::map<int,std::vector<int> > modindexhit2D; //key = layer, mapped value = module index of hits in that layer:
  std::map<int,std::vector<int> > clustindexhit2D; //key = layer, mapped value = index of hits within 2D cluster array of module in question
  std::map<int,std::vector<bool> > hitused2D; //flag to tell each track-finding iteration whether hit was already used in a previous track

  //////////////////// "Free hit list" arrays used on individual track-finding iterations: /////////////////////////////
  std::map<int,int> Nfreehits_layer; //key = layer, mapped value = number of unused hits available:
  std::set<int> layerswithfreehits; //list of layers with at least one unused hit
  std::map<int,std::vector<int> > freehitlist_layer; //list of unused hits mapped by layer: index in the unchanging array defined above (clustindex2D)
  std::map<int,int> freehitcounter; //When using a "brute force" track-finding algorithm, this counter is used for looping over hit combinations ("odometer" algorithm)
  
  
  //Hit lists mapped by grid bin:
  std::map<int, std::vector<int> > Nfreehits_binxy_layer; //number of free hits by grid bin in each layer
  std::map<int, std::vector<std::vector<int> > > freehitlist_binxy_layer; //list of free hits by layer and 2D grid bin; again, the "hit list" contains the index in the unchanging array clustindex2D
  
  
  //////////////////// Tracking results: //////////////////////////////
  
  int fNtracks_found;
  std::vector<int> fNhitsOnTrack; //number of hits on track:
  std::vector<std::vector<int> > fModListTrack; //list of modules containing hits on fitted tracks
  std::vector<std::vector<int> > fHitListTrack; //list of hits on fitted tracks: NOTE--the "hit list" of the track refers to the index in the 2D cluster array. To locate the hit and its properties you need the module index and the hit index, i.e., fModules[fModListTrack[ihit]]->
  std::vector<std::vector<double> > fresidu_hits; //inclusive residuals: track - hit along direction measured by u strips
  std::vector<std::vector<double> > fresidv_hits; //inclusive residuals: track - hit along direction measured by v strips
  std::vector<std::vector<double> > feresidu_hits; //exclusive residuals: track - hit along direction measured by u strips
  std::vector<std::vector<double> > feresidv_hits; //exclusive residuals: track - hit along direction measured by v strips

  //Fitted track parameters: coordinates at Z = 0 and slopes
  std::vector<double> fXtrack;
  std::vector<double> fYtrack;
  std::vector<double> fXptrack;
  std::vector<double> fYptrack;
  std::vector<double> fChi2Track; //chi2/ndf
  
  
};

#endif
