#ifndef SBSGEMTRACKERBASE_H
#define SBSGEMTRACKERBASE_H

#include <vector>
#include <map>
#include <set>
#include "SBSGEMModule.h"
#include "TVector3.h"
//#include <THaTrackingDetector.h>

//class THaRunBase;
//class THaApparatus;
//class THaEvData;
class SBSGEMModule;
//class THaCrateMap;

//This class is not going to inherit from THaAnything or from TObject.
//Instead, this class is only going to contain the common data members and methods needed by SBSGEMSpectrometerTracker and SBSGEMPolarimeterTracker, largely following the stand-alone clustering and track finding codes. The database reading and initialization will be taken care of by the derived classes: 



class SBSGEMTrackerBase {
public:
  
  ;
protected:
SBSGEMTrackerBase(); //only derived classes can construct me.
  std::vector <SBSGEMModule *> fModules; //array of SBSGEMModules:
  
  bool fIsMC;

  int fNmodules; //Total number of modules
  int fNlayers;  //total number of tracking layers

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
  double fGridBinWidthX, fGridBinWidthY; //Default bin size = 10 mm for both.
  int fGridNbinsX, fGridNbinsY;          //In the standalone code, these are typically derived from the grid bin size and the layer active area dimensions.

  double fTrackChi2Cut; //chi2/NDF cut for track validity

  // "Constraint points" to restrict the search region for track-finding:
  TVector3 fConstraintPoint_Front;
  TVector3 fConstraintPoint_Back;
  
};

#endif
