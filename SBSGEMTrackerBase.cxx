#include "SBSGEMTrackerBase.h"
#include "SBSGEMModule.h"
#include "TRotation.h"

SBSGEMTrackerBase::SBSGEMTrackerBase(){ //Set default values of important parameters: 
  Clear();

  fIsMC = false;
  fNmodules = 0;
  fNlayers = 0;
  fTrackingAlgorithmFlag = 2;

  fMinHitsOnTrack = 3;

  fOnlinePedestalSubtraction = true;
  fZeroSuppress = true;
  fZeroSuppressRMS = 10.0; //10 ADC channels. We are free to define how this is actually used;

  fGridBinWidthX = 0.01; //1 cm = 10 mm;
  fGridBinWidthY = 0.01; //1 cm = 10 mm;

  fTrackChi2Cut = 100.0; //Max. chi2/ndf for a combination of hits to form a track

  // set defaults for constraint points and constraint widths:
  fConstraintPoint_Front.SetXYZ(0,0,0);
  fConstraintPoint_Back.SetXYZ(0,0,10.0);

  //wide-open constraints for now:
  fConstraintWidth_Front.Set( 1.5, 0.5 );
  fConstraintWidth_Back.Set( 1.5, 0.5 ); 
}

void SBSGEMTrackerBase::~SBSGEMTrackerBase(){
  //for now, do nothing; let the derived classes handle the clearing out of the modules
  
}

void SBSGEMTrackerBase::Clear(){ //Clear out any event-specific stuff
  fNtracks_found = 0;
  fNhitsOnTrack.clear();
  fModListTrack.clear();
  fHitListTrack.clear();
  fresidu_hits.clear();
  fresidv_hits.clear();
  feresidu_hits.clear();
  feresidv_hits.clear();

  fXtrack.clear();
  fYtrack.clear();
  fXptrack.clear();
  fYptrack.clear();
  fChi2Track.clear();
  
}
