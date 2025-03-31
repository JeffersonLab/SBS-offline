#ifndef SBSGEMPOLARIMETERTRACKER_H
#define SBSGEMPOLARIMETERTRACKER_H 1
#include <vector>
#include <THaTrackingDetector.h>
#include "SBSGEMTrackerBase.h"

class THaRunBase;
class THaApparatus;
class THaEvData;
class SBSGEMPlane;
class THaCrateMap;
class THaTrack;
class TClonesArray;
//class THaSpectrometer;

class SBSGEMPolarimeterTracker : public THaNonTrackingDetector, public SBSGEMTrackerBase {
 public:
  explicit SBSGEMPolarimeterTracker( const char *name, const char *description = "",
				     THaApparatus *app = nullptr );

  virtual ~SBSGEMPolarimeterTracker();

  virtual void    Clear( Option_t* opt="" );
  virtual Int_t   Decode( const THaEvData& );
  virtual EStatus Init( const TDatime& date );

  virtual Int_t   ReadDatabase( const TDatime& date );
  // We're going to need to override ReadGeometry for the GEM tracker classes since our definition of the module orientation angles differs from the standard definition
  // in THaDetectorBase:
  //virtual Int_t   ReadGeometry( FILE *file, const TDatime &date, Bool_t required = false );

  virtual Int_t   CoarseProcess( TClonesArray& tracks );
  virtual Int_t   FineProcess( TClonesArray& tracks );
  virtual Int_t   DefineVariables( EMode mode = kDefine );
  virtual void    Print(const Option_t* opt) const;
  virtual void    SetDebug( Int_t level );

  virtual Int_t   Begin( THaRunBase* r=0 );
  virtual Int_t   End( THaRunBase* r=0 );

  //  virtual bool PassedOpticsConstraint( TVector3 track_origin, TVector3 track_direction, bool coarsecheck=false );

  //Loop on all found tracks and calculate sclose, zclose, theta, phi:
  void CalcScatteringParameters();

  void SetFrontTrack( TVector3 track_origin, TVector3 track_direction );
  void SetFrontTrack( double x, double y, double theta, double phi );

  bool HasFrontTrack() const { return fFrontTrackIsSet; };
  
 private:
  // std::vector <SBSGEMModule *> fPlanes; storing the modules moved to SBSGEMTrackerBase

  //Not sure if we'll need this "test track" array for the polarimeter context.
  //bool fFrontTrackInitialized; 
  //bool fUseFrontTrackConstraint;
  
  //"Front tracks" array for polarimeter tracking:
  //TClonesArray *fFrontTracks; 
  //bool fIsMC; moved to SBSGEMTrackerBase

  //Add additional track properties we want to store:

  
  
  //For the time being, we will only consider one possible "front track":
  double fFrontTrackX;
  double fFrontTrackY;
  double fFrontTrackXp;
  double fFrontTrackYp;

  bool fFrontTrackIsSet;
  
  std::vector<double> fTrackTheta; //Polar scattering angle relative to front track
  std::vector<double> fTrackPhi; //Azimuthal scattering angle relative to front track
  std::vector<double> fTrackSClose; //distance of closest approach relative to front track
  std::vector<double> fTrackZClose; //Z of point of closest approach relative to front track.
  
  //THaCrateMap *fCrateMap; //Does this do anything? Not as far as I can tell. I wish someone would have commented about why they added this. AJRP
  ClassDef(SBSGEMPolarimeterTracker, 0);

};


#endif
