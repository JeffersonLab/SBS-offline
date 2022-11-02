#ifndef SBSGEMSPECTROMETERTRACKER_H
#define SBSGEMSPECTROMETERTRACKER_H 1
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

class SBSGEMSpectrometerTracker : public THaTrackingDetector, public SBSGEMTrackerBase {
 public:
  explicit SBSGEMSpectrometerTracker( const char *name, const char *description = "",
                                      THaApparatus *app = nullptr );

  virtual ~SBSGEMSpectrometerTracker();

  virtual void    Clear( Option_t* opt="" );
  virtual Int_t   Decode( const THaEvData& );
  virtual EStatus Init( const TDatime& date );

  virtual Int_t   ReadDatabase( const TDatime& date );
  // We're going to need to override ReadGeometry for the GEM tracker classes since our definition of the module orientation angles differs from the standard definition
  // in THaDetectorBase:
  //virtual Int_t   ReadGeometry( FILE *file, const TDatime &date, Bool_t required = false );

  virtual Int_t   CoarseTrack( TClonesArray& tracks );
  virtual Int_t   FineTrack( TClonesArray& tracks );
  virtual Int_t   DefineVariables( EMode mode = kDefine );
  virtual void    Print(const Option_t* opt) const;
  virtual void    SetDebug( Int_t level );

  virtual Int_t   Begin( THaRunBase* r=0 );
  virtual Int_t   End( THaRunBase* r=0 );

  virtual bool PassedOpticsConstraint( TVector3 track_origin, TVector3 track_direction, bool coarsecheck=false );

  
  
 private:
  // std::vector <SBSGEMModule *> fPlanes; storing the modules moved to SBSGEMTrackerBase

  bool fTestTrackInitialized;
  
  TClonesArray *fTestTracks; 
  //bool fIsMC; moved to SBSGEMTrackerBase
	
  //THaCrateMap *fCrateMap; //Does this do anything? Not as far as I can tell. I wish someone would have commented about why they added this. AJRP
  ClassDef(SBSGEMSpectrometerTracker, 0);

};


#endif
