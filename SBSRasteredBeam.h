#ifndef SBSRasteredBeam_H 
#define SBSRasteredBeam_H

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// SBSRasteredBeam                                                           //
//   Apparatus describing a rastered beam                                     //
//   Adapted from Bodo Reitz (April 2003)                                     //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "THaBeam.h"
#include "VarDef.h"

#include "TMath.h"
#include "TDatime.h"
#include "TList.h"

#include "SBSRaster.h"
#include "SBSBPM.h"

#include <deque>

class SBSRasteredBeam : public THaBeam {

 public:
  SBSRasteredBeam( const char* name, const char* description ) ;
  virtual ~SBSRasteredBeam() {}
  virtual Int_t  ReadDatabase( const TDatime& date );
  virtual Int_t Reconstruct() ;
  virtual Int_t CoarseReconstruct() ;

  TVector3 GetPositionAvg() const {return fBPM_rollingavg;}
  void SetBeamPosition(double xbeam, double ybeam, double zbeam){
    fbeam_position.SetXYZ(xbeam,ybeam,zbeam);
  }
  TVector3 GetBeamPosition() const {return fbeam_position;}

  double fBPM_L;
  double fBPMA_tg;
  double fRasterx_cen, fRastery_cen;
  double fRasterx_scale, fRastery_scale;
  Int_t fRaster_flag;
  TVector3 fbeam_position;

  void UpdateRollingAverage(TVector3 BPM, std::deque<TVector3> &BPM_container, Int_t &Nevents);

  std::deque<TVector3> fBPM_container_rollingavg;
  TVector3 fBPM_rollingavg;
  Int_t fNevents_rollingavg;
  Int_t fNevents_rollingmax;

protected:
  ClassDef(SBSRasteredBeam,0)    // A beam with rastered beam, analyzed event by event using raster currents

};

#endif

