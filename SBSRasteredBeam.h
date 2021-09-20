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

class SBSRasteredBeam : public THaBeam {

public:
  SBSRasteredBeam( const char* name, const char* description ) ;
  virtual ~SBSRasteredBeam() {}
  virtual Int_t Reconstruct() ;

protected:
  ClassDef(SBSRasteredBeam,0)    // A beam with rastered beam, analyzed event by event using raster currents

};

#endif

