#ifndef SBSHCal_h
#define SBSHCal_h

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// SBSHCal                                                                   //
//                                                                           //
// A sub-class of SBSCalorimeter that has Multi-valued ADC + TDC             //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "SBSCalorimeter.h"

class SBSHCal : protected SBSCalorimeter {
public:
  SBSHCal( const char* name, const char* description = "",
      THaApparatus* a = NULL);
  virtual ~SBSHCal();


  ClassDef(SBSHCal,0)     // HCal detector class
};
#endif // SBSHCal_h
