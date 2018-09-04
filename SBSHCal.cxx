///////////////////////////////////////////////////////////////////////////////
//
// SBSHCal
//
///////////////////////////////////////////////////////////////////////////////
#include "SBSHCal.h"

ClassImp(SBSHCal);

/*
 * SBSHCal constructor.
 *
 * Specify SBSCalorimeter to use both TDC and ADC Multi-samples
 */
SBSHCal::SBSHCal( const char* name, const char* description,
    THaApparatus* apparatus ) : SBSCalorimeter(name,description,apparatus)
{
  SetWithADCSamples(true);
  SetWithTDC(true);
}

/*
 * Generic SBSHCal destructor
 */
SBSHCal::~SBSHCal()
{
}
