//////////////////////////////////////////////////////////////////////////
//
// Bare-bones SBB BigBite spectrometer class 
// Used for testing.
//
//////////////////////////////////////////////////////////////////////////

#include "SBSBigBite.h"

using namespace std;

ClassImp(SBSBigBite)


//_____________________________________________________________________________
SBSBigBite::SBSBigBite( const char* name, const char* description ) :
  THaSpectrometer( name, description )
{
  // Constructor. Defines standard detectors

}

//_____________________________________________________________________________
SBSBigBite::~SBSBigBite()
{
  // Destructor
}

//_____________________________________________________________________________
  Int_t SBSBigBite::FindVertices( TClonesArray& /* tracks */ )
{
  // Reconstruct target coordinates for all tracks found.

  // TODO

  return 0;
}

//_____________________________________________________________________________
Int_t SBSBigBite::TrackCalc()
{
  // Additioal track calculations

  // TODO

  return 0;
}

//_____________________________________________________________________________

