//////////////////////////////////////////////////////////////////////////
//
// Bare-bones SBB BigBite spectrometer class 
// Used for testing.
//
//////////////////////////////////////////////////////////////////////////

#include "SBSEArm.h"

using namespace std;

ClassImp(SBSEArm)


//_____________________________________________________________________________
SBSEArm::SBSEArm( const char* name, const char* description ) :
  THaSpectrometer( name, description )
{
  // Constructor. Defines standard detectors

}

//_____________________________________________________________________________
SBSEArm::~SBSEArm()
{
  // Destructor
}

//_____________________________________________________________________________
  Int_t SBSEArm::FindVertices( TClonesArray& /* tracks */ )
{
  // Reconstruct target coordinates for all tracks found.

  // TODO

  return 0;
}

//_____________________________________________________________________________
Int_t SBSEArm::TrackCalc()
{
  // Additioal track calculations

  // TODO

  return 0;
}

//_____________________________________________________________________________

