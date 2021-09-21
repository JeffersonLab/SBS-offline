#include "SBSRasteredBeam.h"
//_____________________________________________________________________________
SBSRasteredBeam::SBSRasteredBeam( const char* name, const char* description ) :
    THaBeam( name, description ) 
{
  AddDetector( new SBSRaster("Raster2","downstream raster") );
  AddDetector( new SBSRaster("Raster","upstream raster") );
  AddDetector( new SBSBPM("BPMA","1st BPM") );
  AddDetector( new SBSBPM("BPMB","2nd BPM") );

}
//_____________________________________________________________________________
Int_t SBSRasteredBeam::Reconstruct()
{

  TIter nextDet( fDetectors ); 

  nextDet.Reset();

  // This apparatus assumes that there is only one detector 
  // in the list. If someone adds detectors by hand, the first 
  // detector in the list will be used to get the beam position
  // the others will be processed

  // This is the target position traditionally

  if (THaBeamDet* theBeamDet=
      static_cast<THaBeamDet*>( nextDet() )) {
    theBeamDet->Process();
    fPosition = theBeamDet->GetPosition();
    fDirection = theBeamDet->GetDirection();
  }
  else {
    Error( Here("Reconstruct"), 
	   "Beamline Detectors Missing in Detector List" );
  }


  // Process any other detectors that may have been added (by default none)
  while (THaBeamDet * theBeamDet=
	 static_cast<THaBeamDet*>( nextDet() )) {
    theBeamDet->Process();
  }

  Update();

  return 0;

}
//______________________________________________________________________________
ClassImp(SBSRasteredBeam)
