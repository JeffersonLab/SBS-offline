#include "SBSRasteredBeam.h"
//_____________________________________________________________________________
SBSRasteredBeam::SBSRasteredBeam( const char* name, const char* description ) :
    THaBeam( name, description ) 
{
  AddDetector( new SBSRaster("Raster2","downstream raster") );
  AddDetector( new SBSRaster("Raster","upstream raster") );
  AddDetector( new SBSBPM("BPMA","1st BPM") );
  AddDetector( new SBSBPM("BPMB","2nd BPM") );

  fBPM_L = 5.15;
  fBPMA_tg = 7.53;
  fNevents_rollingmax = 500;
  fRasterx_cen = 40000;
  fRastery_cen = 40000;
  fRasterx_scale = 1.6e-4;
  fRastery_scale = 1.6e-4;

}
//_____________________________________________________________________________
Int_t SBSRasteredBeam::ReadDatabase( const TDatime& date )
{

  const char* const here = "ReadDatabase";

  vector<Double_t> Rasterx_range, Rastery_range;
  Double_t Raster_size;

  FILE* file = OpenFile( date );
  
  if(file){
    DBRequest calib_request[] = {
      { "BPM_separation", &fBPM_L, kDouble, 0, 1, 0 },
      { "BPMA_dist_to_targ", &fBPMA_tg, kDouble, 0, 1, 0 },
      { "Nevents_rollingmax"  ,   &fNevents_rollingmax,   kInt, 0, 1, 1 },
      { "Raster_curr_x_range", &Rasterx_range, kDoubleV, 0, 1, 1},
      { "Raster_curr_y_range", &Rastery_range, kDoubleV, 0, 1, 1},
      { "Raster_size", &Raster_size, kDouble, 0, 1, 1 },
      { 0 }
    };
    LoadDB( file, date, calib_request, fPrefix );

    fclose(file);

    fRasterx_cen = (Rasterx_range[1] + Rasterx_range[0]) / 2;
    fRastery_cen = (Rastery_range[1] + Rastery_range[0]) / 2;
  
    fRasterx_scale = Raster_size*1.0/1000 / (Rasterx_range[1] - Rasterx_range[0]);
    fRastery_scale = Raster_size/1000 / (Rastery_range[1] - Rastery_range[0]);
  }
  
  fNevents_rollingavg = 0;
  fBPM_container_rollingavg.resize(fNevents_rollingmax);
 
  
  return 0;
}
//_____________________________________________________________________________
Int_t SBSRasteredBeam::Reconstruct()
{

  return 0;
}
//_____________________________________________________________________________
Int_t SBSRasteredBeam::CoarseReconstruct()
{
  
  TIter nextDet( fDetectors ); 

  nextDet.Reset();

  // This apparatus assumes that there is only one detector 
  // in the list. If someone adds detectors by hand, the first 
  // detector in the list will be used to get the beam position
  // the others will be processed

  // This is the target position traditionally

  TVector3 BPMA;
  TVector3 BPMB;
  TVector3 Raster;

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

    if(theBeamDet->GetName() == std::string("Raster")){
      SBSRaster* RasterDet = reinterpret_cast<SBSRaster*>(theBeamDet);
      Raster.SetXYZ(RasterDet->GetRawPosX(),RasterDet->GetRawPosY(),0); 
    }

    if(theBeamDet->GetName() == std::string("BPMA"))
      BPMA = theBeamDet->GetPosition();
      
    if(theBeamDet->GetName() == std::string("BPMB"))
      BPMB = theBeamDet->GetPosition();
    
  }
  
  Update();

  TVector3 BPM_at_tgt = (BPMB-BPMA) * (fBPMA_tg/fBPM_L) + BPMA;

  UpdateRollingAverage(BPM_at_tgt, fBPM_container_rollingavg, fNevents_rollingavg);

  double xbeam = GetPositionAvg().X() + (Raster.X() - fRasterx_cen)*fRasterx_scale;
  double ybeam = GetPositionAvg().Y() + (Raster.Y() - fRastery_cen)*fRastery_scale;
  
  SetBeamPosition(xbeam,ybeam,0);
  
  return 0;

}
//______________________________________________________________________________
void SBSRasteredBeam::UpdateRollingAverage(TVector3 BPM, deque<TVector3> &BPM_container, Int_t &Nevents){
    
  TVector3 sum;
  
  if(Nevents < fNevents_rollingmax){
    BPM_container[Nevents] = BPM;

    if(Nevents == 0){     
      fBPM_rollingavg = BPM;
      sum = BPM;
    } else {

      TVector3 oldavg = fBPM_rollingavg;
      
      sum = Nevents * oldavg + BPM;
      
      fBPM_rollingavg = sum * (1.0/(Nevents + 1));
     
    }
    Nevents++;
  } else {
    
    TVector3 oldfirst = BPM_container.front();
    BPM_container.pop_front();
    BPM_container.push_back(BPM);
    
    TVector3 oldavg = fBPM_rollingavg;
    TVector3 oldsum = fBPM_rollingavg * fNevents_rollingmax;
    TVector3 newsum = oldsum - oldfirst + BPM;

    fBPM_rollingavg = newsum *(1.0 / fNevents_rollingmax);
  }
   
}
//____________________________________________________________________
ClassImp(SBSRasteredBeam)
