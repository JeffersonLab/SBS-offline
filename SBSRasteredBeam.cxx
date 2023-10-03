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
  fRaster_flag = 0;
}
//_____________________________________________________________________________
Int_t SBSRasteredBeam::ReadDatabase( const TDatime& date )
{

  const char* const here = "ReadDatabase";

  vector<Double_t> Rasterx_range, Rastery_range, Raster2x_range, Raster2y_range;
  Double_t Raster_size;

  FILE* file = OpenFile( date );
  
  if(file){
    DBRequest calib_request[] = {
      { "BPM_separation", &fBPM_L, kDouble, 0, 1, 0 },
      { "BPMA_dist_to_targ", &fBPMA_tg, kDouble, 0, 1, 0 },
      { "Nevents_rollingmax"  ,   &fNevents_rollingmax,   kInt, 0, 1, 1 },
      { "Raster_curr_x_range", &Rasterx_range, kDoubleV, 0, 1, 1},
      { "Raster_curr_y_range", &Rastery_range, kDoubleV, 0, 1, 1},
      { "Raster2_curr_x_range", &Raster2x_range, kDoubleV, 0, 1, 1},
      { "Raster2_curr_y_range", &Raster2y_range, kDoubleV, 0, 1, 1},
      { "Raster_size", &Raster_size, kDouble, 0, 1, 1 },
      { "Raster_type", &fRaster_flag, kInt, 0, 1, 1 },
      { 0 }
    };
    LoadDB( file, date, calib_request, fPrefix );

    fclose(file);

    double Rasterx_bounds[2], Rastery_bounds[2];
    
    // If we are using 2 rasters we take the difference
    if(fRaster_flag == 2 && Rasterx_range.size() == 2 && Rastery_range.size() == 2 && Raster2x_range.size() == 2 && Raster2y_range.size() == 2){
      Rasterx_bounds[0] = Rasterx_range[0] - Raster2x_range[1];
      Rasterx_bounds[1] = Rasterx_range[1] - Raster2x_range[0];
      Rastery_bounds[0] = Rastery_range[0] - Raster2y_range[1];
      Rastery_bounds[1] = Rastery_range[1] - Raster2y_range[0];
      
      // Calculate the raster center
      fRasterx_cen = (Rasterx_bounds[1] + Rasterx_bounds[0]) / 2;
      fRastery_cen = (Rastery_bounds[1] + Rastery_bounds[0]) / 2;

      // Calculate the scale factor. 1000 factor is used to convert mm to m
      fRasterx_scale = Raster_size*1.0/1000 / (Rasterx_bounds[1] - Rasterx_bounds[0]);
      fRastery_scale = Raster_size*1.0/1000 / (Rastery_bounds[1] - Rastery_bounds[0]);
    } else if (fRaster_flag == 1 && Rasterx_range.size() == 2 && Rastery_range.size() == 2) { //In every other case we just use upstream raster
      Rasterx_bounds[0] = Rasterx_range[0];
      Rasterx_bounds[1] = Rasterx_range[1];
      Rastery_bounds[0] = Rastery_range[0];
      Rastery_bounds[1] = Rastery_range[1];

      // Calculate the raster center
      fRasterx_cen = (Rasterx_bounds[1] + Rasterx_bounds[0]) / 2;
      fRastery_cen = (Rastery_bounds[1] + Rastery_bounds[0]) / 2;

      // Calculate the scale factor. 1000 factor is used to convert mm to m
      fRasterx_scale = Raster_size*1.0/1000 / (Rasterx_bounds[1] - Rasterx_bounds[0]);
      fRastery_scale = Raster_size*1.0/1000 / (Rastery_bounds[1] - Rastery_bounds[0]);
    }
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

  TVector3 BPMA, BPMB;
  TVector3 Raster, Raster2;

  if (THaBeamDet* theBeamDet=
      static_cast<THaBeamDet*>( nextDet() )) {
    theBeamDet->Process();
    fPosition = theBeamDet->GetPosition();
    fDirection = theBeamDet->GetDirection();
    
    //Retrieve the Downstream raster X and Y currents
    if(theBeamDet->GetName() == std::string("Raster2")){
      SBSRaster* RasterDet = reinterpret_cast<SBSRaster*>(theBeamDet);
      Raster2.SetXYZ(RasterDet->GetRawPosX(),RasterDet->GetRawPosY(),0); 
    }
  }
  else {
    Error( Here("Reconstruct"), 
	   "Beamline Detectors Missing in Detector List" );
  }
  
  // Process any other detectors that may have been added (by default none)
  while (THaBeamDet * theBeamDet=
	 static_cast<THaBeamDet*>( nextDet() )) {
    theBeamDet->Process();

    // Retrieve the upstream raster X and Y currents
    if(theBeamDet->GetName() == std::string("Raster")){
      SBSRaster* RasterDet = reinterpret_cast<SBSRaster*>(theBeamDet);
      Raster.SetXYZ(RasterDet->GetRawPosX(),RasterDet->GetRawPosY(),0); 
    }

    // Get the BPM A and B positions
    if(theBeamDet->GetName() == std::string("BPMA"))
      BPMA = theBeamDet->GetPosition();
      
    if(theBeamDet->GetName() == std::string("BPMB"))
      BPMB = theBeamDet->GetPosition();
    
  }
  
  Update();

  // Project the BPM positions to the target
  TVector3 BPM_at_tgt = (BPMB-BPMA) * (fBPMA_tg/fBPM_L) + BPMA;

  // Update the rolling average of the target position
  UpdateRollingAverage(BPM_at_tgt, fBPM_container_rollingavg, fNevents_rollingavg);

  double xbeam, ybeam;

  ////////////////////////////////////////////////////////////////////////////////
  // Calculate the final beam position using the BPM projection to the target
  // and the raster value converted to mm
  // 
  // x/y = x/y_avg + (Raster_x/y - Raster_center_x/y) * Raster_scale_x/y
  // 
  // Flag = 2: Two rasters so Raster_x/y is the difference between the upstream
  //           and downstream raster.
  // Flag = 1: One raster so Raster_x/y is the upstream raster.
  // Flag = 0: Unrastered beam so the position is just teh BPM projection average
  ////////////////////////////////////////////////////////////////////////////////
  if(fRaster_flag == 2){
    //cout<<fRasterx_scale<<" "<<fRastery_scale<<endl;
    xbeam = GetPositionAvg().X() + (Raster.X() - Raster2.X() - fRasterx_cen)*fRasterx_scale;
    ybeam = GetPositionAvg().Y() + (Raster.Y() - Raster2.Y() - fRastery_cen)*fRastery_scale;
  } else if (fRaster_flag == 1) {
    xbeam = GetPositionAvg().X() + (Raster.X() - fRasterx_cen)*fRasterx_scale;
    ybeam = GetPositionAvg().Y() + (Raster.Y() - fRastery_cen)*fRastery_scale;
  } else {
    xbeam = GetPositionAvg().X();
    ybeam = GetPositionAvg().Y();
  }
  
  SetBeamPosition(xbeam,ybeam,0);
  
  return 0;

}
//______________________________________________________________________________
// This function takes the BPM value for one event and updates to average BPM
// We are defaulted to take a running average of 500 events
void SBSRasteredBeam::UpdateRollingAverage(TVector3 BPM, deque<TVector3> &BPM_container, Int_t &Nevents){
    
  TVector3 sum;
  
  // Less than max events we add everything to the list of BPM values
  if(Nevents < fNevents_rollingmax){
    BPM_container[Nevents] = BPM;

    if(Nevents == 0){ // Initialize the first event
      fBPM_rollingavg = BPM;
      sum = BPM;
    } else {

      TVector3 oldavg = fBPM_rollingavg;
      
      sum = Nevents * oldavg + BPM;
      
      fBPM_rollingavg = sum * (1.0/(Nevents + 1)); // Update the average
     
    }
    Nevents++;
  } else { // We have reached the maximum and do not want the list to get longer
    
    // Remove the oldest event info and add the current event to the end of the list
    TVector3 oldfirst = BPM_container.front();
    BPM_container.pop_front();
    BPM_container.push_back(BPM);
    
    TVector3 oldavg = fBPM_rollingavg;
    TVector3 oldsum = fBPM_rollingavg * fNevents_rollingmax;
    TVector3 newsum = oldsum - oldfirst + BPM;

    // Calculate the new average
    fBPM_rollingavg = newsum *(1.0 / fNevents_rollingmax);
  }
   
}
//____________________________________________________________________
ClassImp(SBSRasteredBeam)
