#include <iostream>

#include "SBSGEMModule.h"
#include "TDatime.h"
#include "THaEvData.h"

//using namespace SBSGEMModule;

//This should not be hard-coded, I think, but read in from the database (or perhaps not, if it never changes? For now we keep it hard-coded)
const int APVMAP[128] = {1, 33, 65, 97, 9, 41, 73, 105, 17, 49, 81, 113, 25, 57, 89, 121, 3, 35, 67, 99, 11, 43, 75, 107, 19, 51, 83, 115, 27, 59, 91, 123, 5, 37, 69, 101, 13, 45, 77, 109, 21, 53, 85, 117, 29, 61, 93, 125, 7, 39, 71, 103, 15, 47, 79, 111, 23, 55, 87, 119, 31, 63, 95, 127, 0, 32, 64, 96, 8, 40, 72, 104, 16, 48, 80, 112, 24, 56, 88, 120, 2, 34, 66, 98, 10, 42, 74, 106, 18, 50, 82, 114, 26, 58, 90, 122, 4, 36, 68, 100, 12, 44, 76, 108, 20, 52, 84, 116, 28, 60, 92, 124, 6, 38, 70, 102, 14, 46, 78, 110, 22, 54, 86, 118, 30, 62, 94, 126};


SBSGEMModule::SBSGEMModule( const char *name, const char *description,
			    THaDetectorBase* parent ):
  THaSubDetector(name,description,parent)
{
  // FIXME:  To database
  //Set Default values for fZeroSuppress and fZeroSuppressRMS:
  fZeroSuppress    = kTRUE;
  fZeroSuppressRMS = 5.0;

  //Default online zero suppression to FALSE: actually I wonder if it would be better to use this in 
  // Moved this to MPDModule, since this should be done during the decoding of the raw APV data:
  //fOnlineZeroSuppression = kFALSE;
  
  //Set default values for decode map parameters:
  fN_APV25_CHAN = 128;
  fN_MPD_TIME_SAMP = 6;
  fMPDMAP_ROW_SIZE = 9;

  fNumberOfChannelInFrame = 129;

  fSamplePeriod = 25.0; //nanoseconds:
  
  // for( Int_t i = 0; i < N_MPD_TIME_SAMP; i++ ){
  //     fadc[i] = NULL;
  // }

  return;
}

SBSGEMModule::~SBSGEMModule() {
  // if( fStrip ){
  //     fadc0 = NULL;
  //     fadc1 = NULL;
  //     fadc2 = NULL;
  //     fadc3 = NULL;
  //     fadc4 = NULL;
  //     fadc5 = NULL;

  //     for( Int_t i = 0; i < N_MPD_TIME_SAMP; i++ ){
  //         delete fadc[i];
  //         fadc[i] = NULL;
  //     }
  //     delete fPedestal;
  //     fPedestal = NULL;
  //     delete fStrip;
  //     fStrip = NULL;
  // }

  return;
}

Int_t SBSGEMModule::ReadDatabase( const TDatime& date ){
  std::cout << "[SBSGEMModule::ReadDatabase]" << std::endl;

  Int_t status;

  FILE* file = OpenFile( date );
  if( !file ) return kFileError;

  std::vector<Double_t> rawped;
  std::vector<Double_t> rawrms;

  //UShort_t layer;
  
  //I think we can set up the entire database parsing in one shot here (AJRP). Together with the call to ReadGeometry, this should define basically everything we need.
  //After loading, we will run some checks on the loaded information:
  //copying commented structure of DBRequest here for reference (see VarDef.h):
  // struct DBRequest {
  //   const char*      name;     // Key name
  //   void*            var;      // Pointer to data (default to Double*)
  //   VarType          type;     // (opt) data type (see VarType.h, default Double_t)
  //   UInt_t           nelem;    // (opt) number of array elements (0/1 = 1 or auto)
  //   Bool_t           optional; // (opt) If true, missing key is ok
  //   Int_t            search;   // (opt) Search for key along name tree
  //   const char*      descript; // (opt) Key description (if 0, same as name)
  // };
  
  const DBRequest request[] = {
    { "chanmap",        &fChanMapData,        kIntV, 0, 0, 0}, // mandatory: decode map info
    { "ped",            &rawped,        kDoubleV, 0, 1, 1}, // optional raw pedestal info (not clear how we want to use or organize yet)
    { "rms",            &rawrms,        kDoubleV, 0, 1, 1}, // optional pedestal rms info (not clear how we want to use or organize yet)
    { "layer",          &fLayer,         kUShort, 0, 0, 0}, // mandatory: logical tracking layer must be specified for every module:
    { "nstripsu",       &fNstripsU,     kUInt, 0, 0, 1}, //mandatory: number of strips in module along U axis
    { "nstripsv",       &fNstripsV,     kUint, 0, 0, 1}, //mandatory: number of strips in module along V axis
    { "uangle",         &fUAngle,       kDouble, 0, 0, 1}, //mandatory: Angle of "U" strips wrt X axis
    { "vangle",         &fVangle,       kDouble, 0, 0, 1}, //mandatory: Angle of "V" strips wrt X axis
    { "upitch",         &fUStripPitch,  kDouble, 0, 0, 1}, //mandatory: Pitch of U strips
    { "vpitch",         &fVStripPitch,  kDouble, 0, 0, 1}, //mandatory: Pitch of V strips
    { "ugain",          &fUgain,        kDoubleV, 0, 1, 0}, //(optional): Gain of U strips by APV card (ordered by strip position, NOT by order of appearance in decode map)
    { "vgain",          &fVgain,        kDoubleV, 0, 1, 0}, //(optional): Gain of V strips by APV card (ordered by strip position, NOT by order of appearance in decode map)
    { "threshold_sample",  &fThresholdSample, kDouble, 0, 1, 1}, //(optional): threshold on max. ADC sample to keep strip (baseline-subtracted)
    { "threshold_stripsum", &fThresholdStripSum, kDouble, 0, 1, 1}, //(optional): threshold on sum of ADC samples on a strip (baseline-subtracted)
    { "threshold_clustersum", &fThresholdClusterSum, kDouble, 0, 1, 1}, //(optional): threshold on sum of all ADCs over all strips in a cluster (baseline-subtracted)
    { "ADCasym_cut", &fADCasymCut, kDouble, 0, 1, 1}, //(optional): filter 2D hits by ADC asymmetry, |Asym| < cut
    { "deltat_cut", &fTimeCutUVdiff, kDouble, 0, 1, 1}, //(optional): filter 2D hits by U/V time difference
    { "peakprominence_minsigma", &fThresh_2ndMax_nsigma, kDouble, 0, 1, 1}, //(optional): reject overlapping clusters with peak prominence less than this number of sigmas
    { "peakprominence_minfraction", &fThresh_2ndMax_fraction, kDouble, 0, 1, 1}, //(optional): reject overlapping clusters with peak prominence less than this fraction of height of nearby higher peak
    { "maxnu_charge", &fMaxNeighborsU_totalcharge, kUShort, 0, 1, 1}, //(optional): cluster size restriction along U for total charge calculation
    { "maxnv_charge", &fMaxNeighbotsV_totalcharge, kUShort, 0, 1, 1}, //(optional): cluster size restriction along V for total charge calculation
    { "maxnu_pos", &fMaxNeighborsU_hitpos, kUShort, 0, 1, 1}, //(optional): cluster size restriction for position reconstruction
    { "maxnv_pos", &fMaxNeighborsV_hitpos, kUShort, 0, 1, 1}, //(optional): cluster size restriction for position reconstruction
    {0}
  };
  status = LoadDB( file, date, request, fPrefix, 1 ); //The "1" after fPrefix means search up the tree

  if( status != 0 ){
    fclose(file);
    return status;
  }
  
  Int_t nentry = fChanMapData.size()/fMPDMAP_ROW_SIZE;
  for( Int_t mapline = 0; mapline < nentry; mapline++ ){
    mpdmap_t thisdata;
    thisdata.crate  = fChanMapData[0+mapline*MPDMAP_ROW_SIZE];
    thisdata.slot   = fChanMapData[1+mapline*MPDMAP_ROW_SIZE];
    thisdata.mpd_id = fChanMapData[2+mapline*MPDMAP_ROW_SIZE];
    thisdata.gem_id = fChanMapData[3+mapline*MPDMAP_ROW_SIZE];
    thisdata.adc_id = fChanMapData[4+mapline*MPDMAP_ROW_SIZE];
    thisdata.i2c    = fChanMapData[5+mapline*MPDMAP_ROW_SIZE];
    thisdata.pos    = fChanMapData[6+mapline*MPDMAP_ROW_SIZE];
    thisdata.invert = fChanMapData[7+mapline*MPDMAP_ROW_SIZE];
    thisdata.axis   = fChanMapData[8+mapline*MPDMAP_ROW_SIZE];
    fMPDmap.push_back(thisdata);
  }

  std::cout << fName << " mapped to " << nentry << " APV25 chips" << std::endl;

  //Geometry info is required to be present in the database for each module:
  Int_t err = ReadGeometry( file, date, true );
  if( err ) {
    fclose(file);
    return err;
  }

  //Initialize all pedestals to zero, RMS values to default:
  fPedestalU.clear();
  fPedestalU.resize( fNstripsU ); 

  fPedRMSU.clear();
  fPedRMSU.resize( fNstripsU );
  
  for ( UInt_t istrip=0; istrip<fNstripsU; istrip++ ){
    fPedestalU[istrip] = 0.0;
    fPedRMSU[istrip] = fZeroSuppressRMS; //placeholder to be replaced by value from database
  }

  //Initialize all pedestals to zero, RMS values to default:
  fPedestalV.clear();
  fPedestalV.resize( fNstripsV ); 

  fPedRMSV.clear();
  fPedRMSV.resize( fNstripsV );
  
  for( UInt_t istrip=0; istrip<fNstripsV; istrip++ ){
    fPedestalV[istrip] = 0.0;
    fPedRMSV[istrip] = fZeroSuppressRMS;
  }

  // for( UInt_t i = 0; i < rawped.size(); i++ ){
  //   if( (i % 2) == 1 ) continue;
  //   int idx = (int) rawped[i];
	
  //   if( idx < N_APV25_CHAN*nentry ){
  //     fPedestal[idx] = rawped[i+1];
  //   } else {
		
  //     std::cout << "[SBSGEMModule::ReadDatabase]  WARNING: " << " strip " << idx  << " listed but not enough strips in cratemap" << std::endl;
  //   }
  // }

  // for( UInt_t i = 0; i < rawrms.size(); i++ ){
  //   if( (i % 2) == 1 ) continue;
  //   int idx = (int) rawrms[i];
  //   if( idx < N_APV25_CHAN*nentry ){
  //     fRMS[idx] = rawrms[i+1];
  //   } else {
  //     std::cout << "[SBSGEMModule::ReadDatabase]  WARNING: " << " strip " << idx  << " listed but not enough strips in cratemap" << std::endl;
  //   }
  // }

  fclose(file);
  
  return 0;
}

Int_t SBSGEMModule::ReadGeometry( FILE *file, const TDaTime &date, Bool_t required ){ //We start with a copy of THaDetectorBase::ReadGeometry and modify accordingly:
  // Read this detector's basic geometry information from the database.
  // Derived classes may override to read more advanced data.

  const char* const here = "ReadGeometry";

  vector<double> position, size, angles;
  Bool_t optional = !required;
  DBRequest request[] = {
    { "position", &position, kDoubleV, 0, optional, 0,
      "\"position\" (detector position [m])" },
    { "size",     &size,     kDoubleV, 0, optional, 0,
      "\"size\" (detector size [m])" },
    { "angle",    &angles,   kDoubleV, 0, true, 0,
      "\"angle\" (detector angles(s) [deg]" },
    { nullptr }
  };
  Int_t err = LoadDB( file, date, request );
  if( err )
    return kInitError;

  if( !position.empty() ) {
    if( position.size() != 3 ) {
      Error( Here(here), "Incorrect number of values = %u for "
	     "detector position. Must be exactly 3. Fix database.",
	     static_cast<unsigned int>(position.size()) );
      return 1;
    }
    fOrigin.SetXYZ( position[0], position[1], position[2] );
  }
  else
    fOrigin.SetXYZ(0,0,0);

  if( !size.empty() ) {
    if( size.size() != 3 ) {
      Error( Here(here), "Incorrect number of values = %u for "
	     "detector size. Must be exactly 3. Fix database.",
	     static_cast<unsigned int>(size.size()) );
      return 2;
    }
    if( size[0] == 0 || size[1] == 0 || size[2] == 0 ) {
      Error( Here(here), "Illegal zero detector dimension. Fix database." );
      return 3;
    }
    if( size[0] < 0 || size[1] < 0 || size[2] < 0 ) {
      Warning( Here(here), "Illegal negative value for detector dimension. "
	       "Taking absolute. Check database." );
    }
    fSize[0] = 0.5 * TMath::Abs(size[0]);
    fSize[1] = 0.5 * TMath::Abs(size[1]);
    fSize[2] = TMath::Abs(size[2]);
  }
  else
    fSize[0] = fSize[1] = fSize[2] = kBig;

  if( !angles.empty() ) {
    if( angles.size() != 1 && angles.size() != 3 ) {
      Error( Here(here), "Incorrect number of values = %u for "
	     "detector angle(s). Must be either 1 or 3. Fix database.",
	     static_cast<unsigned int>(angles.size()) );
      return 4;
    }
    // If one angle is given, it indicates a rotation about y, as before.
    // If three angles are given, they are interpreted as rotations about the X, Y, and Z axes, respectively:
    // 
    if( angles.size() == 1 ) {
      DefineAxes( angles[0] * TMath::DegToRad() );
    }
    else {
      TRotation RotTemp;

      // So let's review how to define the detector axes correctly.

      // THaDetectorBase::DetToTrackCoord(TVector3 p) returns returns p.X * fXax() + p.Y * fYax() + p.Z * fZax() + fOrigin
      // In the standalone code, we do Rot * (p) + fOrigin (essentially):
      // So in matrix form, when we do TRotation::RotateX(alpha), we get:

      // RotTemp * Point =  |  1    0            0         |    |  p.X()  |
      //                    |  0   cos(alpha) -sin(alpha)  | *  |  p.Y()  |
      //                    |  0   sin(alpha)  cos(alpha)  |    |  p.Z()  |
      // 
      // This definition ***appears**** to be consistent with the "sense" of the rotation as applied by the standalone code.
      // The detector axes are defined as the ROWS of the rotation matrix. We will have to test that it is working correctly, however:
   
      RotTemp.RotateX( angles[0] );
      RotTemp.RotateY( angles[1] );
      RotTemp.RotateZ( angles[2] );
      
      fXax.SetXYZ( RotTemp.XX(), RotTemp.XY(), RotTemp.XZ() );
      fYax.SetXYZ( RotTemp.YX(), RotTemp.YY(), RotTemp.YZ() );
      fZax.SetXYZ( RotTemp.ZX(), RotTemp.ZY(), RotTemp.ZZ() );
    }
  } else
    DefineAxes(0);

  return 0;
}

    
Int_t SBSGEMModule::DefineVariables( EMode mode ) {
  if( mode == kDefine and fIsSetup ) return kOK;
  fIsSetup = ( mode == kDefine );

  RVarDef vars[] = {
    { "nch",   "Number of channels",   "fNch" },
    { "strip", "Strip number mapping", "fStrip" },
    { "adc0", "Strip number mapping", "fadc0" },
    { "adc1", "Strip number mapping", "fadc1" },
    { "adc2", "Strip number mapping", "fadc2" },
    { "adc3", "Strip number mapping", "fadc3" },
    { "adc4", "Strip number mapping", "fadc4" },
    { "adc5", "Strip number mapping", "fadc5" },
    { 0 },
  };


  Int_t ret = DefineVarsFromList( vars, mode );

  if( ret != kOK )
    return ret;
    
  return kOK;
    
}

void    SBSGEMModule::Clear( Option_t* opt){
  fNch = 0;
  fIsDecoded = false;
  
  return;
}

Int_t   SBSGEMModule::Decode( const THaEvData& evdata ){
  //    std::cout << "[SBSGEMModule::Decode " << fName << "]" << std::endl;

  int i;

  //fNch = 0;
  fNstrips_hit = 0;
  //This could be written much more efficiently, in principle. However, it's not yet clear it's a speed bottleneck, so for now let's not worry about it:
  for (std::vector<mpdmap_t>::iterator it = fMPDmap.begin() ; it != fMPDmap.end(); ++it){
    //loop over all decode map entries associated with this module (each decode map entry is one APV card)
    Int_t effChan = it->mpd_id << 8 | it->adc_id; //left-shift mpd id by 8 bits and take the bitwise OR with ADC_id to uniquely identify the APV card.
    //mpd_id is not necessarily equal to slot, but that seems to be the convention in many cases
    // Find channel for this crate/slot

    Int_t nchan = evdata.GetNumChan( it->crate, it->slot );

    GEMaxis_t axis = it->axis == 0 ? SBSGEMModule::kUaxis : SBSGEMModule::kVaxis; 
    
    //        printf("nchan = %d\n", nchan );

    for( Int_t ichan = 0; ichan < nchan; ++ichan ) { //this is looping over all the "channels" (APV cards) in the crate and slot containing this decode map entry/APV card:
      Int_t chan = evdata.GetNextChan( it->crate, it->slot, ichan ); //"chan" here refers to one APV card 
      if( chan != effChan ) continue; // 


      Int_t nsamp = evdata.GetNumHits( it->crate, it->slot, chan );
      assert(nsamp%fN_MPD_TIME_SAMP==0); //this is making sure that the number of samples is equal to an integer multiple of the number of time samples per strip
      Int_t nstrips = nsamp/fN_MPD_TIME_SAMP; //number of strips fired on this APV card (should be exactly 128 if online zero suppression is NOT used):

      // Loop over all the strips and samples in the data
      //Int_t ihit = 0;
      for( Int_t istrip = 0; istrip < nstrips; ++istrip ) {
	//Temporary vector to hold ped-subtracted ADC samples for this strip:
	std::vector<double> ADCtemp;

	//sums over time samples
	double ADCsum_temp = 0.0;
	double maxADC = 0.0;
	UShort_t iSampMax = -1;
	
	//crude timing calculations:
	double Tsum = 0.0;
	double T2sum = 0.0;

	//At this point, "strip" equals APV channel number:
	Int_t strip = evdata.GetRawData(it->crate, it->slot, chan, fN_MPD_TIME_SAMP*istrip );

	assert(strip>=0&&strip<128);
	// Madness....   copy pasted from stand alone decoder
	// I bet there's a more elegant way to express this
	//Int_t RstripNb= 32*(strip%4)+ 8*(int)(strip/4)- 31*(int)(strip/16);
	//RstripNb=RstripNb+1+RstripNb%4-5*(((int)(RstripNb/4))%2);
	// New: Use a pre-computed array from Danning to skip the above
	// two steps.
	Int_t RstripNb = APVMAP[strip];
	RstripNb=RstripNb+(127-2*RstripNb)*it->invert;
	Int_t RstripPos = RstripNb + 128*it->pos;
	strip = RstripPos; //At this point, "strip" should correspond to increasing order of position along the U or V axis; i.e., what we think it should!

	//NOTE that we are replacing the value of "strip" with the line above!

	//Now loop over the time samples:
	for( Int_t adc_samp = 0; adc_samp < fN_MPD_TIME_SAMP; adc_samp++ ){

	  Int_t ihit = adc_samp + fN_MPD_TIME_SAMP * istrip; //index in the "hit" array for this APV card:
	  assert(ihit<nsamp);
	  
	  // Grab appropriate pedestal based on axis: existing code seems to assume that pedestal is specific to an individual strip, but does not vary
	  // sample-to-sample: When operating without online zero suppression, these should all probably be set to zero, since we will
	  // generally do offline common-mode calculation and subtraction in that case:
	  double pedtemp = ( axis == SBSGEMModule::kUaxis ) ? fPedestalU[strip] : fPedestalV[strip];
	  double rmstemp = ( axis == SBSGEMModule::kUaxis ) ? fPedRMSU[strip] : fPedRMSV[strip];
	  
	  double ADCvalue = double( evdata.GetData(it->crate, it->slot, chan, ihit) ) - pedtemp;
	  
	  //subtract "pedestal" from raw ADC value:
	  ADCtemp.push_back( ADCvalue );
	  // fadc[adc_samp][fNch] =  evdata.GetData(it->crate, it->slot,
	  // 					 chan, isamp++) - fPedestal[strip];

	  ADCsum_temp += ADCvalue;
	  
	  if( iSampMax < 0 || ADCvalue > maxADC ){
	    maxADC = ADCvalue;
	    iSampMax = adc_samp;
	  }

	  //for crude strip timing, just take simple time bins at the center of each sample (we'll worry about trigger time words later):
	  double Tsamp = fSamplePeriod * ( adc_samp + 0.5 );
	  
	  Tsum += Tsamp * ADCvalue;
	  T2sum += pow(Tsamp,2) * ADCvalue;
	  
	  //assert( ((UInt_t) fNch) < fMPDmap.size()*fN_APV25_CHAN );
	  //assert( fNstrips_hit < fMPDmap.size()*fN_APV25_CHAN );
	}
	assert(strip>=0); // Make sure we don't end up with negative strip numbers!

	// Zero suppression based on third time sample only?
	//Maybe better to do based on max ADC sample:
	if(!fZeroSuppress ||
	   ( fRMS[strip] > 0.0 && fabs(maxADC)/
	     fRMS[strip] > fZeroSuppressRMS ) ){ //Default threshold is 5-sigma!
	  //Increment hit count and populate decoded data structures:
	  fNstrips_hit++;
	  fStrip.push_back( strip );
	  fAxis.push_back( axis );
	  fADCsamples.push_back( ADCtemp ); //pedestal-subtracted
	  fADCsums.push_back( ADCsum_temp ); //sum of all (pedestal-subtracted) samples
	  fKeepStrip.push_back( true ); //keep all strips by default
	  fMaxSamp.push_back( iSampMax );
	  fADCmax.push_back( maxADC );
	  fTmean.push_back( Tsum/ADCsum_temp );
	  fTsigma.push_back( sqrt( T2sum/ADCsum_temp - pow( fTmean.back(), 2 ) ) );
	  fTcorr.push_back( fTmean.back() ); //don't apply any corrections for now
	}
      }
    }

  }

  //    std::cout << fName << " channels found  " << fNch << std::endl;

  fIsDecoded = true;
  
  return 0;
}

void SBSGEMModule::find_2Dhits(){ //version with no arguments calls 1D cluster finding with default (wide-open) track search constraints
  //these functions will fill the 1D cluster arrays:
  find_clusters_1D(SBSGEMModule::kUaxis); //u strips
  find_clusters_1D(SBSGEMModule::kVaxis); //v strips

  //Now make 2D clusters:

  fxcmin = -1.e12;
  fxcmax = 1.e12;
  fycmin = -1.e12;
  fycmax = 1.e12;
  
  fill_2D_hit_arrays();
  
}

void SBSGEMModule::find_2Dhits(TVector2 constraint_center, TVector2 constraint_width ){
  //constraint center and constraint width are assumed given in meters in "local" detector x,y coordinates.

  double ucenter = constraint_center.X() * fPxU + constraint_center.Y() * fPyU;
  double vcenter = constraint_center.X() * fPxV + constraint_center.Y() * fPyV;

  //To determine the constraint width along u/v, we need to transform the X/Y constraint widths, which define a rectangular region,
  //into U and V limits for 1D clustering:

  double umin,umax,vmin,vmax;

  double xmin = constraint_center.X() - constraint_width.X();
  double xmax = constraint_center.X() + constraint_width.X();
  double ymin = constraint_center.Y() - constraint_width.Y();
  double ymax = constraint_center.Y() + constraint_width.Y();

  //store these for later use:
  fxcmin = xmin;
  fxcmax = xmax;
  fycmin = ymin;
  fycmax = ymax;
  
  //check the four corners of the rectangle and compute the maximum values of u and v occuring at the four corners of the rectangular region:
  // NOTE: we will ALSO enforce the 2D search region in X and Y when we combine 1D U/V clusters into 2D X/Y hits, which, depending on the U/V strip orientation
  // can exclude some 2D hits that would have passed the U/V constraints defined by the corners of the X/Y rectangle, but been outside the X/Y constraint rectangle

  double u00 = xmin * fPxU + ymin * fPyU;
  double u01 = xmin * fPxU + ymax * fPyU;
  double u10 = xmax * fPxU + ymin * fPyU;
  double u11 = xmax * fPxU + ymax * fPyU;

  //this is some elegant-looking (compact) code, but perhaps algorithmically clunky:
  umin = std::min( u00, std::min(u01, std::min(u10, u11) ) );
  umax = std::max( u00, std::max(u01, std::max(u10, u11) ) );

  double v00 = xmin * fPxV + ymin * fPyV;
  double v01 = xmin * fPxV + ymax * fPyV;
  double v10 = xmax * fPxV + ymin * fPyV;
  double v11 = xmax * fPxV + ymax * fPyV;
  
  vmin = std::min( v00, std::min(v01, std::min(v10, v11) ) );
  vmax = std::max( v00, std::max(v01, std::max(v10, v11) ) );

  //The following routines will loop on all the strips and populate the 1D "cluster list" (vector<sbsgemcluster_t> )
  //Taking half the difference between max and min as the "width" is consistent with how it is used in
  // find_clusters_1D; i.e., the peak is required to be within |peak position - constraint center| <= constraint width
  find_clusters_1D(SBSGEMModule::kUaxis, ucenter, (umax-umin)/2.0 ); //U clusters
  find_clusters_1D(SBSGEMModule::kVaxis, vcenter, (vmax-vmin)/2.0 );  //V clusters

  //Now make 2D clusters
  
  fill_2D_hit_arrays();
}

void SBSGEMModule::find_clusters_1D( SBSGEMModule::GEMaxis_t axis, Double_t constraint_center, Double_t constraint_width ){

  //Constraint center and constraint width are assumed to be given in "standard" Hall A units (meters) in module-local coordinates
  // (SPECIFICALLY: constraint center and constraint width are assumed to refer to the direction measured by the strips being considered here)
  //This method will generally only be called by the reconstruction methods of SBSGEMTrackerBase
  
  if( !fIsDecoded ){
    cout << "find_clusters invoked before decoding for GEM Module " << GetName() << ", doing nothing" << endl;
    return;
  }

  UShort_t maxsep;
  UShort_t maxsepcoord; 
  UInt_t Nstrips;
  Double_t pitch;

  if( axis == SBSGEMModule::kUaxis ){
    maxsep = fMaxNeighborsU_totalcharge;
    maxsepcoord = fMaxNeighborsU_hitpos; 
    Nstrips = fNstripsU;
    pitch = fUStripPitch;
  } else { //V axis, no need to compare axis to kVaxis:
    maxsep = fMaxNeighborsV_totalcharge;
    maxsepcoord = fMaxNeighborsV_hitpos; 
    Nstrips = fNstripsV;
    pitch = fVStripPitch;
  }
  
  std::set<UShort_t> striplist;  //sorted list of strips for 1D clustering
  std::map<UShort_t, UInt_t> hitindex; //key = strip ID, mapped value = index in decoded hit array, needed to access the other information efficiently:

  for( int ihit=0; ihit<fNstrips_hit; ihit++ ){
    if( fAxis[ihit] == axis && fKeepStrip[ihit] ){
      
      bool newstrip = (striplist.insert( fStrip[ihit] ) ).second;

      if( newstrip ){ //should always be true:
	hitindex[fStrip[ihit]] = ihit;
      }
    }
  }

  std::set<UShort_t> localmaxima;
  std::map<UShort_t,bool> islocalmax;
  //std::map<UShort_t,bool> passed_constraint;
  
  
  for( std::set<UShort_t>::iterator i=striplist.begin(); i != striplist.end(); ++i ){
    int strip = *i;
    //int hitidx = hitindex[strip];
    islocalmax[strip] = false;

    double sumstrip = fADCsums[hitindex[strip]];
    double sumleft = 0.0;
    double sumright = 0.0;

    
    if( striplist.find( strip - 1 ) != striplist.end() ){
      sumleft = fADCsums[hitindex[strip-1]]; //if strip - 1 is found in strip list, hitindex is guaranteed to have been initialized above
    }
    if( striplist.find( strip + 1 ) != striplist.end() ){
      sumright = fADCsums[hitindex[strip+1]];
    }

    if( sumstrip >= sumleft && sumstrip >= sumright ){ //new local max:
      islocalmax[strip] = true;
      localmaxima.insert( strip );
    } 
  } // end loop over list of strips along this axis:

  vector<int> peakstoerase; 

  //now calculate "prominence" for all peaks and erase "insignificant" peaks:

  for( std::set<int>::iterator i=localmaxima.begin(); i != localmaxima.end(); ++i ){
    int stripmax = *i;

    double ADCmax = fADCsums[hitindex[stripmax]];
    double prominence = ADCmax;

    int striplo = stripmax, striphi = stripmax;
    double ADCminright=ADCmax, ADCminleft=ADCmax;

    bool higherpeakright=false,higherpeakleft=false;
    int peakright = -1, peakleft = -1;

    while( striplist.find( striphi+1 ) != striplist.end() ){
      striphi++;

      Double_t ADCtest = fADCsums[hitindex[striphi]];
      
      if( ADCtest < ADCminright && !higherpeakright ){ //as long as we haven't yet found a higher peak to the right, this is the lowest point between the current maximum and the next higher peak to the right:
	ADCminright = ADCtest;
      }

      if( islocalmax[striphi] && ADCtest > ADCmax ){ //then this peak is in a contiguous group with another higher peak to the right:
	higherpeakright = true;
	peakright = striphi;
      }
    }

    while( striplist.find(striplo-1) != striplist.end() ){
      striplo--;
      Double_t ADCtest = fADCsums[hitindex[striplo]];
      if( ADCtest < ADCminleft && !higherpeakleft ){ //as long as we haven't yet found a higher peak to the left, this is the lowest point between the current maximum and the next higher peak to the left:
	ADCminleft = ADCtest;
      }

      if( islocalmax[striplo] && ADCtest > ADCmax ){ //then this peak is in a contiguous group with another higher peak to the left:
	higherpeakleft = true;
	peakleft = striplo;
      }
    }

    double sigma_sum = sqrt( double(fN_MPD_TIME_SAMP) )*fZeroSuppressRMS; //~25-50 ADC

    bool peak_close = false;
    if( !higherpeakleft ) ADCminleft = 0.0;
    if( !higherpeakright ) ADCminright = 0.0;

    if( higherpeakright || higherpeakleft ){ //this peak is contiguous with higher peaks on either the left or right or both:
      prominence = ADCmax - std::max( ADCminleft, ADCminright );

      if( higherpeakleft && std::abs( peakleft - stripmax ) <= 2*maxsep ) peak_close = true;
      if( higherpeakright && std::abs( peakright - stripmax ) <= 2*maxsep ) peak_close = true;

      if( peak_close && (prominence < fThreshold_2ndMax_nsigma * sigma_sum ||
			 prominance/ADCmax < fThresh_2ndMax_fraction ) ){
	peakstoerase.push_back( stripmax );
      }
    }
  }

  //Erase "insignificant" peaks (those in contiguous grouping with higher peak with prominence below thresholds):
  for( int ipeak=0; ipeak<peakstoerase.size(); ipeak++ ){
    localmaxima.erase( peakstoerase[ipeak] );
    islocalmax[peakstoerase[ipeak]] = false;
  }

 
  //Cluster formation and cluster splitting from remaining local maxima:
  for( std::set<int>::iterator i = localmax.begin(); i != localmax.end(); ++i ){
    int stripmax = *i;
    int striplo = stripmax;
    int striphi = stripmax;
    double ADCmax = fADCsums[hitindex[stripmax]];

    while( striplist.find( striplo-1 ) != striplist.end() &&
	   stripmax - striplo < maxsep ){
      striplo--;
    }

    while( striplist.find( striphi+1 ) != striplist.end() &&
	   striphi - stripmax < maxsep ){
      striphi++;
    }

    int nstrips = striphi-striplo+1;

    double sumx = 0.0, sumx2 = 0.0, sumADC = 0.0, sumt = 0.0, sumt2 = 0.0;
    double sumwx = 0.0;

    map<int,double> splitfraction;
    vector<double> stripADCsum(nstrips);

    double maxpos = (stripmax + 0.5 - 0.5*Nstrips) * pitch;

    //If peak position falls inside the "track search region" constraint, add a new cluster: 
    if( fabs( maxpos - constraint_center ) <= constraint_width ){
    
      //create a cluster, but don't add it to the 1D cluster array unless it passes the track search region constraint:
      sbsgemcluster_t clusttemp;
      clusttemp.nstrips = nstrips;
      clusttemp.istriplo = striplo;
      clusttemp.istriphi = striphi;
      clusttemp.istripmax = stripmax;
      clusttemp.ADCsamples.resize(fN_MPD_TIME_SAMP);
      clusttemp.stripADCsum.clear();
      clusttemp.hitindex.clear();
      
      for( int isamp=0; isamp<fN_MPD_TIME_SAMP; isamp++ ){ //initialize cluster-summed ADC samples to zero:
	clusttemp.ADCsamples[isamp] = 0.0;
      }
      
      for( int istrip=striplo; istrip<=striphi; istrip++ ){
	int nmax_strip = 1;
	double sumweight = ADCmax/(1.0 + pow( (stripmax-istrip)*pitch/fSigma_hitshape, 2 ) );
	double maxweight = sumweight;
	for( int jstrip=istrip-maxsep; jstrip<=istrip+maxsep; jstrip++ ){
	  if( localmaxima.find( jstrip ) != localmaxima.end() && jstrip != stripmax ){
	    sumweight += fADCsums[hitindex[jstrip]]/( 1.0 + pow( (jstrip-istrip)*pitch/fSigma_hitshape, 2 ) );
	  }
	}
   
	splitfraction[istrip] = maxweight/sumweight;

	double hitpos = (istrip + 0.5 - 0.5*Nstrips) * pitch; //local hit position along direction measured by these strips
	double ADCstrip = fADCsums[hitindex[istrip]] * splitfraction[istrip];
	double tstrip = fTmean[hitindex[istrip]];

	for( int isamp=0; isamp<fN_MPD_TIME_SAMP; isamp++ ){
	  clusttemp.ADCsamples[isamp] += fADCsamples[hitindex[istrip]]*splitfraction[istrip];
	}

	clusttemp.stripADCsum.push_back( ADCstrip );

	clusttemp.hitindex.push_back( hitindex[istrip] );
	
	sumADC += ADCstrip;
	
	if( std::abs( istrip - stripmax ) <= std::max(1,std::min(fMaxNeighborsU_hitpos,fMaxNeighborsU_totalcharge)) ){ 
	  sumx += hitpos * ADCstrip;
	  sumx2 += pow(hitpos,2) * ADCstrip;
	  sumwx += ADCstrip;
	  //use same strip cuts for cluster timing determination as for position reconstruction: may revisit later:
	  sumt += tstrip * ADCstrip;
	  sumt2 += pow(tstrip,2) * ADCstrip;
	} 
      }

      clusttemp.hitpos_mean = sumx / sumwx;
      clusttemp.hitpos_sigma = sqrt( sumx2/sumwx - pow(clusttemp.hitpos_mean,2) );
      clusttemp.clusterADCsum = sumADC;
      clusttemp.t_mean = sumt / sumwx;
      clusttemp.t_sigma = sqrt( sumt2 / sumwx - pow(clusttemp.t_mean,2) );

      if( axis == SBSGEMModule::kVaxis ){
	fVclusters.push_back( clusttemp );
      } else {
	fUclusters.push_back( clusttemp );
      }
    } //Check if peak is inside track search region constraint
  } //end loop on local maxima
}

void SBSGEMModule::fill_2D_hit_arrays(){
  //Clear out the 2D hit array to get rid of any leftover junk from prior events:
  fHits.clear();
  
  //This routine is simple: just form all possible 2D hits from combining one "U" cluster with one "V" cluster. Here we assume that find_clusters_1D has already
  //been called, if that is NOT the case, then this routine will just do nothing:
  for( int iu=0; iu<fUclusters.size(); iu++ ){
    for( int iv=0; iv<fVclusters.size(); iv++ ){
      //Initialize sums for computing cluster and strip correlation coefficients:
      sbsgemhit_t hittemp; // declare a temporary "hit" object:

      //copying overhead, might be inefficient:
      //sbsgemcluster_t uclusttemp = fUclusters[iu];
      //sbsgemcluster_t vclusttemp = fVclusters[iv];
      
      //Initialize "keep" to true:
      hittemp.keep = true;
      hittemp.ontrack = false;
      hittemp.trackidx = -1;
      hittemp.iuclust = iu;
      hittemp.ivclust = iv;
      
      hittemp.uhit = fUclusters[iu].hitpos_mean;
      hittemp.vhit = fVclusters[iv].hitpos_mean;

      double pos_maxstripu = ( fUclusters[iu].istripmax + 0.5 - 0.5*fNstripsU ) * fUStripPitch;
      double pos_maxstripv = ( fVclusters[iv].istripmax + 0.5 - 0.5*fNstripsV ) * fVStripPitch;

      //"Cluster moments" defined as difference between reconstructed hit position and center of strip with max. signal in the cluster:
      hittemp.umom = (hittemp.uhit - pos_maxstripu)/fUStripPitch;
      hittemp.vmom = (hittemp.vhit - pos_maxstripv)/fVStripPitch;
      
      TVector2 UVtemp(hittemp.uhit,hittemp.vhit);
      TVector2 XYtemp = UVtoXY( UVtemp );

      hittemp.xhit = XYtemp.X();
      hittemp.yhit = XYtemp.Y();

      //Check if candidate 2D hit is inside the constraint region before doing anything else:
      if( fxcmin <= hittemp.xhit && hittemp.xhit <= fxcmax &&
	  fycmin <= hittemp.yhit && hittemp.yhit <= fycmax ){
    
	hittemp.thit = 0.5*(fUclusters[iu].t_mean + fVclusters[iv].t_mean);
	hittemp.Ehit = 0.5*(fUclusters[iu].clusterADCsum + fVclusters[iv].clusterADCsum);
	
	hittemp.thitcorr = hittemp.thit; //don't apply any corrections on thit yet
	
	//Next up is to calculate "global" hit coordinates (actually coordinates in "tracker-local" system)
	//DetToTrackCoord is a utility function defined in THaDetectorBase
	TVector3 hitpos_global = DetToTrackCoord( hittemp.xhit, hittemp.yhit );
	
	// Unclear whether it is actually necessary to store these variables, but we also probably want to avoid 
	// repeated calls to THaDetectorBase::DetToTrackCoord, so let's keep these for now:
	hittemp.xghit = hitpos_global.X();
	hittemp.yghit = hitpos_global.Y();
	hittemp.zghit = hitpos_global.Z();
	
	hittemp.ADCasym = ( fUclusters[iu].clusterADCsum - fVclusters[iv].clusterADCsum )/( 2.0*hittemp.Ehit );
	hittemp.tdiff = fUclusters[iu].t_mean - fVclusters[iv].t_mean;
	
	//Calculate correlation coefficients:
	hittemp.corrcoeff_clust = CorrCoeff( fN_MPD_TIME_SAMP, fUclusters[iu].ADCsamples, fVclusters[iv].ADCsamples );
	
	//compute index of strip with max ADC sum within cluster strip array:
	int ustripidx = fUclusters[iu].istripmax-fUclusters[iu].istriplo; //
	int vstripidx = fVclusters[iv].istripmax-fVclusters[iv].istriplo; // 
	
	//compute index of strip with max ADC sum within decoded hit array:
	int uhitidx = fUclusters[iu].hitindex[ustripidx]; 
	int vhitidx = fVclusters[iv].hitindex[vstripidx];
	
	hittemp.corrcoeff_strip = CorrCoeff( fN_MPD_TIME_SAMP, fADCsamples[uhitidx], fADCsamples[vhitidx] );
	
	//Okay, that should be everything. Now add it to the 2D hit array:
	fHits.push_back( hittemp );
	
      } //end check that 2D point is inside track search region
    } //end loop over "V" clusters
  } //end loop over "U" clusters
  
}

void    SBSGEMModule::Print( Option_t* opt) const{
  return;
}

Int_t   SBSGEMModule::Begin( THaRunBase* r){
  return 0;
}

Int_t   SBSGEMModule::End( THaRunBase* r){
  return 0;
}

//utility method to calculate correlation coefficient of U and V samples: 
Double_t CorrCoeff( int nsamples, std::vector<double> Usamples, std::vector<double> Vsamples ){
  double sumu=0.0, sumv=0.0, sumu2=0.0, sumv2=0.0, sumuv=0.0;

  if ( Usamples.size() < nsamples || Vsamples.size() < nsamples ){
    return -10.0; //nonsense value, correlation coefficient by definition is -1 < c < 1
  }
  
  for( int isamp=0; isamp<nsamples; isamp++ ){
    sumu += Usamples[isamp];
    sumv += Vsamples[isamp];
    sumu2 += pow(Usamples[isamp],2);
    sumv2 += pow(Vsamples[isamp],2);
    sumuv += Usamples[isamp]*Vsamples[isamp];
  }

  double nSAMP = double(nsamples);
  double mu = sumu/nSAMP;
  double mv = sumv/nSAMP;
  double varu = sumu2/nSAMP - pow(mu,2);
  double varv = sumv2/nSAMP - pow(mv,2);
  double sigu = sqrt(varu);
  double sigv = sqrt(varv);

  return (sumuv - nSAMP*mu*mv)/(nSAMP*sigu*sigv);
  
}

TVector2 SBSGEMModule::UVtoXY( TVector2 UV ){
  double det = fPxU*fPyV - fPyU*fPxV;

  double Utemp = UV.X();
  double Vtemp = UV.Y();
  
  double Xtemp = (fPyV*Utemp - fPyU*Vtemp)/det;
  double Ytemp = (fPxU*Vtemp - fPxV*Utemp)/det;

  return TVector2(Xtemp,Ytemp);
}

TVector2 SBSGEMModule::XYtoUV( TVector2 XY ){
  double Xtemp = XY.X();
  double Ytemp = XY.Y();

  double Utemp = Xtemp*fPxU + Ytemp*fPyU;
  double Vtemp = Xtemp*fPxV + Ytemp*fPyV;

  return TVector2(Utemp,Vtemp);
}
