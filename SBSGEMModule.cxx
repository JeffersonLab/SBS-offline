#include <iostream>

#include "SBSGEMModule.h"
#include "TDatime.h"
#include "THaEvData.h"
#include "THaApparatus.h"
#include "TRotation.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TClonesArray.h"
#include <algorithm>

using namespace std;
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

  fPedestalMode = kFALSE;
  
  //Default online zero suppression to FALSE: actually I wonder if it would be better to use this in 
  // Moved this to MPDModule, since this should be done during the decoding of the raw APV data:
  fOnlineZeroSuppression = kFALSE;
  
  //Set default values for decode map parameters:
  fN_APV25_CHAN = 128;
  fN_MPD_TIME_SAMP = 6;
  fMPDMAP_ROW_SIZE = 9;

  fNumberOfChannelInFrame = 129;

  fSamplePeriod = 25.0; //nanoseconds:

  fSigma_hitshape = 0.0004; //0.4 mm
  // for( Int_t i = 0; i < N_MPD_TIME_SAMP; i++ ){
  //     fadc[i] = NULL;
  // }

  fMakeEfficiencyPlots = true;
  
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

  std::vector<Double_t> rawpedu,rawpedv;
  std::vector<Double_t> rawrmsu,rawrmsv;

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
    { "pedu",           &rawpedu,        kDoubleV, 0, 0, 0}, // optional raw pedestal info (u strips)
    { "pedv",           &rawpedv,        kDoubleV, 0, 0, 0}, // optional raw pedestal info (v strips)
    { "rmsu",           &rawrmsu,        kDoubleV, 0, 0, 0}, // optional pedestal rms info (u strips)
    { "rmsv",           &rawrmsv,        kDoubleV, 0, 0, 0}, // optional pedestal rms info (v strips)
    { "layer",          &fLayer,         kUShort, 0, 0, 0}, // mandatory: logical tracking layer must be specified for every module:
    { "nstripsu",       &fNstripsU,     kUInt, 0, 0, 1}, //mandatory: number of strips in module along U axis
    { "nstripsv",       &fNstripsV,     kUInt, 0, 0, 1}, //mandatory: number of strips in module along V axis
    { "uangle",         &fUAngle,       kDouble, 0, 0, 1}, //mandatory: Angle of "U" strips wrt X axis
    { "vangle",         &fVAngle,       kDouble, 0, 0, 1}, //mandatory: Angle of "V" strips wrt X axis
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
    { "maxnv_charge", &fMaxNeighborsV_totalcharge, kUShort, 0, 1, 1}, //(optional): cluster size restriction along V for total charge calculation
    { "maxnu_pos", &fMaxNeighborsU_hitpos, kUShort, 0, 1, 1}, //(optional): cluster size restriction for position reconstruction
    { "maxnv_pos", &fMaxNeighborsV_hitpos, kUShort, 0, 1, 1}, //(optional): cluster size restriction for position reconstruction
    { "sigmahitshape", &fSigma_hitshape, kDouble, 0, 1, 1}, //(optional): width parameter for cluster-splitting algorithm
    { "zerosuppress_nsigma", &fZeroSuppressRMS, kDouble, 0, 1, 1}, //(optional): 
    {0}
  };
  status = LoadDB( file, date, request, fPrefix, 1 ); //The "1" after fPrefix means search up the tree

  
  
  
  //std::cout << GetName() << " fThresholdStripSum " << fThresholdStripSum 
  //<< " fThresholdSample " << fThresholdSample << std::endl;
  
  if( status != 0 ){
    fclose(file);
    return status;
  }

  fPxU = cos( fUAngle * TMath::DegToRad() );
  fPyU = sin( fUAngle * TMath::DegToRad() );
  fPxV = cos( fVAngle * TMath::DegToRad() );
  fPyV = sin( fVAngle * TMath::DegToRad() );

  fAPVch_by_Ustrip.clear();
  fAPVch_by_Vstrip.clear();
  fMPDID_by_Ustrip.clear();
  fMPDID_by_Vstrip.clear();
  fADCch_by_Ustrip.clear();
  fADCch_by_Vstrip.clear();
  
  Int_t nentry = fChanMapData.size()/fMPDMAP_ROW_SIZE;
  for( Int_t mapline = 0; mapline < nentry; mapline++ ){
    mpdmap_t thisdata;
    thisdata.crate  = fChanMapData[0+mapline*fMPDMAP_ROW_SIZE];
    thisdata.slot   = fChanMapData[1+mapline*fMPDMAP_ROW_SIZE];
    thisdata.mpd_id = fChanMapData[2+mapline*fMPDMAP_ROW_SIZE];
    thisdata.gem_id = fChanMapData[3+mapline*fMPDMAP_ROW_SIZE];
    thisdata.adc_id = fChanMapData[4+mapline*fMPDMAP_ROW_SIZE];
    thisdata.i2c    = fChanMapData[5+mapline*fMPDMAP_ROW_SIZE];
    thisdata.pos    = fChanMapData[6+mapline*fMPDMAP_ROW_SIZE];
    thisdata.invert = fChanMapData[7+mapline*fMPDMAP_ROW_SIZE];
    thisdata.axis   = fChanMapData[8+mapline*fMPDMAP_ROW_SIZE];

    //Populate relevant quantities mapped by strip index:
    for( int ich=0; ich<fN_APV25_CHAN; ich++ ){
      int strip = GetStripNumber( ich, thisdata.pos, thisdata.invert );
      if( thisdata.axis == SBSGEM::kUaxis ){
	fAPVch_by_Ustrip[strip] = ich;
	fMPDID_by_Ustrip[strip] = thisdata.mpd_id;
	fADCch_by_Ustrip[strip] = thisdata.adc_id;
      } else {
	fAPVch_by_Vstrip[strip] = ich;
	fMPDID_by_Vstrip[strip] = thisdata.mpd_id;
	fADCch_by_Vstrip[strip] = thisdata.adc_id;
      }
    }
    
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
  
  // std::cout << "got " << rawpedu.size() << " u pedestal mean values and " << rawrmsu.size() << " u pedestal rms values" << std::endl;
  // std::cout << "got " << rawpedv.size() << " v pedestal mean values and " << rawrmsv.size() << " v pedestal rms values" << std::endl;

  // for( int i=0; i<rawpedu.size(); i++ ){
  //   cout << i << ", " << rawpedu[i] << endl;
  // }
  
  for ( UInt_t istrip=0; istrip<fNstripsU; istrip++ ){
    fPedestalU[istrip] = 0.0;
    fPedRMSU[istrip] = 10.0; //placeholder to be replaced by value from database
    
    if( rawpedu.size() == fNstripsU ){
      fPedestalU[istrip] = rawpedu[istrip];
    }else if( rawpedu.size() == fNstripsU/128 ){
      fPedestalU[istrip] = rawpedu[istrip/128];
    }else if(rawpedu.size()){ 
      fPedestalU[istrip] = rawpedu[0];
    }

    if( rawrmsu.size() == fNstripsU ){
      fPedRMSU[istrip] = rawrmsu[istrip];
    }else if( rawrmsu.size() == fNstripsU/128 ){
      fPedRMSU[istrip] = rawrmsu[istrip/128];
    }else if(rawrmsu.size()){ 
      fPedRMSU[istrip] = rawrmsu[0];
    }
 
  }

  //Initialize all pedestals to zero, RMS values to default:
  fPedestalV.clear();
  fPedestalV.resize( fNstripsV ); 

  fPedRMSV.clear();
  fPedRMSV.resize( fNstripsV );
  
  for( UInt_t istrip=0; istrip<fNstripsV; istrip++ ){
    fPedestalV[istrip] = 0.0;
    fPedRMSV[istrip] = 10.0;

    if( rawpedv.size() == fNstripsV ){
      fPedestalV[istrip] = rawpedv[istrip];
    }else if( rawpedv.size() == fNstripsV/128 ){
      fPedestalV[istrip] = rawpedv[istrip/128];
    }else if(rawpedv.size()){ 
      fPedestalV[istrip] = rawpedv[0];
    }

    if( rawrmsv.size() == fNstripsV ){
      fPedRMSV[istrip] = rawrmsv[istrip];
    }else if( rawrmsv.size() == fNstripsV/128 ){
      fPedRMSV[istrip] = rawrmsv[istrip/128];
    }else if(rawrmsv.size()){ 
      fPedRMSV[istrip] = rawrmsv[0];
    } 
  }

  //default all gains to 1 if they were not loaded from the DB:
  if( fUgain.size() != fNstripsU/128 ){
    fUgain.resize(fNstripsU/128);
    for( int iAPV=0; iAPV<fNstripsU/128; iAPV++ ){
      fUgain[iAPV] = 1.0;
    }
  }

  if( fVgain.size() != fNstripsV/128 ){
    fVgain.resize(fNstripsV/128);
    for( int iAPV=0; iAPV<fNstripsV/128; iAPV++ ){
      fVgain[iAPV] = 1.0;
    }
  }

  if( fPedestalMode ){
    fZeroSuppress = false;
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

Int_t SBSGEMModule::ReadGeometry( FILE *file, const TDatime &date, Bool_t required ){ //We start with a copy of THaDetectorBase::ReadGeometry and modify accordingly:
  // Read this detector's basic geometry information from the database.
  // Derived classes may override to read more advanced data.

  const char* const here = "ReadGeometry";

  vector<double> position, size, angles;
  Bool_t optional = !required;
  DBRequest request[] = {
    { "position", &position, kDoubleV, 0, optional, 0,
      "\"position\" (detector position [m])" },
    { "size",     &size,     kDoubleV, 0, optional, 1,
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

      
      
      RotTemp.RotateX( -angles[0] * TMath::DegToRad() );
      RotTemp.RotateY( -angles[1] * TMath::DegToRad() );
      RotTemp.RotateZ( -angles[2] * TMath::DegToRad() );
      
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
  
  //Raw strip info:
  RVarDef varstrip[] = {
    { "nstripsfired",   "Number of strips fired",   "fNstrips_hit" },
    { "strip", "Strip index", "fStrip" },
    { "stripIsU", "U strip?", "fStripIsU"},
    { "stripIsV", "V strip?", "fStripIsV"},
    { "stripADCsamples", "ADC samples (index = isamp+Nsamples*istrip)", "fADCsamples1D" },
    { "striprawADCsamples", "raw ADC samples (no baseline subtraction)", "fRawADCsamples1D" },
    { "stripADCsum", "Sum of ADC samples on a strip", "fADCsums" },
    { "stripisampmax", "sample in which max ADC occurred on a strip", "fMaxSamp" },
    { "stripADCmax", "Value of max ADC sample on a strip", "fADCmax" },
    { "stripTmean", "ADC-weighted mean strip time", "fTmean" },
    { "stripTsigma", "ADC-weighted rms strip time", "fTsigma" },
    { "stripTcorr", "Corrected strip time", "fTcorr" },
    { "stripItrack", "Index of track containing this strip (-1 if not on any track)", "fStripTrackIndex" },
    { "stripADCavg", "average of ADC samples on a strip", "fStripADCavg" },
    { nullptr },
  };


  Int_t ret = DefineVarsFromList( varstrip, mode );

  if( ret != kOK )
    return ret;
  

  RVarDef varclust[] = {
    { "nclustu",   "Number of clusters in u",   "fNclustU" },
    { "clustu_strips",   "u clusters strip multiplicity",   "fUclusters.nstrips" },
    { "clustu_pos",   "u clusters position",   "fUclusters.hitpos_mean" },
    { "clustu_adc",   "u clusters adc sum",   "fUclusters.clusterADCsum" },
    { "clustu_time",   "u clusters time",   "fUclusters.t_mean" },
    { "nclustv",   "Number of clusters in v",   "fNclustV" },
    { "clustv_strips",   "v clusters strip multiplicity",   "fVclusters.nstrips" },
    { "clustv_pos",   "v clusters position",   "fVclusters.hitpos_mean" },
    { "clustv_adc",   "v clusters adc sum",   "fVclusters.clusterADCsum" },
    { "clustv_time",   "v clusters time",   "fVclusters.t_mean" },
    { nullptr },
  };

  ret = DefineVarsFromList( varclust, mode );

  if( ret != kOK )
    return ret;

  RVarDef varhits[] = {
    { "nhits2d",   "Number of 2d hits",   "fN2Dhits" },
    { "hitx",   "local X coordinate of hit",   "fHits.xhit" },
    { "hity",   "local Y coordinate of hit",   "fHits.yhit" },
    { "hitxg",   "transport X coordinate of hit",   "fHits.xghit" },
    { "hityg",   "transport Y coordinate of hit",   "fHits.yghit" },
    { "hitADCasym",   "hit ADC asymmetry",   "fHits.ADCasym" },
    { "hitTdiff",   "hit time difference (u-v)",   "fHits.tdiff" },
    { nullptr },
  };

  ret = DefineVarsFromList( varhits, mode );

  if( ret != kOK )
    return ret;

  return kOK;
    
}

void SBSGEMModule::Clear( Option_t* opt){ //we will want to clear out many more things too
  fNstrips_hit = 0;
  fIsDecoded = false;

  fNclustU = 0;
  fNclustV = 0;
  fUclusters.clear();
  fVclusters.clear();
  fN2Dhits = 0;
  fHits.clear();

  //fStripAxis.clear();
  fADCsamples1D.clear();
  fStripTrackIndex.clear();
  fRawADCsamples1D.clear();

  fStripIsU.clear();
  fStripIsV.clear();
  fStripADCavg.clear();
  
  fUstripIndex.clear();
  fVstripIndex.clear();
  fStrip.clear();
  fAxis.clear();
  fADCsamples.clear();
  fRawADCsamples.clear();
  fADCsums.clear();
  fKeepStrip.clear();
  fMaxSamp.clear();
  fADCmax.clear();
  fTmean.clear();
  fTsigma.clear();
  fTcorr.clear();
  
  return;
}

Int_t   SBSGEMModule::Decode( const THaEvData& evdata ){
  //std::cout << "[SBSGEMModule::Decode " << fName << "]" << std::endl;

  fNstrips_hit = 0;
  //initialize "U" and "V" strip counters to zero:
  fNstrips_hitU = 0;
  fNstrips_hitV = 0;
  
  fUstripIndex.clear();
  fVstripIndex.clear();
  //This could be written more efficiently, in principle. However, it's not yet clear it's a speed bottleneck, so for now let's not worry about it too much:
  for (std::vector<mpdmap_t>::iterator it = fMPDmap.begin() ; it != fMPDmap.end(); ++it){
    //loop over all decode map entries associated with this module (each decode map entry is one APV card)
    Int_t effChan = it->mpd_id << 4 | it->adc_id; //left-shift mpd id by 4 bits and take the bitwise OR with ADC_id to uniquely identify the APV card.
    //mpd_id is not necessarily equal to slot, but that seems to be the convention in many cases
    // Find channel for this crate/slot
    Int_t nchan = evdata.GetNumChan( it->crate, it->slot );

    SBSGEM::GEMaxis_t axis = it->axis == 0 ? SBSGEM::kUaxis : SBSGEM::kVaxis; 
    
    //printf("nchan = %d\n", nchan );

    for( Int_t ichan = 0; ichan < nchan; ++ichan ) { //this is looping over all the "channels" (APV cards) in the crate and slot containing this decode map entry/APV card:
      Int_t chan = evdata.GetNextChan( it->crate, it->slot, ichan ); //"chan" here refers to one APV card 
      //std::cout << it->crate << " " << it->slot << " mpd_id ??? " << it->mpd_id << " " << chan << " " << effChan << std::endl;
      
      if( chan != effChan ) continue; // 


      Int_t nsamp = evdata.GetNumHits( it->crate, it->slot, chan );
      assert(nsamp%fN_MPD_TIME_SAMP==0); //this is making sure that the number of samples is equal to an integer multiple of the number of time samples per strip
      Int_t nstrips = nsamp/fN_MPD_TIME_SAMP; //number of strips fired on this APV card (should be exactly 128 if online zero suppression is NOT used):

      //std::cout << "MPD ID, ADC channel, number of strips fired = " << it->mpd_id << ", "
      //	<< it->adc_id << ", " << nstrips << std::endl;
      
      //declare temporary array to hold common mode values for this APV card and, if necessary, calculate them:
      double commonMode[fN_MPD_TIME_SAMP];

      for( int isamp=0; isamp<fN_MPD_TIME_SAMP; isamp++ ){
	commonMode[isamp] = 0.0;
      }

      vector<int> Strip( nsamp );
      vector<int> rawStrip( nsamp );
      vector<int> rawADC( nsamp );
      vector<double> pedsubADC( nsamp ); //ped-subtracted, not necessarily common-mode subtracted
      vector<double> commonModeSubtractedADC( nsamp );

      //First loop over the hits: populate strip, raw strip, raw ADC, ped sub ADC and common-mode-subtracted aDC:
      for( int iraw=0; iraw<nsamp; iraw++ ){ //NOTE: iraw = isamp + fN_MPD_TIME_SAMP * istrip
	int strip = evdata.GetRawData( it->crate, it->slot, chan, iraw );
	int ADC = evdata.GetData( it->crate, it->slot, chan, iraw );

	rawStrip[iraw] = strip;
	Strip[iraw] = GetStripNumber( strip, it->pos, it->invert );

	double ped = (axis == SBSGEM::kUaxis ) ? fPedestalU[Strip[iraw]] : fPedestalV[Strip[iraw]];

	//not sure if this is correct, ask Ben:
	if( fOnlineZeroSuppression ) ped = 0.0;
	
	rawADC[iraw] = ADC;
	pedsubADC[iraw] = double(ADC) - ped;
	commonModeSubtractedADC[iraw] = double(ADC) - ped; 

	//the calculation of common mode in pedestal mode analysis differs from the
	// offline or online zero suppression analysis; here we use a simple average of all 128 channels:
	if( fPedestalMode && nstrips == fN_APV25_CHAN ){
	  //do simple common-mode calculation involving the simple average of all 128 (ped-subtracted) ADC
	  //values
	  int isamp = iraw%fN_MPD_TIME_SAMP;
	  commonMode[isamp] += (double(ADC)-ped)/double(nstrips);
	}
      }
      
      //(OPTIONAL) second loop over the hits to calculate and apply common-mode correction (sorting method)
      if( !fOnlineZeroSuppression && nstrips == fN_APV25_CHAN && !fPedestalMode ){ //need to calculate common mode:
	for( int isamp=0; isamp<fN_MPD_TIME_SAMP; isamp++ ){
	  vector<double> sortedADCs(fN_APV25_CHAN);
	  
	  for( int ihit=0; ihit<fN_APV25_CHAN; ihit++ ){
	    int iraw = isamp + fN_MPD_TIME_SAMP * ihit;	    
	    sortedADCs[ihit] = pedsubADC[ iraw ];
	  }
	  
	  std::sort( sortedADCs.begin(), sortedADCs.end() );

	  commonMode[isamp] = 0.0;
	  
	  for( int k=28; k<100; k++ ){
	    commonMode[isamp] += sortedADCs[k];
	  }
	  commonMode[isamp] /= 72.0;
	  //moved this to the next loop for efficiency and code re-use:
	  // for( int ihit=0; ihit<fN_APV25_CHAN; ihit++ ){
	  //   int iraw = isamp + fN_MPD_TIME_SAMP * ihit;	
	  //   commonModeSubtractedADC[ iraw ] = pedsubADC[ iraw ] - commonMode[isamp];
	  // }
	  
	}
      }
      
      //std::cout << "finished common mode " << std::endl;
      // Last loop over all the strips and samples in the data and populate/calculate global variables that are passed to track-finding:
      //Int_t ihit = 0;
      for( Int_t istrip = 0; istrip < nstrips; ++istrip ) {
	//Temporary vector to hold ped-subtracted ADC samples for this strip:
	std::vector<double> ADCtemp;
	std::vector<int> rawADCtemp;
	
	//sums over time samples
	double ADCsum_temp = 0.0;
	double maxADC = 0.0;
	UShort_t iSampMax = -1;
	
	//crude timing calculations:
	double Tsum = 0.0;
	double T2sum = 0.0;

	//grab decoded strip number directly:
	int strip = Strip[fN_MPD_TIME_SAMP * istrip];
	
	//Pedestal has already been subtracted by the time we get herre, but let's grab anyway in case it's needed:
	
	double pedtemp = ( axis == SBSGEM::kUaxis ) ? fPedestalU[strip] : fPedestalV[strip];
	double rmstemp = ( axis == SBSGEM::kUaxis ) ? fPedRMSU[strip] : fPedRMSV[strip];
	double gaintemp = ( axis == SBSGEM::kUaxis ) ? fUgain[strip/128] : fVgain[strip/128]; //should probably not hard-code 128 here
	
	// std::cout << "pedestal temp, rms temp, nsigma cut, threshold, zero suppress, pedestal mode = " << pedtemp << ", " << rmstemp << ", " << fZeroSuppressRMS
	//   	  << ", " << fZeroSuppressRMS * rmstemp << ", " << fZeroSuppress << ", " << fPedestalMode << std::endl;
	
	//Now loop over the time samples:
	for( Int_t adc_samp = 0; adc_samp < fN_MPD_TIME_SAMP; adc_samp++ ){

	  int iraw = adc_samp + fN_MPD_TIME_SAMP * istrip;

	  //If applicable, subtract common-mode here:
	  if( (fPedestalMode || !fOnlineZeroSuppression) && nstrips == fN_APV25_CHAN ){
	    commonModeSubtractedADC[ iraw ] = pedsubADC[ iraw ] - commonMode[adc_samp];
	  }
	  
	  // Int_t ihit = adc_samp + fN_MPD_TIME_SAMP * istrip; //index in the "hit" array for this APV card:
	  // assert(ihit<nsamp);
	  
	  
	  //Int_t rawADC = evdata.GetData(it->crate, it->slot, chan, ihit);
	  Int_t RawADC = rawADC[iraw]; //this value has no corrections applied:
         
	  //cout << adc_samp << " " << istrip << " " << rawADC << " ";// << endl;
	  
	  rawADCtemp.push_back( RawADC );

	  //The following value already has pedestal and common-mode subtracted:
	  double ADCvalue = commonModeSubtractedADC[iraw]; //zero-suppress BEFORE we apply gain correction

	  // if( fPedestalMode ){ //If we are analyzing pedestal data, DON'T substract the pedestal
	  //   ADCvalue = commonModeSubtractedADC[adc_samp][istrip];
	  // }
	  //pedestal-subtracted ADC values:
	  ADCtemp.push_back( ADCvalue );
	  // fadc[adc_samp][fNch] =  evdata.GetData(it->crate, it->slot,
	  // 					 chan, isamp++) - fPedestal[strip];

	  ADCsum_temp += ADCvalue;
	  //cout << ADCvalue << " "<< endl;
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

	//cout << ADCsum_temp << " >? " << fThresholdStripSum << "; " 
	//   << ADCsum_temp/double(fN_MPD_TIME_SAMP) << " >? " << fZeroSuppressRMS*rmstemp << "; " 
	//   << maxADC << " >? " << fZeroSuppressRMS*rmstemp 
	//   << "; " << fThresholdSample << endl;
	
	//the ROOTgui multicrate uses a threshold on the AVERAGE ADC sample (not the MAX). To be consistent
	// with how the "Hit" root files are produced, let's use the same threshold;
	// this amounts to using a higher effective threshold than cutting on the max ADC sample would have been:
	//fThresholdStripSum is in many respects redundant with fZeroSuppressRMS
	if(!fZeroSuppress ||
	   ( ADCsum_temp/double(fN_MPD_TIME_SAMP) > fZeroSuppressRMS*rmstemp &&
	     maxADC > fThresholdSample && ADCsum_temp > fThresholdStripSum ) ){ //Default threshold is 5-sigma!
	  //Increment hit count and populate decoded data structures:
	  //Theshold on the max. sample is also not used in the standalone decoder, but for sufficiently low values should be redundant with the
	  //threshold on the average ADC
	  
	  //NOW apply gain correction:
	  //Don't apply any gain correction if we are doing pedestal mode analysis:
	  for( Int_t isamp=0; isamp<fN_MPD_TIME_SAMP; isamp++ ){
	    if( !fPedestalMode ) ADCtemp[isamp] /= gaintemp;
	  }
	  
	  //  ADCsum_temp /= gaintemp;
	
	  
	  fStrip.push_back( strip );
	  fAxis.push_back( axis );
	  //std::cout << "axis, Int_t(axis) = " << axis << ", " << Int_t(axis) << std::endl;
	  //fStripAxis.push_back( Int_t(axis) );
	  fADCsamples.push_back( ADCtemp ); //pedestal-subtracted
	  fRawADCsamples.push_back( rawADCtemp ); //Raw
	  for( Int_t isamp=0; isamp<fN_MPD_TIME_SAMP; isamp++ ){
	    fADCsamples1D.push_back( ADCtemp[isamp] );
	    fRawADCsamples1D.push_back( rawADCtemp[isamp] );
	  }
	  fStripTrackIndex.push_back( -1 ); //This could be modified later based on tracking results
	  
	  fKeepStrip.push_back( true ); //keep all strips by default
	  fMaxSamp.push_back( iSampMax );

	  if( !fPedestalMode ) maxADC /= gaintemp;
	  
	  fADCmax.push_back( maxADC );
	  fTmean.push_back( Tsum/ADCsum_temp );
	  fTsigma.push_back( sqrt( T2sum/ADCsum_temp - pow( fTmean.back(), 2 ) ) );
	  fTcorr.push_back( fTmean.back() ); //don't apply any corrections for now

	  if( !fPedestalMode ) ADCsum_temp /= gaintemp;

	  fADCsums.push_back( ADCsum_temp ); //sum of all (pedestal-subtracted) samples

	  fStripADCavg.push_back( ADCsum_temp/double(fN_MPD_TIME_SAMP) );

	  UInt_t isU = (axis == SBSGEM::kUaxis) ? 1 : 0;
	  UInt_t isV = (axis == SBSGEM::kVaxis) ? 1 : 0;
	  fStripIsU.push_back( isU );
	  fStripIsV.push_back( isV );
	  
	  if( axis == SBSGEM::kUaxis ) fUstripIndex[strip] = fNstrips_hit;
	  if( axis == SBSGEM::kVaxis ) fVstripIndex[strip] = fNstrips_hit;

	  //std::cout << "starting pedestal histograms..." << std::endl;
	  
	  if( fPedestalMode ){ //fill pedestal histograms:
	    int iAPV = strip/fN_APV25_CHAN;
	    
	    if( axis == SBSGEM::kUaxis ){
	      for( int isamp=0; isamp<fN_MPD_TIME_SAMP; isamp++ ){
		// std::cout << "U axis: isamp, strip, rawADC, ADC, ped, commonmode = " << isamp << ", " << strip << ", "
		// 	  << rawADCtemp[isamp] << ", " << ADCtemp[isamp] << ", " << pedtemp << ", " << commonMode[isamp] << std::endl;
		
		hrawADCs_by_stripU->Fill( strip, rawADCtemp[isamp] );
		hpedestal_subtracted_ADCs_by_stripU->Fill( strip, ADCtemp[isamp] ); //common-mode AND ped-subtracted
		hcommonmode_subtracted_ADCs_by_stripU->Fill( strip, ADCtemp[isamp] + pedtemp ); //common-mode subtraction only, no ped:
		hpedestal_subtracted_rawADCs_by_stripU->Fill( strip, ADCtemp[isamp] + commonMode[isamp] ); //pedestal subtraction only, no common-mode

		hpedestal_subtracted_rawADCsU->Fill( ADCtemp[isamp] + commonMode[isamp] ); //1D distribution of ped-subtracted ADCs w/o common-mode subtraction
		hpedestal_subtracted_ADCsU->Fill( ADCtemp[isamp] ); //1D distribution of ped-and-common-mode subtracted ADCs

		hcommonmode_mean_by_APV_U->Fill( iAPV, commonMode[isamp] );
		// ( (TH2F*) (*hrawADCs_by_strip_sampleU)[isamp] )->Fill( strip, rawADCtemp[isamp] );
		// //for this one, we add back in the pedestal:
		// ( (TH2F*) (*hcommonmode_subtracted_ADCs_by_strip_sampleU)[isamp] )->Fill( strip, ADCtemp[isamp] + pedtemp );
		// ( (TH2F*) (*hpedestal_subtracted_ADCs_by_strip_sampleU)[isamp] )->Fill( strip, ADCtemp[isamp] );
	      }
	    } else {
	      for( int isamp=0; isamp<fN_MPD_TIME_SAMP; isamp++ ){
		// std::cout << "V axis: isamp, strip, rawADC, ADC, ped, commonmode = " << isamp << ", " << strip << ", "
		// 	  << rawADCtemp[isamp] << ", " << ADCtemp[isamp] << ", " << pedtemp << ", " << commonMode[isamp] << std::endl;
		hrawADCs_by_stripV->Fill( strip, rawADCtemp[isamp] );
		hpedestal_subtracted_ADCs_by_stripV->Fill( strip, ADCtemp[isamp] );
		hcommonmode_subtracted_ADCs_by_stripV->Fill( strip, ADCtemp[isamp] + pedtemp );
		hpedestal_subtracted_rawADCs_by_stripV->Fill( strip, ADCtemp[isamp] + commonMode[isamp] );

		hpedestal_subtracted_rawADCsV->Fill( ADCtemp[isamp] + commonMode[isamp] );
		hpedestal_subtracted_ADCsV->Fill( ADCtemp[isamp] );

		
		hcommonmode_mean_by_APV_V->Fill( iAPV, commonMode[isamp] );
		
		// ( (TH2F*) (*hrawADCs_by_strip_sampleV)[isamp] )->Fill( strip, rawADCtemp[isamp] );
		// ( (TH2F*) (*hcommonmode_subtracted_ADCs_by_strip_sampleV)[isamp] )->Fill( strip, ADCtemp[isamp] );
		// ( (TH2F*) (*hpedestal_subtracted_ADCs_by_strip_sampleV)[isamp] )->Fill( strip, ADCtemp[isamp] - pedtemp );
	      }
	    }

	    // std::cout << "finished pedestal histograms..." << std::endl;
	    
	  }
	  
	  fNstrips_hit++;
	}
      } //end loop over strips on this APV card
    } //end loop over APV cards with hits
  } //end loop on decode map entries for this module

  fNstrips_hitU = fUstripIndex.size();
  fNstrips_hitV = fVstripIndex.size();
  
  //std::cout << fName << " decoded, number of strips fired = " << fNstrips_hit << std::endl;

  fIsDecoded = true;
  
  return 0;
}

void SBSGEMModule::find_2Dhits(){ //version with no arguments calls 1D cluster finding with default (wide-open) track search constraints
  //these functions will fill the 1D cluster arrays:
  
  find_clusters_1D(SBSGEM::kUaxis); //u strips
  find_clusters_1D(SBSGEM::kVaxis); //v strips

  //Now make 2D clusters:

  if( fNclustU > 0 && fNclustV > 0 ){
  
    fxcmin = -1.e12;
    fxcmax = 1.e12;
    fycmin = -1.e12;
    fycmax = 1.e12;
  
    fill_2D_hit_arrays();

    //std::cout << "Module " << fName << ", found " << fHits.size() << " 2D cluster candidates" << std::endl;
  }
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
  find_clusters_1D(SBSGEM::kUaxis, ucenter, (umax-umin)/2.0 ); //U clusters
  find_clusters_1D(SBSGEM::kVaxis, vcenter, (vmax-vmin)/2.0 );  //V clusters

  //Now make 2D clusters
  
  fill_2D_hit_arrays();
}

void SBSGEMModule::find_clusters_1D( SBSGEM::GEMaxis_t axis, Double_t constraint_center, Double_t constraint_width ){

  //Constraint center and constraint width are assumed to be given in "standard" Hall A units (meters) in module-local coordinates
  // (SPECIFICALLY: constraint center and constraint width are assumed to refer to the direction measured by the strips being considered here)
  //This method will generally only be called by the reconstruction methods of SBSGEMTrackerBase

  // cout << "constraint center, constraint width = " << constraint_center << ", " << constraint_width << endl;
  
  if( !fIsDecoded ){
    cout << "find_clusters invoked before decoding for GEM Module " << GetName() << ", doing nothing" << endl;
    return;
  }

  UShort_t maxsep;
  UShort_t maxsepcoord; 
  UInt_t Nstrips;
  Double_t pitch;

  if( axis == SBSGEM::kUaxis ){
    maxsep = fMaxNeighborsU_totalcharge;
    maxsepcoord = fMaxNeighborsU_hitpos; 
    Nstrips = fNstripsU;
    pitch = fUStripPitch;
    fNclustU = 0;
    fUclusters.clear();
  } else { //V axis, no need to compare axis to kVaxis:
    maxsep = fMaxNeighborsV_totalcharge;
    maxsepcoord = fMaxNeighborsV_hitpos; 
    Nstrips = fNstripsV;
    pitch = fVStripPitch;
    fNclustV = 0;
    fVclusters.clear();
  }
  
  std::set<UShort_t> striplist;  //sorted list of strips for 1D clustering
  std::map<UShort_t, UInt_t> hitindex; //key = strip ID, mapped value = index in decoded hit array, needed to access the other information efficiently:
  std::map<UShort_t, Double_t> pedrms_strip;
  
  for( int ihit=0; ihit<fNstrips_hit; ihit++ ){
    if( fAxis[ihit] == axis && fKeepStrip[ihit] ){
      
      bool newstrip = (striplist.insert( fStrip[ihit] ) ).second;

      if( newstrip ){ //should always be true:
	hitindex[fStrip[ihit]] = ihit;
	if( axis == SBSGEM::kUaxis ){
	  pedrms_strip[fStrip[ihit]] = fPedRMSU[fStrip[ihit]];
	} else {
	  pedrms_strip[fStrip[ihit]] = fPedRMSV[fStrip[ihit]];
	}
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

  //cout << "before peak erasing, n local maxima = " << localmaxima.size() << endl;
  
  vector<int> peakstoerase; 

  //now calculate "prominence" for all peaks and erase "insignificant" peaks:

  for( std::set<UShort_t>::iterator i=localmaxima.begin(); i != localmaxima.end(); ++i ){
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

    //    double sigma_sum = sqrt( double(fN_MPD_TIME_SAMP) )*fZeroSuppressRMS; //~25-50 ADC

    // the above calculation is not consistent with the intended usages of the above variables.
    // We should instead use the pedestal RMS value for the strip in question. The strip RMS values
    // represent the RMS of the AVERAGE of the samples. So the prominence threshold should be expressed in terms of the same thing to be consistent:
    double sigma_sum = double(fN_MPD_TIME_SAMP)*pedrms_strip[stripmax]; //Typically 60-70 ADC. Since pedrms_strip represents the rms of the sum of the ADC samples divided by the number of samples, to get the RMS of the sum, we need only multiply by the number of samples.
    //A 5-sigma threshold on this quantity would typically be about 300-350 ADC.
    
    bool peak_close = false;
    if( !higherpeakleft ) ADCminleft = 0.0;
    if( !higherpeakright ) ADCminright = 0.0;

    if( higherpeakright || higherpeakleft ){ //this peak is contiguous with higher peaks on either the left or right or both:
      prominence = ADCmax - std::max( ADCminleft, ADCminright ); //subtract the higher of the two valleys to get the prominence

      if( higherpeakleft && std::abs( peakleft - stripmax ) <= 2*maxsep ) peak_close = true;
      if( higherpeakright && std::abs( peakright - stripmax ) <= 2*maxsep ) peak_close = true;

      if( peak_close && (prominence < fThresh_2ndMax_nsigma * sigma_sum ||
			 prominence/ADCmax < fThresh_2ndMax_fraction ) ){
	peakstoerase.push_back( stripmax );
      }
    }
  }

  //Erase "insignificant" peaks (those in contiguous grouping with higher peak with prominence below thresholds):
  for( int ipeak=0; ipeak<peakstoerase.size(); ipeak++ ){
    localmaxima.erase( peakstoerase[ipeak] );
    islocalmax[peakstoerase[ipeak]] = false;
  }

  //cout << "After peak erasing, n local maxima = " << localmaxima.size() << endl;
  
  
  //Cluster formation and cluster splitting from remaining local maxima:
  for( std::set<UShort_t>::iterator i = localmaxima.begin(); i != localmaxima.end(); ++i ){
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
	//int nmax_strip = 1;
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
	  clusttemp.ADCsamples[isamp] += fADCsamples[hitindex[istrip]][isamp]*splitfraction[istrip];
	}

	clusttemp.stripADCsum.push_back( ADCstrip );

	clusttemp.hitindex.push_back( hitindex[istrip] ); //do we use this anywhere? Yes, it is good to keep track of this if we want to access raw strip info later on 
	
	sumADC += ADCstrip;
	
	if( std::abs( istrip - stripmax ) <= std::max(UShort_t(1),std::min(maxsepcoord,maxsep)) ){ 
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

      // In the standalone we don't apply an independent threshold on the cluster sum in the context of 1D cluster-finding:
      // if( sumADC >= fThresholdClusterSum ){
      if( axis == SBSGEM::kVaxis ){
	fVclusters.push_back( clusttemp );
	fNclustV++;
      } else {
	fUclusters.push_back( clusttemp );
	fNclustU++;
      }

	// std::cout << "found cluster, (axis, istripmax, nstrips, ADCsum, hit pos (mm) )=(" << axis << ", " << clusttemp.istripmax << ", "
	// 	  << clusttemp.nstrips << ", " << clusttemp.clusterADCsum
	// 	  << ", " << clusttemp.hitpos_mean*1000.0 << ")" << std::endl;
	
      //}
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

      //"Cluster moments" defined as differences between reconstructed hit position and center of strip with max. signal in the cluster:
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
	fN2Dhits++;
      } //end check that 2D point is inside track search region
    } //end loop over "V" clusters
  } //end loop over "U" clusters

  
  
}

void    SBSGEMModule::Print( Option_t* opt) const{
  return;
}

Int_t   SBSGEMModule::Begin( THaRunBase* r){ //Does nothing
  //Here we can create some histograms that will be written to the ROOT file:
  //This is a natural place to do the hit maps/efficiency maps:

  //pulled these lines out of the if-block below to avoid code duplication:
  TString histname,histtitle;
  TString detname = GetParent()->GetName();
  detname.ReplaceAll(".","_");
  detname += "_";
  detname += GetName();
  
  if( fMakeEfficiencyPlots ){
    
  
    fhdidhitx = new TH1F( histname.Format( "hdidhitx_%s", detname.Data() ), "local x coordinate of hits on good tracks (m)", 100, -0.51*GetXSize(), 0.51*GetXSize() );
    fhdidhity = new TH1F( histname.Format( "hdidhity_%s", detname.Data() ), "local y coordinate of hits on good tracks (m)", 100, -0.51*GetYSize(), 0.51*GetYSize() );
    fhdidhitxy = new TH2F( histname.Format( "hdidhitxy_%s", detname.Data() ), "x vs y of hits on good tracks (m)",
			   100, -0.51*GetYSize(), 0.51*GetYSize(),
			   100, -0.51*GetXSize(), 0.51*GetXSize() );

    fhshouldhitx = new TH1F( histname.Format( "hshouldhitx_%s", detname.Data() ), "x of good track passing through (m)", 100, -0.51*GetXSize(), 0.51*GetXSize() );
    fhshouldhity = new TH1F( histname.Format( "hshouldhity_%s", detname.Data() ), "y of good track passing through (m)", 100, -0.51*GetYSize(), 0.51*GetYSize() );
    fhshouldhitxy = new TH2F( histname.Format( "hshouldhitxy_%s", detname.Data() ), "x vs y of good track passing through (m)",
			      100, -0.51*GetYSize(), 0.51*GetYSize(),
			      100, -0.51*GetXSize(), 0.51*GetXSize() );			
  }

  if( fPedestalMode ){ //make pedestal histograms:

    //Procedure:
    // 1. Analyze pedestal data with no pedestals supplied from the database.
    // 2. Extract "coarse" strip offsets from raw ADC histograms (i.e., common-mode means)
    // 3. Extract "fine" strip offsets from common-mode-subtracted histograms 
    // 3. Add "coarse" strip offsets (common-mode means) back into "common-mode subtracted" ADC 
    
    hrawADCs_by_stripU = new TH2F( histname.Format( "hrawADCs_by_stripU_%s", detname.Data() ), "Raw ADCs by U strip number, no corrections",
				   fNstripsU, -0.5, fNstripsU-0.5,
				   2048, -0.5, 4095.5 );
    hrawADCs_by_stripV = new TH2F( histname.Format( "hrawADCs_by_stripV_%s", detname.Data() ), "Raw ADCs by V strip number, no corrections",
				   fNstripsV, -0.5, fNstripsV-0.5,
				   2048, -0.5, 4095.5 );

    hcommonmode_subtracted_ADCs_by_stripU = new TH2F( histname.Format( "hpedestalU_%s", detname.Data() ), "ADCs by U strip number, w/common mode correction, no ped. subtraction",
						      fNstripsU, -0.5, fNstripsU-0.5,
						      2500, -500.0, 4500.0 );
    hcommonmode_subtracted_ADCs_by_stripV = new TH2F( histname.Format( "hpedestalV_%s", detname.Data() ), "ADCs by V strip number, w/common mode correction, no ped. subtraction",
						      fNstripsV, -0.5, fNstripsV-0.5,
						      2500, -500.0, 4500.0 );

    hpedestal_subtracted_ADCs_by_stripU = new TH2F( histname.Format( "hADCpedsubU_%s", detname.Data() ), "Pedestal and common-mode subtracted ADCs by U strip number",
						    fNstripsU, -0.5, fNstripsU-0.5,
						    1000,-500.,500. );
    hpedestal_subtracted_ADCs_by_stripV = new TH2F( histname.Format( "hADCpedsubV_%s", detname.Data() ), "Pedestal and common-mode subtracted ADCs by V strip number",
						    fNstripsV, -0.5, fNstripsV-0.5,
						    1000,-500.,500. );

    hpedestal_subtracted_rawADCs_by_stripU = new TH2F( histname.Format( "hrawADCpedsubU_%s", detname.Data() ), "ADCs by U strip, ped-subtracted, no common-mode correction",
						       fNstripsU, -0.5, fNstripsU-0.5,
						       2500,-500.,4500. );
    hpedestal_subtracted_rawADCs_by_stripV = new TH2F( histname.Format( "hrawADCpedsubV_%s", detname.Data() ), "ADCs by V strip, ped-subtracted, no common-mode correction",
						       fNstripsV, -0.5, fNstripsV-0.5,
						       2500,-500.,4500. );

    hpedestal_subtracted_rawADCsU = new TH1F( histname.Format( "hrawADCpedsubU_allstrips_%s", detname.Data() ), "distribution of ped-subtracted U strip ADCs w/o common-mode correction",
					      2500, -500.,4500. );
    hpedestal_subtracted_rawADCsV = new TH1F( histname.Format( "hrawADCpedsubV_allstrips_%s", detname.Data() ), "distribution of ped-subtracted V strip ADCs w/o common-mode correction",
					      2500, -500.,4500. );

    hpedestal_subtracted_ADCsU = new TH1F( histname.Format( "hADCpedsubU_allstrips_%s", detname.Data() ), "distribution of ped-subtracted U strip ADCs w/common-mode correction",
					      1000, -500.,500. );
    hpedestal_subtracted_ADCsV = new TH1F( histname.Format( "hADCpedsubV_allstrips_%s", detname.Data() ), "distribution of ped-subtracted V strip ADCs w/common-mode correction",
					      1000, -500.,500. );

    int nAPVs_U = fNstripsU/fN_APV25_CHAN;
    hcommonmode_mean_by_APV_U = new TH2F( histname.Format( "hCommonModeMean_by_APV_U_%s", detname.Data() ), "distribution of common-mode means for U strip pedestal data",
					  nAPVs_U, -0.5, nAPVs_U-0.5,  
					  2048, -0.5, 4095.5 );
    int nAPVs_V = fNstripsV/fN_APV25_CHAN;
    hcommonmode_mean_by_APV_V = new TH2F( histname.Format( "hCommonModeMean_by_APV_V_%s", detname.Data() ), "distribution of common-mode means for V strip pedestal data",
					  nAPVs_V, -0.5, nAPVs_V-0.5,
					  2048, -0.5, 4095.5 );
    
    // Uncomment these later if you want them:
    // hrawADCs_by_strip_sampleU = new TClonesArray( "TH2F", fN_MPD_TIME_SAMP );
    // hrawADCs_by_strip_sampleV = new TClonesArray( "TH2F", fN_MPD_TIME_SAMP );

    // hcommonmode_subtracted_ADCs_by_strip_sampleU = new TClonesArray( "TH2F", fN_MPD_TIME_SAMP );
    // hcommonmode_subtracted_ADCs_by_strip_sampleV = new TClonesArray( "TH2F", fN_MPD_TIME_SAMP );

    // hpedestal_subtracted_ADCs_by_strip_sampleU = new TClonesArray( "TH2F", fN_MPD_TIME_SAMP );
    // hpedestal_subtracted_ADCs_by_strip_sampleV = new TClonesArray( "TH2F", fN_MPD_TIME_SAMP );
    
    // for( int isamp = 0; isamp<fN_MPD_TIME_SAMP; isamp++ ){
    //   new( (*hrawADCs_by_strip_sampleU)[isamp] ) TH2F( histname.Format( "hrawADCU_%s_sample%d", detname.Data(), isamp ), "Raw U ADCs by strip and sample",
    // 						       fNstripsU, -0.5, fNstripsU-0.5,
    // 						       1024, -0.5, 4095.5 );
    //   new( (*hrawADCs_by_strip_sampleV)[isamp] ) TH2F( histname.Format( "hrawADCV_%s_sample%d", detname.Data(), isamp ), "Raw V ADCs by strip and sample",
    // 						       fNstripsV, -0.5, fNstripsV-0.5,
    // 						       1024, -0.5, 4095.5 );

    //   new( (*hcommonmode_subtracted_ADCs_by_strip_sampleU)[isamp] ) TH2F( histname.Format( "hpedestalU_%s_sample%d", detname.Data(), isamp ), "Pedestals by strip and sample",
    // 									  fNstripsU, -0.5, fNstripsU-0.5,
    // 									  1500, -500.0, 1000.0 );
    //   new( (*hcommonmode_subtracted_ADCs_by_strip_sampleV)[isamp] ) TH2F( histname.Format( "hpedestalV_%s_sample%d", detname.Data(), isamp ), "Pedestals by strip and sample",
    // 									  fNstripsV, -0.5, fNstripsV-0.5,
    // 									  1500, -500.0, 1000.0 );

    //   new( (*hpedestal_subtracted_ADCs_by_strip_sampleU)[isamp] ) TH2F( histname.Format( "hADCpedsubU_%s_sample%d", detname.Data(), isamp ), "Pedestal-subtracted ADCs by strip and sample",
    // 									fNstripsU, -0.5, fNstripsU-0.5,
    // 									1000, -500.0, 500.0 );
    //   new( (*hpedestal_subtracted_ADCs_by_strip_sampleV)[isamp] ) TH2F( histname.Format( "hADCpedsubV_%s_sample%d", detname.Data(), isamp ), "Pedestal-subtracted ADCs by strip and sample",
    // 									fNstripsV, -0.5, fNstripsV-0.5,
    // 									1000, -500.0, 500.0 );
      
    // }
    
  }
  
    
  return 0;
}

void SBSGEMModule::PrintPedestals( std::ofstream &dbfile, std::ofstream &daqfile, std::ofstream &daqfile_cmr ){
  //The first argument is a file in the format expected by the database,
  //The second argument is a file in the format expected by the DAQ:

  //Step 1: we need to extract pedestal mean and rms by channel, and also grab
  // "common-mode min" and "common-mode max" values:
  // we can simply use the existing arrays to store the values:

  std::cout << "[SBSGEMModule::PrintPedestals]: module " << GetName() << std::endl;
  
  //start with pedestals;
  for( int iu = 0; iu<fNstripsU; iu++ ){
    TH1D *htemp = hcommonmode_subtracted_ADCs_by_stripU->ProjectionY( "htemp", iu+1,iu+1 );
    
    fPedestalU[iu] = htemp->GetMean();

    //htemp here represents the individual sample noise width, but the threshold is applied on the average of the six (or other number) time samples:
    // sigma(average) = sigma(1 sample)/sqrt(number of samples):
    fPedRMSU[iu] = htemp->GetRMS() / sqrt( double( fN_MPD_TIME_SAMP ) ); 

    htemp->Delete();
  }

  for( int iv = 0; iv<fNstripsV; iv++ ){
    TH1D *htemp = hcommonmode_subtracted_ADCs_by_stripV->ProjectionY( "htemp", iv+1, iv+1 );
    fPedestalV[iv] = htemp->GetMean();
    fPedRMSV[iv] = htemp->GetRMS() / sqrt( double( fN_MPD_TIME_SAMP ) );

    htemp->Delete();
  }

  
  
  TString appname = static_cast<THaDetector *>(GetParent())->GetApparatus()->GetName();
  TString detname = GetParent()->GetName();
  TString modname = GetName();
  
  TString header;
  header.Form( "%s.%s.%s.pedu = ", appname.Data(), detname.Data(), modname.Data() );
  dbfile << std::endl << header << std::endl;
  for( int iu=0; iu<fNstripsU; iu++ ){
    TString sentry;
    sentry.Form( "  %15.5g ", fPedestalU[iu] );
    dbfile << sentry;
    if( (iu+1)%16 == 0 ) dbfile << std::endl;
  }

  dbfile << std::endl;
  
  header.Form( "%s.%s.%s.rmsu = ", appname.Data(), detname.Data(), modname.Data() );
  dbfile << std::endl << header << std::endl;

  for( int iu=0; iu<fNstripsU; iu++ ){
    TString sentry;
    sentry.Form( "  %15.5g ", fPedRMSU[iu] );
    dbfile << sentry;
    if( (iu+1)%16 == 0 ) dbfile << std::endl;
  }

  dbfile << std::endl;

  header.Form( "%s.%s.%s.pedv = ", appname.Data(), detname.Data(), modname.Data() );
  dbfile << std::endl << header << std::endl;
  for( int iv=0; iv<fNstripsV; iv++ ){
    TString sentry;
    sentry.Form( "  %15.5g ", fPedestalV[iv] );
    dbfile << sentry;
    if( (iv+1)%16 == 0 ) dbfile << std::endl;
  }

  dbfile << std::endl;
  
  header.Form( "%s.%s.%s.rmsv = ", appname.Data(), detname.Data(), modname.Data() );
  dbfile << std::endl << header << std::endl;

  for( int iv=0; iv<fNstripsV; iv++ ){
    TString sentry;
    sentry.Form( "  %15.5g ", fPedRMSV[iv] );
    dbfile << sentry;
    if( (iv+1)%16 == 0 ) dbfile << std::endl;
  }

  dbfile << std::endl;
  
  //TO DO: common-mode mean, min, and max by APV card:
  int nAPVsU = fNstripsU/fN_APV25_CHAN;
  int nAPVsV = fNstripsV/fN_APV25_CHAN;
  std::vector<double> commonmode_meanU(nAPVsU), commonmode_rmsU(nAPVsU);
  std::vector<double> commonmode_meanV(nAPVsV), commonmode_rmsV(nAPVsV);

  for( int iAPV = 0; iAPV<nAPVsU; iAPV++ ){
    TH1D *htemp = hcommonmode_mean_by_APV_U->ProjectionY("htemp", iAPV+1, iAPV+1 );

    commonmode_meanU[iAPV] = htemp->GetMean();
    commonmode_rmsU[iAPV] = htemp->GetRMS();
  }

  for( int iAPV = 0; iAPV<nAPVsV; iAPV++ ){
    TH1D *htemp = hcommonmode_mean_by_APV_V->ProjectionY("htemp", iAPV+1, iAPV+1 );

    commonmode_meanV[iAPV] = htemp->GetMean();
    commonmode_rmsV[iAPV] = htemp->GetRMS();
  }
  
  //That takes care of the database file. For the "DAQ" file, we need to organize things by APV card. For this we can loop over the MPDmap:
  
  for( auto iapv = fMPDmap.begin(); iapv != fMPDmap.end(); iapv++ ){
    int crate = iapv->crate;
    int mpd = iapv->mpd_id;
    int adc_ch = iapv->adc_id;
    int pos = iapv->pos;
    int invert = iapv->invert;
    int axis = iapv->axis;

    daqfile << "APV "
	    << std::setw(16) << crate
	    << std::setw(16) << mpd
	    << std::setw(16) << adc_ch
	    << std::endl;
    
    //Then loop over the strips:
    for( int ich=0; ich<fN_APV25_CHAN; ich++ ){
      int strip = GetStripNumber( ich, pos, invert );
      double pedmean = (axis == SBSGEM::kUaxis ) ? fPedestalU[strip] : fPedestalV[strip];
      double pedrms = (axis == SBSGEM::kUaxis ) ? fPedRMSU[strip] : fPedRMSV[strip];

      daqfile << std::setw(16) << ich
	      << std::setw(16) << std::setprecision(4) << pedmean
	      << std::setw(16) << std::setprecision(4) << pedrms
	      << std::endl;
    }

    double cm_mean = (axis == SBSGEM::kUaxis ) ? commonmode_meanU[ pos ] : commonmode_meanV[ pos ];
    double cm_rms = (axis == SBSGEM::kUaxis ) ? commonmode_rmsU[ pos ] : commonmode_rmsV[ pos ];
    
    double cm_min = cm_mean - fZeroSuppressRMS * cm_rms;
    double cm_max = cm_mean + fZeroSuppressRMS * cm_rms;
    
    daqfile_cmr << std::setw(12) << crate
		<< std::setw(12) << mpd
		<< std::setw(12) << adc_ch
		<< std::setw(12) << int( cm_min )
		<< std::setw(12) << int( cm_max )
		<< std::endl;
		
    
  }

  

  
  
}

Int_t   SBSGEMModule::End( THaRunBase* r){ //Calculates efficiencies and writes hit maps and efficiency histograms and/or pedestal info to ROOT file:

  if( fMakeEfficiencyPlots ){
    //Create the track-based efficiency histograms at the end of the run:
    TString histname;
    TString detname = GetParent()->GetName();
    detname.ReplaceAll(".","_");
    detname += "_";
    detname += GetName();
  
    if( fhdidhitx != NULL && fhshouldhitx != NULL ){ //Create efficiency histograms and write to the ROOT file:
      TH1F *hefficiency_vs_x = new TH1F(*fhdidhitx);
      hefficiency_vs_x->SetName( histname.Format( "hefficiency_vs_x_%s", detname.Data() ) );
      hefficiency_vs_x->SetTitle( histname.Format( "Track-based efficiency vs x, module %s", GetName() ) );
      hefficiency_vs_x->Divide( fhshouldhitx );
      hefficiency_vs_x->Write();
    }

    if( fhdidhity != NULL && fhshouldhity != NULL ){ //Create efficiency histograms and write to the ROOT file:
      TH1F *hefficiency_vs_y = new TH1F(*fhdidhity);
      hefficiency_vs_y->SetName( histname.Format( "hefficiency_vs_y_%s", detname.Data() ) );
      hefficiency_vs_y->SetTitle( histname.Format( "Track-based efficiency vs y, module %s", GetName() ) );
      hefficiency_vs_y->Divide( fhshouldhity );
      hefficiency_vs_y->Write();
    }

    if( fhdidhitxy != NULL && fhshouldhitxy != NULL ){ //Create efficiency histograms and write to the ROOT file:
      TH2F *hefficiency_vs_xy = new TH2F(*fhdidhitxy);
      hefficiency_vs_xy->SetName( histname.Format( "hefficiency_vs_xy_%s", detname.Data() ) );
      hefficiency_vs_xy->SetTitle( histname.Format( "Track-based efficiency vs x and y, module %s", GetName() ) );
      hefficiency_vs_xy->Divide( fhshouldhitxy );
      hefficiency_vs_xy->Write();
    }
  
    if( fhdidhitx != NULL  ) fhdidhitx->Write();
    if( fhdidhity != NULL  ) fhdidhity->Write();
    if( fhdidhitxy != NULL  ) fhdidhitxy->Write();

    if( fhshouldhitx != NULL  ) fhshouldhitx->Write();
    if( fhshouldhity != NULL  ) fhshouldhity->Write();
    if( fhshouldhitxy != NULL  ) fhshouldhitxy->Write();
  }

  if( fPedestalMode ){ //write out pedestal histograms, print out pedestals in the format needed for both database and DAQ:

    //The channel map and pedestal information are specific to a GEM module.
    //But we want ONE database file and ONE DAQ file for the entire tracker.
    //So probably the best way to proceed is to have the  
    
    hrawADCs_by_stripU->Write();
    hrawADCs_by_stripV->Write();
    hcommonmode_subtracted_ADCs_by_stripU->Write();
    hcommonmode_subtracted_ADCs_by_stripV->Write();
    hpedestal_subtracted_ADCs_by_stripU->Write();
    hpedestal_subtracted_ADCs_by_stripV->Write();
    hpedestal_subtracted_rawADCs_by_stripU->Write();
    hpedestal_subtracted_rawADCs_by_stripV->Write();

    hpedestal_subtracted_rawADCsU->Write();
    hpedestal_subtracted_rawADCsV->Write();
    hpedestal_subtracted_ADCsU->Write();
    hpedestal_subtracted_ADCsV->Write();

    hcommonmode_mean_by_APV_U->Write();
    hcommonmode_mean_by_APV_V->Write();

    // hrawADCs_by_strip_sampleU->Write();
    // hrawADCs_by_strip_sampleV->Write();
    // hcommonmode_subtracted_ADCs_by_strip_sampleU->Write();
    // hcommonmode_subtracted_ADCs_by_strip_sampleV->Write();
    // hpedestal_subtracted_ADCs_by_strip_sampleU->Write();
    // hpedestal_subtracted_ADCs_by_strip_sampleV->Write();
  }
  
    
  return 0;
}

//utility method to calculate correlation coefficient of U and V samples: 
Double_t SBSGEMModule::CorrCoeff( int nsamples, std::vector<double> Usamples, std::vector<double> Vsamples ){
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

Int_t SBSGEMModule::GetStripNumber( UInt_t rawstrip, UInt_t pos, UInt_t invert ){
  Int_t RstripNb = APVMAP[rawstrip];
  RstripNb = RstripNb + (127-2*RstripNb)*invert;
  Int_t RstripPos = RstripNb + 128*pos;

  if( fIsMC ){
    return rawstrip + 128*pos;
  }
  
  return RstripPos;
}
