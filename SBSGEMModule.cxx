#include <iostream>

#include "SBSGEMModule.h"
#include "TDatime.h"
#include "THaEvData.h"
#include "THaApparatus.h"
#include "TRotation.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TClonesArray.h"
#include <algorithm>

using namespace std;
//using namespace SBSGEMModule;

//This should not be hard-coded, I think, but read in from the database (or perhaps not, if it never changes? For now we keep it hard-coded)
// const int APVMAP[128] = {1, 33, 65, 97, 9, 41, 73, 105, 17, 49, 81, 113, 25, 57, 89, 121, 3, 35, 67, 99, 11, 43, 75, 107, 19, 51, 83, 115, 27, 59, 91, 123, 5, 37, 69, 101, 13, 45, 77, 109, 21, 53, 85, 117, 29, 61, 93, 125, 7, 39, 71, 103, 15, 47, 79, 111, 23, 55, 87, 119, 31, 63, 95, 127, 0, 32, 64, 96, 8, 40, 72, 104, 16, 48, 80, 112, 24, 56, 88, 120, 2, 34, 66, 98, 10, 42, 74, 106, 18, 50, 82, 114, 26, 58, 90, 122, 4, 36, 68, 100, 12, 44, 76, 108, 20, 52, 84, 116, 28, 60, 92, 124, 6, 38, 70, 102, 14, 46, 78, 110, 22, 54, 86, 118, 30, 62, 94, 126};


SBSGEMModule::SBSGEMModule( const char *name, const char *description,
			    THaDetectorBase* parent ):
  THaSubDetector(name,description,parent)
{
  // FIXME:  To database
  //Set Default values for fZeroSuppress and fZeroSuppressRMS:
  fZeroSuppress    = kTRUE;
  fZeroSuppressRMS = 5.0; //threshold in units of RMS:

  fPedestalMode = kFALSE;
  fPedHistosInitialized = kFALSE;
  
  //Default online zero suppression to FALSE: actually I wonder if it would be better to use this in 
  // Moved this to MPDModule, since this should be done during the decoding of the raw APV data:
  fOnlineZeroSuppression = kFALSE;

  fCommonModeFlag = 0; //"sorting" method
  //Default: discard highest and lowest 28 strips for "sorting method" common-mode calculation:
  fCommonModeNstripRejectHigh = 28; 
  fCommonModeNstripRejectLow = 28;

  fCommonModeNumIterations = 3;
  fCommonModeMinStripsInRange = 10;
  fMakeCommonModePlots = false;
  fCommonModePlotsInitialized = false;
  
  fPedSubFlag = 1; //default to online ped subtraction, as that is the mode we will run in most of the time

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

  //Default clustering parameters:
  fThresholdSample = 50.0;
  fThresholdStripSum = 250.0;
  fThresholdClusterSum = 500.0;

  fADCasymCut = 1.1;
  fTimeCutUVdiff = 30.0;
  fCorrCoeffCut = 0.5;

  fFiltering_flag1D = 0; //"soft" cuts
  fFiltering_flag2D = 0; //"soft" cuts
  
  // default these offsets to zero: 
  fUStripOffset = 0.0;
  fVStripOffset = 0.0;
  
  fMakeEfficiencyPlots = true;
  fEfficiencyInitialized = false;

  fChan_CM_flags = 512; //default to 512:
  fChan_TimeStamp_low = 513;
  fChan_TimeStamp_high = 514;
  fChan_MPD_EventCount = 515;

  
  UInt_t MAXNSAMP_PER_APV = fN_APV25_CHAN * fN_MPD_TIME_SAMP;
  //arrays to hold raw data from one APV card:
  fStripAPV.resize( MAXNSAMP_PER_APV );
  fRawStripAPV.resize( MAXNSAMP_PER_APV );
  fRawADC_APV.resize( MAXNSAMP_PER_APV );
  fPedSubADC_APV.resize( MAXNSAMP_PER_APV );
  fCommonModeSubtractedADC_APV.resize( MAXNSAMP_PER_APV );

  //default to 
  fMAX2DHITS = 250000;

  fRMS_ConversionFactor = sqrt(fN_MPD_TIME_SAMP); //=2.45

  fIsMC = false; //need to set default value!

  fAPVmapping = SBSGEM::kUVA_XY; //default to UVA X/Y style APV mapping, but require this in the database::

  InitAPVMAP();
  
  //  std::cout << "SBSGEMModule constructor invoked, name = " << name << std::endl;
  
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

  // default these offsets to zero: 
  fUStripOffset = 0.0;
  fVStripOffset = 0.0;

  //std::cout << "Before loadDB, fCommonModePlotsInitialized = " << fCommonModePlotsInitialized << std::endl;

  int cmplots_flag = fMakeCommonModePlots ? 1 : 0;
  int zerosuppress_flag = fZeroSuppress ? 1 : 0;
  int onlinezerosuppress_flag = fOnlineZeroSuppression ? 1 : 0;
  

  const DBRequest request[] = {
    { "chanmap",        &fChanMapData,        kIntV, 0, 0, 0}, // mandatory: decode map info
    { "apvmap",         &fAPVmapping,    kUInt, 0, 1, 1}, //optional, allow search up the tree if all modules in a setup have the same APV mapping
    { "pedu",           &rawpedu,        kDoubleV, 0, 1, 0}, // optional raw pedestal info (u strips)
    { "pedv",           &rawpedv,        kDoubleV, 0, 1, 0}, // optional raw pedestal info (v strips)
    { "rmsu",           &rawrmsu,        kDoubleV, 0, 1, 0}, // optional pedestal rms info (u strips)
    { "rmsv",           &rawrmsv,        kDoubleV, 0, 1, 0}, // optional pedestal rms info (v strips)
    { "layer",          &fLayer,         kUShort, 0, 0, 0}, // mandatory: logical tracking layer must be specified for every module:
    { "nstripsu",       &fNstripsU,     kUInt, 0, 0, 1}, //mandatory: number of strips in module along U axis
    { "nstripsv",       &fNstripsV,     kUInt, 0, 0, 1}, //mandatory: number of strips in module along V axis
    { "uangle",         &fUAngle,       kDouble, 0, 0, 1}, //mandatory: Angle of "U" strips wrt X axis
    { "vangle",         &fVAngle,       kDouble, 0, 0, 1}, //mandatory: Angle of "V" strips wrt X axis
    { "uoffset",        &fUStripOffset, kDouble, 0, 1, 1}, //optional: position of first U strip
    { "voffset",        &fVStripOffset, kDouble, 0, 1, 1}, //optional: position of first V strip
    { "upitch",         &fUStripPitch,  kDouble, 0, 0, 1}, //mandatory: Pitch of U strips
    { "vpitch",         &fVStripPitch,  kDouble, 0, 0, 1}, //mandatory: Pitch of V strips
    { "ugain",          &fUgain,        kDoubleV, 0, 1, 0}, //(optional): Gain of U strips by APV card (ordered by strip position, NOT by order of appearance in decode map)
    { "vgain",          &fVgain,        kDoubleV, 0, 1, 0}, //(optional): Gain of V strips by APV card (ordered by strip position, NOT by order of appearance in decode map)
    { "threshold_sample",  &fThresholdSample, kDouble, 0, 1, 1}, //(optional): threshold on max. ADC sample to keep strip (baseline-subtracted)
    { "threshold_stripsum", &fThresholdStripSum, kDouble, 0, 1, 1}, //(optional): threshold on sum of ADC samples on a strip (baseline-subtracted)
    { "threshold_clustersum", &fThresholdClusterSum, kDouble, 0, 1, 1}, //(optional): threshold on sum of all ADCs over all strips in a cluster (baseline-subtracted)
    { "ADCasym_cut", &fADCasymCut, kDouble, 0, 1, 1}, //(optional): filter 2D hits by ADC asymmetry, |Asym| < cut
    { "deltat_cut", &fTimeCutUVdiff, kDouble, 0, 1, 1}, //(optional): filter 2D hits by U/V time difference
    { "corrcoeff_cut", &fCorrCoeffCut, kDouble, 0, 1, 1},
    { "filterflag1D", &fFiltering_flag1D, kInt, 0, 1, 1},
    { "filterflag2D", &fFiltering_flag2D, kInt, 0, 1, 1},
    { "peakprominence_minsigma", &fThresh_2ndMax_nsigma, kDouble, 0, 1, 1}, //(optional): reject overlapping clusters with peak prominence less than this number of sigmas
    { "peakprominence_minfraction", &fThresh_2ndMax_fraction, kDouble, 0, 1, 1}, //(optional): reject overlapping clusters with peak prominence less than this fraction of height of nearby higher peak
    { "maxnu_charge", &fMaxNeighborsU_totalcharge, kUShort, 0, 1, 1}, //(optional): cluster size restriction along U for total charge calculation
    { "maxnv_charge", &fMaxNeighborsV_totalcharge, kUShort, 0, 1, 1}, //(optional): cluster size restriction along V for total charge calculation
    { "maxnu_pos", &fMaxNeighborsU_hitpos, kUShort, 0, 1, 1}, //(optional): cluster size restriction for position reconstruction
    { "maxnv_pos", &fMaxNeighborsV_hitpos, kUShort, 0, 1, 1}, //(optional): cluster size restriction for position reconstruction
    { "sigmahitshape", &fSigma_hitshape, kDouble, 0, 1, 1}, //(optional): width parameter for cluster-splitting algorithm
    { "zerosuppress", &zerosuppress_flag, kUInt, 0, 1, 1}, //(optional, search): toggle offline zero suppression (default = true).
    { "zerosuppress_nsigma", &fZeroSuppressRMS, kDouble, 0, 1, 1}, //(optional, search):
    { "onlinezerosuppress", &onlinezerosuppress_flag, kUInt, 0, 1, 1}, //(optional, search)
    { "commonmode_meanU", &fCommonModeMeanU, kDoubleV, 0, 1, 0}, //(optional, don't search)
    { "commonmode_meanV", &fCommonModeMeanV, kDoubleV, 0, 1, 0}, //(optional, don't search)
    { "commonmode_rmsU", &fCommonModeRMSU, kDoubleV, 0, 1, 0}, //(optional, don't search)
    { "commonmode_rmsV", &fCommonModeRMSV, kDoubleV, 0, 1, 0}, //(optional, don't search)
    { "commonmode_flag", &fCommonModeFlag, kInt, 0, 1, 1}, //optional, search up the tree
    { "commonmode_nstriplo", &fCommonModeNstripRejectLow, kInt, 0, 1, 1}, //optional, search up the tree:
    { "commonmode_nstriphi", &fCommonModeNstripRejectHigh, kInt, 0, 1, 1}, //optional, search:
    { "commonmode_niter", &fCommonModeNumIterations, kInt, 0, 1, 1},
    { "commonmode_minstrips", &fCommonModeMinStripsInRange, kInt, 0, 1, 1},
    { "plot_common_mode", &cmplots_flag, kInt, 0, 1, 1},
    { "chan_cm_flags", &fChan_CM_flags, kUInt, 0, 1, 1}, //optional, search up the tree: must match the value in crate map!
    { "chan_timestamp_low", &fChan_TimeStamp_low, kUInt, 0, 1, 1},
    { "chan_timestamp_high", &fChan_TimeStamp_high, kUInt, 0, 1, 1},
    { "chan_event_count", &fChan_MPD_EventCount, kUInt, 0, 1, 1},
    { "pedsub_online", &fPedSubFlag, kInt, 0, 1, 1},
    { "max2Dhits", &fMAX2DHITS, kUInt, 0, 1, 1}, //optional, search up tree
    {0}
  };
  status = LoadDB( file, date, request, fPrefix, 1 ); //The "1" after fPrefix means search up the tree

  if( status != 0 ){
    fclose(file);
    return status;
  }

  fMakeCommonModePlots = cmplots_flag != 0;
  fZeroSuppress = zerosuppress_flag != 0;
  fOnlineZeroSuppression = onlinezerosuppress_flag != 0;

  //  std::cout << "After loadDB, fCommonModePlotsInitialized = " << fCommonModePlotsInitialized << std::endl;
 
  if( fAPVmapping < SBSGEM::kINFN || fAPVmapping > SBSGEM::kMC ) {
    std::cout << "Warning in SBSGEMModule::Decode for module " << GetParent()->GetName() << "." << GetName() << ": invalid APV mapping choice, defaulting to UVA X/Y." << std::endl
	      << " Analysis results may be incorrect" << std::endl;
    fAPVmapping = SBSGEM::kUVA_XY;
  }
  
  //prevent the user from defining something silly for the common-mode stuff:
  fCommonModeNstripRejectLow = std::min( 50, std::max( 0, fCommonModeNstripRejectLow ) );
  fCommonModeNstripRejectHigh = std::min( 50, std::max( 0, fCommonModeNstripRejectHigh ) );
  fCommonModeNumIterations = std::min( 10, std::max( 2, fCommonModeNumIterations ) );
  fCommonModeMinStripsInRange = std::min( fN_APV25_CHAN-25, std::max(1, fCommonModeMinStripsInRange ) );

  
  //std::cout << GetName() << " fThresholdStripSum " << fThresholdStripSum 
  //<< " fThresholdSample " << fThresholdSample << std::endl;
  
  if( fIsMC ){
    fCommonModeFlag = -1;
    fPedestalMode = false;
    fOnlineZeroSuppression = true;
    fAPVmapping = SBSGEM::kMC;
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

  //Count APV cards by axis. Each APV card must have one decode map entry:
  fNAPVs_U = 0;
  fNAPVs_V = 0;

  fMPDmap.clear();
  
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

    if( thisdata.axis == SBSGEM::kUaxis ){
      fNAPVs_U++;
    } else {
      fNAPVs_V++;
    }
    
    fMPDmap.push_back(thisdata);

    fT0_by_APV.push_back( 0 );
    fTcoarse_by_APV.push_back( 0 );
    fTfine_by_APV.push_back( 0 );
  }

  //if a different number of decode map entries is counted than the expectation based on the number of strips, 
  //e.g., because we commented out one or more decode map entries, then we go with the larger of the two numbers. 
  //This isn't perfectly idiot-proof, but prevents the kind of undefined behavior we want to avoid:
  //if( fNAPVs_U != fNstripsU/fN_APV25_CHAN ){
  fNAPVs_U = std::max( fNAPVs_U, fNstripsU/fN_APV25_CHAN );
    //}
  fNAPVs_V = std::max( fNAPVs_V, fNstripsV/fN_APV25_CHAN );


  //resize vectors that hold APV-card specific parameters:
  // fUgain.resize( fNAPVs_U );
  // fVgain.resize( fNAPVs_V );
  // fCommonModeMeanU.resize( fNAPVs_U );
  // fCommonModeMeanV.resize( fNAPVs_V );
  // fCommonModeRMSU.resize( fNAPVs_U );
  // fCommonModeRMSV.resize( fNAPVs_V );
  
  std::cout << fName << " mapped to " << nentry << " APV25 chips" << std::endl;

  // fT0_by_APV.resize( nentry );
  // fTcoarse_by_APV.resize( nentry );
  // fTfine_by_APV.resize( nentry );
  
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

  // //resize all the "decoded strip" arrays to their maximum possible values for this module:
  UInt_t nstripsmax = fNstripsU + fNstripsV;
  
  fStrip.resize( nstripsmax );
  fAxis.resize( nstripsmax );
  fADCsamples.resize( nstripsmax );
  fRawADCsamples.resize( nstripsmax );
  //The lines below are problematic and unnecessary
  for( int istrip=0; istrip<nstripsmax; istrip++ ){
    fADCsamples[istrip].resize( fN_MPD_TIME_SAMP );
    fRawADCsamples[istrip].resize( fN_MPD_TIME_SAMP );
  }
  fADCsums.resize( nstripsmax );
  fStripADCavg.resize( nstripsmax );
  fStripIsU.resize( nstripsmax );
  fStripIsV.resize( nstripsmax );
  fKeepStrip.resize( nstripsmax );
  fMaxSamp.resize( nstripsmax );
  fADCmax.resize( nstripsmax );
  fTmean.resize( nstripsmax );
  fTsigma.resize( nstripsmax );
  fTcorr.resize( nstripsmax );

  fADCsamples1D.resize( nstripsmax * fN_MPD_TIME_SAMP );
  fRawADCsamples1D.resize( nstripsmax * fN_MPD_TIME_SAMP );
  fStripTrackIndex.resize( nstripsmax );
  
  
  //default all common-mode mean and RMS values to 0 and 10 respectively if they were
  // NOT loaded from the DB and/or they are loaded with the wrong size:
  if( fCommonModeMeanU.size() != fNAPVs_U ){
    fCommonModeMeanU.resize( fNAPVs_U );
    for( int iAPV=0; iAPV<fNAPVs_U; iAPV++ ){
      fCommonModeMeanU[iAPV] = 0.0;
    }
  }

  if( fCommonModeRMSU.size() != fNAPVs_U ){
    fCommonModeRMSU.resize( fNAPVs_U );
    for( int iAPV=0; iAPV<fNAPVs_U; iAPV++ ){
      fCommonModeRMSU[iAPV] = 10.0;
    }
  }

  //default all common-mode mean and RMS values to 0 and 10 respectively if they were
  // NOT loaded from the DB and/or they were loaded with the wrong size:
  if( fCommonModeMeanV.size() != fNAPVs_V ){
    fCommonModeMeanV.resize( fNAPVs_V );
    for( int iAPV=0; iAPV<fNAPVs_V; iAPV++ ){
      fCommonModeMeanV[iAPV] = 0.0;
    }
  }

  if( fCommonModeRMSV.size() != fNAPVs_V ){
    fCommonModeRMSV.resize( fNAPVs_V );
    for( int iAPV=0; iAPV<fNAPVs_V; iAPV++ ){
      fCommonModeRMSV[iAPV] = 10.0;
    }
  }
  
  //default all gains to 1 if they were not loaded from the DB and/or if they were loaded with the 
  //wrong size: 
  if( fUgain.size() != fNAPVs_U ){
    fUgain.resize(fNAPVs_U);
    for( int iAPV=0; iAPV<fNAPVs_U; iAPV++ ){
      fUgain[iAPV] = 1.0;
    }
  }

  if( fVgain.size() != fNAPVs_V ){
    fVgain.resize(fNAPVs_V);
    for( int iAPV=0; iAPV<fNAPVs_V; iAPV++ ){
      fVgain[iAPV] = 1.0;
    }
  }

  if( fPedestalMode ){
    fZeroSuppress = false;
    fOnlineZeroSuppression = false;
    //fPedSubFlag = 0;
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
      // The detector axes are defined as the columns of the rotation matrix. We will have to test that it is working correctly, however:
      
      RotTemp.RotateX( angles[0] * TMath::DegToRad() );
      RotTemp.RotateY( angles[1] * TMath::DegToRad() );
      RotTemp.RotateZ( angles[2] * TMath::DegToRad() );
      
      fXax.SetXYZ( RotTemp.XX(), RotTemp.YX(), RotTemp.ZX() );
      fYax.SetXYZ( RotTemp.XY(), RotTemp.YY(), RotTemp.ZY() );
      fZax.SetXYZ( RotTemp.XZ(), RotTemp.YZ(), RotTemp.ZZ() );
    }
  } else
    DefineAxes(0);

  return 0;
}

    
Int_t SBSGEMModule::DefineVariables( EMode mode ) {
  if( mode == kDefine and fIsSetup ) return kOK;
  fIsSetup = ( mode == kDefine );

  //Int_t nstripsmax = fNstripsU + fNstripsV;
  
  
  VarDef varstrip[] = {
    {"strip.nstripsfired", "Number of strips fired", kUInt, 0, &fNstrips_hit },
    {"strip.istrip", "strip index", kUInt, 0, &(fStrip[0]), &fNstrips_hit },
    {"strip.IsU", "U strip?", kUInt, 0, &(fStripIsU[0]), &fNstrips_hit },
    {"strip.IsV", "V strip?", kUInt, 0, &(fStripIsV[0]), &fNstrips_hit },
    {"strip.ADCsamples", "ADC samples (index = isamp+Nsamples*istrip)", kDouble, 0, &(fADCsamples1D[0]), &fNdecoded_ADCsamples },
    { "strip.rawADCsamples", "raw ADC samples (no baseline subtraction)", kInt, 0, &(fRawADCsamples1D[0]), &fNdecoded_ADCsamples },
    { "strip.ADCsum", "Sum of ADC samples on a strip", kDouble, 0, &(fADCsums[0]), &fNstrips_hit },
    { "strip.isampmax", "sample in which max ADC occurred on a strip", kUInt, 0, &(fMaxSamp[0]), &fNstrips_hit },
    { "strip.ADCmax", "Value of max ADC sample on a strip", kDouble, 0, &(fADCmax[0]), &fNstrips_hit },
    { "strip.Tmean", "ADC-weighted mean strip time", kDouble, 0, &(fTmean[0]), &fNstrips_hit },
    { "strip.Tsigma", "ADC-weighted rms strip time", kDouble, 0, &(fTsigma[0]), &fNstrips_hit },
    { "strip.Tcorr", "Corrected strip time", kDouble, 0, &(fTcorr[0]), &fNstrips_hit },
    { "strip.itrack", "Index of track containing this strip (-1 if not on any track)", kInt, 0, &(fStripTrackIndex[0]), &fNstrips_hit },
    { "strip.ADCavg", "average of ADC samples on a strip", kDouble, 0, &(fStripADCavg[0]), &fNstrips_hit },
    { nullptr },
  };
  
  // //Raw strip info:
  // RVarDef varstrip[] = {
  //   { "nstripsfired",   "Number of strips fired",   "fNstrips_hit" },
  //   { "strip", "Strip index", "fStrip" },
  //   { "stripIsU", "U strip?", "fStripIsU"},
  //   { "stripIsV", "V strip?", "fStripIsV"},
  //   { "stripADCsamples", "ADC samples (index = isamp+Nsamples*istrip)", "fADCsamples1D" },
  //   { "striprawADCsamples", "raw ADC samples (no baseline subtraction)", "fRawADCsamples1D" },
  //   { "stripADCsum", "Sum of ADC samples on a strip", "fADCsums" },
  //   { "stripisampmax", "sample in which max ADC occurred on a strip", "fMaxSamp" },
  //   { "stripADCmax", "Value of max ADC sample on a strip", "fADCmax" },
  //   { "stripTmean", "ADC-weighted mean strip time", "fTmean" },
  //   { "stripTsigma", "ADC-weighted rms strip time", "fTsigma" },
  //   { "stripTcorr", "Corrected strip time", "fTcorr" },
  //   { "stripItrack", "Index of track containing this strip (-1 if not on any track)", "fStripTrackIndex" },
  //   { "stripADCavg", "average of ADC samples on a strip", "fStripADCavg" },
  //   { nullptr },
  // };


  Int_t ret = DefineVarsFromList( varstrip, mode );

  if( ret != kOK )
    return ret;
  

  RVarDef varclust[] = {
    { "clust.nclustu",   "Number of clusters in u",   "fNclustU" },
    { "clust.clustu_strips",   "u clusters strip multiplicity",   "fUclusters.nstrips" },
    { "clust.clustu_pos",   "u clusters position",   "fUclusters.hitpos_mean" },
    { "clust.clustu_adc",   "u clusters adc sum",   "fUclusters.clusterADCsum" },
    { "clust.clustu_time",   "u clusters time",   "fUclusters.t_mean" },
    { "clust.nclustv",   "Number of clusters in v",   "fNclustV" },
    { "clust.clustv_strips",   "v clusters strip multiplicity",   "fVclusters.nstrips" },
    { "clust.clustv_pos",   "v clusters position",   "fVclusters.hitpos_mean" },
    { "clust.clustv_adc",   "v clusters adc sum",   "fVclusters.clusterADCsum" },
    { "clust.clustv_time",   "v clusters time",   "fVclusters.t_mean" },
    { nullptr },
  };

  ret = DefineVarsFromList( varclust, mode );

  if( ret != kOK )
    return ret;

  RVarDef varhits[] = {
    { "hit.nhits2d",   "Number of 2d hits",   "fN2Dhits" },
    { "hit.hitx",   "local X coordinate of hit",   "fHits.xhit" },
    { "hit.hity",   "local Y coordinate of hit",   "fHits.yhit" },
    { "hit.hitxg",   "transport X coordinate of hit",   "fHits.xghit" },
    { "hit.hityg",   "transport Y coordinate of hit",   "fHits.yghit" },
    { "hit.hitADCasym",   "hit ADC asymmetry (ADCU-ADCV)/2",   "fHits.ADCasym" },
    { "hit.hitADCavg",  "(ADCU+ADCV)/2", "fHits.Ehit" },
    { "hit.hitTdiff",   "hit time difference (u-v)",   "fHits.tdiff" },
    { "hit.hitTavg",   "average time of 2D hit", "fHits.thitcorr" },
    { "hit.hit_iuclust", "index in u cluster array", "fHits.iuclust" },
    { "hit.hit_ivclust", "index in v cluster array", "fHits.ivclust" },
    { nullptr },
  };

  ret = DefineVarsFromList( varhits, mode );

  RVarDef vartiming[] = {
    { "time.T0_by_APV", "Coarse MPD timestamp of first event", "fT0_by_APV" },
    { "time.Tref_coarse", "Reference coarse MPD time stamp for this event", "fTref_coarse" },
    { "time.Tcoarse_by_APV", "Coarse MPD timestamp by APV relative to Tref_coarse", "fTcoarse_by_APV" },
    { "time.Tfine_by_APV", "Fine MPD timestamp by APV", "fTfine_by_APV" },
    { "time.EventCount_by_APV", "MPD event counter by APV (these should all agree in any one event)", "fEventCount_by_APV" },
    { "time.T_ns_by_APV", "Time stamp in ns relative to coarse T_ref", "fTimeStamp_ns_by_APV" },
    { nullptr },
  };

  ret = DefineVarsFromList( vartiming, mode );
  
  if( ret != kOK )
    return ret;

  return kOK;
    
}

void SBSGEMModule::Clear( Option_t* opt){ //we will want to clear out many more things too
  // Modify this a little bit so we only clear out the "hit counters", not necessarily the
  // arrays themselves, to make the decoding more efficient:
  fNstrips_hit = 0;
  fNstrips_hitU = 0;
  fNstrips_hitV = 0;
  fNdecoded_ADCsamples = 0;
  fIsDecoded = false;

  fNclustU = 0;
  fNclustV = 0;
  //later we may need to check whether this is a performance bottleneck:
  fUclusters.clear();
  fVclusters.clear();
  fN2Dhits = 0;
  //similar here:
  fHits.clear();

  //fStripAxis.clear();
  // fADCsamples1D.clear();
  // fStripTrackIndex.clear();
  // fRawADCsamples1D.clear();

  // fStripIsU.clear();
  // fStripIsV.clear();
  // fStripADCavg.clear();
  
  // fUstripIndex.clear();
  // fVstripIndex.clear();
  // fStrip.clear();
  // fAxis.clear();
  // fADCsamples.clear();
  // fRawADCsamples.clear();
  // fADCsums.clear();
  // fKeepStrip.clear();
  // fMaxSamp.clear();
  // fADCmax.clear();
  // fTmean.clear();
  // fTsigma.clear();
  // fTcorr.clear();
  
  return;
}

Int_t   SBSGEMModule::Decode( const THaEvData& evdata ){
  //std::cout << "[SBSGEMModule::Decode " << fName << "]" << std::endl;

  //initialize generic "strip" counter to zero:
  fNstrips_hit = 0;
  //initialize "U" and "V" strip counters to zero:
  fNstrips_hitU = 0;
  fNstrips_hitV = 0;

  //UInt_t MAXNSAMP_PER_APV = fN_APV25_CHAN * fN_MPD_TIME_SAMP;

  //std::cout << "MAXNSAMP_PER_APV = " << MAXNSAMP_PER_APV << std::endl;

  //we could save some time on these allocations by making these data members of SBSGEMModule: these are probably expensive:

  //to avoid rewriting the other code below, declare references to the fixed-size arrays:
  vector<UInt_t> &Strip = fStripAPV;
  vector<UInt_t> &rawStrip = fRawStripAPV;
  vector<Int_t> &rawADC = fRawADC_APV;
  vector<Double_t> &pedsubADC = fPedSubADC_APV; //ped-subtracted, not necessarily common-mode subtracted
  vector<Double_t> &commonModeSubtractedADC = fCommonModeSubtractedADC_APV;
  
  //resize all the "decoded strip" arrays to their maximum possible values for this module:
  //we need to do this event-by-event, because we shrink the size of the arrays to fNstrips_hit after decoding to prevent
  //enormous ROOT output:
  //UInt_t nstripsmax = fNstripsU + fNstripsV;
  
  // fStrip.resize( nstripsmax );
  // fAxis.resize( nstripsmax );
  // fADCsamples.resize( nstripsmax );
  // fRawADCsamples.resize( nstripsmax );
  // for( int istrip=0; istrip<nstripsmax; istrip++ ){
  //   fADCsamples[istrip].resize( fN_MPD_TIME_SAMP );
  //   fRawADCsamples[istrip].resize( fN_MPD_TIME_SAMP );
  // }
  // fADCsums.resize( nstripsmax );
  // fStripADCavg.resize( nstripsmax );
  // fStripIsU.resize( nstripsmax );
  // fStripIsV.resize( nstripsmax );
  // fKeepStrip.resize( nstripsmax );
  // fMaxSamp.resize( nstripsmax );
  // fADCmax.resize( nstripsmax );
  // fTmean.resize( nstripsmax );
  // fTsigma.resize( nstripsmax );
  // fTcorr.resize( nstripsmax );
  
  // fADCsamples1D.resize( nstripsmax * fN_MPD_TIME_SAMP );
  // fRawADCsamples1D.resize( nstripsmax * fN_MPD_TIME_SAMP );
  // fStripTrackIndex.resize( nstripsmax );
  
  // fUstripIndex.clear();
  // fVstripIndex.clear();
  //This could be written more efficiently, in principle. However, it's not yet clear it's a speed bottleneck, so for now let's not worry about it too much:

  //Do we need to loop on all APV cards? maybe not,
  int apvcounter=0;
  
  for (std::vector<mpdmap_t>::iterator it = fMPDmap.begin() ; it != fMPDmap.end(); ++it){
    //loop over all decode map entries associated with this module (each decode map entry is one APV card)
    Int_t effChan = it->mpd_id << 4 | it->adc_id; //left-shift mpd id by 4 bits and take the bitwise OR with ADC_id to uniquely identify the APV card.
    //mpd_id is not necessarily equal to slot, but that seems to be the convention in many cases
    // Find channel for this crate/slot

    //First get time stamp info:
    UInt_t nhits_timestamp_low = evdata.GetNumHits( it->crate, it->slot, fChan_TimeStamp_low );
    UInt_t nhits_timestamp_high = evdata.GetNumHits( it->crate, it->slot, fChan_TimeStamp_high );
    UInt_t nhits_event_count = evdata.GetNumHits( it->crate, it->slot, fChan_MPD_EventCount );

    if( nhits_timestamp_low > 0 && nhits_timestamp_high == nhits_timestamp_low && nhits_event_count == nhits_timestamp_low ){
      for( int ihit=0; ihit<nhits_timestamp_low; ihit++ ){
	int fiber = evdata.GetRawData( it->crate, it->slot, fChan_TimeStamp_low, ihit );
	if( fiber == it->mpd_id ){ //this is the channel we want: 
	  UInt_t Tlow = evdata.GetData( it->crate, it->slot, fChan_TimeStamp_low, ihit );
	  UInt_t Thigh = evdata.GetData( it->crate, it->slot, fChan_TimeStamp_high, ihit );
	  UInt_t EvCnt = evdata.GetData( it->crate, it->slot, fChan_MPD_EventCount, ihit );

	  fEventCount_by_APV[apvcounter] = EvCnt;
	  
	  // Fine time stamp is in the first 8 bits of Tlow;
	  fTfine_by_APV[apvcounter] = Tlow & 0x000000FF;

	  Long64_t Tcoarse = Thigh << 16 | ( Tlow << 8 ); 
	  if( EvCnt == 0 ) fT0_by_APV[apvcounter] = Tcoarse;

	  //T ref is the coarse time stamp of the reference APV (the first one, in this case)
	  if( apvcounter == 0 ) fTref_coarse = Tcoarse - fT0_by_APV[apvcounter];

	  //This SHOULD make fTcoarse_by_APV the Tcoarse RELATIVE to the
	  // "reference" APV
	  fTcoarse_by_APV[apvcounter] = Tcoarse - fT0_by_APV[apvcounter] - fTref_coarse;

	  //We probably don't want to hard-code 24 ns and 4 ns here for the units of
	  //Tcoarse and Tfine, but this should be fine for initial checkout of decoding:
	  fTimeStamp_ns_by_APV[apvcounter] = 24.0 * fTcoarse_by_APV[apvcounter] + 4.0 * (fTfine_by_APV[apvcounter] % 6);
	  
	  break;
	}	
      }
    }
    
    // Get common-mode flags, if applicable:
    // Default to the values from the database (or the default values):

    Bool_t CM_ENABLED = fCommonModeFlag != 0 && fCommonModeFlag != 1 && !fPedestalMode;
    Bool_t BUILD_ALL_SAMPLES = !fOnlineZeroSuppression;
    Bool_t CM_OUT_OF_RANGE = false;
    
    UInt_t cm_flags=4*CM_OUT_OF_RANGE + 2*CM_ENABLED + BUILD_ALL_SAMPLES;
    UInt_t nhits_cm_flag=evdata.GetNumHits( it->crate, it->slot, fChan_CM_flags );
    
    bool cm_flags_found = false;
											       
    if( nhits_cm_flag > 0 ){
      
      // If applicable, find the common-mode/zero-suppression settings loaded from the raw data for this APV:
      // In principle in the SSP/VTP event format, there should be exactly one "hit" per APV in this "channel":
      for( int ihit=0; ihit<nhits_cm_flag; ihit++ ){ 
	int chan_temp = evdata.GetRawData( it->crate, it->slot, fChan_CM_flags, ihit );
	if( chan_temp == effChan ){ //assume that this is only filled once per MPD per event, and exit the loop when we find this MPD:
	  // std::cout << "Before decoding cm flags, CM_ENABLED, BUILD_ALL_SAMPLES = " << CM_ENABLED << ", "
	  // 	    << BUILD_ALL_SAMPLES << std::endl;
	  cm_flags = evdata.GetData( it->crate, it->slot, fChan_CM_flags, ihit );
	  cm_flags_found = true;
	  break;
	}
      }
    }

    //The proper logic of common-mode calculation/subtraction and zero suppression is as follows:
    // 1. If CM_ENABLED is true, we never calculate the common-mode ourselves, it has already been subtracted from the data:
    // 2. If BUILD_ALL_SAMPLES is false, then online zero suppression is enabled. We can, in addition, apply our own higher thresholds if we want:
    // 3. If CM_ENABLED is true, the pedestal has also been subtracted, so we don't subtract it again.
    // 4. If CM_ENABLED is false, we need to subtract the pedestals AND calculate and subtract the common-mode:
    // 5. If BUILD_ALL_SAMPLES is false then CM_ENABLED had better be true!
    // 6. If CM_OUT_OF_RANGE is true then BUILD_ALL_SAMPLES must be true and CM_ENABLED
    //    must be false!
    
    CM_OUT_OF_RANGE = cm_flags/4;
    CM_ENABLED = cm_flags/2;
    BUILD_ALL_SAMPLES = cm_flags%2;
    
    // if( cm_flags_found ){
    //   std::cout << "cm flag defaults overridden by raw data, effChan = " << effChan << ", CM_ENABLED = " << CM_ENABLED << ", BUILD_ALL_SAMPLES = " << BUILD_ALL_SAMPLES << std::endl;
    // }
    
    //fOnlineZeroSuppression = !BUILD_ALL_SAMPLES;
    if( !BUILD_ALL_SAMPLES && !CM_ENABLED ) { //This should never happen: skip this APV card
      continue;
    }
    if( CM_OUT_OF_RANGE && !BUILD_ALL_SAMPLES ){ // this should also never happen: skip this APV card:
      continue;
    }

    if( CM_OUT_OF_RANGE && CM_ENABLED ){ //force CM_ENABLED to false if CM_OUT_OF_RANGE is true;
      //This should have already been done online:
      CM_ENABLED = false; 
    }
    
    //Int_t nchan = evdata.GetNumChan( it->crate, it->slot ); //this could be made faster

    SBSGEM::GEMaxis_t axis = it->axis == 0 ? SBSGEM::kUaxis : SBSGEM::kVaxis; 
    
    //printf("nchan = %d\n", nchan );

    //std::cout << "crate, slot, nchan = " << it->crate << ", " << it->slot << ", " << nchan << std::endl;

    //this is looping on all the 
    //for( Int_t ichan = 0; ichan < nchan; ++ichan ) { //this is looping over all the "channels" (APV cards) in the crate and slot containing this decode map entry/APV card:
    //Int_t chan = evdata.GetNextChan( it->crate, it->slot, ichan ); //"chan" here refers to one APV card 
      //std::cout << it->crate << " " << it->slot << " mpd_id ??? " << it->mpd_id << " " << chan << " " << effChan << std::endl;

      //if( chan != effChan ) continue; // 

    Int_t nsamp = evdata.GetNumHits( it->crate, it->slot, effChan );

    
    if( nsamp > 0 ){
      
      assert(nsamp%fN_MPD_TIME_SAMP==0); //this is making sure that the number of samples is equal to an integer multiple of the number of time samples per strip
      Int_t nstrips = nsamp/fN_MPD_TIME_SAMP; //number of strips fired on this APV card (should be exactly 128 if online zero suppression is NOT used):
      
      // std::cout << "MPD ID, ADC channel, number of strips fired = " << it->mpd_id << ", "
      // 		<< it->adc_id << ", " << nstrips << std::endl;

      double commonMode[fN_MPD_TIME_SAMP];
	
      for( int isamp=0; isamp<fN_MPD_TIME_SAMP; isamp++ ){
	commonMode[isamp] = 0.0;
      }
      
      if( !CM_ENABLED && BUILD_ALL_SAMPLES && nstrips == fN_APV25_CHAN ){ //then two loops over the data are necessary, first one to calculate common-mode:
	//declare temporary array to hold common mode values for this APV card and, if necessary, calculate them:

	//std::cout << "Common-mode calculation: " << std::endl;
	
	//First loop over the hits: populate strip, raw strip, raw ADC, ped sub ADC and common-mode-subtracted aDC:
	for( int iraw=0; iraw<nsamp; iraw++ ){ //NOTE: iraw = isamp + fN_MPD_TIME_SAMP * istrip
	  int strip = evdata.GetRawData( it->crate, it->slot, effChan, iraw );
	  UInt_t decoded_rawADC = evdata.GetData( it->crate, it->slot, effChan, iraw );

	  Int_t ADC = Int_t( decoded_rawADC );
	  
	  rawStrip[iraw] = strip;
	  Strip[iraw] = GetStripNumber( strip, it->pos, it->invert );

	  double ped = (axis == SBSGEM::kUaxis ) ? fPedestalU[Strip[iraw]] : fPedestalV[Strip[iraw]];

	  // If pedestal subtraction was done online, don't do it again:
	  if( fPedSubFlag != 0 ) ped = 0.0;
	
	  rawADC[iraw] = ADC;
	  pedsubADC[iraw] = double(ADC) - ped;
	  commonModeSubtractedADC[iraw] = double(ADC) - ped; 

	  //the calculation of common mode in pedestal mode analysis differs from the
	  // offline or online zero suppression analysis; here we use a simple average of all 128 channels:
	  if( fPedestalMode ){
	    //do simple common-mode calculation involving the simple average of all 128 (ped-subtracted) ADC
	    //values
	    int isamp = iraw%fN_MPD_TIME_SAMP;
	    commonMode[isamp] += pedsubADC[iraw]/double(fN_APV25_CHAN);
	  }
	}
	
	// second loop over the hits to calculate and apply common-mode correction (sorting method)
	//if( !fPedestalMode ){ //need to calculate common mode:
	if( fMakeCommonModePlots || !fPedestalMode ){ // calculate both ways:
	  for( int isamp=0; isamp<fN_MPD_TIME_SAMP; isamp++ ){
	    
	    //moved common-mode calculation to its own function:
	    
	    // Now: a question is should we modify the behavior/do added checks if
	    // CM_OUT_OF_RANGE is set? 
	  
	    if( fMakeCommonModePlots ){
	      double cm_danning = GetCommonMode( isamp, 1, *it );
	      double cm_sorting = GetCommonMode( isamp, 0, *it );

	      //std::cout << "cm danning, sorting = " << cm_danning << ", " << cm_sorting << std::endl;
	      
	      if( !fPedestalMode ) commonMode[isamp] = fCommonModeFlag == 0 ? cm_sorting : cm_danning;

	      double cm_mean;
	      
	      UInt_t iAPV = it->pos;
	      
	      // std::cout << "Filling common-mode histograms..." << std::endl;
	      
	      // std::cout << "iAPV, nAPVsU, nAPVsV, axis = " << iAPV << ", " << fNAPVs_U << ", "
	      // 		<< fNAPVs_V << ", " << axis << std::endl;
	      
	      if( axis == SBSGEM::kUaxis ){
		cm_mean = fCommonModeMeanU[iAPV];
		
		fCommonModeDistU->Fill( iAPV, commonMode[isamp] - cm_mean );
		fCommonModeDistU_Sorting->Fill( iAPV, cm_sorting - cm_mean );
		fCommonModeDistU_Danning->Fill( iAPV, cm_danning - cm_mean );
		fCommonModeDiffU->Fill( iAPV, cm_sorting - cm_danning );
	      } else {
		cm_mean = fCommonModeMeanV[iAPV];
		
		fCommonModeDistV->Fill( iAPV, commonMode[isamp] - cm_mean );
		fCommonModeDistV_Sorting->Fill( iAPV, cm_sorting - cm_mean );
		fCommonModeDistV_Danning->Fill( iAPV, cm_danning - cm_mean );
		fCommonModeDiffV->Fill( iAPV, cm_sorting - cm_danning );
	      }

	      //std::cout << "Done..." << std::endl;
	      
	    } else if( !fPedestalMode ) { //if not doing diagnostic plots, just calculate whichever way the user wanted:
	      
	      commonMode[isamp] = GetCommonMode( isamp, fCommonModeFlag, *it );
	      
	    }
	    //std::cout << "effChan, isamp, Common-mode = " << effChan << ", " << isamp << ", " << commonMode[isamp] << std::endl;
	    
	  } //loop over time samples
	} //check if conditions are satisfied to require offline common-mode calculation
      } //End check !CM_ENABLED && BUILD_ALL_SAMPLES

      
      //std::cout << "finished common mode " << std::endl;
      // Last loop over all the strips and samples in the data and populate/calculate global variables that are passed to track-finding:
      //Int_t ihit = 0;
            
      for( Int_t istrip = 0; istrip < nstrips; ++istrip ) {
	if( CM_ENABLED ){
	  //then we skipped the first loop over the data; need to grab the actual data:
	  for( int isamp=0; isamp<fN_MPD_TIME_SAMP; isamp++ ){
	    int iraw = isamp + fN_MPD_TIME_SAMP * istrip;
	    int strip = evdata.GetRawData( it->crate, it->slot, effChan, iraw );
	    int ADC = evdata.GetData( it->crate, it->slot, effChan, iraw );
	    rawStrip[iraw] = strip;
	    Strip[iraw] = GetStripNumber( strip, it->pos, it->invert );
	    //no need to grab pedestal if CM_ENABLED is true:
	    
	    double ped = 0;
	    if(fIsMC)ped = (axis == SBSGEM::kUaxis ) ? fPedestalU[Strip[iraw]] : fPedestalV[Strip[iraw]];
	    
	    rawADC[iraw] = Int_t(ADC);
	    pedsubADC[iraw] = double(ADC) - ped;
	    commonModeSubtractedADC[iraw] = double(ADC) - ped;
	  }
	}
	
	//Temporary vector to hold ped-subtracted ADC samples for this strip:
	std::vector<double> ADCtemp(fN_MPD_TIME_SAMP);
	std::vector<int> rawADCtemp(fN_MPD_TIME_SAMP);
	
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
	
	//"pedtemp" is only used to fill pedestal histograms as of now:
	double pedtemp = ( axis == SBSGEM::kUaxis ) ? fPedestalU[strip] : fPedestalV[strip];

	if( fPedSubFlag != 0 && !fIsMC ) pedtemp = 0.0;

	double rmstemp = ( axis == SBSGEM::kUaxis ) ? fPedRMSU[strip] : fPedRMSV[strip];
	double gaintemp = ( axis == SBSGEM::kUaxis ) ? fUgain[strip/fN_APV25_CHAN] : fVgain[strip/fN_APV25_CHAN]; //should probably not hard-code 128 here
	
	// std::cout << "pedestal temp, rms temp, nsigma cut, threshold, zero suppress, pedestal mode = " << pedtemp << ", " << rmstemp << ", " << fZeroSuppressRMS
	//   	  << ", " << fZeroSuppressRMS * rmstemp << ", " << fZeroSuppress << ", " << fPedestalMode << std::endl;
	
	//Now loop over the time samples:
	for( Int_t adc_samp = 0; adc_samp < fN_MPD_TIME_SAMP; adc_samp++ ){

	  int iraw = adc_samp + fN_MPD_TIME_SAMP * istrip;

	  //If applicable, subtract common-mode here:
	  //if( (fPedestalMode || !fOnlineZeroSuppression) && nstrips == fN_APV25_CHAN ){
	  //We need to subtract the common-mode if it was calculated:
	  if( !CM_ENABLED && BUILD_ALL_SAMPLES && nstrips == fN_APV25_CHAN ){
	    
	    // std::cout << "isamp, commonMode = " << adc_samp << ", " << commonMode[adc_samp]
	    // 	      << std::endl;
	    commonModeSubtractedADC[ iraw ] = pedsubADC[ iraw ] - commonMode[adc_samp];
	  }
	  
	  // Int_t ihit = adc_samp + fN_MPD_TIME_SAMP * istrip; //index in the "hit" array for this APV card:
	  // assert(ihit<nsamp);
	  
	  
	  //Int_t rawADC = evdata.GetData(it->crate, it->slot, chan, ihit);
	  Int_t RawADC = rawADC[iraw]; //this value has no corrections applied:
         
	  //cout << adc_samp << " " << istrip << " " << rawADC << " ";// << endl;
	  
	  //rawADCtemp.push_back( RawADC );
	  rawADCtemp[adc_samp] = RawADC;
	  
	  //The following value already has pedestal and common-mode subtracted (if applicable):
	  double ADCvalue = commonModeSubtractedADC[iraw]; //zero-suppress BEFORE we apply gain correction

	  // if( fPedestalMode ){ //If we are analyzing pedestal data, DON'T substract the pedestal
	  //   ADCvalue = commonModeSubtractedADC[adc_samp][istrip];
	  // }
	  //pedestal-subtracted ADC values:
	  //ADCtemp.push_back( ADCvalue );
	  ADCtemp[adc_samp] = ADCvalue;
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

	// if( !CM_ENABLED && BUILD_ALL_SAMPLES ){
	//   std::cout << " before zero suppression: common-mode and pedestal-subtracted ADC sum: effChan, istrip, ADC sum = " 
	// 	    << effChan << ", " << istrip << ", " << ADCsum_temp << std::endl;
	// }

	//cout << ADCsum_temp << " >? " << fThresholdStripSum << "; " 
	//   << ADCsum_temp/double(fN_MPD_TIME_SAMP) << " >? " << fZeroSuppressRMS*rmstemp << "; " 
	//   << maxADC << " >? " << fZeroSuppressRMS*rmstemp 
	//   << "; " << fThresholdSample << endl;

	//fill pedestal diagnostic histograms if and only if we are in pedestal mode or plot common mode 
	// AND the CM_ENABLED is not set, meaning we did cm and ped subtraction offline
	if( (fPedestalMode || fMakeCommonModePlots) && !CM_ENABLED ){ 
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
	      // ( (TH2D*) (*hrawADCs_by_strip_sampleU)[isamp] )->Fill( strip, rawADCtemp[isamp] );
	      // //for this one, we add back in the pedestal:
	      // ( (TH2D*) (*hcommonmode_subtracted_ADCs_by_strip_sampleU)[isamp] )->Fill( strip, ADCtemp[isamp] + pedtemp );
	      // ( (TH2D*) (*hpedestal_subtracted_ADCs_by_strip_sampleU)[isamp] )->Fill( strip, ADCtemp[isamp] );
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
		
	      // ( (TH2D*) (*hrawADCs_by_strip_sampleV)[isamp] )->Fill( strip, rawADCtemp[isamp] );
	      // ( (TH2D*) (*hcommonmode_subtracted_ADCs_by_strip_sampleV)[isamp] )->Fill( strip, ADCtemp[isamp] );
	      // ( (TH2D*) (*hpedestal_subtracted_ADCs_by_strip_sampleV)[isamp] )->Fill( strip, ADCtemp[isamp] - pedtemp );
	    }
	  }

	  // std::cout << "finished pedestal histograms..." << std::endl;
	    
	}
	
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
	
	  
	  //fStrip.push_back( strip );
	  //fAxis.push_back( axis );

	  //std::cout << "strip, axis = " << strip << ", " << axis << std::endl;
	  
	  fStrip[fNstrips_hit] = strip;
	  fAxis[fNstrips_hit] = axis;
	  
	  //std::cout << "axis, Int_t(axis) = " << axis << ", " << Int_t(axis) << std::endl;
	  //fStripAxis.push_back( Int_t(axis) );
	  // fADCsamples.push_back( ADCtemp ); //pedestal-subtracted
	  // fRawADCsamples.push_back( rawADCtemp ); //Raw
	  for( Int_t isamp=0; isamp<fN_MPD_TIME_SAMP; isamp++ ){
	    fADCsamples[fNstrips_hit][isamp] = ADCtemp[isamp];
	    fRawADCsamples[fNstrips_hit][isamp] = rawADCtemp[isamp];
	    
	    //fADCsamples1D.push_back( ADCtemp[isamp] );
	    //fRawADCsamples1D.push_back( rawADCtemp[isamp] );
	    fADCsamples1D[isamp + fN_MPD_TIME_SAMP * fNstrips_hit ] = ADCtemp[isamp];
	    fRawADCsamples1D[isamp + fN_MPD_TIME_SAMP * fNstrips_hit ] = rawADCtemp[isamp];
	  }
	  //fStripTrackIndex.push_back( -1 ); //This could be modified later based on tracking results
	  fStripTrackIndex[fNstrips_hit] = -1;
	  
	  //	  fKeepStrip.push_back( true ); //keep all strips by default
	  fKeepStrip[fNstrips_hit] = true;
	  //	  fMaxSamp.push_back( iSampMax );
	  fMaxSamp[fNstrips_hit] = iSampMax;
	  
	  if( !fPedestalMode ) maxADC /= gaintemp;
	  
	  //	  fADCmax.push_back( maxADC );
	  fADCmax[fNstrips_hit] = maxADC; 
	  //	  fTmean.push_back( Tsum/ADCsum_temp );
	  fTmean[fNstrips_hit] = Tsum/ADCsum_temp; 
	  //  fTsigma.push_back( sqrt( T2sum/ADCsum_temp - pow( fTmean.back(), 2 ) ) );
	  fTsigma[fNstrips_hit] = sqrt( T2sum/ADCsum_temp - pow( fTmean[fNstrips_hit], 2 ) );
	  //fTcorr.push_back( fTmean.back() ); //don't apply any corrections for now
	  fTcorr[fNstrips_hit] = fTmean[fNstrips_hit];
	  
	  if( !fPedestalMode ) ADCsum_temp /= gaintemp;

	  //  fADCsums.push_back( ADCsum_temp ); //sum of all (pedestal-subtracted) samples
	  fADCsums[fNstrips_hit] = ADCsum_temp;
	  
	  //  fStripADCavg.push_back( ADCsum_temp/double(fN_MPD_TIME_SAMP) );
	  fStripADCavg[fNstrips_hit] = ADCsum_temp/double(fN_MPD_TIME_SAMP);
	  
	  UInt_t isU = (axis == SBSGEM::kUaxis) ? 1 : 0;
	  UInt_t isV = (axis == SBSGEM::kVaxis) ? 1 : 0;
	  //	  fStripIsU.push_back( isU );
	  //      fStripIsV.push_back( isV );
	  fStripIsU[fNstrips_hit] = isU;
	  fStripIsV[fNstrips_hit] = isV;

	  fNstrips_hitU += isU;
	  fNstrips_hitV += isV;
	  
	  //	  if( axis == SBSGEM::kUaxis ) fUstripIndex[strip] = fNstrips_hit;
	  //      if( axis == SBSGEM::kVaxis ) fVstripIndex[strip] = fNstrips_hit;

	  //std::cout << "starting pedestal histograms..." << std::endl;
	  
	  
	  
	  fNstrips_hit++;
	} //check if passed zero suppression cuts
      } //end loop over strips on this APV card
    } //end if( nsamp > 0 )

    apvcounter++;
  } //end loop on decode map entries for this module

  fNdecoded_ADCsamples = fNstrips_hit * fN_MPD_TIME_SAMP;
  
  //We will want to resize the 

  //resize all the "decoded strip" arrays to the actual number of fired strips:
  
  // fStrip.resize( fNstrips_hit );
  // fAxis.resize( fNstrips_hit );
  // fADCsamples.resize( fNstrips_hit );
  // fRawADCsamples.resize( fNstrips_hit );
  // for( int istrip=0; istrip<fNstrips_hit; istrip++ ){
  //   fADCsamples[istrip].resize( fN_MPD_TIME_SAMP );
  //   fRawADCsamples[istrip].resize( fN_MPD_TIME_SAMP );
  // }
  // fADCsums.resize( fNstrips_hit );
  // fStripADCavg.resize( fNstrips_hit );
  // fStripIsU.resize( fNstrips_hit );
  // fStripIsV.resize( fNstrips_hit );
  // fKeepStrip.resize( fNstrips_hit );
  // fMaxSamp.resize( fNstrips_hit );
  // fADCmax.resize( fNstrips_hit );
  // fTmean.resize( fNstrips_hit );
  // fTsigma.resize( fNstrips_hit );
  // fTcorr.resize( fNstrips_hit );

  // fADCsamples1D.resize( fNstrips_hit * fN_MPD_TIME_SAMP );
  // fRawADCsamples1D.resize( fNstrips_hit * fN_MPD_TIME_SAMP );
  // fStripTrackIndex.resize( fNstrips_hit );
  
  // No longer necessary:
  //fNstrips_hitU = fUstripIndex.size();
  //fNstrips_hitV = fVstripIndex.size();
  
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
  Double_t offset = (axis == SBSGEM::kUaxis) ? fUStripOffset : fVStripOffset;
  
  //hopefully this compiles and works correctly:
  std::vector<sbsgemcluster_t> &clusters = (axis == SBSGEM::kUaxis ) ? fUclusters : fVclusters;

  UInt_t &nclust = ( axis == SBSGEM::kUaxis ) ? fNclustU : fNclustV; 

  nclust = 0;

  clusters.clear();
  
  if( axis == SBSGEM::kUaxis ){
    maxsep = fMaxNeighborsU_totalcharge;
    maxsepcoord = fMaxNeighborsU_hitpos; 
    Nstrips = fNstripsU;
    pitch = fUStripPitch;
    // fNclustU = 0;
    // fUclusters.clear();
  } else { //V axis, no need to compare axis to kVaxis:
    maxsep = fMaxNeighborsV_totalcharge;
    maxsepcoord = fMaxNeighborsV_hitpos; 
    Nstrips = fNstripsV;
    pitch = fVStripPitch;
    // fNclustV = 0;
    // fVclusters.clear();
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
  //for speed/efficiency, resize the cluster array to the theoretical maximum
  //and use operator[] rather than push_back:
  clusters.resize( localmaxima.size() );
  
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

    double maxpos = (stripmax + 0.5 - 0.5*Nstrips) * pitch + offset;

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

	double hitpos = (istrip + 0.5 - 0.5*Nstrips) * pitch + offset; //local hit position along direction measured by these strips
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

      //initialize "keep" flag for all 1D clusters to true:
      clusttemp.keep = true;
      
      // In the standalone we don't apply an independent threshold on the cluster sum in the context of 1D cluster-finding:
      // if( sumADC >= fThresholdClusterSum ){
      // if( axis == SBSGEM::kVaxis ){
      // 	fVclusters.push_back( clusttemp );
      // 	fNclustV++;
      // } else {
      // 	fUclusters.push_back( clusttemp );
      // 	fNclustU++;
      // }

      //Hopefully this works correctly:
      clusters[nclust] = clusttemp;
      nclust++;
      

	// std::cout << "found cluster, (axis, istripmax, nstrips, ADCsum, hit pos (mm) )=(" << axis << ", " << clusttemp.istripmax << ", "
	// 	  << clusttemp.nstrips << ", " << clusttemp.clusterADCsum
	// 	  << ", " << clusttemp.hitpos_mean*1000.0 << ")" << std::endl;
	
      //}
    } //Check if peak is inside track search region constraint
  } //end loop on local maxima

  clusters.resize(nclust); //just to make sure no pathological behavior later on
  
  filter_1Dhits(SBSGEM::kUaxis);
  filter_1Dhits(SBSGEM::kVaxis);
}

void SBSGEMModule::fill_2D_hit_arrays(){
  //Clear out the 2D hit array to get rid of any leftover junk from prior events:
  //fHits.clear();
  fN2Dhits = 0;

  //fHits.resize( fNclustU * fNclustV );
  //fHits.clear();
  fHits.resize( std::min( fNclustU*fNclustV, fMAX2DHITS ) );

  // if( fNclustU * fNclustV > fMAX2DHITS ){
  //   std::cout << "Warning in SBSGEMModule::fill_2D_hit_arrays(): 
  // }
  
  //This routine is simple: just form all possible 2D hits from combining one "U" cluster with one "V" cluster. Here we assume that find_clusters_1D has already
  //been called, if that is NOT the case, then this routine will just do nothing:
  for( int iu=0; iu<fNclustU; iu++ ){
    for( int iv=0; iv<fNclustV; iv++ ){

      if( fUclusters[iu].keep && fVclusters[iv].keep ){
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

	double pos_maxstripu = ( fUclusters[iu].istripmax + 0.5 - 0.5*fNstripsU ) * fUStripPitch + fUStripOffset;
	double pos_maxstripv = ( fVclusters[iv].istripmax + 0.5 - 0.5*fNstripsV ) * fVStripPitch + fVStripOffset;

	//"Cluster moments" defined as differences between reconstructed hit position and center of strip with max. signal in the cluster:
	hittemp.umom = (hittemp.uhit - pos_maxstripu)/fUStripPitch;
	hittemp.vmom = (hittemp.vhit - pos_maxstripv)/fVStripPitch;
      
	TVector2 UVtemp(hittemp.uhit,hittemp.vhit);
	TVector2 XYtemp = UVtoXY( UVtemp );

	hittemp.xhit = XYtemp.X();
	hittemp.yhit = XYtemp.Y();

	//Check if candidate 2D hit is inside the constraint region before doing anything else:
	if( fxcmin <= hittemp.xhit && hittemp.xhit <= fxcmax &&
	    fycmin <= hittemp.yhit && hittemp.yhit <= fycmax &&
	    IsInActiveArea( hittemp.xhit, hittemp.yhit ) &&
	    fN2Dhits < fMAX2DHITS ){
    
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
	  //fHits.push_back( hittemp );
	  fHits[fN2Dhits++] = hittemp; //should be faster than push_back();
	  //fN2Dhits++;
	} //end check that 2D point is inside track search region
      } //end check that both U and V clusters passed filtering criteria:
    } //end loop over "V" clusters
  } //end loop over "U" clusters

  fHits.resize( fN2Dhits );
  
  filter_2Dhits();
  
}

void    SBSGEMModule::Print( Option_t* opt) const{
  return;
}

Int_t   SBSGEMModule::Begin( THaRunBase* r){ //Does nothing
  //Here we can create some histograms that will be written to the ROOT file:
  //This is a natural place to do the hit maps/efficiency maps:
  fZeroSuppress = fZeroSuppress && !fPedestalMode; 
  
  //pulled these lines out of the if-block below to avoid code duplication:
  TString histname,histtitle;
  TString appname = (static_cast<THaDetector *>(GetParent()) )->GetApparatus()->GetName();
  appname.ReplaceAll(".","_");
  appname += "_";
  TString detname = GetParent()->GetName();
  detname.Prepend(appname);
  detname.ReplaceAll(".","_");
  detname += "_";
  detname += GetName();
  
  
  if( fMakeEfficiencyPlots && !fEfficiencyInitialized ){
    fEfficiencyInitialized = true;

    int nbinsx1D = int( round( 1.02*GetXSize() /fBinSize_efficiency1D ) );
    int nbinsy1D = int( round( 1.02*GetYSize() /fBinSize_efficiency1D ) );
    int nbinsx2D = int( round( 1.02*GetXSize() /fBinSize_efficiency2D ) );
    int nbinsy2D = int( round( 1.02*GetYSize() /fBinSize_efficiency2D ) );
    
    fhdidhitx = new TH1D( histname.Format( "hdidhitx_%s", detname.Data() ), "local x coordinate of hits on good tracks (m)", nbinsx1D, -0.51*GetXSize(), 0.51*GetXSize() );
    fhdidhity = new TH1D( histname.Format( "hdidhity_%s", detname.Data() ), "local y coordinate of hits on good tracks (m)", nbinsy1D, -0.51*GetYSize(), 0.51*GetYSize() );
    fhdidhitxy = new TH2D( histname.Format( "hdidhitxy_%s", detname.Data() ), "x vs y of hits on good tracks (m)",
			   nbinsy2D, -0.51*GetYSize(), 0.51*GetYSize(),
			   nbinsx2D, -0.51*GetXSize(), 0.51*GetXSize() );

    fhshouldhitx = new TH1D( histname.Format( "hshouldhitx_%s", detname.Data() ), "x of good track passing through (m)", nbinsx1D, -0.51*GetXSize(), 0.51*GetXSize() );
    fhshouldhity = new TH1D( histname.Format( "hshouldhity_%s", detname.Data() ), "y of good track passing through (m)", nbinsy1D, -0.51*GetYSize(), 0.51*GetYSize() );
    fhshouldhitxy = new TH2D( histname.Format( "hshouldhitxy_%s", detname.Data() ), "x vs y of good track passing through (m)",
			      nbinsy2D, -0.51*GetYSize(), 0.51*GetYSize(),
			      nbinsx2D, -0.51*GetXSize(), 0.51*GetXSize() );

    fEfficiencyInitialized = true;
  }

  if( (fPedestalMode || fMakeCommonModePlots) && !fPedHistosInitialized ){ //make pedestal histograms:

    //Procedure:
    // 1. Analyze pedestal data with no pedestals supplied from the database.
    // 2. Extract "coarse" strip offsets from raw ADC histograms (i.e., common-mode means)
    // 3. Extract "fine" strip offsets from common-mode-subtracted histograms 
    // 3. Add "coarse" strip offsets (common-mode means) back into "common-mode subtracted" ADC 

    //U strips:
    hpedrmsU_distribution = new TH1D( histname.Format( "hpedrmsU_distribution_%s", detname.Data() ), "Pedestal RMS distribution, U strips", 500, 0, 100 );
    hpedmeanU_distribution = new TH1D( histname.Format( "hpedmeanU_distribution_%s", detname.Data() ), "Pedestal mean distribution, U strips", 500, -500, 500 );
    hpedrmsU_by_strip = new TH1D( histname.Format( "hpedrmsU_by_strip_%s", detname.Data() ), "Pedestal rms U by strip", fNstripsU, -0.5, fNstripsU - 0.5 );
    hpedmeanU_by_strip = new TH1D( histname.Format( "hpedmeanU_by_strip_%s", detname.Data() ), "Pedestal mean U by strip", fNstripsU, -0.5, fNstripsU - 0.5 );

    //V strips:
    hpedrmsV_distribution = new TH1D( histname.Format( "hpedrmsV_distribution_%s", detname.Data() ), "Pedestal RMS distribution, V strips", 500, 0, 100 );
    hpedmeanV_distribution = new TH1D( histname.Format( "hpedmeanV_distribution_%s", detname.Data() ), "Pedestal mean distribution, V strips", 500, -500, 500 );
    hpedrmsV_by_strip = new TH1D( histname.Format( "hpedrmsV_by_strip_%s", detname.Data() ), "Pedestal rms V by strip", fNstripsV, -0.5, fNstripsV - 0.5 );
    hpedmeanV_by_strip = new TH1D( histname.Format( "hpedmeanV_by_strip_%s", detname.Data() ), "Pedestal mean V by strip", fNstripsV, -0.5, fNstripsV - 0.5 );

    if( !fPedestalMode ){ //fill the above histograms with the pedestal info loaded from the database:
      for( int istrip = 0; istrip<fNstripsU; istrip++ ){
	hpedmeanU_distribution->Fill( fPedestalU[istrip] );
	hpedrmsU_distribution->Fill( fPedRMSU[istrip] );
	hpedmeanU_by_strip->SetBinContent( istrip+1, fPedestalU[istrip] );
	hpedrmsU_by_strip->SetBinContent( istrip+1, fPedRMSU[istrip] );	
      }

      for( int istrip=0; istrip<fNstripsV; istrip++ ){
	hpedmeanV_distribution->Fill( fPedestalV[istrip] );
	hpedrmsV_distribution->Fill( fPedRMSV[istrip] );
	hpedmeanV_by_strip->SetBinContent( istrip+1, fPedestalV[istrip] );
	hpedrmsV_by_strip->SetBinContent( istrip+1, fPedRMSV[istrip] );
      }
      
    }
    
    
    hrawADCs_by_stripU = new TH2D( histname.Format( "hrawADCs_by_stripU_%s", detname.Data() ), "Raw ADCs by U strip number, no corrections",
				   fNstripsU, -0.5, fNstripsU-0.5,
				   2048, -0.5, 4095.5 );
    hrawADCs_by_stripV = new TH2D( histname.Format( "hrawADCs_by_stripV_%s", detname.Data() ), "Raw ADCs by V strip number, no corrections",
				   fNstripsV, -0.5, fNstripsV-0.5,
				   2048, -0.5, 4095.5 );

    hcommonmode_subtracted_ADCs_by_stripU = new TH2D( histname.Format( "hpedestalU_%s", detname.Data() ), "ADCs by U strip number, w/common mode correction, no ped. subtraction",
						      fNstripsU, -0.5, fNstripsU-0.5,
						      2500, -500.0, 4500.0 );
    hcommonmode_subtracted_ADCs_by_stripV = new TH2D( histname.Format( "hpedestalV_%s", detname.Data() ), "ADCs by V strip number, w/common mode correction, no ped. subtraction",
						      fNstripsV, -0.5, fNstripsV-0.5,
						      2500, -500.0, 4500.0 );

    hpedestal_subtracted_ADCs_by_stripU = new TH2D( histname.Format( "hADCpedsubU_%s", detname.Data() ), "Pedestal and common-mode subtracted ADCs by U strip number",
						    fNstripsU, -0.5, fNstripsU-0.5,
						    1000,-500.,500. );
    hpedestal_subtracted_ADCs_by_stripV = new TH2D( histname.Format( "hADCpedsubV_%s", detname.Data() ), "Pedestal and common-mode subtracted ADCs by V strip number",
						    fNstripsV, -0.5, fNstripsV-0.5,
						    1000,-500.,500. );

    hpedestal_subtracted_rawADCs_by_stripU = new TH2D( histname.Format( "hrawADCpedsubU_%s", detname.Data() ), "ADCs by U strip, ped-subtracted, no common-mode correction",
						       fNstripsU, -0.5, fNstripsU-0.5,
						       2500,-500.,4500. );
    hpedestal_subtracted_rawADCs_by_stripV = new TH2D( histname.Format( "hrawADCpedsubV_%s", detname.Data() ), "ADCs by V strip, ped-subtracted, no common-mode correction",
						       fNstripsV, -0.5, fNstripsV-0.5,
						       2500,-500.,4500. );

    hpedestal_subtracted_rawADCsU = new TH1D( histname.Format( "hrawADCpedsubU_allstrips_%s", detname.Data() ), "distribution of ped-subtracted U strip ADCs w/o common-mode correction",
					      2500, -500.,4500. );
    hpedestal_subtracted_rawADCsV = new TH1D( histname.Format( "hrawADCpedsubV_allstrips_%s", detname.Data() ), "distribution of ped-subtracted V strip ADCs w/o common-mode correction",
					      2500, -500.,4500. );

    hpedestal_subtracted_ADCsU = new TH1D( histname.Format( "hADCpedsubU_allstrips_%s", detname.Data() ), "distribution of ped-subtracted U strip ADCs w/common-mode correction",
					      1000, -500.,500. );
    hpedestal_subtracted_ADCsV = new TH1D( histname.Format( "hADCpedsubV_allstrips_%s", detname.Data() ), "distribution of ped-subtracted V strip ADCs w/common-mode correction",
					      1000, -500.,500. );

    int nAPVs_U = fNstripsU/fN_APV25_CHAN;
    hcommonmode_mean_by_APV_U = new TH2D( histname.Format( "hCommonModeMean_by_APV_U_%s", detname.Data() ), "distribution of common-mode means for U strip pedestal data",
					  nAPVs_U, -0.5, nAPVs_U-0.5,  
					  2048, -0.5, 4095.5 );
    int nAPVs_V = fNstripsV/fN_APV25_CHAN;
    hcommonmode_mean_by_APV_V = new TH2D( histname.Format( "hCommonModeMean_by_APV_V_%s", detname.Data() ), "distribution of common-mode means for V strip pedestal data",
					  nAPVs_V, -0.5, nAPVs_V-0.5,
					  2048, -0.5, 4095.5 );

    fPedHistosInitialized = true;
    
    // Uncomment these later if you want them:
    // hrawADCs_by_strip_sampleU = new TClonesArray( "TH2D", fN_MPD_TIME_SAMP );
    // hrawADCs_by_strip_sampleV = new TClonesArray( "TH2D", fN_MPD_TIME_SAMP );

    // hcommonmode_subtracted_ADCs_by_strip_sampleU = new TClonesArray( "TH2D", fN_MPD_TIME_SAMP );
    // hcommonmode_subtracted_ADCs_by_strip_sampleV = new TClonesArray( "TH2D", fN_MPD_TIME_SAMP );

    // hpedestal_subtracted_ADCs_by_strip_sampleU = new TClonesArray( "TH2D", fN_MPD_TIME_SAMP );
    // hpedestal_subtracted_ADCs_by_strip_sampleV = new TClonesArray( "TH2D", fN_MPD_TIME_SAMP );
    
    // for( int isamp = 0; isamp<fN_MPD_TIME_SAMP; isamp++ ){
    //   new( (*hrawADCs_by_strip_sampleU)[isamp] ) TH2D( histname.Format( "hrawADCU_%s_sample%d", detname.Data(), isamp ), "Raw U ADCs by strip and sample",
    // 						       fNstripsU, -0.5, fNstripsU-0.5,
    // 						       1024, -0.5, 4095.5 );
    //   new( (*hrawADCs_by_strip_sampleV)[isamp] ) TH2D( histname.Format( "hrawADCV_%s_sample%d", detname.Data(), isamp ), "Raw V ADCs by strip and sample",
    // 						       fNstripsV, -0.5, fNstripsV-0.5,
    // 						       1024, -0.5, 4095.5 );

    //   new( (*hcommonmode_subtracted_ADCs_by_strip_sampleU)[isamp] ) TH2D( histname.Format( "hpedestalU_%s_sample%d", detname.Data(), isamp ), "Pedestals by strip and sample",
    // 									  fNstripsU, -0.5, fNstripsU-0.5,
    // 									  1500, -500.0, 1000.0 );
    //   new( (*hcommonmode_subtracted_ADCs_by_strip_sampleV)[isamp] ) TH2D( histname.Format( "hpedestalV_%s_sample%d", detname.Data(), isamp ), "Pedestals by strip and sample",
    // 									  fNstripsV, -0.5, fNstripsV-0.5,
    // 									  1500, -500.0, 1000.0 );

    //   new( (*hpedestal_subtracted_ADCs_by_strip_sampleU)[isamp] ) TH2D( histname.Format( "hADCpedsubU_%s_sample%d", detname.Data(), isamp ), "Pedestal-subtracted ADCs by strip and sample",
    // 									fNstripsU, -0.5, fNstripsU-0.5,
    // 									1000, -500.0, 500.0 );
    //   new( (*hpedestal_subtracted_ADCs_by_strip_sampleV)[isamp] ) TH2D( histname.Format( "hADCpedsubV_%s_sample%d", detname.Data(), isamp ), "Pedestal-subtracted ADCs by strip and sample",
    // 									fNstripsV, -0.5, fNstripsV-0.5,
    // 									1000, -500.0, 500.0 );
      
    // }
    
  }

  //std::cout << "fCommonModePlotsInitialized = " << fCommonModePlotsInitialized << std::endl;
  if( fMakeCommonModePlots && !fCommonModePlotsInitialized ){

    //std::cout << "SBSGEMModule::Begin: making common-mode histograms for module " << GetName() << std::endl;
    //U strips:
   
    fCommonModeDistU = new TH2D( histname.Format( "hcommonmodeU_%s", detname.Data() ), "U strips: Common mode - common-mode mean(iAPV) vs APV", fNAPVs_U, -0.5, fNAPVs_U-0.5, 500, -1000.0,1000.0 );
    
    
    fCommonModeDistU_Sorting = new TH2D( histname.Format( "hcommonmodeU_sorting_%s", detname.Data() ), "U strips: Common mode - common-mode mean(iAPV) vs APV card, Sorting Method", fNAPVs_U, -0.5, fNAPVs_U-0.5, 500, -1000.0,1000.0 );
    
    
    fCommonModeDistU_Danning = new TH2D( histname.Format( "hcommonmodeU_danning_%s", detname.Data() ), "U strips: Common mode - common-mode mean(iAPV) vs APV card, Danning Method", fNAPVs_U, -0.5, fNAPVs_U-0.5, 500, -1000.0,1000.0 );
    
   
    fCommonModeDiffU = new TH2D( histname.Format( "hcommonmodeU_diff_%s", detname.Data() ), "U strips: Common mode (Sorting) - Common mode (Danning) vs APV card", fNAPVs_U, -0.5, fNAPVs_U-0.5, 250, -25.0, 25.0 );
    

    //V strips:
    
    fCommonModeDistV = new TH2D( histname.Format( "hcommonmodeV_%s", detname.Data() ), "V strips: Common mode - common-mode mean(iAPV) vs APV card", fNAPVs_V, -0.5, fNAPVs_V-0.5, 500, -1000.0,1000.0 );
    
    
    fCommonModeDistV_Sorting = new TH2D( histname.Format( "hcommonmodeV_sorting_%s", detname.Data() ), "V strips: Common mode - common-mode mean(iAPV) vs APV card, Sorting Method", fNAPVs_V, -0.5, fNAPVs_V-0.5, 500, -1000.0,1000.0 );
    
    
    fCommonModeDistV_Danning = new TH2D( histname.Format( "hcommonmodeV_danning_%s", detname.Data() ), "V strips: Common mode - common-mode mean(iAPV) vs APV card, Danning Method", fNAPVs_V, -0.5, fNAPVs_V-0.5, 500, -1000.0,1000.0 );
    
    
    fCommonModeDiffV = new TH2D( histname.Format( "hcommonmodeV_diff_%s", detname.Data() ), "V strips: Common mode (Sorting) - Common mode (Danning) vs APV card", fNAPVs_V, -0.5, fNAPVs_V-0.5, 250, -25.0,25.0 );
    

    // fCommonModeDistU->Print();
    // fCommonModeDistU_Sorting->Print();
    // fCommonModeDistU_Danning->Print();
    // fCommonModeDiffU->Print();

    // fCommonModeDistV->Print();
    // fCommonModeDistV_Sorting->Print();
    // fCommonModeDistV_Danning->Print();
    // fCommonModeDiffV->Print();
    fCommonModePlotsInitialized = true;
  }

  //  std::cout << "fCommonModePlotsInitialized = " << fCommonModePlotsInitialized << std::endl;
    
  return 0;
}

void SBSGEMModule::PrintPedestals( std::ofstream &dbfile, std::ofstream &dbfile_CM, std::ofstream &daqfile, std::ofstream &daqfile_cmr ){
  //The first argument is a file in the format expected by the database,
  //The second argument is a file in the format expected by the DAQ:

  //Step 1: we need to extract pedestal mean and rms by channel, and also grab
  // "common-mode min" and "common-mode max" values:
  // we can simply use the existing arrays to store the values:

  std::cout << "[SBSGEMModule::PrintPedestals]: module " << GetName() << std::endl;

  hpedmeanU_distribution->Reset();
  hpedmeanV_distribution->Reset();
  hpedrmsU_distribution->Reset();
  hpedrmsV_distribution->Reset();

  hpedmeanU_by_strip->Reset();
  hpedmeanV_by_strip->Reset();
  hpedrmsU_by_strip->Reset();
  hpedrmsV_by_strip->Reset();
  
  //start with pedestals;
  for( int iu = 0; iu<fNstripsU; iu++ ){
    TH1D *htemp = hcommonmode_subtracted_ADCs_by_stripU->ProjectionY( "htemp", iu+1,iu+1 );
    
    fPedestalU[iu] = htemp->GetMean();

    //htemp here represents the individual sample noise width, but the threshold is applied on the average of the six (or other number) time samples:
    // sigma(average) = sigma(1 sample)/sqrt(number of samples):
    fPedRMSU[iu] = htemp->GetRMS() / sqrt( double( fN_MPD_TIME_SAMP ) ); 

    hpedmeanU_distribution->Fill( fPedestalU[iu] );
    hpedrmsU_distribution->Fill( fPedRMSU[iu] );

    hpedmeanU_by_strip->SetBinContent( iu+1, fPedestalU[iu] );
    hpedrmsU_by_strip->SetBinContent( iu+1, fPedRMSU[iu] );
    
    htemp->Delete();
  }

  for( int iv = 0; iv<fNstripsV; iv++ ){
    TH1D *htemp = hcommonmode_subtracted_ADCs_by_stripV->ProjectionY( "htemp", iv+1, iv+1 );
    fPedestalV[iv] = htemp->GetMean();
    fPedRMSV[iv] = htemp->GetRMS() / sqrt( double( fN_MPD_TIME_SAMP ) );

    hpedmeanV_distribution->Fill( fPedestalV[iv] );
    hpedrmsV_distribution->Fill( fPedRMSV[iv] );

    hpedmeanV_by_strip->SetBinContent( iv+1, fPedestalV[iv] );
    hpedrmsV_by_strip->SetBinContent( iv+1, fPedRMSV[iv] );
    
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

  header.Form( "%s.%s.%s.commonmode_meanU = ", appname.Data(), detname.Data(), modname.Data() );

  dbfile_CM << std::endl << header << std::endl;
  for( int iAPV = 0; iAPV<nAPVsU; iAPV++ ){
    TString sentry;
    sentry.Form( "  %15.5g ", commonmode_meanU[iAPV] );
    dbfile_CM << sentry;

    if( (iAPV+1) % 16 == 0 ) dbfile_CM << std::endl;
  }

  header.Form( "%s.%s.%s.commonmode_rmsU = ", appname.Data(), detname.Data(), modname.Data() );

  dbfile_CM << std::endl << header << std::endl;
  for( int iAPV = 0; iAPV<nAPVsU; iAPV++ ){
    TString sentry;
    sentry.Form( "  %15.5g ", commonmode_rmsU[iAPV] );
    dbfile_CM << sentry;

    if( (iAPV+1) % 16 == 0 ) dbfile_CM << std::endl;
  }
  
  
  for( int iAPV = 0; iAPV<nAPVsV; iAPV++ ){
    TH1D *htemp = hcommonmode_mean_by_APV_V->ProjectionY("htemp", iAPV+1, iAPV+1 );

    commonmode_meanV[iAPV] = htemp->GetMean();
    commonmode_rmsV[iAPV] = htemp->GetRMS();
  }

  header.Form( "%s.%s.%s.commonmode_meanV = ", appname.Data(), detname.Data(), modname.Data() );

  dbfile_CM << std::endl << header << std::endl;
  for( int iAPV = 0; iAPV<nAPVsV; iAPV++ ){
    TString sentry;
    sentry.Form( "  %15.5g ", commonmode_meanV[iAPV] );
    dbfile_CM << sentry;

    if( (iAPV+1) % 16 == 0 ) dbfile_CM << std::endl;
  }

  header.Form( "%s.%s.%s.commonmode_rmsV = ", appname.Data(), detname.Data(), modname.Data() );

  dbfile_CM << std::endl << header << std::endl;
  for( int iAPV = 0; iAPV<nAPVsV; iAPV++ ){
    TString sentry;
    sentry.Form( "  %15.5g ", commonmode_rmsV[iAPV] );
    dbfile_CM << sentry;

    if( (iAPV+1) % 16 == 0 ) dbfile_CM << std::endl;
  }
  

  dbfile_CM << std::endl << std::endl;
  
  //That takes care of the database file. For the "DAQ" file, we need to organize things by APV card. For this we can loop over the MPDmap:
  
  for( auto iapv = fMPDmap.begin(); iapv != fMPDmap.end(); iapv++ ){
    int crate = iapv->crate;
    int slot = iapv->slot;
    int mpd = iapv->mpd_id;
    int adc_ch = iapv->adc_id;
    int pos = iapv->pos;
    int invert = iapv->invert;
    int axis = iapv->axis;

    
    daqfile << "APV "
	    << std::setw(16) << crate 
	    << std::setw(16) << slot 
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
		<< std::setw(12) << slot
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
    TString appname = (static_cast<THaDetector *>(GetParent()) )->GetApparatus()->GetName();
    appname.ReplaceAll(".","_");
    appname += "_";
    TString detname = GetParent()->GetName();
    detname.Prepend(appname);
    detname.ReplaceAll(".","_");
    detname += "_";
    detname += GetName();
  
    if( fhdidhitx != NULL && fhshouldhitx != NULL ){ //Create efficiency histograms and write to the ROOT file:
      TH1D *hefficiency_vs_x = new TH1D(*fhdidhitx);
      hefficiency_vs_x->SetName( histname.Format( "hefficiency_vs_x_%s", detname.Data() ) );
      hefficiency_vs_x->SetTitle( histname.Format( "Track-based efficiency vs x, module %s", GetName() ) );
      hefficiency_vs_x->Divide( fhshouldhitx );
      hefficiency_vs_x->Write( 0, kOverwrite );
      hefficiency_vs_x->Delete();
    }

    if( fhdidhity != NULL && fhshouldhity != NULL ){ //Create efficiency histograms and write to the ROOT file:
      TH1D *hefficiency_vs_y = new TH1D(*fhdidhity);
      hefficiency_vs_y->SetName( histname.Format( "hefficiency_vs_y_%s", detname.Data() ) );
      hefficiency_vs_y->SetTitle( histname.Format( "Track-based efficiency vs y, module %s", GetName() ) );
      hefficiency_vs_y->Divide( fhshouldhity );
      hefficiency_vs_y->Write( 0, kOverwrite);
      hefficiency_vs_y->Delete();
    }

    if( fhdidhitxy != NULL && fhshouldhitxy != NULL ){ //Create efficiency histograms and write to the ROOT file:
      TH2D *hefficiency_vs_xy = new TH2D(*fhdidhitxy);
      hefficiency_vs_xy->SetName( histname.Format( "hefficiency_vs_xy_%s", detname.Data() ) );
      hefficiency_vs_xy->SetTitle( histname.Format( "Track-based efficiency vs x and y, module %s", GetName() ) );
      hefficiency_vs_xy->Divide( fhshouldhitxy );
      hefficiency_vs_xy->Write( 0, kOverwrite );
      hefficiency_vs_xy->Delete();
    }
  
    if( fhdidhitx != NULL  ) fhdidhitx->Write(fhdidhitx->GetName(), kOverwrite );
    if( fhdidhity != NULL  ) fhdidhity->Write(fhdidhity->GetName(), kOverwrite );
    if( fhdidhitxy != NULL  ) fhdidhitxy->Write(fhdidhitxy->GetName(), kOverwrite );

    if( fhshouldhitx != NULL  ) fhshouldhitx->Write(fhshouldhitx->GetName(), kOverwrite );
    if( fhshouldhity != NULL  ) fhshouldhity->Write(fhshouldhity->GetName(), kOverwrite );
    if( fhshouldhitxy != NULL  ) fhshouldhitxy->Write(fhshouldhitxy->GetName(), kOverwrite );
  }

  if( fPedestalMode || fMakeCommonModePlots ){ //write out pedestal histograms, print out pedestals in the format needed for both database and DAQ:

    //The channel map and pedestal information are specific to a GEM module.
    //But we want ONE database file and ONE DAQ file for the entire tracker.
    //So probably the best way to proceed is to have the  
    hpedrmsU_distribution->Write(0,kOverwrite);
    hpedmeanU_distribution->Write(0,kOverwrite);
    hpedrmsV_distribution->Write(0,kOverwrite);
    hpedmeanV_distribution->Write(0,kOverwrite);

    hpedrmsU_by_strip->Write(0,kOverwrite);
    hpedmeanU_by_strip->Write(0,kOverwrite);
    hpedrmsV_by_strip->Write(0,kOverwrite);
    hpedmeanV_by_strip->Write(0,kOverwrite);
    
    hrawADCs_by_stripU->Write(0,kOverwrite);
    hrawADCs_by_stripV->Write(0,kOverwrite);
    hcommonmode_subtracted_ADCs_by_stripU->Write(0,kOverwrite);
    hcommonmode_subtracted_ADCs_by_stripV->Write(0,kOverwrite);
    hpedestal_subtracted_ADCs_by_stripU->Write(0,kOverwrite);
    hpedestal_subtracted_ADCs_by_stripV->Write(0,kOverwrite);
    hpedestal_subtracted_rawADCs_by_stripU->Write(0,kOverwrite);
    hpedestal_subtracted_rawADCs_by_stripV->Write(0,kOverwrite);

    hpedestal_subtracted_rawADCsU->Write(0,kOverwrite);
    hpedestal_subtracted_rawADCsV->Write(0,kOverwrite);
    hpedestal_subtracted_ADCsU->Write(0,kOverwrite);
    hpedestal_subtracted_ADCsV->Write(0,kOverwrite);

    hcommonmode_mean_by_APV_U->Write(0,kOverwrite);
    hcommonmode_mean_by_APV_V->Write(0,kOverwrite);

    // hrawADCs_by_strip_sampleU->Write();
    // hrawADCs_by_strip_sampleV->Write();
    // hcommonmode_subtracted_ADCs_by_strip_sampleU->Write();
    // hcommonmode_subtracted_ADCs_by_strip_sampleV->Write();
    // hpedestal_subtracted_ADCs_by_strip_sampleU->Write();
    // hpedestal_subtracted_ADCs_by_strip_sampleV->Write();
  }

  if ( fMakeCommonModePlots ){
    fCommonModeDistU->Write(0,kOverwrite);
    fCommonModeDistU_Sorting->Write(0,kOverwrite);
    fCommonModeDistU_Danning->Write(0,kOverwrite);
    fCommonModeDiffU->Write(0,kOverwrite);

    fCommonModeDistV->Write(0,kOverwrite);
    fCommonModeDistV_Sorting->Write(0,kOverwrite);
    fCommonModeDistV_Danning->Write(0,kOverwrite);
    fCommonModeDiffV->Write(0,kOverwrite);
  }
    
  return 0;
}

//utility method to calculate correlation coefficient of U and V samples: 
Double_t SBSGEMModule::CorrCoeff( int nsamples, std::vector<double> Usamples, std::vector<double> Vsamples ){
  Double_t sumu=0.0, sumv=0.0, sumu2=0.0, sumv2=0.0, sumuv=0.0;

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
  Int_t RstripNb = APVMAP[fAPVmapping][rawstrip];
  RstripNb = RstripNb + (127-2*RstripNb)*invert;
  Int_t RstripPos = RstripNb + 128*pos;

  if( fIsMC ){
    return rawstrip + 128*pos;
  }
  
  return RstripPos;
}

void SBSGEMModule::filter_1Dhits(SBSGEM::GEMaxis_t axis){

  if( fFiltering_flag1D < 0 ) return; //flag < 0 means don't filter 1D clusters at all
  // flag = 0 means use a "soft" filter (only reject if at least one other cluster passed)
  // flag > 0 means use a "hard" filter (reject failing clusters no matter what)
  //First filter on cluster ADC sum:
  int ngood = 0;
  int nclust = (axis == SBSGEM::kUaxis) ? fNclustU : fNclustV; 

  std::vector<sbsgemcluster_t> &clusters =  (axis == SBSGEM::kUaxis) ? fUclusters : fVclusters;

  bool passed[nclust];
  
  for( int ipass=0; ipass<2; ipass++ ){
    if( ipass == 0 ) ngood = 0;
    for( int icl=0; icl<nclust; icl++ ){

      //On the first pass, determine which clusters passed the criterion and count the number of good clusters:
      if( ipass == 0 ){
	passed[icl] = clusters[icl].keep && clusters[icl].clusterADCsum >= fThresholdClusterSum;
	if( passed[icl] ) ngood++;
      }

      //on the second pass, if at least one good cluster was found passing the criterion, we set "keep" for all others to false:
      if( ipass == 1 && !passed[icl] && (ngood > 0 || fFiltering_flag1D > 0 ) ){
	clusters[icl].keep = false;
      }
    }
  }

  //Second, filter on cluster size (we may not actually want to filter on cluster size at this stage):
  ngood = 0;

  for( int ipass=0; ipass<2; ipass++ ){
    if( ipass == 0 ) ngood = 0;

    for( int icl=0; icl<nclust; icl++ ){
      if( ipass == 0 ){
	passed[icl] = clusters[icl].keep && clusters[icl].nstrips >= 2;
	if( passed[icl] ) ngood++;
      }

      if( ipass == 1 && !passed[icl] && (ngood > 0 || fFiltering_flag1D > 0 ) ){
	clusters[icl].keep = false;
      }
    }
  }
    
  
  
  
}

void SBSGEMModule::filter_2Dhits(){
  //Here we will initially filter only based on time U/V time difference, ADC asymmetry, and perhaps correlation coefficient:

  if( fFiltering_flag2D < 0 ) return;
  
  //First U/V time difference:
  bool passed[fN2Dhits];
  int ngood = 0;
  for( int ipass=0; ipass<2; ipass++ ){
    if( ipass == 0 ) ngood = 0;
    for( int ihit=0; ihit<fN2Dhits; ihit++ ){
      if( ipass == 0 ){
	passed[ihit] = fHits[ihit].keep && fabs( fHits[ihit].tdiff ) <= fTimeCutUVdiff;
	if( passed[ihit] ) ngood++;
      }

      if( ipass == 1 && !passed[ihit] && ( ngood > 0 || fFiltering_flag2D > 0 ) ){
	fHits[ihit].keep = false;
      }
    }
  }

  //Second: Cluster Correlation Coefficient:
  ngood = 0;
  for( int ipass=0; ipass<2; ipass++ ){
    if( ipass == 0 ) ngood = 0;
    for( int ihit=0; ihit<fN2Dhits; ihit++ ){
      if( ipass == 0 ){
	passed[ihit] = fHits[ihit].keep &&  fHits[ihit].corrcoeff_clust >= fCorrCoeffCut;
	if( passed[ihit] ) ngood++;
      }

      if( ipass == 1 && !passed[ihit] && ( ngood > 0 || fFiltering_flag2D > 0 ) ){
	fHits[ihit].keep = false;
      }
    }
  }
  
  //Third: ADC asymmetry:
  ngood = 0;
  for( int ipass=0; ipass<2; ipass++ ){
    if( ipass == 0 ) ngood = 0;
    for( int ihit=0; ihit<fN2Dhits; ihit++ ){
      if( ipass == 0 ){
	passed[ihit] = fHits[ihit].keep && fabs( fHits[ihit].ADCasym ) <= fADCasymCut;
	if( passed[ihit] ) ngood++;
      }

      if( ipass == 1 && !passed[ihit] && ( ngood > 0 || fFiltering_flag2D > 0 ) ){
	fHits[ihit].keep = false;
      }
    }
  }

  //other criteria could include correlation coefficient, time sample peaking, etc. 
  
}

double SBSGEMModule::GetCommonMode( UInt_t isamp, Int_t flag, const mpdmap_t &apvinfo ){
  if( isamp > fN_MPD_TIME_SAMP ) return 0;

  if( flag == 0 ){ //Sorting method: doesn't actually use the apv info: 
    vector<double> sortedADCs(fN_APV25_CHAN);
	    
    for( int ihit=0; ihit<fN_APV25_CHAN; ihit++ ){
      int iraw = isamp + fN_MPD_TIME_SAMP * ihit;	    
      sortedADCs[ihit] = fPedSubADC_APV[ iraw ];
    }
	    
    std::sort( sortedADCs.begin(), sortedADCs.end() );
	    
    //   commonMode[isamp] = 0.0;
    double cm_temp = 0.0;
    int stripcount=0;

    
    for( int k=fCommonModeNstripRejectLow; k<fN_APV25_CHAN-fCommonModeNstripRejectHigh; k++ ){
      cm_temp += sortedADCs[k];
      stripcount++;
    }
    return  cm_temp/double(stripcount);
  } else { //Danning method: requires apv info for cm-mean and cm-rms values:
    int iAPV = apvinfo.pos;
    double cm_mean = ( apvinfo.axis == SBSGEM::kUaxis ) ? fCommonModeMeanU[iAPV] : fCommonModeMeanV[iAPV];
    double cm_rms = ( apvinfo.axis == SBSGEM::kUaxis ) ? fCommonModeRMSU[iAPV] : fCommonModeRMSV[iAPV];

    //TODO: allow to use a different parameter than the one used for
    // zero-suppression:
    double cm_min = cm_mean - fZeroSuppressRMS*cm_rms;
    double cm_max = cm_mean + fZeroSuppressRMS*cm_rms;

    double cm_temp = 0.0;
    
    //for now, hard-code 3 iterations:
    for( int iter=0; iter<fCommonModeNumIterations; iter++ ){
      int nstripsinrange=0;
      double sumADCinrange=0.0;
      //double sum2ADCinrange=0.0;
      for( int ihit=0; ihit<fN_APV25_CHAN; ihit++ ){
	int iraw=isamp + fN_MPD_TIME_SAMP * ihit;
	
	double ADCtemp = fPedSubADC_APV[iraw];
	
	//on iterations after the first iteration, reject strips with signals above nsigma * pedrms:
	double rmstemp = ( apvinfo.axis == SBSGEM::kUaxis ) ? fPedRMSU[fStripAPV[iraw]] : fPedRMSV[fStripAPV[iraw]];

	double mintemp = cm_min;
	double maxtemp = cm_max;
	
	if( iter > 0 ) {
	  maxtemp = cm_temp + fZeroSuppressRMS*rmstemp*fRMS_ConversionFactor; //2.45 = sqrt(6), don't want to calculate sqrt every time
	  mintemp = 0.0;
	}
	
	if( ADCtemp >= mintemp && ADCtemp <= maxtemp ){
	  nstripsinrange++;
	  sumADCinrange += ADCtemp;
	  //sum2ADCinrange += pow(ADCtemp,2);
	}
      }
      
      //double rmsinrange = cm_rms;
      
      //TO-DO: don't hard-code the minimum strip count. 
      
      if( nstripsinrange >= fCommonModeMinStripsInRange ){ //require minimum 10 strips in range:
	cm_temp = sumADCinrange/double(nstripsinrange);
	//rmsinrange = sqrt( sum2ADCinrange/double(nstripsinrange) - pow(commonMode[isamp],2) );
	//cm_max = commonMode[isamp] + fZeroSuppressRMS*std::min( rmsinrange, cm_rms );
      } else if( iter==0 ){ //not enough strips on FIRST iteration, use mean from sorting-method:
	// std::cout << "Warning: fewer than ten strips in range on first iteration for common-mode calculation, crate, slot, MPD, ADC, cm_min, cm_max, n strips in range, cm_mean = "
	// 	  << apvinfo.crate << ", " << apvinfo.slot << ", " << apvinfo.mpd_id << ", " << apvinfo.adc_id << ", "
	// 	  << cm_min << ", " << cm_max << ", " << nstripsinrange << ", "
	// 	  << cm_mean << ", defaulting to sorting method" << std::endl;

	return GetCommonMode( isamp, 0, apvinfo );
      }

      // std::cout << "iteration " << iter << ": isamp, crate, slot, mpd, adc, cm_mean, cm_min, cm_max, nstrips in range, cm calc, cm sorting method = " << isamp << ", "
      // 		<< apvinfo.crate << ", " << apvinfo.slot << ", " << apvinfo.mpd_id << ", " << apvinfo.adc_id << ", "
      // 		<< cm_mean << ", " << cm_min << ", " << cm_max << ", " << nstripsinrange << ", " << cm_temp << ", "
      // 		<< GetCommonMode( isamp, 0, apvinfo ) << std::endl;
    } //loop over iterations for "Danning method" CM calculation

    //std::cout << std::endl;
    return cm_temp;
  }
}

void SBSGEMModule::InitAPVMAP(){
  APVMAP.clear();

  APVMAP[SBSGEM::kINFN].resize(fN_APV25_CHAN);
  APVMAP[SBSGEM::kUVA_XY].resize(fN_APV25_CHAN);
  APVMAP[SBSGEM::kUVA_UV].resize(fN_APV25_CHAN);
  APVMAP[SBSGEM::kMC].resize(fN_APV25_CHAN);

  for( UInt_t i=0; i<fN_APV25_CHAN; i++ ){
    Int_t strip1 = 32*(i%4) + 8*(i/4) - 31*(i/16);
    Int_t strip2 = strip1 + 1 + strip1 % 4 - 5 * ( ( strip1/4 ) % 2 );
    Int_t strip3 = ( strip2 % 2 == 0 ) ? strip2/2 + 32 : ( (strip2<64) ? (63 - strip2)/2 : 127 + (65-strip2)/2 ); 
    APVMAP[SBSGEM::kINFN][i] = strip1; 
    APVMAP[SBSGEM::kUVA_XY][i] = strip2;
    APVMAP[SBSGEM::kUVA_UV][i] = strip3;
    APVMAP[SBSGEM::kMC][i] = i;
  }

  //Print out to test:
  // //std::cout << "INFN X/Y APV mapping: " << std::endl;
  // for( UInt_t i=0; i<fN_APV25_CHAN; i++ ){
  //   std::cout << std::setw(8) << APVMAP[SBSGEM::kINFN][i] << ", ";
  //   if( (i+1) % 10 == 0 ) std::cout << endl;
  // }
  // std::cout << endl;


  // //Print out to test:
  // std::cout << "UVA X/Y APV mapping: " << std::endl;
  // for( UInt_t i=0; i<fN_APV25_CHAN; i++ ){
  //   std::cout << std::setw(8) << APVMAP[SBSGEM::kUVA_XY][i] << ", ";
  //   if( (i+1) % 10 == 0 ) std::cout << endl;
  // }
  // std::cout << endl;

  // //Print out to test:
  // std::cout << "UVA U/V APV mapping: " << std::endl;
  // for( UInt_t i=0; i<fN_APV25_CHAN; i++ ){
  //   std::cout << std::setw(8) << APVMAP[SBSGEM::kUVA_UV][i] << ", ";
  //   if( (i+1) % 10 == 0 ) std::cout << endl;
  // }
  // std::cout << endl;
  
}
