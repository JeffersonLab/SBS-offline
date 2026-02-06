///////////////////////////////////////////////////////////////////////////////
//
// SBSHCal
//
///////////////////////////////////////////////////////////////////////////////
#include "SBSHCal.h"
#include <iostream>
#include "THaEvData.h"

ClassImp(SBSHCal);

/*
 * SBSHCal constructor.
 *
 * Specify SBSCalorimeter to use both TDC and ADC Multi-samples
 */
SBSHCal::SBSHCal( const char* name, const char* description,
    THaApparatus* apparatus ) : SBSCalorimeter(name,description,apparatus)
{
  SetModeADC(SBSModeADC::kWaveform);
  SetModeTDC(SBSModeTDC::kTDCSimple);
  SetDisableRefTDC(true);
  fWithLED = true;

  //Default values for time-based cuts for best cluster selection:
  fRequireTDCGoodCluster = false;
  fAtimeMinGoodCluster = -10000.;
  fAtimeMaxGoodCluster = 10000.;
  fRefADCtimeGoodCluster = 0.0;
}

///////////////////////////////////////////////////////////////////////////////
/// Read SBSHCal Database
Int_t SBSHCal::ReadDatabase( const TDatime& date )
{
  // We can use this name here for logs
  static const char* const here = "ReadDatabase()";

  FILE* file = OpenFile( date );
  if( !file ) return kFileError;
  Int_t err = kOK;
  
  if(fWithLED) {
   
  /*  
    // Read in required geometry variables, which include fOrigin and fSize
    Int_t err = ReadGeometry( file, date, true );
    if( err ) {
      fclose(file);
      return err;
    }
    */
  
    std::vector<Int_t> ledmap;
    DBRequest led_request[] = {
      { "ledmap", &ledmap, kIntV, 2, false }, ///< ledmap of LED
      {0}
    };
    err = LoadDB( file, date, led_request, fPrefix );
    //fclose(file);
    if(err) {
	return err;
    }
    if(ledmap.size()<2) {
      Error(Here(here), "Need crate slot for LED");
      return kInitError;
    }
    fLEDCrate = ledmap[0];
    fLEDSlot  = ledmap[1];
  }

  Int_t tdc_flag = fRequireTDCGoodCluster ? 1 : 0;
  
  DBRequest thisreq[] = {
    { "adctime_min_clusterselect", &fAtimeMinGoodCluster, kDouble, 0, 1},
    { "adctime_max_clusterselect", &fAtimeMaxGoodCluster, kDouble, 0, 1},
    { "requireTDC_clusterselect", &tdc_flag, kInt, 0, 1},
    { 0 }
  };

  err = LoadDB( file, date, thisreq, fPrefix );
  
  fRequireTDCGoodCluster = tdc_flag != 0 ? true : false;
  
  if( fAtimeMinGoodCluster > fAtimeMaxGoodCluster ){ // wrong order:
    double min = fAtimeMaxGoodCluster;
    fAtimeMaxGoodCluster = fAtimeMinGoodCluster;
    fAtimeMinGoodCluster = min;
  }

  fclose(file);
  
  if(err) {  
    return err;
  }
  
  return SBSCalorimeter::ReadDatabase(date);

}


//_____________________________________________________________________________
Int_t SBSHCal::Decode( const THaEvData& evdata )
{
  Int_t err = SBSCalorimeter::Decode(evdata);
  if(fWithLED) {
    UInt_t ihit = evdata.GetNumChan(fLEDCrate,fLEDSlot);
    if(ihit!=2 ) {
      //std::cerr << "ihit=" << ihit << std::endl;
      return 0;
    }
    fLEDBit = evdata.GetData(fLEDCrate,fLEDSlot,1,0);
    fLEDCount = evdata.GetData(fLEDCrate,fLEDSlot,2,0);
  }
  return err;
}
//
Int_t SBSHCal::CoarseProcess(TClonesArray& tracks)
{
  Int_t err = SBSCalorimeter::CoarseProcess(tracks);
  if(err) {
    return err;
  }
  /*Int_t BlockSize = */ SBSCalorimeter::MakeGoodBlocks();

  Int_t ClusSize = SBSCalorimeter::FindClusters();

  return ClusSize;
}
//_____________________________________________________________________________
Int_t SBSHCal::DefineVariables( EMode mode )
{
  // Initialize global variables
  Int_t err = SBSCalorimeter::DefineVariables(mode);
  if(err) {
    return err;
  }

  RVarDef vars[] = {
    { "ledbit", "LEDBit",  "fLEDBit" },
    { "ledcount", "LEDCount",  "fLEDCount" },
    { 0 }
  };

  err = DefineVarsFromList( vars, mode );
  return err;
}

void SBSHCal::Clear( Option_t* opt )
{
  fLEDBit = -1;
  fLEDCount = 0;
  SBSCalorimeter::Clear(opt);
  fRefADCtimeGoodCluster = 0.0;
}
/*
 * Generic SBSHCal destructor
 */
SBSHCal::~SBSHCal()
{
}

Int_t SBSHCal::SelectBestCluster(){ //Default is just highest-energy cluster regardless of timing:
  //This will implement the equivalent of the "highest-energy in-time" algorithm for HCAL best cluster selection:
  if( fNclus <= 0 ) return -1;

  // std::cout << "Calling SBSHCal::SelectBestCluster(), atime (min,max,ref,ibest)=("
  //  	    << fAtimeMinGoodCluster << ", " << fAtimeMaxGoodCluster << ", " << fRefADCtimeGoodCluster
  // 	    << ", " << fBestClusterIndex << ")" << std::endl;
  
  int oldindex = fBestClusterIndex;

  //we're not implementing the filtering algorithm below unless and until it is fully debugged/understood.
  //if( true ) return fBestClusterIndex;
  
  int best=-1;
  std::vector<Bool_t> keep(fNclus,true);

  //"Filter" clusters based on two criteria, always keeping at least one!
  // 1) ADC time of highest-energy block within "good cluster" limits:
  // 2) at least one good TDC hit
  //First: ADC time check:
  int ngood = 0;

  double Told = fClusters[oldindex]->GetAtime();
  double Eold = fClusters[oldindex]->GetE();
  double Tbest = Told;
  double Ebest = Eold;
  
  for( int ipass=0; ipass<2; ipass++ ){
    for( int iclus=0; iclus<fNclus; iclus++ ){
      double Tcheck = fClusters[iclus]->GetAtime();
      bool goodADCtime = fAtimeMinGoodCluster < Tcheck && Tcheck < fAtimeMaxGoodCluster;
      if( ipass == 0 && keep[iclus] && goodADCtime ) ngood++;
      if( ipass > 0 && !goodADCtime && ngood > 0 ) keep[iclus] = false; //If at least one cluster with ADC time in "good" window, reject other clusters
    }
  }

  if( fRequireTDCGoodCluster ){ //Second: filter on presence of at least one TDC hit in the cluster
    ngood = 0;
    for( int ipass=0; ipass<2; ipass++ ){
      for( int iclus=0; iclus<fNclus; iclus++ ){
	bool goodTDC = fClusters[iclus]->GetNgoodTDChits() > 0;
	if( ipass == 0 && keep[iclus] && goodTDC ) ngood++;
	if( ipass > 0 && !goodTDC && ngood > 0 ) keep[iclus] = false; //If at least one cluster with a good TDC hit (and ADC time), reject other clusters
      }
    }
  }

  double Emax = 0.0;
  
  //Finally, keep the highest-energy cluster passing all applicable filtering criteria:
  for( int iclus=0; iclus<fNclus; iclus++ ){
    if( keep[iclus] && fClusters[iclus]->GetE() > Emax ){
      Emax = fClusters[iclus]->GetE();
      best = iclus;
    }
  }

  Tbest = fClusters[best]->GetAtime();
  Ebest = fClusters[best]->GetE();
  
  //if( best != oldindex ) {
    // std::cout << "Changed best cluster index from " << oldindex << " to " << best
    // 	      << ", (Told,Tnew)=(" << Told << ", " << Tbest << "), (Eold, Enew)=("
    // 	      << Eold << ", " << Ebest << ")" << std::endl;

  //}
    
  if( best >= 0 ) fBestClusterIndex = best;
  return fBestClusterIndex;
}
