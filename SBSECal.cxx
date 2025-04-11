///////////////////////////////////////////////////////////////////////////////
//
// SBSECal
//
///////////////////////////////////////////////////////////////////////////////
#include "SBSECal.h"
#include <iostream>
#include "THaEvData.h"
#include "THaApparatus.h"

ClassImp(SBSECal);

/*
 * SBSECal constructor.
 *
 * Specify SBSCalorimeter to use both TDC and ADC Multi-samples
 */
SBSECal::SBSECal( const char* name, const char* description,
    THaApparatus* apparatus ) : SBSCalorimeter(name,description,apparatus)
{
  SetModeADC(SBSModeADC::kWaveform);
  // SetModeTDC(SBSModeTDC::kTDCSimple);
  // SetDisableRefTDC(true);
  // fWithLED = true;

  //Default values for time-based cuts for best cluster selection:
  // fRequireTDCGoodCluster = false;
  fAtimeMinGoodCluster = -10000.;
  fAtimeMaxGoodCluster = 10000.;
  fRefADCtimeGoodCluster = 0.0;

  //Defaults for shower shape corrections:
  fUseShowerShapeCorr = false;

  //Ecal block sizes (m, approx):
  fShowerShapeLxProf = 0.043;
  fShowerShapeLyProf = 0.043;
  //Default binnings for ECAL and profiles:
  fShowerShapeNbinsX = 1;
  fShowerShapeNbinsY = 1;
  fShowerShapeNbinsProf = 100;
  fShowerShapeMomMin = -1.0;
  fShowerShapeMomMax = 1.0;

  fShowerShapeXminProf = -1.5;
  fShowerShapeXmaxProf = 1.5;
  fShowerShapeYminProf = -0.65;
  fShowerShapeYmaxProf = 0.65;

  fShowerProfilesX.clear();
  fShowerProfilesY.clear();
  
}

///////////////////////////////////////////////////////////////////////////////
/// Read SBSECal Database
Int_t SBSECal::ReadDatabase( const TDatime& date )
{
  // We can use this name here for logs
  static const char* const here = "ReadDatabase()";

  FILE* file = OpenFile( date );
  if( !file ) return kFileError;
  Int_t err = kOK;
  
  // if(fWithLED) {
   
  // /*  
  //   // Read in required geometry variables, which include fOrigin and fSize
  //   Int_t err = ReadGeometry( file, date, true );
  //   if( err ) {
  //     fclose(file);
  //     return err;
  //   }
  //   */
  
  //   std::vector<Int_t> ledmap;
  //   DBRequest led_request[] = {
  //     { "ledmap", &ledmap, kIntV, 2, false }, ///< ledmap of LED
  //     {0}
  //   };
  //   err = LoadDB( file, date, led_request, fPrefix );
  //   //fclose(file);
  //   if(err) {
  // 	return err;
  //   }
  //   if(ledmap.size()<2) {
  //     Error(Here(here), "Need crate slot for LED");
  //     return kInitError;
  //   }
  //   fLEDCrate = ledmap[0];
  //   fLEDSlot  = ledmap[1];
  // }

  // Int_t tdc_flag = fRequireTDCGoodCluster ? 1 : 0;
  Int_t useshowershapeflag = fUseShowerShapeCorr ? 1 : 0;

  std::vector<Double_t> xproftemp;
  std::vector<Double_t> yproftemp; 
  
  DBRequest thisreq[] = {
    { "adctime_min_clusterselect", &fAtimeMinGoodCluster, kDouble, 0, 1},
    { "adctime_max_clusterselect", &fAtimeMaxGoodCluster, kDouble, 0, 1},
    { "nbinsxprof", &fShowerShapeNbinsX, kInt, 0, 1, 1 },
    { "nbinsyprof", &fShowerShapeNbinsY, kInt, 0, 1, 1 },
    { "nbinsprof", &fShowerShapeNbinsProf, kInt, 0, 1, 1 },
    { "mom_min", &fShowerShapeMomMin, kDouble, 0, 1, 1 },
    { "mom_max", &fShowerShapeMomMax, kDouble, 0, 1, 1 },
    { "xminprof", &fShowerShapeXminProf, kDouble, 0, 1, 1 },
    { "xmaxprof", &fShowerShapeXmaxProf, kDouble, 0, 1, 1 },
    { "yminprof", &fShowerShapeYminProf, kDouble, 0, 1, 1 },
    { "ymaxprof", &fShowerShapeYmaxProf, kDouble, 0, 1, 1 },
    { "useshowershapecorr", &useshowershapeflag, kInt, 0, 1, 1 },
    { "Lxprof", &fShowerShapeLxProf, kDouble, 0, 1, 1 },
    { "Lyprof", &fShowerShapeLyProf, kDouble, 0, 1, 1 },
    { "xprofiles", &xproftemp, kDoubleV, 0, 1, 1 },
    { "yprofiles", &yproftemp, kDoubleV, 0, 1, 1 },
    // { "requireTDC_clusterselect", &tdc_flag, kInt, 0, 1},
    { 0 }
  };

  err = LoadDB( file, date, thisreq, fPrefix );
  
  // fRequireTDCGoodCluster = tdc_flag != 0 ? true : false;
  
  if( fAtimeMinGoodCluster > fAtimeMaxGoodCluster ){ // wrong order:
    double min = fAtimeMaxGoodCluster;
    fAtimeMaxGoodCluster = fAtimeMinGoodCluster;
    fAtimeMinGoodCluster = min;
  }

  fUseShowerShapeCorr = (useshowershapeflag > 0);

  if( fUseShowerShapeCorr ){ //Then check that the xprofile and yprofile data provided have the correct size:
    if( xproftemp.size() != fShowerShapeNbinsX*fShowerShapeNbinsY*fShowerShapeNbinsProf ||
	yproftemp.size() != fShowerShapeNbinsX*fShowerShapeNbinsY*fShowerShapeNbinsProf ){
      std::cout << "Warning in SBSECal::ReadDatabase for detector " << GetApparatus()->GetName() << "." << GetName() << ": shower x and/or y profiles incorrect size; disabling shower shape corrections (fix database if you want these)!" << std::endl;
      fUseShowerShapeCorr = false; 
    } else {
      // Set up the shower profiles:
      fShowerProfilesX.resize( fShowerShapeNbinsX*fShowerShapeNbinsY );
      fShowerProfilesY.resize( fShowerShapeNbinsX*fShowerShapeNbinsY );
      //ORDER of the shower profiles is y bin first, then x bin
      int i=0;
      for( int binx=0; binx<fShowerShapeNbinsX; binx++ ){
	for( int biny=0; biny<fShowerShapeNbinsY; biny++ ){
	  int bin = biny + fShowerShapeNbinsY*binx;
	  fShowerProfilesX[bin].resize( fShowerShapeNbinsProf );
	  fShowerProfilesY[bin].resize( fShowerShapeNbinsProf );
	  for( int binfrac=0; binfrac<fShowerShapeNbinsProf; binfrac++ ){
	    fShowerProfilesX[bin][binfrac] = xproftemp[i];
	    fShowerProfilesY[bin][binfrac] = yproftemp[i];
	    i++;
	  }
	}
      }
    }
  }
  
  fclose(file);
  
  if(err) {  
    return err;
  }
  
  return SBSCalorimeter::ReadDatabase(date);

}


//_____________________________________________________________________________
Int_t SBSECal::Decode( const THaEvData& evdata )
{
  Int_t err = SBSCalorimeter::Decode(evdata);
  // if(fWithLED) {
  //   UInt_t ihit = evdata.GetNumChan(fLEDCrate,fLEDSlot);
  //   if(ihit!=2 ) {
  //     //std::cerr << "ihit=" << ihit << std::endl;
  //     return 0;
  //   }
  //   fLEDBit = evdata.GetData(fLEDCrate,fLEDSlot,1,0);
  //   fLEDCount = evdata.GetData(fLEDCrate,fLEDSlot,2,0);
  // }
  return err;
}
//
Int_t SBSECal::CoarseProcess(TClonesArray& tracks)
{
  Int_t err = SBSCalorimeter::CoarseProcess(tracks);
  if(err) {
    return err;
  }
  SBSCalorimeter::MakeGoodBlocks();

  Int_t ClusSize = SBSCalorimeter::FindClusters();

  // Modify shower coordinates using shape corrections, if provided.
  if( fUseShowerShapeCorr ){
    CalcShowerCoord();
  }
  
  return ClusSize;
}
//_____________________________________________________________________________
Int_t SBSECal::DefineVariables( EMode mode )
{
  // Initialize global variables
  Int_t err = SBSCalorimeter::DefineVariables(mode);
  if(err) {
    return err;
  }

  // RVarDef vars[] = {
  //   { "ledbit", "LEDBit",  "fLEDBit" },
  //   { "ledcount", "LEDCount",  "fLEDCount" },
  //   { 0 }
  // };

  // err = DefineVarsFromList( vars, mode );
  return err;
}

void SBSECal::Clear( Option_t* opt )
{
  // fLEDBit = -1;
  // fLEDCount = 0;
  SBSCalorimeter::Clear(opt);
}
/*
 * Generic SBSECal destructor
 */
SBSECal::~SBSECal()
{
}

Int_t SBSECal::SelectBestCluster(){ //Default is just highest-energy cluster regardless of timing:
  //This will implement the equivalent of the "highest-energy in-time" algorithm for ECAL best cluster selection:
  if( fNclus <= 0 ) return -1;

  // std::cout << "Calling SBSECal::SelectBestCluster(), atime (min,max,ref,ibest)=("
  //  	    << fAtimeMinGoodCluster << ", " << fAtimeMaxGoodCluster << ", " << fRefADCtimeGoodCluster
  // 	    << ", " << fBestClusterIndex << ")" << std::endl;
  
  int oldindex = fBestClusterIndex;

  //we're not implementing the filtering algorithm below unless and until it is fully debugged/understood.
  //if( true ) return fBestClusterIndex;
  
  int best=-1;
  std::vector<Bool_t> keep(fNclus,true);

  //"Filter" clusters based on two criteria, always keeping at least one!
  // 1) ADC time of highest-energy block within "good cluster" limits:
  // 2) at least one good TDC hit -- Not being used for ECAL [P. Datta - 12/09/24]
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

  // if( fRequireTDCGoodCluster ){ //Second: filter on presence of at least one TDC hit in the cluster
  //   ngood = 0;
  //   for( int ipass=0; ipass<2; ipass++ ){
  //     for( int iclus=0; iclus<fNclus; iclus++ ){
  // 	bool goodTDC = fClusters[iclus]->GetNgoodTDChits() > 0;
  // 	if( ipass == 0 && keep[iclus] && goodTDC ) ngood++;
  // 	if( ipass > 0 && !goodTDC && ngood > 0 ) keep[iclus] = false; //If at least one cluster with a good TDC hit (and ADC time), reject other clusters
  //     }
  //   }
  // }

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

void SBSECal::CalcShowerCoord(){
  if( !fUseShowerShapeCorr ) return; //This flag will be forced to false if
  // the database is not correctly parsed with all required parameters.

  //  std::cout << "SBSECal::CalcShowerCoord(): calculating shower shape corrections" << std::endl;
  
  //loop on all clusters
  int nclust = GetNclust();

  //declare a reference so we don't need to make a local copy (even though it's just a vector of pointers):
  std::vector<SBSCalorimeterCluster*> &clusters = GetClusters();
  
  for( int i=0; i<nclust; i++ ){
    Double_t xsum = 0.0, ysum = 0.0, esum = 0.0;
    int nblk = 0;

    int jmax=-1;
    double Emax = 0.0; //while the first block "should" always be the highest-energy block, let's not ASSUME that!
    double xmax = 0.0, ymax = 0.0;
    //loop on all the blocks and calculate the shower "moments":
    for( int j=0; j<clusters[i]->GetSize(); j++ ){
      double xblk = clusters[i]->GetElement(j)->GetX();
      double yblk = clusters[i]->GetElement(j)->GetY();
      double Eblk = clusters[i]->GetElement(j)->GetE();

      xsum += xblk*Eblk;
      ysum += yblk*Eblk;
      esum += Eblk;

      if( j==0 || Eblk > Emax ){
	jmax = j;
	Emax = Eblk;
	ymax = yblk;
	xmax = xblk;
      }
      nblk++;
    }

    double xmom = (xsum / esum - xmax)/fShowerShapeLxProf;
    double ymom = (ysum / esum - ymax)/fShowerShapeLyProf;

    int binx = int( (xmax - fShowerShapeXminProf)/(fShowerShapeXmaxProf-fShowerShapeXminProf)*double(fShowerShapeNbinsX) );
    int biny = int( (ymax - fShowerShapeYminProf)/(fShowerShapeXmaxProf-fShowerShapeXminProf)*double(fShowerShapeNbinsY) );
    if( binx >= 0 && binx < fShowerShapeNbinsX &&
	biny >= 0 && biny < fShowerShapeNbinsY ){
      int bin = biny + fShowerShapeNbinsY * binx;

      std::vector<Double_t> &fracx = fShowerProfilesX[bin];
      std::vector<Double_t> &fracy = fShowerProfilesY[bin];

      double binwidth_frac = (fShowerShapeMomMax-fShowerShapeMomMin)/double(fShowerShapeNbinsProf);
      
      int binfracx = int( (xmom-fShowerShapeMomMin)/binwidth_frac );
      int binfracy = int( (ymom-fShowerShapeMomMin)/binwidth_frac );

      //For now we'll just linearly interpolate between bins. Later we might want to get fancy and do cubic spline interpolation...

      //In fact, the way the profiles are calculated, the value in each bin is the fraction below the HIGH edge of that bin.
      //So this calculation needs to be adjusted:
      
      if( binfracx >= 0 && binfracx < fShowerShapeNbinsProf &&
	  binfracy >= 0 && binfracy < fShowerShapeNbinsProf ){

	double fracxhigh = fracx[binfracx];
	double fracxlow = fracxhigh;
	if( binfracx > 0 ){
	  fracxlow = fracx[binfracx-1];
	}
	double fracyhigh = fracy[binfracy];
	double fracylow = fracyhigh;
	if( binfracy > 0 ){
	  fracylow = fracy[binfracy-1];
	}

	double lowedgex = fShowerShapeMomMin + binfracx*binwidth_frac;
	double lowedgey = fShowerShapeMomMin + binfracy*binwidth_frac;

	double xinterpfrac = (xmom-lowedgex)/binwidth_frac;
	double yinterpfrac = (ymom-lowedgey)/binwidth_frac;

	double fracxfinal = fracxlow*(1.0-xinterpfrac)+fracxhigh*xinterpfrac;
	double fracyfinal = fracylow*(1.0-yinterpfrac)+fracyhigh*yinterpfrac;

	double xcorrected = (xmax + (fracxfinal - 0.5) * fShowerShapeLxProf);
	double ycorrected = (ymax + (fracyfinal - 0.5) * fShowerShapeLyProf);

	double xold = clusters[i]->GetX();
	double yold = clusters[i]->GetY();
	
	clusters[i]->SetX( xcorrected );
	clusters[i]->SetY( ycorrected );

	// std::cout << "updating cluster coordinates, (iclust,xfrac,yfrac,xold,yold,xnew,ynew)=(" << i << ", "
	// 	  << fracxfinal << ", " << fracyfinal << ", " << xold << ", " << yold
	// 	  << ", " << clusters[i]->GetX() << ", " << clusters[i]->GetY() << ")" << std::endl;

	if( i == fBestClusterIndex ){ //update fMainClus position variables:
	  //This check is in principle unnecessary but do it anyway to be safe:
	  if( fMainclus.x.size() > 0 && fMainclus.y.size() > 0 ){
	    fMainclus.x[0] = xcorrected;
	    fMainclus.y[0] = ycorrected;
	  }
	}
	
      } //end check of xmom and ymom within limits (otherwise we leave the shower coordinates alone)
      
    } //end check that position of cluster seed lies within limits

  } //end loop over clusters

  return;
}
