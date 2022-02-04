///////////////////////////////////////////////////////////////////////////////
//
// SBSBBShower
//
///////////////////////////////////////////////////////////////////////////////
#include "SBSBBShower.h"
#include "SBSCalorimeter.h"
#include <iostream>
#include "THaEvData.h"
#include <iomanip>

using namespace std;
ClassImp(SBSBBShower);

/*
 * SBSBBShower constructor.
 *
 * Specify SBSCalorimeter to use both TDC and ADC Multi-samples
 */
SBSBBShower::SBSBBShower( const char* name, const char* description,
    THaApparatus* apparatus ) : SBSCalorimeter(name,description,apparatus),
  fSearchRegion(0), fSearchRowmin(0), fSearchRowmax(0), fSearchColmin(0),
  fSearchColmax(0)
{
  SetModeADC(SBSModeADC::kWaveform); //< Multi-function ADC
  SetModeTDC(SBSModeTDC::kNone); //< No TDC information
}

//_____________________________________________________________________________
Int_t SBSBBShower::ReadDatabase( const TDatime& date )
{
  cout << "******** Detector " << GetName() << " ReadDatabase ********" << endl;
  //static const char* const here = "ReadDatabase()";
  // Call the parent class ReadDatabase first
  Int_t err = SBSCalorimeter::ReadDatabase(date);
  if(err) {
    return err;
  }
  fIsInit = false;

  FILE* file = OpenFile( date );
  if( !file ) return kFileError;

  std::vector<Double_t> dxyz;
  // Readout components needed by BBShower
  DBRequest config_request[] = {
    { "thr_adc",      &fThrADC,     kDouble,  0, true },
    { "clus_rad",     &fClusRadius, kDouble,   0, true },
    { "mc_data",      &fMCdata,     kInt,     0, true },// flag for MC data
    { "dxdydz",         &dxyz,         kDoubleV, 3 },  // dx and dy block spacings
    { 0 } ///< Request must end in a NULL
  };
  err = LoadDB( file, date, config_request, fPrefix );
  fclose(file);
  if(err) {
    return err;
  }
  
  fClusBlockRadX = Int_t(fClusRadius/dxyz[0]);
  fClusBlockRadY = Int_t(fClusRadius/dxyz[1]);

  if(fMaxNclus>1)fMultClus = true;

  fIsInit = true;
  return 0;
}


Int_t SBSBBShower::DefineVariables( EMode mode )
{
  if( mode == kDefine && fIsSetup ) return kOK;
  // Initialize parent variables first
  Int_t err = SBSCalorimeter::DefineVariables(mode);
  std::cout << " SBS BBShower define variables " << err << std::endl;

  if(err)
    return err;

  // Register variables in global list
  
  if(fMCdata){
    RVarDef varsmc[] = {
      { "e_m_res", "Energy resolution of main cluster",    "fEres" },
      { "x_m_res", "x-position resolution (m) of main cluster", "fXres" },
      { "y_m_res", "y-position resolution  (m) of main cluster", "fYres" },
      { "e_res", "Energy resolution of all clusters",    "fE_cl_res" },
      { "x_res", "x-position resolution (m) of all clusters", "fX_cl_res" },
      { "y_res", "y-position resolution  (m) of all clusters", "fY_cl_res" },
      { 0 }
    };

    err = DefineVarsFromList( varsmc, mode );
    if( err != kOK )
      return err;
  }
  

  return err;
};


//_____________________________________________________________________________
Int_t SBSBBShower::CoarseProcess(TClonesArray& tracks) 
{
  //  std::cout << "******** Detector " << GetName() << "BBshower  Coarse process = " << fCoarseProcessed << std::endl;
  
  // if(fCoarseProcessed)    return 0;

  // Call the parent's parent class coarse process to start filling out output variables
  //std::cout << "SBSGen  Coarse process " << std::endl;
  SBSGenericDetector::CoarseProcess(tracks);
  fCoarseProcessed = 1;
  return 0;
}

//_____________________________________________________________________________
Int_t SBSBBShower::FineProcess(TClonesArray& tracks)
{
  // Fine Shower processing.
  // Call parent's parent class to prepare any other variables
  SBSCalorimeter::FineProcess(tracks);


  // The parent class already sorted by energy, and the first cluster is the
  // one with the highest energy.
  // This function now needs to store the MCdata
  
  for (size_t i=0;i<fClusters.size();i++) {
     if(fDebug){
       cout << " cluster " << i << " E = " << fClusters[i]->GetE() << " " 
        << fClusters[i]->GetX() << " " << fClusters[i]->GetY() 
        << " " << fClusters[i]->GetMult()  << endl; 
    }

    if(fMCdata){
      fE_cl_res.push_back(1.0 - fClusters[i]->GetE());
      fX_cl_res.push_back(fClusters[i]->GetX());
      fY_cl_res.push_back(fClusters[i]->GetY());
    }
  }
  
  fFineProcessed = 1;
  return 0;

}



SBSBBShower::~SBSBBShower()
{
}

void SBSBBShower::LoadMCHitAt( Double_t x, Double_t y, Double_t E )
{  
  ClearEvent();
  SBSCalorimeterCluster *cluster = new SBSCalorimeterCluster(fNclublk);
  cluster->SetE(E);
  cluster->SetX(x);
  cluster->SetY(y);
  cluster->SetMult(0);
  fClusters.push_back(cluster);
  fNclus=fClusters.size();
}
 
 

void SBSBBShower::MakeCluster(Int_t nblk_size, SBSElement* blk) 
{
	     SBSCalorimeterCluster* cluster = new SBSCalorimeterCluster(nblk_size,blk);
	     fClusters.push_back(cluster);
}

void SBSBBShower::MakeCluster(Int_t nblk_size) 
{
	     SBSCalorimeterCluster* cluster = new SBSCalorimeterCluster(nblk_size);
	     fClusters.push_back(cluster);
}

void SBSBBShower::AddToCluster(Int_t nc,SBSElement* blk) 
{
  if (nc < (int)fClusters.size()) fClusters[nc]->AddElement(blk);
}

void SBSBBShower::MakeMainCluster() 
{
  if(!fClusters.empty()) {
    SBSCalorimeterCluster *clus = fClusters[0];
    fMainclus.e.push_back(clus->GetE());
    fMainclus.atime.push_back(clus->GetAtime());
    fMainclus.tdctime.push_back(clus->GetTDCtime());
    fMainclus.e_c.push_back(clus->GetE()*(fConst + fSlope*fAccCharge));
    fMainclus.x.push_back(clus->GetX());
    fMainclus.y.push_back(clus->GetY());
    fMainclus.n.push_back(clus->GetMult());
    fMainclus.blk_e.push_back(clus->GetEblk());
    fMainclus.blk_e_c.push_back(clus->GetEblk()*(fConst + fSlope*fAccCharge));
    fMainclus.id.push_back(clus->GetElemID());
    fMainclus.row.push_back(clus->GetRow());
    fMainclus.col.push_back(clus->GetCol());
  }
  //
  fNclus=0;
  for( const auto& cluster: fClusters ) {
    if( cluster->GetMult() > 0 ) fNclus++;
    }
  //
}

void SBSBBShower::SetSearchRegion(int rowmin, int rowmax, int colmin, int colmax)
{
  fSearchRowmin = rowmin;
  fSearchRowmax = rowmax;
  fSearchColmin = colmin;
  fSearchColmax = colmax;

  fSearchRegion = true;
  fMultClus = false;
}

void SBSBBShower::ClearEvent()
{
  SBSCalorimeter::ClearEvent();

  fEres = fXres = fYres = 0.0;
  fE_cl_res.clear();
  fX_cl_res.clear();
  fY_cl_res.clear();
}



//_____________________________________________________________
SBSElement* SBSBBShower::GetElement(UInt_t i)
{
  SBSElement* blk=nullptr;
  if(i < fElements.size()) blk = fElements[i];
  return blk;
}
