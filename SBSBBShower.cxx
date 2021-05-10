///////////////////////////////////////////////////////////////////////////////
//
// SBSBBShower
//
///////////////////////////////////////////////////////////////////////////////
#include "SBSBBShower.h"
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
  SBSGenericDetector::SetDisableRefADC(kTRUE);
  // Call the parent class ReadDatabase first
  Int_t err = SBSCalorimeter::ReadDatabase(date);
  std::cout << " return from SBSCal " << std::endl;
  if(err) {
    return err;
  }

  FILE* file = OpenFile( date );
  if( !file ) return kFileError;

  std::vector<Double_t> dxyz;
  // Readout components needed by BBShower
  DBRequest config_request[] = {
    { "thr_adc",      &fThrADC,     kDouble,  0, true },
    { "clus_rad",     &fClusRadius, kFloat,   0, true },
    { "mc_data",      &fMCdata,     kInt,     0, true },// flag for MC data
    { "dxdydz",         &dxyz,         kDoubleV, 3 },  // dx and dy block spacings
    { 0 } ///< Request must end in a NULL
  };
  std::cout << " loading DB  " << fPrefix << std::endl;
  err = LoadDB( file, date, config_request, fPrefix );
  std::cout << " " << dxyz[0]<< " " << dxyz[1]<< " " << dxyz[2] << " " << std::endl;
  if(err) {
    return err;
  }
  fClusBlockRadX = Int_t(fClusRadius/dxyz[0]);
  fClusBlockRadY = Int_t(fClusRadius/dxyz[1]);

  if(fMaxNclus>1)fMultClus = true;

  return 0;
}


Int_t SBSBBShower::DefineVariables( EMode mode )
{
  if( mode == kDefine && fIsSetup ) return kOK;
  // Initialize parent variables first
  Int_t err = SBSCalorimeter::DefineVariables(mode);

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
Int_t SBSBBShower::FindGoodHit(SBSElement *blk) {
  Int_t GoodHit=0;
  if (blk->ADC() && blk->HasData() && fModeADC != SBSModeADC::kWaveform) {
    //std::cout << " findGOOHIt Row = " << blk->GetRow() << " " << " " << blk->GetCol() << " " << std::endl;
    // std::cout << " Num hits = " << blk->ADC()->GetNHits() << std::endl;
    Int_t bnhits = blk->ADC()->GetNHits();
    UInt_t GoodHitIndex = 999;
    Float_t CentTime = 300. ;
    Float_t WidthTime = 50. ;    
    for (Int_t ih=0;ih<bnhits;ih++) {
          const SBSData::PulseADCData &hit = blk->ADC()->GetHit(ih);
	  //	  std::cout << "ih = "<< ih << " "  << abs(hit.time.val- CentTime) << " " << hit.time.val << std::endl;
	  if ( abs(hit.time.val- CentTime) < WidthTime) GoodHitIndex=ih;
    }
    blk->ADC()->SetGoodHit(0);
    GoodHit=1;
  }
  //  if (blk->ADC()&& fModeADC == SBSModeADC::kWaveform) {
        SBSData::Waveform *wave = blk->Waveform();
	if (wave->GetTime().val > 0) GoodHit=1;
	//	std::cout << " findg = " << GoodHit << std::endl;
	// }
  return GoodHit;
}
//_____________________________________________________________________________
Int_t SBSBBShower::CoarseProcess(TClonesArray& tracks) 
{
  // Someone already called us for this event
  // std::cout << "BBshower  Coarse process = " << fCoarseProcessed << std::endl;
  
  // if(fCoarseProcessed)    return 0;

  // Call the parent's parent class coarse process to start filling out output variables
  //std::cout << "SBSGen  Coarse process " << std::endl;
  SBSGenericDetector::CoarseProcess(tracks);
  Int_t col, row;
  Double_t  energy_max = 0.0;

  std::set<Double_t> locenergy_max;
  std::map< Double_t, std::pair<Int_t, Int_t> > locrowcol_max;

# if not defined(_WIN32)//Win32 compiler do not support variable as array size
  Double_t energyDep[fNcolsMax][fNrows];
# else
  Double_t energyDep[100][100];
# endif

  Int_t Colblk = 0;
  Int_t Rowblk = 0;

  Double_t energyX = 0.0;
  Double_t energyY = 0.0;
  Double_t energyTotal = 0.0;
  Double_t energyClusterTotal = 0.0;
  Double_t X, Y;
  SBSCalorimeterCluster *cluster = new SBSCalorimeterCluster(fNclublk);

  // much simpler if a search region is defined:
  // check the blocks in the search region. 
  // If there is at least one above the defined threshold, build the cluster.
  // otherwise, return false
  if(fSearchRegion){
    if(fDebug)cout << " fSearchColmin "  << fSearchColmin 
      << " fSearchColmax "  << fSearchColmax 
        << " fSearchRowmin "  << fSearchRowmin 
        << " fSearchRowmax "  << fSearchRowmax 
        << endl;

    for(col = fSearchColmin; col <= fSearchColmax; col++){
      for(row = fSearchRowmin; row <= fSearchRowmax; row++){
        SBSElement *blk = fElementGrid[row][col][0];
        if(fDebug)cout << " col " << col << " row " << row
          << " blk->GetE()"<< blk->GetE() << endl;
        cluster->AddElement(blk); //< Add element to cluster
      }
    }

    if(cluster->GetMaxE()>=fEmin){
      fClusters.push_back(cluster);
      fNclus = fClusters.size();
    }

     fCoarseProcessed = 1;
    return 0; 
  }

  // Otherwise look over entire calorimeter
  for( row = 0; row < fNrows; row++ ){
    for( col = 0; col < fNcols[row]; col++ ){
      SBSElement *blk = fElementGrid[row][col][0];
      energyDep[col][row] = blk->GetE();

      //cout << " energy dep ? " << energyDep[col][row] << endl;

      if( energyDep[col][row] < 0.0 ) 
        energyDep[col][row] = 0.0;
      energyTotal += energyDep[col][row];

      //cluster seeds
      if(energyDep[col][row]>fEmin){
        if(energyDep[col][row]>energy_max){
          energy_max=energyDep[col][row];
          Colblk = col;
          Rowblk = row;
        }

        if(fDebug)
          cout << " col " << col << " row " << row << " block ? " << blk->GetID()
            << " Edep ? " << energyDep[col][row] << endl; 

        if(fMultClus){
          locenergy_max.insert(energyDep[col][row]);
          locrowcol_max[energyDep[col][row]] = std::make_pair(row, col);
        }
      }//end      
    }
    //      cout << endl;
  }

  if(energy_max < fEmin){
    fCoarseProcessed = 1;
    return 0;
  }

  if(fMultClus){
    for(std::set<Double_t>::iterator it = locenergy_max.begin(); it!=locenergy_max.end(); ++it){
      Double_t emax_i = *it;
      if(fDebug)
        cout << emax_i << " " << locrowcol_max[emax_i].first << " " << locrowcol_max[emax_i].second << endl;
      for(std::set<Double_t>::iterator jt = locenergy_max.begin(); jt!=it; ++jt){
        Double_t emax_j = *jt;

        if(abs(locrowcol_max[emax_i].first-locrowcol_max[emax_j].first)<=fClusBlockRadX &&
            abs(locrowcol_max[emax_i].second-locrowcol_max[emax_j].second)<=fClusBlockRadY){
          if(emax_i<emax_j){
            locenergy_max.erase(it);
          }else{
            locenergy_max.erase(jt);
          }
        }
      }
    }

    if(fDebug){
      cout << "after cleaning neighbors" << endl;//after cleaning
      for(std::set<Double_t>::iterator it = locenergy_max.begin(); it!=locenergy_max.end(); ++it){
        Double_t emax_i = *it;
        cout << emax_i << " " << locrowcol_max[emax_i].first << " " << locrowcol_max[emax_i].second << endl;
      }
    }

    while(int(locenergy_max.size())>fMaxNclus){
      std::set<Double_t>::iterator it = locenergy_max.begin();
      locenergy_max.erase(it);
    }

    if(fDebug){
      cout << "after cutting tails" << endl;//after cleaning
      for(std::set<Double_t>::iterator it = locenergy_max.begin(); it!=locenergy_max.end(); ++it){
        Double_t emax_i = *it;
        cout << emax_i << " " << locrowcol_max[emax_i].first << " " << locrowcol_max[emax_i].second << endl;
      }

      cout << " col max ? " << Colblk << " row max ? " << Rowblk << " Emax ? " << energy_max << " Emin = " << fEmin << endl;
      cout << " secondary clusters seeds size " << locrowcol_max.size() << " " << locenergy_max.size() << endl;
    }
  }else{//end if(fMultClus)
    locenergy_max.insert(energy_max);
    locrowcol_max[energy_max] = std::make_pair(Rowblk, Colblk);
  }

  //  Double_t energyClusterGreatest = 0.0;

  Int_t i, j;//, k=0;
  Int_t mnrow, mxrow, mncol, mxcol;

  //for(size_t cls = 0; cls<locenergy_max.size(); cls++){
  for(std::set<Double_t>::iterator it = locenergy_max.begin(); it!=locenergy_max.end(); ++it){
    Double_t emax_i = *it;
    energyClusterTotal = 0.0;

    /*
       mnrow=TMath::Max(locrowmax[cls]-fClusBlockRadX,0);
       mxrow=TMath::Min(locrowmax[cls]+fClusBlockRadX,fNrows-1);
       mncol=TMath::Max(loccolmax[cls]-fClusBlockRadY,0);
       mxcol=TMath::Min(loccolmax[cls]+fClusBlockRadY,fNcols-1);
       */

    mnrow=TMath::Max(locrowcol_max[emax_i].first-fClusBlockRadX,0);
    mxrow=TMath::Min(locrowcol_max[emax_i].first+fClusBlockRadX,fNrows-1);
    mncol=TMath::Max(locrowcol_max[emax_i].second-fClusBlockRadY,0);
    mxcol=TMath::Min(locrowcol_max[emax_i].second+fClusBlockRadY,fNcolsMax-1);

    if(fDebug){
      cout << " Cluster: mnrow " << mnrow << " mxrow " << mxrow 
        << " mncol " << mncol << " mxcol " << mxcol << endl;
    }
    for( i = mnrow; i <= mxrow; i++ ){
      for( j = mncol; j <= mxcol; j++){
        energyClusterTotal += energyDep[j][i];
        //fEblk[k] = energyDep[j][i];
        //k++;
      }
    }

    //Double_t energyCluster = energyClusterTotal;
    //if( energyCluster < 0.0 ) return 0;

    X = fElementGrid[Rowblk][Colblk][0]->GetX();
    Y = fElementGrid[Rowblk][Colblk][0]->GetY();

    if(fDebug)cout << "Got a cluster! E = " << energyClusterTotal << " X = " << X << " Y = " << Y << endl;


    Int_t  blockcounter = 0;
    for( i = mnrow; i <= mxrow; i++ ){
      for( j = mncol; j <= mxcol; j++ ){
        if( (i >= 0 && i < fNrows ) && ( j >=0 && j < fNcols[i] ) ){
          SBSElement *blk = fElementGrid[i][j][0];
          energyX += energyDep[j][i]*blk->GetX();
          energyY += energyDep[j][i]*blk->GetY();

          /*
             if(i!=fBlocks[BlockColRowToNumber(j,i)]->GetRow()){
             cout << "row " << i << " " << fBlocks[BlockColRowToNumber(j,i)]->GetRow() << endl;
             }
             if(j!=fBlocks[BlockColRowToNumber(j,i)]->GetCol()){
             cout << "col " << j << " " << fBlocks[BlockColRowToNumber(j,i)]->GetCol() << endl;
             }
             */
          cluster->AddElement( blk );
          if(fDebug){
            cout << "Cluster " << &cluster << " Adding block row " 
              << blk->GetRow() 
              << " col " << blk->GetCol() 
              << " E " << blk->GetE() 
              << " new size " << cluster->GetSize() << " block " << blockcounter 
              << " address " << cluster->GetElement(blockcounter) << endl;
          } 
          blockcounter++;
        }
      }
    }

    X = energyX/energyClusterTotal;
    Y = energyY/energyClusterTotal;

    if(fDebug){
      cout << energyClusterTotal << " " << X << " " << Y 
        << " " << fOrigin.X() << " " << fOrigin.Y() << " " << cluster->GetSize() << endl;
    }

    cluster->SetE( energyClusterTotal );
    cluster->SetX( X+fOrigin.X() );
    cluster->SetY( Y+fOrigin.Y() );
    //cluster.SetX( X );
    //cluster.SetY( Y );

    if(fMCdata){
      fEres = cluster->GetE();//-
      fXres = cluster->GetX();//-
      fYres = cluster->GetY();//-
    }

    fClusters.push_back(cluster);
    fNclus = fClusters.size();
    if(fDebug)
      cout << "Added - we now have " << fNclus << endl;
  }

  // locrowmax.clear();  
  // loccolmax.clear();  
  locenergy_max.clear();  
  locrowcol_max.clear();  
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
    SBSCalorimeterCluster *cluster = fClusters[i];
    if(fDebug){
      cout << fClusters[i] << " " << cluster->GetE() << " " 
        << cluster->GetX() << " " << cluster->GetY() 
        << " " << cluster->GetMult()  << endl; 
    }

    if(fMCdata){
      fE_cl_res.push_back(1.0 - cluster->GetE());
      fX_cl_res.push_back(cluster->GetX());
      fY_cl_res.push_back(cluster->GetY());
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

