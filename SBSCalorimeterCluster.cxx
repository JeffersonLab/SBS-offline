///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// SBSCalorimeterCluster                                                     //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include "SBSCalorimeterCluster.h"
#include <iostream>
#include <DataType.h>

ClassImp(SBSCalorimeterCluster)   // Generic shower cluster class

using namespace std;


//_____________________________________________________________
SBSCalorimeterCluster::SBSCalorimeterCluster(Int_t nmaxblk, SBSElement* block) 
  : fNMaxElements(nmaxblk), fMaxElement(0)
{
  fElements.clear();
  fElements.reserve(fNMaxElements);
  fElements.push_back(block);
  fX = block->GetX();
  fY = block->GetY();
  fE = block->GetE();
  fEblk = block->GetE();
  fMult = 1;
  fAgain = block->GetAgain();
  fAtime  = block->GetAtime();
  fAtimeMean = fAtime;
  
  fTDCtime  = -1000.0; //Assign a default of -1000 ns for TDC time unless we have actual data!
  fTDCtimeTW  = -1000.0; //Assign a default of -1000 ns for TDC time unless we have actual data!
    
  fRow  = block->GetRow();
  fCol  = block->GetCol();
  fElemID  = block->GetID();

  fMaxElement = block; 
  
  fNgoodTDChits = 0;
  fTDCtimeMean = 0.0;
  fTDCtimeMeanTW = 0.0;
  
  fEblk_GoodTDC = 0.0;
  fE_GoodTDC = 0.0;

  fRowGoodTDC = -1;
  fColGoodTDC = -1;
  fElemIDGoodTDC = -1;
  
  if( block->HasTDCData() ){
    fNgoodTDChits = 1;
    fTDCtime = block->GetTDCtime();
    fTDCtimeTW = block->GetTDCtimeTW();
    fTDCtimeMean = block->GetTDCtime();
    fTDCtimeMeanTW = block->GetTDCtimeTW();
    fEblk_GoodTDC = block->GetE();
    fE_GoodTDC = block->GetE();
    fRowGoodTDC = block->GetRow();
    fColGoodTDC = block->GetCol();
    fElemIDGoodTDC = block->GetID();
  }
  
}


//_____________________________________________________________
SBSCalorimeterCluster::SBSCalorimeterCluster(Int_t nmaxblk) 
  : fNMaxElements(nmaxblk)
{
  fElements.clear();
  fElements.reserve(fNMaxElements);
  fX = 0;
  fY = 0;
  fE = 0;
  fEblk = 0;
  fMult = 0;
  fAgain = 0.;
  fAtime = 0;
  fTDCtime = 0;
  fTDCtimeTW = 0;
  fRow  = -1;
  fCol  = -1;
  fElemID  = -1;
  fMaxElement = 0;
  
  fNgoodTDChits = 0;
  fTDCtimeMean = 0.0;
  fTDCtimeMeanTW = 0.0;
  fAtimeMean = 0.0;
  fEblk_GoodTDC = 0.0;
  fE_GoodTDC = 0.0;

  fRowGoodTDC = -1;
  fColGoodTDC = -1;
  fElemIDGoodTDC = -1;
}


//_____________________________________________________________
SBSCalorimeterCluster::SBSCalorimeterCluster() {
  fX = kBig;
  fY = kBig;
  fE = kBig;
  fEblk = 0;
  fMult = 0;
  fAgain = 0.;
  fAtime = 0;
  fTDCtime = 0;
  fTDCtimeTW = 0;
  fRow  = -1;
  fCol  = -1;
  fElemID  = -1;
  fMaxElement = 0;

  fElements.clear();
  
  fNgoodTDChits = 0;
  fTDCtimeMean = 0.0;
  fTDCtimeMeanTW = 0.0;
  fAtimeMean = 0.0;
  fEblk_GoodTDC = 0.0;
  fE_GoodTDC = 0.0;
  
  fRowGoodTDC = -1;
  fColGoodTDC = -1;
  fElemIDGoodTDC = -1;
}

//_____________________________________________________________
SBSCalorimeterCluster::~SBSCalorimeterCluster()
{ 
}

//_____________________________________________________________
void SBSCalorimeterCluster::AddElement(SBSElement* block) {
  if (fMult<fNMaxElements) {
    fElements.push_back(block);
    fMult = fElements.size();
    block->SetStat(1);
    fX = (fX*fE + block->GetX()*block->GetE()) / (fE+block->GetE());
    fY = (fY*fE + block->GetY()*block->GetE()) / (fE+block->GetE());
    fAtimeMean = (fAtimeMean*fE + block->GetAtime()*block->GetE())/(fE+block->GetE());
	
    fE += block->GetE();
    if(block->GetE() > fEblk) {
      fEblk = block->GetE();
      fAtime = block->GetAtime();
      fAgain = block->GetAgain();
      fRow = block->GetRow();
      fCol = block->GetCol();
      fElemID = block->GetID();
    }

    if( block->HasTDCData() ){
      fTDCtimeMean = (fTDCtimeMean*fE_GoodTDC + block->GetTDCtime() * block->GetE())/(fE_GoodTDC + block->GetE());
      fTDCtimeMeanTW = (fTDCtimeMeanTW*fE_GoodTDC + block->GetTDCtimeTW() * block->GetE())/(fE_GoodTDC + block->GetE());
      fE_GoodTDC += block->GetE();
      if( block->GetE() > fEblk_GoodTDC ){
	fTDCtime = block->GetTDCtime();
	fTDCtimeTW = block->GetTDCtimeTW();
	fEblk_GoodTDC = block->GetE();
	fRowGoodTDC = block->GetRow();
	fColGoodTDC = block->GetCol();
	fElemIDGoodTDC = block->GetID();
      }
      
      fNgoodTDChits++;
    }
    
    // Keep a pointer to the element with the highest energy
    if(!fMaxElement) {
      fMaxElement = block;
    } else if ( block->GetE() > fMaxElement->GetE() ) {
      fMaxElement = block;
    }
  }

}

//_____________________________________________________________
void SBSCalorimeterCluster::Clear( Option_t* opt ) {
  fMult=0;fX=fY=fE=0.;
  fEblk=0;
  fAgain=0.;
  fAtime=0;
  fTDCtime=0;
  fTDCtimeTW=0;
  fRow=fCol=fElemID-1;
  fMaxElement = 0;
  fElements.clear();
  fNgoodTDChits = 0;
  fAtimeMean = 0.0;
  fTDCtimeMean = 0.0;
  fTDCtimeMeanTW = 0.0;
  fEblk_GoodTDC = 0.0;
  fE_GoodTDC = 0.0;
  fRowGoodTDC=fColGoodTDC=fElemIDGoodTDC=-1;
}

//_____________________________________________________________
SBSElement* SBSCalorimeterCluster::GetElement(UInt_t i)
{
  SBSElement* blk=0;
  if(i < fElements.size()) blk = fElements[i];
  return blk;
}



ClassImp(SBSCalorimeterCluster)
