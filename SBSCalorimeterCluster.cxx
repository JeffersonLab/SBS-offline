///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// SBSCalorimeterCluster                                                     //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include "SBSCalorimeterCluster.h"
#include <iostream>

ClassImp(SBSCalorimeterCluster)   // Generic shower cluster class

using namespace std;

const Float_t SBSCalorimeterCluster::kBig =(Float_t)1e15;

//_____________________________________________________________
SBSCalorimeterCluster::SBSCalorimeterCluster(Int_t nmaxblk, SBSCalorimeterBlock* block) 
: fNMaxBlocks(nmaxblk)
{
    fBlocks = new SBSCalorimeterBlock*[fNMaxBlocks];
    fBlocks[0] = block;
    fX = block->GetX();
    fY = block->GetY();
    fE = block->GetE();
    fEblk = block->GetE();
    fMult = 1;
    fRow  = block->GetRow();
    fCol  = block->GetCol();
}


//_____________________________________________________________
SBSCalorimeterCluster::SBSCalorimeterCluster(Int_t nmaxblk) 
: fNMaxBlocks(nmaxblk)
{
    fBlocks = new SBSCalorimeterBlock*[fNMaxBlocks];
    fX = 0;
    fY = 0;
    fE = 0;
    fEblk = 0;
    fMult = 0;
    fRow  = -1;
    fCol  = -1;
}


//_____________________________________________________________
SBSCalorimeterCluster::SBSCalorimeterCluster() {
    fX = kBig;
    fY = kBig;
    fE = kBig;
    fEblk = 0;
    fMult = 0;
    fRow  = -1;
    fCol  = -1;
}

//_____________________________________________________________
SBSCalorimeterCluster::~SBSCalorimeterCluster() {

    DeleteArrays();

}

//_____________________________________________________________
void SBSCalorimeterCluster::AddBlock(SBSCalorimeterBlock* block) {

    if (fMult<fNMaxBlocks) {
        fBlocks[fMult] = block;
        block->SetStat(1);
        fMult++;
        //     cout << fX << " " << fE << " " << fMult << " " 
        // 	 << block->GetX() << " " << block->GetE() << endl;
        //     cout << fX*fE*(fMult-1) << " " << block->GetX()*block->GetE() << endl;

        // Why is this multiplied by fMult-1 in the numerator, but
        // no fMult-1 in the denominator? This will inherently break
        // things as X and Y will increasingly become larger and larger.
        //fX = (fX*fE*(fMult-1) + block->GetX()*block->GetE()) / 
        //    (fE+block->GetE());
        //fY = (fY*fE*(fMult-1) + block->GetY()*block->GetE()) / 
        //    (fE+block->GetE());
        fX = (fX*fE + block->GetX()*block->GetE()) / (fE+block->GetE());
        fY = (fY*fE + block->GetY()*block->GetE()) / (fE+block->GetE());
        fE += block->GetE();
        if(block->GetE() > fEblk) {
          fEblk = block->GetE();
          fRow = block->GetRow();
          fCol = block->GetCol();
        }
    }

}

//_____________________________________________________________
void SBSCalorimeterCluster::ClearEvent() {
    fMult=0;fX=fY=fE=0.;
    fEblk=0;
    fRow=fCol=-1;
    DeleteArrays();
    fBlocks = new SBSCalorimeterBlock*[fNMaxBlocks];
}

//_____________________________________________________________
void SBSCalorimeterCluster::DeleteArrays() {
    delete [] fBlocks; fBlocks = 0;
}

ClassImp(SBSCalorimeterCluster)


