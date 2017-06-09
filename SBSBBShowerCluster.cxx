///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// THaBBShowerCluster                                                        //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include "THaBBShowerCluster.h"
#include <iostream>

using namespace std;

const Float_t THaBBShowerCluster::kBig =(Float_t)1e15;

//_____________________________________________________________
THaBBShowerCluster::THaBBShowerCluster(Int_t nmaxblk, THaBBShowerBlock* block) 
: fNMaxBlocks(nmaxblk)
{
    fBlocks = new THaBBShowerBlock*[fNMaxBlocks];
    fBlocks[0] = block;
    fX = block->GetX();
    fY = block->GetY();
    fE = block->GetE();
    fMult = 1;
}


//_____________________________________________________________
THaBBShowerCluster::THaBBShowerCluster(Int_t nmaxblk) 
: fNMaxBlocks(nmaxblk)
{
    fBlocks = new THaBBShowerBlock*[fNMaxBlocks];
    fX = 0;
    fY = 0;
    fE = 0;
    fMult = 0;
}


//_____________________________________________________________
THaBBShowerCluster::THaBBShowerCluster() {
    fX = kBig;
    fY = kBig;
    fE = kBig;
    fMult = 0;
}

//_____________________________________________________________
THaBBShowerCluster::~THaBBShowerCluster() {

    DeleteArrays();

}

//_____________________________________________________________
void THaBBShowerCluster::AddBlock(THaBBShowerBlock* block) {

    if (fMult<fNMaxBlocks) {
        fBlocks[fMult] = block;
        block->SetStat(1);
        fMult++;
        //     cout << fX << " " << fE << " " << fMult << " " 
        // 	 << block->GetX() << " " << block->GetE() << endl;
        //     cout << fX*fE*(fMult-1) << " " << block->GetX()*block->GetE() << endl;
        fX = (fX*fE*(fMult-1) + block->GetX()*block->GetE()) / 
            (fE+block->GetE());
        fY = (fY*fE*(fMult-1) + block->GetY()*block->GetE()) / 
            (fE+block->GetE());
        fE += block->GetE();
    }

}

//_____________________________________________________________
void THaBBShowerCluster::ClearEvent() {
    fMult=0;fX=fY=fE=0.;
    DeleteArrays();
    fBlocks = new THaBBShowerBlock*[fNMaxBlocks];
}

//_____________________________________________________________
void THaBBShowerCluster::DeleteArrays() {
    delete [] fBlocks; fBlocks = 0;
}

ClassImp(THaBBShowerCluster)

