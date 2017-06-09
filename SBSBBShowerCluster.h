#ifndef ROOT_THaBBShowerCluster
#define ROOT_THaBBShowerCluster

////////////////////////////////////////////////////////////////////////////
//                                                                        //
// THaBBShowerCluster                                                      //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include "TObject.h"
#include "THaBBShowerBlock.h"

class THaBBShowerCluster : public TObject {

public:

    THaBBShowerCluster();
    THaBBShowerCluster(Int_t nmaxblk);
    THaBBShowerCluster(Int_t nmaxblk, THaBBShowerBlock* block);
    virtual ~THaBBShowerCluster();

    Float_t GetX() const {return fX;}
    Float_t GetY() const {return fY;}
    Float_t GetE() const {return fE;}
    Int_t   GetMult() const {return fMult;}

    Int_t GetNMaxBlocks() const {return fNMaxBlocks;}

    void SetX(Float_t var) {fX=var;}
    void SetY(Float_t var) {fY=var;}
    void SetE(Float_t var) {fE=var;}
    void SetMult(Int_t var) {fMult=var;}

    THaBBShowerBlock** GetBlocks() {return fBlocks;}

    Int_t GetSize() {return fMult;}

    void AddBlock(THaBBShowerBlock* block);

    void ClearEvent();
    void DeleteArrays();

private:

    static const Float_t kBig; // = 1e15;

    Float_t fX;       // x position of the center
    Float_t fY;       // y position of the center

    Float_t fE;       // Energy deposit in block

    Int_t fMult;      // Number of blocks in the cluster
    Int_t fNMaxBlocks;// Max number of blocks


    THaBBShowerBlock** fBlocks; //[fNMaxBlocks] List of blocks in cluster

    ClassDef(THaBBShowerCluster,0)   // Generic shower cluster class
};

#endif
