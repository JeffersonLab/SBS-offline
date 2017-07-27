#ifndef ROOT_SBSECalCluster
#define ROOT_SBSECalCluster

////////////////////////////////////////////////////////////////////////////
//                                                                        //
// SBSECalCluster                                                      //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include "TObject.h"
#include "SBSShowerBlock.h"

class SBSECalCluster : public TObject {

public:

    SBSECalCluster();
    SBSECalCluster(Int_t nmaxblk);
    SBSECalCluster(Int_t nmaxblk, SBSShowerBlock* block);
    virtual ~SBSECalCluster();

    Float_t GetX() const {return fX;}
    Float_t GetY() const {return fY;}
    Float_t GetE() const {return fE;}
    Int_t   GetMult() const {return fMult;}

    Int_t GetNMaxBlocks() const {return fNMaxBlocks;}

    void SetX(Float_t var) {fX=var;}
    void SetY(Float_t var) {fY=var;}
    void SetE(Float_t var) {fE=var;}
    void SetMult(Int_t var) {fMult=var;}

    SBSShowerBlock** GetBlocks() {return fBlocks;}

    Int_t GetSize() {return fMult;}

    void AddBlock(SBSShowerBlock* block);

    void ClearEvent();
    void DeleteArrays();

private:

    static const Float_t kBig; // = 1e15;

    Float_t fX;       // x position of the center
    Float_t fY;       // y position of the center

    Float_t fE;       // Energy deposit in block

    Int_t fMult;      // Number of blocks in the cluster
    Int_t fNMaxBlocks;// Max number of blocks


    SBSShowerBlock** fBlocks; //[fNMaxBlocks] List of blocks in cluster

    ClassDef(SBSECalCluster,0)   // Generic shower cluster class
};

#endif
