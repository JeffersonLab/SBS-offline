#ifndef ROOT_SBSCalorimeterCluster
#define ROOT_SBSCalorimeterCluster

////////////////////////////////////////////////////////////////////////////
//                                                                        //
// SBSCalorimeterCluster                                                  //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include "TObject.h"
#include "SBSElement.h"

class SBSCalorimeterCluster : public TObject {

public:

    SBSCalorimeterCluster();
    SBSCalorimeterCluster(Int_t nmaxblk);
    SBSCalorimeterCluster(Int_t nmaxblk, SBSElement* block);
    virtual ~SBSCalorimeterCluster();

    Float_t GetX() const {return fX;}
    Float_t GetY() const {return fY;}
    Float_t GetE() const {return fE;}
    Float_t GetEblk() const {return fEblk;}
    Int_t   GetMult() const {return fMult;}
    Int_t   GetRow()  const {return fRow; }
    Int_t   GetCol()  const {return fCol; }
    Int_t   GetElemID()  const {return fElemID; }
    SBSElement* GetMaxElement() { return fMaxElement; }
    Float_t GetMaxE() const {if(fMaxElement) return fMaxElement->GetE(); return 0; }

    Int_t GetNMaxElements() const {return fNMaxElements;}

    void SetX(Float_t var) {fX=var;}
    void SetY(Float_t var) {fY=var;}
    void SetE(Float_t var) {fE=var;}
    void SetEblk(Float_t var) {fEblk=var;}
    void SetMult(Int_t var) {fMult=var;}
    void SetRow(Int_t var) { fRow=var; }
    void SetCol(Int_t var) { fCol=var; }
    void SetElemID(Int_t var) { fElemID=var; }

    SBSElement* GetElement(UInt_t i);
    std::vector<SBSElement*>& GetElements() {return fElements;}

    Int_t GetSize() {return fMult;}

    void AddElement(SBSElement* block);

    void ClearEvent();

private:

    Float_t fX;       // x position of the center
    Float_t fY;       // y position of the center

    Float_t fE;       // Energy deposit in block
    Float_t fEblk;    // Energy of block with highest E
    Int_t   fRow;     // Row of block with highest E
    Int_t   fCol;     // Row of block with highest E
    Int_t   fElemID;     // ElemID of block with highest E

    Int_t fMult;      // Number of blocks in the cluster
    Int_t fNMaxElements;// Max number of blocks
    SBSElement* fMaxElement;  // Element with maximum in the cluster

    std::vector<SBSElement*> fElements; //

    ClassDef(SBSCalorimeterCluster,1)   // Generic shower cluster class
};

struct SBSCalorimeterClusterCompare {
  bool operator() (const SBSCalorimeterCluster *l,
      const SBSCalorimeterCluster *r) {
    if(!l || !r)
      return false;
    return l->GetE() < r->GetE();
  }
};


#endif
