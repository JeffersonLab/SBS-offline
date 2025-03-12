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
  explicit SBSCalorimeterCluster(Int_t nmaxblk);
  SBSCalorimeterCluster(Int_t nmaxblk, SBSElement* block);
  virtual ~SBSCalorimeterCluster();

  Double_t GetX() const {return fX;}
  Double_t GetY() const {return fY;}
  Double_t GetE() const {return fE;}
  Double_t GetAgain() const {return fAgain;}
  Double_t GetAtime() const {return fAtime;}
  Double_t GetTDCtime() const {return fTDCtime;}
  Double_t GetTDCtimeTW() const {return fTDCtimeTW;}
  Double_t GetEblk() const {return fEblk;}
  Int_t   GetMult() const {return fMult;}
  Int_t   GetRow()  const {return fRow; }
  Int_t   GetCol()  const {return fCol; }
  Int_t   GetElemID()  const {return fElemID; }
  SBSElement* GetMaxElement() { return fMaxElement; }
  Double_t GetMaxE() const {if(fMaxElement) return fMaxElement->GetE(); return 0; }

  Int_t GetNMaxElements() const {return fNMaxElements;}

  //New Get methods for Blocks with good TDC time:
  Double_t GetAtimeMean() const { return fAtimeMean; }
  Double_t GetE_GoodTDC() const { return fE_GoodTDC; }
  Double_t GetTDCtimeMean() const { return fTDCtimeMean; }
  Double_t GetTDCtimeMeanTW() const { return fTDCtimeMeanTW; }
  Double_t GetEblk_GoodTDC() const { return fEblk_GoodTDC; }
  Int_t GetNgoodTDChits() const { return fNgoodTDChits; }
  Int_t GetRowGoodTDC() const { return fRowGoodTDC; }
  Int_t GetColGoodTDC() const { return fColGoodTDC; }
  Int_t GetElemIDGoodTDC() const { return fElemIDGoodTDC; }
  //We could also write "Set" methods for these new variables, but at least for now, we
  //shouldn't need to. 
  
  void SetX(Double_t var) {fX=var;}
  void SetY(Double_t var) {fY=var;}
  void SetE(Double_t var) {fE=var;}
  void SetAgain(Double_t var) {fAgain=var;}
  void SetAtime(Double_t var) {fAtime=var;}
  void SetTDCtime(Double_t var) {fTDCtime=var;}
  void SetEblk(Double_t var) {fEblk=var;}
  void SetMult(Int_t var) {fMult=var;}
  void SetRow(Int_t var) { fRow=var; }
  void SetCol(Int_t var) { fCol=var; }
  void SetElemID(Int_t var) { fElemID=var; }

  SBSElement* GetElement(UInt_t i);
  std::vector<SBSElement*>& GetElements() {return fElements;}

  
  Int_t GetSize() const {return fMult;}

  void AddElement(SBSElement* block);

  virtual void Clear( Option_t* opt="" );

private:
  
  Double_t fX;       // x position of the center
  Double_t fY;       // y position of the center
  Double_t fE;       // sum of all Energy deposits in blocks
  Double_t fE_GoodTDC; //sum of energy deposits in blocks with good TDC hits
    
  Double_t fAgain;   // ADC gain coefficient (GeV/pC)
  Double_t fAtime;       // ADC time  of block with highest E
  Double_t fTDCtime;       // TDC time  of block with highest E
  Double_t fTDCtimeTW;    // Walk corrected TDC time of block with highest E
  Double_t fEblk;    // Energy of block with highest E
  Double_t fAtimeMean; //Energy-weighted mean ADC time of all blocks
  Double_t fTDCtimeMean; //Energy-weighted mean TDC time of all blocks in the cluster with good TDC hits
  Double_t fTDCtimeMeanTW;    // Energy-weighted mean TDC time of all Walk corrected blocks in the cluster with good TDC hits
  Double_t fEblk_GoodTDC; //Energy of block with highest E subject to the requirement of a good TDC hit;
  Int_t   fRow;     // Row of block with highest E
  Int_t   fCol;     // Row of block with highest E
  Int_t   fElemID;     // ElemID of block with highest E

  //Store row, col, and element ID for highest-energy block with good TDC hit:
  Int_t fRowGoodTDC; 
  Int_t fColGoodTDC;
  Int_t fElemIDGoodTDC; 

  Int_t fNgoodTDChits;
  
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
