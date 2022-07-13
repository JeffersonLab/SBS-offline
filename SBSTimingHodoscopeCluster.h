#ifndef ROOT_SBSTimingHodoscopeCluster
#define ROOT_SBSTimingHodoscopeCluster

////////////////////////////////////////////////////////////////////////////
//                                                                        //
// SBSTimingHodoscopeCluster                                              //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include "TObject.h"
#include "SBSTimingHodoscopeBar.h"

class SBSTimingHodoscopeCluster : public TObject {

public:

    SBSTimingHodoscopeCluster();
    explicit SBSTimingHodoscopeCluster(Int_t nmaxblk);
    SBSTimingHodoscopeCluster(Int_t nmaxblk, SBSTimingHodoscopeBar* bar);
    virtual ~SBSTimingHodoscopeCluster();

    Double_t GetXmean() const {return fXmean;}
    Double_t GetYmean() const {return fYmean;}
    Double_t GetTmean() const {return fTmean;}
    Double_t GetToTmean() const {return fToTmean;}
    Double_t GetTdiff() const {return fMaxElement->GetTimeDiff(); }
    
    Int_t GetNMaxElements() const {return fNMaxElements;}
    Int_t GetMaxBarID() const { return fMaxElement->GetBarNum(); }
    
    void SetXmean(Double_t var) {fXmean=var;}
    void SetYmean(Double_t var) {fYmean=var;}
    void SetTmean(Double_t var) {fTmean=var;}
    void SetToTmean(Double_t var) {fToTmean=var;}
    
    SBSTimingHodoscopeBar* GetElement(UInt_t i);
    std::vector<SBSTimingHodoscopeBar*>& GetElements() {return fElements;}

    Int_t GetSize() const {return fMult;}

    Bool_t AddElement(SBSTimingHodoscopeBar* bar);

    virtual void Clear( Option_t* opt="" );

private:

    Double_t fXmean;       // x position of the center
    Double_t fYmean;       // y position of the center
    Double_t fTmean;       // Energy deposit in block
    Double_t fToTmean;       // Energy deposit in block
    
    Int_t fMult;      // Number of bars in the cluster
    Int_t fNMaxElements;// Max number of blocks
    SBSTimingHodoscopeBar* fMaxElement;  // Element with maximum ToT in the cluster

    std::vector<SBSTimingHodoscopeBar*> fElements; //

    ClassDef(SBSTimingHodoscopeCluster,1)   // Generic shower cluster class
};

struct SBSTimingHodoscopeClusterCompare {
  bool operator() (const SBSTimingHodoscopeCluster *l,
      const SBSTimingHodoscopeCluster *r) {
    if(!l || !r)
      return false;
    return l->GetTmean() < r->GetTmean();
  }
};


#endif
