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
    SBSTimingHodoscopeCluster(Int_t nmaxblk);
    //SBSTimingHodoscopeCluster(Int_t nmaxblk, SBSTimingHodoscopeBar* bar);
    virtual ~SBSTimingHodoscopeCluster();

    Float_t GetXmean() const {return fXmean;}
    Float_t GetYmean() const {return fYmean;}
    Float_t GetTmean() const {return fTmean;}

    Int_t GetNMaxElements() const {return fNMaxElements;}

    void SetXmean(Float_t var) {fXmean=var;}
    void SetYmean(Float_t var) {fYmean=var;}
    void SetTmean(Float_t var) {fTmean=var;}

    SBSTimingHodoscopeBar* GetElement(UInt_t i);
    std::vector<SBSTimingHodoscopeBar*>& GetElements() {return fElements;}

    Int_t GetSize() {return fMult;}

    void AddElement(SBSTimingHodoscopeBar* bar);

    void ClearEvent();

private:

    Float_t fXmean;       // x position of the center
    Float_t fYmean;       // y position of the center
    Float_t fTmean;       // Energy deposit in block

    Int_t fMult;      // Number of bars in the cluster
    Int_t fNMaxElements;// Max number of blocks
    //SBSElement* fMaxElement;  // Element with maximum in the cluster

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
