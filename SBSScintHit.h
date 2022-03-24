#ifndef ROOT_SBSScintHit
#define ROOT_SBSScintHit

/////////////////////////////////////////////////////////////////////////////
//                                                                         //
// SBSScintHit                                                             //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

#include "TObject.h"
#include "TRef.h"
#include "SBSScintBar.h"
#include <cstdio>

class SBSScintHit : public TObject {

 public:
 
    
  SBSScintHit(const SBSScintBar* bar,Int_t planenum, Int_t barnum, 
	      Double_t ypos, Double_t Tof, Double_t HitEnergy, Double_t Tdiff);   
  explicit SBSScintHit(const SBSScintHit* pScHit = nullptr);
  SBSScintHit(const SBSScintHit* pScHit, Int_t clusternum);
  SBSScintHit(const SBSScintHit* pScHit, Int_t planenum, Int_t barnum_nd);

  virtual ~SBSScintHit();

  SBSScintBar* GetScintBar() const { return (SBSScintBar*)fScBar.GetObject(); }
  Int_t GetPlaneNum() const {return fPlaneNum;}
  Int_t GetBarNum() const {return fBarNum;}
  Int_t GetBarNum_nd() const {return fBarNum_nd;}
  
  Double_t GetHitXPos() const {return fHitXPos;}
  Double_t GetHitYPos() const {return fHitYPos;}
  Double_t GetHitZPos() const {return fHitZPos;}
  
  Double_t GetHitTOF() const {return fHitTOF;}
  Double_t GetHitEdep() const {return fHitEdep;}
  Double_t GetHitTdiff() const {return fTdiff;}
  Int_t GetOrder() const {return fOrder;}
  Int_t GetClusterNum() const {return fClusterNum;}

  void SetScintBar(SBSScintBar* bar) {fScBar=bar;}
  void SetPlaneNum(Int_t planenum) {fPlaneNum=planenum;}
  void SetBarNum(Int_t barnum) {fBarNum=barnum;}
  void SetBarNum_nd(Int_t barnum_nd) {fBarNum_nd=barnum_nd;}
  void SetYHitPos(Double_t ypos) {fHitYPos = ypos;}
  void SetHitTOF(Double_t Tof) {fHitTOF = Tof;}
  void SetHitEdep(Double_t HitEnergy) {fHitEdep = HitEnergy;}
  void SetHitTdiff(Double_t Tdiff) {fTdiff = Tdiff;}
  void SetHitOrder(Int_t order) {fOrder = order;}
  void SetClusterNum(Int_t clusternum) {fClusterNum = clusternum;}

  void AddEnergy(Double_t HitEnergy) {fHitEdep += HitEnergy;}

  Int_t CopyScintHit(const SBSScintHit* pScHit);

  void Clear(Option_t *s="");
  
  // need to be able to sort the hits by bar number

  Bool_t IsSortable() const { return kTRUE; }
  Int_t  Compare(const TObject* obj) const;
  
 private:

  TRef fScBar;
  Int_t fPlaneNum;   // this or first plane number
  Int_t fBarNum;     // 
  Int_t fBarNum_nd;  // 
  Double_t fHitXPos; // 
  Double_t fHitYPos; // 
  Double_t fHitZPos; // 
  Double_t fHitTOF;  //   time of hit
  Double_t fTdiff;   //   left-to-right PMT time difference
  Double_t fHitEdep; //   total energy deposited
  Int_t fOrder;
  Int_t fClusterNum;
  
 public:

  ClassDef(SBSScintHit,2)  // Reconstructed hit-information for a complete hit
};
/////////////////////////////////////////////////////////////////
#endif
