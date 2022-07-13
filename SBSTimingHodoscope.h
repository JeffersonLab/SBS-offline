//////////////////////////////////////////////////////////////////////////
//
// SBSTimingHodoscope
//
// General detector with TDC information
//
//////////////////////////////////////////////////////////////////////////

#ifndef SBSTimingHodoscope_h
#define SBSTimingHodoscope_h

#include "TClonesArray.h"
#include "THaTrack.h"

#include "SBSGenericDetector.h"
#include "SBSTimingHodoscopePMT.h"
#include "SBSTimingHodoscopeBar.h"
#include "SBSTimingHodoscopeCluster.h"

struct SBSTimingHodoscopeOutput {
  std::vector<Int_t>   n;   // Number of elements
  std::vector<Int_t>   id;   // ID
  std::vector<Double_t> x;   //< []
  std::vector<Double_t> y;   //< []
  std::vector<Double_t> t;   //< []
  std::vector<Double_t> tot;   //< []
  std::vector<Double_t> tdiff;   //< []
  std::vector<Int_t> trackindex; 
};

class SBSTimingHodoscope : public SBSGenericDetector {
public:
  SBSTimingHodoscope( const char* name, const char* description, 
      THaApparatus* apparatus=0);

  virtual ~SBSTimingHodoscope();
  
  // Overwritten derived functions
  virtual Int_t   ReadDatabase( const TDatime& date );
  virtual Int_t   DefineVariables( EMode mode = kDefine );
  virtual Int_t   CoarseProcess( TClonesArray& tracks );
  virtual Int_t   FineProcess( TClonesArray& tracks );
  /* virtual Int_t   FindGoodHit(SBSElement *element); */
  virtual void    Clear( Option_t* opt="" );
  // new functions
  Int_t   ConstructHodoscope();
  Double_t TimeWalk(Double_t time, Double_t tot, Double_t timewalk0, Double_t timewalk1);
  
  Int_t GetGoodBarsSize() const {return fGoodBarIDsTDC.size();};
  std::vector<Int_t> GetGoodBarsIDs() const {return fGoodBarIDsTDC;};
  std::vector<Double_t> GetGoodBarsMeanTimes() const {return fGoodBarTDCmean;};
  // time diff pos is hit pos along the bar length from time diff
  std::vector<Double_t> GetGoodBarsTimeDiffPos() const {return fGoodBarTDCpos;};
  std::vector<Double_t> GetGoodBarsVPos() const {return fGoodBarTDCvpos;};
  Double_t GetGoodBarVPosElement(Int_t i) const {
    if((Double_t)i<fGoodBarTDCvpos.size())
      return fGoodBarTDCvpos[i];
    else return -999.0;};
  Double_t GetGoodBarHPosElement(Int_t i) const {
    if((Double_t)i<fGoodBarTDCpos.size())
      return fGoodBarTDCpos[i];
    else return -999.0;};
  Double_t GetGoodBarMeanTimeElement(Int_t i) const {
    if((Double_t)i<fGoodBarTDCmean.size())
      return fGoodBarTDCmean[i];
    else return -999.0;};
  
  Int_t GetNClusters() const {return fClusters.size();};
  SBSTimingHodoscopeCluster* GetCluster(int i);
  
  Int_t GetID();
  Int_t GetSize();
  Double_t GetX();
  Double_t GetY();
  Double_t GetT();
  Double_t GetToT();

  Double_t GetVVal(std::vector<Double_t> &v, UInt_t i = 0 );
  Int_t GetVVal(std::vector<Int_t> &v, UInt_t i = 0 );
  
  void SetDataOutputLevel(int var) { fDataOutputLevel = var; }
    
 protected:
  
  Int_t DoClustering();
  Int_t MatchTrack(THaTrack*);
  
  void ClearHodoOutput(SBSTimingHodoscopeOutput &var);
  
  Double_t fHorizPosBarCut;
  Double_t fTimeRef;//? really useful?
  Double_t fTimeBarCut;
  
  int fClusMaxSize;
  Double_t fMaxYposDiffCluster;// maximum Y pos difference to incorporate a new bar in a cluster 
  Double_t fMaxTimeDiffCluster;// maximum time difference to incorporate a new bar in a cluster 
  
  Double_t AttLength = 1.0/3.8;//380cm for ej200. units per m?
  /* speed of light in scint? what is n for ej200 */
  // n = 1.58 => n=c/v =>v=c/n
  //  Double_t n=1.58;

  std::vector<Double_t> fvScint; 
  std::vector<Double_t> ftDiff0;

  //  Double_t fvScint; //default to 0.454c, later we might want to define separately for different bars?
  //  Double_t ftDiff0; //offset of time difference to align horizontal position from time difference with horizontal projection of tracks
  Double_t fTrackMatchCutX;
  Double_t fTrackMatchCutY;
  
  /* std::vector<SBSTimingHodoscopePMT> fPMTMap; */
  std::vector<SBSTimingHodoscopePMT*> fPMTMapL;
  std::vector<SBSTimingHodoscopePMT*> fPMTMapR;
  std::vector<SBSTimingHodoscopeBar*> fBars;
  Int_t fTDCBarOffset;
  Int_t fADCBarOffset;
  Int_t fTDCRefLeL;
  Int_t fTDCRefLeR;
  Double_t fTDCWinMin;
  Double_t fTDCWinMax;
  Double_t fTotMin;
  Double_t fTotMax;
  
  std::vector<Int_t>   fGoodBarIDsTDC;
  std::vector<Double_t> fGoodBarTDCmean;
  std::vector<Double_t> fGoodBarTDCdiff;
  std::vector<Double_t> fGoodBarTDCvpos;
  std::vector<Double_t> fGoodBarTDCpos;
  std::vector<Double_t> fGoodBarTDCLle;
  std::vector<Double_t> fGoodBarTDCLleW;
  std::vector<Double_t> fGoodBarTDCLte;
  std::vector<Double_t> fGoodBarTDCLteW;
  std::vector<Double_t> fGoodBarTDCLtot;
  std::vector<Double_t> fGoodBarTDCLtotW;
  std::vector<Double_t> fGoodBarTDCRle;
  std::vector<Double_t> fGoodBarTDCRleW;
  std::vector<Double_t> fGoodBarTDCRte;
  std::vector<Double_t> fGoodBarTDCRteW;
  std::vector<Double_t> fGoodBarTDCRtot;
  std::vector<Double_t> fGoodBarTDCRtotW;

  std::vector<Int_t>   fGoodBarIDsADC;
  std::vector<Double_t> fGoodBarADCmean;
  std::vector<Double_t> fGoodBarADCLa;
  std::vector<Double_t> fGoodBarADCLap;
  std::vector<Double_t> fGoodBarADCLac;
  std::vector<Double_t> fGoodBarADCRa;
  std::vector<Double_t> fGoodBarADCRap;
  std::vector<Double_t> fGoodBarADCRac;
  
  //cluster output only:
  /*
  std::vector<Int_t> fClusterMult;
  std::vector<Double_t> fClusterXmean;
  std::vector<Double_t> fClusterYmean;
  std::vector<Double_t> fClusterTmean;
  std::vector<Double_t> fClusterToTmean;
  */
  
  // Mapping (see also fDetMap)
  // maps of time walk parameters - in row, col, lay
  std::vector<std::vector<std::vector<Double_t>>> fTimeWalkPar0;
  std::vector<std::vector<std::vector<Double_t>>> fTimeWalkPar1;
  
  SBSTimingHodoscopeOutput fMainClus;
  SBSTimingHodoscopeOutput fMainClusBars;
  SBSTimingHodoscopeOutput fOutClus;
  std::vector<SBSTimingHodoscopeCluster*> fClusters; 
  
  Int_t fDataOutputLevel;   //0 (default): only main cluster; 1: (0)+ main cluster bars; 2: (1)+all clusters; ": (2)+all bars
  
  ClassDef(SBSTimingHodoscope,5)  // Describes scintillator plane with F1TDC as a detector
};

inline Double_t SBSTimingHodoscope::GetVVal(std::vector<Double_t> &v, UInt_t i)
{
  if (v.size() > i) {
    return v[i];
  }
  return 0.0;
}

inline Int_t SBSTimingHodoscope::GetVVal(std::vector<Int_t> &v, UInt_t i)
{
  if (v.size() > i) {
    return v[i];
  }
  return 0.0;
}

inline Int_t SBSTimingHodoscope::GetID(){
  return GetVVal(fMainClus.id);
}

inline Int_t SBSTimingHodoscope::GetSize(){
  return GetVVal(fMainClus.n);
}

inline Double_t SBSTimingHodoscope::GetX(){
  return GetVVal(fMainClus.x);
}

inline Double_t SBSTimingHodoscope::GetY(){
  return GetVVal(fMainClus.y);
}

inline Double_t SBSTimingHodoscope::GetT(){
  return GetVVal(fMainClus.t);
}
  
inline Double_t SBSTimingHodoscope::GetToT(){
  return GetVVal(fMainClus.tot);
}




#endif
