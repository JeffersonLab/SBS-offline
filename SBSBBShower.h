#ifndef ROOT_SBSBBShower
#define ROOT_SBSBBShower

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// SBSBBShower                                                               //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "THaPidDetector.h"
#include "THaShower.h"
#include "SBSBBShowerCluster.h"
#include "TRotation.h"
#include "TVector3.h"

class SBSBBShower : public THaShower {//THaPidDetector {

 public:
  SBSBBShower( const char* name, const char* description = "",
	     THaApparatus* a = NULL );
  virtual ~SBSBBShower();

  virtual Int_t      Decode( const THaEvData& );
  virtual Int_t      CoarseProcess(TClonesArray& tracks);
  virtual Int_t      FineProcess(TClonesArray& tracks);
  
  /* Int_t      GetNclust() const { return fNclust; } */
  /* Int_t      GetNhits() const  { return fNhits; } */
  /* Float_t    GetE(int i) const { return fE[i]; } */
  /* Float_t    GetX(int i) const { return fX[i]; } */
  /* Float_t    GetY(int i) const { return fY[i]; } */
  
  Int_t      GetNBlocks() { return (fNrows * fNcols);}
  Float_t    GetBlockX( Int_t i )  { if(i < fNrows*fNcols) return fBlocks[i]->GetX(); else return 0.0;}
  Float_t    GetBlockY( Int_t i )  { if(i < fNrows*fNcols) return fBlocks[i]->GetY(); else  return 0.0;}
  
  Float_t    GetBlockdX()  {return fdX;}
  Float_t    GetBlockdY()  {return fdY;}
  Float_t    GetBlockdZ()  {return fdZ;}
  
  
  Float_t    GetBlockA_c( Int_t i ) const { return fA_c[i]; }
  
  Int_t    GetNRows() {return fNrows;}
  Int_t    GetNCols() {return fNcols;}
  Int_t    BlockColRowToNumber( Int_t col, Int_t row );
  
  // Blocks should have a Z!!!
  
  SBSBBShowerCluster* GetClust(Int_t i) { return fClusters[i]; }
  
  void       AddCluster(SBSBBShowerCluster* clus);
  void       RemoveCluster(int i);
  void       AddCluster(SBSBBShowerCluster& clus);
  
  void       LoadMCHitAt( Double_t x, Double_t y, Double_t E );
  
 protected:
  
  Bool_t fCoarseProcessed;
  Bool_t fFineProcessed;

  // Maximum number of clusters
  Int_t fMaxNClust;
  //Int_t fNclusters;
  // Mapping (see also fDetMap)
  //UShort_t*  fNChan;     // Number of channels for each module
  //std::vector< std::vector<UShort_t> > fChanMap; 
  // Logical channel   UShort_t** fChanMap;   // Logical channel numbers 
  Int_t chanmap_start;
  // Configuration
  //Int_t      fNclublk;   // Max. number of blocks composing a cluster
  //Int_t      fNrows;     // Number of rows
  Int_t      fNcols;     // Number of columns
  
  // Geometry
  //Float_t*   fBlockX;    // [fNelem] x positions (cm) of block centers
  //Float_t*   fBlockY;    // [fNelem] y positions (cm) of block centers

  // Calibration
  //Float_t*   fPed;       // [fNelem] Pedestals for each block
  //Float_t*   fGain;      // [fNelem] Gains for each block

  //Gain correction 
  Float_t   gconst;     // const from gain correction 
  Float_t   gslope;     // slope for gain correction
  Float_t   acc_charge; // accumulated charge
  
  // Cluster parameters Other parameters
  Float_t fClusRadius; // radius for cluster search around the max 
  Int_t fClusBlockRadX;
  Int_t fClusBlockRadY;
  
  // Per-event data
  //Int_t      fNhits;     // Number of hits
  //Float_t*   fA;         // [fNelem] Array of ADC amplitudes of blocks
  //Float_t*   fA_p;       // [fNelem] Array of ADC minus pedestal values of blocks
  //Float_t*   fA_c;       // [fNelem] Array of corrected ADC amplitudes of blocks
  //Float_t    fAsum_p;    // Sum of blocks ADC minus pedestal values
  //Float_t    fAsum_c;    // Sum of blocks corrected ADC amplitudes
  //Int_t      fNclust;    // Number of clusters
  Float_t*   fE_cl;        // [fNClust] Energy (MeV) of clusters
  Float_t*   fX_cl;        // [fNClust] x position (m) of clusters
  Float_t*   fY_cl;        // [fNClust] y position (m) of clusters 
  //Float_t*   fXtarg;     // [fNClust] x position (m) of clusters in target coords
  //Float_t*   fYtarg;     // [fNClust] y position (m) of clusters in target coords
  //Float_t*   fZtarg;     // [fNClust] z position (m) of clusters in target coords
  Int_t*     fMult_cl;     // [fNClust]  Number of blocks in all clusters
  Int_t**    fNblk_cl;     // [fNclublk] Block numbers composing main cluster
  Float_t**  fEblk_cl;     // [fNclublk] Block energies composing main cluster
  
  Float_t    fTRX;       // x position of track cross point
  Float_t    fTRY;       // y position of track cross point

  Float_t    fE_corr;        // [fNClust] Corrected Energy (MeV) of clusters
  Float_t*   fE_cl_corr;        // [fNClust] Corrected Energy (MeV) of clusters
 
  Double_t   fdX;
  Double_t   fdY;
  Double_t   fdZ;
  
  SBSShowerBlock** fBlocks; //[fNelem] Array of blocks
  SBSBBShowerCluster** fClusters; //[fMaxNClust] 
  //SBSShowerBlock*** fBlkGrid; //[fNrows]

  Float_t   fEres;        // [fNClust] Energy resolution of main cluster
  Float_t   fXres;        // [fNClust] x position (m) of main cluster
  Float_t   fYres;        // [fNClust] y position (m) of main cluster
  
  Float_t*   fE_cl_res;        // [fNClust] Energy resolution of clusters
  Float_t*   fX_cl_res;        // [fNClust] x position (m) of clusters
  Float_t*   fY_cl_res;        // [fNClust] y position (m) of clusters 
  //TRotation  fDetToTarg;
  //TVector3   fDetOffset;

  // Useful derived quantities for internal usage.
  Double_t fThrADC;
  
  Double_t tan_angle, sin_angle, cos_angle;
  
  bool fMultClus;// allow multiple clustering
  bool fMCdata;// easy way to enable/disable the use of MC data.
  
  //void           ClearEvent();
  virtual void   Clear( Option_t* opt="" );
  void           DeleteArrays();
  virtual Int_t  ReadDatabase( const TDatime& date );
  virtual Int_t  DefineVariables( EMode mode = kDefine );
  
  ClassDef(SBSBBShower,0)     //Generic shower detector class
};

////////////////////////////////////////////////////////////////////////////////

#endif
