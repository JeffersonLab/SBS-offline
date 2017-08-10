#ifndef ROOT_SBSHCal
#define ROOT_SBSHCal

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// SBSHCal                                                               //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "THaPidDetector.h"
#include "SBSShowerBlock.h"
//2//#include "SBSBBShowerCluster.h"
#include "TRotation.h"
#include "TVector3.h"

class SBSHCal : public THaPidDetector {

 public:
  SBSHCal( const char* name, const char* description = "",
	     THaApparatus* a = NULL );
  virtual ~SBSHCal();

  virtual Int_t      Decode( const THaEvData& );
  virtual Int_t      CoarseProcess(TClonesArray& tracks);
  virtual Int_t      FineProcess(TClonesArray& tracks);
  //virtual Int_t      DecodeFADC( const THaEvData& );
  
          Int_t      GetNclust() const { return fNclust; }
          Int_t      GetNhits() const  { return fNhits; }
          Float_t    GetE(int i) const { return fE[i]; }
          Float_t    GetX(int i) const { return fX[i]; }
          Float_t    GetY(int i) const { return fY[i]; }

	  Int_t      GetNBlocks() { return (fNrows * fNcols);}
	  //2//Float_t    GetBlockX( Int_t i )  { if(i < fNrows*fNcols) return fModules[i]->GetX(); else return 0.0;}
	  //2//Float_t    GetBlockY( Int_t i )  { if(i < fNrows*fNcols) return fModules[i]->GetY(); else  return 0.0;}
	  
	  Float_t    GetBlockdX()  {return fdX;}
	  Float_t    GetBlockdY()  {return fdY;}
	  Float_t    GetBlockdZ()  {return fdZ;}
	  

	  //Float_t    GetBlockA_c( Int_t i ) const { return fA_c[i]; }

	  Int_t    GetNRows() {return fNrows;}
	  Int_t    GetNCols() {return fNcols;}
	  Int_t    BlockColRowToNumber( Int_t col, Int_t row );

	  // Blocks should have a Z!!!
	  
	  //2//SBSBBShowerCluster* GetClust(Int_t i) { return fClusters[i]; }
  
	  //2//void       AddCluster(SBSBBShowerCluster* clus);
	  void       RemoveCluster(int i);
	  //2//void       AddCluster(SBSBBShowerCluster& clus);

	  void       LoadMCHitAt( Double_t x, Double_t y, Double_t E );
	  
 protected:
	  
  Bool_t fCoarseProcessed;
  Bool_t fFineProcessed;

  // Maximum number of clusters
  Int_t fMaxNClust;

  // Mapping (see also fDetMap)
  UShort_t*  fNChan;     // Number of channels for each module
  std::vector< std::vector<UShort_t> > fChanMap; // Logical channel numbers
                                                 // for each detector map module

  // Configuration
  Int_t      fNclublk;   // Max. number of blocks composing a cluster
  Int_t      fNrows;     // Number of rows
  Int_t      fNcols;     // Number of columns

  //Gain correction 
   Float_t   gconst;     // const from gain correction 
   Float_t   gslope;     // slope for gain correction
   Float_t   acc_charge; // accumulated charge

  // Other parameters
  Float_t    fEmin;      // Minimum energy for a cluster center

  // Per-event data
  Int_t      fNhits;     // Number of hits

  std::vector<Int_t>   fNumSamples; // [fNelem] Number of samples in each ADC/module
  std::vector<std::vector<Double_t> > fASamples; // [fNumSamples] Raw ADC samples
  std::vector<std::vector<Double_t> > fASamplesPed; // [fNumSamples] Ped correct ADC samples
  std::vector<std::vector<Double_t> > fASamplesCal; // [fNumSamples] Calibrated ADC samples
  std::vector<Float_t> fPed; // [fNelem] Pedestal for each module
  std::vector<Float_t> fGain; // [fNelem] Gain for each module

  Float_t    fAsum_p;    // Sum of blocks ADC minus pedestal values
  Float_t    fAsum_c;    // Sum of blocks corrected ADC amplitudes
  Int_t      fNclust;    // Number of clusters
  Float_t*   fE;        // [fNClust] Energy (MeV) of clusters
  Float_t*   fX;        // [fNClust] x position (m) of clusters
  Float_t*   fY;        // [fNClust] y position (m) of clusters 
  //Float_t*   fXtarg;     // [fNClust] x position (m) of clusters in target coords
  //Float_t*   fYtarg;     // [fNClust] y position (m) of clusters in target coords
  //Float_t*   fZtarg;     // [fNClust] z position (m) of clusters in target coords
  Int_t*     fMult;      // [fNClust]  Number of blocks in main cluster
  Int_t*     fNblk;      // [fNclublk] Numbers of blocks composing main cluster
  Float_t*   fEblk;      // [fNclublk] Energies of blocks composing main cluster
  Float_t    fTRX;       // x position of track cross point
  Float_t    fTRY;       // y position of track cross point

  Float_t*   fE_c;        // [fNClust] Corrected Energy (MeV) of clusters
 
  Double_t   fdX;
  Double_t   fdY;
  Double_t   fdZ;

  std::vector<SBSShowerBlock*> fBlocks; //[fNelem] Array of modules
  //SBSBBShowerCluster** fClusters; //[fMaxNClust] 
  std::vector<std::vector<SBSShowerBlock*> > fBlkGrid; //[fNrows]

  //TRotation  fDetToTarg;
  //TVector3   fDetOffset;

  // Useful derived quantities for internal usage.

  Double_t tan_angle, sin_angle, cos_angle;

  void           ClearEvent();
  void           DeleteArrays();
  virtual Int_t  ReadDatabase( const TDatime& date );
  virtual Int_t  DefineVariables( EMode mode = kDefine );

 private:
  void SetBlockADCSamples(Int_t block, std::vector<UInt_t> samples);

  // Simple and quick routine to clear most vectors (of integers, floats,
  // doubles, etc...)
  template<class T>
  void ResetVector(std::vector<T> &vec, T val = 0) {
    for(size_t i = 0; i < vec.size(); i++) {
      vec[i] = val;
    }
  }

  // Reset 2D vector
  template<class T>
  void ResetVector(std::vector<std::vector<T> > &vec, T val = 0) {
    for(size_t i = 0; i < vec.size(); i++) {
      ResetVector(vec[i],val);
    }
  }

  ClassDef(SBSHCal,0)     //Generic shower detector class
};

////////////////////////////////////////////////////////////////////////////////

#endif
