#ifndef SBSCalorimeter_h
#define SBSCalorimeter_h

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// SBSCalorimeter                                                            //
//                                                                           //
// A generic calorimeter class which could contain the following types       //
// of data:                                                                  //
//  - Single-valued ADC                                                      //
//  - Single-valued ADC + TDC                                                //
//  - Multi-valued ADC                                                       //
//  - Multi-valued ADC + TDC                                                 //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "THaNonTrackingDetector.h"
//#include "SBSCalorimeterCluster.h"
#include "SBSCalorimeterBlock.h"
#include "TRotation.h"
#include "TVector3.h"
#include "THaDetMap.h"


class SBSCalorimeter : public THaNonTrackingDetector {

public:
  SBSCalorimeter( const char* name, const char* description = "",
      THaApparatus* a = NULL);
  virtual ~SBSCalorimeter();

  void           ClearEvent();
  void SetWithADCSamples(Bool_t var) { fWithADCSamples = var; }
  void SetWithTDC(Bool_t var)        { fWithTDC = var; }

  // Standard apparatus re-implemented functions
  virtual Int_t      Decode( const THaEvData& );
  virtual Int_t      CoarseProcess(TClonesArray& tracks);
  virtual Int_t      FineProcess(TClonesArray& tracks);

  virtual Int_t      DecodeADC( const THaEvData&, SBSCalorimeterBlock *blk,
      THaDetMap::Module *d, Int_t chan);
  virtual Int_t      DecodeTDC( const THaEvData&, SBSCalorimeterBlock *blk,
      THaDetMap::Module *d, Int_t chan);

  // A structure to hold the output that will go to the tree
  struct OutputData {
    std::vector<Int_t> fRow;
    std::vector<Int_t> fCol;
    std::vector<Int_t> fLayer;
    // ADC Integral
    std::vector<Double_t> fA;
    std::vector<Double_t> fA_p;
    std::vector<Double_t> fA_c;
    // ADC Samples
    std::vector<Int_t> fSampsIdx;
    std::vector<Int_t> fNSamps;
    std::vector<Double_t> fSamps;
    std::vector<Double_t> fSamps_p;
    std::vector<Double_t> fSamps_c;
    // TDC values
    std::vector<Double_t> fTDC;
    std::vector<Double_t> fTDC_c;
    void ClearEvent();
  };

  // Utility functions
  // Get index of block
  Int_t blkidx(Int_t row, Int_t col, Int_t layer = 0);
  // Get corresponding row, col, layer of block from index
  void blkrcl(Int_t index, Int_t &row, Int_t &col, Int_t &layer);

protected:

  virtual Int_t  ReadDatabase( const TDatime& date );
  virtual Int_t  DefineVariables( EMode mode = kDefine );

  // Configuration
  Int_t  fNclublk;      ///< Max number of blocks composing a cluster
  Int_t  fNrows;        ///< Number of rows
  Int_t  fNcols;        ///< Number of columns
  Int_t  fNlayers;      ///< Number of layers (in z-direction)
  Bool_t fWithTDC;      ///< Does this calorimeter have TDC readout?
  Bool_t fWithADCSamples; ///< Does this calorimeter have multi-valued ADC readout?

  // Mapping (see also fDetMap)
  UShort_t   fChanMapStart; ///< Starting number for block number (i.e. 0 or 1)
  std::vector<std::vector<UShort_t> > fChanMap;
  //std::map<UShort_t,std::map<UShort_t,std::vector<UShort_t> > > fChanMap;

  // Output vectors
  OutputData fDataOut;

  // Blocks
  std::vector<SBSCalorimeterBlock*> fBlocks;


  // Per event data
  Int_t      fNhits;     ///< Number of hits in event

  /*
     Int_t      GetNclust() const { return fNclust; }
     Int_t      GetNhits() const  { return fNhits; }
     Float_t    GetE(int i) const { return fE[i]; }
     Float_t    GetX(int i) const { return fX[i]; }
     Float_t    GetY(int i) const { return fY[i]; }

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

  SBSCalorimeterCluster* GetClust(Int_t i) { return fClusters[i]; }

  void       AddCluster(SBSCalorimeterCluster* clus);
  void       RemoveCluster(int i);
  void       AddCluster(SBSCalorimeterCluster& clus);

  void       LoadMCHitAt( Double_t x, Double_t y, Double_t E );

  protected:

  Bool_t fCoarseProcessed;
  Bool_t fFineProcessed;

  // Maximum number of clusters
  Int_t fMaxNClust;

  // Mapping (see also fDetMap)
  UShort_t*  fNChan;     // Number of channels for each module
  UShort_t** fChanMap;   // Logical channel numbers 

  // Configuration
  Int_t      fNclublk;   // Max. number of blocks composing a cluster
  Int_t      fNrows;     // Number of rows
  Int_t      fNcols;     // Number of columns

  // Geometry
  Float_t*   fBlockX;    // [fNelem] x positions (cm) of block centers
  Float_t*   fBlockY;    // [fNelem] y positions (cm) of block centers

  // Calibration
  Float_t*   fPed;       // [fNelem] Pedestals for each block
  Float_t*   fGain;      // [fNelem] Gains for each block

  //Gain correction 
  Float_t   gconst;     // const from gain correction 
  Float_t   gslope;     // slope for gain correction
  Float_t   acc_charge; // accumulated charge

  // Other parameters
  Float_t    fEmin;      // Minimum energy for a cluster center

  // Per-event data
  Float_t*   fA;         // [fNelem] Array of ADC amplitudes of blocks
  Float_t*   fA_p;       // [fNelem] Array of ADC minus pedestal values of blocks
  Float_t*   fA_c;       // [fNelem] Array of corrected ADC amplitudes of blocks
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

  SBSShowerBlock** fBlocks; //[fNelem] Array of blocks
  SBSCalorimeterCluster** fClusters; //[fMaxNClust] 
  SBSShowerBlock*** fBlkGrid; //[fNrows]

  //TRotation  fDetToTarg;
  //TVector3   fDetOffset;

  // Useful derived quantities for internal usage.

  Double_t tan_angle, sin_angle, cos_angle;

  void           DeleteArrays();
  */

private:
  // Simple and quick routine to init and clear most vectors
  // (of integers, floats, doubles, etc...)
  // Reset/Init 1D vector
  template<class T>
  void InitVector(std::vector<T> &vec, T val = 0, size_t n = 0) {
    vec.resize(n);
    ResetVector(vec,val);
  }

  template<class T>
  void ResetVector(std::vector<T> &vec, T val = 0, size_t n = 0) {
    if(n > 0) {
      vec.clear();
      vec.resize(n);
    }
    for(size_t i = 0; i < vec.size(); i++) {
      vec[i] = val;
    }
  }

  // Reset 2D vector
  template<class T>
  void ResetVector(std::vector<std::vector<T> > &vec, T val = 0,
      size_t nr = 0, size_t nc = 0) {
    if(nr > 0) {
      vec.clear();
      vec.resize(nr);
    }
    for(size_t i = 0; i < vec.size(); i++) {
      ResetVector(vec[i],val,nc);
    }
  }

  ClassDef(SBSCalorimeter,0)     //Generic shower detector class
};

inline Int_t SBSCalorimeter::blkidx(Int_t row, Int_t col, Int_t layer)
{
  return fNlayers*(fNcols*row + col) + layer;
}

inline void SBSCalorimeter::blkrcl(Int_t index, Int_t &row, Int_t &col,
    Int_t &layer)
{
  row = index/(fNlayers*fNcols);
  index -= row*fNlayers*fNcols;
  col = index/fNlayers;
  layer = index%fNlayers;
}

////////////////////////////////////////////////////////////////////////////////

#endif // SBSCalorimeter_h
