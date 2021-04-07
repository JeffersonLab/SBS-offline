#ifndef SBSGenericDetector_h
#define SBSGenericDetector_h

////////////////////////////////////////////////////////////////////////////////
//
//SBSGenericDetector
//
// A generic detector class which could contain the following types
// of data:
// - ADC
//  - Can be single valued OR
//  - Can contain pulse information such as integral, amplitude and time
//  - Can be using Waveform data
// - TDC
//  - Can be single valued (leading edge) OR:
//  - Can contain both leading and trailing edge and Time-Over-Threshold
//
//  Organized into row, col, and layer structure (but makes no assumption about
//  number of rows, cols or layers are constant throughout).
//
////////////////////////////////////////////////////////////////////////////////

#include "THaNonTrackingDetector.h"
#include "TRotation.h"
#include "TVector3.h"
#include "THaDetMap.h"
#include "SBSElement.h"

namespace SBSModeADC {
  enum Mode{
    kNone,
    kADC,       //< Contains pulse information also
    kADCSimple, //< Does not contain pulse information (nor reference info)
    kWaveform,  //< Contains waveform data
  };
}

namespace SBSModeTDC {
  enum Mode{
    kNone,
    kIgnore,        //< Useful to preserve DB but ignore ADCs otherwise
    kTDC,       //< Includes Trailing edge (and ToT)
    kTDCSimple,  //< Does not contain trailing edge, and hence no ToT info
  };
}

class SBSGenericDetector : public THaNonTrackingDetector {
  //class SBSGenericDetector : public THaShower {

public:
  SBSGenericDetector( const char* name, const char* description = "",
      THaApparatus* a = NULL);
  virtual ~SBSGenericDetector();

  virtual void ClearEvent();
  virtual void ClearOutputVariables();

  void SetModeADC(SBSModeADC::Mode mode);
  void SetModeTDC(SBSModeTDC::Mode mode) { fModeTDC = mode; }
  void SetDisableRefADC(Bool_t b) { fDisableRefADC = b; }
  void SetDisableRefTDC(Bool_t b) { fDisableRefTDC = b; }
  void SetStoreRawData(Bool_t var) { fStoreRawData = var; }
  void SetStoreEmptyElements(Bool_t b) { fStoreEmptyElements = b; }

  Bool_t WithTDC() { return fModeTDC != SBSModeTDC::kNone; };
  Bool_t WithADC() { return fModeADC != SBSModeADC::kNone; };

  // Standard apparatus re-implemented functions
  virtual Int_t      Decode( const THaEvData& );
  virtual Int_t      CoarseProcess(TClonesArray& tracks);
  virtual Int_t      FineProcess(TClonesArray& tracks);

  virtual Int_t      DecodeADC( const THaEvData&, SBSElement *blk,
      THaDetMap::Module *d, Int_t chan);
  virtual Int_t      DecodeTDC( const THaEvData&, SBSElement *blk,
      THaDetMap::Module *d, Int_t chan);

  // Utility functions
  // Get index of element
  //Int_t blkidx(Int_t row, Int_t col, Int_t layer = 0);
  // Get corresponding row, col, layer of element from index
  //void blkrcl(Int_t index, Int_t &row, Int_t &col, Int_t &layer);

protected:

  virtual Int_t  ReadDatabase( const TDatime& date );
  virtual Int_t  DefineVariables( EMode mode = kDefine );

  // Configuration
  Int_t  fNrows;        ///< Number of rows
  std::vector<Int_t>  fNcols; ///< Number of columns per row
  Int_t  fNlayers;      ///< Number of layers (in z-direction)
  SBSModeADC::Mode fModeADC;      //< ADC Mode
  SBSModeTDC::Mode fModeTDC;      //< TDC Mode
  Bool_t fDisableRefADC; //< Reference ADC may be optionally disabled
  Bool_t fDisableRefTDC; //< Reference TDC may be optionally disabled
  Bool_t fStoreEmptyElements; //< Do not store data for empty elements in rootfile

  // Mapping (see also fDetMap)
  UShort_t   fChanMapStart; ///< Starting number for element number (i.e. 0 or 1)
  std::vector<std::vector<Int_t> > fChanMap; //< Maps modules in THaDetMap to calorimeter element number

  // Output variables
  // Note [] means it is a vector with variable width
  std::vector<Int_t> fRow;        //< [] row number for given element
  std::vector<Int_t> fCol;        //< [] col number for given element
  std::vector<Int_t> fLayer;      //< [] layer number for given element
  // ADC Integral
  std::vector<SBSData::ADC*> fA_all;        //< [] ADC structures
  std::vector<Float_t> fA;      //< [] ADC pulse integral
  std::vector<Float_t> fA_amp;  //< [] ADC pulse amplitude (i.e. peak value)
  std::vector<Float_t> fA_time; //< [] ADC pulse time
  // ADC Samples
  std::vector<Int_t> fNsamps;     //< [] Number of samples for each element
  std::vector<Int_t> fSampsIdx;   //< []*fNsamps Index of the corresponding element
  std::vector<Float_t> fSamps;    //< []*fNsamps ADC samples for each element
  std::vector<Float_t> fSamps_c;  //< []*fNsamps ADC calibrated samples for each element
  // TDC values
  std::vector<Float_t> fTDC;    //< [] TDC good hit times
  std::vector<Float_t> fTDC_te; //< [] TDC good hit  trailing edge times
  std::vector<Float_t> fTDC_ToT; //< [] TDC good hit Time-Over-Threshold

  // Blocks, where the grid is just for easy axis to the elements by row,col,layer
  std::vector<SBSElement*> fElements;
  std::vector<SBSElement*> fRefElements; //< Reference elements (for TDCs)
  std::vector<std::vector<std::vector<SBSElement*> > > fElementGrid;
  // Clusters for this event

  // Other parameters
  Bool_t fCoarseProcessed;  //< Was CourseProcessed already called?
  Bool_t fFineProcessed;    //< Was fine processed already called

  //Gain correction 
  Float_t   fConst;     // const from gain correction 
  Float_t   fSlope;     // slope for gain correction
  Float_t   fAccCharge; // accumulated charge

  // Per event data
  Int_t      fNhits;     ///< Number of hits in event

  // Flags for enabling and disabling various features
  Bool_t    fStoreRawData; ///< Store the raw data in the root tree?

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

  ClassDef(SBSGenericDetector,0)     //Generic shower detector class
};

/*inline Int_t SBSGenericDetector::blkidx(Int_t row, Int_t col, Int_t layer)
{
  return fNlayers*(fNcols[row]*row + col) + layer;
}

inline void SBSGenericDetector::blkrcl(Int_t index, Int_t &row, Int_t &col,
    Int_t &layer)
{
  row = index/(fNlayers*fNcols);
  index -= row*fNlayers*fNcols;
  col = index/fNlayers;
  layer = index%fNlayers;
}
*/

////////////////////////////////////////////////////////////////////////////////
// Specify some default sub-classes available to the user

#endif // SBSGenericDetector_h
