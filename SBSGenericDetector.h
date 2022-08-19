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
#include "Helper.h"

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
// structure for sorting TDC hits
struct TDCHits {
  UInt_t edge;
  UInt_t rawtime;
};

// This structure has output data when the user wants every hit to be stored
// in the rootfile.
struct SBSGenericOutputData {
  // Note [] means it can be variable sized data per event/module
  // Module info
  std::vector<Int_t> TDCrow;       //< [] row
  std::vector<Int_t> TDCcol;         //< [] col
  std::vector<Int_t> TDClayer;       //< [] layer
  std::vector<Int_t> TDCelemID;      //< [] element ID
  std::vector<Int_t> ADCrow;         //< [] row
  std::vector<Int_t> ADCcol;         //< [] col
  std::vector<Int_t> ADClayer;       //< [] layer
  std::vector<Int_t> ADCelemID;      //< [] element ID
  // ADC variables
  std::vector<Int_t> ped;         //< [] pedestal
  std::vector<Int_t> a_mult;         //< [] ADC # of hits per channel
  std::vector<Double_t> a;         //< [] ADC integral
  std::vector<Double_t> a_p;         //< [] ADC integral -pedestal
  std::vector<Double_t> a_c;         //< [] (ADC integral -pedestal)*calib
  std::vector<Double_t> a_amp;     //< [] ADC pulse amplitude
  std::vector<Double_t> a_amp_p;     //< [] ADC pulse amplitude -pedestal
  std::vector<Double_t> a_amp_c;     //< [] ADC pulse amplitude -pedestal
  std::vector<Double_t> a_amptrig_p;     //< [] ADC pulse amplitude -pedestal
  std::vector<Double_t> a_amptrig_c;     //< [] ADC pulse amplitude -pedestal
  std::vector<Double_t> a_time;    //< [] ADC pulse time
  // TDC variables
  std::vector<Int_t> t_mult;         //< [] TDC # of hits per channel
  std::vector<Double_t> t;         //< [] TDC (leading edge) time
  std::vector<Double_t> t_te;      //< [] TDC trailing edge time
  std::vector<Double_t> t_ToT;     //< [] TDC Time-Over-Threshold
  // Waveform variables
  std::vector<Int_t> samps_elemID;      //< [] Element ID of samples
  std::vector<Int_t> nsamps;      //< [] Number of ADC samples
  std::vector<Int_t> sidx;        //< [] Index of start of ADC samples in for this row-col-layer
  std::vector<Double_t> samps;     //< []*nsamps ADC samples

  // Quick clear class
  void clear() {
    TDCrow.clear();
    TDCcol.clear();
    TDCelemID.clear();
    TDClayer.clear();
    ADCrow.clear();
    ADCcol.clear();
    ADCelemID.clear();
    ADClayer.clear();
    ped.clear();
    a_mult.clear();
    a.clear();
    a_p.clear();
    a_c.clear();
    a_amp.clear();
    a_amp_p.clear();
    a_amp_c.clear();
    a_amptrig_p.clear();
    a_amptrig_c.clear();
    a_time.clear();
    t.clear();
    t_mult.clear();
    t_te.clear();
    t_ToT.clear();
    nsamps.clear();
    sidx.clear();
    samps.clear();
    samps_elemID.clear();
  }
};

class SBSGenericDetector : public THaNonTrackingDetector {
  //class SBSGenericDetector : public THaShower {

public:
  explicit SBSGenericDetector( const char* name, const char* description = "",
                               THaApparatus* a = nullptr);
  virtual ~SBSGenericDetector();

  virtual void Clear( Option_t* opt="" );

  void SetModeADC(SBSModeADC::Mode mode);
  void SetModeTDC(SBSModeTDC::Mode mode) { fModeTDC = mode; }
  void SetDisableRefADC(Bool_t b) { fDisableRefADC = b; }
  void SetDisableRefTDC(Bool_t b) { fDisableRefTDC = b; }
  void SetStoreRawHits(Bool_t var) { fStoreRawHits = var; }
  void SetStoreEmptyElements(Bool_t b) { fStoreEmptyElements = b; }

  Bool_t WithTDC() { return fModeTDC != SBSModeTDC::kNone; };
  Bool_t WithADC() { return fModeADC != SBSModeADC::kNone; };

  // Standard apparatus re-implemented functions
  virtual Int_t      Decode( const THaEvData& );
  virtual Int_t      CoarseProcess(TClonesArray& tracks);
  virtual Int_t      FineProcess(TClonesArray& tracks);

  virtual Int_t      DecodeADC( const THaEvData&, SBSElement *blk,
				THaDetMap::Module *d, Int_t chan, Bool_t IsRef);
  virtual Int_t      DecodeTDC( const THaEvData&, SBSElement *blk,
      THaDetMap::Module *d, Int_t chan, Bool_t IsRef);

  // Utility functions
  // Can be re-implemented by other classes to specify a different
  // SBSElement sub-class (i.e. useful when one wants to chang  the logic
  // in SBSElement::CoarseProcess()
  virtual SBSElement* MakeElement(Double_t x, Double_t y, Double_t z, Int_t row,
      Int_t col, Int_t layer, Int_t id = 0);
  
  Double_t SizeRow() const { return fSizeRow; };
  Double_t SizeCol() const { return fSizeCol; };
  
protected:

  virtual Int_t  ReadDatabase( const TDatime& date );
  virtual Int_t  DefineVariables( EMode mode = kDefine );
  virtual Int_t  FindGoodHit(SBSElement *); // 

  // Configuration
  Int_t  fNRefElem;        ///< Number of Ref Time elements
  Int_t  fNrows;        ///< Number of rows
  std::vector<Int_t>  fNcols; ///< Number of columns per row
  Double_t fSizeRow;
  Double_t fSizeCol;
  Int_t fNcolsMax;      ///< Max number of columns out of all rows
  Int_t  fNlayers;      ///< Number of layers (in z-direction)
  SBSModeADC::Mode fModeADC;      //< ADC Mode
  SBSModeTDC::Mode fModeTDC;      //< TDC Mode
  Bool_t fDisableRefADC; //< Reference ADC may be optionally disabled
  Bool_t fDisableRefTDC; //< Reference TDC may be optionally disabled
  Bool_t fStoreEmptyElements; //< Do not store data for empty elements in rootfile
  Bool_t fIsMC; // flag to indicate if data are simulated;
  Int_t fF1_RollOver;
  Int_t fF1_TimeWindow;
  
  // Mapping (see also fDetMap)
  UShort_t   fChanMapStart; ///< Starting number for element number (i.e. 0 or 1)
  std::vector<std::vector<Int_t> > fChanMap; //< Maps modules in THaDetMap to element number
  std::vector<std::vector<Int_t> > fRefChanMap; //< Maps modules in THaDetMap to Reference time
  std::vector<Bool_t> fModuleRefTimeFlag; //< If module has Reftime set to true
  std::vector<Int_t> fRefChanLo; //< Module Reftime Low channel number
  std::vector<Int_t> fRefChanHi; //< Module Reftime High channel number

  // Output variable containers
  SBSGenericOutputData fGood;     //< Good data output
  SBSGenericOutputData fRaw;      //< All hits
  SBSGenericOutputData fRefGood;     //< Good Ref time data output
  SBSGenericOutputData fRefRaw;      //< All Ref time hits

  // Blocks, where the grid is just for easy access to the elements by row,col,layer
  std::vector<SBSElement*> fElements;
  std::vector<SBSElement*> fRefElements; //< Reference elements (for TDCs and multi-function ADCs)
  std::vector<std::vector<std::vector<SBSElement*> > > fElementGrid;
  // Clusters for this event

  // Other parameters
  Bool_t fCoarseProcessed;  //< Was CourseProcessed already called?
  Bool_t fFineProcessed;    //< Was fine processed already called

  //Gain correction 
  Double_t   fConst;     // const from gain correction 
  Double_t   fSlope;     // slope for gain correction
  Double_t   fAccCharge; // accumulated charge

  // Per event data
  Int_t      fNhits;     ///< Number of hits in event
  Int_t      fNRefhits;     ///< Number of reference hits in event
  Int_t      fNGoodTDChits;     ///< Number of good TDC hits in event
  Int_t      fNGoodADChits;     ///< Number of good ADC hits in event

  // Flags for enabling and disabling various features
  Bool_t    fStoreRawHits; ///< Store the raw data in the root tree?

private:
  void ClearOutputVariables();

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
