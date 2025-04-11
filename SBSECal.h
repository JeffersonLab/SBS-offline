#ifndef SBSECal_h
#define SBSECal_h

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// SBSECal                                                                   //
//                                                                           //
// A derived class of SBSCalorimeter that has Multi-valued ADC + TDC             //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "SBSCalorimeter.h"

class SBSECal : public SBSCalorimeter {
public:
  explicit SBSECal( const char* name, const char* description = "",
                    THaApparatus* a = nullptr);
  virtual ~SBSECal();
  virtual Int_t  ReadDatabase( const TDatime& date );
  virtual Int_t  Decode( const THaEvData& evdata );
  virtual Int_t  CoarseProcess(TClonesArray& tracks );
  //virtual Int_t  FineProcess(TClonesArray &tracks );
  virtual Int_t  DefineVariables( EMode mode = kDefine );
  virtual void   Clear( Option_t* opt="" );

  virtual Int_t SelectBestCluster();

  inline void SetRefADCtimeGoodCluster(Double_t tval){ fRefADCtimeGoodCluster = tval; }

  void CalcShowerCoord(); //Calculate shower coordinates
  
protected:
  // Bool_t fWithLED;
  // Int_t fLEDCrate;
  // Int_t fLEDSlot;
  // Int_t fLEDBit;
  // Int_t fLEDCount;

  //ADC time window for "good" cluster selection
  Double_t fRefADCtimeGoodCluster; // default to zero. This can be set based on the electron arm information using an InterStage module (forthcoming)
  Double_t fAtimeMinGoodCluster; //default to wide open limits t0 +/- width
  Double_t fAtimeMaxGoodCluster; //default to wide optn limits 

  // //Flag to filter clusters based on presence of at least one good TDC hit?
  // Bool_t fRequireTDCGoodCluster;

  //We'll need to add some output variables
  //Actually no we won't. 
  // We just need a public method to calculate and apply the shower shape corrections
  // to the cluster position.
  
  // Adding parameters for shower shape corrections:

  Bool_t fUseShowerShapeCorr; //Flag to correct ECAL cluster positions using "shower shape" profiles

  Int_t fShowerShapeNbinsX;
  Int_t fShowerShapeNbinsY;
  Int_t fShowerShapeNbinsProf;
  Double_t fShowerShapeMomMin; //default to -1
  Double_t fShowerShapeMomMax; //default to +1

  Double_t fShowerShapeXminProf;
  Double_t fShowerShapeXmaxProf;
  Double_t fShowerShapeYminProf;
  Double_t fShowerShapeYmaxProf;
  
  Double_t fShowerShapeLxProf;
  Double_t fShowerShapeLyProf;

  std::vector<std::vector<Double_t> > fShowerProfilesX;
  std::vector<std::vector<Double_t> > fShowerProfilesY;

  
  ClassDef(SBSECal,0)     // ECal detector class
};
#endif // SBSECal_h
