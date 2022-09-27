#ifndef ROOT_TreeSearch_SBSEArm
#define ROOT_TreeSearch_SBSEArm

#include "THaSpectrometer.h"

class SBSEArm : public THaSpectrometer {

    public:
    SBSEArm( const char *name, const char *description );
    virtual ~SBSEArm();

    virtual Int_t FindVertices( TClonesArray& tracks );
    virtual Int_t TrackCalc();
    virtual Int_t CoarseReconstruct();
    virtual void  Clear( Option_t* opt="");

    protected:
    virtual Int_t ReadDatabase( const TDatime& date );
    virtual Int_t DefineVariables( EMode mode = kDefine );

    Double_t fFrontConstraintWidthX;
    Double_t fFrontConstraintWidthY;
    Double_t fBackConstraintWidthX;
    Double_t fBackConstraintWidthY;
    Double_t fFrontConstraintX0;
    Double_t fFrontConstraintY0;
    Double_t fBackConstraintX0; 
    Double_t fBackConstraintY0;


    std::vector<double> fFrontConstraintX;
  std::vector<double> fFrontConstraintY;
  std::vector<double> fFrontConstraintZ;
  std::vector<double> fBackConstraintX;
  std::vector<double> fBackConstraintY;
  std::vector<double> fBackConstraintZ;

    TVector3 fGEMorigin;  //Absolute position of GEM origin relative to target center, in TARGET transport coordinates



    ClassDef(SBSEArm,0) // BigBite spectrometer

};


#endif//ROOT_TreeSearch_SBSEArm

