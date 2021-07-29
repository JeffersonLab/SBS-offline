//////////////////////////////////////////////////////////////////////////
//
// SBSTimingHodoscope class implementation
//
//////////////////////////////////////////////////////////////////////////

#include "SBSTimingHodoscope.h"

ClassImp(SBSTimingHodoscope);

/*
 * SBSTimingHodoscope constructor.
 *
 * Use a TDC with trailing edge info, default is no ADC, but available for
 * commissioning only
 */
SBSTimingHodoscope::SBSTimingHodoscope( const char* name, const char* description,
    THaApparatus* apparatus ) : SBSGenericDetector(name,description,apparatus)
{
  SetModeTDC(SBSModeTDC::kTDC); //  A TDC with leading & trailing edge info
  SetModeADC(SBSModeADC::kNone); // Default is No ADC, but can be re-enabled later
}

///////////////////////////////////////////////////////////////////////////////
/// Read SBSTimingHodoscope Database
Int_t SBSTimingHodoscope::ReadDatabase( const TDatime& date )
{
  // We can use this name here for logs
  //static const char* const here = "ReadDatabase()";

  // If we want to add any new variables, uncomment the following and add
  // the new variables we want to read from the database
  //FILE* file = OpenFile( date );
  //if( !file ) return kFileError;
  //Int_t err;

  std::cout << "******** Detector " << GetName() << " ReadDatabase ********" << std::endl;

  // Make sure to call parent class so that the generic variables can be read
  return SBSGenericDetector::ReadDatabase(date);
}

//_____________________________________________________________________________
Int_t SBSTimingHodoscope::DefineVariables( EMode mode )
{
  // Initialize global variables
  Int_t err = SBSGenericDetector::DefineVariables(mode);
  if(err) {
    return err;
  }


  //RVarDef vars[] = {
  //  { "example", "Example variable",  "fExampleVariable" },
  //  { 0 }
  //};
  //err = DefineVarsFromList( vars, mode );
  
  // Finally go back
  return err;
}

/*
 * ClearEvent()
 * called at the end of every event
 */
void SBSTimingHodoscope::ClearEvent()
{
  // If we defined any new variables that we need to clear prior to the next event
  // clear them here:
  // fExample = 0.0;

  // Make sure to call parent class's ClearEvent() also!
  SBSGenericDetector::ClearEvent();
}

/*
 * FindGoodHit()
 */
Int_t SBSTimingHodoscope::FindGoodHit(SBSElement *blk)
{
  Int_t GoodHit=0;  
  if (blk->TDC()&& blk->HasData()) {
    blk->TDC()->SetGoodHit(0);
    GoodHit=1;
  }
  return GoodHit;
}

Int_t SBSTimingHodoscope::CoarseProcess( TClonesArray& tracks )
{
  if(fCoarseProcessed)
    return 0;

  // Call the parent class so that it can prepare the data structure on the
  // event it just read from file
  SBSGenericDetector::CoarseProcess(tracks);

  // All good hits now defined.  Now determine the position based on
  // time differences between paddles
  // For example:
  //
  //for(int row = 0; row < fNrows; row++) {
  //{
  //  SBSData::TDCHit lPMT = fElements[row][0][0]->TDC()->GetGoodHit();
  //  SBSData::TDCHit rPMT = fElements[row][1][0]->TDC()->GetGoodHit();
  //  diff = lPMT.le.val - rPMT.le.val; // Leading edge difference
  //}

  fCoarseProcessed = 1;
  return 0;
}

Int_t SBSTimingHodoscope::FineProcess( TClonesArray& tracks )
{

  if(fFineProcessed)
    return 0;

  // Do more detailed processing here.  Parent class does nothing, so no need
  // to call it.
  // We can prepare more detailed output if we want.

  fFineProcessed = 1;
  return 0;
}



/*
 * Generic SBSTimingHodoscope destructor
 */
SBSTimingHodoscope::~SBSTimingHodoscope()
{
  // Delete any new objects/instances created here
}
