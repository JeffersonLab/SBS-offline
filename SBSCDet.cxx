//////////////////////////////////////////////////////////////////////////
//
// SBSCDet class implementation
//
//////////////////////////////////////////////////////////////////////////

#include "SBSCDet.h"

ClassImp(SBSCDet);

/*
 * SBSCDet constructor.
 *
 * Use a TDC with trailing edge info, default is no ADC, but available for
 * commissioning only
 */
SBSCDet::SBSCDet( const char* name, const char* description,
    THaApparatus* apparatus ) : SBSGenericDetector(name,description,apparatus)
{
  SetModeTDC(SBSModeTDC::kTDC); //  A TDC with leading & trailing edge info
  SetModeADC(SBSModeADC::kNone); // Default is No ADC, but can be re-enabled later
}

///////////////////////////////////////////////////////////////////////////////
/// Read SBSCDet Database
Int_t SBSCDet::ReadDatabase( const TDatime& date )
{
  // We can use this name here for logs
  //static const char* const here = "ReadDatabase()";

  // If we want to add any new variables, uncomment the following and add
  // the new variables we want to read from the database
  FILE* file = OpenFile( date );
  if( !file ) return kFileError;
  //Int_t err;

  //std::cout<<"ReadDatabase method"<<std::endl;

  // Make sure to call parent class so that the generic variables can be read
  return SBSGenericDetector::ReadDatabase(date);
}

//_____________________________________________________________________________
Int_t SBSCDet::DefineVariables( EMode mode )
{
  // Initialize global variables
  Int_t err = SBSGenericDetector::DefineVariables(mode);
  if(err) {
    return err;
  }

  // Uncomment the following to add a ny new variables we want to define
  // as the output

  RVarDef vars[] = {
   //{ "nhits",       " number of PMT hits", "GetNumHits()"  },
   //{ "hit.pmtnum",  " Hit PMT num",        "fHits.SBSCDet.GetPMTNum()"},
   //{ "hit.xhit",    " PMT hit X",          "fHits.SBSCDet.GetX()"     },
   //{ "hit.yhit",    " PMT hit Y",          "fHits.SBSCDet.GetY()"     },
   //{ "hit.row",     " PMT hit row",        "fHits.SBSCDet.GetRow()"   },
   //{ "hit.col",     " PMT hit column",     "fHits.SBSCDet.GetCol()"   },
   //{ "hit.adc_r",   " PMT hit ADC right",  "fHits.SBSCDet.GetADC_r()" },
   //{ "hit.adc_l",   " PMT hit ADC left",   "fHits.SBSCDet.GetADC_l()" },
   //{ "hit.tdc_r",   " PMT hit TDC right",  "fHits.SBSCDet.GetTDC_r()" },
   //{ "hit.tdc_l",   " PMT hit TDC left",   "fHits.SBSCDet.GetTDC_l()" },
   { 0 }
  };
  err = DefineVarsFromList( vars, mode );
  
  // Finally go back
  return err;
}

//_____________________________________________________________________________
Int_t SBSCDet::Decode( const THaEvData& evdata )
{
  //std::cout << "SBSCDet::Decode" << std::endl;
  Int_t err = SBSGenericDetector::Decode(evdata);
  return err;
}
//


/*
 * Clear()
 * called at the end of every event
 */
void SBSCDet::Clear( Option_t* opt )
{
  // If we defined any new variables that we need to clear prior to the next event
  // clear them here:
  // fExample = 0.0;

  // Make sure to call parent class's Clear() also!
  SBSGenericDetector::Clear(opt);
}

/*
 * FindGoodHit()
 */
Int_t SBSCDet::FindGoodHit(SBSElement *)
{
  // Documentation for Hodoscope - not for CDet (yet)
  // The variable passed defines is one CDet paddle PMT
  // We can use it alone to find the good hits in that paddle or since
  // we know the row and column of that paddle, we can find it's corresponding
  // pair in fElementGrid[row][col][0]  (the last [0] is for the layer, we always
  // use only one layer

  // TODO: Implement logic here to determine good TDC Hit

  return 0;
}

Int_t SBSCDet::CoarseProcess( TClonesArray& tracks )
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

Int_t SBSCDet::FineProcess( TClonesArray& tracks )
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
 * Generic SBSCDet destructor
 */
SBSCDet::~SBSCDet()
{
  // Delete any new objects/instances created here
}
