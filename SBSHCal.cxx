///////////////////////////////////////////////////////////////////////////////
//
// SBSHCal
//
///////////////////////////////////////////////////////////////////////////////
#include "SBSHCal.h"
#include <iostream>
#include "THaEvData.h"

ClassImp(SBSHCal);

/*
 * SBSHCal constructor.
 *
 * Specify SBSCalorimeter to use both TDC and ADC Multi-samples
 */
SBSHCal::SBSHCal( const char* name, const char* description,
    THaApparatus* apparatus ) : SBSCalorimeter(name,description,apparatus)
{
  SetWithADCSamples(true);
  SetWithTDC(true);
  fWithLED = true;
}

///////////////////////////////////////////////////////////////////////////////
/// Read SBSHCal Database
Int_t SBSHCal::ReadDatabase( const TDatime& date )
{
  // We can use this name here for logs
  static const char* const here = "ReadDatabase()";

  if(fWithLED) {
    FILE* file = OpenFile( date );
    if( !file ) return kFileError;
  
    // Read in required geometry variables, which include fOrigin and fSize
    Int_t err = ReadGeometry( file, date, true );
    if( err ) {
      fclose(file);
      return err;
    }
  
    std::vector<Int_t> ledmap;
    DBRequest led_request[] = {
      { "ledmap", &ledmap, kIntV, 2, false }, ///< ledmap of LED
      {0}
    };
    err = LoadDB( file, date, led_request, fPrefix );
    if(err) {
	return err;
    }
    if(ledmap.size()<2) {
      Error(Here(here), "Need crate slot for LED");
      return kInitError;
    }
    fLEDCrate = ledmap[0];
    fLEDSlot  = ledmap[1];
  }

  return SBSCalorimeter::ReadDatabase(date);

}


//_____________________________________________________________________________
Int_t SBSHCal::Decode( const THaEvData& evdata )
{
  Int_t err = SBSCalorimeter::Decode(evdata);
  if(fWithLED) {
    Int_t ihit = evdata.GetNumChan(fLEDCrate,fLEDSlot);
    if(ihit!=2 ) {
      //std::cerr << "ihit=" << ihit << std::endl;
      return 0;
    }
    //assert(ihit==2);
    fLEDBit = evdata.GetData(fLEDCrate,fLEDSlot,1,0);
    fLEDCount = evdata.GetData(fLEDCrate,fLEDSlot,2,0);
    //std::cerr << "ihit LED: " << ihit << ", ledbit: " << fLEDBit << std::endl;
    if(fLEDBit==0) {
    SBSCalorimeterBlock *blk = 0;
    for(Int_t r = 0; r < fNrows; r++) {
      for(Int_t c = 0; c < fNcols; c++) { 
      blk = fBlocksGrid[r][c][0];
        if(blk->Samples()->GetDataSumRaw()>0&&false) {
          std::cerr << "[ " << (r+1) << ", " << (c+1) << "], Sum: " << blk->Samples()->GetDataSumRaw() << ", ADC:";
	  std::vector<Float_t> &s_r = blk->Samples()->GetDataRaw();
          size_t nsamples = s_r.size();
          for(size_t s = 0; s < nsamples; s++) {
            std::cerr << " " << s_r[s];
          }
          std::cerr << std::endl;
        }
      }
    }
    }
  }
  return err;
}


//_____________________________________________________________________________
Int_t SBSHCal::DefineVariables( EMode mode )
{
  // Initialize global variables
  Int_t err = SBSCalorimeter::DefineVariables(mode);
  if(err) {
    return err;
  }

  RVarDef vars[] = {
    { "ledbit", "LEDBit",  "fLEDBit" },
    { "ledcount", "LEDCount",  "fLEDCount" },
    { 0 }
  };

  err = DefineVarsFromList( vars, mode );
  return err;
}

void SBSHCal::ClearEvent()
{
  fLEDBit = -1;
  fLEDCount = 0;
  SBSCalorimeter::ClearEvent();
}
/*
 * Generic SBSHCal destructor
 */
SBSHCal::~SBSHCal()
{
}
