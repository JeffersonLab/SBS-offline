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

  fHits             = new TClonesArray("SBSCDet_Hit",200);
  fHit_tmin 	    = -1000;
  fHit_tmax	    = 1000000; // wide open for now!

  Clear();
}

SBSCDet::~SBSCDet()
{
  // Destructor. Remove variables from global list and free up the memory
  // allocated by us.
  Clear();// so the prgram doesn't complain when deleting clusters
  RemoveVariables();
  delete fHits;
}


///////////////////////////////////////////////////////////////////////////////
/// Read SBSCDet Database
Int_t SBSCDet::ReadDatabase( const TDatime& date )
{
  // We can use this name here for logs
  //static const char* const here = "ReadDatabase()";

  // If we want to add any new variables, uncomment the following and add
  // the new variables we want to read from the database
  FILE* fi = OpenFile( date );
  if( !fi ) return kFileError;
  //Int_t err;

  std::cout<<"SBSCDet::ReadDatabase method"<<std::endl;
  
  Int_t err = SBSGenericDetector::ReadDatabase(date);
  if(err) {
    return err;
  }
  fIsInit = false;

  std::vector<Double_t> xpos,ypos,zpos;

  DBRequest config_request[] = {
    { "xpos", &xpos,    kDoubleV, 0, 1 },
    { "ypos", &ypos,    kDoubleV, 0, 1 },
    { "zpos", &zpos,    kDoubleV, 0, 1 },
    { 0 } ///< Request must end in a NULL
  };
  err = LoadDB( fi, date, config_request, fPrefix );

  if (!xpos.empty()) {
    if ((int)xpos.size() == fNelem) {
      for (Int_t ne=0;ne<fNelem;ne++) {
        fElements[ne]->SetX(xpos[ne]);
	//std::cout << "ne = " << ne << " xpos = " << xpos[ne] << std::endl;
      }
    } else {
      std::cout << "  vector too small " << xpos.size() << " # of elements =" << fNelem << std::endl;
    }
  }

  if (!ypos.empty()) {
    if ((int)ypos.size() == fNelem) {
      for (Int_t ne=0;ne<fNelem;ne++) {
        fElements[ne]->SetY(ypos[ne]);
	//std::cout << "ne = " << ne << " ypos = " << ypos[ne] << std::endl;
      }
    } else {
      std::cout << " ypos vector too small " << ypos.size() << " # of elements =" << fNelem << std::endl;
    }
  }
  
  if (!zpos.empty()) {
    if ((int)zpos.size() == fNelem) {
      for (Int_t ne=0;ne<fNelem;ne++) {
        fElements[ne]->SetZ(zpos[ne]);
	//std::cout << "ne = " << ne << " zpos = " << zpos[ne] << std::endl;
      }
    } else {
      std::cout << " zpos vector too small " << zpos.size() << " # of elements =" << fNelem << std::endl;
    }
  }

  fIsInit = true;

  fclose(fi);
  return kOK;

  // Make sure to call parent class so that the generic variables can be read
  //return SBSGenericDetector::ReadDatabase(date);

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
   { "ngoodhits",       " number of Good PMT hits", "GetNumHits()"  },
   { "hit.pmtnum",  " Hit PMT num",        "fHits.SBSCDet_Hit.GetPMTNum()"},
   { "hit.row",     " PMT hit row",        "fHits.SBSCDet_Hit.GetRow()"   },
   { "hit.col",     " PMT hit column",     "fHits.SBSCDet_Hit.GetCol()"   },
   { "hit.layer",     " PMT hit layer",     "fHits.SBSCDet_Hit.GetLayer()"   },
   { "hit.xhit",    " PMT hit X",          "fHits.SBSCDet_Hit.GetX()"     },
   { "hit.yhit",    " PMT hit Y",          "fHits.SBSCDet_Hit.GetY()"     },
   { "hit.zhit",    " PMT hit Z",          "fHits.SBSCDet_Hit.GetZ()"     },
   { "hit.tdc_le",   " PMT hit TDC LE",  "fHits.SBSCDet_Hit.GetTDC_LE()" },
   { "hit.tdc_te",   " PMT hit TDC TE",   "fHits.SBSCDet_Hit.GetTDC_TE()" },
   { "hit.tdc_tot",   " PMT hit TDC TOT",   "fHits.SBSCDet_Hit.GetToT()" },
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
  fHits->Clear("C");
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

  //std::cout << "SBSCDet::CoarseProcess ... " << std::endl;

  // Call the parent class so that it can prepare the data structure on the
  // event it just read from file

  //std::cout << "fNGoodTDChits = " << fNGoodTDChits << std::endl;
  SBSGenericDetector::CoarseProcess(tracks);
  //std::cout << "Back SBSGenericDetector::CoarseProcess ... fNGoodTDChits = " << fNGoodTDChits << std::endl;

  double x, y, z;
  //double tmin, tmax;

  Int_t nHit = 0;
  SBSCDet_Hit* the_hit = nullptr;

  //fNtrackMatch = 0;

  for(int k = 0; k<fNGoodTDChits; k++){
    //tmin = -fElements[fGood.TDCelemID[k]]->TDC()->GetGoodTimeCut();
    //tmax = +fElements[fGood.TDCelemID[k]]->TDC()->GetGoodTimeCut();

    //double t0 = fElements[fGood.TDCelemID[k]]->TDC()->GetGoodTimeCut();

    //    if(tmin<=fGood.t[k] && fGood.t[k]<=tmax){
    if (fGood.TDCelemID[k] >= 2688) {
	std::cout << "Processing good hit " << k << " fHit_tmin = " << fHit_tmin << " fHit_tmax = " << 
    	fHit_tmax << " le time = " << fGood.t[k] << " te time = " << fGood.t_te[k] << " PMT = " << fGood.TDCelemID[k] << std::endl; 
    }
    if( fHit_tmin <= fGood.t[k] && fGood.t[k] <= fHit_tmax ){
      the_hit = new( (*fHits)[nHit++] ) SBSCDet_Hit();

      the_hit->SetPMTNum(fGood.TDCelemID[k]);
      the_hit->SetRow(fGood.TDCrow[k]);
      the_hit->SetCol(fGood.TDCcol[k]);
      the_hit->SetLayer(fGood.TDClayer[k]);
      the_hit->SetTDC_LE(fGood.t[k]);
      the_hit->SetTDC_TE(fGood.t_te[k]);
      the_hit->SetToT(fGood.t_ToT[k]);

      x = (fElements[fGood.TDCelemID[k]])->GetX();
      y = (fElements[fGood.TDCelemID[k]])->GetY();
      z = (fElements[fGood.TDCelemID[k]])->GetZ();
      //std::cout << "X = " << x << " Y = " << y << " Z = " << z << std::endl;

      the_hit->SetX(x);
      the_hit->SetY(y);
      the_hit->SetZ(z);
    }
  }
  //clustering to be done by dereived class...




// EJB ---- old code below - March 29, 2025 --------------------------------
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

  //std::cout << "End Coarse Process "  << std::endl;
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



