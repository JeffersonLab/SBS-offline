//*-- Author :    Ole Hansen   26 March 2001; completed by Guido Urciuoli;
//*-- last change: 9 January 2004.

//////////////////////////////////////////////////////////////////////////
//
// SBSCDet_Hit
//
//////////////////////////////////////////////////////////////////////////


#include "SBSCDet_Hit.h"
#include "TClonesArray.h"
#include "TMath.h"
#include <DataType.h>

#include <fstream>
#include <cstdio>
#include <cstring>
#include <cmath>

using namespace std;

//_____________________________________________________________________________
SBSCDet_Hit::SBSCDet_Hit(): 
  fPMTNum(-1), fRow(-1), fCol(-1), fLayer(-1),
  fX(kBig), fY(kBig), fZ(kBig), fTDC_LE(kBig), fTDC_TE(kBig), fToT(kBig)
  //fFlag(0), fVeto(0) 
{
  //fClustIndex = -1;
  fTrackIndex = -1;
} 

//_____________________________________________________________________________
SBSCDet_Hit::SBSCDet_Hit( Int_t pmtnum, Int_t i, Int_t j, Int_t k,
				    Double_t x, Double_t y, Double_t z, Double_t le, 
				    Double_t te, Double_t tot ):
  fPMTNum(pmtnum), fRow(i), fCol(j), fLayer(k),
  fX(x), fY(y), fZ(z), fTDC_LE(le), fTDC_TE(te), fToT(tot)
  //fFlag(0), fVeto(0), tdcr_set(false), tdcf_set(false) 
{
  //fClustIndex = -1;
  fTrackIndex = -1;
}

//_____________________________________________________________________________

Int_t SBSCDet_Hit::Compare( const TObject* theOtherHit ) const
{
  if (fTDC_LE < static_cast<const SBSCDet_Hit*>( theOtherHit )->fTDC_LE ||
      fTDC_TE < static_cast<const SBSCDet_Hit*>( theOtherHit )->fTDC_TE)
    return -1;
  if (fTDC_LE > static_cast<const SBSCDet_Hit*>( theOtherHit )->fTDC_LE ||
      fTDC_TE > static_cast<const SBSCDet_Hit*>( theOtherHit )->fTDC_TE)
    return +1;
  else
    return 0;
}

void SBSCDet_Hit::Clear( Option_t *opt ){ 
  fPMTNum = -1;
  fRow = -1;
  fCol = -1;
  fLayer = -1;
  //fClustIndex = -1;
  fTrackIndex = -1;
  fX = kBig;
  fY = kBig;
  fZ = kBig;
  //
  fTDC_LE = kBig;
  fTDC_TE = kBig;
  fToT = kBig;
}

/*
//_____________________________________________________________________________
void SBSCDet_Hit::Show(FILE * fout1, FILE* fout2) 
{
  // FIXME comment on fprintf(fout1...) should be changed accordingly 
  // when ntuple will be created 
  fprintf(fout2," PMT num ");
  fprintf(fout2,"%4d",fPMTNum);
  fprintf(fout2," Row, Col ");
  fprintf(fout2,"%4d,%4d",fRow,fCol);
  fprintf(fout2," ; X, Y ");
  fprintf(fout2,"%4f,%4f",fX,fY);
  fprintf(fout2,"; fTDC, fToT ");
  fprintf(fout2,"%4d,%4d",fTDC, fToT);
  fprintf(fout2,"; fADC = ");
  fprintf(fout2,"%4d",fADC);
  fprintf(fout2,"\n");
  //  Show(fout1);
}

void SBSCDet_Hit::Show(FILE * fout1) 
{
  // FIXME comment on fprintf(fout1...) should be changed accordingly 
  // when ntuple will be created
  fprintf(fout1," %4d",fPMTNum);
  fprintf(fout1," %4d %4d",fRow,fCol);
  fprintf(fout1,"% 4f %4f",fX,fY);
  fprintf(fout1," %4d %4d \n",fTime, fAmp);
  //fprintf(fout1," %4d \n",fADC);
}
*/

