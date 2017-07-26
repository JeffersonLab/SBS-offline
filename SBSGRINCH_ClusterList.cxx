//*-- Author :    Ole Hansen   26 March 2001; completed by Guido Urciuoli;
//*-- last change: 9 January 2004.

//////////////////////////////////////////////////////////////////////////
//
// SBSGRINCH_ClusterList
//
//////////////////////////////////////////////////////////////////////////


#include "SBSGRINCH_ClusterList.h"
#include "TClonesArray.h"
#include "TMath.h"

#include <fstream>
#include <cstdio>
#include <cstring>
#include <cmath>

using namespace std;

//_____________________________________________________________________________
Int_t SBSGRINCH_Hit::Compare( const TObject* theOtherHit ) const
{
  if (fADC < static_cast<const SBSGRINCH_Hit*>( theOtherHit )->fADC)
    return -1;
  if (fADC > static_cast<const SBSGRINCH_Hit*>( theOtherHit )->fADC)
    return +1;
  else
    return 0;
}

//_____________________________________________________________________________
void SBSGRINCH_Hit::Show(FILE * fout1, FILE* fout2) 
{
  // FIXME comment on fprintf(fout1...) should be changed accordingly 
  // when ntuple will be created 
  fprintf(fout2," I, J ");
  fprintf(fout2,"%4d,%4d",fI,fJ);
  fprintf(fout2," ; X, Y ");
  fprintf(fout2,"%4f,%4f",fX,fY);
  fprintf(fout2,"; fADC = ");
  fprintf(fout2,"%4d",fADC);
  fprintf(fout2,"\n");
  //  Show(fout1);
}

void SBSGRINCH_Hit::Show(FILE * fout1) 
{
  // FIXME comment on fprintf(fout1...) should be changed accordingly 
  // when ntuple will be created
  fprintf(fout1," %4d %4d",fI,fJ);
  fprintf(fout1,"% 4f %4f",fX,fY);
  fprintf(fout1," %4d \n",fADC);
}


//=============================================================================
SBSGRINCH_Cluster::SBSGRINCH_Cluster() :
  fLocalMaximumNumber(1), fMIPflag(kFALSE), fFictious_Mip_Flag(0),
  fPionChi2AnalysisFlag(kFALSE), fKaonChi2AnalysisFlag(kFALSE), 
  fProtonChi2AnalysisFlag(kFALSE), fXcenter(0.0), fYcenter(0.0), fCharge(0.0),
  fTheta_photon(0.0), fPhi_photon(0.0), fAngle(0.0), fMIP(0), fTrack(0), 
  fnoise_cut_success(0), fResolved_noise_cut_success(0) 
{
  for(int i=0; i<3; i++) {
    N_Photon[i] = 0; 
    angle[i] = 0.0;
    angle_corrected[i] = 0.0;
    ResolvedN_Photon[i]= 0;
    Resolvedangle[i] = 0.0;
    Resolvedangle_corrected[i] = 0.0;
    N_chi2_Photon[i] = 0;
    chi2[i] = 0.0;
    chi2_prob[i] = 0.0;
    N_chi2_corrected_Photon[i] = 0;
    chi2_corrected[i]    = 0.0;
    chi2_corrected_prob[i] = 0.0;
    N_MaximumLikelihood_Photon[i] = 0;
    MaximumLikelihood[i] = 0.0;
    ResolvedN_chi2_Photon[i] = 0;
    Resolvedchi2[i] = 0.0;
    Resolvedchi2_prob[i] = 0.0;
    ResolvedN_chi2_corrected_Photon[i] = 0;
    Resolvedchi2_corrected[i] = 0.;
    Resolvedchi2_corrected_prob[i] = 0.;
    ResolvedN_MaximumLikelihood_Photon[i] = 0;
    ResolvedMaximumLikelihood[i] = 0.0;
  }
  fHitList = new TList; 
}

//_____________________________________________________________________________
SBSGRINCH_Cluster::SBSGRINCH_Cluster( const SBSGRINCH_Cluster& rhs ) :
  fLocalMaximumNumber(rhs.fLocalMaximumNumber), fMIPflag(rhs.fMIPflag),
  fFictious_Mip_Flag(rhs.fFictious_Mip_Flag),
  fPionChi2AnalysisFlag(rhs.fPionChi2AnalysisFlag),
  fKaonChi2AnalysisFlag(rhs.fKaonChi2AnalysisFlag),
  fProtonChi2AnalysisFlag(rhs.fProtonChi2AnalysisFlag),
  fXcenter(rhs.fXcenter), fYcenter(rhs.fYcenter), fCharge(rhs.fCharge),
  fTheta_photon(rhs.fTheta_photon), fPhi_photon(rhs.fPhi_photon),
  fAngle(rhs.fAngle), fMIP(rhs.fMIP), fTrack(rhs.fTrack),
  fnoise_cut_success(rhs.fnoise_cut_success), 
  fResolved_noise_cut_success(rhs.fResolved_noise_cut_success)
{
  // Copy constructor
  for(int i=0; i<3; i++) {
    N_Photon[i] = rhs.N_Photon[i]; 
    angle[i]    = rhs.angle[i];
    angle_corrected[i]    = rhs.angle_corrected[i];
    N_chi2_Photon[i] = rhs.N_chi2_Photon[i]; 
    chi2[i]    = rhs.chi2[i];
    chi2_prob[i]    = rhs.chi2_prob[i];
    N_chi2_corrected_Photon[i] = rhs.N_chi2_corrected_Photon[i];
    chi2_corrected[i]    = rhs.chi2_corrected[i]; 
    chi2_corrected_prob[i]    = rhs.chi2_corrected_prob[i];
    N_MaximumLikelihood_Photon[i] =  rhs.N_MaximumLikelihood_Photon[i];
    MaximumLikelihood[i] =  rhs.MaximumLikelihood[i];
    ResolvedN_Photon[i] = rhs.ResolvedN_Photon[i];
    Resolvedangle[i] = rhs.Resolvedangle[i];
    Resolvedangle_corrected[i] = rhs.Resolvedangle_corrected[i];
    ResolvedN_chi2_Photon[i] = rhs.ResolvedN_chi2_Photon[i];
    Resolvedchi2[i] = rhs.Resolvedchi2[i];
    Resolvedchi2_prob[i] = rhs.Resolvedchi2_prob[i];
    ResolvedN_chi2_corrected_Photon[i] =rhs.ResolvedN_chi2_corrected_Photon[i];
    Resolvedchi2_corrected[i] = rhs.Resolvedchi2_corrected[i];
    Resolvedchi2_corrected_prob[i] = rhs.Resolvedchi2_corrected_prob[i];
    ResolvedN_MaximumLikelihood_Photon[i] =  
      rhs.ResolvedN_MaximumLikelihood_Photon[i];
    ResolvedMaximumLikelihood[i] = rhs.ResolvedMaximumLikelihood[i];
  }
  fHitList = new TList; 
  if( rhs.fHitList && (rhs.fHitList->GetSize() > 0 )) {
    TIter next( rhs.fHitList );
    while( SBSGRINCH_Hit* pHit = static_cast<SBSGRINCH_Hit*>( next() ))
      fHitList->AddLast(pHit);
  }
}

//_____________________________________________________________________________
SBSGRINCH_Cluster& SBSGRINCH_Cluster::operator=( const SBSGRINCH_Cluster& rhs )
{
  // Assignment operator

  if( this != &rhs ) {
    fLocalMaximumNumber = rhs.fLocalMaximumNumber;
    fMIPflag = rhs.fMIPflag;
    fPionChi2AnalysisFlag = rhs.fPionChi2AnalysisFlag;
    fKaonChi2AnalysisFlag = rhs.fKaonChi2AnalysisFlag;
    fProtonChi2AnalysisFlag = rhs.fProtonChi2AnalysisFlag;
    fXcenter = rhs.fXcenter;
    fYcenter = rhs.fYcenter;
    fCharge = rhs.fCharge;
    fTheta_photon = rhs.fTheta_photon;
    fPhi_photon = rhs.fPhi_photon;
    fAngle = rhs.fAngle;
    fMIP = rhs.fMIP;
    fTrack = rhs.fTrack;
    fnoise_cut_success = rhs.fnoise_cut_success;
    fResolved_noise_cut_success = rhs.fResolved_noise_cut_success;
    fFictious_Mip_Flag = rhs.fFictious_Mip_Flag;
    for(int i=0; i<3; i++) {
      N_Photon[i] = rhs.N_Photon[i]; 
      angle[i]    = rhs.angle[i];
      angle_corrected[i]    = rhs.angle_corrected[i];
      N_chi2_Photon[i] = rhs.N_chi2_Photon[i]; 
      chi2[i]    = rhs.chi2[i];
      chi2_prob[i]    = rhs.chi2_prob[i];
      N_chi2_corrected_Photon[i] = rhs.N_chi2_corrected_Photon[i]; 
      chi2_corrected[i]    = rhs.chi2_corrected[i];
      chi2_corrected_prob[i]    = rhs.chi2_corrected_prob[i];
      N_MaximumLikelihood_Photon[i] =  rhs.N_MaximumLikelihood_Photon[i];
      MaximumLikelihood[i] = rhs.MaximumLikelihood[i];
      ResolvedN_Photon[i] = rhs.ResolvedN_Photon[i];
      Resolvedangle[i] =  rhs.Resolvedangle[i];
      Resolvedangle_corrected[i] =  rhs.Resolvedangle_corrected[i];
      ResolvedN_chi2_Photon[i] = rhs.ResolvedN_chi2_Photon[i];
      Resolvedchi2[i] = rhs.Resolvedchi2[i];
      Resolvedchi2_prob[i] = rhs.Resolvedchi2_prob[i];
      ResolvedN_chi2_corrected_Photon[i]=rhs.ResolvedN_chi2_corrected_Photon[i];
      Resolvedchi2_corrected[i] = rhs.Resolvedchi2_corrected[i];
      Resolvedchi2_corrected_prob[i] = rhs.Resolvedchi2_corrected_prob[i];
      ResolvedN_MaximumLikelihood_Photon[i] =  
	rhs.ResolvedN_MaximumLikelihood_Photon[i];
      ResolvedMaximumLikelihood[i] = rhs.ResolvedMaximumLikelihood[i];
    }
    if( !fHitList )
      fHitList = new TList;
    else
      fHitList->Clear("nodelete");
    if( rhs.fHitList && (rhs.fHitList->GetSize() > 0 )) {
      TIter next( rhs.fHitList ); 
      while( SBSGRINCH_Hit* pHit = static_cast<SBSGRINCH_Hit*>( next() ))
	fHitList->AddLast(pHit);
    }
 }
  return *this;
}

//_____________________________________________________________________________
void SBSGRINCH_Cluster::Clear( Option_t* opt )
{
  //Reset the cluster to an empty state.

  //Don't delete the hits, just clear the internal list.
  if( fHitList ) 
    fHitList->Clear("nodelete");

  // Full clear
  if( opt && opt[0] == 'F' ) {
    fMIPflag = kFALSE;
    fPionChi2AnalysisFlag = kFALSE;
    fKaonChi2AnalysisFlag = kFALSE;
    fProtonChi2AnalysisFlag = kFALSE;
    fXcenter = 0.0;
    fYcenter = 0.0;
    fCharge  = 0.0;
    fAngle   = 0.0;
    fTrack   = 0;
    fMIP     = 0;
    fTheta_photon = 0.;
    fPhi_photon = 0;
    fnoise_cut_success = 0;
    fResolved_noise_cut_success = 0;
    fFictious_Mip_Flag = 0;
    for (Int_t i = 0; i < 3; i++) { 
      N_Photon[i] = 0;
      angle[i] = 0.0;
      angle_corrected[i] = 0.0;
      N_chi2_Photon[i] = 0;
      chi2[i]    = 0.0;
      chi2_prob[i]    = 0.0;
      N_chi2_corrected_Photon[i] = 0;
      chi2_corrected[i]    = 0.0;
      chi2_corrected_prob[i]    = 0.0;
      N_MaximumLikelihood_Photon[i] = 0;
      MaximumLikelihood[i] = 0.0;
      ResolvedN_Photon[i] = 0;
      Resolvedangle[i] = 0.0;
      Resolvedangle_corrected[i] = 0.0;
      ResolvedN_chi2_Photon[i] = 0;
      Resolvedchi2[i] = 0.0;
      Resolvedchi2_prob[i] = 0.0;
      ResolvedN_chi2_corrected_Photon[i] = 0;
      Resolvedchi2_corrected[i] = 0.;
      Resolvedchi2_corrected_prob[i] = 0.;
      N_MaximumLikelihood_Photon[i] = 0;
      ResolvedMaximumLikelihood[i] = 0.0;
    }
  } else {
    // Fast clear for clearing TClonesArrays of clusters
    // Needs to be deleted since a TList allocates memory
    // FIXME: performance issue?
    delete fHitList;
    fHitList = 0;
  }  
}

//_____________________________________________________________________________
Float_t SBSGRINCH_Cluster::Dist( const SBSGRINCH_Cluster* c ) const
{
  // Calculate distance of this cluster to another cluster

  Float_t dx = GetXcenter() - c->GetXcenter();
  Float_t dy = GetYcenter() - c->GetYcenter();
  return sqrt( dx*dx + dy*dy );
}

//_____________________________________________________________________________
void SBSGRINCH_Cluster::Insert( SBSGRINCH_Hit* theHit )
{
  //Add a hit to the cluster

  if( !fHitList ) fHitList = new TList;
  fHitList->AddLast( theHit );
  fXcenter = fXcenter*fCharge;
  fYcenter = fYcenter*fCharge;
  fCharge  += theHit->GetADC();
  fXcenter += theHit->GetADC()*theHit->GetX();
  fYcenter += theHit->GetADC()*theHit->GetY();
  fXcenter = fXcenter/fCharge;
  fYcenter = fYcenter/fCharge;
}

//_____________________________________________________________________________
void SBSGRINCH_Cluster::Insert( SBSGRINCH_Hit* theHit, Float_t factor )
{
  //Add a hit to the cluster with its charge weighted by a factor.

  if( !fHitList ) fHitList = new TList;
  fHitList->AddLast( theHit );
  fXcenter = fXcenter*fCharge;
  fYcenter = fYcenter*fCharge;
  fCharge  += theHit->GetADC()*factor;
  fXcenter += theHit->GetADC()*factor*theHit->GetX();
  fYcenter += theHit->GetADC()*factor*theHit->GetY();
  fXcenter = fXcenter/fCharge;
  fYcenter = fYcenter/fCharge;
}

//_____________________________________________________________________________
void SBSGRINCH_Cluster::Insert_Photon(Int_t flag, Float_t angle_to_add,
				    Int_t ResolvedFlag)
{
  // add a photon to the MIP

  if(flag == 0) {
    angle[flag] = angle[flag]*N_Photon[flag];
    angle_corrected[flag] = angle_corrected[flag]*N_Photon[flag];
    N_Photon[flag]++;
    angle[flag] = angle[flag] + angle_to_add;
    angle_corrected[flag] = angle_corrected[flag] + angle_to_add;
    angle[flag] = angle[flag]/N_Photon[flag];
    angle_corrected[flag] = angle_corrected[flag]/N_Photon[flag];

  } else {
    Resolvedangle[flag] = Resolvedangle[flag]*ResolvedN_Photon[flag];
    Resolvedangle_corrected[flag] = 
      Resolvedangle_corrected[flag]*ResolvedN_Photon[flag];
    ResolvedN_Photon[flag]++;
    Resolvedangle_corrected[flag] = Resolvedangle_corrected[flag] + 
      angle_to_add;
    Resolvedangle[flag] = Resolvedangle[flag]/ResolvedN_Photon[flag];
    Resolvedangle_corrected[flag] = 
      Resolvedangle_corrected[flag]/ResolvedN_Photon[flag];
  }
}
//_____________________________________________________________________________
void SBSGRINCH_Cluster::Insert_Photon(Int_t flag, Float_t angle_to_add,
				    Int_t ResolvedFlag, Float_t expected_angle,
				    Float_t central_momentum_expected_angle)
{
  // add a photon to the MIP

  //cout << "ResolvedFlag " << ResolvedFlag << endl;
  if(ResolvedFlag == 0) {
    angle[flag] = angle[flag]*N_Photon[flag];
    angle_corrected[flag] = angle_corrected[flag]*N_Photon[flag];
    N_Photon[flag]++;
    angle[flag] = angle[flag] + angle_to_add;
    angle_corrected[flag] = angle_corrected[flag] + angle_to_add 
      - expected_angle + central_momentum_expected_angle;
    angle[flag] = angle[flag]/N_Photon[flag];
    angle_corrected[flag] = angle_corrected[flag]/N_Photon[flag];

  } else {
    Resolvedangle[flag] = Resolvedangle[flag]*ResolvedN_Photon[flag];
    Resolvedangle_corrected[flag] = 
      Resolvedangle_corrected[flag]*ResolvedN_Photon[flag];
    ResolvedN_Photon[flag]++;
    Resolvedangle[flag] = Resolvedangle[flag] + angle_to_add;
    Resolvedangle_corrected[flag] = Resolvedangle_corrected[flag] + 
      angle_to_add - expected_angle + central_momentum_expected_angle;
    Resolvedangle[flag] = Resolvedangle[flag]/ResolvedN_Photon[flag];
    Resolvedangle_corrected[flag] = 
      Resolvedangle_corrected[flag]/ResolvedN_Photon[flag];
    //cout << "Resolvedangle[flag]" << Resolvedangle[flag] << endl;
  }
}
//____________________________________________________________________________
void SBSGRINCH_Cluster::Insert_chi2(Int_t flag, Float_t chi2_to_add)
{
  // only for simple clusters, it calculates the single contribution 
  // of a cluster to the chi square value
    N_chi2_Photon[flag]++; 
    chi2[flag] = chi2[flag] + chi2_to_add; 
} 
//____________________________________________________________________________
void SBSGRINCH_Cluster::Insert_chi2(Int_t flag, Float_t chi2_to_add,
				    Int_t ResolvedFlag)
{
  // add a photon to the MIP for chi square calculation 
  if(ResolvedFlag == 0) {
    N_chi2_Photon[flag]++; 
    chi2[flag] = chi2[flag] + chi2_to_add; 
  } else {
    ResolvedN_chi2_Photon[flag]++;
    Resolvedchi2[flag] = Resolvedchi2[flag] + chi2_to_add;
    Resolvedchi2_prob[flag] = TMath::Prob(Resolvedchi2[flag],
					  ResolvedN_chi2_Photon[flag]);
  }
}
//____________________________________________________________________________
void SBSGRINCH_Cluster::Insert_chi2_corrected(Int_t flag, Float_t chi2,
				    Int_t ResolvedFlag)
{
  // add a photon to the MIP

  if(ResolvedFlag == 0) { 
    chi2_corrected[flag] = chi2;
  } else {
    Resolvedchi2_corrected[flag] = chi2;
  }
}
//____________________________________________________________________________
void SBSGRINCH_Cluster::Insert_N_chi2_corrected_Photon(Int_t flag, Int_t N_Photon, Int_t ResolvedFlag)
{
  // add a photon to the MIP

  if(ResolvedFlag == 0) { 
    N_chi2_corrected_Photon[flag] = N_Photon; 
  } else {
    ResolvedN_chi2_corrected_Photon[flag] = N_Photon;
  }
}
//____________________________________________________________________________
void SBSGRINCH_Cluster::Insert_MaximumLikelihood(Int_t flag, 
					     Float_t MaximumLikelihood_to_add,
					     Int_t ResolvedFlag)
{
  // add a photon to the MIP for Maximum Likelihood calculation
  if(ResolvedFlag == 0) {
    N_MaximumLikelihood_Photon[flag]++; 
    MaximumLikelihood[flag] = MaximumLikelihood[flag] + 
      MaximumLikelihood_to_add; 
  } else {
    ResolvedN_MaximumLikelihood_Photon[flag]++;
    ResolvedMaximumLikelihood[flag] = ResolvedMaximumLikelihood[flag] + 
      MaximumLikelihood_to_add;
  }
}
//____________________________________________________________________________
  void SBSGRINCH_Cluster::Setchi2_prob(Int_t flag, Float_t chi2_value, 
				       Int_t N_Photon, Int_t ResolvedFlag)
{
  // calculate chi2 probability

  if(ResolvedFlag == 0) { 
    chi2_prob[flag] = TMath::Prob(chi2_value, N_Photon); 
  } else {
    Resolvedchi2_prob[flag] = TMath::Prob(chi2_value, N_Photon);
  }
}
//____________________________________________________________________________
  void SBSGRINCH_Cluster::Setchi2_corrected_prob(Int_t flag, Float_t chi2_value,
					       Int_t N_Photon, 
					       Int_t ResolvedFlag)
{
  // calculate chi2 probability

  if(ResolvedFlag == 0) { 
    chi2_corrected_prob[flag] = TMath::Prob(chi2_value, N_Photon); 
  } else {
    Resolvedchi2_corrected_prob[flag] = TMath::Prob(chi2_value, N_Photon);
  }
}
//_____________________________________________________________________________
void SBSGRINCH_Cluster::Setnoise_cut_success(Int_t value, Int_t ResolvedFlag)
{
  if(ResolvedFlag == 0)
    {
      fnoise_cut_success = value;
    }
  else
    {
      fResolved_noise_cut_success = value;
    }
}
//_____________________________________________________________________________
void SBSGRINCH_Cluster::SetN_Photon(Int_t i, Int_t Value, Int_t ResolvedFlag)
{
  if(ResolvedFlag == 0)
    {
      N_Photon[i] = Value;
    }
  else
    {
      ResolvedN_Photon[i] = Value;
    }
}
//_____________________________________________________________________________
void SBSGRINCH_Cluster::Setangle(Int_t i, Float_t Value, Int_t ResolvedFlag)
{
  if(ResolvedFlag == 0)
    {
      angle[i] = Value;
    }
  else
    {
      Resolvedangle[i] = Value;
    }
}
//_____________________________________________________________________________
void SBSGRINCH_Cluster::Setangle_corrected(Int_t i, Float_t Value, 
					 Int_t ResolvedFlag)
{
  if(ResolvedFlag == 0)
    {
      angle_corrected[i] = Value;
    }
  else
    {
      Resolvedangle_corrected[i] = Value;
    }
}
//_____________________________________________________________________________
Int_t SBSGRINCH_Cluster::GetN_Photon(Int_t i, Int_t ResolvedFlag) const     
{
  if(ResolvedFlag == 0)
    {
      return N_Photon[i];
    }
  else
    {
      return ResolvedN_Photon[i];
    }
}
//_____________________________________________________________________________
Int_t SBSGRINCH_Cluster::GetNphot_pi(Int_t ResolvedFlag) const             
{
  if(ResolvedFlag == 0)
    {
      return N_Photon[0];
    }
  else
    {
      return ResolvedN_Photon[0];
    }
}
//_____________________________________________________________________________
Int_t SBSGRINCH_Cluster::GetNphot_k(Int_t ResolvedFlag) const             
{
  if(ResolvedFlag == 0)
    {
      return N_Photon[1];
    }
  else
    {
      return ResolvedN_Photon[1];
    }
}
//_____________________________________________________________________________
Int_t SBSGRINCH_Cluster::GetNphot_p(Int_t ResolvedFlag) const             
{
  if(ResolvedFlag == 0)
    {
      return N_Photon[2];
    }
  else
    {
      return ResolvedN_Photon[2];
    }
}
//_____________________________________________________________________________
Float_t SBSGRINCH_Cluster::Getangle(Int_t i, Int_t ResolvedFlag) const       
{
  if(ResolvedFlag == 0)
    {
      return angle[i];
    }
  else
    {
      return Resolvedangle[i];
    }
}
//_____________________________________________________________________________
Float_t SBSGRINCH_Cluster::Getangle_pi(Int_t ResolvedFlag) const             
{
  if(ResolvedFlag == 0)
    {
      return angle[0];
    }
  else
    {
      return Resolvedangle[0];
    }
}
//_____________________________________________________________________________
Float_t SBSGRINCH_Cluster::Getangle_k(Int_t ResolvedFlag) const             
{
  if(ResolvedFlag == 0)
    {
      return angle[1];
    }
  else
    {
      return Resolvedangle[1];
    }
}
//_____________________________________________________________________________
Float_t SBSGRINCH_Cluster::Getangle_p(Int_t ResolvedFlag) const             
{
  if(ResolvedFlag == 0)
    {
      return angle[2];
    }
  else
    {
      return Resolvedangle[2];
    }
}
//_____________________________________________________________________________
Float_t SBSGRINCH_Cluster::Getangle_corrected(Int_t i, Int_t ResolvedFlag) const       
{
  if(ResolvedFlag == 0)
    {
      return angle_corrected[i];
    }
  else
    {
      return Resolvedangle_corrected[i];
    }
}
//_____________________________________________________________________________
Float_t SBSGRINCH_Cluster::Getangle_corrected_pi(Int_t ResolvedFlag) const             
{
  if(ResolvedFlag == 0)
    {
      return angle_corrected[0];
    }
  else
    {
      return Resolvedangle_corrected[0];
    }
}
//_____________________________________________________________________________
Float_t SBSGRINCH_Cluster::Getangle_corrected_k(Int_t ResolvedFlag) const             
{
  if(ResolvedFlag == 0)
    {
      return angle_corrected[1];
    }
  else
    {
      return Resolvedangle_corrected[1];
    }
}
//_____________________________________________________________________________
Float_t SBSGRINCH_Cluster::Getangle_corrected_p(Int_t ResolvedFlag) const             
{
  if(ResolvedFlag == 0)
    {
      return angle_corrected[2];
    }
  else
    {
      return Resolvedangle_corrected[2];
    }
}
//____________________________________________________________________________
 Int_t SBSGRINCH_Cluster::GetN_chi2_phot_pi(Int_t ResolvedFlag) const 
{
  if(ResolvedFlag == 0)
    {
      return N_chi2_Photon[0];
    }
  else
    {
      return ResolvedN_chi2_Photon[0];
    }
}
//____________________________________________________________________________
 Int_t SBSGRINCH_Cluster::GetN_chi2_phot_k(Int_t ResolvedFlag) const 
{
  if(ResolvedFlag == 0)
    {
      return N_chi2_Photon[1];
    }
  else
    {
      return ResolvedN_chi2_Photon[1];
    }
}
//____________________________________________________________________________
 Int_t SBSGRINCH_Cluster::GetN_chi2_phot_p(Int_t ResolvedFlag) const 
{
  if(ResolvedFlag == 0)
    {
      return N_chi2_Photon[2];
    }
  else
    {
      return ResolvedN_chi2_Photon[2];
    }
}
//____________________________________________________________________________
Float_t    SBSGRINCH_Cluster::Getchi2_pi(Int_t ResolvedFlag) const
{
  if(ResolvedFlag == 0)
    {
      return chi2[0];
    }
  else
    {
      return Resolvedchi2[0];
    }
}
//____________________________________________________________________________
Float_t    SBSGRINCH_Cluster::Getchi2_k(Int_t ResolvedFlag) const
{
  if(ResolvedFlag == 0)
    {
      return chi2[1];
    }
  else
    {
      return Resolvedchi2[1];
    }
}
//____________________________________________________________________________
Float_t    SBSGRINCH_Cluster::Getchi2_p(Int_t ResolvedFlag) const
{
  if(ResolvedFlag == 0)
    {
      return chi2[2];
    }
  else
    {
      return Resolvedchi2[2];
    }
}
//____________________________________________________________________________
Float_t    SBSGRINCH_Cluster::Getchi2_prob_pi(Int_t ResolvedFlag) const
{
  if(ResolvedFlag == 0)
    {
      return chi2_prob[0];
    }
  else
    {
      return Resolvedchi2_prob[0];
    }
}
//____________________________________________________________________________
Float_t    SBSGRINCH_Cluster::Getchi2_prob_k(Int_t ResolvedFlag) const
{
  if(ResolvedFlag == 0)
    {
      return chi2_prob[1];
    }
  else
    {
      return Resolvedchi2_prob[1];
    }
}
//____________________________________________________________________________
Float_t    SBSGRINCH_Cluster::Getchi2_prob_p(Int_t ResolvedFlag) const
{
  if(ResolvedFlag == 0)
    {
      return chi2_prob[2];
    }
  else
    {
      return Resolvedchi2_prob[2];
    }
}
//____________________________________________________________________________
Float_t    SBSGRINCH_Cluster::Getchi2_prob(Int_t flag, Int_t ResolvedFlag) const
{
  if(ResolvedFlag == 0)
    {
      return chi2_prob[flag];
    }
  else
    {
      return Resolvedchi2_prob[flag];
    }
}
//____________________________________________________________________________
Int_t SBSGRINCH_Cluster::GetN_chi2_corrected_phot_pi(Int_t ResolvedFlag) const 
{
  if(ResolvedFlag == 0)
    {
      return N_chi2_corrected_Photon[0];
    }
  else
    {
       return ResolvedN_chi2_corrected_Photon[0];
    }
}
//____________________________________________________________________________
Int_t SBSGRINCH_Cluster::GetN_chi2_corrected_phot_k(Int_t ResolvedFlag) const 
{
  if(ResolvedFlag == 0)
    {
      return N_chi2_corrected_Photon[1];
    }
  else
    {
       return ResolvedN_chi2_corrected_Photon[1];
    }
}
//____________________________________________________________________________
Int_t SBSGRINCH_Cluster::GetN_chi2_corrected_phot_p(Int_t ResolvedFlag) const
{
  if(ResolvedFlag == 0)
    {
      return N_chi2_corrected_Photon[2];
    }
  else
    {
       return ResolvedN_chi2_corrected_Photon[2];
    }
}
//____________________________________________________________________________
Float_t SBSGRINCH_Cluster::Getchi2_corrected_pi(Int_t ResolvedFlag) const
{
  if(ResolvedFlag == 0)
    {
      return chi2_corrected[0];
    }
  else
    {
       return Resolvedchi2_corrected[0];
    }
}
//____________________________________________________________________________
Float_t SBSGRINCH_Cluster::Getchi2_corrected_k(Int_t ResolvedFlag) const 
{
  if(ResolvedFlag == 0)
    {
      return chi2_corrected[1];
    }
  else
    {
       return Resolvedchi2_corrected[1];
    }
}
//____________________________________________________________________________
Float_t SBSGRINCH_Cluster::Getchi2_corrected_p(Int_t ResolvedFlag) const 
{
  if(ResolvedFlag == 0)
    {
      return chi2_corrected[2];
    }
  else
    {
       return Resolvedchi2_corrected[2];
    }
}
//____________________________________________________________________________
Float_t SBSGRINCH_Cluster::Getchi2_corrected_prob_pi(Int_t ResolvedFlag) const
{
  if(ResolvedFlag == 0)
    {
      return chi2_corrected_prob[0];
    }
  else
    {
       return Resolvedchi2_corrected_prob[0];
    }
}
//____________________________________________________________________________
Float_t SBSGRINCH_Cluster::Getchi2_corrected_prob_k(Int_t ResolvedFlag) const 
{
  if(ResolvedFlag == 0)
    {
      return chi2_corrected_prob[1];
    }
  else
    {
       return Resolvedchi2_corrected_prob[1];
    }
}
//____________________________________________________________________________
Float_t SBSGRINCH_Cluster::Getchi2_corrected_prob_p(Int_t ResolvedFlag) const 
{
  if(ResolvedFlag == 0)
    {
      return chi2_corrected_prob[2];
    }
  else
    {
       return Resolvedchi2_corrected_prob[2];
    }
}
//____________________________________________________________________________
Float_t    SBSGRINCH_Cluster::Getchi2_corrected_prob(Int_t flag, 
						   Int_t ResolvedFlag) const
{
  if(ResolvedFlag == 0)
    {
      return chi2_corrected_prob[flag];
    }
  else
    {
      return Resolvedchi2_corrected_prob[flag];
    }
}
//____________________________________________________________________________
Int_t  SBSGRINCH_Cluster::Getnoise_cut_success(Int_t ResolvedFlag) const
{
  if(ResolvedFlag == 0)
    {
      return fnoise_cut_success;
    }
  else
    {
      return fResolved_noise_cut_success;
    }
}
//_____________________________________________________________________________
Int_t SBSGRINCH_Cluster::Test( const SBSGRINCH_Hit* theHit, Float_t par1, 
			     Float_t par2, Float_t par3 ) const
{
  // Check if theHit is sufficiently close to any of the hits of the
  // cluster to be considered part of this cluster
  //
  // Parameters: par1
  //             par2
  //             par3
  //

  Float_t dx,dy,dist;
  if( !fHitList ) return 0;
  TIter next( fHitList );

  while( SBSGRINCH_Hit* pHit = static_cast<SBSGRINCH_Hit*>( next() )) {
    dx   = theHit->GetX() - pHit->GetX();
    dy   = theHit->GetY() - pHit->GetY();
    dist = sqrt( dx*dx + dy*dy );
    if( dist>par1 && TMath::Abs(dx)<par2 && TMath::Abs(dy)<par3 )
      return 1;
  }
  return 0;
}
//_____________________________________________________________________________
Int_t SBSGRINCH_Cluster::FindLocalMaximumNumber( )
{
  // determine how many local maxima are in the cluster
  //

  if( !fHitList ) return 0;

  fLocalMaximumNumber = 0;

  TIter next( fHitList );
  while( SBSGRINCH_Hit* pHit = static_cast<SBSGRINCH_Hit*>( next() ))
    pHit->SetVeto(0); // all veto flags zeroed

  next.Reset();
  while( SBSGRINCH_Hit* pHitSave = static_cast<SBSGRINCH_Hit*>( next() )) {
    Int_t I_Position = pHitSave->GetI();
    Int_t J_Position = pHitSave->GetJ();

    TIter next1( next );
    while( SBSGRINCH_Hit* pHit = static_cast<SBSGRINCH_Hit*>( next1() )) {
      Int_t intdist = TMath::Abs(I_Position - pHit->GetI()) + 
	TMath::Abs(J_Position - pHit->GetJ());
      if(intdist == 1) {
	//set the flag that forbids the Hit to be a local maximum
	if (pHit->Compare(pHitSave) == -1)
	  pHit->SetVeto(1); 
	else
	  pHitSave->SetVeto(1); 
      }
    }
  }
    
    //--old--
#if 0    
  Int_t HitNumber = fHitList->GetSize(); 
  for (i = 0; i < HitNumber; i++) {
    iteration = 0;

    next.Reset();
    while( SBSGRINCH_Hit* pHit = static_cast<SBSGRINCH_Hit*>( next() )) {
      //FIXME: urgh
      if(iteration < i) {
	iteration++; 
	continue;
      }
      if(iteration == i) {
	I_Position = pHit->GetI();
	J_Position = pHit->GetJ();
	pHitSave = pHit; 
	iteration++;
	continue;
      }
      intdist = TMath::Abs(I_Position - pHit->GetI()) + 
	TMath::Abs(J_Position - pHit->GetJ());
      if(intdist == 1) {
	//set the flag that forbids the Hit to be a local maximum
	if (pHit->Compare(pHitSave) == -1)
	  pHit->SetVeto(1); 
	else
	  pHitSave->SetVeto(1); 
      }
      iteration++;
    }
  }
#endif

  next.Reset();
  while( SBSGRINCH_Hit* pHit = static_cast<SBSGRINCH_Hit*>( next() )) {
    if(!pHit->GetVeto())
      fLocalMaximumNumber++;
  }

  return fLocalMaximumNumber;
}

//---------------------------------------------------------------------------
Int_t SBSGRINCH_Cluster::FindResolvedClusterElements
  ( const SBSGRINCH_Hit* theLocalMaximum, SBSGRINCH_Cluster* resolvedCluster,
    TClonesArray* resolvedHits )

  // It finds the hits of a new resolved cluster from the present one. 
  // Arguments:
  //   theLocalMaximum: local maximum hit to use as the center of the
  //                    new resolved cluster
  //   resolvedCluster: the new resolved cluster
  //   resolvedHits: array of new resolved hits. The new hits of the
  //                 new resolved cluster are appended to this array.
  //
  // resolvedCluster and resolvedHits are allocated by the caller
  
  // factors: factors the charge of the corresponding hit should be scaled
  //          (the charge of the hit could be fragmented when the same 
  //           Hit is attributed to several resolved clusters).
{
  if( !theLocalMaximum || ! resolvedCluster || !resolvedHits ) {
    Error("FindResolvedClusterElements", "internal error: one or more "
	  "NULL arguments!");
    return -1;
  }

  Int_t fIRef = 0, fJRef = 0;
  Float_t ChargeofMaximum = 0.;
  Int_t isaved = -1;
  Int_t fISaved[4], fJSaved[4];
  Float_t fChargeSaved[4];

  for (Int_t i = 0; i < 4 ; i++) {
    fISaved[i] = 0;
    fJSaved[i] = 0;
    fChargeSaved[i] = 0;
  }

  Int_t nResolvedHits = resolvedHits->GetLast()+1;

  const SBSGRINCH_Hit* pHit = theLocalMaximum;
  SBSGRINCH_Hit* newHit = new( (*resolvedHits)[nResolvedHits++] ) 
    SBSGRINCH_Hit( *pHit );
  resolvedCluster->Insert( newHit );
  fIRef = pHit->GetI();
  fJRef = pHit->GetJ();
  ChargeofMaximum = pHit->GetADC();

  // take as the elements of the new clusters the elements next to the maximum.

  TIter next( fHitList );
  while( (pHit = static_cast<SBSGRINCH_Hit*>( next() ))) {
    if( ( TMath::Abs(pHit->GetI() - fIRef) + 
	  TMath::Abs(pHit->GetJ() - fJRef) ) == 1) {
      Int_t fIofHit = pHit->GetI();
      Int_t fJofHit = pHit->GetJ();
      Float_t SumofCharges = 0;
      TIter next2( fHitList );
      while( SBSGRINCH_Hit* pHit1 = static_cast<SBSGRINCH_Hit*>( next2() )) {
	// check if the Hit belongs to other local maxima.
	// Make the sum of the charges of the local maximum the Hit
	// belongs to. 
	if( (!pHit1->GetVeto()) &&
	    ( TMath::Abs(pHit1->GetI() - fIofHit) + 
	      TMath::Abs(pHit1->GetJ() - fJofHit)  == 1) ) {
	  SumofCharges = SumofCharges + pHit1->GetADC();
	}
	if( (!pHit1->GetVeto()) &&
	    ( TMath::Abs(pHit1->GetI() - fIofHit) + 
	      TMath::Abs(pHit1->GetJ() - fJofHit)  == 2) ) {
	  // for the posibility of more distant Hits to belong 
	  // to a local maximum see below.  
	  Int_t fItocheck = pHit1->GetI();
	  Int_t fJtocheck = pHit1->GetJ();
	  Int_t BelongingFlag = 1;
	  TIter next5( fHitList );
	  while( SBSGRINCH_Hit* pHit2 = 
		 static_cast<SBSGRINCH_Hit*>( next5() )) {
	    if(( TMath::Abs(pHit2->GetI() - fItocheck) + 
		 TMath::Abs(pHit2->GetJ() - fJtocheck)  == 1) 
	       && (pHit->GetADC() > pHit2->GetADC()) ) {
	      BelongingFlag = 0;
	    }
	  }
	  if(BelongingFlag)
	    SumofCharges = SumofCharges + pHit1->GetADC();
	}
      }
      // If the Hit belong to other local maximum divide 
      // its charge accordingly.   
      isaved++;
      fISaved[isaved] = fIofHit;
      fJSaved[isaved] = fJofHit;
      fChargeSaved[isaved] = pHit->GetADC();
      Float_t factor = ChargeofMaximum/SumofCharges;
      newHit = new( (*resolvedHits)[nResolvedHits++] ) SBSGRINCH_Hit( *pHit );
      newHit->SetADC(int(pHit->GetADC()*factor));
      resolvedCluster->Insert( newHit );
    }
  }

  // check for Hits more distant from the new cluster center;
  next.Reset();
  while( (pHit = static_cast<SBSGRINCH_Hit*>( next() ))) {
    if(!pHit->GetVeto()) continue; 
    // do not analyze local maximum 
    if( ( TMath::Abs(pHit->GetI() - fIRef) 
	  + TMath::Abs(pHit->GetJ() - fJRef) ) == 2) {
      Int_t fIofHit = pHit->GetI();
      Int_t fJofHit = pHit->GetJ();
      Float_t CheckedCharge = pHit->GetADC();
      Int_t ClusterFlag = 1;
      for(Int_t icheck = 0; icheck <= isaved; icheck++) {
	if( (TMath::Abs(fISaved[icheck] - fIofHit) == 1) || 
	    (TMath::Abs(fJSaved[icheck] - fJofHit) == 1) ) {
	  if (CheckedCharge > fChargeSaved[icheck]) 
	    ClusterFlag = 0;
	}
      }
      if(ClusterFlag) {
	// check if the Hit belongs to other local maximum too.
	Float_t SumofCharges = 0;
	TIter next4( fHitList );
	while( SBSGRINCH_Hit* pHit1 = static_cast<SBSGRINCH_Hit*>( next4() )){
	  if( (!pHit1->GetVeto()) &&
	      ( TMath::Abs(pHit1->GetI() - fIofHit) + 
		TMath::Abs(pHit1->GetJ() - fJofHit)  == 1) ){
	    SumofCharges = SumofCharges + pHit1->GetADC();
	  }
	  if( (!pHit1->GetVeto()) &&
	      ( TMath::Abs(pHit1->GetI() - fIofHit) + 
		TMath::Abs(pHit1->GetJ() - fJofHit)  == 2) ) {
	    Int_t fItocheck = pHit1->GetI();
	    Int_t fJtocheck = pHit1->GetJ();
	    Int_t BelongingFlag = 1;
	    TIter next5( fHitList );
	    while( SBSGRINCH_Hit* pHit2 = 
		   static_cast<SBSGRINCH_Hit*>( next5() )) {
	      if(( TMath::Abs(pHit2->GetI() - fItocheck) + 
		   TMath::Abs(pHit2->GetJ() - fJtocheck)  == 1) 
		 && (CheckedCharge > pHit2->GetADC()) ){
		BelongingFlag = 0;
	      }
	    }
	    if(BelongingFlag) 
	      SumofCharges = SumofCharges + pHit1->GetADC();
	  }
	}		
	Float_t factor;
	if (SumofCharges == 0){
	  factor  = 0;
	  // In the above algorithm if a charge is located 
	  // two pads distance from a local maximum it could be
	  // assigned to its cluster even if there is a not hit pad
	  // between them (and hence it does not really belongs 
	  // to that cluster). In this case SumofCharges = 0.
	  // The Hit will continue stay in this cluster but with 
	  // charge 0 and hence with no effect in the analysis.
	  // This check however does not avoid completely this 
	  // problem FIX ME. 
	}
	else
	  factor = ChargeofMaximum/SumofCharges;

	newHit = new( (*resolvedHits)[nResolvedHits++] ) SBSGRINCH_Hit( *pHit );
	newHit->SetADC(int(pHit->GetADC()*factor));
	resolvedCluster->Insert( newHit );
      }
    }
  }
  return 0;
}

//_____________________________________________________________________________
void SBSGRINCH_Cluster::Show( FILE* fout1, FILE* fout2 ) const
{
  //Print info about cluster statistics and all the hits in the cluster.

  fprintf(fout2," X ");
  fprintf(fout1," %4f",fXcenter);
  fprintf(fout2,"%4f",fXcenter);
  fprintf(fout2," Y ");
  fprintf(fout1," %4f",fYcenter);
  fprintf(fout2,"%4f",fYcenter);
  fprintf(fout2," Charge ");
  fprintf(fout1," %4f",fCharge);
  fprintf(fout2,"%4f",fCharge); 
  fprintf(fout2," Element Number ");
  fprintf(fout1," %4d", GetNHits() ); 
  fprintf(fout2," %4d", GetNHits() );
  fprintf(fout2,"\n");
  fprintf(fout2," Number of local maximum");
  fprintf(fout2," %4d", fLocalMaximumNumber);
  fprintf(fout2,"\n");
  if(IsMIP())
    {
      Float_t Zero = 0.; 
      fprintf(fout2," Theta Angle ");
      fprintf(fout1," %4f", Zero);
      fprintf(fout2,"%4f", Zero);
      fprintf(fout2," Phi Angle ");
      fprintf(fout1,"% 4f", Zero);
      fprintf(fout2,"%4f", Zero);
      fprintf(fout2," Cherenkov Angle ");
      fprintf(fout1," %4f", Zero);
      fprintf(fout2,"%4f", Zero);
      fprintf(fout2,"\n"); 
    }
  else
    {
      fprintf(fout2," Theta Angle ");
      fprintf(fout1," %4f", fTheta_photon);
      fprintf(fout2,"%4f", fTheta_photon);
      fprintf(fout2," Phi Angle ");
      fprintf(fout1,"% 4f", fPhi_photon);
      fprintf(fout2,"%4f", fPhi_photon);
      fprintf(fout2," Cherenkov Angle ");
      fprintf(fout1," %4f",fAngle);
      fprintf(fout2,"%4f",fAngle);
  fprintf(fout2,"\n"); 
    }
  fprintf(fout2," Cluster element list: ");
  fprintf(fout2," \n");

  // fHitList->Print();
  
  if( fHitList ) {
    TIter next(fHitList);
    while (SBSGRINCH_Hit* d = (SBSGRINCH_Hit*)next()) {   
      d->Show(fout1, fout2);
    }
  }

  return;
}

//_____________________________________________________________________________
void SBSGRINCH_Cluster::Show( FILE* fout1) const
{
  //Print info about cluster statistics and all the hits in the cluster.

  fprintf(fout1,"% 4f",fXcenter);
  fprintf(fout1,"% 4f",fYcenter);
  fprintf(fout1,"% 4f",fCharge);
  fprintf(fout1," %4d", GetNHits() ); 
  if(IsMIP())
    {
      Float_t Zero = 0.;
      fprintf(fout1," %4f", Zero);
      fprintf(fout1," %4f", Zero);
      fprintf(fout1," %4f", Zero);
    }
  else
    {
      fprintf(fout1," %4f", fTheta_photon);
      fprintf(fout1," %4f", fPhi_photon);
      fprintf(fout1," %4f",fAngle);
    }

  return;
}


//_____________________________________________________________________________
void SBSGRINCH_Cluster::ShowElements(FILE * fout)  const
{
  fprintf(fout,"\n");
  fprintf(fout, "Cluster Number Elements ");
  fprintf(fout,"%4d", GetNHits());
  fprintf(fout,"\n");
  fprintf(fout," X, Y ");
  fprintf(fout,"%4f,%4f",fXcenter,fYcenter);
  fprintf(fout,"; Charge  = ");
  fprintf(fout,"%4f",fCharge);
  fprintf(fout,"\n");
  fprintf(fout,"Cluster elements: ");
  fprintf(fout,"\n");
  Show(fout);
}


//_____________________________________________________________________________

ClassImp(SBSGRINCH_Hit)
ClassImp(SBSGRINCH_Cluster)

