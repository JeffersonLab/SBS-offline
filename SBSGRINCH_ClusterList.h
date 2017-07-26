#ifndef ROOT_SBSGRINCH_ClusterList
#define ROOT_SBSGRINCH_ClusterList

//****************************************************************
// The header file for the "ClusterRICH.C"
//****************************************************************

#include "TObject.h"
#include "TList.h"

class THaTrack;
class TClonesArray;

#include <iostream>
#include <stdio.h>

// Data class contains the data of a single pad (of a cluster).
// (i,j) coordinate of a pad; (x,y) coordinated of a pad; ADC value
// and a Veto to point it cannot be a local minimum inside the cluster
// useful to tell overlapping clusters apart). 
// Above all it contains a module to compare adc values. This allows 
// the ordering according the adc values itself of the components of a cluster.

class SBSGRINCH_Hit : public TObject {

public:
  SBSGRINCH_Hit() 
    : fFlag(0), fVeto(0) {}
  SBSGRINCH_Hit( Int_t number, Int_t ADC, Int_t i, Int_t j, 
	       Float_t x, Float_t y ) :
    fNumber(number), fADC(ADC), fI(i), fJ(j), fX(x), fY(y), fFlag(0), 
    fVeto(0) {}
  virtual ~SBSGRINCH_Hit() {}

  void       Show(FILE * fout1);
  void       Show(FILE * fout1, FILE * fout2);

  Int_t      GetNumber()   const {return fNumber;}
  Float_t    GetX()        const {return fX;}
  Float_t    GetY()        const {return fY;}
  Int_t      GetI()        const {return fI;}
  Int_t      GetJ()        const {return fJ;}
  Int_t      GetADC()      const {return fADC;}
  Int_t      GetFlag()     const {return fFlag;}
  Int_t      GetVeto()     const {return fVeto;} 
  void       SetNumber( Int_t number ) {fNumber = number;}
  void       SetX( Float_t x )         {fX = x;}
  void       SetY( Float_t y )         {fY = y;}
  void       SetI( Int_t i )           {fI = i;}
  void       SetJ( Int_t j )           {fJ = j;}
  void       SetADC( Int_t ADC )       {fADC = ADC;}
  void       SetFlag( Int_t Flag )     {fFlag = Flag;}
  void       SetVeto( Int_t Veto )     {fVeto = Veto;}

  virtual Int_t   Compare( const TObject* ) const;
  virtual Bool_t  IsSortable() const { return kTRUE; }

private:
  Int_t     fNumber;
  Int_t     fADC;
  Int_t     fI;
  Int_t     fJ;
  Float_t   fX;
  Float_t   fY;
  Int_t     fFlag;
  Int_t     fVeto;

  ClassDef(SBSGRINCH_Hit,0)   //A hit in the RICH
};


// --------------------------------------------------------------

// ClusterElement: class of elements that make up a cluster. They are orderd
// according to the ADC value read in the corresponding Pad (from the highest
// value to the smallest).
// ClusterCompon is the general (ADT) class from which all the other classes
// are derived.
// HeadClusterElement and TailClusterElement are just null elements, 
// that just contain the pointer to the first and the last ClusterElement 
// respectively.
// InternalClusterElement contain a pointer to the next element in the list 
// of the cluster components. It contains also a pointer to the Class Hit
// that contains the information (ADC (X,Y) etc. relative).
// The definition of the classes above allows:
// 1) to avoid creation of arrays of predefined sizes of Hit. 
//         Only the elements  needed will be created to make up the cluster.
// 2) To just accomplish further algorithm  of ordering cluster element 
//  just changing module Hit::Compare above.


// A cluster of hits in the RICH.
// Calculates charge-weighted center position of all hits.
// Hits belonging to the cluster are kept in a TList.
// Hits must be TObjects.
// If you want a sorted list, just replace TList with TSortedList

class SBSGRINCH_Cluster : public TObject {

public:
  SBSGRINCH_Cluster();
  SBSGRINCH_Cluster( const SBSGRINCH_Cluster& rhs );
  SBSGRINCH_Cluster& operator=( const SBSGRINCH_Cluster& rhs );

  ~SBSGRINCH_Cluster() { delete fHitList; }


  void       Clear( Option_t* opt="" );
  Float_t    Dist( const SBSGRINCH_Cluster* c ) const;
  void       Insert( SBSGRINCH_Hit* theHit );
  void       Insert( SBSGRINCH_Hit* theHit, Float_t factor);
  void       Insert_Photon(Int_t flag, Float_t angle, Int_t ResolvedFlag);
  void       Insert_Photon(Int_t flag, Float_t angle, Int_t ResolvedFlag, 
			   Float_t expected_angle,
			   Float_t central_momentum_expected_angle);
  void       Insert_chi2(Int_t flag, Float_t chi2);
  void       Insert_chi2(Int_t flag, Float_t chi2, Int_t ResolvedFlag);
  void       Insert_chi2_corrected (Int_t flag, Float_t chi2, 
				    Int_t ResolvedFlag);
  void       Insert_N_chi2_corrected_Photon (Int_t flag, Int_t N_Photon, 
				    Int_t ResolvedFlag);
  void       Insert_MaximumLikelihood (Int_t flag, Float_t MaximumLikelihood, 
				     Int_t ResolvedFlag);
  Bool_t     IsMIP()      const              { return fMIPflag; }
  Bool_t     GetPionChi2AnalysisFlag() const { return fPionChi2AnalysisFlag; }
  Bool_t     GetKaonChi2AnalysisFlag() const { return fKaonChi2AnalysisFlag; }
  Bool_t    GetProtonChi2AnalysisFlag() const {return fProtonChi2AnalysisFlag;}
  Int_t      GetNHits()   const              { return fHitList->GetSize(); }
  TList*     GetHitList()                    { return fHitList; }
  Int_t      GetLocalMaximumNumber()         { return fLocalMaximumNumber;} 
  Float_t    GetXcenter() const              { return fXcenter; }
  Float_t    GetYcenter() const              { return fYcenter; }
  Float_t    GetCharge()  const              { return fCharge; }
  void       SetXcenter(Float_t value)       { fXcenter = value; }
  void       SetYcenter(Float_t value)       { fYcenter = value; }
  void       SetCharge(Float_t value)        { fCharge = value; }
  void       SetMIPflag(Bool_t value)        { fMIPflag = value; }
  void       SetPionChi2AnalysisFlag(Bool_t value) 
  { fPionChi2AnalysisFlag = value; }
  void       SetKaonChi2AnalysisFlag(Bool_t value) 
  { fKaonChi2AnalysisFlag = value; }
  void       SetProtonChi2AnalysisFlag(Bool_t value) 
  { fProtonChi2AnalysisFlag = value; }
  void       Setchi2_prob(Int_t flag, Float_t chi2_value, Int_t N_Photon, 
			  Int_t ResolvedFlag);
  void       Setchi2_corrected_prob(Int_t flag, Float_t chi2_value, 
				    Int_t N_Photon, Int_t ResolvedFlag);
  void       Setnoise_cut_success(Int_t value, Int_t ResolvedFlag);
  Int_t      FindLocalMaximumNumber();
  Int_t FindResolvedClusterElements (const SBSGRINCH_Hit* localMaximum,
				     SBSGRINCH_Cluster* resolvedCluster,
				     TClonesArray* resolvedHits );
  Float_t    GetTheta_photon() const         { return fTheta_photon; }
  Float_t    GetPhi_photon() const           { return fPhi_photon; }
  Float_t    GetAngle()   const              { return fAngle; }
  THaTrack*  GetTrack()   const              { return fTrack; }
  void       MakeMIP( Bool_t flag = kTRUE );
  void       SetMIP( SBSGRINCH_Cluster* mip )  { fMIP = mip; }
  void       SetFictious_MIP_Flag(Int_t value) { fFictious_Mip_Flag = value; }
  void       SetTheta_photon ( Float_t Theta_photon) 
    { fTheta_photon = Theta_photon; }
  void       SetPhi_photon (Float_t Phi_photon) 
    { fPhi_photon = Phi_photon; }
  void       SetAngle( Float_t angle )       { fAngle = angle; }
  void       SetTrack( THaTrack* track )     { fTrack = track; }
   // to be applied on Mips unless otherwise specified
  void       SetN_Photon(Int_t i, Int_t Value, Int_t ResolvedFlag);
  void       Setangle(Int_t i, Float_t Value, Int_t ResolvedFlag);
  void       Setangle_corrected(Int_t i, Float_t Value, Int_t ResolvedFlag);
  Int_t      GetN_Photon(Int_t i, Int_t ResolvedFlag) const; 
  Float_t    Getangle(Int_t i,Int_t ResolvedFlag ) const;
  Float_t    Getangle_corrected(Int_t i,Int_t ResolvedFlag ) const;
  // Convenience functions for access via global variables
  Int_t      GetNphot_pi(Int_t ResolvedFlag) const;
  Int_t      GetNphot_k(Int_t ResolvedFlag)  const; 
  Int_t      GetNphot_p(Int_t ResolvedFlag)  const; 
  Float_t    Getangle_pi(Int_t ResolvedFlag) const;
  Float_t    Getangle_k(Int_t ResolvedFlag)  const;
  Float_t    Getangle_p(Int_t ResolvedFlag)  const;
  Float_t    Getangle_corrected_pi(Int_t ResolvedFlag) const;
  Float_t    Getangle_corrected_k(Int_t ResolvedFlag)  const;
  Float_t    Getangle_corrected_p(Int_t ResolvedFlag)  const;
  Int_t      GetN_chi2_phot_pi(Int_t ResolvedFlag) const;
  Int_t      GetN_chi2_phot_k(Int_t ResolvedFlag)  const;
  Int_t      GetN_chi2_phot_p(Int_t ResolvedFlag)  const;
  Float_t    Getchi2_pi() const { return chi2[0]; } // only for single clysters
  Float_t    Getchi2_k() const  { return chi2[1]; } // only for single clysters
  Float_t    Getchi2_p() const  { return chi2[2]; } // only for single clysters
  Float_t    Getchi2_pi(Int_t ResolvedFlag) const;
  Float_t    Getchi2_k(Int_t ResolvedFlag)  const;
  Float_t    Getchi2_p(Int_t ResolvedFlag)  const;
  Float_t    Getchi2_prob_pi(Int_t ResolvedFlag) const;
  Float_t    Getchi2_prob_k(Int_t ResolvedFlag)  const;
  Float_t    Getchi2_prob_p(Int_t ResolvedFlag)  const;
  Float_t    Getchi2_prob(Int_t flag, Int_t ResolvedFlag)  const;
  Int_t      GetN_chi2_corrected_phot_pi(Int_t ResolvedFlag) const;
  Int_t      GetN_chi2_corrected_phot_k(Int_t ResolvedFlag)  const; 
  Int_t      GetN_chi2_corrected_phot_p(Int_t ResolvedFlag)  const; 
  Float_t    Getchi2_corrected_pi(Int_t ResolvedFlag) const;    
  Float_t    Getchi2_corrected_k(Int_t ResolvedFlag)  const;
  Float_t    Getchi2_corrected_p(Int_t ResolvedFlag)  const;
  Float_t    Getchi2_corrected_prob_pi(Int_t ResolvedFlag) const;    
  Float_t    Getchi2_corrected_prob_k(Int_t ResolvedFlag)  const;
  Float_t    Getchi2_corrected_prob_p(Int_t ResolvedFlag)  const;
  Float_t    Getchi2_corrected_prob(Int_t flag, Int_t ResolvedFlag)  const;
  Int_t      Getnoise_cut_success(Int_t ResolvedFlag) const;
  Int_t      Getunresolved_noise_cut_success() const
  { return fnoise_cut_success; };
  Int_t      Getresolved_noise_cut_success() const
  { return fResolved_noise_cut_success; };
  Int_t      GetFictious_Mip_Flag() const { return fFictious_Mip_Flag; };
 
  void       Show( FILE* fout1) const;
  void       Show( FILE* fout1, FILE* fout2 ) const;
  void       ShowElements( FILE* fout ) const;
  Int_t      Test( const SBSGRINCH_Hit* theHit,  Float_t par1, 
		   Float_t par2,  Float_t par3 ) const;

private:
  TList*     fHitList;   //List of hits belonging to this cluster
  Int_t fLocalMaximumNumber; // number of local maxima.
  Bool_t     fMIPflag;   //True if this cluster is a MIP
  Int_t fFictious_Mip_Flag; 
                         // only for MIPs; flag that is equal to:
                         // 3           when the analysis (that is 
                         //             MIP_through_interception value
                         //             in the data base) forced the MIP 
                         //             to be the interception between the
                         //             track and the PAD Plane regardless
                         //             of the cluster pattern in the Pad 
                         //             plane.
                         // 2           when the MIP was forced (by the 
                         //             MIP_through_interception value in the 
                         //             data base) to be the interception 
                         //             between the track and the PAD plane
                         //             because a cluster inside the MIP 
                         //             search radius has not been found. 
                         // 1           when the MIP was forced (by the 
                         //             MIP_through_interception value in the 
                         //             data base) to be the interception 
                         //             between the track and the PAD plane
                         //             because both of this conditions were
                         //             fulfilled: 
                         //                       a) a cluster inside the MIP 
                         //                          search radius was not 
                         //                          found and 
                         //                       b) no MIP spot was not 
                         //                          expected on the pad plane
                         //                          because the track 
                         //                          interception with it
                         //                          felt in a not sensible 
                         //                          region of the rich
                         // 0           the MIP is the maximum charge cluster
                         //             inside the MIP search radius.
  Bool_t fPionChi2AnalysisFlag; // true if the cluster is employed in the 
                               // chi square calculation performed under the
                               // hypothesis the particle crossing the Rich
                               // is a pion
  Bool_t fKaonChi2AnalysisFlag; // true if the cluster is employed in the 
                               // chi square calculation performed under the
                               // hypothesis the particle crossing the Rich
                               // is a kaon
  Bool_t fProtonChi2AnalysisFlag; // true if the cluster is employed in the 
                               // chi square calculation performed under the
                               // hypothesis the particle crossing the Rich
                               // is a proton
  Float_t    fXcenter;   // (Sum of x*adc)/(sum adc)  of all hits in the list
  Float_t    fYcenter;   // (Sum of y*adc)/sum(adc) of all hits in the list
  Float_t    fCharge;    //Sum of adc values of all hits
  Float_t fTheta_photon; // Theta angle in the RICH system of the photon 
                         // associated with the cluster
  Float_t fPhi_photon;   // Phi angle in the RICH system of the photon 
                         // associated with the cluster.
  Float_t    fAngle;     //Calculated angle wrt particle ray (not used yet)
  SBSGRINCH_Cluster* fMIP; //Pointer to MIP cluster belonging to this cluster
  THaTrack*  fTrack;     //Track associated with this cluster (only for MIPs)
  Int_t N_Photon[3];     //Only for MIPs; number of clusters whose angles 
                         // with respect to the MIP considered are in the 
                         // fiducial region in case the MIP is generated by:
                         // a pion:    N_Photon[0]
                         // a kaon:    N_Photon[1]
                         // a proton:  N_Photon[2]
                         // respectively.
  Float_t angle[3];      // only for MIPs; average of the angles of the 
                         // clusters in the fiducial region of the MIP 
                         // considered. Three cases are considered as 
                         // possible:
                         // the MIP is generated by a pion (angle[0] = 
                         // average of the N_Photon[0] clusters);
                         // the MIp is generated by a kaon (angle[1]); 
                         // The MIP is generated by a proton(angle[2]).
  Float_t angle_corrected[3]; // only for MIPs;  
                         // it is equal to:
                         //
                         // angle[i] - expected_angle[i] + 
                         //           central_momentum_angle[i]
                         //
                         // where:
                         //
                         //  - angle[i] is the variable angle[i] described
                         //             in the item above;
                         //  - expected_angle[i] is the expected cherenkov
                         //             photon angle considering the particle
                         //             momentum;
                         //  - central_momentum_angle[i] is the expected
                         //             cherenkov photon angle for particles
                         //             having the momentum just equal to 
                         //             HRS right arm central momentum.
                         //
                         // Unlike the variable angle[i], the distribution of 
                         // "angle_corrected[i]" is not affected by the 
                         // smearing caused by the particle momentum variation 
                         // inside HRS acceptance. 
  Int_t N_chi2_Photon[3]; // Only for MIPS; number of photon the chi2 test
                         // is performed (see chi2[3] definition).
                         //The meaning of the index of the array is the same 
                         // of N_Photon[3]. 
  Float_t chi2[3];       // only for MIPS; chi square value calculated
                         // as Sum[(angle-expected_angle)/sigma)]**2
                         // where angle is the cherenkov angles of the single
                         // clusters and expected_angle is the expected 
                         // cherenkov angle if the mip is a pion (chi2[0]),
                         // a kaon (chi2[1]) and a proton (chi2[2]). Sigma
                         // is the variange of the single cluster distributon
                         // (from the database).
  Float_t chi2_prob[3];  // only for MIPs; probability of chi2.
  Int_t N_chi2_corrected_Photon[3]; // Only for MIPs; number of photon 
                                    // the "corrected" chi2 test is performed 
                                    // (see chi2_corrected[3] definition).
                                    //The meaning of the index of the array 
                                    // is the same of N_Photon[3].
  Float_t chi2_corrected[3];        // only for MIPs, chi2 values when 
                                    // clusters  originated by noise are cut 
                                    // away.
  Float_t chi2_corrected_prob[3];   // only for MIPs, probability of 
                                    // chi2_corrected.
  Int_t fnoise_cut_success;          // only for mips. It is equal to 1 when 
                                    // the algorithm to throw away the noise 
                                    // was succesful; 0 otherwise (see 
                                    // fineprocess and clearnoise functions in
                                    // ThaRICH.cxx
  Float_t MaximumLikelihood[3];      // only for MIPS; Maximum Likeliood 
                                    // Algorithm defined as:
                                    // 
                           // (1/N)*Sum(log{1+(1/(sqrt(2*pi)*sigma*epsilon))*
                           // *exp[(angle-expected_angle)""2/(2*sigma**2)]})
                                    //
                                    // where 
                                    // angle            is the cherenkov angle 
                                    //                  of the single clusters 
                                    // expected_angle   is the expected 
                                    //                  cherenkov angle if the
                                    //                  mip is a pion 
                                    //                  (MaximumLikelihood[0]),
                                    //                  a kaon 
                                    //                  (MaximumLikelihood[1]),
                                    //                  and a proton 
                                    //                  (MaximumLikelihood[2]).
                                    // Sigma            is the variane of the 
                                    //                  single cluster 
                                    //                  distributon
                                    //                  (from the database).
                                    // epsilon          is an effecttive 
                                    //                  cut parameter
                                    //                  (from the database).
                                    // N                is equal to 
                                    //                  log{[1+(1/(sqrt(2*pi)*
                                    //                  sigma*epsilon)]
  Int_t N_MaximumLikelihood_Photon[3]; // Only for MIPS; number of photon the 
                                    // Maximum Likelihood calculation
  Int_t ResolvedN_Photon[3];        // The equivalent of N_Photon when 
                                    // considering resolved clusters.
  Float_t Resolvedangle[3];         // The equivalent of angle calculated by 
                                    // resolved clusters.
  Float_t Resolvedangle_corrected[3]; // The equivalent of angle_corrected 
                                    // calculated by resolved clusters.
  Int_t ResolvedN_chi2_Photon[3];   // The equivalent of N_chi2_Photon 
                                    // when considering resolved cluster.
  Float_t Resolvedchi2[3];          // The equivalent of chi2 when considering
                                    // resolved cluster.
  Float_t Resolvedchi2_prob[3];     // The equivalent of chi2_prob when
                                    // considering resolved cluster.
  Int_t ResolvedN_chi2_corrected_Photon[3]; // The equivalent of 
                                            // N_chi2_corrected_Photon when
                                            // considering resolve3d cluster
  Float_t Resolvedchi2_corrected[3];        // The equivalent of 
                                            // chi2_corrected when considering
                                            // resolved clusters.
  Float_t Resolvedchi2_corrected_prob[3];   // The equivalent of 
                                            // chi2_corrected_prob when
                                            // considering resolved clusters.
  Int_t fResolved_noise_cut_success;         // The equivalent of 
                                            // fnoise_cut_success for resolved
                                            // clusters
  Float_t ResolvedMaximumLikelihood[3];     // The equivalent of 
                                            // MaximumLikelihood for resolved 
                                            // clusters.
  Int_t ResolvedN_MaximumLikelihood_Photon[3]; // The equivalent of 
                                            // MaximumLikelihood_Photon for 
                                            // resolved clusters.

  ClassDef(SBSGRINCH_Cluster,0)  //A cluster of hits in the RICH
};


// ---------------- inlines -------------------------------------

inline
void SBSGRINCH_Cluster::MakeMIP( Bool_t flag )
{
  fMIPflag = flag;
  if( flag )
    fMIP = this;
  else
    fMIP = NULL;
}

#endif


















