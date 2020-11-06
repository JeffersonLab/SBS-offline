//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Jan  7 11:54:23 2016 by ROOT version 5.34/32
// from TTree T/Geant4 SBS Simulation
// found on file: gep_spin_transport_Sx.root
//////////////////////////////////////////////////////////

#ifndef __G4SBS_TREE_H
#define __G4SBS_TREE_H

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TString.h>

// ----------------------------
// This class is useful to activate *all* variables from the tree 
// stored in g4sbs output root files. 
// In libsolgem (libsbsgem) it is used by the class TSBSGeant4File.
// Note it can perfectly be used in standalone.
// 
// It is more particularly dedicated to unfold files obtained with GEp setup.
// It includes information of CDET, GEp ECal, FT, FPP1&2, HCal.
// For more info check the following link: 
// https://hallaweb.jlab.org/wiki/index.php/Documentation_of_g4sbs#ROOT_Tree_Structure 
// 

// Header file for the classes stored in the TTree if any.
#include <vector>
//#include "sbstypes.hh"
#include "g4sbs_data.h"
//#include "g4sbs_types.h"
// Fixed size dimensions of array or collections stored in the TTree if any.

/*
struct gem_branch {
   int id;
   const char* name;
   TSBSGeant4::GEMData_t tree;
   const char* Track_name;
   TSBSGeant4::TrackerData_t Track_tree;
   gem_branch(int det_id = 0, const char *b_name = 0, const char *b_Track_name = 0)
     : id(det_id), name(b_name), Track_name(b_Track_name)
   {
   }
};
*/
#define treeName "T"

using namespace std;

class g4sbs_tree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain
   //Exp_t           fExpt;    // Choose experiment type: defined in "g4sbs_types.h"
   //bool            fPythia;// needed to turn on/off the reading of the pythia variables
   // EPAF: those were inherited from an attempt of "standard tree" 
   // for the "raw" analysis of g4sbs files
   // these have nothing to do with digitization
   // bool            fEcalBox;// needed to turn on/off the reading of the ECAL_box data
   // bool            fHcalBox;// needed to turn on/off the reading of the HCAL_box data

   // Declaration of leaf types

   // Event variables
   Double_t        ev_count;
   Double_t        ev_rate;
   Double_t        ev_solang;
   Double_t        ev_sigma;
   Double_t        ev_W2;
   Double_t        ev_xbj;
   Double_t        ev_Q2;
   Double_t        ev_th;
   Double_t        ev_ph;
   Double_t        ev_Aperp;
   Double_t        ev_Apar;
   Double_t        ev_Pt;
   Double_t        ev_Pl;
   Double_t        ev_vx;
   Double_t        ev_vy;
   Double_t        ev_vz;
   Double_t        ev_ep;
   Double_t        ev_np;
   Double_t        ev_epx;
   Double_t        ev_epy;
   Double_t        ev_epz;
   Double_t        ev_npx;
   Double_t        ev_npy;
   Double_t        ev_npz;
   Double_t        ev_nth;
   Double_t        ev_nph;
   Double_t        ev_pmperp;
   Double_t        ev_pmpar;
   Double_t        ev_pmparsm;
   Double_t        ev_z;
   Double_t        ev_phperp;
   Double_t        ev_phih;
   Double_t        ev_phiS;
   Double_t        ev_MX2;
   Double_t        ev_Sx;
   Double_t        ev_Sy;
   Double_t        ev_Sz;
   Int_t           ev_nucl;
   Int_t           ev_fnucl;
   Int_t           ev_hadr;
   Int_t           ev_earmaccept;
   Int_t           ev_harmaccept;
   /**/
   // TODO: do some cleaning in here: I don't think we want ecal box data structures.
   // GEM variables
   // std::vector<gem_branch> GEMs;

   //BB GEMs variables
   TSBSGeant4::GEMData_t Earm_BBGEM;
   TSBSGeant4::DigGEMData_t Earm_BBGEM_Dig;
   TSBSGeant4::TrackerData_t Earm_BBGEM_Track;

   // BB timing hodoscope
   TSBSGeant4::CalData_t Earm_BBHodoScint;
   TSBSGeant4::DigTimingData_t Earm_BBHodo_Dig;

   // GRINCH variables
   TSBSGeant4::RICHData_t Earm_GRINCH;
   TSBSGeant4::DigTimingData_t Earm_GRINCH_Dig;

   //BB ECal variables
   //TSBSGeant4::ECalData_t Earm_BBPS;
   TSBSGeant4::CalData_t Earm_BBPSTF1;
   TSBSGeant4::DigCalData_t Earm_BBPS_Dig;
   //TSBSGeant4::ECalData_t Earm_BBSH;
   TSBSGeant4::CalData_t Earm_BBSHTF1;
   TSBSGeant4::DigCalData_t Earm_BBSH_Dig;

   // Coordinate detector hits
   //TSBSGeant4::ECalData_t Earm_CDET;
   //TSBSGeant4::CalData_t  Earm_CDET_Scint;
   //TSBSGeant4::ECalData_t Harm_CDET;
   //TSBSGeant4::CalData_t  Harm_CDET_Scint;
   //exception: the only det that can be used on either side
   TSBSGeant4::CalData_t CDET_Scint;
   TSBSGeant4::DigTimingData_t CDET_Dig;
  
   // GEp Electromagnetic calorimeter hits
   // TSBSGeant4::CalData_t  Earm_ECAL_box;
   //TSBSGeant4::ECalData_t Earm_ECAL;
   TSBSGeant4::CalData_t Earm_ECalTF1;
   TSBSGeant4::DigCalData_t Earm_ECal_Dig;
   
   // Focal Plane Polarimeter 1 hits
   TSBSGeant4::GEMData_t Harm_FPP1;
   TSBSGeant4::DigGEMData_t Harm_FPP1_Dig;
   TSBSGeant4::TrackerData_t Harm_FPP1_Track;
   TSBSGeant4::GEMData_t Harm_FPP2;
   TSBSGeant4::DigGEMData_t Harm_FPP2_Dig;
   TSBSGeant4::TrackerData_t Harm_FPP2_Track;
   TSBSGeant4::GEMData_t Harm_FT;
   TSBSGeant4::DigGEMData_t Harm_FT_Dig;
   TSBSGeant4::TrackerData_t Harm_FT_Track;
   
   // Hadronic calorimeter hits
   // An example for how to simplify tree objects
   // TODO: Don't hard code detectors here, but rather read them in
   // through a database if possible
   // TSBSGeant4::CalData_t      hcalbox;
   TSBSGeant4::CalData_t Harm_HCalScint;
   TSBSGeant4::DigSampCalData_t Harm_HCal_Dig;
   //TSBSGeant4::ECalData_t     hcal;
   //TSBSGeant4::ECalPartData_t hcalpart;

   // GEn-RP Active analyzer hits
   TSBSGeant4::CalData_t Harm_ActAnScint;
   TSBSGeant4::DigTimingData_t Harm_ActAn_Dig;
   
   // GEn-RP PR polarimeter Scintillators hits;
   TSBSGeant4::CalData_t Harm_PRPolScintBeamSide;
   TSBSGeant4::DigTimingData_t Harm_PRPolScintBeamSide_Dig;
   TSBSGeant4::CalData_t Harm_PRPolScintFarSide;
   TSBSGeant4::DigTimingData_t Harm_PRPolScintFarSide_Dig;
   
   // GEn-RP PR polarimeter Scintillators hits;
   TSBSGeant4::GEMData_t Harm_CEPolFront;
   TSBSGeant4::DigGEMData_t Harm_CEPolFront_Dig;
   TSBSGeant4::TrackerData_t Harm_CEPolFront_Track;
   TSBSGeant4::GEMData_t Harm_CEPolRear;
   TSBSGeant4::DigGEMData_t Harm_CEPolRear_Dig;
   TSBSGeant4::TrackerData_t Harm_CEPolRear_Track;
   
   TSBSGeant4::GEMData_t Harm_PrPolGEMBeamSide;
   TSBSGeant4::DigGEMData_t Harm_PrPolGEMBeamSide_Dig;
   TSBSGeant4::TrackerData_t Harm_PrPolGEMBeamSide_Track;
   TSBSGeant4::GEMData_t Harm_PrPolGEMFarSide;
   TSBSGeant4::DigGEMData_t Harm_PrPolGEMFarSide_Dig;
   TSBSGeant4::TrackerData_t Harm_PrPolGEMFarSide_Track;
   
   //SBS GEMs variables
   TSBSGeant4::GEMData_t Harm_SBSGEM;
   TSBSGeant4::DigGEMData_t Harm_SBSGEM_Dig;
   TSBSGeant4::TrackerData_t Harm_SBSGEM_Track;

   // RICH variables
   TSBSGeant4::RICHData_t Harm_RICH;
   TSBSGeant4::DigTimingData_t Harm_RICH_Dig;

   //Pythia variables
   Double_t              primaries_Sigma;
   Double_t              primaries_Ebeam;
   Double_t              primaries_Eprime;
   Double_t              primaries_theta_e;
   Double_t              primaries_phi_e;
   Double_t              primaries_px_e;
   Double_t              primaries_py_e;
   Double_t              primaries_pz_e;
   Double_t              primaries_vx_e;
   Double_t              primaries_vy_e;
   Double_t              primaries_vz_e;
   Double_t              primaries_Egamma;
   Double_t              primaries_theta_gamma;
   Double_t              primaries_phi_gamma;
   Double_t              primaries_px_gamma;
   Double_t              primaries_py_gamma;
   Double_t              primaries_pz_gamma;
   Double_t              primaries_vx_gamma;
   Double_t              primaries_vy_gamma;
   Double_t              primaries_vz_gamma;
   
   Int_t                 Primaries_Nprimaries;
   std::vector<int>     *Primaries_PID;
   std::vector<int>     *Primaries_genflag;
   std::vector<double>  *Primaries_Px;
   std::vector<double>  *Primaries_Py;
   std::vector<double>  *Primaries_Pz;
   std::vector<double>  *Primaries_vx;
   std::vector<double>  *Primaries_vy;
   std::vector<double>  *Primaries_vz;
   std::vector<double>  *Primaries_M;
   std::vector<double>  *Primaries_E;
   std::vector<double>  *Primaries_P;
   std::vector<double>  *Primaries_t;
   std::vector<double>  *Primaries_theta;
   std::vector<double>  *Primaries_phi;
   
   // List of branches

   TBranch        *b_ev;   //!
   //TBranch        *b_gen;   //!
   /*
   TBranch        *b_primaries_Sigma;   //!
   TBranch        *b_primaries_Ebeam;   //!
   TBranch        *b_primaries_Eprime;   //!
   TBranch        *b_primaries_theta_e;   //!
   TBranch        *b_primaries_phi_e;   //!
   TBranch        *b_primaries_px_e;   //!
   TBranch        *b_primaries_py_e;   //!
   TBranch        *b_primaries_pz_e;   //!
   TBranch        *b_primaries_vx_e;   //!
   TBranch        *b_primaries_vy_e;   //!
   TBranch        *b_primaries_vz_e;   //!
   TBranch        *b_primaries_Egamma;   //!
   TBranch        *b_primaries_theta_gamma;   //!
   TBranch        *b_primaries_phi_gamma;   //!
   TBranch        *b_primaries_px_gamma;   //!
   TBranch        *b_primaries_py_gamma;   //!
   TBranch        *b_primaries_pz_gamma;   //!
   TBranch        *b_primaries_vx_gamma;   //!
   TBranch        *b_primaries_vy_gamma;   //!
   TBranch        *b_primaries_vz_gamma;   //!
   
   TBranch        *b_Primaries_Nprimaries;   //!
   TBranch        *b_Primaries_PID;   //!
   TBranch        *b_Primaries_genflag;   //!
   TBranch        *b_Primaries_Px;   //!
   TBranch        *b_Primaries_Py;   //!
   TBranch        *b_Primaries_Pz;   //!
   TBranch        *b_Primaries_vx;   //!
   TBranch        *b_Primaries_vy;   //!
   TBranch        *b_Primaries_vz;   //!
   TBranch        *b_Primaries_M;   //!
   TBranch        *b_Primaries_E;   //!
   TBranch        *b_Primaries_P;   //!
   TBranch        *b_Primaries_t;   //!
   TBranch        *b_Primaries_theta;   //!
   TBranch        *b_Primaries_phi;   //!
   */

   g4sbs_tree(TTree *tree);//, std::vector<TString> det_list);
   //g4sbs_tree(TTree *tree=0, Exp_t expt = kGMN, bool pythia = false);
   //, bool ecalbox = false, bool hcalbox = false);
   // EPAF: We need to clean this. 
   virtual ~g4sbs_tree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntries();
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);//, std::vector<TString> det_list);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   
   // void ClearDigBranches();
   // void FillDigBranches();
   
protected:
   void SetupDetBranch(TSBSGeant4::VDetData_t &det, const char* prefix);
};

#endif



