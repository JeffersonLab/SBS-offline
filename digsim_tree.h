//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Feb 19 07:46:51 2020 by ROOT version 6.14/04
// from TTree digtree/
// found on file: simdig_test.root
//////////////////////////////////////////////////////////

#ifndef digsim_tree_h
#define digsim_tree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <digsim_data.h>
#include <map>

#define treeName "digtree"

// Header file for the classes stored in the TTree if any.
using namespace SBSDigSim;

class digsim_tree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           RunID;
   Int_t           EvtID;
   Double_t        Weight;
   Int_t           NSignal;
   // List of branches
   TBranch        *b_RunID;   //!
   TBranch        *b_EvtID;   //!
   TBranch        *b_Weight;   //!
   TBranch        *b_NSignal;   //!
   
   std::map<std::string, MCTrack_t*> MCTrack;
   std::map<std::string, TrackMCHit_t*> TrackMCHitDet;
   std::map<std::string, PMTSimHit_t*> PMTSimHitDet;
   std::map<std::string, GEMSimHit_t*> GEMSimHitDet;
   
   //std::map<std::string, HitData_t*> HitDataDet;
   //std::map<std::string, SampHitData_t*> SampHitDataDet;
   std::map<std::string, UHitData_t*> HitDataDet;
   
   
   //Need to declare and fill standard containers
   //std::map<string, vector<> >

   digsim_tree(TTree *tree=0);
   virtual ~digsim_tree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

protected:
   //void SetupDetBranch(SBSDigSim::VDetData_t &det, const char* prefix);
   void SetupDetBranch(SBSDigSim::VDetData_t* det, const char* prefix);
};

#endif
