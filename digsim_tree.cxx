//#define digsim_tree_cxx
#include "digsim_tree.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

digsim_tree::digsim_tree(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("simdig_test.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("simdig_test.root");
      }
      f->GetObject("digtree",tree);

   }
   Init(tree);
}

digsim_tree::~digsim_tree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t digsim_tree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}

Long64_t digsim_tree::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void digsim_tree::Init(TTree *tree)
{
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the branch addresses and branch
  // pointers of the tree will be set.
  // It is normally not necessary to make changes to the generated
  // code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF
  // (once per file to be processed).
  
  // Set branch addresses and branch pointers
  if (!tree) return;
  fChain = tree;
  fCurrent = -1;
  
  fChain->SetBranchAddress("RunID", &RunID, &b_RunID);
  fChain->SetBranchAddress("EvtID", &EvtID, &b_EvtID);
  fChain->SetBranchAddress("Weight", &Weight, &b_Weight);
  fChain->SetBranchAddress("NSignal", &NSignal, &b_NSignal);
  
  // HCal
  // sbs_hcal_simhits = new PMTSimHit_t(true, true);
  // SetupDetBranch(sbs_hcal_simhits, "sbs.hcal.simhit");
  // sbs_hcal_hits = new SampHitData_t(true);
  // SetupDetBranch(sbs_hcal_hits, "sbs.hcal.hit");
  PMTSimHitDet["sbs.hcal"] = new PMTSimHit_t(true, true);
  SetupDetBranch(PMTSimHitDet["sbs.hcal"], "sbs.hcal.simhit");
  SampHitDataDet["sbs.hcal"] = new SampHitData_t(true);
  SetupDetBranch(SampHitDataDet["sbs.hcal"], "sbs.hcal.hit");
  

  // PS/SH
  // bb_sh_simhits = new PMTSimHit_t(false, true);
  // SetupDetBranch(bb_sh_simhits, "bb.sh.simhit");
  // bb_sh_hits = new HitData_t(true, false);
  // SetupDetBranch(bb_sh_hits, "bb.sh.hit");
  PMTSimHitDet["bb.sh"] = new PMTSimHit_t(false, true);
  SetupDetBranch(PMTSimHitDet["bb.sh"], "bb.sh.simhit");
  HitDataDet["bb.sh"] = new HitData_t(true, false);
  SetupDetBranch(SampHitDataDet["bb.sh"], "bb.sh.hit");
  // bb_ps_simhits = new PMTSimHit_t(false, true);
  // SetupDetBranch(bb_ps_simhits, "bb.ps.simhit");
  // bb_ps_hits = new HitData_t(true, false);
  // SetupDetBranch(bb_ps_hits, "bb.ps.hit");
  PMTSimHitDet["bb.ps"] = new PMTSimHit_t(false, true);
  SetupDetBranch(PMTSimHitDet["bb.ps"], "bb.ps.simhit");
  HitDataDet["bb.ps"] = new HitData_t(true, false);
  SetupDetBranch(SampHitDataDet["bb.ps"], "bb.ps.hit");
  
  // Hodoscope
  // bb_hodo_simhits = new PMTSimHit_t(true, true);
  // SetupDetBranch(bb_hodo_simhits, "bb.hodo.simhit");
  // bb_hodo_hits = new HitData_t(false, true);
  // SetupDetBranch(bb_hodo_hits, "bb.hodo.hit");
  PMTSimHitDet["bb.hodo"] = new PMTSimHit_t(true, true);
  SetupDetBranch(PMTSimHitDet["bb.hodo"], "bb.hodo.simhit");
  HitDataDet["bb.hodo"] = new HitData_t(false, true);
  SetupDetBranch(SampHitDataDet["bb.hodo"], "bb.hodo.hit");
  
  // Grinch
  // bb_grinch_simhits = new PMTSimHit_t(true, false);
  // SetupDetBranch(bb_hodo_simhits, "bb.grinch.simhit");
  // bb_grinch_hits = new HitData_t(false, true);
  // SetupDetBranch(bb_hodo_hits, "bb.grinch.hit");
  PMTSimHitDet["bb.grinch"] = new PMTSimHit_t(true, false);
  SetupDetBranch(PMTSimHitDet["bb.grinch"], "bb.grinch.simhit");
  HitDataDet["bb.grinch"] = new HitData_t(false, true);
  SetupDetBranch(SampHitDataDet["bb.grinch"], "bb.grinch.hit");
  
  // GEMs
  // bb_gem_simhits = new GEMSimHit_t();
  // SetupDetBranch(bb_hodo_simhits, "bb.gem.simhit");
  // bb_gem_hits = new GEMData_t();
  // SetupDetBranch(bb_hodo_hits, "bb.gem.hit");
  GEMSimHitDet["bb.gem"] = new GEMSimHit_t();
  SetupDetBranch(GEMSimHitDet["bb.gem"], "bb.gem.simhit");
  GEMDataDet["bb.gem"] = new GEMData_t();
  SetupDetBranch(GEMDataDet["bb.gem"], "bb.gem.hit");
  
}

Bool_t digsim_tree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void digsim_tree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t digsim_tree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

void digsim_tree::Loop()
{
//   In a ROOT session, you can do:
//      root> .L digsim_tree.C
//      root> digsim_tree t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
   }
}

//void digsim_tree::SetupDetBranch(SBSDigSim::VDetData_t &det, const char *prefix)
void digsim_tree::SetupDetBranch(SBSDigSim::VDetData_t* det, const char *prefix)
{
  det->SetupBranches(fChain,prefix);
}
