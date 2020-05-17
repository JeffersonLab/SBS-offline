//#define digsim_tree_cxx
#include "digsim_tree.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#define DEBUG 0

digsim_tree::digsim_tree(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("simdig_test.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("simdig_test.root");
      }
      f->GetObject(treeName,tree);

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
  // TrackMCHitDet["sbs.hcal"] = new TrackMCHit_t();
  // SetupDetBranch(TrackMCHitDet["sbs.hcal"], "sbs.hcal.trackmchit");
  PMTSimHitDet["sbs.hcal"] = new PMTSimHit_t(true, true);
  SetupDetBranch(PMTSimHitDet["sbs.hcal"], "sbs.hcal.simhit");
  //SampHitDataDet["sbs.hcal"] = new SampHitData_t(true);
  //SetupDetBranch(SampHitDataDet["sbs.hcal"], "sbs.hcal.hit");
  HitDataDet["sbs.hcal"] = new UHitData_t(true, true, true);
  SetupDetBranch(HitDataDet["sbs.hcal"], "sbs.hcal.hit");
  
  
  // PS/SH
  // TrackMCHitDet["bb.sh"] = new TrackMCHit_t();
  // SetupDetBranch(TrackMCHitDet["bb.sh"], "bb.sh.trackmchit");
  PMTSimHitDet["bb.sh"] = new PMTSimHit_t(false, true);
  SetupDetBranch(PMTSimHitDet["bb.sh"], "bb.sh.simhit");
  //HitDataDet["bb.sh"] = new HitData_t(true, false);
  HitDataDet["bb.sh"] = new UHitData_t(true, false, false);
  SetupDetBranch(HitDataDet["bb.sh"], "bb.sh.hit");
  
  // TrackMCHitDet["bb.ps"] = new TrackMCHit_t();
  // SetupDetBranch(TrackMCHitDet["bb.ps"], "bb.ps.trackmchit");
  PMTSimHitDet["bb.ps"] = new PMTSimHit_t(false, true);
  SetupDetBranch(PMTSimHitDet["bb.ps"], "bb.ps.simhit");
  //HitDataDet["bb.ps"] = new HitData_t(true, false);
  HitDataDet["bb.ps"] = new UHitData_t(true, false, false);
  SetupDetBranch(HitDataDet["bb.ps"], "bb.ps.hit");
  
  // Hodoscope
  // TrackMCHitDet["bb.hodo"] = new TrackMCHit_t();
  // SetupDetBranch(TrackMCHitDet["bb.hodo"], "bb.hodo.trackmchit");
  PMTSimHitDet["bb.hodo"] = new PMTSimHit_t(true, true);
  SetupDetBranch(PMTSimHitDet["bb.hodo"], "bb.hodo.simhit");
  //HitDataDet["bb.hodo"] = new HitData_t(false, true);
  HitDataDet["bb.hodo"] = new UHitData_t(false, true, false);
  SetupDetBranch(HitDataDet["bb.hodo"], "bb.hodo.hit");
  
  // Grinch
  // TrackMCHitDet["bb.grinch"] = new TrackMCHit_t();
  // SetupDetBranch(TrackMCHitDet["bb.grinch"], "bb.grinch.trackmchit");
  PMTSimHitDet["bb.grinch"] = new PMTSimHit_t(true, false);
  SetupDetBranch(PMTSimHitDet["bb.grinch"], "bb.grinch.simhit");
  //HitDataDet["bb.grinch"] = new HitData_t(false, true);
  HitDataDet["bb.grinch"] = new UHitData_t(false, true, false);
  SetupDetBranch(HitDataDet["bb.grinch"], "bb.grinch.hit");
  
  // GEMs
  // MCTrack["bb.gem"] = new MCTrack_t();
  // SetupDetBranch(MCTrack["bb.gem"], "bb.gem.mctrack");
  GEMSimHitDet["bb.gem"] = new GEMSimHit_t();
  SetupDetBranch(GEMSimHitDet["bb.gem"], "bb.gem.simhit");
  std::string fullgemname;
  for(int ipl = 0; ipl<5; ipl++){
    // int nmod = 3;
    // if(ipl==4)nmod = 4;
    // for(int imod = 0; imod<nmod; imod++){
    for(int ipr = 0; ipr<2; ipr++){
      //fullgemname = Form("bb.gem.p%d.m%d.%s", 
      //		   ipl+1, imod+1, kProj_str[ipr].c_str());
      fullgemname = Form("bb.gem.%d.%s", ipl+1, kProj_str[ipr].c_str());
      //SampHitDataDet[fullgemname] = new SampHitData_t(false);
      //SetupDetBranch(SampHitDataDet[fullgemname], 
      HitDataDet[fullgemname] = new UHitData_t(true, false, true);
      SetupDetBranch(HitDataDet[fullgemname], 
		     Form("%s.hit", fullgemname.c_str()));
    }
    //}
  }
  
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
#if DEBUG>0
  printf("SetupDetBranch %s\n", prefix);
#endif
  det->SetupBranches(fChain,prefix);
}

