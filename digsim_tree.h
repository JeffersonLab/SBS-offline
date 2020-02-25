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
   
   //just GMn... 
   // we'll see later to make it configurable
   PMTSimHit_t *sbs_hcal_simhits;
   PMTSimHit_t *bb_sh_simhits;
   PMTSimHit_t *bb_ps_simhits;
   PMTSimHit_t *bb_hodo_simhits;
   PMTSimHit_t *bb_grinch_simhits;
   GEMSimHit_t *bb_gem_simhits;
   
   SampHitData_t *sbs_hcal_hits;
   HitData_t *bb_sh_hits;
   HitData_t *bb_ps_hits;
   HitData_t *bb_hodo_hits;
   HitData_t *bb_grinch_hits;
   GEMData_t *bb_gem_hits;
   
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

/*

#ifdef digsim_tree_cxx
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

   // Set object pointer
   sbs_hcal_simhit_src = 0;
   sbs_hcal_simhit_trid = 0;
   sbs_hcal_simhit_pid = 0;
   sbs_hcal_simhit_chan = 0;
   sbs_hcal_simhit_edep = 0;
   sbs_hcal_simhit_npe = 0;
   sbs_hcal_simhit_time = 0;
   sbs_hcal_simhit_t_lead = 0;
   sbs_hcal_simhit_t_trail = 0;
   sbs_hcal_hit_chan = 0;
   sbs_hcal_hit_nwords = 0;
   sbs_hcal_hit_adcsum = 0;
   sbs_hcal_hit_samps_adc = 0;
   sbs_hcal_hit_samps_datawords = 0;
   sbs_hcal_hit_tdc_l = 0;
   sbs_hcal_hit_tdc_t = 0;
   bb_sh_simhit_src = 0;
   bb_sh_simhit_trid = 0;
   bb_sh_simhit_pid = 0;
   bb_sh_simhit_chan = 0;
   bb_sh_simhit_edep = 0;
   bb_sh_simhit_npe = 0;
   bb_sh_simhit_time = 0;
   bb_sh_hit_chan = 0;
   bb_sh_hit_dataword = 0;
   bb_sh_hit_adc = 0;
   bb_ps_simhit_src = 0;
   bb_ps_simhit_trid = 0;
   bb_ps_simhit_pid = 0;
   bb_ps_simhit_chan = 0;
   bb_ps_simhit_edep = 0;
   bb_ps_simhit_npe = 0;
   bb_ps_simhit_time = 0;
   bb_ps_hit_chan = 0;
   bb_ps_hit_dataword = 0;
   bb_ps_hit_adc = 0;
   bb_hodo_simhit_src = 0;
   bb_hodo_simhit_trid = 0;
   bb_hodo_simhit_pid = 0;
   bb_hodo_simhit_chan = 0;
   bb_hodo_simhit_edep = 0;
   bb_hodo_simhit_npe = 0;
   bb_hodo_simhit_time = 0;
   bb_hodo_simhit_t_lead = 0;
   bb_hodo_simhit_t_trail = 0;
   bb_hodo_hit_chan = 0;
   bb_hodo_hit_dataword = 0;
   bb_hodo_hit_tdc_l = 0;
   bb_hodo_hit_tdc_t = 0;
   bb_grinch_simhit_src = 0;
   bb_grinch_simhit_trid = 0;
   bb_grinch_simhit_pid = 0;
   bb_grinch_simhit_chan = 0;
   bb_grinch_simhit_npe = 0;
   bb_grinch_simhit_time = 0;
   bb_grinch_simhit_t_lead = 0;
   bb_grinch_simhit_t_trail = 0;
   bb_grinch_hit_chan = 0;
   bb_grinch_hit_dataword = 0;
   bb_grinch_hit_tdc_l = 0;
   bb_grinch_hit_tdc_t = 0;
   bb_gem_simhit_src = 0;
   bb_gem_simhit_trid = 0;
   bb_gem_simhit_pid = 0;
   bb_gem_simhit_plane = 0;
   bb_gem_simhit_module = 0;
   bb_gem_simhit_edep = 0;
   bb_gem_simhit_time = 0;
   bb_gem_simhit_xpos = 0;
   bb_gem_simhit_ypos = 0;
   bb_gem_simhit_px = 0;
   bb_gem_simhit_py = 0;
   bb_gem_simhit_pz = 0;
   bb_gem_simhit_sizex = 0;
   bb_gem_simhit_sizey = 0;
   bb_gem_simhit_startx = 0;
   bb_gem_simhit_starty = 0;
   bb_gem_hit_plane = 0;
   bb_gem_hit_module = 0;
   bb_gem_hit_proj = 0;
   bb_gem_hit_nwords = 0;
   bb_gem_hit_strip = 0;
   bb_gem_hit_samp = 0;
   bb_gem_hit_samps_adc = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("RunID", &RunID, &b_RunID);
   fChain->SetBranchAddress("EvtID", &EvtID, &b_EvtID);
   fChain->SetBranchAddress("Weight", &Weight, &b_Weight);
   fChain->SetBranchAddress("NSignal", &NSignal, &b_NSignal);
   fChain->SetBranchAddress("sbs.hcal.simhit.nhits", &sbs_hcal_simhit_nhits, &b_sbs_hcal_simhit_nhits);
   fChain->SetBranchAddress("sbs.hcal.simhit.src", &sbs_hcal_simhit_src, &b_sbs_hcal_simhit_src);
   fChain->SetBranchAddress("sbs.hcal.simhit.trid", &sbs_hcal_simhit_trid, &b_sbs_hcal_simhit_trid);
   fChain->SetBranchAddress("sbs.hcal.simhit.pid", &sbs_hcal_simhit_pid, &b_sbs_hcal_simhit_pid);
   fChain->SetBranchAddress("sbs.hcal.simhit.chan", &sbs_hcal_simhit_chan, &b_sbs_hcal_simhit_chan);
   fChain->SetBranchAddress("sbs.hcal.simhit.edep", &sbs_hcal_simhit_edep, &b_sbs_hcal_simhit_edep);
   fChain->SetBranchAddress("sbs.hcal.simhit.npe", &sbs_hcal_simhit_npe, &b_sbs_hcal_simhit_npe);
   fChain->SetBranchAddress("sbs.hcal.simhit.time", &sbs_hcal_simhit_time, &b_sbs_hcal_simhit_time);
   fChain->SetBranchAddress("sbs.hcal.simhit.t_lead", &sbs_hcal_simhit_t_lead, &b_sbs_hcal_simhit_t_lead);
   fChain->SetBranchAddress("sbs.hcal.simhit.t_trail", &sbs_hcal_simhit_t_trail, &b_sbs_hcal_simhit_t_trail);
   fChain->SetBranchAddress("sbs.hcal.hit.nhits", &sbs_hcal_hit_nhits, &b_sbs_hcal_hit_nhits);
   fChain->SetBranchAddress("sbs.hcal.hit.chan", &sbs_hcal_hit_chan, &b_sbs_hcal_hit_chan);
   fChain->SetBranchAddress("sbs.hcal.hit.nwords", &sbs_hcal_hit_nwords, &b_sbs_hcal_hit_nwords);
   fChain->SetBranchAddress("sbs.hcal.hit.adcsum", &sbs_hcal_hit_adcsum, &b_sbs_hcal_hit_adcsum);
   fChain->SetBranchAddress("sbs.hcal.hit.samps_adc", &sbs_hcal_hit_samps_adc, &b_sbs_hcal_hit_samps_adc);
   fChain->SetBranchAddress("sbs.hcal.hit.samps_datawords", &sbs_hcal_hit_samps_datawords, &b_sbs_hcal_hit_samps_datawords);
   fChain->SetBranchAddress("sbs.hcal.hit.tdc_l", &sbs_hcal_hit_tdc_l, &b_sbs_hcal_hit_tdc_l);
   fChain->SetBranchAddress("sbs.hcal.hit.tdc_t", &sbs_hcal_hit_tdc_t, &b_sbs_hcal_hit_tdc_t);
   fChain->SetBranchAddress("bb.sh.simhit.nhits", &bb_sh_simhit_nhits, &b_bb_sh_simhit_nhits);
   fChain->SetBranchAddress("bb.sh.simhit.src", &bb_sh_simhit_src, &b_bb_sh_simhit_src);
   fChain->SetBranchAddress("bb.sh.simhit.trid", &bb_sh_simhit_trid, &b_bb_sh_simhit_trid);
   fChain->SetBranchAddress("bb.sh.simhit.pid", &bb_sh_simhit_pid, &b_bb_sh_simhit_pid);
   fChain->SetBranchAddress("bb.sh.simhit.chan", &bb_sh_simhit_chan, &b_bb_sh_simhit_chan);
   fChain->SetBranchAddress("bb.sh.simhit.edep", &bb_sh_simhit_edep, &b_bb_sh_simhit_edep);
   fChain->SetBranchAddress("bb.sh.simhit.npe", &bb_sh_simhit_npe, &b_bb_sh_simhit_npe);
   fChain->SetBranchAddress("bb.sh.simhit.time", &bb_sh_simhit_time, &b_bb_sh_simhit_time);
   fChain->SetBranchAddress("bb.sh.hit.nhits", &bb_sh_hit_nhits, &b_bb_sh_hit_nhits);
   fChain->SetBranchAddress("bb.sh.hit.chan", &bb_sh_hit_chan, &b_bb_sh_hit_chan);
   fChain->SetBranchAddress("bb.sh.hit.dataword", &bb_sh_hit_dataword, &b_bb_sh_hit_dataword);
   fChain->SetBranchAddress("bb.sh.hit.adc", &bb_sh_hit_adc, &b_bb_sh_hit_adc);
   fChain->SetBranchAddress("bb.ps.simhit.nhits", &bb_ps_simhit_nhits, &b_bb_ps_simhit_nhits);
   fChain->SetBranchAddress("bb.ps.simhit.src", &bb_ps_simhit_src, &b_bb_ps_simhit_src);
   fChain->SetBranchAddress("bb.ps.simhit.trid", &bb_ps_simhit_trid, &b_bb_ps_simhit_trid);
   fChain->SetBranchAddress("bb.ps.simhit.pid", &bb_ps_simhit_pid, &b_bb_ps_simhit_pid);
   fChain->SetBranchAddress("bb.ps.simhit.chan", &bb_ps_simhit_chan, &b_bb_ps_simhit_chan);
   fChain->SetBranchAddress("bb.ps.simhit.edep", &bb_ps_simhit_edep, &b_bb_ps_simhit_edep);
   fChain->SetBranchAddress("bb.ps.simhit.npe", &bb_ps_simhit_npe, &b_bb_ps_simhit_npe);
   fChain->SetBranchAddress("bb.ps.simhit.time", &bb_ps_simhit_time, &b_bb_ps_simhit_time);
   fChain->SetBranchAddress("bb.ps.hit.nhits", &bb_ps_hit_nhits, &b_bb_ps_hit_nhits);
   fChain->SetBranchAddress("bb.ps.hit.chan", &bb_ps_hit_chan, &b_bb_ps_hit_chan);
   fChain->SetBranchAddress("bb.ps.hit.dataword", &bb_ps_hit_dataword, &b_bb_ps_hit_dataword);
   fChain->SetBranchAddress("bb.ps.hit.adc", &bb_ps_hit_adc, &b_bb_ps_hit_adc);
   fChain->SetBranchAddress("bb.hodo.simhit.nhits", &bb_hodo_simhit_nhits, &b_bb_hodo_simhit_nhits);
   fChain->SetBranchAddress("bb.hodo.simhit.src", &bb_hodo_simhit_src, &b_bb_hodo_simhit_src);
   fChain->SetBranchAddress("bb.hodo.simhit.trid", &bb_hodo_simhit_trid, &b_bb_hodo_simhit_trid);
   fChain->SetBranchAddress("bb.hodo.simhit.pid", &bb_hodo_simhit_pid, &b_bb_hodo_simhit_pid);
   fChain->SetBranchAddress("bb.hodo.simhit.chan", &bb_hodo_simhit_chan, &b_bb_hodo_simhit_chan);
   fChain->SetBranchAddress("bb.hodo.simhit.edep", &bb_hodo_simhit_edep, &b_bb_hodo_simhit_edep);
   fChain->SetBranchAddress("bb.hodo.simhit.npe", &bb_hodo_simhit_npe, &b_bb_hodo_simhit_npe);
   fChain->SetBranchAddress("bb.hodo.simhit.time", &bb_hodo_simhit_time, &b_bb_hodo_simhit_time);
   fChain->SetBranchAddress("bb.hodo.simhit.t_lead", &bb_hodo_simhit_t_lead, &b_bb_hodo_simhit_t_lead);
   fChain->SetBranchAddress("bb.hodo.simhit.t_trail", &bb_hodo_simhit_t_trail, &b_bb_hodo_simhit_t_trail);
   fChain->SetBranchAddress("bb.hodo.hit.nhits", &bb_hodo_hit_nhits, &b_bb_hodo_hit_nhits);
   fChain->SetBranchAddress("bb.hodo.hit.chan", &bb_hodo_hit_chan, &b_bb_hodo_hit_chan);
   fChain->SetBranchAddress("bb.hodo.hit.dataword", &bb_hodo_hit_dataword, &b_bb_hodo_hit_dataword);
   fChain->SetBranchAddress("bb.hodo.hit.tdc_l", &bb_hodo_hit_tdc_l, &b_bb_hodo_hit_tdc_l);
   fChain->SetBranchAddress("bb.hodo.hit.tdc_t", &bb_hodo_hit_tdc_t, &b_bb_hodo_hit_tdc_t);
   fChain->SetBranchAddress("bb.grinch.simhit.nhits", &bb_grinch_simhit_nhits, &b_bb_grinch_simhit_nhits);
   fChain->SetBranchAddress("bb.grinch.simhit.src", &bb_grinch_simhit_src, &b_bb_grinch_simhit_src);
   fChain->SetBranchAddress("bb.grinch.simhit.trid", &bb_grinch_simhit_trid, &b_bb_grinch_simhit_trid);
   fChain->SetBranchAddress("bb.grinch.simhit.pid", &bb_grinch_simhit_pid, &b_bb_grinch_simhit_pid);
   fChain->SetBranchAddress("bb.grinch.simhit.chan", &bb_grinch_simhit_chan, &b_bb_grinch_simhit_chan);
   fChain->SetBranchAddress("bb.grinch.simhit.npe", &bb_grinch_simhit_npe, &b_bb_grinch_simhit_npe);
   fChain->SetBranchAddress("bb.grinch.simhit.time", &bb_grinch_simhit_time, &b_bb_grinch_simhit_time);
   fChain->SetBranchAddress("bb.grinch.simhit.t_lead", &bb_grinch_simhit_t_lead, &b_bb_grinch_simhit_t_lead);
   fChain->SetBranchAddress("bb.grinch.simhit.t_trail", &bb_grinch_simhit_t_trail, &b_bb_grinch_simhit_t_trail);
   fChain->SetBranchAddress("bb.grinch.hit.nhits", &bb_grinch_hit_nhits, &b_bb_grinch_hit_nhits);
   fChain->SetBranchAddress("bb.grinch.hit.chan", &bb_grinch_hit_chan, &b_bb_grinch_hit_chan);
   fChain->SetBranchAddress("bb.grinch.hit.dataword", &bb_grinch_hit_dataword, &b_bb_grinch_hit_dataword);
   fChain->SetBranchAddress("bb.grinch.hit.tdc_l", &bb_grinch_hit_tdc_l, &b_bb_grinch_hit_tdc_l);
   fChain->SetBranchAddress("bb.grinch.hit.tdc_t", &bb_grinch_hit_tdc_t, &b_bb_grinch_hit_tdc_t);
   fChain->SetBranchAddress("bb.gem.simhit.nhits", &bb_gem_simhit_nhits, &b_bb_gem_simhit_nhits);
   fChain->SetBranchAddress("bb.gem.simhit.src", &bb_gem_simhit_src, &b_bb_gem_simhit_src);
   fChain->SetBranchAddress("bb.gem.simhit.trid", &bb_gem_simhit_trid, &b_bb_gem_simhit_trid);
   fChain->SetBranchAddress("bb.gem.simhit.pid", &bb_gem_simhit_pid, &b_bb_gem_simhit_pid);
   fChain->SetBranchAddress("bb.gem.simhit.plane", &bb_gem_simhit_plane, &b_bb_gem_simhit_plane);
   fChain->SetBranchAddress("bb.gem.simhit.module", &bb_gem_simhit_module, &b_bb_gem_simhit_module);
   fChain->SetBranchAddress("bb.gem.simhit.edep", &bb_gem_simhit_edep, &b_bb_gem_simhit_edep);
   fChain->SetBranchAddress("bb.gem.simhit.time", &bb_gem_simhit_time, &b_bb_gem_simhit_time);
   fChain->SetBranchAddress("bb.gem.simhit.xpos", &bb_gem_simhit_xpos, &b_bb_gem_simhit_xpos);
   fChain->SetBranchAddress("bb.gem.simhit.ypos", &bb_gem_simhit_ypos, &b_bb_gem_simhit_ypos);
   fChain->SetBranchAddress("bb.gem.simhit.px", &bb_gem_simhit_px, &b_bb_gem_simhit_px);
   fChain->SetBranchAddress("bb.gem.simhit.py", &bb_gem_simhit_py, &b_bb_gem_simhit_py);
   fChain->SetBranchAddress("bb.gem.simhit.pz", &bb_gem_simhit_pz, &b_bb_gem_simhit_pz);
   fChain->SetBranchAddress("bb.gem.simhit.sizex", &bb_gem_simhit_sizex, &b_bb_gem_simhit_sizex);
   fChain->SetBranchAddress("bb.gem.simhit.sizey", &bb_gem_simhit_sizey, &b_bb_gem_simhit_sizey);
   fChain->SetBranchAddress("bb.gem.simhit.startx", &bb_gem_simhit_startx, &b_bb_gem_simhit_startx);
   fChain->SetBranchAddress("bb.gem.simhit.starty", &bb_gem_simhit_starty, &b_bb_gem_simhit_starty);
   fChain->SetBranchAddress("bb.gem.hit.nhits", &bb_gem_hit_nhits, &b_bb_gem_hit_nhits);
   fChain->SetBranchAddress("bb.gem.hit.plane", &bb_gem_hit_plane, &b_bb_gem_hit_plane);
   fChain->SetBranchAddress("bb.gem.hit.module", &bb_gem_hit_module, &b_bb_gem_hit_module);
   fChain->SetBranchAddress("bb.gem.hit.proj", &bb_gem_hit_proj, &b_bb_gem_hit_proj);
   fChain->SetBranchAddress("bb.gem.hit.nwords", &bb_gem_hit_nwords, &b_bb_gem_hit_nwords);
   fChain->SetBranchAddress("bb.gem.hit.strip", &bb_gem_hit_strip, &b_bb_gem_hit_strip);
   fChain->SetBranchAddress("bb.gem.hit.samp", &bb_gem_hit_samp, &b_bb_gem_hit_samp);
   fChain->SetBranchAddress("bb.gem.hit.samps_adc", &bb_gem_hit_samps_adc, &b_bb_gem_hit_samps_adc);
   Notify();
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
#endif // #ifdef digsim_tree_cxx
*/

   /*
   UInt_t          sbs_hcal_simhit_nhits;
   std::vector<short>   *sbs_hcal_simhit_src;
   std::vector<short>   *sbs_hcal_simhit_trid;
   std::vector<int>     *sbs_hcal_simhit_pid;
   std::vector<short>   *sbs_hcal_simhit_chan;
   std::vector<double>  *sbs_hcal_simhit_edep;
   std::vector<int>     *sbs_hcal_simhit_npe;
   std::vector<double>  *sbs_hcal_simhit_time;
   std::vector<double>  *sbs_hcal_simhit_t_lead;
   std::vector<double>  *sbs_hcal_simhit_t_trail;
   UInt_t          sbs_hcal_hit_nhits;
   std::vector<short>   *sbs_hcal_hit_chan;
   std::vector<unsigned int> *sbs_hcal_hit_nwords;
   std::vector<int>     *sbs_hcal_hit_adcsum;
   std::vector<std::vector<int> > *sbs_hcal_hit_samps_adc;
   std::vector<std::vector<unsigned int> > *sbs_hcal_hit_samps_datawords;
   std::vector<int>     *sbs_hcal_hit_tdc_l;
   std::vector<int>     *sbs_hcal_hit_tdc_t;
   UInt_t          bb_sh_simhit_nhits;
   std::vector<short>   *bb_sh_simhit_src;
   std::vector<short>   *bb_sh_simhit_trid;
   std::vector<int>     *bb_sh_simhit_pid;
   std::vector<short>   *bb_sh_simhit_chan;
   std::vector<double>  *bb_sh_simhit_edep;
   std::vector<int>     *bb_sh_simhit_npe;
   std::vector<double>  *bb_sh_simhit_time;
   UInt_t          bb_sh_hit_nhits;
   std::vector<short>   *bb_sh_hit_chan;
   std::vector<unsigned int> *bb_sh_hit_dataword;
   std::vector<int>     *bb_sh_hit_adc;
   UInt_t          bb_ps_simhit_nhits;
   std::vector<short>   *bb_ps_simhit_src;
   std::vector<short>   *bb_ps_simhit_trid;
   std::vector<int>     *bb_ps_simhit_pid;
   std::vector<short>   *bb_ps_simhit_chan;
   std::vector<double>  *bb_ps_simhit_edep;
   std::vector<int>     *bb_ps_simhit_npe;
   std::vector<double>  *bb_ps_simhit_time;
   UInt_t          bb_ps_hit_nhits;
   std::vector<short>   *bb_ps_hit_chan;
   std::vector<unsigned int> *bb_ps_hit_dataword;
   std::vector<int>     *bb_ps_hit_adc;
   UInt_t          bb_hodo_simhit_nhits;
   std::vector<short>   *bb_hodo_simhit_src;
   std::vector<short>   *bb_hodo_simhit_trid;
   std::vector<int>     *bb_hodo_simhit_pid;
   std::vector<short>   *bb_hodo_simhit_chan;
   std::vector<double>  *bb_hodo_simhit_edep;
   std::vector<int>     *bb_hodo_simhit_npe;
   std::vector<double>  *bb_hodo_simhit_time;
   std::vector<double>  *bb_hodo_simhit_t_lead;
   std::vector<double>  *bb_hodo_simhit_t_trail;
   UInt_t          bb_hodo_hit_nhits;
   std::vector<short>   *bb_hodo_hit_chan;
   std::vector<unsigned int> *bb_hodo_hit_dataword;
   std::vector<int>     *bb_hodo_hit_tdc_l;
   std::vector<int>     *bb_hodo_hit_tdc_t;
   UInt_t          bb_grinch_simhit_nhits;
   std::vector<short>   *bb_grinch_simhit_src;
   std::vector<short>   *bb_grinch_simhit_trid;
   std::vector<int>     *bb_grinch_simhit_pid;
   std::vector<short>   *bb_grinch_simhit_chan;
   std::vector<int>     *bb_grinch_simhit_npe;
   std::vector<double>  *bb_grinch_simhit_time;
   std::vector<double>  *bb_grinch_simhit_t_lead;
   std::vector<double>  *bb_grinch_simhit_t_trail;
   UInt_t          bb_grinch_hit_nhits;
   std::vector<short>   *bb_grinch_hit_chan;
   std::vector<unsigned int> *bb_grinch_hit_dataword;
   std::vector<int>     *bb_grinch_hit_tdc_l;
   std::vector<int>     *bb_grinch_hit_tdc_t;
   UInt_t          bb_gem_simhit_nhits;
   std::vector<short>   *bb_gem_simhit_src;
   std::vector<short>   *bb_gem_simhit_trid;
   std::vector<int>     *bb_gem_simhit_pid;
   std::vector<short>   *bb_gem_simhit_plane;
   std::vector<short>   *bb_gem_simhit_module;
   std::vector<double>  *bb_gem_simhit_edep;
   std::vector<double>  *bb_gem_simhit_time;
   std::vector<double>  *bb_gem_simhit_xpos;
   std::vector<double>  *bb_gem_simhit_ypos;
   std::vector<double>  *bb_gem_simhit_px;
   std::vector<double>  *bb_gem_simhit_py;
   std::vector<double>  *bb_gem_simhit_pz;
   std::vector<short>   *bb_gem_simhit_sizex;
   std::vector<short>   *bb_gem_simhit_sizey;
   std::vector<short>   *bb_gem_simhit_startx;
   std::vector<short>   *bb_gem_simhit_starty;
   UInt_t          bb_gem_hit_nhits;
   std::vector<short>   *bb_gem_hit_plane;
   std::vector<short>   *bb_gem_hit_module;
   std::vector<short>   *bb_gem_hit_proj;
   std::vector<unsigned int> *bb_gem_hit_nwords;
   std::vector<std::vector<short> > *bb_gem_hit_strip;
   std::vector<std::vector<short> > *bb_gem_hit_samp;
   std::vector<std::vector<int> > *bb_gem_hit_samps_adc;

   // List of branches
   TBranch        *b_RunID;   //!
   TBranch        *b_EvtID;   //!
   TBranch        *b_Weight;   //!
   TBranch        *b_NSignal;   //!
   TBranch        *b_sbs_hcal_simhit_nhits;   //!
   TBranch        *b_sbs_hcal_simhit_src;   //!
   TBranch        *b_sbs_hcal_simhit_trid;   //!
   TBranch        *b_sbs_hcal_simhit_pid;   //!
   TBranch        *b_sbs_hcal_simhit_chan;   //!
   TBranch        *b_sbs_hcal_simhit_edep;   //!
   TBranch        *b_sbs_hcal_simhit_npe;   //!
   TBranch        *b_sbs_hcal_simhit_time;   //!
   TBranch        *b_sbs_hcal_simhit_t_lead;   //!
   TBranch        *b_sbs_hcal_simhit_t_trail;   //!
   TBranch        *b_sbs_hcal_hit_nhits;   //!
   TBranch        *b_sbs_hcal_hit_chan;   //!
   TBranch        *b_sbs_hcal_hit_nwords;   //!
   TBranch        *b_sbs_hcal_hit_adcsum;   //!
   TBranch        *b_sbs_hcal_hit_samps_adc;   //!
   TBranch        *b_sbs_hcal_hit_samps_datawords;   //!
   TBranch        *b_sbs_hcal_hit_tdc_l;   //!
   TBranch        *b_sbs_hcal_hit_tdc_t;   //!
   TBranch        *b_bb_sh_simhit_nhits;   //!
   TBranch        *b_bb_sh_simhit_src;   //!
   TBranch        *b_bb_sh_simhit_trid;   //!
   TBranch        *b_bb_sh_simhit_pid;   //!
   TBranch        *b_bb_sh_simhit_chan;   //!
   TBranch        *b_bb_sh_simhit_edep;   //!
   TBranch        *b_bb_sh_simhit_npe;   //!
   TBranch        *b_bb_sh_simhit_time;   //!
   TBranch        *b_bb_sh_hit_nhits;   //!
   TBranch        *b_bb_sh_hit_chan;   //!
   TBranch        *b_bb_sh_hit_dataword;   //!
   TBranch        *b_bb_sh_hit_adc;   //!
   TBranch        *b_bb_ps_simhit_nhits;   //!
   TBranch        *b_bb_ps_simhit_src;   //!
   TBranch        *b_bb_ps_simhit_trid;   //!
   TBranch        *b_bb_ps_simhit_pid;   //!
   TBranch        *b_bb_ps_simhit_chan;   //!
   TBranch        *b_bb_ps_simhit_edep;   //!
   TBranch        *b_bb_ps_simhit_npe;   //!
   TBranch        *b_bb_ps_simhit_time;   //!
   TBranch        *b_bb_ps_hit_nhits;   //!
   TBranch        *b_bb_ps_hit_chan;   //!
   TBranch        *b_bb_ps_hit_dataword;   //!
   TBranch        *b_bb_ps_hit_adc;   //!
   TBranch        *b_bb_hodo_simhit_nhits;   //!
   TBranch        *b_bb_hodo_simhit_src;   //!
   TBranch        *b_bb_hodo_simhit_trid;   //!
   TBranch        *b_bb_hodo_simhit_pid;   //!
   TBranch        *b_bb_hodo_simhit_chan;   //!
   TBranch        *b_bb_hodo_simhit_edep;   //!
   TBranch        *b_bb_hodo_simhit_npe;   //!
   TBranch        *b_bb_hodo_simhit_time;   //!
   TBranch        *b_bb_hodo_simhit_t_lead;   //!
   TBranch        *b_bb_hodo_simhit_t_trail;   //!
   TBranch        *b_bb_hodo_hit_nhits;   //!
   TBranch        *b_bb_hodo_hit_chan;   //!
   TBranch        *b_bb_hodo_hit_dataword;   //!
   TBranch        *b_bb_hodo_hit_tdc_l;   //!
   TBranch        *b_bb_hodo_hit_tdc_t;   //!
   TBranch        *b_bb_grinch_simhit_nhits;   //!
   TBranch        *b_bb_grinch_simhit_src;   //!
   TBranch        *b_bb_grinch_simhit_trid;   //!
   TBranch        *b_bb_grinch_simhit_pid;   //!
   TBranch        *b_bb_grinch_simhit_chan;   //!
   TBranch        *b_bb_grinch_simhit_npe;   //!
   TBranch        *b_bb_grinch_simhit_time;   //!
   TBranch        *b_bb_grinch_simhit_t_lead;   //!
   TBranch        *b_bb_grinch_simhit_t_trail;   //!
   TBranch        *b_bb_grinch_hit_nhits;   //!
   TBranch        *b_bb_grinch_hit_chan;   //!
   TBranch        *b_bb_grinch_hit_dataword;   //!
   TBranch        *b_bb_grinch_hit_tdc_l;   //!
   TBranch        *b_bb_grinch_hit_tdc_t;   //!
   TBranch        *b_bb_gem_simhit_nhits;   //!
   TBranch        *b_bb_gem_simhit_src;   //!
   TBranch        *b_bb_gem_simhit_trid;   //!
   TBranch        *b_bb_gem_simhit_pid;   //!
   TBranch        *b_bb_gem_simhit_plane;   //!
   TBranch        *b_bb_gem_simhit_module;   //!
   TBranch        *b_bb_gem_simhit_edep;   //!
   TBranch        *b_bb_gem_simhit_time;   //!
   TBranch        *b_bb_gem_simhit_xpos;   //!
   TBranch        *b_bb_gem_simhit_ypos;   //!
   TBranch        *b_bb_gem_simhit_px;   //!
   TBranch        *b_bb_gem_simhit_py;   //!
   TBranch        *b_bb_gem_simhit_pz;   //!
   TBranch        *b_bb_gem_simhit_sizex;   //!
   TBranch        *b_bb_gem_simhit_sizey;   //!
   TBranch        *b_bb_gem_simhit_startx;   //!
   TBranch        *b_bb_gem_simhit_starty;   //!
   TBranch        *b_bb_gem_hit_nhits;   //!
   TBranch        *b_bb_gem_hit_plane;   //!
   TBranch        *b_bb_gem_hit_module;   //!
   TBranch        *b_bb_gem_hit_proj;   //!
   TBranch        *b_bb_gem_hit_nwords;   //!
   TBranch        *b_bb_gem_hit_strip;   //!
   TBranch        *b_bb_gem_hit_samp;   //!
   TBranch        *b_bb_gem_hit_samps_adc;   //!
   */
