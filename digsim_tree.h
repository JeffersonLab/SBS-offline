//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Feb 12 21:41:18 2020 by ROOT version 6.14/04
// from TTree digtree/
// found on file: simdig_test.root
//////////////////////////////////////////////////////////

#ifndef data_digtree_h
#define data_digtree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"
#include "vector"
#include "vector"

class data_digtree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           RunID;
   Int_t           EvtID;
   Double_t        Weight;
   Int_t           NSignal;
   UInt_t          sbs_hcal_Nsimhits;
   vector<short>   *sbs_hcal_simhit_src;
   vector<short>   *sbs_hcal_simhit_trid;
   vector<int>     *sbs_hcal_simhit_pid;
   vector<short>   *sbs_hcal_simhit_chan;
   vector<double>  *sbs_hcal_simhit_Edep;
   vector<int>     *sbs_hcal_simhit_npe;
   vector<double>  *sbs_hcal_simhit_time;
   vector<double>  *sbs_hcal_simhit_t_lead;
   vector<double>  *sbs_hcal_simhit_t_trail;
   UInt_t          sbs_hcal_Nhits;
   vector<short>   *sbs_hcal_hit_chan;
   vector<unsigned int> *sbs_hcal_hit_nwords;
   vector<int>     *sbs_hcal_hit_adcsum;
   vector<vector<int> > *sbs_hcal_hit_samps_adc;
   vector<vector<unsigned int> > *sbs_hcal_hit_samps_datawords;
   vector<int>     *sbs_hcal_hit_tdc_l;
   vector<int>     *sbs_hcal_hit_tdc_t;
   UInt_t          sbs_cdet_Nsimhits;
   vector<short>   *sbs_cdet_simhit_src;
   vector<short>   *sbs_cdet_simhit_trid;
   vector<int>     *sbs_cdet_simhit_pid;
   vector<short>   *sbs_cdet_simhit_chan;
   vector<double>  *sbs_cdet_simhit_Edep;
   vector<int>     *sbs_cdet_simhit_npe;
   vector<double>  *sbs_cdet_simhit_time;
   vector<double>  *sbs_cdet_simhit_t_lead;
   vector<double>  *sbs_cdet_simhit_t_trail;
   UInt_t          sbs_cdet_Nhits;
   vector<short>   *sbs_cdet_hit_chan;
   vector<unsigned int> *sbs_cdet_hit_dataword;
   vector<int>     *sbs_cdet_hit_tdc_l;
   vector<int>     *sbs_cdet_hit_tdc_t;
   UInt_t          bb_sh_Nsimhits;
   vector<short>   *bb_sh_simhit_src;
   vector<short>   *bb_sh_simhit_trid;
   vector<int>     *bb_sh_simhit_pid;
   vector<short>   *bb_sh_simhit_chan;
   vector<double>  *bb_sh_simhit_Edep;
   vector<int>     *bb_sh_simhit_npe;
   vector<double>  *bb_sh_simhit_time;
   UInt_t          bb_sh_Nhits;
   vector<short>   *bb_sh_hit_chan;
   vector<unsigned int> *bb_sh_hit_dataword;
   vector<int>     *bb_sh_hit_adc;
   UInt_t          bb_ps_Nsimhits;
   vector<short>   *bb_ps_simhit_src;
   vector<short>   *bb_ps_simhit_trid;
   vector<int>     *bb_ps_simhit_pid;
   vector<short>   *bb_ps_simhit_chan;
   vector<double>  *bb_ps_simhit_Edep;
   vector<int>     *bb_ps_simhit_npe;
   vector<double>  *bb_ps_simhit_time;
   UInt_t          bb_ps_Nhits;
   vector<short>   *bb_ps_hit_chan;
   vector<unsigned int> *bb_ps_hit_dataword;
   vector<int>     *bb_ps_hit_adc;
   UInt_t          bb_hodo_Nsimhits;
   vector<short>   *bb_hodo_simhit_src;
   vector<short>   *bb_hodo_simhit_trid;
   vector<int>     *bb_hodo_simhit_pid;
   vector<short>   *bb_hodo_simhit_chan;
   vector<double>  *bb_hodo_simhit_Edep;
   vector<int>     *bb_hodo_simhit_npe;
   vector<double>  *bb_hodo_simhit_time;
   vector<double>  *bb_hodo_simhit_t_lead;
   vector<double>  *bb_hodo_simhit_t_trail;
   UInt_t          bb_hodo_Nhits;
   vector<short>   *bb_hodo_hit_chan;
   vector<unsigned int> *bb_hodo_hit_dataword;
   vector<int>     *bb_hodo_hit_tdc_l;
   vector<int>     *bb_hodo_hit_tdc_t;
   UInt_t          bb_grinch_Nsimhits;
   vector<short>   *bb_grinch_simhit_src;
   vector<short>   *bb_grinch_simhit_trid;
   vector<int>     *bb_grinch_simhit_pid;
   vector<short>   *bb_grinch_simhit_chan;
   vector<int>     *bb_grinch_simhit_npe;
   vector<double>  *bb_grinch_simhit_time;
   vector<double>  *bb_grinch_simhit_t_lead;
   vector<double>  *bb_grinch_simhit_t_trail;
   UInt_t          bb_grinch_Nhits;
   vector<short>   *bb_grinch_hit_chan;
   vector<unsigned int> *bb_grinch_hit_dataword;
   vector<int>     *bb_grinch_hit_tdc_l;
   vector<int>     *bb_grinch_hit_tdc_t;
   UInt_t          bb_gem_Nsimhits;
   vector<short>   *bb_gem_simhit_src;
   vector<short>   *bb_gem_simhit_trid;
   vector<int>     *bb_gem_simhit_pid;
   vector<short>   *bb_gem_simhit_plane;
   vector<short>   *bb_gem_simhit_module;
   vector<double>  *bb_gem_simhit_Edep;
   vector<double>  *bb_gem_simhit_time;
   vector<double>  *bb_gem_simhit_xpos;
   vector<double>  *bb_gem_simhit_ypos;
   vector<double>  *bb_gem_simhit_px;
   vector<double>  *bb_gem_simhit_py;
   vector<double>  *bb_gem_simhit_pz;
   vector<short>   *bb_gem_simhit_sizex;
   vector<short>   *bb_gem_simhit_sizey;
   vector<short>   *bb_gem_simhit_startx;
   vector<short>   *bb_gem_simhit_starty;
   UInt_t          bb_gem_NHits;
   vector<short>   *bb_gem_Plane;
   vector<short>   *bb_gem_Module;
   vector<short>   *bb_gem_Proj;
   vector<unsigned int> *bb_gem_nwords;
   vector<vector<short> > *bb_gem_Strip;
   vector<vector<short> > *bb_gem_Samp;
   vector<vector<int> > *bb_gem_hit_samps_adc;

   // List of branches
   TBranch        *b_RunID;   //!
   TBranch        *b_EvtID;   //!
   TBranch        *b_Weight;   //!
   TBranch        *b_NSignal;   //!
   TBranch        *b_sbs_hcal_Nsimhits;   //!
   TBranch        *b_sbs_hcal_simhit_src;   //!
   TBranch        *b_sbs_hcal_simhit_trid;   //!
   TBranch        *b_sbs_hcal_simhit_pid;   //!
   TBranch        *b_sbs_hcal_simhit_chan;   //!
   TBranch        *b_sbs_hcal_simhit_Edep;   //!
   TBranch        *b_sbs_hcal_simhit_npe;   //!
   TBranch        *b_sbs_hcal_simhit_time;   //!
   TBranch        *b_sbs_hcal_simhit_t_lead;   //!
   TBranch        *b_sbs_hcal_simhit_t_trail;   //!
   TBranch        *b_sbs_hcal_Nhits;   //!
   TBranch        *b_sbs_hcal_hit_chan;   //!
   TBranch        *b_sbs_hcal_hit_nwords;   //!
   TBranch        *b_sbs_hcal_hit_adcsum;   //!
   TBranch        *b_sbs_hcal_hit_samps_adc;   //!
   TBranch        *b_sbs_hcal_hit_samps_datawords;   //!
   TBranch        *b_sbs_hcal_hit_tdc_l;   //!
   TBranch        *b_sbs_hcal_hit_tdc_t;   //!
   TBranch        *b_sbs_cdet_Nsimhits;   //!
   TBranch        *b_sbs_cdet_simhit_src;   //!
   TBranch        *b_sbs_cdet_simhit_trid;   //!
   TBranch        *b_sbs_cdet_simhit_pid;   //!
   TBranch        *b_sbs_cdet_simhit_chan;   //!
   TBranch        *b_sbs_cdet_simhit_Edep;   //!
   TBranch        *b_sbs_cdet_simhit_npe;   //!
   TBranch        *b_sbs_cdet_simhit_time;   //!
   TBranch        *b_sbs_cdet_simhit_t_lead;   //!
   TBranch        *b_sbs_cdet_simhit_t_trail;   //!
   TBranch        *b_sbs_cdet_Nhits;   //!
   TBranch        *b_sbs_cdet_hit_chan;   //!
   TBranch        *b_sbs_cdet_hit_dataword;   //!
   TBranch        *b_sbs_cdet_hit_tdc_l;   //!
   TBranch        *b_sbs_cdet_hit_tdc_t;   //!
   TBranch        *b_bb_sh_Nsimhits;   //!
   TBranch        *b_bb_sh_simhit_src;   //!
   TBranch        *b_bb_sh_simhit_trid;   //!
   TBranch        *b_bb_sh_simhit_pid;   //!
   TBranch        *b_bb_sh_simhit_chan;   //!
   TBranch        *b_bb_sh_simhit_Edep;   //!
   TBranch        *b_bb_sh_simhit_npe;   //!
   TBranch        *b_bb_sh_simhit_time;   //!
   TBranch        *b_bb_sh_Nhits;   //!
   TBranch        *b_bb_sh_hit_chan;   //!
   TBranch        *b_bb_sh_hit_dataword;   //!
   TBranch        *b_bb_sh_hit_adc;   //!
   TBranch        *b_bb_ps_Nsimhits;   //!
   TBranch        *b_bb_ps_simhit_src;   //!
   TBranch        *b_bb_ps_simhit_trid;   //!
   TBranch        *b_bb_ps_simhit_pid;   //!
   TBranch        *b_bb_ps_simhit_chan;   //!
   TBranch        *b_bb_ps_simhit_Edep;   //!
   TBranch        *b_bb_ps_simhit_npe;   //!
   TBranch        *b_bb_ps_simhit_time;   //!
   TBranch        *b_bb_ps_Nhits;   //!
   TBranch        *b_bb_ps_hit_chan;   //!
   TBranch        *b_bb_ps_hit_dataword;   //!
   TBranch        *b_bb_ps_hit_adc;   //!
   TBranch        *b_bb_hodo_Nsimhits;   //!
   TBranch        *b_bb_hodo_simhit_src;   //!
   TBranch        *b_bb_hodo_simhit_trid;   //!
   TBranch        *b_bb_hodo_simhit_pid;   //!
   TBranch        *b_bb_hodo_simhit_chan;   //!
   TBranch        *b_bb_hodo_simhit_Edep;   //!
   TBranch        *b_bb_hodo_simhit_npe;   //!
   TBranch        *b_bb_hodo_simhit_time;   //!
   TBranch        *b_bb_hodo_simhit_t_lead;   //!
   TBranch        *b_bb_hodo_simhit_t_trail;   //!
   TBranch        *b_bb_hodo_Nhits;   //!
   TBranch        *b_bb_hodo_hit_chan;   //!
   TBranch        *b_bb_hodo_hit_dataword;   //!
   TBranch        *b_bb_hodo_hit_tdc_l;   //!
   TBranch        *b_bb_hodo_hit_tdc_t;   //!
   TBranch        *b_bb_grinch_Nsimhits;   //!
   TBranch        *b_bb_grinch_simhit_src;   //!
   TBranch        *b_bb_grinch_simhit_trid;   //!
   TBranch        *b_bb_grinch_simhit_pid;   //!
   TBranch        *b_bb_grinch_simhit_chan;   //!
   TBranch        *b_bb_grinch_simhit_npe;   //!
   TBranch        *b_bb_grinch_simhit_time;   //!
   TBranch        *b_bb_grinch_simhit_t_lead;   //!
   TBranch        *b_bb_grinch_simhit_t_trail;   //!
   TBranch        *b_bb_grinch_Nhits;   //!
   TBranch        *b_bb_grinch_hit_chan;   //!
   TBranch        *b_bb_grinch_hit_dataword;   //!
   TBranch        *b_bb_grinch_hit_tdc_l;   //!
   TBranch        *b_bb_grinch_hit_tdc_t;   //!
   TBranch        *b_bb_gem_Nsimhits;   //!
   TBranch        *b_bb_gem_simhit_src;   //!
   TBranch        *b_bb_gem_simhit_trid;   //!
   TBranch        *b_bb_gem_simhit_pid;   //!
   TBranch        *b_bb_gem_simhit_plane;   //!
   TBranch        *b_bb_gem_simhit_module;   //!
   TBranch        *b_bb_gem_simhit_Edep;   //!
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
   TBranch        *b_bb_gem_NHits;   //!
   TBranch        *b_bb_gem_Plane;   //!
   TBranch        *b_bb_gem_Module;   //!
   TBranch        *b_bb_gem_Proj;   //!
   TBranch        *b_bb_gem_nwords;   //!
   TBranch        *b_bb_gem_Strip;   //!
   TBranch        *b_bb_gem_Samp;   //!
   TBranch        *b_bb_gem_hit_samps_adc;   //!

   data_digtree(TTree *tree=0);
   virtual ~data_digtree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef data_digtree_cxx
data_digtree::data_digtree(TTree *tree) : fChain(0) 
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

data_digtree::~data_digtree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t data_digtree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t data_digtree::LoadTree(Long64_t entry)
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

void data_digtree::Init(TTree *tree)
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
   sbs_hcal_simhit_Edep = 0;
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
   sbs_cdet_simhit_src = 0;
   sbs_cdet_simhit_trid = 0;
   sbs_cdet_simhit_pid = 0;
   sbs_cdet_simhit_chan = 0;
   sbs_cdet_simhit_Edep = 0;
   sbs_cdet_simhit_npe = 0;
   sbs_cdet_simhit_time = 0;
   sbs_cdet_simhit_t_lead = 0;
   sbs_cdet_simhit_t_trail = 0;
   sbs_cdet_hit_chan = 0;
   sbs_cdet_hit_dataword = 0;
   sbs_cdet_hit_tdc_l = 0;
   sbs_cdet_hit_tdc_t = 0;
   bb_sh_simhit_src = 0;
   bb_sh_simhit_trid = 0;
   bb_sh_simhit_pid = 0;
   bb_sh_simhit_chan = 0;
   bb_sh_simhit_Edep = 0;
   bb_sh_simhit_npe = 0;
   bb_sh_simhit_time = 0;
   bb_sh_hit_chan = 0;
   bb_sh_hit_dataword = 0;
   bb_sh_hit_adc = 0;
   bb_ps_simhit_src = 0;
   bb_ps_simhit_trid = 0;
   bb_ps_simhit_pid = 0;
   bb_ps_simhit_chan = 0;
   bb_ps_simhit_Edep = 0;
   bb_ps_simhit_npe = 0;
   bb_ps_simhit_time = 0;
   bb_ps_hit_chan = 0;
   bb_ps_hit_dataword = 0;
   bb_ps_hit_adc = 0;
   bb_hodo_simhit_src = 0;
   bb_hodo_simhit_trid = 0;
   bb_hodo_simhit_pid = 0;
   bb_hodo_simhit_chan = 0;
   bb_hodo_simhit_Edep = 0;
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
   bb_gem_simhit_Edep = 0;
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
   bb_gem_Plane = 0;
   bb_gem_Module = 0;
   bb_gem_Proj = 0;
   bb_gem_nwords = 0;
   bb_gem_Strip = 0;
   bb_gem_Samp = 0;
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
   fChain->SetBranchAddress("sbs.hcal_Nsimhits", &sbs_hcal_Nsimhits, &b_sbs_hcal_Nsimhits);
   fChain->SetBranchAddress("sbs.hcal_simhit_src", &sbs_hcal_simhit_src, &b_sbs_hcal_simhit_src);
   fChain->SetBranchAddress("sbs.hcal_simhit_trid", &sbs_hcal_simhit_trid, &b_sbs_hcal_simhit_trid);
   fChain->SetBranchAddress("sbs.hcal_simhit_pid", &sbs_hcal_simhit_pid, &b_sbs_hcal_simhit_pid);
   fChain->SetBranchAddress("sbs.hcal_simhit_chan", &sbs_hcal_simhit_chan, &b_sbs_hcal_simhit_chan);
   fChain->SetBranchAddress("sbs.hcal_simhit_Edep", &sbs_hcal_simhit_Edep, &b_sbs_hcal_simhit_Edep);
   fChain->SetBranchAddress("sbs.hcal_simhit_npe", &sbs_hcal_simhit_npe, &b_sbs_hcal_simhit_npe);
   fChain->SetBranchAddress("sbs.hcal_simhit_time", &sbs_hcal_simhit_time, &b_sbs_hcal_simhit_time);
   fChain->SetBranchAddress("sbs.hcal_simhit_t_lead", &sbs_hcal_simhit_t_lead, &b_sbs_hcal_simhit_t_lead);
   fChain->SetBranchAddress("sbs.hcal_simhit_t_trail", &sbs_hcal_simhit_t_trail, &b_sbs_hcal_simhit_t_trail);
   fChain->SetBranchAddress("sbs.hcal_Nhits", &sbs_hcal_Nhits, &b_sbs_hcal_Nhits);
   fChain->SetBranchAddress("sbs.hcal_hit_chan", &sbs_hcal_hit_chan, &b_sbs_hcal_hit_chan);
   fChain->SetBranchAddress("sbs.hcal_hit_nwords", &sbs_hcal_hit_nwords, &b_sbs_hcal_hit_nwords);
   fChain->SetBranchAddress("sbs.hcal_hit_adcsum", &sbs_hcal_hit_adcsum, &b_sbs_hcal_hit_adcsum);
   fChain->SetBranchAddress("sbs.hcal_hit_samps_adc", &sbs_hcal_hit_samps_adc, &b_sbs_hcal_hit_samps_adc);
   fChain->SetBranchAddress("sbs.hcal_hit_samps_datawords", &sbs_hcal_hit_samps_datawords, &b_sbs_hcal_hit_samps_datawords);
   fChain->SetBranchAddress("sbs.hcal_hit_tdc_l", &sbs_hcal_hit_tdc_l, &b_sbs_hcal_hit_tdc_l);
   fChain->SetBranchAddress("sbs.hcal_hit_tdc_t", &sbs_hcal_hit_tdc_t, &b_sbs_hcal_hit_tdc_t);
   fChain->SetBranchAddress("sbs.cdet_Nsimhits", &sbs_cdet_Nsimhits, &b_sbs_cdet_Nsimhits);
   fChain->SetBranchAddress("sbs.cdet_simhit_src", &sbs_cdet_simhit_src, &b_sbs_cdet_simhit_src);
   fChain->SetBranchAddress("sbs.cdet_simhit_trid", &sbs_cdet_simhit_trid, &b_sbs_cdet_simhit_trid);
   fChain->SetBranchAddress("sbs.cdet_simhit_pid", &sbs_cdet_simhit_pid, &b_sbs_cdet_simhit_pid);
   fChain->SetBranchAddress("sbs.cdet_simhit_chan", &sbs_cdet_simhit_chan, &b_sbs_cdet_simhit_chan);
   fChain->SetBranchAddress("sbs.cdet_simhit_Edep", &sbs_cdet_simhit_Edep, &b_sbs_cdet_simhit_Edep);
   fChain->SetBranchAddress("sbs.cdet_simhit_npe", &sbs_cdet_simhit_npe, &b_sbs_cdet_simhit_npe);
   fChain->SetBranchAddress("sbs.cdet_simhit_time", &sbs_cdet_simhit_time, &b_sbs_cdet_simhit_time);
   fChain->SetBranchAddress("sbs.cdet_simhit_t_lead", &sbs_cdet_simhit_t_lead, &b_sbs_cdet_simhit_t_lead);
   fChain->SetBranchAddress("sbs.cdet_simhit_t_trail", &sbs_cdet_simhit_t_trail, &b_sbs_cdet_simhit_t_trail);
   fChain->SetBranchAddress("sbs.cdet_Nhits", &sbs_cdet_Nhits, &b_sbs_cdet_Nhits);
   fChain->SetBranchAddress("sbs.cdet_hit_chan", &sbs_cdet_hit_chan, &b_sbs_cdet_hit_chan);
   fChain->SetBranchAddress("sbs.cdet_hit_dataword", &sbs_cdet_hit_dataword, &b_sbs_cdet_hit_dataword);
   fChain->SetBranchAddress("sbs.cdet_hit_tdc_l", &sbs_cdet_hit_tdc_l, &b_sbs_cdet_hit_tdc_l);
   fChain->SetBranchAddress("sbs.cdet_hit_tdc_t", &sbs_cdet_hit_tdc_t, &b_sbs_cdet_hit_tdc_t);
   fChain->SetBranchAddress("bb.sh_Nsimhits", &bb_sh_Nsimhits, &b_bb_sh_Nsimhits);
   fChain->SetBranchAddress("bb.sh_simhit_src", &bb_sh_simhit_src, &b_bb_sh_simhit_src);
   fChain->SetBranchAddress("bb.sh_simhit_trid", &bb_sh_simhit_trid, &b_bb_sh_simhit_trid);
   fChain->SetBranchAddress("bb.sh_simhit_pid", &bb_sh_simhit_pid, &b_bb_sh_simhit_pid);
   fChain->SetBranchAddress("bb.sh_simhit_chan", &bb_sh_simhit_chan, &b_bb_sh_simhit_chan);
   fChain->SetBranchAddress("bb.sh_simhit_Edep", &bb_sh_simhit_Edep, &b_bb_sh_simhit_Edep);
   fChain->SetBranchAddress("bb.sh_simhit_npe", &bb_sh_simhit_npe, &b_bb_sh_simhit_npe);
   fChain->SetBranchAddress("bb.sh_simhit_time", &bb_sh_simhit_time, &b_bb_sh_simhit_time);
   fChain->SetBranchAddress("bb.sh_Nhits", &bb_sh_Nhits, &b_bb_sh_Nhits);
   fChain->SetBranchAddress("bb.sh_hit_chan", &bb_sh_hit_chan, &b_bb_sh_hit_chan);
   fChain->SetBranchAddress("bb.sh_hit_dataword", &bb_sh_hit_dataword, &b_bb_sh_hit_dataword);
   fChain->SetBranchAddress("bb.sh_hit_adc", &bb_sh_hit_adc, &b_bb_sh_hit_adc);
   fChain->SetBranchAddress("bb.ps_Nsimhits", &bb_ps_Nsimhits, &b_bb_ps_Nsimhits);
   fChain->SetBranchAddress("bb.ps_simhit_src", &bb_ps_simhit_src, &b_bb_ps_simhit_src);
   fChain->SetBranchAddress("bb.ps_simhit_trid", &bb_ps_simhit_trid, &b_bb_ps_simhit_trid);
   fChain->SetBranchAddress("bb.ps_simhit_pid", &bb_ps_simhit_pid, &b_bb_ps_simhit_pid);
   fChain->SetBranchAddress("bb.ps_simhit_chan", &bb_ps_simhit_chan, &b_bb_ps_simhit_chan);
   fChain->SetBranchAddress("bb.ps_simhit_Edep", &bb_ps_simhit_Edep, &b_bb_ps_simhit_Edep);
   fChain->SetBranchAddress("bb.ps_simhit_npe", &bb_ps_simhit_npe, &b_bb_ps_simhit_npe);
   fChain->SetBranchAddress("bb.ps_simhit_time", &bb_ps_simhit_time, &b_bb_ps_simhit_time);
   fChain->SetBranchAddress("bb.ps_Nhits", &bb_ps_Nhits, &b_bb_ps_Nhits);
   fChain->SetBranchAddress("bb.ps_hit_chan", &bb_ps_hit_chan, &b_bb_ps_hit_chan);
   fChain->SetBranchAddress("bb.ps_hit_dataword", &bb_ps_hit_dataword, &b_bb_ps_hit_dataword);
   fChain->SetBranchAddress("bb.ps_hit_adc", &bb_ps_hit_adc, &b_bb_ps_hit_adc);
   fChain->SetBranchAddress("bb.hodo_Nsimhits", &bb_hodo_Nsimhits, &b_bb_hodo_Nsimhits);
   fChain->SetBranchAddress("bb.hodo_simhit_src", &bb_hodo_simhit_src, &b_bb_hodo_simhit_src);
   fChain->SetBranchAddress("bb.hodo_simhit_trid", &bb_hodo_simhit_trid, &b_bb_hodo_simhit_trid);
   fChain->SetBranchAddress("bb.hodo_simhit_pid", &bb_hodo_simhit_pid, &b_bb_hodo_simhit_pid);
   fChain->SetBranchAddress("bb.hodo_simhit_chan", &bb_hodo_simhit_chan, &b_bb_hodo_simhit_chan);
   fChain->SetBranchAddress("bb.hodo_simhit_Edep", &bb_hodo_simhit_Edep, &b_bb_hodo_simhit_Edep);
   fChain->SetBranchAddress("bb.hodo_simhit_npe", &bb_hodo_simhit_npe, &b_bb_hodo_simhit_npe);
   fChain->SetBranchAddress("bb.hodo_simhit_time", &bb_hodo_simhit_time, &b_bb_hodo_simhit_time);
   fChain->SetBranchAddress("bb.hodo_simhit_t_lead", &bb_hodo_simhit_t_lead, &b_bb_hodo_simhit_t_lead);
   fChain->SetBranchAddress("bb.hodo_simhit_t_trail", &bb_hodo_simhit_t_trail, &b_bb_hodo_simhit_t_trail);
   fChain->SetBranchAddress("bb.hodo_Nhits", &bb_hodo_Nhits, &b_bb_hodo_Nhits);
   fChain->SetBranchAddress("bb.hodo_hit_chan", &bb_hodo_hit_chan, &b_bb_hodo_hit_chan);
   fChain->SetBranchAddress("bb.hodo_hit_dataword", &bb_hodo_hit_dataword, &b_bb_hodo_hit_dataword);
   fChain->SetBranchAddress("bb.hodo_hit_tdc_l", &bb_hodo_hit_tdc_l, &b_bb_hodo_hit_tdc_l);
   fChain->SetBranchAddress("bb.hodo_hit_tdc_t", &bb_hodo_hit_tdc_t, &b_bb_hodo_hit_tdc_t);
   fChain->SetBranchAddress("bb.grinch_Nsimhits", &bb_grinch_Nsimhits, &b_bb_grinch_Nsimhits);
   fChain->SetBranchAddress("bb.grinch_simhit_src", &bb_grinch_simhit_src, &b_bb_grinch_simhit_src);
   fChain->SetBranchAddress("bb.grinch_simhit_trid", &bb_grinch_simhit_trid, &b_bb_grinch_simhit_trid);
   fChain->SetBranchAddress("bb.grinch_simhit_pid", &bb_grinch_simhit_pid, &b_bb_grinch_simhit_pid);
   fChain->SetBranchAddress("bb.grinch_simhit_chan", &bb_grinch_simhit_chan, &b_bb_grinch_simhit_chan);
   fChain->SetBranchAddress("bb.grinch_simhit_npe", &bb_grinch_simhit_npe, &b_bb_grinch_simhit_npe);
   fChain->SetBranchAddress("bb.grinch_simhit_time", &bb_grinch_simhit_time, &b_bb_grinch_simhit_time);
   fChain->SetBranchAddress("bb.grinch_simhit_t_lead", &bb_grinch_simhit_t_lead, &b_bb_grinch_simhit_t_lead);
   fChain->SetBranchAddress("bb.grinch_simhit_t_trail", &bb_grinch_simhit_t_trail, &b_bb_grinch_simhit_t_trail);
   fChain->SetBranchAddress("bb.grinch_Nhits", &bb_grinch_Nhits, &b_bb_grinch_Nhits);
   fChain->SetBranchAddress("bb.grinch_hit_chan", &bb_grinch_hit_chan, &b_bb_grinch_hit_chan);
   fChain->SetBranchAddress("bb.grinch_hit_dataword", &bb_grinch_hit_dataword, &b_bb_grinch_hit_dataword);
   fChain->SetBranchAddress("bb.grinch_hit_tdc_l", &bb_grinch_hit_tdc_l, &b_bb_grinch_hit_tdc_l);
   fChain->SetBranchAddress("bb.grinch_hit_tdc_t", &bb_grinch_hit_tdc_t, &b_bb_grinch_hit_tdc_t);
   fChain->SetBranchAddress("bb.gem_Nsimhits", &bb_gem_Nsimhits, &b_bb_gem_Nsimhits);
   fChain->SetBranchAddress("bb.gem_simhit_src", &bb_gem_simhit_src, &b_bb_gem_simhit_src);
   fChain->SetBranchAddress("bb.gem_simhit_trid", &bb_gem_simhit_trid, &b_bb_gem_simhit_trid);
   fChain->SetBranchAddress("bb.gem_simhit_pid", &bb_gem_simhit_pid, &b_bb_gem_simhit_pid);
   fChain->SetBranchAddress("bb.gem_simhit_plane", &bb_gem_simhit_plane, &b_bb_gem_simhit_plane);
   fChain->SetBranchAddress("bb.gem_simhit_module", &bb_gem_simhit_module, &b_bb_gem_simhit_module);
   fChain->SetBranchAddress("bb.gem_simhit_Edep", &bb_gem_simhit_Edep, &b_bb_gem_simhit_Edep);
   fChain->SetBranchAddress("bb.gem_simhit_time", &bb_gem_simhit_time, &b_bb_gem_simhit_time);
   fChain->SetBranchAddress("bb.gem_simhit_xpos", &bb_gem_simhit_xpos, &b_bb_gem_simhit_xpos);
   fChain->SetBranchAddress("bb.gem_simhit_ypos", &bb_gem_simhit_ypos, &b_bb_gem_simhit_ypos);
   fChain->SetBranchAddress("bb.gem_simhit_px", &bb_gem_simhit_px, &b_bb_gem_simhit_px);
   fChain->SetBranchAddress("bb.gem_simhit_py", &bb_gem_simhit_py, &b_bb_gem_simhit_py);
   fChain->SetBranchAddress("bb.gem_simhit_pz", &bb_gem_simhit_pz, &b_bb_gem_simhit_pz);
   fChain->SetBranchAddress("bb.gem_simhit_sizex", &bb_gem_simhit_sizex, &b_bb_gem_simhit_sizex);
   fChain->SetBranchAddress("bb.gem_simhit_sizey", &bb_gem_simhit_sizey, &b_bb_gem_simhit_sizey);
   fChain->SetBranchAddress("bb.gem_simhit_startx", &bb_gem_simhit_startx, &b_bb_gem_simhit_startx);
   fChain->SetBranchAddress("bb.gem_simhit_starty", &bb_gem_simhit_starty, &b_bb_gem_simhit_starty);
   fChain->SetBranchAddress("bb.gem_NHits", &bb_gem_NHits, &b_bb_gem_NHits);
   fChain->SetBranchAddress("bb.gem_Plane", &bb_gem_Plane, &b_bb_gem_Plane);
   fChain->SetBranchAddress("bb.gem_Module", &bb_gem_Module, &b_bb_gem_Module);
   fChain->SetBranchAddress("bb.gem_Proj", &bb_gem_Proj, &b_bb_gem_Proj);
   fChain->SetBranchAddress("bb.gem_nwords", &bb_gem_nwords, &b_bb_gem_nwords);
   fChain->SetBranchAddress("bb.gem_Strip", &bb_gem_Strip, &b_bb_gem_Strip);
   fChain->SetBranchAddress("bb.gem_Samp", &bb_gem_Samp, &b_bb_gem_Samp);
   fChain->SetBranchAddress("bb.gem_hit_samps_adc", &bb_gem_hit_samps_adc, &b_bb_gem_hit_samps_adc);
   Notify();
}

Bool_t data_digtree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void data_digtree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t data_digtree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef data_digtree_cxx
