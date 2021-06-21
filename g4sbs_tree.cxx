//#define g4sbs_tree_cxx
#include "g4sbs_tree.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>

//to read the detector list
#include "THaGlobals.h"
#include "THaApparatus.h"
#include "SBSBBShower.h"
#include "SBSBBTotalShower.h"

// g4sbs_tree constructor: the tree will be the 
// the boolean is a flag to consider(true) or ignore(false) the ECal_box and HCal_box data
//g4sbs_tree::g4sbs_tree(TTree *tree, Exp_t expt, bool pythia)
//, bool ecalbox, bool have_hcalbox) 
g4sbs_tree::g4sbs_tree(TTree *tree)//, std::vector<TString> det_list)
  : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
  if (tree == 0) {
    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/volatile/halla/sbs/efuchey/gmn13.5_elastic_20200228_17/elastic_0.root");
    if (!f || !f->IsOpen()) {
      f = new TFile("/volatile/halla/sbs/efuchey/gmn13.5_elastic_20200228_17/elastic_0.root");
    }
    f->GetObject("T",tree);
  }
  // fExpt = expt;
  // fPythia = pythia;
  // fEcalBox = ecalbox;
  // fHcalBox = have_hcalbox;
  Init(tree);
  //Init(tree, det_list);
}

//default destructor
g4sbs_tree::~g4sbs_tree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

//overload of the TTree::GetEntries() function
Int_t g4sbs_tree::GetEntries()
{
  // Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntries();
}
//overload of the TTree::GetEntry(Long64_t) function
Int_t g4sbs_tree::GetEntry(Long64_t entry)
{
  // Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}
Long64_t g4sbs_tree::LoadTree(Long64_t entry)
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

void g4sbs_tree::Init(TTree *tree)//, std::vector<TString> det_list)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).
  cout << "initialize g4sbs_tree" << endl;

  std::vector<TString> det_list;
  det_list.clear();
  
  TIter aiter(gHaApps);
  THaApparatus* app = 0;
  while( (app=(THaApparatus*)aiter()) ){
    TList* listdet = app->GetDetectors();
    TIter diter(listdet);
    TObject* det = 0;
    while( (det=(TObject*)diter()) ){
      det_list.push_back(det->GetName() );
      if(strcmp(app->GetDetector(det->GetName())->GetClassName(),"SBSBBTotalShower")==0){
	SBSBBTotalShower* TS = (SBSBBTotalShower*)app->GetDetector(det->GetName());
	det_list.push_back(TS->GetShower()->GetName());
	det_list.push_back(TS->GetPreShower()->GetName());
	
      }
    }
  }

   // Set object pointer

   Primaries_PID = 0;
   Primaries_genflag = 0;
   Primaries_Px = 0;
   Primaries_Py = 0;
   Primaries_Pz = 0;
   Primaries_vx = 0;
   Primaries_vy = 0;
   Primaries_vz = 0;
   Primaries_M = 0;
   Primaries_E = 0;
   Primaries_P = 0;
   Primaries_t = 0;
   Primaries_theta = 0;
   Primaries_phi = 0;
    
   
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   // (jc2): Why do we want to make a class again??
   // I disabled this so that we can force the tree to check each
   // SetBranchStatus for matches.
   //fChain->SetMakeClass(1);
   
   //fChain->Print();
   
   //Setup "Event branch": can be useful? But is it???
   //fChain->SetBranchAddress("ev", &ev_count, &b_ev);
   
   for(int k = 0; k<det_list.size(); k++){
     //GMN/GEN
     if(det_list[k]=="ps"){
       //printf(" ps  branches set up! \n");
       SetupDetBranch(Earm_BBPSTF1, "Earm.BBPSTF1.hit");
       SetupDetBranch(Earm_BBPS_Dig, "Earm.BBPS.dighit");
     }
     if(det_list[k]=="sh"){
       //printf(" sh  branches set up! \n");
       SetupDetBranch(Earm_BBSHTF1, "Earm.BBSHTF1.hit");
       SetupDetBranch(Earm_BBSH_Dig, "Earm.BBSH.dighit");
     }
     if(det_list[k]=="grinch"){
       //printf(" grinch  branches set up! \n");
       SetupDetBranch(Earm_GRINCH, "Earm.GRINCH.hit");
       SetupDetBranch(Earm_GRINCH_Dig, "Earm.GRINCH.dighit");
     }
     if(det_list[k]=="hodo"){
       //printf(" hodo  branches set up! \n");
       SetupDetBranch(Earm_BBHodoScint, "Earm.BBHodoScint.hit");
       SetupDetBranch(Earm_BBHodo_Dig, "Earm.BBHodo.dighit");
     }
     if(det_list[k]=="gem"){
       printf(" bbgem branches set up! \n");
       SetupDetBranch(Earm_BBGEM, "Earm.BBGEM.hit");
       SetupDetBranch(Earm_BBGEM_Track, "Earm.BBGEM.Track");
       SetupDetBranch(Earm_BBGEM_Dig, "Earm.BBGEM.dighit");
     }
     if(det_list[k]=="hcal"){
       //printf(" hcal branches set up! \n");
       SetupDetBranch(Harm_HCalScint,"Harm.HCalScint.hit");
       SetupDetBranch(Harm_HCal_Dig, "Harm.HCal.dighit");
     }
     //GENRP
     if(det_list[k]=="h_cdet"){
       SetupDetBranch(CDET_Scint,"Harm.CDET_Scint.hit");
       SetupDetBranch(CDET_Dig, "Harm.CDET.dighit");
     }
     if(det_list[k]=="activeana"){
       SetupDetBranch(Harm_ActAnScint, "Harm.ActAnScint.hit");
       SetupDetBranch(Harm_ActAn_Dig, "Harm.ActAn.dighit");
     }
     if(det_list[k]=="prpolscint_bs"){
       SetupDetBranch(Harm_PRPolScintBeamSide, "Harm.PRPolScintBeamSide.hit");
       SetupDetBranch(Harm_PRPolScintBeamSide_Dig, "Harm.PRPolScintBeamSide.dighit");
     }
     if(det_list[k]=="prpolscint_fs"){
       SetupDetBranch(Harm_PRPolScintFarSide, "Harm.PRPolScintFarSide.hit");
       SetupDetBranch(Harm_PRPolScintFarSide_Dig, "Harm.PRPolScintFarSide.dighit");
     }
     if(det_list[k]=="cepol_front"){
       SetupDetBranch(Harm_CEPolFront, "Harm.CEPolFront.hit");
       SetupDetBranch(Harm_CEPolFront_Dig, "Harm.CEPolFront.dighit");
       SetupDetBranch(Harm_CEPolFront_Track, "Harm.CEPolFront.Track");
     }
     if(det_list[k]=="cepol_rear"){
       SetupDetBranch(Harm_CEPolRear, "Harm.CEPolRear.hit");
       SetupDetBranch(Harm_CEPolRear_Dig, "Harm.CEPolRear.dighit");
       SetupDetBranch(Harm_CEPolRear_Track, "Harm.CEPolFront.Track");
     }
     if(det_list[k]=="prpolgem_bs"){
       SetupDetBranch(Harm_PrPolGEMBeamSide, "Harm.PRPolGEMBeamSide.hit");
       SetupDetBranch(Harm_PrPolGEMBeamSide_Dig, "Harm.PRPolGEMBeamSide.dighit");
       SetupDetBranch(Harm_PrPolGEMBeamSide_Track, "Harm.PRPolGEMBeamSide.Track");
     }
     if(det_list[k]=="prpolgem_fs"){
       SetupDetBranch(Harm_PrPolGEMFarSide, "Harm.PRPolGEMFarSide.hit");
       SetupDetBranch(Harm_PrPolGEMFarSide_Dig, "Harm.PRPolGEMFarSide.dighit");
       SetupDetBranch(Harm_PrPolGEMFarSide_Track, "Harm.PRPolGEMFarSide.Track");
     }
     //GEP
     if(det_list[k]=="e_cdet"){
       SetupDetBranch(CDET_Scint,"Earm.CDET_Scint.hit");
       SetupDetBranch(CDET_Dig,"Earm.CDET.dighit");
     }
     if(det_list[k]=="ecal"){
       SetupDetBranch(Earm_ECalTF1, "Earm.ECalTF1.hit");
       SetupDetBranch(Earm_ECal_Dig,"Earm.ECal.dighit");
     }
     if(det_list[k]=="ft"){
       SetupDetBranch(Harm_FT, "Harm.FT.hit");
       SetupDetBranch(Harm_FT_Dig, "Harm.FT.dighit");
       SetupDetBranch(Harm_FT_Track, "Harm.FT.Track");
     }
     if(det_list[k]=="fpp1"){
       SetupDetBranch(Harm_FPP1, "Harm.FPP1.hit");
       SetupDetBranch(Harm_FPP1_Dig, "Harm.FPP1.dighit");
       SetupDetBranch(Harm_FT_Track, "Harm.FPP1.Track");
     }
     if(det_list[k]=="fpp2"){
       SetupDetBranch(Harm_FPP2, "Harm.FPP2.hit");
       SetupDetBranch(Harm_FPP2_Dig, "Harm.FPP2.dighit");
       SetupDetBranch(Harm_FPP2_Track, "Harm.FT.Track");
     }
     if(det_list[k]=="sbsgem"){
       SetupDetBranch(Harm_SBSGEM, "Harm.SBSGEM.hit");
       SetupDetBranch(Harm_SBSGEM_Dig, "Harm.SBSGEM.dighit");
       SetupDetBranch(Harm_SBSGEM_Track, "Harm.SBSGEM.Track");
     }
     if(det_list[k]=="rich"){
       SetupDetBranch(Harm_RICH,"Harm.RICH.hit");
       SetupDetBranch(Harm_RICH_Dig,"Harm.RICH.dighit");
     }
   }

   /*
   //BigBite detector package: all expts except GEp
   if(fExpt==kGMN || fExpt==kGEN || fExpt==kGEnRP || fExpt==kSIDISExp || fExpt==kA1n){
     //gem_branch gem(BBGEM_UNIQUE_DETID,"Earm.BBGEM.hit","Earm.BBGEM.Track");
     //GEMs.push_back(gem);
     // if(fEcalBox){
     //   SetupDetBranch(Earm_ECAL_box, "Earm.BBCal.hit.nhits");
     // }else{
     // SetupDetBranch(Earm_BBPS,    "Earm.BBPS.hit");
     // SetupDetBranch(Earm_BBSH,    "Earm.BBSH.hit");
     //}
     
   }
   
   if(fExpt==kGEp){
     //SetupDetBranch(Earm_CDET,"Earm.CDET.hit");
     
     // if(fEcalBox){
     //   SetupDetBranch(Earm_ECAL_box, "Earm.ECAL_box.hit.nhits");
     // }else{
     //   SetupDetBranch(Earm_ECAL, "Earm.ECAL.hit");
     // }
     

     // Focal plane polarimeters
     
     // gem_branch gem_ft(FT_UNIQUE_DETID,"Harm.FT.hit","Harm.FT.Track");
     // GEMs.push_back(gem_ft);
     // gem_branch gem_fpp1(FPP1_UNIQUE_DETID,"Harm.FPP1.hit","Harm.FPP1.Track");
     // GEMs.push_back(gem_fpp1);
     // gem_branch gem_fpp2(FPP2_UNIQUE_DETID,"Harm.FPP2.hit","Harm.FPP2.Track");
     // GEMs.push_back(gem_fpp2);
   }

   if(fExpt==kGEnRP){
     //SetupDetBranch(Harm_CDET,"Harm.CDET.hit");
     //TODO: Add polarimeter GEMs
     //gem_branch gem_cefront(CEPOL_GEMFRONT_UNIQUE_DETID,"Harm.CEPolFront.hit","Harm.CEPolFront.Track");
     //GEMs.push_back(gem_cefront);
     //gem_branch gem_cerear(CEPOL_GEMREAR_UNIQUE_DETID,"Harm.CEPolRear.hit","Harm.CEPolRear.Track");
     //GEMs.push_back(gem_cerear);
     
     //gem_branch gem_prbs(PRPOLBS_GEM_UNIQUE_DETID,"Harm.PRPolGEMBeamSide.hit","Harm.PRPolGEMBeamSide.Track");
     //GEMs.push_back(gem_prbs);
     //gem_branch gem_prfs(PRPOLFS_GEM_UNIQUE_DETID,"Harm.PRPolGEMFarSide.hit","Harm.PRPolGEMFarSide.Track");
     //GEMs.push_back(gem_prfs);
   }
   if(fExpt!=kTDIS && fExpt!=kNDVCS){
     // Example of simplified HCAL branch setup
     // if(fHcalBox){
     //   SetupDetBranch(hcalbox,"Harm.HCAL_box.hit");
     // }else{
     //   SetupDetBranch(hcal,"Harm.HCal.hit");
     // }
     // EPAF: for the time being, remove.
     // we might reestablish it to turn on some test mode.
     // SetupDetBranch(hcalpart,"Harm.HCal");
   }
   
   if(fExpt==kSIDISExp || fExpt==kA1n || fExpt==kTDIS || fExpt==kNDVCS){
     //gem_branch gem(SBSGEM_UNIQUE_DETID,"Harm.SBSGEM.hit","Harm.SBSGEM.Track");
     //GEMs.push_back(gem);
     
   }
   
   if(fPythia){
     fChain->SetBranchAddress("primaries.Sigma", &primaries_Sigma, &b_primaries_Sigma);
     fChain->SetBranchAddress("primaries.Ebeam", &primaries_Ebeam, &b_primaries_Ebeam);
     fChain->SetBranchAddress("primaries.Eprime", &primaries_Eprime, &b_primaries_Eprime);
     fChain->SetBranchAddress("primaries.theta_e", &primaries_theta_e, &b_primaries_theta_e);
     fChain->SetBranchAddress("primaries.phi_e", &primaries_phi_e, &b_primaries_phi_e);
     fChain->SetBranchAddress("primaries.px_e", &primaries_px_e, &b_primaries_px_e);
     fChain->SetBranchAddress("primaries.py_e", &primaries_py_e, &b_primaries_py_e);
     fChain->SetBranchAddress("primaries.pz_e", &primaries_pz_e, &b_primaries_pz_e);
     fChain->SetBranchAddress("primaries.vx_e", &primaries_vx_e, &b_primaries_vx_e);
     fChain->SetBranchAddress("primaries.vy_e", &primaries_vy_e, &b_primaries_vy_e);
     fChain->SetBranchAddress("primaries.vz_e", &primaries_vz_e, &b_primaries_vz_e);
     fChain->SetBranchAddress("primaries.Egamma", &primaries_Egamma, &b_primaries_Egamma);
     fChain->SetBranchAddress("primaries.theta_gamma", &primaries_theta_gamma, &b_primaries_theta_gamma);
     fChain->SetBranchAddress("primaries.phi_gamma", &primaries_phi_gamma, &b_primaries_phi_gamma);
     fChain->SetBranchAddress("primaries.px_gamma", &primaries_px_gamma, &b_primaries_px_gamma);
     fChain->SetBranchAddress("primaries.py_gamma", &primaries_py_gamma, &b_primaries_py_gamma);
     fChain->SetBranchAddress("primaries.pz_gamma", &primaries_pz_gamma, &b_primaries_pz_gamma);
     fChain->SetBranchAddress("primaries.vx_gamma", &primaries_vx_gamma, &b_primaries_vx_gamma);
     fChain->SetBranchAddress("primaries.vy_gamma", &primaries_vy_gamma, &b_primaries_vy_gamma);
     fChain->SetBranchAddress("primaries.vz_gamma", &primaries_vz_gamma, &b_primaries_vz_gamma);
     
     fChain->SetBranchAddress("Primaries.Nprimaries", &Primaries_Nprimaries, &b_Primaries_Nprimaries);
     fChain->SetBranchAddress("Primaries.PID", &Primaries_PID, &b_Primaries_PID);
     fChain->SetBranchAddress("Primaries.genflag", &Primaries_genflag, &b_Primaries_genflag);
     fChain->SetBranchAddress("Primaries.Px", &Primaries_Px, &b_Primaries_Px);
     fChain->SetBranchAddress("Primaries.Py", &Primaries_Py, &b_Primaries_Py);
     fChain->SetBranchAddress("Primaries.Pz", &Primaries_Pz, &b_Primaries_Pz);
     fChain->SetBranchAddress("Primaries.vx", &Primaries_vx, &b_Primaries_vx);
     fChain->SetBranchAddress("Primaries.vy", &Primaries_vy, &b_Primaries_vy);
     fChain->SetBranchAddress("Primaries.vz", &Primaries_vz, &b_Primaries_vz);
     fChain->SetBranchAddress("Primaries.M", &Primaries_M, &b_Primaries_M);
     fChain->SetBranchAddress("Primaries.E", &Primaries_E, &b_Primaries_E);
     fChain->SetBranchAddress("Primaries.P", &Primaries_P, &b_Primaries_P);
     fChain->SetBranchAddress("Primaries.t", &Primaries_t, &b_Primaries_t);
     fChain->SetBranchAddress("Primaries.theta", &Primaries_theta, &b_Primaries_theta);
     fChain->SetBranchAddress("Primaries.phi", &Primaries_phi, &b_Primaries_phi);
   }

   // Now config all the GEM branches and trees

   for(std::vector<gem_branch>::iterator it = GEMs.begin();
       it != GEMs.end(); it++) {
     SetupDetBranch(it->tree,it->name);
     SetupDetBranch(it->Track_tree,it->Track_name);
   }
   */
   Notify();
}


Bool_t g4sbs_tree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void g4sbs_tree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}

Int_t g4sbs_tree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

void g4sbs_tree::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L gep_tree_with_spin.C
//      Root > gep_tree_with_spin t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
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


void g4sbs_tree::SetupDetBranch(TSBSGeant4::VDetData_t &det, const char *prefix)
{
  det.SetupBranches(fChain,prefix);
}

/*
void g4sbs_tree::ClearDigBranches()
{
  Earm_BBGEM_Dig.ClearBranches();
  Earm_BBHodo_Dig.ClearBranches();
  Earm_GRINCH_Dig.ClearBranches();
  Earm_BBPS_Dig.ClearBranches();
  Earm_BBSH_Dig.ClearBranches();
  CDET_Dig.ClearBranches();
  Earm_ECal_Dig.ClearBranches();
  Harm_FPP1_Dig.ClearBranches();
  Harm_FPP2_Dig.ClearBranches();
  Harm_FT_Dig.ClearBranches();
  Harm_HCal_Dig.ClearBranches();
  Harm_ActAn_Dig.ClearBranches();
  Harm_PRPolScintBeamSide_Dig.ClearBranches();
  Harm_PRPolScintFarSide_Dig.ClearBranches();
  Harm_CEPolFront_Dig.ClearBranches();
  Harm_CEPolRear_Dig.ClearBranches();
  Harm_PrPolGEMBeamSide_Dig.ClearBranches();
  Harm_PrPolGEMFarSide_Dig.ClearBranches();
  Harm_SBSGEM_Dig.ClearBranches();
  Harm_RICH_Dig.ClearBranches();
}

void g4sbs_tree::FillDigBranches()
{
  Earm_BBGEM_Dig.FillBranches();
  Earm_BBHodo_Dig.FillBranches();
  Earm_GRINCH_Dig.FillBranches();
  Earm_BBPS_Dig.FillBranches();
  Earm_BBSH_Dig.FillBranches();
  CDET_Dig.FillBranches();
  Earm_ECal_Dig.FillBranches();
  Harm_FPP1_Dig.FillBranches();
  Harm_FPP2_Dig.FillBranches();
  Harm_FT_Dig.FillBranches();
  Harm_HCal_Dig.FillBranches();
  Harm_ActAn_Dig.FillBranches();
  Harm_PRPolScintBeamSide_Dig.FillBranches();
  Harm_PRPolScintFarSide_Dig.FillBranches();
  Harm_CEPolFront_Dig.FillBranches();
  Harm_CEPolRear_Dig.FillBranches();
  Harm_PrPolGEMBeamSide_Dig.FillBranches();
  Harm_PrPolGEMFarSide_Dig.FillBranches();
  Harm_SBSGEM_Dig.FillBranches();
  Harm_RICH_Dig.FillBranches();
}
*/
