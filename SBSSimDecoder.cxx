//*-- Author :    Ole Hansen (ole@jlab.org)    9-Dec-2011

/////////////////////////////////////////////////////////////////////
//
//   SBSSimDecoder
//
//   Decoder for SoLID simulation data
//
//   Interprets event buffer from input as SBSSimEvent objects
//   (containing digitized simulation data) and unpacks them into
//   crateslot arrays for low-level decoding by detectors.
//
/////////////////////////////////////////////////////////////////////

#include "SBSSimDecoder.h"
#include "SBSSimDataDecoder.h"
#include "THaCrateMap.h"
#include "THaBenchmark.h"
#include "VarDef.h"
#include "THaSlotData.h"

#include "TError.h"
#include "TSystem.h"
#include "TMath.h"
#include "TDatabasePDG.h"
#include "TRandom.h"
#include "THaVarList.h"
#include "THaDetMap.h"
#include "THaDetector.h"
#include "SBSBBShower.h"
#include "SBSBBTotalShower.h"
#include "THaCrateMap.h"
#include "Textvars.h"
//#include "THaAnalysisObject.h"

//#include <SBSSimFadc250Module.h>// we need not to need this
#include "TList.h"
#include "TObject.h"

#include <cstdlib>
#include <iostream>
#include <string>
#include <fstream>
#include <utility>
#include <stdexcept>

using namespace std;
using namespace Podd;

class THaAnalysisObject;

ClassImp(SBSSimDecoder) // Implements SBSSimDecoder


//static const Int_t kPrimaryType = 1, kPrimarySource = 0;
// Projection types must match the definitions in SBS-offline
//enum EProjType { kUPlane = 0, kVPlane =1, kXPlane = 2, kYPlane = 3};
//typedef vector<int>::size_type vsiz_t;

//-----------------------------------------------------------------------------
SBSSimDecoder::SBSSimDecoder()// : fCheckedForEnabledDetectors(false), fTreeIsSet(false)
{
  // Constructor
  DefineVariables();
  fDetectors.clear();
  //fTree = 0;
  // Load detectors: rely on gHaApps (please tell me it works!!!)
  //cout << " Calling SBSSimDecoder! "<< endl;
  //cout << " Make sure you have already declared your apparatuses and detectors, and added these to gHaApps" << endl;
  //SetDetectors();
  
  // h1_sizeHCal = new TH1D("h1_sizeHCal", "", 500, 0, 5000);
  // h1_sizeGEMs = new TH1D("h1_sizeGEMs", "", 500, 0, 5000);
  
  gSystem->Load("libEG.so");  // for TDatabasePDG
  // Get MPD encoder for GEMs
  // FIXME: a bit of a kludge... 
  // we shouldn't have to do that to initialize all encoders... shall we?
  fDecoderMPD = dynamic_cast<SBSSimSADCEncoder*>
    (SBSSimDataDecoder::GetEncoderByName("mpd"));
  
  fIsInit = false;
  
}

//-----------------------------------------------------------------------------
SBSSimDecoder::~SBSSimDecoder() {
  // h1_sizeHCal->Write();
  // h1_sizeGEMs->Write();
  //DefineVariables( THaAnalysisObject::kDelete );
  // h1_sizeHCal->Delete();
  // h1_sizeGEMs->Delete();
}

Int_t SBSSimDecoder::Init(){ 

  Int_t status = THaEvData::Init();

  fDetectors.clear();

  //status += DefineVariables();
  SetDetectors();

  fIsInit = true;

  return status;

}

//-----------------------------------------------------------------------------
Int_t SBSSimDecoder::DefineVariables( THaAnalysisObject::EMode mode )
{
  // Define global variables for the MC quantities. Extends the base
  // class method.
  
  const char* const here = "SBSSimDecoder::DefineVariables";
  
  if( mode == THaAnalysisObject::kDefine && fIsSetup )
    return THaAnalysisObject::kOK;
  
  SimDecoder::DefineVariables( mode );
  
  cout << "Read SBSSimDecoder variables " << endl;

  RVarDef vars[] = {
    //simc variables
    {"simc_sigma",    "MC cross section from SIMC gen.",   "fSigma_simc"},
    {"simc_Weight",   "MC cross section weight from SIMC gen.",   "fWeight_simc"},
    {"simc_Q2",       "MC Q2 from SIMC gen.",   "fQ2_simc"},
    {"simc_xbj",      "MC xbj from SIMC gen.",   "fXbj_simc"},
    {"simc_nu",       "MC nu from SIMC gen.",   "fNu_simc"},
    {"simc_W",        "MC W from SIMC gen.",   "fW_simc"},
    {"simc_epsilon",  "MC epsilon from SIMC gen.",   "fEpsilon_simc"},
    {"simc_Ebeam",    "MC Ebeam from SIMC gen.",   "fEbeam_simc"},
    {"simc_p_e",      "MC e momentum from SIMC gen.",   "fEp_simc"},
    {"simc_theta_e",  "MC e polar angle from SIMC gen.",   "fEtheta_simc"},
    {"simc_phi_e",    "MC e azimuthal angle from SIMC gen.",   "fEphi_simc"},
    {"simc_px_e",     "MC e mom. x componant from SIMC gen.",   "fEPx_simc"},
    {"simc_py_e",     "MC e mom. y componant from SIMC gen.",   "fEPy_simc"},
    {"simc_pz_e",     "MC e mom. z componant from SIMC gen.",   "fEPz_simc"},
    {"simc_fnucl",    "MC final-state nucleon type from SIMC gen.",   "fFnucl_simc"},
    {"simc_p_n",      "MC nucleon mom. from SIMC gen.",   "fNp_simc"},
    {"simc_theta_n",  "MC nucleon polar angle from SIMC gen.",   "fEtheta_simc"},
    {"simc_phi_n",    "MC nucleon azimuthal angle from SIMC gen.",   "fNphi_simc"},
    {"simc_px_n",     "MC nucleon mom. x componant from SIMC gen.",   "fNPx_simc"},
    {"simc_py_n",     "MC nucleon mom. y componant from SIMC gen.",   "fNPy_simc"},
    {"simc_pz_n",     "MC nucleon mom. z componant from SIMC gen.",   "fNPz_simc"},
    {"simc_vx",       "MC vertex x co-ordinate from SIMC gen.",   "fVx_simc"},
    {"simc_vy",       "MC vertex y co-ordinate from SIMC gen.",   "fVy_simc"},
    {"simc_vz",       "MC vertex z co-ordinate from SIMC gen.",   "fVz_simc"},
    {"simc_veE",      "MC scattered e- energy at vertex from SIMC gen.",   "fVeE_simc"},
    {"simc_vetheta",  "MC scattered e- theta at vertex from SIMC gen.",   "fVetheta_simc"},
    // ** ^^ **
    {"mc_sigma",   "MC cross section",   "fSigma"},
    {"mc_omega",   "MC phase spece generation",   "fOmega"},
    {"mc_epx",     "MC electron momentum x",   "fEPx"},
    {"mc_epy",     "MC electron momentum y",   "fEPy"},
    {"mc_epz",     "MC electron momentum z",   "fEPz"},
    {"mc_npx",     "MC nucleon momentum x",   "fNPx"},
    {"mc_npy",     "MC nucleon momentum y",   "fNPy"},
    {"mc_npz",     "MC nucleon momentum z",   "fNPz"},
    {"mc_vx",      "MC vertex x",   "fVx"},
    {"mc_vy",      "MC vertex y",   "fVy"},
    {"mc_vz",      "MC vertex z",   "fVz"},
    {"mc_ep",      "MC Initial momentum of the final state electron in GeV",   "fEp"},
    {"mc_np",      "MC Initial momentum of the final state nucleon in GeV",   "fNp"},
    {"mc_nucl",    "MC Initial (struck) nucleon type: 1 = proton, 0 = neutron",   "fNucl"},
    {"mc_fnucl",   "MC Final-state (detected) nucleon type: 1 = proton, 0 = neutron",   "fFnucl"},
    {"nbbtracks",   "number of BB MC tracks",   "fNBBtracks"},
    {"bbtrack_nhits",   "BB MC track hit mult",   "fBBtrack_Nhits"},
    {"bbtrack_tid",   "BB MC track TID",   "fBBtrack_TID"},
    {"bbtrack_pid",   "BB MC track PID",   "fBBtrack_PID"},
    {"bbtrack_mid",   "BB MC track MID",   "fBBtrack_MID"},
    {"bbtrack_p",   "BB MC track momentum",   "fBBtrack_P"},
    {"bbtrack_x",   "BB MC track transport X position",   "fBBtrack_X"},
    {"bbtrack_y",   "BB MC track transport Y position",   "fBBtrack_Y"},
    {"bbtrack_dx",   "BB MC track transport dX slope",   "fBBtrack_dX"},
    {"bbtrack_dy",   "BB MC track transport dY slope",   "fBBtrack_dY"},
    {"nbbgemhits",   "number of BBGEM MC hits",   "fNBBGEMhits"},
    {"bbgemhit_plane",   "BBGEM MC hit plane",   "fBBGEMhit_plane"},
    {"bbgemhit_tid",   "BBGEM MC hit TID",   "fBBGEMhit_TID"},
    {"bbgemhit_pid",   "BBGEM MC hit PID",   "fBBGEMhit_PID"},
    {"bbgemhit_mid",   "BBGEM MC hit MID",   "fBBGEMhit_MID"},
    {"bbgemhit_edep",   "BBGEM MC hit edep",   "fBBGEMhit_edep"},
    {"bbgemhit_x",   "BBGEM MC hit transport X",   "fBBGEMhit_x"},
    {"bbgemhit_y",   "BBGEM MC hit transport Y",   "fBBGEMhit_y"},
    {"bbps_esum",   "BBPS total energy sum",   "fBBPS_esum"},
    {"bbsh_esum",   "BBSH total energy sum",   "fBBSH_esum"},
    {"bbgemhit_ptridx",    "Primary track index for BBGEM SD",        "fBBGEMhit_ptridx"},
    {"bbgemhit_sdtridx",   "SD track index for BBGEM SD",             "fBBGEMhit_sdtridx"},
    {"bbgemtrack_ptridx",  "Primary track index for BBGEM Track SD",  "fBBGEMtrack_ptridx"},
    {"bbgemtrack_sdtridx", "SD track index for BBGEM Track SD",       "fBBGEMtrack_sdtridx"},
    {"bbhodohit_ptridx",   "Primary track index for BBHodo SD",       "fBBHODOhit_ptridx"},
    {"bbhodohit_sdtridx",  "SD track index for BBHodo SD",            "fBBHODOhit_sdtridx"},
    {"bbpshit_ptridx",   "Primary track index for BBPSTF1 SD",        "fBBPSTF1hit_ptridx"},
    {"bbpshit_sdtridx",  "SD track index for BBPSTF1 SD",             "fBBPSTF1hit_sdtridx"},
    {"bbshhit_ptridx",   "Primary track index for BBSHTF1 SD",        "fBBSHTF1hit_ptridx"},
    {"bbshhit_sdtridx",  "SD track index for BBSHTF1 SD",             "fBBSHTF1hit_sdtridx"},
    {"hcalhit_ptridx",   "Primary track index for HCalScint SD",      "fHCALhit_ptridx"},
    {"hcalhit_sdtridx",  "SD track index for HCalScint SD",           "fHCALhit_sdtridx"},
    {"ptrack_ntracks",   "Primary track ntracks", "fPTrack_ntacks"},
    {"ptrack_tid",       "Primary track TID",     "fPTrack_TID"},
    {"ptrack_pid",       "Primary track PID",     "fPTrack_PID"},
    {"ptrack_posx",      "Primary track posx",    "fPTrack_posx"},
    {"ptrack_posy",      "Primary track posy",    "fPTrack_posy"},
    {"ptrack_posz",      "Primary track posz",    "fPTrack_posz"},
    {"ptrack_momx",      "Primary track momx",    "fPTrack_momx"},
    {"ptrack_momy",      "Primary track momy",    "fPTrack_momy"},
    {"ptrack_momz",      "Primary track momz",    "fPTrack_momz"},
    {"ptrack_polx",      "Primary track polx",    "fPTrack_polx"},
    {"ptrack_poly",      "Primary track poly",    "fPTrack_poly"},
    {"ptrack_polz",      "Primary track polz",    "fPTrack_polz"},
    {"ptrack_etot",      "Primary track Etot",    "fPTrack_Etot"},
    {"ptrack_t",         "Primary track T",       "fPTrack_T"},
    {"sdtrack_ntracks",  "SD track ntracks",      "fSDTrack_ntacks"},
    {"sdtrack_tid",      "SD track TID",          "fSDTrack_TID"},
    {"sdtrack_mid",      "SD track MID",          "fSDTrack_MID"},
    {"sdtrack_pid",      "SD track PID",          "fSDTrack_PID"},
    {"sdtrack_posx",     "SD track posx",    "fSDTrack_posx"},
    {"sdtrack_posy",     "SD track posy",    "fSDTrack_posy"},
    {"sdtrack_posz",     "SD track posz",    "fSDTrack_posz"},
    {"sdtrack_momx",     "SD track momx",    "fSDTrack_momx"},
    {"sdtrack_momy",     "SD track momy",    "fSDTrack_momy"},
    {"sdtrack_momz",     "SD track momz",    "fSDTrack_momz"},
    {"sdtrack_polx",     "SD track polx",    "fSDTrack_polx"},
    {"sdtrack_poly",     "SD track poly",    "fSDTrack_poly"},
    {"sdtrack_polz",     "SD track polz",    "fSDTrack_polz"},
    {"sdtrack_etot",     "SD track Etot",    "fSDTrack_Etot"},
    {"sdtrack_t",        "SD track T",       "fSDTrack_T"},
    {"sdtrack_vx",       "SD track vx",    "fSDTrack_vx"},
    {"sdtrack_vy",       "SD track vy",    "fSDTrack_vy"},
    {"sdtrack_vz",       "SD track vz",    "fSDTrack_vz"},
    {"sdtrack_vnx",      "SD track vnx",    "fSDTrack_vnx"},
    {"sdtrack_vny",      "SD track vny",    "fSDTrack_vny"},
    {"sdtrack_vnz",      "SD track vnz",    "fSDTrack_vnz"},
    {"sdtrack_vEkin",    "SD track vEkin",  "fSDTrack_vEkin"},
    { 0 }
  };

  return THaAnalysisObject::
    DefineVarsFromList( vars, THaAnalysisObject::kRVarDef,
			mode, "", this, Podd::MC_PREFIX, here );
}

//-----------------------------------------------------------------------------
void SBSSimDecoder::Clear( Option_t* opt )
{
  // Clear track and plane data

  SimDecoder::Clear(opt);   // clears fMCCherHits, fMCCherClus
  
  //fPMTMap.clear(); 
}

//-----------------------------------------------------------------------------
#if ANALYZER_VERSION_CODE >= ANALYZER_VERSION(1,6,0)
int SBSSimDecoder::LoadEvent(const UInt_t* evbuffer )
#else
int SBSSimDecoder::LoadEvent(const Int_t* evbuffer )
#endif
{
  // Wrapper around DoLoadEvent so we can conveniently stop the benchmark
  // counter in case of errors

  if( !fIsInit ) Init();
  
  int ret = -1;
  if(sizeof(evbuffer)!=0){
    ret = DoLoadEvent( evbuffer );
  }
  
  if( fDoBench ) fBench->Stop("physics_decode");
  
  return ret;
}

//-----------------------------------------------------------------------------
#if ANALYZER_VERSION_CODE >= ANALYZER_VERSION(1,6,0)
Int_t SBSSimDecoder::DoLoadEvent(const UInt_t* evbuffer )
#else
Int_t SBSSimDecoder::DoLoadEvent(const Int_t* evbuffer )
#endif
{
  // Uncommenting the three lines below, 
  // SBS-offline processes 5000 simulated GMn events with no background within ~60s
  // instead of the ~85s it takes without.
  // commenting "event type = 1;", those 5000 events take ~3s
  /* 
  event_type = 1;
  event_num++;
  return HED_OK;
  */
  // Fill crateslot structures with Monte Carlo event data in 'evbuffer'
  static const char* const here = "SBSSimDecoder::LoadEvent";

#if ANALYZER_VERSION_CODE < ANALYZER_VERSION(1,6,0)
  Bool_t fNeedInit = fgNeedInit;
#endif
  assert( fMap || fNeedInit );
  
  // Local copy of evbuffer pointer, used in GetMCHitInfo
  buffer = evbuffer;
  
  // if(!fTreeIsSet){
  //   std::cerr << "SBSSimDecoder Tree not initialized correctly - exiting" << std::endl;
  //   return HED_FATAL;
  // }
  //fTree->GetEntry(GetEvNum());
  // Cast the evbuffer pointer back to exactly the event type that is present
  // in the input file (in SBSSimFile). The pointer-to-unsigned integer is
  // needed compatibility with the standard decoder.
  if(fDebug>2)std::cout << "Processing " << here << std::endl;
  
  const SBSSimEvent* simEvent = reinterpret_cast<const SBSSimEvent*>(buffer);
  // add a check here!!!
  
  //simc variables
  fSigma_simc = simEvent->Tgmn->simc_sigma;
  fWeight_simc = simEvent->Tgmn->simc_Weight;
  fQ2_simc = simEvent->Tgmn->simc_Q2;
  fXbj_simc = simEvent->Tgmn->simc_xbj;
  fNu_simc = simEvent->Tgmn->simc_nu;
  fW_simc = simEvent->Tgmn->simc_W;
  fEpsilon_simc = simEvent->Tgmn->simc_epsilon;
  fEbeam_simc = simEvent->Tgmn->simc_Ebeam;
  fEp_simc = simEvent->Tgmn->simc_p_e;
  fEtheta_simc = simEvent->Tgmn->simc_theta_e;
  fEphi_simc = simEvent->Tgmn->simc_phi_e;
  fEPx_simc = simEvent->Tgmn->simc_px_e;
  fEPy_simc = simEvent->Tgmn->simc_py_e;
  fEPz_simc = simEvent->Tgmn->simc_pz_e;
  fFnucl_simc = simEvent->Tgmn->simc_fnucl;
  fNp_simc = simEvent->Tgmn->simc_p_n;
  fNtheta_simc = simEvent->Tgmn->simc_theta_n;
  fNphi_simc = simEvent->Tgmn->simc_phi_n;
  fNPx_simc = simEvent->Tgmn->simc_px_n;
  fNPy_simc = simEvent->Tgmn->simc_py_n;
  fNPz_simc = simEvent->Tgmn->simc_pz_n;
  fVx_simc = simEvent->Tgmn->simc_vx;
  fVy_simc = simEvent->Tgmn->simc_vy;
  fVz_simc = simEvent->Tgmn->simc_vz;
  fVeE_simc = simEvent->Tgmn->simc_veE;
  fVetheta_simc = simEvent->Tgmn->simc_vetheta;
  //g4sbs variables
  fSigma = simEvent->Tgmn->ev_sigma;
  fOmega = simEvent->Tgmn->ev_solang;
  fEPx = simEvent->Tgmn->ev_epx;
  fEPy = simEvent->Tgmn->ev_epy;
  fEPz = simEvent->Tgmn->ev_epz;
  fNPx = simEvent->Tgmn->ev_npx;
  fNPy = simEvent->Tgmn->ev_npy;
  fNPz = simEvent->Tgmn->ev_npz;
  fVx = simEvent->Tgmn->ev_vx;
  fVy = simEvent->Tgmn->ev_vy;
  fVz = simEvent->Tgmn->ev_vz;
  fEp = simEvent->Tgmn->ev_ep;
  fNp = simEvent->Tgmn->ev_np;
  fNucl = simEvent->Tgmn->ev_nucl;
  fFnucl = simEvent->Tgmn->ev_fnucl;
  fNBBtracks = simEvent->Tgmn->Earm_BBGEM_Track_ntracks;
  fBBtrack_Nhits = *(simEvent->Tgmn->Earm_BBGEM_Track_NumHits);
  fBBtrack_TID = *(simEvent->Tgmn->Earm_BBGEM_Track_TID);
  fBBtrack_PID = *(simEvent->Tgmn->Earm_BBGEM_Track_PID);
  fBBtrack_MID = *(simEvent->Tgmn->Earm_BBGEM_Track_MID);
  fBBtrack_P = *(simEvent->Tgmn->Earm_BBGEM_Track_P);
  fBBtrack_X = *(simEvent->Tgmn->Earm_BBGEM_Track_X);
  fBBtrack_Y = *(simEvent->Tgmn->Earm_BBGEM_Track_Y);
  fBBtrack_dX = *(simEvent->Tgmn->Earm_BBGEM_Track_Xp);
  fBBtrack_dY = *(simEvent->Tgmn->Earm_BBGEM_Track_Yp);
  fNBBGEMhits = simEvent->Tgmn->Earm_BBGEM_hit_nhits;
  fBBGEMhit_plane = *(simEvent->Tgmn->Earm_BBGEM_hit_plane);
  fBBGEMhit_TID = *(simEvent->Tgmn->Earm_BBGEM_hit_trid);
  fBBGEMhit_PID = *(simEvent->Tgmn->Earm_BBGEM_hit_pid);
  fBBGEMhit_MID = *(simEvent->Tgmn->Earm_BBGEM_hit_mid);
  fBBGEMhit_edep = *(simEvent->Tgmn->Earm_BBGEM_hit_edep);
  fBBGEMhit_x = *(simEvent->Tgmn->Earm_BBGEM_hit_tx);
  fBBGEMhit_y = *(simEvent->Tgmn->Earm_BBGEM_hit_ty);
  fBBPS_esum = simEvent->Tgmn->Earm_BBPSTF1_det_esum;
  fBBSH_esum = simEvent->Tgmn->Earm_BBSHTF1_det_esum;
  fBBGEMhit_ptridx = *(simEvent->Tgmn->Earm_BBGEM_hit_ptridx);
  fBBGEMhit_sdtridx = *(simEvent->Tgmn->Earm_BBGEM_hit_sdtridx);
  fBBGEMtrack_ptridx = *(simEvent->Tgmn->Earm_BBGEM_Track_ptridx);
  fBBGEMtrack_sdtridx = *(simEvent->Tgmn->Earm_BBGEM_Track_sdtridx);
  fBBHODOhit_ptridx = *(simEvent->Tgmn->Earm_BBHodoScint_hit_ptridx);
  fBBHODOhit_sdtridx = *(simEvent->Tgmn->Earm_BBHodoScint_hit_sdtridx);
  fBBPSTF1hit_ptridx = *(simEvent->Tgmn->Earm_BBPSTF1_hit_ptridx);
  fBBPSTF1hit_sdtridx = *(simEvent->Tgmn->Earm_BBPSTF1_hit_sdtridx);
  fBBSHTF1hit_ptridx = *(simEvent->Tgmn->Earm_BBSHTF1_hit_ptridx);
  fBBSHTF1hit_sdtridx = *(simEvent->Tgmn->Earm_BBSHTF1_hit_sdtridx);
  fHCALhit_ptridx = *(simEvent->Tgmn->Harm_HCalScint_hit_ptridx);
  fHCALhit_sdtridx = *(simEvent->Tgmn->Harm_HCalScint_hit_sdtridx);
  fPTrack_ntracks = simEvent->Tgmn->PTrack_ntracks;
  fPTrack_TID = *(simEvent->Tgmn->PTrack_TID);
  fPTrack_PID = *(simEvent->Tgmn->PTrack_PID);
  fPTrack_posx = *(simEvent->Tgmn->PTrack_posx);
  fPTrack_posy = *(simEvent->Tgmn->PTrack_posy);
  fPTrack_posz = *(simEvent->Tgmn->PTrack_posz);
  fPTrack_momx = *(simEvent->Tgmn->PTrack_momx);
  fPTrack_momy = *(simEvent->Tgmn->PTrack_momy);
  fPTrack_momz = *(simEvent->Tgmn->PTrack_momz);
  fPTrack_polx = *(simEvent->Tgmn->PTrack_polx);
  fPTrack_poly = *(simEvent->Tgmn->PTrack_poly);
  fPTrack_polz = *(simEvent->Tgmn->PTrack_polz);
  fPTrack_Etot = *(simEvent->Tgmn->PTrack_Etot);
  fPTrack_T = *(simEvent->Tgmn->PTrack_T);
  fSDTrack_ntracks = simEvent->Tgmn->SDTrack_ntracks;
  fSDTrack_TID = *(simEvent->Tgmn->SDTrack_TID);
  fSDTrack_MID = *(simEvent->Tgmn->SDTrack_MID);
  fSDTrack_PID = *(simEvent->Tgmn->SDTrack_PID);
  fSDTrack_posx = *(simEvent->Tgmn->SDTrack_posx);
  fSDTrack_posy = *(simEvent->Tgmn->SDTrack_posy);
  fSDTrack_posz = *(simEvent->Tgmn->SDTrack_posz);
  fSDTrack_momx = *(simEvent->Tgmn->SDTrack_momx);
  fSDTrack_momy = *(simEvent->Tgmn->SDTrack_momy);
  fSDTrack_momz = *(simEvent->Tgmn->SDTrack_momz);
  fSDTrack_polx = *(simEvent->Tgmn->SDTrack_polx);
  fSDTrack_poly = *(simEvent->Tgmn->SDTrack_poly);
  fSDTrack_polz = *(simEvent->Tgmn->SDTrack_polz);
  fSDTrack_Etot = *(simEvent->Tgmn->SDTrack_Etot);
  fSDTrack_T = *(simEvent->Tgmn->SDTrack_T);
  fSDTrack_vx = *(simEvent->Tgmn->SDTrack_vx);
  fSDTrack_vy = *(simEvent->Tgmn->SDTrack_vy);
  fSDTrack_vz = *(simEvent->Tgmn->SDTrack_vz);
  fSDTrack_vnx = *(simEvent->Tgmn->SDTrack_vnx);
  fSDTrack_vny = *(simEvent->Tgmn->SDTrack_vny);
  fSDTrack_vnz = *(simEvent->Tgmn->SDTrack_vnz);
  fSDTrack_vEkin = *(simEvent->Tgmn->SDTrack_vEkin);
  
  Int_t ret = HED_OK;
  if (first_decode || fNeedInit) {
    //fMap->print();
    if( (ret = init_cmap()) != HED_OK )
      return ret;
#if ANALYZER_VERSION_CODE < ANALYZER_VERSION(1,6,0)
    if( (ret = init_slotdata(fMap)) != HED_OK)
#else
    if( (ret = init_slotdata()) != HED_OK)
#endif
      return ret;
    first_decode = false;
  }

  if( fDoBench ) fBench->Begin("clearEvent");
  Clear();
  for(unsigned short i : fSlotClear)
    crateslot[i]->clearEvent();
  if( fDoBench ) fBench->Stop("clearEvent");

  // FIXME: needed?
  evscaler = 0;
  event_length = 0;
  
  //event_type = 1;//not smart to set event_type to 1 automatically...
  event_type = 0;//event_type set to 0 by default 
  // only set it to 1 if there is some signal in at least one detector...
  event_num = simEvent->EvtID;//++;
  //int recent_event = event_num; // no longer used

  // Event weight
  fWeight = simEvent->Tgmn->ev_sigma*simEvent->Tgmn->ev_solang;

  //
  if( fDoBench ) fBench->Begin("physics_decode");
  
  //Bool_t newclus;
  //Int_t crate, slot, chan,lchan;
  
  std::vector<std::map<Decoder::THaSlotData*, std::vector<UInt_t> > > detmaps;
  detmaps.resize(fDetectors.size());
  
  for(size_t d = 0; d<fDetectors.size(); d++){
    if(fDebug>2)cout << fDetectors[d] << endl;
    //SBSDigSim::UHitData_t* HitData_Det = simEvent->HitDataDet.at(fDetectors[d]);
    LoadDetector(detmaps[d], fDetectors[d], simEvent);
  }
  
  // Now call LoadSlot for the different detectors
  for(size_t d = 0; d < fDetectors.size(); d++) {
    //int size_det = 0;
    if(fDebug>2)
      cout << " " << fDetectors[d] << endl;
    for( std::map<Decoder::THaSlotData*, std::vector<UInt_t> >::iterator it =
	   detmaps[d].begin(); it != detmaps[d].end(); ++it) {
      if(it->first->GetModule()==0) {
        if(fDebug>2) {
	  std::cout << "No data available for detector "
		    << fDetectors[d] << std::endl;
        }
      } else {
	event_type = 1;
	//if there is data in at least one detector, event_type set to 1 
	if(fDebug>2){
	  std::cout << "load crate/slot: " << it->first->getCrate() << "/" << it->first->getSlot() << " it->second = {";
	  for(size_t k = 0; k<it->second.size(); k++)std::cout << it->second[k] << " ; ";
	  std::cout << " } " << std::endl;
	}
	it->first->GetModule()->LoadSlot(it->first,
					 it->second.data(),0,it->second.size() );
      }
      //cout << fDetectors[d].c_str() << " " << it->first->getCrate() << " " << it->first->getSlot() << " " << it->second.size() << endl;
      //size_det+=it->second.size();
    }
    //if(strcmp(fDetectors[d].c_str(), "sbs.hcal")==0)h1_sizeHCal->Fill(size_det);
    //if(strcmp(fDetectors[d].c_str(), "bb.gem")==0)h1_sizeGEMs->Fill(size_det);
  }
  
  return HED_OK;
}


//Utilities
/*
Int_t SBSSimDecoder::RetrieveDetMapParam(const char* detname, 
					  int& chanperslot, int& slotpercrate, 
					  int& firstcrate, int& firstslot)
{
  // chanperslot = ((TDetInfo &)fManager->GetDetInfo("hcal")).ChanPerSlot();
  // slotpercrate = ((TDetInfo &)fManager->GetDetInfo("hcal")).SlotPerCrate();
  // firstslot = ((TDetInfo &)fManager->GetDetInfo("hcal")).FirstSlot();
  // firstcrate = ((TDetInfo &)fManager->GetDetInfo("hcal")).FirstCrate();
  TDetInfo detinfo = fManager->GetDetInfo(detname);
  chanperslot = detinfo.ChanPerSlot();
  slotpercrate = detinfo.SlotPerCrate();
  firstslot = detinfo.FirstSlot();
  firstcrate = detinfo.FirstCrate();
}
*/


Int_t SBSSimDecoder::LoadDetector( std::map<Decoder::THaSlotData*,
				   std::vector<UInt_t> > &map,
				   const std::string& detname,
				   const SBSSimEvent* simev)
{
  if(fDebug>1)std::cout << "SBSSimDecoder::LoadDectector(" << detname << ")" << std::endl;
  //int detid = detinfo.DetUniqueId();
  Int_t crate, slot;
  //unsigned int nwords = 0;
  unsigned short chan = 0;//, data_type = 0, chan_mult = 0;
  int lchan;
  int mod, apvnum;
  //SimEncoder::mpd_data tmp_mpd;
  //UInt_t* mpd_hdr = new UInt_t[2];
  std::vector<UInt_t> strips;
  std::vector<UInt_t> samps;
  std::vector<UInt_t> times;

  bool loadevt = false;
  //int cur_apv = -1;
  
  Decoder::THaSlotData *sldat = 0;
  //This should be *general* and work for *every* subsystem
  // Loop over all raw data in this event
  //UInt_t j = 0;
  //FIXME: we don't want that, I just set it up this way for the sake of going forward
  //Simple fix (might not be ideal): do "if(detname=="xyz")"
  //cout << detname.c_str() << endl;
  int row, col;
  
  if(strcmp(detname.c_str(), "bb.ps")==0){
    //cout << " ouh " << detname.c_str() << " " << simev->Tgmn->Earm_BBPSTF1_hit_nhits << " " << simev->Tgmn->Earm_BBPS_dighit_nchan << endl;
    samps.clear();
    assert(simev->Tgmn->b_Earm_BBPS_dighit_nchan);
    for(int j = 0; j<simev->Tgmn->Earm_BBPS_dighit_nchan; j++){
      loadevt = false;
      //if(simev->Tgmn->Earm_BBPS_dighit_samp->at(j)==0)cout << "SBSSimDecoder, BBPS " << simev->Tgmn->Earm_BBPS_dighit_chan->at(j);// << endl;
      lchan = simev->Tgmn->Earm_BBPS_dighit_chan->at(j);
      //if(simev->Tgmn->Earm_BBPS_dighit_samp->at(j)==0)
      
      if(simev->Tgmn->Earm_BBPS_dighit_samp->at(j)>=0){
	samps.push_back(simev->Tgmn->Earm_BBPS_dighit_adc->at(j));
      }
      
      if(j==simev->Tgmn->Earm_BBPS_dighit_nchan-1){
	loadevt = true;
      }else if(simev->Tgmn->Earm_BBPS_dighit_chan->at(j+1)!=lchan){
	loadevt = true;
      }
      
      if(loadevt){
	/*
      for(int k = 0; k<simev->Tgmn->Earm_BBPSTF1_hit_nhits;k++){
	if(simev->Tgmn->Earm_BBPSTF1_hit_cell->at(k)==simev->Tgmn->Earm_BBPS_dighit_chan->at(j)){
	  cout << "/" << simev->Tgmn->Earm_BBPSTF1_hit_cell->at(k) << " " << simev->Tgmn->Earm_BBPSTF1_hit_row->at(k) << " " << simev->Tgmn->Earm_BBPSTF1_hit_col->at(k) << " " << simev->Tgmn->Earm_BBPSTF1_hit_xcell->at(k) << " " << simev->Tgmn->Earm_BBPSTF1_hit_ycell->at(k);// << endl;
	  break;
	}
      }
	*/
	//That stuff below is confusing... let's stick to the use of the block position in the DB!
	row = lchan%26;
	col = (lchan-row)/26;
	lchan = row*2+col;
	//row = 25-row;
	//lchan = col*26+row;
	//cout << " => " << row << ", " << col << " new lchan = " << lchan << endl;
	//ADC
	ChanToROC(detname, lchan, crate, slot, chan);
	
	if( crate >= 0 || slot >=  0 ) {
	  sldat = crateslot[idx(crate,slot)].get();
	}
	std::vector<UInt_t> *myev = &(map[sldat]);
	
	// cout << detname.c_str() << " det channel " << lchan << ", crate " << crate 
	//      << ", slot " << slot << " chan " << chan << " size " << samps.size() << endl;
	if(!samps.empty()){
	  myev->push_back(SBSSimDataDecoder::EncodeHeader(5, chan, samps.size()));
	  for(unsigned int samp : samps){
	    myev->push_back(samp);
	    //cout << " " << samps[k];
	  }
	}
	//cout << endl;
	
	samps.clear();
      }
      
      /*
      //cout << j << " " << simev->Tgmn->Earm_BBPS_dighit_chan->at(j) << " " << simev->Tgmn->Earm_BBPS_dighit_adc->at(j) << endl;
      lchan = simev->Tgmn->Earm_BBPS_dighit_chan->at(j);
      ChanToROC(detname, lchan, crate, slot, chan);
      
      if( crate >= 0 || slot >=  0 ) {
	sldat = crateslot[idx(crate,slot)].get();
      }
      std::vector<UInt_t> *myev = &(map[sldat]);
      
      myev->push_back(SBSSimDataDecoder::EncodeHeader(6, chan, 1));
   
      myev->push_back(simev->Tgmn->Earm_BBPS_dighit_adc->at(j));
      
      if(fDebug>2){
	std::cout << " j = " << j << " my ev = {";
	for(size_t k = 0; k<myev->size(); k++)std::cout << myev->at(k) << " ; ";
	std::cout << " } " << std::endl;
      }
      */
    }
  }
  if(strcmp(detname.c_str(), "bb.sh")==0){
    //cout << " ouh " << detname.c_str() << " " << simev->Tgmn->Earm_BBSHTF1_hit_nhits << " " << simev->Tgmn->Earm_BBSH_dighit_nchan << endl;
    samps.clear();
    assert(simev->Tgmn->b_Earm_BBSH_dighit_nchan);
    for(int j = 0; j<simev->Tgmn->Earm_BBSH_dighit_nchan; j++){
      loadevt = false;
      //if(simev->Tgmn->Earm_BBSH_dighit_samp->at(j)==0)cout << "SBSSimDecoder, BBSH " << simev->Tgmn->Earm_BBSH_dighit_chan->at(j);// << endl;
      lchan = simev->Tgmn->Earm_BBSH_dighit_chan->at(j);
      if(simev->Tgmn->Earm_BBSH_dighit_samp->at(j)>=0){
	samps.push_back(simev->Tgmn->Earm_BBSH_dighit_adc->at(j));
      }
      
      if(j==simev->Tgmn->Earm_BBSH_dighit_nchan-1){
	loadevt = true;
      }else if(simev->Tgmn->Earm_BBSH_dighit_chan->at(j+1)!=lchan){
	loadevt = true;
      }
     
      if(loadevt){
	/*
      for(int k = 0; k<simev->Tgmn->Earm_BBSHTF1_hit_nhits;k++){
	if(simev->Tgmn->Earm_BBSHTF1_hit_cell->at(k)==simev->Tgmn->Earm_BBSH_dighit_chan->at(j)){
	  cout << " " << simev->Tgmn->Earm_BBSHTF1_hit_cell->at(k) << " " << simev->Tgmn->Earm_BBSHTF1_hit_row->at(k) << " " << simev->Tgmn->Earm_BBSHTF1_hit_col->at(k) << " " << simev->Tgmn->Earm_BBSHTF1_hit_xcell->at(k) << " " << simev->Tgmn->Earm_BBSHTF1_hit_ycell->at(k);// << endl;
	  break;
	}
      }
	*/
	row = lchan%27;
	col = (lchan-row)/27;
      // row = 26-row;
      // col = 6-col;
	lchan = row*7+col;
      //cout << " => " << row << ", " << col << " new lchan = " << lchan << endl;
	//ADC
	ChanToROC(detname, lchan, crate, slot, chan);
	
	if( crate >= 0 || slot >=  0 ) {
	  sldat = crateslot[idx(crate,slot)].get();
	}
	std::vector<UInt_t> *myev = &(map[sldat]);
	
	// cout << detname.c_str() << " det channel " << lchan << ", crate " << crate 
	//      << ", slot " << slot << " chan " << chan << " size " << samps.size() << endl;
	if(!samps.empty()){
	  myev->push_back(SBSSimDataDecoder::EncodeHeader(5, chan, samps.size()));
	  for(unsigned int samp : samps){
	    myev->push_back(samp);
	    //cout << " " << samps[k];
	  }
	}
	//cout << endl;
	
	samps.clear();
      }
      
      /*
      //cout << j << " " << simev->Tgmn->Earm_BBSH_dighit_chan->at(j) << " " << simev->Tgmn->Earm_BBSH_dighit_adc->at(j) << endl;
      lchan = simev->Tgmn->Earm_BBSH_dighit_chan->at(j);
      ChanToROC(detname, lchan, crate, slot, chan);
      
      if( crate >= 0 || slot >=  0 ) {
	sldat = crateslot[idx(crate,slot)].get();
      }
      std::vector<UInt_t> *myev = &(map[sldat]);
      
      myev->push_back(SBSSimDataDecoder::EncodeHeader(6, chan, 1));
   
      myev->push_back(simev->Tgmn->Earm_BBSH_dighit_adc->at(j));
      
      if(fDebug>2){
	std::cout << " j = " << j << " my ev = {";
	for(size_t k = 0; k<myev->size(); k++)std::cout << myev->at(k) << " ; ";
	std::cout << " } " << std::endl;
      }
      */
    }
  }
  if(strcmp(detname.c_str(), "bb.hodo")==0){
    //cout << " ouh " << detname.c_str() << " " << simev->Tgmn->Earm_BBHodoScint_hit_nhits << " " << simev->Tgmn->Earm_BBHodo_dighit_nchan << endl;
    // cout << simev->Tgmn->Earm_BBHodo_dighit_chan->size() << " " 
    // 	 << simev->Tgmn->Earm_BBHodo_dighit_adc->size() << " " 
    // 	 << simev->Tgmn->Earm_BBHodo_dighit_tdc_l->size() << " " 
    // 	 << simev->Tgmn->Earm_BBHodo_dighit_tdc_t->size() << endl; 
    /*
    ChanToROC(detname, 180, crate, slot, chan);
    cout << crate << " " << slot << " " << chan << endl;
    if( crate >= 0 || slot >=  0 ) {
      sldat = crateslot[idx(crate,slot)].get();
    }
    std::vector<UInt_t> *myev = &(map[sldat]);
    myev->push_back(SBSSimDataDecoder::EncodeHeader(1, chan, 2));
    myev->push_back(0);
    */
    int ntdc = 0;
    assert(simev->Tgmn->b_Earm_BBHodo_dighit_nchan);
    for(int j = 0; j<simev->Tgmn->Earm_BBHodo_dighit_nchan; j++){
      ntdc = 0;
      lchan = simev->Tgmn->Earm_BBHodo_dighit_chan->at(j);
      col = lchan%2;
      row = (lchan-col)/2;
      lchan = col*90+row;
      ChanToROC(detname, lchan, crate, slot, chan);
      //cout << detname << " " << simev->Tgmn->Earm_BBHodo_dighit_chan->at(j) << " " << lchan << " " << crate << " " << slot << " " << chan << endl;
      //cout << j << " " << simev->Tgmn->Earm_BBHodo_dighit_chan->at(j) << " " << simev->Tgmn->Earm_BBHodo_dighit_adc->at(j) << " " << simev->Tgmn->Earm_BBHodo_dighit_tdc_l->at(j) << " " << simev->Tgmn->Earm_BBHodo_dighit_tdc_t->at(j) << endl;
      if( crate >= 0 || slot >=  0 ) {
	sldat = crateslot[idx(crate,slot)].get();
      }
      if(simev->Tgmn->Earm_BBHodo_dighit_tdc_l->at(j)>-1000000)ntdc++;
      if(simev->Tgmn->Earm_BBHodo_dighit_tdc_t->at(j)>-1000000)ntdc++;
      
      if(ntdc){
	std::vector<UInt_t> *myev = &(map[sldat]);
	myev->push_back(SBSSimDataDecoder::EncodeHeader(1, chan, ntdc));
	
	if(simev->Tgmn->Earm_BBHodo_dighit_tdc_l->at(j)>-1000000)myev->push_back(simev->Tgmn->Earm_BBHodo_dighit_tdc_l->at(j));
	if(simev->Tgmn->Earm_BBHodo_dighit_tdc_t->at(j)>-1000000){
	  uint tdc =  simev->Tgmn->Earm_BBHodo_dighit_tdc_t->at(j)|(1<<31);
	  //cout << tdc << endl;
	  myev->push_back( tdc );
	}
      /*
      ChanToROC(detname, lchan, crate, slot, chan);//+91 ??? that might be the trick
      if( crate >= 0 || slot >=  0 ) {
	sldat = crateslot[idx(crate,slot)].get();
      }
      myev = &(map[sldat]);
      
      myev->push_back(SBSSimDataDecoder::EncodeHeader(8, chan, 1));
      myev->push_back(simev->Tgmn->Earm_BBHodo_dighit_adc->at(j));
      */
	if(fDebug>2){
	  std::cout << " j = " << j << " my ev = {";
	  for(size_t k = 0; k<myev->size(); k++)std::cout << myev->at(k) << " ; ";
	  std::cout << " } " << std::endl;
	}
      }
    }
  }
  if(strcmp(detname.c_str(), "bb.grinch_tdc")==0){
    int ntdc = 0;
    //if(simev->Tgmn->b_Earm_GRINCH_dighit_nchan==0)
    //cout << "*** Warning: your GRINCH variables are probably missing in the tree you are analyzing. " << endl << " consider using another file or removing the grinch for your analysis " << endl;
    //cout << " ouh " << detname.c_str() << " " << simev->Tgmn->Earm_GRINCH_hit_nhits << " " << simev->Tgmn->b_Earm_GRINCH_dighit_nchan << " " << simev->Tgmn->Earm_GRINCH_dighit_nchan << endl;
    assert(simev->Tgmn->b_Earm_GRINCH_dighit_nchan);
    for(int j = 0; j<simev->Tgmn->Earm_GRINCH_dighit_nchan; j++){
      ntdc = 0;
      //cout << j << " " << simev->Tgmn->Earm_GRINCH_dighit_chan->at(j) << " " << simev->Tgmn->Earm_GRINCH_dighit_adc->at(j) << " " << simev->Tgmn->Earm_GRINCH_dighit_tdc_l->at(j) << " " << simev->Tgmn->Earm_GRINCH_dighit_tdc_t->at(j) << endl;
      lchan = simev->Tgmn->Earm_GRINCH_dighit_chan->at(j);
      ChanToROC(detname, lchan, crate, slot, chan);
      
      if( crate >= 0 || slot >=  0 ) {
	sldat = crateslot[idx(crate,slot)].get();
      }

      if(simev->Tgmn->Earm_GRINCH_dighit_tdc_l->at(j)>-1000000)ntdc++;
      if(simev->Tgmn->Earm_GRINCH_dighit_tdc_t->at(j)>-1000000)ntdc++;

      if(ntdc){
	std::vector<UInt_t> *myev = &(map[sldat]);
	
	myev->push_back(SBSSimDataDecoder::EncodeHeader(1, chan, ntdc));
	
	if(simev->Tgmn->Earm_GRINCH_dighit_tdc_l->at(j)>-1000000)myev->push_back(simev->Tgmn->Earm_GRINCH_dighit_tdc_l->at(j));
	if(simev->Tgmn->Earm_GRINCH_dighit_tdc_t->at(j)>-1000000){
	  uint tdc =  simev->Tgmn->Earm_GRINCH_dighit_tdc_t->at(j)|(1<<31);
	  //cout << tdc << endl;
	  myev->push_back( tdc );
	}
      /*
      ChanToROC(detname, lchan, crate, slot, chan);//+288 ??? that might be the trick
      if( crate >= 0 || slot >=  0 ) {
	sldat = crateslot[idx(crate,slot)].get();
      }
      myev = &(map[sldat]);
      
      myev->push_back(SBSSimDataDecoder::EncodeHeader(8, chan, 1));
      myev->push_back(simev->Tgmn->Earm_GRINCH_dighit_adc->at(j));
      */
	if(fDebug>2){
	  std::cout << " j = " << j << " my ev = {";
	  for(size_t k = 0; k<myev->size(); k++)std::cout << myev->at(k) << " ; ";
	  std::cout << " } " << std::endl;
	}
      }
    }
  }
  
  if(strcmp(detname.c_str(), "bb.gem")==0){
    //cout << fPx << " " << fPy << " " << fPz << "   " << fVz << endl;
    samps.clear();  
    strips.clear();  
    //cout << " ouh " << detname.c_str() << " " << simev->Tgmn->Earm_BBGEM_dighit_nstrips << endl;
    assert(simev->Tgmn->b_Earm_BBGEM_dighit_nstrips);
    for(int j = 0; j<simev->Tgmn->Earm_BBGEM_dighit_nstrips; j++){
      loadevt = false;
      mod = simev->Tgmn->Earm_BBGEM_dighit_module->at(j);
      lchan = simev->Tgmn->Earm_BBGEM_dighit_strip->at(j);
      apvnum = APVnum(detname, mod, lchan, crate, slot, chan);
      
      if(simev->Tgmn->Earm_BBGEM_dighit_samp->at(j)>=0){
	strips.push_back(chan);
	samps.push_back(simev->Tgmn->Earm_BBGEM_dighit_adc->at(j));
      }
      
      if(fDebug>3)
	cout << " mod " << mod << " lchan " << lchan << " crate " << crate << " slot " << slot << " apvnum " << apvnum << " chan " << chan << " samp " << simev->Tgmn->Earm_BBGEM_dighit_samp->at(j)  << " adc " << simev->Tgmn->Earm_BBGEM_dighit_adc->at(j) << endl;
      //if(mod>=26 && simev->Tgmn->Earm_BBGEM_dighit_samp->at(j)==5)cout << mod << " " << lchan << " " << apvnum << endl;
      
      if(j==simev->Tgmn->Earm_BBGEM_dighit_nstrips-1){
	loadevt = true;
      }else if(mod!=simev->Tgmn->Earm_BBGEM_dighit_module->at(j+1) ||
	       //fabs(lchan-simev->Tgmn->Earm_BBGEM_dighit_strip->at(j+1))>=128
	       floor(simev->Tgmn->Earm_BBGEM_dighit_strip->at(j+1)/128)!=floor(lchan/128)
	       ){
	loadevt = true;
      }
	
      if(loadevt){
	if( crate >= 0 || slot >=  0 ) {
	  sldat = crateslot[idx(crate,slot)].get();
	}
	std::vector<UInt_t> *myev = &(map[sldat]);
	
	if(!samps.empty()){
	  //myev->push_back(SBSSimDataDecoder::EncodeHeader(5, apvnum, samps.size()));
	  //I think I'm onto something here, but I also need to transmit strip num 
	  myev->push_back(SBSSimDataDecoder::EncodeHeader(9, apvnum, samps.size()));
	  for(int k = 0; k<(int)samps.size(); k++){
	    // cout << " " << samps[k];
	    myev->push_back(strips[k]*8192+samps[k]);//strips[k]<< 13 | samps[k]);
	  }
	  //for(int l = 0; l<myev->size();l++)cout << myev->at(l) << " ";
	  //cout << endl;
	}
	//cout << endl;
	
	samps.clear();
	strips.clear();
      }
    }
  }
  
  
  if(strcmp(detname.c_str(), "sbs.hcal")==0){
    //cout << " ouh " << detname.c_str() << " " << simev->Tgmn->Harm_HCalScint_hit_nhits << " " << simev->Tgmn->Harm_HCal_dighit_nchan << endl;
    samps.clear();
    times.clear();
    
    assert(simev->Tgmn->b_Harm_HCal_dighit_nchan);
    for(int j = 0; j<simev->Tgmn->Harm_HCal_dighit_nchan; j++){
      loadevt = false;
      lchan = simev->Tgmn->Harm_HCal_dighit_chan->at(j);
      if(simev->Tgmn->Harm_HCal_dighit_samp->at(j)>=0){
	samps.push_back(simev->Tgmn->Harm_HCal_dighit_adc->at(j));
      }else{
	times.push_back(simev->Tgmn->Harm_HCal_dighit_tdc->at(j));
      }
      
      if(j==simev->Tgmn->Harm_HCal_dighit_nchan-1){
	loadevt = true;
      }else if(simev->Tgmn->Harm_HCal_dighit_chan->at(j+1)!=lchan){
	loadevt = true;
      }

      // In simulation row 0 col 0 block starts at top left corner weheras in real data row 0 col 0 starts at top
      // right corner, while looking at HCAL from front. Lets try the following to eleminate the mismatch:
      col = lchan%12;
      row = (lchan-col)/12; // row in simulation is already same as real data
      col = 12 - 1 - col; // this will fix the mismatch in column numbering
      lchan = row*12 + col; 
      // --
      
      if(loadevt){
	//ADC
	ChanToROC(detname, lchan, crate, slot, chan);
	
	if( crate >= 0 || slot >=  0 ) {
	  sldat = crateslot[idx(crate,slot)].get();
	}
	std::vector<UInt_t> *myev = &(map[sldat]);
	
	// cout << detname.c_str() << " det channel " << lchan << ", crate " << crate 
	//      << ", slot " << slot << " chan " << chan << " size " << samps.size() << endl;
	if(!samps.empty()){
	  myev->push_back(SBSSimDataDecoder::EncodeHeader(5, chan, samps.size()));
	  for(unsigned int samp : samps){
	    myev->push_back(samp);
	    //cout << " " << samps[k];
	  }
	}
	//cout << endl;

	//TDC
	ChanToROC(detname, lchan+288, crate, slot, chan);
	if( crate >= 0 || slot >=  0 ) {
	  sldat = crateslot[idx(crate,slot)].get();
	}
	myev = &(map[sldat]);
	if(!times.empty()){
	  myev->push_back(SBSSimDataDecoder::EncodeHeader(4, chan, times.size()));
	  for(unsigned int time : times){
	    myev->push_back(time);
	  }
	}
	
	samps.clear();
	times.clear();
      }
      
      /*
      ChanToROC(detname, lchan, crate, slot, chan);
      //cout << lchan << " " << crate << " " << slot << " " << chan << endl;
      
      if( crate >= 0 || slot >=  0 ) {
	sldat = crateslot[idx(crate,slot)].get();
      }
      std::vector<UInt_t> *myev = &(map[sldat]);
      
      //cout << SBSSimDataDecoder::EncodeHeader(5, chan, 20) << endl;
      //cout << SBSSimDataDecoder::EncodeHeader(5, chan, 1) << endl;
      myev->push_back(SBSSimDataDecoder::EncodeHeader(5, chan, 20));
      myev->push_back(simev->Tgmn->Harm_HCal_dighit_adc_0->at(j));
      myev->push_back(simev->Tgmn->Harm_HCal_dighit_adc_1->at(j));
      myev->push_back(simev->Tgmn->Harm_HCal_dighit_adc_2->at(j));
      myev->push_back(simev->Tgmn->Harm_HCal_dighit_adc_3->at(j));
      myev->push_back(simev->Tgmn->Harm_HCal_dighit_adc_4->at(j));
      myev->push_back(simev->Tgmn->Harm_HCal_dighit_adc_5->at(j));
      myev->push_back(simev->Tgmn->Harm_HCal_dighit_adc_6->at(j));
      myev->push_back(simev->Tgmn->Harm_HCal_dighit_adc_7->at(j));
      myev->push_back(simev->Tgmn->Harm_HCal_dighit_adc_8->at(j));
      myev->push_back(simev->Tgmn->Harm_HCal_dighit_adc_9->at(j));
      myev->push_back(simev->Tgmn->Harm_HCal_dighit_adc_10->at(j));
      myev->push_back(simev->Tgmn->Harm_HCal_dighit_adc_11->at(j));
      myev->push_back(simev->Tgmn->Harm_HCal_dighit_adc_12->at(j));
      myev->push_back(simev->Tgmn->Harm_HCal_dighit_adc_13->at(j));
      myev->push_back(simev->Tgmn->Harm_HCal_dighit_adc_14->at(j));
      myev->push_back(simev->Tgmn->Harm_HCal_dighit_adc_15->at(j));
      myev->push_back(simev->Tgmn->Harm_HCal_dighit_adc_16->at(j));
      myev->push_back(simev->Tgmn->Harm_HCal_dighit_adc_17->at(j));
      myev->push_back(simev->Tgmn->Harm_HCal_dighit_adc_18->at(j));
      myev->push_back(simev->Tgmn->Harm_HCal_dighit_adc_19->at(j));
      ChanToROC(detname, lchan+288, crate, slot, chan);//+288 ??? that might be the trick

      //cout << lchan+288  << " " << crate << " " << slot << " " << chan << endl;
      if( crate >= 0 || slot >=  0 ) {
	sldat = crateslot[idx(crate,slot)].get();
      }
      myev = &(map[sldat]);
      
      myev->push_back(SBSSimDataDecoder::EncodeHeader(4, chan, 1));
      myev->push_back(simev->Tgmn->Harm_HCal_dighit_tdc->at(j));
      */
    }

  }

  //add here the GEN-RP scintillators
  if(strcmp(detname.c_str(), "sbs.active_ana")==0){
    //cout << " ouh " << detname.c_str() << " " << simev->Tgenrp->Earm_BBSHTF1_hit_nhits << " " << simev->Tgenrp->Earm_BBSH_dighit_nchan << endl;
    samps.clear();
    assert(simev->Tgenrp->b_Harm_ActAn_dighit_nchan);
    for(int j = 0; j<simev->Tgenrp->Harm_ActAn_dighit_nchan; j++){
      loadevt = false;
      //if(simev->Tgenrp->Harm_ActAn_dighit_samp->at(j)==0)cout << "SBSSimDecoder, BBSH " << simev->Tgenrp->Harm_ActAn_dighit_chan->at(j);// << endl;
      lchan = simev->Tgenrp->Harm_ActAn_dighit_chan->at(j);
      if(simev->Tgenrp->Harm_ActAn_dighit_samp->at(j)>=0){
	samps.push_back(simev->Tgenrp->Harm_ActAn_dighit_adc->at(j));
      }
      
      if(j==simev->Tgenrp->Harm_ActAn_dighit_nchan-1){
	loadevt = true;
      }else if(simev->Tgenrp->Harm_ActAn_dighit_chan->at(j+1)!=lchan){
	loadevt = true;
      }
     
      if(loadevt){
	/*
      for(int k = 0; k<simev->Tgenrp->Earm_BBSHTF1_hit_nhits;k++){
	if(simev->Tgenrp->Harm_ActAnScint_hit_cell->at(k)==simev->Tgenrp->Harm_ActAn_dighit_chan->at(j)){
	  cout << " " << simev->Tgenrp->Harm_ActAnScint_hit_cell->at(k) << " " << simev->Tgenrp->Harm_ActAnScint_hit_row->at(k) << " " << simev->Tgenrp->Harm_ActAnScint_hit_col->at(k) << " " << simev->Tgenrp->Harm_ActAnScint_hit_xcell->at(k) << " " << simev->Tgenrp->Harm_ActAnScint_hit_ycell->at(k);// << endl;
	  break;
	}
      }
	*/
	//row = lchan%4;
	//col = (lchan-row)/4;
	//lchan = row*4+col;
	//ADC
	ChanToROC(detname, lchan, crate, slot, chan);
	//if(crate!=9)cout << detname << " " << simev->Tgenrp->Harm_ActAn_dighit_chan->at(j) << " " << lchan << " " << crate << " " << slot << endl;
	
	if( crate >= 0 || slot >=  0 ) {
	  sldat = crateslot[idx(crate,slot)].get();
	}
	std::vector<UInt_t> *myev = &(map[sldat]);
	
	// cout << detname.c_str() << " det channel " << lchan << ", crate " << crate 
	//      << ", slot " << slot << " chan " << chan << " size " << samps.size() << endl;
	if(!samps.empty()){
	  myev->push_back(SBSSimDataDecoder::EncodeHeader(5, chan, samps.size()));
	  for(unsigned int samp : samps){
	    myev->push_back(samp);
	    //cout << " " << samps[k];
	  }
	}
	//cout << endl;
	
	samps.clear();
      }
      
      /*
      //cout << j << " " << simev->Tgenrp->Harm_ActAn_dighit_chan->at(j) << " " << simev->Tgenrp->Harm_ActAn_dighit_adc->at(j) << endl;
      lchan = simev->Tgenrp->Harm_ActAn_dighit_chan->at(j);
      ChanToROC(detname, lchan, crate, slot, chan);
      
      if( crate >= 0 || slot >=  0 ) {
	sldat = crateslot[idx(crate,slot)].get();
      }
      std::vector<UInt_t> *myev = &(map[sldat]);
      
      myev->push_back(SBSSimDataDecoder::EncodeHeader(6, chan, 1));
   
      myev->push_back(simev->Tgenrp->Harm_ActAn_dighit_adc->at(j));
      
      if(fDebug>2){
	std::cout << " j = " << j << " my ev = {";
	for(size_t k = 0; k<myev->size(); k++)std::cout << myev->at(k) << " ; ";
	std::cout << " } " << std::endl;
      }
      */
    }
  }
  if(strcmp(detname.c_str(), "sbs.prpolscint_farside")==0){
    //cout << " ouh " << detname.c_str() << " " << simev->Tgenrp->Earm_BBHodoScint_hit_nhits << " " << simev->Tgenrp->Harm_PRPolScintFarSide_dighit_nchan << endl;
    // cout << simev->Tgenrp->Harm_PRPolScintFarSide_dighit_chan->size() << " " 
    // 	 << simev->Tgenrp->Harm_PRPolScintFarSide_dighit_adc->size() << " " 
    // 	 << simev->Tgenrp->Harm_PRPolScintFarSide_dighit_tdc_l->size() << " " 
    // 	 << simev->Tgenrp->Harm_PRPolScintFarSide_dighit_tdc_t->size() << endl; 
    /*
    ChanToROC(detname, 180, crate, slot, chan);
    cout << crate << " " << slot << " " << chan << endl;
    if( crate >= 0 || slot >=  0 ) {
      sldat = crateslot[idx(crate,slot)].get();
    }
    std::vector<UInt_t> *myev = &(map[sldat]);
    myev->push_back(SBSSimDataDecoder::EncodeHeader(1, chan, 2));
    myev->push_back(0);
    */
    int ntdc = 0;
    assert(simev->Tgenrp->b_Harm_PRPolScintFarSide_dighit_nchan);
    for(int j = 0; j<simev->Tgenrp->Harm_PRPolScintFarSide_dighit_nchan; j++){
      ntdc = 0;
      lchan = simev->Tgenrp->Harm_PRPolScintFarSide_dighit_chan->at(j);
      //do we want that???
      //col = lchan%2;
      //row = (lchan-col)/2;
      //lchan = col*24+row;
      ChanToROC(detname, lchan, crate, slot, chan);
      //if(crate!=9)cout << detname << " " << simev->Tgenrp->Harm_PRPolScintFarSide_dighit_chan->at(j) << " " << lchan << " " << crate << " " << slot << endl;
      if( crate >= 0 || slot >=  0 ) {
	sldat = crateslot[idx(crate,slot)].get();
      }
      if(simev->Tgenrp->Harm_PRPolScintFarSide_dighit_tdc_l->at(j)>-1000000)ntdc++;
      if(simev->Tgenrp->Harm_PRPolScintFarSide_dighit_tdc_t->at(j)>-1000000)ntdc++;
      
      if(ntdc){
	std::vector<UInt_t> *myev = &(map[sldat]);
	myev->push_back(SBSSimDataDecoder::EncodeHeader(1, chan, ntdc));
	
	if(simev->Tgenrp->Harm_PRPolScintFarSide_dighit_tdc_l->at(j)>-1000000)myev->push_back(simev->Tgenrp->Harm_PRPolScintFarSide_dighit_tdc_l->at(j));
	if(simev->Tgenrp->Harm_PRPolScintFarSide_dighit_tdc_t->at(j)>-1000000){
	  uint tdc =  simev->Tgenrp->Harm_PRPolScintFarSide_dighit_tdc_t->at(j)|(1<<31);
	  //cout << tdc << endl;
	  myev->push_back( tdc );
	}
      /*
      ChanToROC(detname, lchan, crate, slot, chan);//+91 ??? that might be the trick
      if( crate >= 0 || slot >=  0 ) {
	sldat = crateslot[idx(crate,slot)].get();
      }
      myev = &(map[sldat]);
      
      myev->push_back(SBSSimDataDecoder::EncodeHeader(8, chan, 1));
      myev->push_back(simev->Tgenrp->Harm_PRPolScintFarSide_dighit_adc->at(j));
      */
	if(fDebug>2){
	  std::cout << " j = " << j << " my ev = {";
	  for(size_t k = 0; k<myev->size(); k++)std::cout << myev->at(k) << " ; ";
	  std::cout << " } " << std::endl;
	}
      }
    }
  }
  
  if(strcmp(detname.c_str(), "sbs.cepolfront_gem")==0){
    //cout << fPx << " " << fPy << " " << fPz << "   " << fVz << endl;
    samps.clear();  
    strips.clear();  
    //cout << " ouh " << detname.c_str() << " " << simev->Tgenrp->Harm_CEPolFront_dighit_nstrips << endl;
    assert(simev->Tgenrp->b_Harm_CEPolFront_dighit_nstrips);
    for(int j = 0; j<simev->Tgenrp->Harm_CEPolFront_dighit_nstrips; j++){
      loadevt = false;
      mod = simev->Tgenrp->Harm_CEPolFront_dighit_module->at(j);
      lchan = simev->Tgenrp->Harm_CEPolFront_dighit_strip->at(j);
      apvnum = APVnum(detname, mod, lchan, crate, slot, chan);
      
      if(simev->Tgenrp->Harm_CEPolFront_dighit_samp->at(j)>=0){
	strips.push_back(chan);
	samps.push_back(simev->Tgenrp->Harm_CEPolFront_dighit_adc->at(j));
      }
      
      if(fDebug>3)
	cout << " mod " << mod << " lchan " << lchan << " crate " << crate << " slot " << slot << " apvnum " << apvnum << " chan " << chan << " samp " << simev->Tgenrp->Harm_CEPolFront_dighit_samp->at(j)  << " adc " << simev->Tgenrp->Harm_CEPolFront_dighit_adc->at(j) << endl;
      //if(mod>=26 && simev->Tgenrp->Harm_CEPolFront_dighit_samp->at(j)==5)cout << mod << " " << lchan << " " << apvnum << endl;
      
      if(j==simev->Tgenrp->Harm_CEPolFront_dighit_nstrips-1){
	loadevt = true;
      }else if(mod!=simev->Tgenrp->Harm_CEPolFront_dighit_module->at(j+1) ||
	       //fabs(lchan-simev->Tgenrp->Harm_CEPolFront_dighit_strip->at(j+1))>=128
	       floor(simev->Tgenrp->Harm_CEPolFront_dighit_strip->at(j+1)/128)!=floor(lchan/128)
	       ){
	loadevt = true;
      }
	
      if(loadevt){
	if( crate >= 0 || slot >=  0 ) {
	  sldat = crateslot[idx(crate,slot)].get();
	}
	std::vector<UInt_t> *myev = &(map[sldat]);
	
	if(!samps.empty()){
	  //myev->push_back(SBSSimDataDecoder::EncodeHeader(5, apvnum, samps.size()));
	  //I think I'm onto something here, but I also need to transmit strip num 
	  myev->push_back(SBSSimDataDecoder::EncodeHeader(9, apvnum, samps.size()));
	  for(int k = 0; k<(int)samps.size(); k++){
	    // cout << " " << samps[k];
	    myev->push_back(strips[k]*8192+samps[k]);//strips[k]<< 13 | samps[k]);
	  }
	  //for(int l = 0; l<myev->size();l++)cout << myev->at(l) << " ";
	  //cout << endl;
	}
	//cout << endl;
	
	samps.clear();
	strips.clear();
      }
    }
  }

  if(strcmp(detname.c_str(), "sbs.cepolrear_gem")==0){
    //cout << fPx << " " << fPy << " " << fPz << "   " << fVz << endl;
    samps.clear();  
    strips.clear();  
    //cout << " ouh " << detname.c_str() << " " << simev->Tgenrp->Harm_CEPolRear_dighit_nstrips << endl;
    assert(simev->Tgenrp->b_Harm_CEPolRear_dighit_nstrips);
    for(int j = 0; j<simev->Tgenrp->Harm_CEPolRear_dighit_nstrips; j++){
      loadevt = false;
      mod = simev->Tgenrp->Harm_CEPolRear_dighit_module->at(j);
      lchan = simev->Tgenrp->Harm_CEPolRear_dighit_strip->at(j);
      apvnum = APVnum(detname, mod, lchan, crate, slot, chan);
      
      if(simev->Tgenrp->Harm_CEPolRear_dighit_samp->at(j)>=0){
	strips.push_back(chan);
	samps.push_back(simev->Tgenrp->Harm_CEPolRear_dighit_adc->at(j));
      }
      
      if(fDebug>3)
	cout << " mod " << mod << " lchan " << lchan << " crate " << crate << " slot " << slot << " apvnum " << apvnum << " chan " << chan << " samp " << simev->Tgenrp->Harm_CEPolRear_dighit_samp->at(j)  << " adc " << simev->Tgenrp->Harm_CEPolRear_dighit_adc->at(j) << endl;
      //if(mod>=26 && simev->Tgenrp->Harm_CEPolRear_dighit_samp->at(j)==5)cout << mod << " " << lchan << " " << apvnum << endl;
      
      if(j==simev->Tgenrp->Harm_CEPolRear_dighit_nstrips-1){
	loadevt = true;
      }else if(mod!=simev->Tgenrp->Harm_CEPolRear_dighit_module->at(j+1) ||
	       //fabs(lchan-simev->Tgenrp->Harm_CEPolRear_dighit_strip->at(j+1))>=128
	       floor(simev->Tgenrp->Harm_CEPolRear_dighit_strip->at(j+1)/128)!=floor(lchan/128)
	       ){
	loadevt = true;
      }
	
      if(loadevt){
	if( crate >= 0 || slot >=  0 ) {
	  sldat = crateslot[idx(crate,slot)].get();
	}
	std::vector<UInt_t> *myev = &(map[sldat]);
	
	if(!samps.empty()){
	  //myev->push_back(SBSSimDataDecoder::EncodeHeader(5, apvnum, samps.size()));
	  //I think I'm onto something here, but I also need to transmit strip num 
	  myev->push_back(SBSSimDataDecoder::EncodeHeader(9, apvnum, samps.size()));
	  for(int k = 0; k<(int)samps.size(); k++){
	    // cout << " " << samps[k];
	    myev->push_back(strips[k]*8192+samps[k]);//strips[k]<< 13 | samps[k]);
	  }
	  //for(int l = 0; l<myev->size();l++)cout << myev->at(l) << " ";
	  //cout << endl;
	}
	//cout << endl;
	
	samps.clear();
	strips.clear();
      }
    }
  }

  if(strcmp(detname.c_str(), "sbs.prpolgem_farside")==0){
    //cout << fPx << " " << fPy << " " << fPz << "   " << fVz << endl;
    samps.clear();  
    strips.clear();  
    //cout << " ouh " << detname.c_str() << " " << simev->Tgenrp->Harm_PRPolGEMFarSide_dighit_nstrips << endl;
    assert(simev->Tgenrp->b_Harm_PRPolGEMFarSide_dighit_nstrips);
    for(int j = 0; j<simev->Tgenrp->Harm_PRPolGEMFarSide_dighit_nstrips; j++){
      loadevt = false;
      mod = simev->Tgenrp->Harm_PRPolGEMFarSide_dighit_module->at(j);
      lchan = simev->Tgenrp->Harm_PRPolGEMFarSide_dighit_strip->at(j);
      apvnum = APVnum(detname, mod, lchan, crate, slot, chan);
      
      if(simev->Tgenrp->Harm_PRPolGEMFarSide_dighit_samp->at(j)>=0){
	strips.push_back(chan);
	samps.push_back(simev->Tgenrp->Harm_PRPolGEMFarSide_dighit_adc->at(j));
      }
      
      if(fDebug>3)
	cout << " mod " << mod << " lchan " << lchan << " crate " << crate << " slot " << slot << " apvnum " << apvnum << " chan " << chan << " samp " << simev->Tgenrp->Harm_PRPolGEMFarSide_dighit_samp->at(j)  << " adc " << simev->Tgenrp->Harm_PRPolGEMFarSide_dighit_adc->at(j) << endl;
      //if(mod>=26 && simev->Tgenrp->Harm_PRPolGEMFarSide_dighit_samp->at(j)==5)cout << mod << " " << lchan << " " << apvnum << endl;
      
      if(j==simev->Tgenrp->Harm_PRPolGEMFarSide_dighit_nstrips-1){
	loadevt = true;
      }else if(mod!=simev->Tgenrp->Harm_PRPolGEMFarSide_dighit_module->at(j+1) ||
	       //fabs(lchan-simev->Tgenrp->Harm_PRPolGEMFarSide_dighit_strip->at(j+1))>=128
	       floor(simev->Tgenrp->Harm_PRPolGEMFarSide_dighit_strip->at(j+1)/128)!=floor(lchan/128)
	       ){
	loadevt = true;
      }
	
      if(loadevt){
	if( crate >= 0 || slot >=  0 ) {
	  sldat = crateslot[idx(crate,slot)].get();
	}
	std::vector<UInt_t> *myev = &(map[sldat]);
	
	if(!samps.empty()){
	  //myev->push_back(SBSSimDataDecoder::EncodeHeader(5, apvnum, samps.size()));
	  //I think I'm onto something here, but I also need to transmit strip num 
	  myev->push_back(SBSSimDataDecoder::EncodeHeader(9, apvnum, samps.size()));
	  for(int k = 0; k<(int)samps.size(); k++){
	    // cout << " " << samps[k];
	    myev->push_back(strips[k]*8192+samps[k]);//strips[k]<< 13 | samps[k]);
	  }
	  //for(int l = 0; l<myev->size();l++)cout << myev->at(l) << " ";
	  //cout << endl;
	}
	//cout << endl;
	
	samps.clear();
	strips.clear();
      }
    }
  }
  
  
  
  /*
  while(j < HitData_Det->nhits){
    //Decode header first
    lchan = 0;
    if(HitData_Det->chan->at(j)<0){
      if(fDebug>2)
	std::cout << "j = " << j << " header = " << HitData_Det->dataword->at(j) << std::endl;
      SBSSimDataDecoder::DecodeHeader(HitData_Det->dataword->at(j),
				       data_type,chan_mult,nwords);
      
      //if header if from GEM detector, also decode the MPD header
      if(detname.find("gem")!=std::string::npos){
	for(uint k = 0; k<(HitData_Det->samps_datawords->at(j)).size();k++){
	  mpd_hdr[k] = (HitData_Det->samps_datawords->at(j)).at(k);
	}
	fEncoderMPD->DecodeMPDHeader(mpd_hdr, tmp_mpd);
	//reencode header for GEMs - not sure why - to set "chan" value ?
      }

      if(nwords>0)j++;
    }
    if(fDebug>2)
      std::cout << "j = " << j << " det chan = " << HitData_Det->chan->at(j) << std::endl;
    //channel should *not* be negative (unless there's a problem with nwords...)
    assert(HitData_Det->chan->at(j)>=0);
    //determine crate/slot
    lchan = (int)HitData_Det->chan->at(j);//+chan_mult*fNChan[detname];
    ChanToROC(detname, lchan, crate, slot, chan);

    if(fDebug>2)
      std::cout << "crate " << crate  << " slot " << slot << " chan " << chan << std::endl;
    if(detname.find("gem")!=std::string::npos){
      fEncoderMPD->EncodeMPDHeader(tmp_mpd, mpd_hdr, chan);
    }

    Decoder::THaSlotData *sldat = 0;
    if( crate >= 0 || slot >=  0 ) {
      sldat = crateslot[idx(crate,slot)].get();
    }
    
    //save the header
    std::vector<UInt_t> *myev = &(map[sldat]);
    myev->push_back(SBSSimDataDecoder::EncodeHeader(data_type,chan,nwords));
    if(detname.find("gem")!=std::string::npos){
      for(int k = 0; k<2;k++){ myev->push_back(mpd_hdr[k]);
      }
    }
    //Then save the hits
    //nwords = n following "hits" for ECal, Cher, Scint;
    //nowrds = n following hits*n data words for HCal, GEMs
    uint i = 0;
    while(i<nwords){
      if(fDebug>2)// || detname.find("grinch")!=std::string::npos 
	std::cout << " i = " << i << " j = " << j << " dataword = " << HitData_Det->dataword->at(j) << std::endl;
      if(detname.find("gem")!=std::string::npos || 
	 detname.find("hcal")!=std::string::npos){
	//if GEM or HCal, loop on "samples datawords" 
	if(HitData_Det->adc->at(j)>-9.e5){
	  // here dataword stores the number of samples datawords
	  for(int k = 0; k<HitData_Det->dataword->at(j);k++, i++){
	    myev->push_back( (HitData_Det->samps_datawords->at(j)).at(k) );
	    if(fDebug>2)
	      std::cout << " samp " << k << " dataword = " << (HitData_Det->samps_datawords->at(j)).at(k) << std::endl;
	  }
	}else{
	  //if adc has dummy value , it is a HCal TDC
	  myev->push_back(HitData_Det->dataword->at(j));
	}
      }else{
	//straightforward for detectors other than GEMs, HCal.
	myev->push_back(HitData_Det->dataword->at(j));
      }
      i++;
      j++;
    }
    if(fDebug>2){
      std::cout << " j = " << j << " my ev = {";
      for(size_t k = 0; k<myev->size(); k++)std::cout << myev->at(k) << " ; ";
      std::cout << " } " << std::endl;
    }
  }//end loop on j
  */
  return HED_OK;
}

/*
void SBSSimDecoder::SetDetMapParam(const std::string detname, int cps, int spc, int fs, int fc)
{
  fChansPerSlotDetMap[detname] = cps;
  fSlotsPerCrateDetMap[detname] = spc;
  fFirstSlotDetMap[detname] = fs;
  fFirstCrateDetMap[detname] = fc;
}
*/

void SBSSimDecoder::CheckForEnabledDetectors()
{
  //fDetectors = fManager->GetAllDetInfo();
  if(fDebug>0) {
    for(size_t i = 0; i < fDetectors.size(); i++) {
      std::cout << "Found detector: " << fDetectors[i].c_str() << endl;
      //<< ", ID: " << fDetectors[i].DetUniqueId() << std::endl;
    }
  }
  fCheckedForEnabledDetectors = true;
}

/*
void SBSSimDecoder::SetTree(TTree *t)
{
  if(t==0)return;
  fTree = new digsim_tree(t);
  if(fTree==0)return;
  fTreeIsSet = true;
}
*/

void SBSSimDecoder::SetDetectors()
{
  //std::cout << "[SBSSimDecoder::SetDetectors()]: rundate = ";

  //TDatime rundate = gHaRun->GetDate(); //will this work? answer appears to be NO
  
  //If the following works, we should be gold:
  TDatime rundate;
  rundate.Set( GetRunTime() ); //GetRunTime() gives the run time as a UNIX time

  rundate.Print();

  TIter aiter(gHaApps);
  THaApparatus* app = 0;
  while( (app=(THaApparatus*)aiter()) ){
    TList* listdet = app->GetDetectors();
    TIter diter(listdet);
    TObject* det = 0;
    while( (det=(TObject*)diter()) ){
      cout << "Setting det " << app->GetName() << "." << det->GetName() 
	   << " into SBSSimDecoder" << endl;
      if(strcmp(app->GetDetector(det->GetName())->GetClassName(),"SBSBBTotalShower")==0){
	SBSBBTotalShower* TS = (SBSBBTotalShower*)app->GetDetector(det->GetName());
	// AddDetector(Form("%s.%s",app->GetName(), TS->GetShower()->GetName()), 
	// 	    (app->GetDetector(det->GetName()))->GetInitDate());
	// AddDetector(Form("%s.%s",app->GetName(), TS->GetPreShower()->GetName()), 
	// 	    (app->GetDetector(det->GetName()))->GetInitDate());

	AddDetector(Form("%s.%s",app->GetName(), TS->GetShower()->GetName()), 
		    rundate);
	AddDetector(Form("%s.%s",app->GetName(), TS->GetPreShower()->GetName()), 
		    rundate);
       }else{
	// AddDetector(Form("%s.%s",app->GetName(), det->GetName()), 
	// 	    (app->GetDetector(det->GetName()))->GetInitDate());

	AddDetector(Form("%s.%s",app->GetName(), det->GetName()), 
		    rundate);
      }
    }
  }
}

Int_t SBSSimDecoder::AddDetector(std::string detname, TDatime date)
{
  fDetectors.push_back(detname);
  return ReadDetectorDB(detname, date);
}

Int_t SBSSimDecoder::ReadDetectorDB(std::string detname, TDatime date)
{
  //EPAF: in here the det name is the "full" det name i.e. including the spectro name
  std::string path = std::string(std::getenv("SBS"))+"/DB/";
  if(std::getenv("DB_DIR")) {
    path = std::string(std::getenv("DB_DIR"))+"/";
  }
  const string& fileName = path+"db_"+detname+".dat";
  
  const string prefix = detname+".";
  // First, open the common db file and parse info there, later, the
  // digitization specific db can be used to override any values
  //FILE* file  = Podd::OpenDBFile(fileName.c_str(), date);
  
  std::cout << "Calling ReadDetectorDB for detector " << detname << ", Date = ";
  date.Print();

  // FILE *file = Podd::OpenDBFile( detname.c_str(), date, "SBSSimDecoder::ReadDetectorDB()", 
  // 				 "r", 2 );

  FILE *file = Podd::OpenDBFile( detname.c_str(), date );
  
  if( !file ) return THaAnalysisObject::kFileError;

  std::vector<int> detmap,chanmap;//, detmap_adc;
  uint nchan, nlogchan = 0, chanmapstart = 0;
  
  //int cps, spc, fs, fc;
  
  bool isgem = (detname.find("gem")!=std::string::npos);
  int apv_num = -1, mpd = -1, mod = 0, axis = -1;
  //int pos = -1;

  DBRequest request[] = {
    {"nchan", &nchan, kInt, 0, false},// 
    {"nlog_chan", &nlogchan, kInt, 0, true},// <- optional
    {"detmap", &detmap, kIntV, 0, false}, //
    {"chanmap", &chanmap, kIntV, 0, true}, // <- optional
    {"chanmap_start", &chanmapstart, kInt, 0, true}, // <- optional
    //{"detmap_adc", &detmap_adc, kIntV, 0, true}, // <- optional
    /*
    {"first_crate", &fc, kInt, 0, true},// <- optional 
    {"first_slot", &fs, kInt, 0, true},//  <- optional
    {"chan_per_slot", &cps, kInt, 0, true},//  <- optional
    {"slot_per_crate", &spc, kInt, 0, true},//  <- optional
    */
    { 0 }
  };
  Int_t err;
  int nparam_mod = 5;
  if(isgem){//gem detectors
    nparam_mod = 9;
  }
  int crate,slot,ch_lo,ch_hi, ch_ref, ch_count = 0, ch_map = 0;
  
  if(isgem){//it's easier if gems are their own thing
    /*
    std::string chambers;
    DBRequest req_chambers[] = {
      {"chambers", &chambers, kString, 0, false}, //
      { 0 }
    };
    err = THaAnalysisObject::LoadDB(file, date, req_chambers, prefix.c_str());
    
    //cout << " prefix " << prefix.c_str() << " err " << err << " chambers " << chambers.c_str() << " size ? " << chambers.size() << endl;
    
    std::vector<std::string> chambers_names;
    if(err==0)chambers_names = vsplit(chambers);
    
    if(!chambers_names.empty()){
      for (std::vector<std::string>::iterator it = chambers_names.begin() ; it != chambers_names.end(); ++it){
    */
    std::string modules;
    //std::string pref_cham = prefix+(*it)+".";
    //std::string pref_cham = prefix;//+".";
    //cout << "prefix chamber "  << pref_cham.c_str() << endl;
    DBRequest req_modules[] = {
      {"modules", &modules, kString, 0, false}, //
      { 0 }
    };
    err = THaAnalysisObject::LoadDB(file, date, req_modules, prefix.c_str());

    if(err)return THaAnalysisObject::kInitError;
    
    //cout << " prefix " << pref_cham.c_str() << " err " << err << " modules " << modules.c_str() << " size ? " << modules.size() << endl;
    
    std::vector<std::string> modules_names;
    if(err==0)modules_names = vsplit(modules);
    if(!modules_names.empty()){
      for (std::vector<std::string>::iterator jt = modules_names.begin() ; jt != modules_names.end(); ++jt){
	std::string pref_mod = prefix+(*jt)+".";
	
	DBRequest request_gem[] = {
	  {"chanmap", &chanmap, kIntV, 0, false}, 
	  { 0 }
	};
	err+= THaAnalysisObject::LoadDB(file, date, request_gem, pref_mod.c_str());

	fInvGEMDetMap[detname].resize(fInvGEMDetMap[detname].size()+2); //increments the size of this container by two. But it never gets initialized prior to now! We have to trust that the size is zero to start with. Is this safe? Probably not... 
	
	//This resizes the fInvGEMDetMap[detname][2*module+axis] to the total size of the decode map
	for(int m = 0; m<2; m++)(fInvGEMDetMap[detname])[mod*2+m].resize(chanmap.size()/nparam_mod);
	
	// std::cout << "(detname, mod, nparam_mod)=(" << detname << ", " << mod << ", " << nparam_mod 
	// 	  << ")" << std::endl;

	//	int nparam_mod = 9;
	int ax_prev = 0;
	int n_ax = 0, n_ax_x = 0, n_ax_y = 0;
	for(size_t k = 0; k < chanmap.size(); k+=nparam_mod) {
	  //for(int m = 0; m<nparam_mod; m++)std::cout << chanmap[k+m] << " ";
	  //std::cout << std::endl;
	  crate  = chanmap[k];
	  slot   = chanmap[k+1];
	  mpd   = chanmap[k+2];
	  apv_num = mpd << 4 | chanmap[k+4];//
	  //pos = chanmap[k+6];
	  axis = chanmap[k+8];
	  if(axis==0)n_ax_x++;
	  if(axis==1)n_ax_y++;
	  if(ax_prev!=axis){
	    n_ax = 0;
	    ax_prev = axis;
	  }
	  ch_lo = 128*n_ax;
	  ch_hi = 128*(n_ax+1)-1;
	  //mod*2+axis???
	  //std::cout << mod << " " << mod*2+axis << " " << fInvGEMDetMap[detname].size() << " " << mpd << " " << chanmap[k+4] << " " << apv_num << " " << n_ax << endl;
	  (fInvGEMDetMap[detname])[mod*2+axis][n_ax]=gemstripinfo(crate, slot, apv_num);
	  n_ax++;
	}
	(fInvGEMDetMap[detname])[mod*2+0].resize(n_ax_x);
	(fInvGEMDetMap[detname])[mod*2+1].resize(n_ax_y);
	/*
	std::string planeconfig;
	//cout << "prefix module "  << pref_mod.c_str() << endl;
	DBRequest req_planeconfig[] = {
	  {"planeconfig", &planeconfig, kString, 0, false}, //
	  { 0 }
	};
	err+= THaAnalysisObject::LoadDB(file, date, req_planeconfig, pref_mod.c_str());
	//cout << " prefix " << pref_mod.c_str() << " err " << err << " planeconfig " << planeconfig.c_str() << " size ? " << planeconfig.size() << endl;
	std::vector<std::string> plane_readouts;
	if(err==0)plane_readouts = vsplit(planeconfig);
	if(!plane_readouts.empty()){
	  for (std::vector<std::string>::iterator kt = plane_readouts.begin() ; kt != plane_readouts.end(); ++kt){
	    //mod++;
	    std::string pref_ro = pref_mod+(*kt)+".";
	    ch_count = 0;
	    fInvGEMDetMap[detname].resize(fInvGEMDetMap[detname].size()+1);
	    
	    if(fDebug>=2)cout << fInvGEMDetMap[detname].size() << " module number " << mod << endl;
	    
	    err+= THaAnalysisObject::LoadDB(file, date, request, pref_ro.c_str());
	    //if(nlogchan==0)
	    nlogchan = nchan;
	    if(fDebug>=2)cout << " prefix " << pref_ro.c_str() << " err " << err << endl;
	    
	    if(err==0)fInvGEMDetMap[detname][mod].resize(nchan);
	    
	    //int nparam_mod = 5;
	    for(size_t k = 0; k < detmap.size(); k+=nparam_mod) {
	      crate  = detmap[k];
	      slot   = detmap[k+1];
	      ch_lo  = detmap[k+2];
	      ch_hi  = detmap[k+3];
	      
	      for(int i = ch_lo; i<=ch_hi; i++, ch_count++){
		if(i%128==0){
		  apv_num++;
		  if(fDebug>=3)cout << crate << " " << slot << " " << i << " " << apv_num << endl;
		}
		if(ch_count>nlogchan){
		  std::cout << " <0> number of channels defined in detmap ( >= " << ch_count << ") exceeds logical number of channels = " << nlogchan << std::endl;
		  return THaAnalysisObject::kInitError;
		}
		(fInvGEMDetMap[detname])[mod][ch_count]=gemstripinfo(crate, slot, i, apv_num);
	      }
	      
	    }
	  }//end loop on kt
	}//end if !plane_readouts
	*/
	mod++;
      }//end loop on jt
    }//end if !modules_names
    //  }//end loop on it
    //}//end if !chambers_names
    
  }else{
   err = THaAnalysisObject::LoadDB(file, date, request, prefix.c_str());
   //}
   // Could close the common file already
   fclose(file);
   if(nlogchan==0)nlogchan = nchan;
   
   
   if(err)return THaAnalysisObject::kInitError;
  
   //fNChanDet[detname] = nchan;
   //fChanMapStartDet[detname] = chanmapstart;
   (fInvDetMap[detname]).resize(nlogchan+1);//for ref
   //if(detmap[4]==-1)nparam_mod = 5;
   for(size_t k = 0; k < detmap.size(); k+=nparam_mod) {
     crate  = detmap[k];
     slot   = detmap[k+1];
     ch_lo  = detmap[k+2];
     ch_hi  = detmap[k+3];
     ch_ref = detmap[k+4];
     if(ch_ref==-1){
       (fInvDetMap[detname])[nlogchan]=detchaninfo(crate, slot, ch_lo);
       continue;
     }
     /*
       if(detname.find("hodo")!=std::string::npos)
       cout << " crate " << crate << " slot " << slot 
       << " ch_lo " << ch_lo << " ch_hi " << ch_hi << endl;
     */
     if(chanmap.empty()){
       for(int i = ch_lo; i<=ch_hi; i++, ch_count++){
	 /*
	 if(isgem && i%128==0){
	   apv_num++;
	   cout << crate << " " << slot << " " << i << " " << apv_num << endl;
	 }
	 */
	 if(ch_count>(int)nlogchan){
	   std::cout << " <1> number of channels defined in detmap ( >= " << ch_count << ") exceeds logical number of channels = " << nlogchan << std::endl;
	   return THaAnalysisObject::kInitError;
	 }
	 
	 (fInvDetMap[detname])[ch_count]=detchaninfo(crate, slot, i);
	 //cout << "ch_count " << ch_count << " crate " << crate << " slot " << slot << " i " << i << " &(fInvDetMap[detname]).at(ch_count) " << &(fInvDetMap[detname]).at(ch_count) << endl ;
	 /*
	   if(detname.find("hodo")!=std::string::npos){
	   cout << " crate " << crate << " slot " << slot 
	   << " i " << i << " ch_count " << ch_count << endl;
	   cout << &(fInvDetMap.at(detname)).at(ch_count) << endl;
	   }
	 */
       }
     }else{
       int chan_offset = 1;
       if(detname.find("sh")!=std::string::npos)chan_offset = 0;
       if(detname.find("hodo")!=std::string::npos)chan_offset = 0;
       if(detname.find("grinch")!=std::string::npos)chan_offset = 0;
       if(detname.find("ps")!=std::string::npos)chan_offset = 0;
       if(detname.find("active_ana")!=std::string::npos)chan_offset = 0;
       if(detname.find("prpolscint_farside")!=std::string::npos)chan_offset = 0;
       for(int i = ch_lo; i<=ch_hi; i++, ch_map++){
	 if(ch_count>(int)nlogchan){
	   std::cout << " <2> number of channels defined in detmap ( >= " << ch_count << ") exceeds logical number of channels = " << nlogchan << std::endl;
	   return THaAnalysisObject::kInitError;
	 }
	 if(fDebug>=2)std::cout << " i = " << i << ", crate = " << crate << ", slot = " << slot <<  ", ch_count = " << ch_count << " chan = " << chanmap[ch_map]-chan_offset << " (+" << nchan << ") " << std::endl;
	 if(chanmap[ch_map]>=0){
	   if(ch_count<(int)nchan){
	     (fInvDetMap[detname])[chanmap[ch_map]-chan_offset]=detchaninfo(crate, slot, i);
	     if(fDebug>=3)std::cout << chanmap[ch_map]-chan_offset << " " << &(fInvDetMap.at(detname)).at(chanmap[ch_map]-chan_offset) << std::endl;
	   }else{
	     (fInvDetMap[detname])[chanmap[ch_map]+nchan-chan_offset]=detchaninfo(crate, slot, i);
	     if(fDebug>=3)std::cout <<&(fInvDetMap.at(detname)).at(chanmap[ch_map]+nchan-chan_offset) << std::endl;
	   }
	   ch_count++;
	 }
       }
     }
   }
  }//end else (if isgem)
  /*
  fChansPerSlotDetMap[detname] = cps;
  fSlotsPerCrateDetMap[detname] = spc;
  fFirstSlotDetMap[detname] = fs;
  fFirstCrateDetMap[detname] = fc;
  */
  
  return(THaAnalysisObject::kOK);
}


//-----------------------------------------------------------------------------
//static inline
void SBSSimDecoder::ChanToROC(const std::string& detname, Int_t h_chan,
			       Int_t& crate, Int_t& slot, UShort_t& chan )const 
{
  // Convert location parameters (row, col, chan) of the given Channel
  // to hardware channel (crate,slot,chan)
  // The (crate,slot,chan) assignment must match the detmap definition in
  // the database!  See TreeSearch/dbconvert.cxx
  // In the case of GRINCH/RICH: 
  // crate = GTP; slot = VETROC; chan = PMT. (NINOs are "transparent", in a similar way to the MPDs)
  
  /*
  int CPS = fChansPerSlotDetMap.at(detname);
  int SPC = fSlotsPerCrateDetMap.at(detname);
  int FS = fFirstSlotDetMap.at(detname);
  int FC = fFirstCrateDetMap.at(detname);
  
  //div_t d = div( h_chan, fManager->GetChanPerSlot() );
  div_t d = div( h_chan, CPS );
  slot = d.quot;
  chan = d.rem;

  d = div( slot, SPC );
  crate = d.quot+FC;
  slot  = d.rem+FS;
  */
  if(h_chan>=fInvDetMap.at(detname).size())std::cout << " " << detname << " "  << h_chan << " " << &fInvDetMap.at(detname) << std::endl;
  assert(h_chan<fInvDetMap.at(detname).size());
  
  if(fDebug>3){
  
    std::cout << &(fInvDetMap.at(detname)).at(h_chan) << std::endl;
  }
  crate = ((fInvDetMap.at(detname)).at(h_chan)).crate;
  slot = ((fInvDetMap.at(detname)).at(h_chan)).slot;
  chan = ((fInvDetMap.at(detname)).at(h_chan)).chan;
  
}

int SBSSimDecoder::APVnum(const std::string& detname, Int_t mod, Int_t h_chan,
			  Int_t &crate, Int_t &slot, UShort_t &chan) const
{
  chan = h_chan%128;
  int n = (h_chan-chan)/128;

  // std::cout << "(detname, mod, h_chan, chan, n )= (" << detname << ", " << mod << ", "
  // 	    << h_chan << ", " << chan << ", " << n << ")" << std::endl;
  
  assert(mod<fInvGEMDetMap.at(detname).size());
  assert(n<(fInvGEMDetMap.at(detname)[mod]).size());

  // if( mod>fInvGEMDetMap.at(detname).size() ){
  //   std::err << "ERROR: map size =  " << " for detector " << detname 
  // 	     << " is smaller than module size =  " << mod << std:: endl;
  //   //return -1;
  // }else if( n>(fInvGEMDetMap.at(detname)[mod]).size() ){
  //   return -1;
  // }
  //if((fInvGEMDetMap.at(detname))[mod][n].chan_lo<=h_chan &&
  // hchan <= (fInvGEMDetMap.at(detname))[mod][n].chan_hi){
  crate = ((fInvGEMDetMap.at(detname))[mod][n]).crate;
  slot = ((fInvGEMDetMap.at(detname))[mod][n]).slot;
  return ((fInvGEMDetMap.at(detname))[mod][n]).apvnum;
  //}else{
  //return -1;
  //}
}

/*
//-----------------------------------------------------------------------------
static inline
Int_t MakeROCKey( Int_t crate, Int_t slot, Int_t chan )
{
  return chan;// +
  //fManager->GetChanPerSlot()*( slot + fManager->GetSlotPerCrate()*crate );
}

//-----------------------------------------------------------------------------
Int_t SBSSimDecoder::ChanFromROC( Int_t crate, Int_t slot, Int_t chan ) const
{
  // Return index of digitized strip correspomding to hardware channel
  // (crate,slot,chan)

  if( fPMTMap.empty() )
    return -1;

  PMTMap_t::const_iterator found = fPMTMap.find( MakeROCKey(crate,slot,chan) );
  if( found == fPMTMap.end() )
    return -1;

  return found->second;
}
*/
