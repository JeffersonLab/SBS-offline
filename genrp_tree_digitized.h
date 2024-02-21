//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Feb 13 16:53:48 2024 by ROOT version 6.26/10
// from TTree T/Geant4 SBS Simulation
// found on file: genrp_elastic_gc_job_0.root
//////////////////////////////////////////////////////////

#ifndef genrp_tree_digitized_h
#define genrp_tree_digitized_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"

class genrp_tree_digitized {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Double_t        TargPol;
   Double_t        TargThetaSpin;
   Double_t        TargPhiSpin;
   Double_t        BeamPol;
   Double_t        BeamThetaSpin;
   Double_t        BeamPhiSpin;
   Double_t        Harm_ActAnScint_det_esum;
   Int_t           Harm_ActAnScint_hit_nhits;
   std::vector<int>     *Harm_ActAnScint_hit_row;
   std::vector<int>     *Harm_ActAnScint_hit_col;
   std::vector<int>     *Harm_ActAnScint_hit_cell;
   std::vector<int>     *Harm_ActAnScint_hit_plane;
   std::vector<int>     *Harm_ActAnScint_hit_wire;
   std::vector<double>  *Harm_ActAnScint_hit_xcell;
   std::vector<double>  *Harm_ActAnScint_hit_ycell;
   std::vector<double>  *Harm_ActAnScint_hit_zcell;
   std::vector<double>  *Harm_ActAnScint_hit_xcellg;
   std::vector<double>  *Harm_ActAnScint_hit_ycellg;
   std::vector<double>  *Harm_ActAnScint_hit_zcellg;
   std::vector<double>  *Harm_ActAnScint_hit_xhit;
   std::vector<double>  *Harm_ActAnScint_hit_yhit;
   std::vector<double>  *Harm_ActAnScint_hit_zhit;
   std::vector<double>  *Harm_ActAnScint_hit_xhitg;
   std::vector<double>  *Harm_ActAnScint_hit_yhitg;
   std::vector<double>  *Harm_ActAnScint_hit_zhitg;
   std::vector<double>  *Harm_ActAnScint_hit_sumedep;
   std::vector<double>  *Harm_ActAnScint_hit_tavg;
   std::vector<double>  *Harm_ActAnScint_hit_trms;
   std::vector<double>  *Harm_ActAnScint_hit_tmin;
   std::vector<double>  *Harm_ActAnScint_hit_tmax;
   std::vector<int>     *Harm_ActAnScint_hit_otridx;
   std::vector<int>     *Harm_ActAnScint_hit_ptridx;
   std::vector<int>     *Harm_ActAnScint_hit_sdtridx;
   Int_t           Harm_CEPolFront_hit_nhits;
   std::vector<int>     *Harm_CEPolFront_hit_plane;
   std::vector<int>     *Harm_CEPolFront_hit_strip;
   std::vector<double>  *Harm_CEPolFront_hit_x;
   std::vector<double>  *Harm_CEPolFront_hit_y;
   std::vector<double>  *Harm_CEPolFront_hit_z;
   std::vector<double>  *Harm_CEPolFront_hit_polx;
   std::vector<double>  *Harm_CEPolFront_hit_poly;
   std::vector<double>  *Harm_CEPolFront_hit_polz;
   std::vector<double>  *Harm_CEPolFront_hit_t;
   std::vector<double>  *Harm_CEPolFront_hit_trms;
   std::vector<double>  *Harm_CEPolFront_hit_tmin;
   std::vector<double>  *Harm_CEPolFront_hit_tmax;
   std::vector<double>  *Harm_CEPolFront_hit_tx;
   std::vector<double>  *Harm_CEPolFront_hit_ty;
   std::vector<double>  *Harm_CEPolFront_hit_xin;
   std::vector<double>  *Harm_CEPolFront_hit_yin;
   std::vector<double>  *Harm_CEPolFront_hit_zin;
   std::vector<double>  *Harm_CEPolFront_hit_xout;
   std::vector<double>  *Harm_CEPolFront_hit_yout;
   std::vector<double>  *Harm_CEPolFront_hit_zout;
   std::vector<double>  *Harm_CEPolFront_hit_txp;
   std::vector<double>  *Harm_CEPolFront_hit_typ;
   std::vector<double>  *Harm_CEPolFront_hit_xg;
   std::vector<double>  *Harm_CEPolFront_hit_yg;
   std::vector<double>  *Harm_CEPolFront_hit_zg;
   std::vector<int>     *Harm_CEPolFront_hit_trid;
   std::vector<int>     *Harm_CEPolFront_hit_mid;
   std::vector<int>     *Harm_CEPolFront_hit_pid;
   std::vector<double>  *Harm_CEPolFront_hit_vx;
   std::vector<double>  *Harm_CEPolFront_hit_vy;
   std::vector<double>  *Harm_CEPolFront_hit_vz;
   std::vector<double>  *Harm_CEPolFront_hit_p;
   std::vector<double>  *Harm_CEPolFront_hit_edep;
   std::vector<double>  *Harm_CEPolFront_hit_beta;
   std::vector<int>     *Harm_CEPolFront_hit_otridx;
   std::vector<int>     *Harm_CEPolFront_hit_ptridx;
   std::vector<int>     *Harm_CEPolFront_hit_sdtridx;
   Int_t           Harm_CEPolFront_Track_ntracks;
   std::vector<int>     *Harm_CEPolFront_Track_TID;
   std::vector<int>     *Harm_CEPolFront_Track_PID;
   std::vector<int>     *Harm_CEPolFront_Track_MID;
   std::vector<int>     *Harm_CEPolFront_Track_NumHits;
   std::vector<int>     *Harm_CEPolFront_Track_NumPlanes;
   std::vector<int>     *Harm_CEPolFront_Track_NDF;
   std::vector<double>  *Harm_CEPolFront_Track_Chi2fit;
   std::vector<double>  *Harm_CEPolFront_Track_Chi2true;
   std::vector<double>  *Harm_CEPolFront_Track_X;
   std::vector<double>  *Harm_CEPolFront_Track_Y;
   std::vector<double>  *Harm_CEPolFront_Track_Xp;
   std::vector<double>  *Harm_CEPolFront_Track_Yp;
   std::vector<double>  *Harm_CEPolFront_Track_T;
   std::vector<double>  *Harm_CEPolFront_Track_P;
   std::vector<double>  *Harm_CEPolFront_Track_Sx;
   std::vector<double>  *Harm_CEPolFront_Track_Sy;
   std::vector<double>  *Harm_CEPolFront_Track_Sz;
   std::vector<double>  *Harm_CEPolFront_Track_Xfit;
   std::vector<double>  *Harm_CEPolFront_Track_Yfit;
   std::vector<double>  *Harm_CEPolFront_Track_Xpfit;
   std::vector<double>  *Harm_CEPolFront_Track_Ypfit;
   std::vector<int>     *Harm_CEPolFront_Track_otridx;
   std::vector<int>     *Harm_CEPolFront_Track_ptridx;
   std::vector<int>     *Harm_CEPolFront_Track_sdtridx;
   Int_t           Harm_CEPolRear_hit_nhits;
   std::vector<int>     *Harm_CEPolRear_hit_plane;
   std::vector<int>     *Harm_CEPolRear_hit_strip;
   std::vector<double>  *Harm_CEPolRear_hit_x;
   std::vector<double>  *Harm_CEPolRear_hit_y;
   std::vector<double>  *Harm_CEPolRear_hit_z;
   std::vector<double>  *Harm_CEPolRear_hit_polx;
   std::vector<double>  *Harm_CEPolRear_hit_poly;
   std::vector<double>  *Harm_CEPolRear_hit_polz;
   std::vector<double>  *Harm_CEPolRear_hit_t;
   std::vector<double>  *Harm_CEPolRear_hit_trms;
   std::vector<double>  *Harm_CEPolRear_hit_tmin;
   std::vector<double>  *Harm_CEPolRear_hit_tmax;
   std::vector<double>  *Harm_CEPolRear_hit_tx;
   std::vector<double>  *Harm_CEPolRear_hit_ty;
   std::vector<double>  *Harm_CEPolRear_hit_xin;
   std::vector<double>  *Harm_CEPolRear_hit_yin;
   std::vector<double>  *Harm_CEPolRear_hit_zin;
   std::vector<double>  *Harm_CEPolRear_hit_xout;
   std::vector<double>  *Harm_CEPolRear_hit_yout;
   std::vector<double>  *Harm_CEPolRear_hit_zout;
   std::vector<double>  *Harm_CEPolRear_hit_txp;
   std::vector<double>  *Harm_CEPolRear_hit_typ;
   std::vector<double>  *Harm_CEPolRear_hit_xg;
   std::vector<double>  *Harm_CEPolRear_hit_yg;
   std::vector<double>  *Harm_CEPolRear_hit_zg;
   std::vector<int>     *Harm_CEPolRear_hit_trid;
   std::vector<int>     *Harm_CEPolRear_hit_mid;
   std::vector<int>     *Harm_CEPolRear_hit_pid;
   std::vector<double>  *Harm_CEPolRear_hit_vx;
   std::vector<double>  *Harm_CEPolRear_hit_vy;
   std::vector<double>  *Harm_CEPolRear_hit_vz;
   std::vector<double>  *Harm_CEPolRear_hit_p;
   std::vector<double>  *Harm_CEPolRear_hit_edep;
   std::vector<double>  *Harm_CEPolRear_hit_beta;
   std::vector<int>     *Harm_CEPolRear_hit_otridx;
   std::vector<int>     *Harm_CEPolRear_hit_ptridx;
   std::vector<int>     *Harm_CEPolRear_hit_sdtridx;
   Int_t           Harm_CEPolRear_Track_ntracks;
   std::vector<int>     *Harm_CEPolRear_Track_TID;
   std::vector<int>     *Harm_CEPolRear_Track_PID;
   std::vector<int>     *Harm_CEPolRear_Track_MID;
   std::vector<int>     *Harm_CEPolRear_Track_NumHits;
   std::vector<int>     *Harm_CEPolRear_Track_NumPlanes;
   std::vector<int>     *Harm_CEPolRear_Track_NDF;
   std::vector<double>  *Harm_CEPolRear_Track_Chi2fit;
   std::vector<double>  *Harm_CEPolRear_Track_Chi2true;
   std::vector<double>  *Harm_CEPolRear_Track_X;
   std::vector<double>  *Harm_CEPolRear_Track_Y;
   std::vector<double>  *Harm_CEPolRear_Track_Xp;
   std::vector<double>  *Harm_CEPolRear_Track_Yp;
   std::vector<double>  *Harm_CEPolRear_Track_T;
   std::vector<double>  *Harm_CEPolRear_Track_P;
   std::vector<double>  *Harm_CEPolRear_Track_Sx;
   std::vector<double>  *Harm_CEPolRear_Track_Sy;
   std::vector<double>  *Harm_CEPolRear_Track_Sz;
   std::vector<double>  *Harm_CEPolRear_Track_Xfit;
   std::vector<double>  *Harm_CEPolRear_Track_Yfit;
   std::vector<double>  *Harm_CEPolRear_Track_Xpfit;
   std::vector<double>  *Harm_CEPolRear_Track_Ypfit;
   std::vector<int>     *Harm_CEPolRear_Track_otridx;
   std::vector<int>     *Harm_CEPolRear_Track_ptridx;
   std::vector<int>     *Harm_CEPolRear_Track_sdtridx;
   Int_t           Harm_PRPolGEMFarSide_hit_nhits;
   std::vector<int>     *Harm_PRPolGEMFarSide_hit_plane;
   std::vector<int>     *Harm_PRPolGEMFarSide_hit_strip;
   std::vector<double>  *Harm_PRPolGEMFarSide_hit_x;
   std::vector<double>  *Harm_PRPolGEMFarSide_hit_y;
   std::vector<double>  *Harm_PRPolGEMFarSide_hit_z;
   std::vector<double>  *Harm_PRPolGEMFarSide_hit_polx;
   std::vector<double>  *Harm_PRPolGEMFarSide_hit_poly;
   std::vector<double>  *Harm_PRPolGEMFarSide_hit_polz;
   std::vector<double>  *Harm_PRPolGEMFarSide_hit_t;
   std::vector<double>  *Harm_PRPolGEMFarSide_hit_trms;
   std::vector<double>  *Harm_PRPolGEMFarSide_hit_tmin;
   std::vector<double>  *Harm_PRPolGEMFarSide_hit_tmax;
   std::vector<double>  *Harm_PRPolGEMFarSide_hit_tx;
   std::vector<double>  *Harm_PRPolGEMFarSide_hit_ty;
   std::vector<double>  *Harm_PRPolGEMFarSide_hit_xin;
   std::vector<double>  *Harm_PRPolGEMFarSide_hit_yin;
   std::vector<double>  *Harm_PRPolGEMFarSide_hit_zin;
   std::vector<double>  *Harm_PRPolGEMFarSide_hit_xout;
   std::vector<double>  *Harm_PRPolGEMFarSide_hit_yout;
   std::vector<double>  *Harm_PRPolGEMFarSide_hit_zout;
   std::vector<double>  *Harm_PRPolGEMFarSide_hit_txp;
   std::vector<double>  *Harm_PRPolGEMFarSide_hit_typ;
   std::vector<double>  *Harm_PRPolGEMFarSide_hit_xg;
   std::vector<double>  *Harm_PRPolGEMFarSide_hit_yg;
   std::vector<double>  *Harm_PRPolGEMFarSide_hit_zg;
   std::vector<int>     *Harm_PRPolGEMFarSide_hit_trid;
   std::vector<int>     *Harm_PRPolGEMFarSide_hit_mid;
   std::vector<int>     *Harm_PRPolGEMFarSide_hit_pid;
   std::vector<double>  *Harm_PRPolGEMFarSide_hit_vx;
   std::vector<double>  *Harm_PRPolGEMFarSide_hit_vy;
   std::vector<double>  *Harm_PRPolGEMFarSide_hit_vz;
   std::vector<double>  *Harm_PRPolGEMFarSide_hit_p;
   std::vector<double>  *Harm_PRPolGEMFarSide_hit_edep;
   std::vector<double>  *Harm_PRPolGEMFarSide_hit_beta;
   std::vector<int>     *Harm_PRPolGEMFarSide_hit_otridx;
   std::vector<int>     *Harm_PRPolGEMFarSide_hit_ptridx;
   std::vector<int>     *Harm_PRPolGEMFarSide_hit_sdtridx;
   Int_t           Harm_PRPolGEMFarSide_Track_ntracks;
   std::vector<int>     *Harm_PRPolGEMFarSide_Track_TID;
   std::vector<int>     *Harm_PRPolGEMFarSide_Track_PID;
   std::vector<int>     *Harm_PRPolGEMFarSide_Track_MID;
   std::vector<int>     *Harm_PRPolGEMFarSide_Track_NumHits;
   std::vector<int>     *Harm_PRPolGEMFarSide_Track_NumPlanes;
   std::vector<int>     *Harm_PRPolGEMFarSide_Track_NDF;
   std::vector<double>  *Harm_PRPolGEMFarSide_Track_Chi2fit;
   std::vector<double>  *Harm_PRPolGEMFarSide_Track_Chi2true;
   std::vector<double>  *Harm_PRPolGEMFarSide_Track_X;
   std::vector<double>  *Harm_PRPolGEMFarSide_Track_Y;
   std::vector<double>  *Harm_PRPolGEMFarSide_Track_Xp;
   std::vector<double>  *Harm_PRPolGEMFarSide_Track_Yp;
   std::vector<double>  *Harm_PRPolGEMFarSide_Track_T;
   std::vector<double>  *Harm_PRPolGEMFarSide_Track_P;
   std::vector<double>  *Harm_PRPolGEMFarSide_Track_Sx;
   std::vector<double>  *Harm_PRPolGEMFarSide_Track_Sy;
   std::vector<double>  *Harm_PRPolGEMFarSide_Track_Sz;
   std::vector<double>  *Harm_PRPolGEMFarSide_Track_Xfit;
   std::vector<double>  *Harm_PRPolGEMFarSide_Track_Yfit;
   std::vector<double>  *Harm_PRPolGEMFarSide_Track_Xpfit;
   std::vector<double>  *Harm_PRPolGEMFarSide_Track_Ypfit;
   std::vector<int>     *Harm_PRPolGEMFarSide_Track_otridx;
   std::vector<int>     *Harm_PRPolGEMFarSide_Track_ptridx;
   std::vector<int>     *Harm_PRPolGEMFarSide_Track_sdtridx;
   Double_t        Harm_PRPolScintFarSide_det_esum;
   Int_t           Harm_PRPolScintFarSide_hit_nhits;
   std::vector<int>     *Harm_PRPolScintFarSide_hit_row;
   std::vector<int>     *Harm_PRPolScintFarSide_hit_col;
   std::vector<int>     *Harm_PRPolScintFarSide_hit_cell;
   std::vector<int>     *Harm_PRPolScintFarSide_hit_plane;
   std::vector<int>     *Harm_PRPolScintFarSide_hit_wire;
   std::vector<double>  *Harm_PRPolScintFarSide_hit_xcell;
   std::vector<double>  *Harm_PRPolScintFarSide_hit_ycell;
   std::vector<double>  *Harm_PRPolScintFarSide_hit_zcell;
   std::vector<double>  *Harm_PRPolScintFarSide_hit_xcellg;
   std::vector<double>  *Harm_PRPolScintFarSide_hit_ycellg;
   std::vector<double>  *Harm_PRPolScintFarSide_hit_zcellg;
   std::vector<double>  *Harm_PRPolScintFarSide_hit_xhit;
   std::vector<double>  *Harm_PRPolScintFarSide_hit_yhit;
   std::vector<double>  *Harm_PRPolScintFarSide_hit_zhit;
   std::vector<double>  *Harm_PRPolScintFarSide_hit_xhitg;
   std::vector<double>  *Harm_PRPolScintFarSide_hit_yhitg;
   std::vector<double>  *Harm_PRPolScintFarSide_hit_zhitg;
   std::vector<double>  *Harm_PRPolScintFarSide_hit_sumedep;
   std::vector<double>  *Harm_PRPolScintFarSide_hit_tavg;
   std::vector<double>  *Harm_PRPolScintFarSide_hit_trms;
   std::vector<double>  *Harm_PRPolScintFarSide_hit_tmin;
   std::vector<double>  *Harm_PRPolScintFarSide_hit_tmax;
   std::vector<int>     *Harm_PRPolScintFarSide_hit_otridx;
   std::vector<int>     *Harm_PRPolScintFarSide_hit_ptridx;
   std::vector<int>     *Harm_PRPolScintFarSide_hit_sdtridx;
   Int_t           Harm_ActAn_dighit_nchan;
   std::vector<int>     *Harm_ActAn_dighit_chan;
   std::vector<int>     *Harm_ActAn_dighit_adc;
   std::vector<int>     *Harm_ActAn_dighit_samp;
   std::vector<int>     *Harm_ActAn_dighit_tdc;
   Int_t           Harm_PRPolScintFarSide_dighit_nchan;
   std::vector<int>     *Harm_PRPolScintFarSide_dighit_chan;
   std::vector<int>     *Harm_PRPolScintFarSide_dighit_adc;
   std::vector<int>     *Harm_PRPolScintFarSide_dighit_tdc_l;
   std::vector<int>     *Harm_PRPolScintFarSide_dighit_tdc_t;
   Int_t           Harm_CEPolFront_dighit_nstrips;
   std::vector<int>     *Harm_CEPolFront_dighit_module;
   std::vector<int>     *Harm_CEPolFront_dighit_strip;
   std::vector<int>     *Harm_CEPolFront_dighit_adc;
   std::vector<int>     *Harm_CEPolFront_dighit_samp;
   Int_t           Harm_CEPolRear_dighit_nstrips;
   std::vector<int>     *Harm_CEPolRear_dighit_module;
   std::vector<int>     *Harm_CEPolRear_dighit_strip;
   std::vector<int>     *Harm_CEPolRear_dighit_adc;
   std::vector<int>     *Harm_CEPolRear_dighit_samp;
   Int_t           Harm_PRPolGEMFarSide_dighit_nstrips;
   std::vector<int>     *Harm_PRPolGEMFarSide_dighit_module;
   std::vector<int>     *Harm_PRPolGEMFarSide_dighit_strip;
   std::vector<int>     *Harm_PRPolGEMFarSide_dighit_adc;
   std::vector<int>     *Harm_PRPolGEMFarSide_dighit_samp;

   // List of branches
   TBranch        *b_TargPol;   //!
   TBranch        *b_TargThetaSpin;   //!
   TBranch        *b_TargPhiSpin;   //!
   TBranch        *b_BeamPol;   //!
   TBranch        *b_BeamThetaSpin;   //!
   TBranch        *b_BeamPhiSpin;   //!
   TBranch        *b_Harm_ActAnScint_det_esum;   //!
   TBranch        *b_Harm_ActAnScint_hit_nhits;   //!
   TBranch        *b_Harm_ActAnScint_hit_row;   //!
   TBranch        *b_Harm_ActAnScint_hit_col;   //!
   TBranch        *b_Harm_ActAnScint_hit_cell;   //!
   TBranch        *b_Harm_ActAnScint_hit_plane;   //!
   TBranch        *b_Harm_ActAnScint_hit_wire;   //!
   TBranch        *b_Harm_ActAnScint_hit_xcell;   //!
   TBranch        *b_Harm_ActAnScint_hit_ycell;   //!
   TBranch        *b_Harm_ActAnScint_hit_zcell;   //!
   TBranch        *b_Harm_ActAnScint_hit_xcellg;   //!
   TBranch        *b_Harm_ActAnScint_hit_ycellg;   //!
   TBranch        *b_Harm_ActAnScint_hit_zcellg;   //!
   TBranch        *b_Harm_ActAnScint_hit_xhit;   //!
   TBranch        *b_Harm_ActAnScint_hit_yhit;   //!
   TBranch        *b_Harm_ActAnScint_hit_zhit;   //!
   TBranch        *b_Harm_ActAnScint_hit_xhitg;   //!
   TBranch        *b_Harm_ActAnScint_hit_yhitg;   //!
   TBranch        *b_Harm_ActAnScint_hit_zhitg;   //!
   TBranch        *b_Harm_ActAnScint_hit_sumedep;   //!
   TBranch        *b_Harm_ActAnScint_hit_tavg;   //!
   TBranch        *b_Harm_ActAnScint_hit_trms;   //!
   TBranch        *b_Harm_ActAnScint_hit_tmin;   //!
   TBranch        *b_Harm_ActAnScint_hit_tmax;   //!
   TBranch        *b_Harm_ActAnScint_hit_otridx;   //!
   TBranch        *b_Harm_ActAnScint_hit_ptridx;   //!
   TBranch        *b_Harm_ActAnScint_hit_sdtridx;   //!
   TBranch        *b_Harm_CEPolFront_hit_nhits;   //!
   TBranch        *b_Harm_CEPolFront_hit_plane;   //!
   TBranch        *b_Harm_CEPolFront_hit_strip;   //!
   TBranch        *b_Harm_CEPolFront_hit_x;   //!
   TBranch        *b_Harm_CEPolFront_hit_y;   //!
   TBranch        *b_Harm_CEPolFront_hit_z;   //!
   TBranch        *b_Harm_CEPolFront_hit_polx;   //!
   TBranch        *b_Harm_CEPolFront_hit_poly;   //!
   TBranch        *b_Harm_CEPolFront_hit_polz;   //!
   TBranch        *b_Harm_CEPolFront_hit_t;   //!
   TBranch        *b_Harm_CEPolFront_hit_trms;   //!
   TBranch        *b_Harm_CEPolFront_hit_tmin;   //!
   TBranch        *b_Harm_CEPolFront_hit_tmax;   //!
   TBranch        *b_Harm_CEPolFront_hit_tx;   //!
   TBranch        *b_Harm_CEPolFront_hit_ty;   //!
   TBranch        *b_Harm_CEPolFront_hit_xin;   //!
   TBranch        *b_Harm_CEPolFront_hit_yin;   //!
   TBranch        *b_Harm_CEPolFront_hit_zin;   //!
   TBranch        *b_Harm_CEPolFront_hit_xout;   //!
   TBranch        *b_Harm_CEPolFront_hit_yout;   //!
   TBranch        *b_Harm_CEPolFront_hit_zout;   //!
   TBranch        *b_Harm_CEPolFront_hit_txp;   //!
   TBranch        *b_Harm_CEPolFront_hit_typ;   //!
   TBranch        *b_Harm_CEPolFront_hit_xg;   //!
   TBranch        *b_Harm_CEPolFront_hit_yg;   //!
   TBranch        *b_Harm_CEPolFront_hit_zg;   //!
   TBranch        *b_Harm_CEPolFront_hit_trid;   //!
   TBranch        *b_Harm_CEPolFront_hit_mid;   //!
   TBranch        *b_Harm_CEPolFront_hit_pid;   //!
   TBranch        *b_Harm_CEPolFront_hit_vx;   //!
   TBranch        *b_Harm_CEPolFront_hit_vy;   //!
   TBranch        *b_Harm_CEPolFront_hit_vz;   //!
   TBranch        *b_Harm_CEPolFront_hit_p;   //!
   TBranch        *b_Harm_CEPolFront_hit_edep;   //!
   TBranch        *b_Harm_CEPolFront_hit_beta;   //!
   TBranch        *b_Harm_CEPolFront_hit_otridx;   //!
   TBranch        *b_Harm_CEPolFront_hit_ptridx;   //!
   TBranch        *b_Harm_CEPolFront_hit_sdtridx;   //!
   TBranch        *b_Harm_CEPolFront_Track_ntracks;   //!
   TBranch        *b_Harm_CEPolFront_Track_TID;   //!
   TBranch        *b_Harm_CEPolFront_Track_PID;   //!
   TBranch        *b_Harm_CEPolFront_Track_MID;   //!
   TBranch        *b_Harm_CEPolFront_Track_NumHits;   //!
   TBranch        *b_Harm_CEPolFront_Track_NumPlanes;   //!
   TBranch        *b_Harm_CEPolFront_Track_NDF;   //!
   TBranch        *b_Harm_CEPolFront_Track_Chi2fit;   //!
   TBranch        *b_Harm_CEPolFront_Track_Chi2true;   //!
   TBranch        *b_Harm_CEPolFront_Track_X;   //!
   TBranch        *b_Harm_CEPolFront_Track_Y;   //!
   TBranch        *b_Harm_CEPolFront_Track_Xp;   //!
   TBranch        *b_Harm_CEPolFront_Track_Yp;   //!
   TBranch        *b_Harm_CEPolFront_Track_T;   //!
   TBranch        *b_Harm_CEPolFront_Track_P;   //!
   TBranch        *b_Harm_CEPolFront_Track_Sx;   //!
   TBranch        *b_Harm_CEPolFront_Track_Sy;   //!
   TBranch        *b_Harm_CEPolFront_Track_Sz;   //!
   TBranch        *b_Harm_CEPolFront_Track_Xfit;   //!
   TBranch        *b_Harm_CEPolFront_Track_Yfit;   //!
   TBranch        *b_Harm_CEPolFront_Track_Xpfit;   //!
   TBranch        *b_Harm_CEPolFront_Track_Ypfit;   //!
   TBranch        *b_Harm_CEPolFront_Track_otridx;   //!
   TBranch        *b_Harm_CEPolFront_Track_ptridx;   //!
   TBranch        *b_Harm_CEPolFront_Track_sdtridx;   //!
   TBranch        *b_Harm_CEPolRear_hit_nhits;   //!
   TBranch        *b_Harm_CEPolRear_hit_plane;   //!
   TBranch        *b_Harm_CEPolRear_hit_strip;   //!
   TBranch        *b_Harm_CEPolRear_hit_x;   //!
   TBranch        *b_Harm_CEPolRear_hit_y;   //!
   TBranch        *b_Harm_CEPolRear_hit_z;   //!
   TBranch        *b_Harm_CEPolRear_hit_polx;   //!
   TBranch        *b_Harm_CEPolRear_hit_poly;   //!
   TBranch        *b_Harm_CEPolRear_hit_polz;   //!
   TBranch        *b_Harm_CEPolRear_hit_t;   //!
   TBranch        *b_Harm_CEPolRear_hit_trms;   //!
   TBranch        *b_Harm_CEPolRear_hit_tmin;   //!
   TBranch        *b_Harm_CEPolRear_hit_tmax;   //!
   TBranch        *b_Harm_CEPolRear_hit_tx;   //!
   TBranch        *b_Harm_CEPolRear_hit_ty;   //!
   TBranch        *b_Harm_CEPolRear_hit_xin;   //!
   TBranch        *b_Harm_CEPolRear_hit_yin;   //!
   TBranch        *b_Harm_CEPolRear_hit_zin;   //!
   TBranch        *b_Harm_CEPolRear_hit_xout;   //!
   TBranch        *b_Harm_CEPolRear_hit_yout;   //!
   TBranch        *b_Harm_CEPolRear_hit_zout;   //!
   TBranch        *b_Harm_CEPolRear_hit_txp;   //!
   TBranch        *b_Harm_CEPolRear_hit_typ;   //!
   TBranch        *b_Harm_CEPolRear_hit_xg;   //!
   TBranch        *b_Harm_CEPolRear_hit_yg;   //!
   TBranch        *b_Harm_CEPolRear_hit_zg;   //!
   TBranch        *b_Harm_CEPolRear_hit_trid;   //!
   TBranch        *b_Harm_CEPolRear_hit_mid;   //!
   TBranch        *b_Harm_CEPolRear_hit_pid;   //!
   TBranch        *b_Harm_CEPolRear_hit_vx;   //!
   TBranch        *b_Harm_CEPolRear_hit_vy;   //!
   TBranch        *b_Harm_CEPolRear_hit_vz;   //!
   TBranch        *b_Harm_CEPolRear_hit_p;   //!
   TBranch        *b_Harm_CEPolRear_hit_edep;   //!
   TBranch        *b_Harm_CEPolRear_hit_beta;   //!
   TBranch        *b_Harm_CEPolRear_hit_otridx;   //!
   TBranch        *b_Harm_CEPolRear_hit_ptridx;   //!
   TBranch        *b_Harm_CEPolRear_hit_sdtridx;   //!
   TBranch        *b_Harm_CEPolRear_Track_ntracks;   //!
   TBranch        *b_Harm_CEPolRear_Track_TID;   //!
   TBranch        *b_Harm_CEPolRear_Track_PID;   //!
   TBranch        *b_Harm_CEPolRear_Track_MID;   //!
   TBranch        *b_Harm_CEPolRear_Track_NumHits;   //!
   TBranch        *b_Harm_CEPolRear_Track_NumPlanes;   //!
   TBranch        *b_Harm_CEPolRear_Track_NDF;   //!
   TBranch        *b_Harm_CEPolRear_Track_Chi2fit;   //!
   TBranch        *b_Harm_CEPolRear_Track_Chi2true;   //!
   TBranch        *b_Harm_CEPolRear_Track_X;   //!
   TBranch        *b_Harm_CEPolRear_Track_Y;   //!
   TBranch        *b_Harm_CEPolRear_Track_Xp;   //!
   TBranch        *b_Harm_CEPolRear_Track_Yp;   //!
   TBranch        *b_Harm_CEPolRear_Track_T;   //!
   TBranch        *b_Harm_CEPolRear_Track_P;   //!
   TBranch        *b_Harm_CEPolRear_Track_Sx;   //!
   TBranch        *b_Harm_CEPolRear_Track_Sy;   //!
   TBranch        *b_Harm_CEPolRear_Track_Sz;   //!
   TBranch        *b_Harm_CEPolRear_Track_Xfit;   //!
   TBranch        *b_Harm_CEPolRear_Track_Yfit;   //!
   TBranch        *b_Harm_CEPolRear_Track_Xpfit;   //!
   TBranch        *b_Harm_CEPolRear_Track_Ypfit;   //!
   TBranch        *b_Harm_CEPolRear_Track_otridx;   //!
   TBranch        *b_Harm_CEPolRear_Track_ptridx;   //!
   TBranch        *b_Harm_CEPolRear_Track_sdtridx;   //!
   TBranch        *b_Harm_PRPolGEMFarSide_hit_nhits;   //!
   TBranch        *b_Harm_PRPolGEMFarSide_hit_plane;   //!
   TBranch        *b_Harm_PRPolGEMFarSide_hit_strip;   //!
   TBranch        *b_Harm_PRPolGEMFarSide_hit_x;   //!
   TBranch        *b_Harm_PRPolGEMFarSide_hit_y;   //!
   TBranch        *b_Harm_PRPolGEMFarSide_hit_z;   //!
   TBranch        *b_Harm_PRPolGEMFarSide_hit_polx;   //!
   TBranch        *b_Harm_PRPolGEMFarSide_hit_poly;   //!
   TBranch        *b_Harm_PRPolGEMFarSide_hit_polz;   //!
   TBranch        *b_Harm_PRPolGEMFarSide_hit_t;   //!
   TBranch        *b_Harm_PRPolGEMFarSide_hit_trms;   //!
   TBranch        *b_Harm_PRPolGEMFarSide_hit_tmin;   //!
   TBranch        *b_Harm_PRPolGEMFarSide_hit_tmax;   //!
   TBranch        *b_Harm_PRPolGEMFarSide_hit_tx;   //!
   TBranch        *b_Harm_PRPolGEMFarSide_hit_ty;   //!
   TBranch        *b_Harm_PRPolGEMFarSide_hit_xin;   //!
   TBranch        *b_Harm_PRPolGEMFarSide_hit_yin;   //!
   TBranch        *b_Harm_PRPolGEMFarSide_hit_zin;   //!
   TBranch        *b_Harm_PRPolGEMFarSide_hit_xout;   //!
   TBranch        *b_Harm_PRPolGEMFarSide_hit_yout;   //!
   TBranch        *b_Harm_PRPolGEMFarSide_hit_zout;   //!
   TBranch        *b_Harm_PRPolGEMFarSide_hit_txp;   //!
   TBranch        *b_Harm_PRPolGEMFarSide_hit_typ;   //!
   TBranch        *b_Harm_PRPolGEMFarSide_hit_xg;   //!
   TBranch        *b_Harm_PRPolGEMFarSide_hit_yg;   //!
   TBranch        *b_Harm_PRPolGEMFarSide_hit_zg;   //!
   TBranch        *b_Harm_PRPolGEMFarSide_hit_trid;   //!
   TBranch        *b_Harm_PRPolGEMFarSide_hit_mid;   //!
   TBranch        *b_Harm_PRPolGEMFarSide_hit_pid;   //!
   TBranch        *b_Harm_PRPolGEMFarSide_hit_vx;   //!
   TBranch        *b_Harm_PRPolGEMFarSide_hit_vy;   //!
   TBranch        *b_Harm_PRPolGEMFarSide_hit_vz;   //!
   TBranch        *b_Harm_PRPolGEMFarSide_hit_p;   //!
   TBranch        *b_Harm_PRPolGEMFarSide_hit_edep;   //!
   TBranch        *b_Harm_PRPolGEMFarSide_hit_beta;   //!
   TBranch        *b_Harm_PRPolGEMFarSide_hit_otridx;   //!
   TBranch        *b_Harm_PRPolGEMFarSide_hit_ptridx;   //!
   TBranch        *b_Harm_PRPolGEMFarSide_hit_sdtridx;   //!
   TBranch        *b_Harm_PRPolGEMFarSide_Track_ntracks;   //!
   TBranch        *b_Harm_PRPolGEMFarSide_Track_TID;   //!
   TBranch        *b_Harm_PRPolGEMFarSide_Track_PID;   //!
   TBranch        *b_Harm_PRPolGEMFarSide_Track_MID;   //!
   TBranch        *b_Harm_PRPolGEMFarSide_Track_NumHits;   //!
   TBranch        *b_Harm_PRPolGEMFarSide_Track_NumPlanes;   //!
   TBranch        *b_Harm_PRPolGEMFarSide_Track_NDF;   //!
   TBranch        *b_Harm_PRPolGEMFarSide_Track_Chi2fit;   //!
   TBranch        *b_Harm_PRPolGEMFarSide_Track_Chi2true;   //!
   TBranch        *b_Harm_PRPolGEMFarSide_Track_X;   //!
   TBranch        *b_Harm_PRPolGEMFarSide_Track_Y;   //!
   TBranch        *b_Harm_PRPolGEMFarSide_Track_Xp;   //!
   TBranch        *b_Harm_PRPolGEMFarSide_Track_Yp;   //!
   TBranch        *b_Harm_PRPolGEMFarSide_Track_T;   //!
   TBranch        *b_Harm_PRPolGEMFarSide_Track_P;   //!
   TBranch        *b_Harm_PRPolGEMFarSide_Track_Sx;   //!
   TBranch        *b_Harm_PRPolGEMFarSide_Track_Sy;   //!
   TBranch        *b_Harm_PRPolGEMFarSide_Track_Sz;   //!
   TBranch        *b_Harm_PRPolGEMFarSide_Track_Xfit;   //!
   TBranch        *b_Harm_PRPolGEMFarSide_Track_Yfit;   //!
   TBranch        *b_Harm_PRPolGEMFarSide_Track_Xpfit;   //!
   TBranch        *b_Harm_PRPolGEMFarSide_Track_Ypfit;   //!
   TBranch        *b_Harm_PRPolGEMFarSide_Track_otridx;   //!
   TBranch        *b_Harm_PRPolGEMFarSide_Track_ptridx;   //!
   TBranch        *b_Harm_PRPolGEMFarSide_Track_sdtridx;   //!
   TBranch        *b_Harm_PRPolScintFarSide_det_esum;   //!
   TBranch        *b_Harm_PRPolScintFarSide_hit_nhits;   //!
   TBranch        *b_Harm_PRPolScintFarSide_hit_row;   //!
   TBranch        *b_Harm_PRPolScintFarSide_hit_col;   //!
   TBranch        *b_Harm_PRPolScintFarSide_hit_cell;   //!
   TBranch        *b_Harm_PRPolScintFarSide_hit_plane;   //!
   TBranch        *b_Harm_PRPolScintFarSide_hit_wire;   //!
   TBranch        *b_Harm_PRPolScintFarSide_hit_xcell;   //!
   TBranch        *b_Harm_PRPolScintFarSide_hit_ycell;   //!
   TBranch        *b_Harm_PRPolScintFarSide_hit_zcell;   //!
   TBranch        *b_Harm_PRPolScintFarSide_hit_xcellg;   //!
   TBranch        *b_Harm_PRPolScintFarSide_hit_ycellg;   //!
   TBranch        *b_Harm_PRPolScintFarSide_hit_zcellg;   //!
   TBranch        *b_Harm_PRPolScintFarSide_hit_xhit;   //!
   TBranch        *b_Harm_PRPolScintFarSide_hit_yhit;   //!
   TBranch        *b_Harm_PRPolScintFarSide_hit_zhit;   //!
   TBranch        *b_Harm_PRPolScintFarSide_hit_xhitg;   //!
   TBranch        *b_Harm_PRPolScintFarSide_hit_yhitg;   //!
   TBranch        *b_Harm_PRPolScintFarSide_hit_zhitg;   //!
   TBranch        *b_Harm_PRPolScintFarSide_hit_sumedep;   //!
   TBranch        *b_Harm_PRPolScintFarSide_hit_tavg;   //!
   TBranch        *b_Harm_PRPolScintFarSide_hit_trms;   //!
   TBranch        *b_Harm_PRPolScintFarSide_hit_tmin;   //!
   TBranch        *b_Harm_PRPolScintFarSide_hit_tmax;   //!
   TBranch        *b_Harm_PRPolScintFarSide_hit_otridx;   //!
   TBranch        *b_Harm_PRPolScintFarSide_hit_ptridx;   //!
   TBranch        *b_Harm_PRPolScintFarSide_hit_sdtridx;   //!
   TBranch        *b_Harm_ActAn_dighit_nchan;   //!
   TBranch        *b_Harm_ActAn_dighit_chan;   //!
   TBranch        *b_Harm_ActAn_dighit_adc;   //!
   TBranch        *b_Harm_ActAn_dighit_samp;   //!
   TBranch        *b_Harm_ActAn_dighit_tdc;   //!
   TBranch        *b_Harm_PRPolScintFarSide_dighit_nchan;   //!
   TBranch        *b_Harm_PRPolScintFarSide_dighit_chan;   //!
   TBranch        *b_Harm_PRPolScintFarSide_dighit_adc;   //!
   TBranch        *b_Harm_PRPolScintFarSide_dighit_tdc_l;   //!
   TBranch        *b_Harm_PRPolScintFarSide_dighit_tdc_t;   //!
   TBranch        *b_Harm_CEPolFront_dighit_nstrips;   //!
   TBranch        *b_Harm_CEPolFront_dighit_module;   //!
   TBranch        *b_Harm_CEPolFront_dighit_strip;   //!
   TBranch        *b_Harm_CEPolFront_dighit_adc;   //!
   TBranch        *b_Harm_CEPolFront_dighit_samp;   //!
   TBranch        *b_Harm_CEPolRear_dighit_nstrips;   //!
   TBranch        *b_Harm_CEPolRear_dighit_module;   //!
   TBranch        *b_Harm_CEPolRear_dighit_strip;   //!
   TBranch        *b_Harm_CEPolRear_dighit_adc;   //!
   TBranch        *b_Harm_CEPolRear_dighit_samp;   //!
   TBranch        *b_Harm_PRPolGEMFarSide_dighit_nstrips;   //!
   TBranch        *b_Harm_PRPolGEMFarSide_dighit_module;   //!
   TBranch        *b_Harm_PRPolGEMFarSide_dighit_strip;   //!
   TBranch        *b_Harm_PRPolGEMFarSide_dighit_adc;   //!
   TBranch        *b_Harm_PRPolGEMFarSide_dighit_samp;   //!

   genrp_tree_digitized(TTree *tree=0);
   virtual ~genrp_tree_digitized();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef genrp_tree_digitized_cxx
genrp_tree_digitized::genrp_tree_digitized(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("genrp_elastic_gc_job_0.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("genrp_elastic_gc_job_0.root");
      }
      f->GetObject("T",tree);

   }
   Init(tree);
}

genrp_tree_digitized::~genrp_tree_digitized()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t genrp_tree_digitized::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t genrp_tree_digitized::LoadTree(Long64_t entry)
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

void genrp_tree_digitized::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   Harm_ActAnScint_hit_row = 0;
   Harm_ActAnScint_hit_col = 0;
   Harm_ActAnScint_hit_cell = 0;
   Harm_ActAnScint_hit_plane = 0;
   Harm_ActAnScint_hit_wire = 0;
   Harm_ActAnScint_hit_xcell = 0;
   Harm_ActAnScint_hit_ycell = 0;
   Harm_ActAnScint_hit_zcell = 0;
   Harm_ActAnScint_hit_xcellg = 0;
   Harm_ActAnScint_hit_ycellg = 0;
   Harm_ActAnScint_hit_zcellg = 0;
   Harm_ActAnScint_hit_xhit = 0;
   Harm_ActAnScint_hit_yhit = 0;
   Harm_ActAnScint_hit_zhit = 0;
   Harm_ActAnScint_hit_xhitg = 0;
   Harm_ActAnScint_hit_yhitg = 0;
   Harm_ActAnScint_hit_zhitg = 0;
   Harm_ActAnScint_hit_sumedep = 0;
   Harm_ActAnScint_hit_tavg = 0;
   Harm_ActAnScint_hit_trms = 0;
   Harm_ActAnScint_hit_tmin = 0;
   Harm_ActAnScint_hit_tmax = 0;
   Harm_ActAnScint_hit_otridx = 0;
   Harm_ActAnScint_hit_ptridx = 0;
   Harm_ActAnScint_hit_sdtridx = 0;
   Harm_CEPolFront_hit_plane = 0;
   Harm_CEPolFront_hit_strip = 0;
   Harm_CEPolFront_hit_x = 0;
   Harm_CEPolFront_hit_y = 0;
   Harm_CEPolFront_hit_z = 0;
   Harm_CEPolFront_hit_polx = 0;
   Harm_CEPolFront_hit_poly = 0;
   Harm_CEPolFront_hit_polz = 0;
   Harm_CEPolFront_hit_t = 0;
   Harm_CEPolFront_hit_trms = 0;
   Harm_CEPolFront_hit_tmin = 0;
   Harm_CEPolFront_hit_tmax = 0;
   Harm_CEPolFront_hit_tx = 0;
   Harm_CEPolFront_hit_ty = 0;
   Harm_CEPolFront_hit_xin = 0;
   Harm_CEPolFront_hit_yin = 0;
   Harm_CEPolFront_hit_zin = 0;
   Harm_CEPolFront_hit_xout = 0;
   Harm_CEPolFront_hit_yout = 0;
   Harm_CEPolFront_hit_zout = 0;
   Harm_CEPolFront_hit_txp = 0;
   Harm_CEPolFront_hit_typ = 0;
   Harm_CEPolFront_hit_xg = 0;
   Harm_CEPolFront_hit_yg = 0;
   Harm_CEPolFront_hit_zg = 0;
   Harm_CEPolFront_hit_trid = 0;
   Harm_CEPolFront_hit_mid = 0;
   Harm_CEPolFront_hit_pid = 0;
   Harm_CEPolFront_hit_vx = 0;
   Harm_CEPolFront_hit_vy = 0;
   Harm_CEPolFront_hit_vz = 0;
   Harm_CEPolFront_hit_p = 0;
   Harm_CEPolFront_hit_edep = 0;
   Harm_CEPolFront_hit_beta = 0;
   Harm_CEPolFront_hit_otridx = 0;
   Harm_CEPolFront_hit_ptridx = 0;
   Harm_CEPolFront_hit_sdtridx = 0;
   Harm_CEPolFront_Track_TID = 0;
   Harm_CEPolFront_Track_PID = 0;
   Harm_CEPolFront_Track_MID = 0;
   Harm_CEPolFront_Track_NumHits = 0;
   Harm_CEPolFront_Track_NumPlanes = 0;
   Harm_CEPolFront_Track_NDF = 0;
   Harm_CEPolFront_Track_Chi2fit = 0;
   Harm_CEPolFront_Track_Chi2true = 0;
   Harm_CEPolFront_Track_X = 0;
   Harm_CEPolFront_Track_Y = 0;
   Harm_CEPolFront_Track_Xp = 0;
   Harm_CEPolFront_Track_Yp = 0;
   Harm_CEPolFront_Track_T = 0;
   Harm_CEPolFront_Track_P = 0;
   Harm_CEPolFront_Track_Sx = 0;
   Harm_CEPolFront_Track_Sy = 0;
   Harm_CEPolFront_Track_Sz = 0;
   Harm_CEPolFront_Track_Xfit = 0;
   Harm_CEPolFront_Track_Yfit = 0;
   Harm_CEPolFront_Track_Xpfit = 0;
   Harm_CEPolFront_Track_Ypfit = 0;
   Harm_CEPolFront_Track_otridx = 0;
   Harm_CEPolFront_Track_ptridx = 0;
   Harm_CEPolFront_Track_sdtridx = 0;
   Harm_CEPolRear_hit_plane = 0;
   Harm_CEPolRear_hit_strip = 0;
   Harm_CEPolRear_hit_x = 0;
   Harm_CEPolRear_hit_y = 0;
   Harm_CEPolRear_hit_z = 0;
   Harm_CEPolRear_hit_polx = 0;
   Harm_CEPolRear_hit_poly = 0;
   Harm_CEPolRear_hit_polz = 0;
   Harm_CEPolRear_hit_t = 0;
   Harm_CEPolRear_hit_trms = 0;
   Harm_CEPolRear_hit_tmin = 0;
   Harm_CEPolRear_hit_tmax = 0;
   Harm_CEPolRear_hit_tx = 0;
   Harm_CEPolRear_hit_ty = 0;
   Harm_CEPolRear_hit_xin = 0;
   Harm_CEPolRear_hit_yin = 0;
   Harm_CEPolRear_hit_zin = 0;
   Harm_CEPolRear_hit_xout = 0;
   Harm_CEPolRear_hit_yout = 0;
   Harm_CEPolRear_hit_zout = 0;
   Harm_CEPolRear_hit_txp = 0;
   Harm_CEPolRear_hit_typ = 0;
   Harm_CEPolRear_hit_xg = 0;
   Harm_CEPolRear_hit_yg = 0;
   Harm_CEPolRear_hit_zg = 0;
   Harm_CEPolRear_hit_trid = 0;
   Harm_CEPolRear_hit_mid = 0;
   Harm_CEPolRear_hit_pid = 0;
   Harm_CEPolRear_hit_vx = 0;
   Harm_CEPolRear_hit_vy = 0;
   Harm_CEPolRear_hit_vz = 0;
   Harm_CEPolRear_hit_p = 0;
   Harm_CEPolRear_hit_edep = 0;
   Harm_CEPolRear_hit_beta = 0;
   Harm_CEPolRear_hit_otridx = 0;
   Harm_CEPolRear_hit_ptridx = 0;
   Harm_CEPolRear_hit_sdtridx = 0;
   Harm_CEPolRear_Track_TID = 0;
   Harm_CEPolRear_Track_PID = 0;
   Harm_CEPolRear_Track_MID = 0;
   Harm_CEPolRear_Track_NumHits = 0;
   Harm_CEPolRear_Track_NumPlanes = 0;
   Harm_CEPolRear_Track_NDF = 0;
   Harm_CEPolRear_Track_Chi2fit = 0;
   Harm_CEPolRear_Track_Chi2true = 0;
   Harm_CEPolRear_Track_X = 0;
   Harm_CEPolRear_Track_Y = 0;
   Harm_CEPolRear_Track_Xp = 0;
   Harm_CEPolRear_Track_Yp = 0;
   Harm_CEPolRear_Track_T = 0;
   Harm_CEPolRear_Track_P = 0;
   Harm_CEPolRear_Track_Sx = 0;
   Harm_CEPolRear_Track_Sy = 0;
   Harm_CEPolRear_Track_Sz = 0;
   Harm_CEPolRear_Track_Xfit = 0;
   Harm_CEPolRear_Track_Yfit = 0;
   Harm_CEPolRear_Track_Xpfit = 0;
   Harm_CEPolRear_Track_Ypfit = 0;
   Harm_CEPolRear_Track_otridx = 0;
   Harm_CEPolRear_Track_ptridx = 0;
   Harm_CEPolRear_Track_sdtridx = 0;
   Harm_PRPolGEMFarSide_hit_plane = 0;
   Harm_PRPolGEMFarSide_hit_strip = 0;
   Harm_PRPolGEMFarSide_hit_x = 0;
   Harm_PRPolGEMFarSide_hit_y = 0;
   Harm_PRPolGEMFarSide_hit_z = 0;
   Harm_PRPolGEMFarSide_hit_polx = 0;
   Harm_PRPolGEMFarSide_hit_poly = 0;
   Harm_PRPolGEMFarSide_hit_polz = 0;
   Harm_PRPolGEMFarSide_hit_t = 0;
   Harm_PRPolGEMFarSide_hit_trms = 0;
   Harm_PRPolGEMFarSide_hit_tmin = 0;
   Harm_PRPolGEMFarSide_hit_tmax = 0;
   Harm_PRPolGEMFarSide_hit_tx = 0;
   Harm_PRPolGEMFarSide_hit_ty = 0;
   Harm_PRPolGEMFarSide_hit_xin = 0;
   Harm_PRPolGEMFarSide_hit_yin = 0;
   Harm_PRPolGEMFarSide_hit_zin = 0;
   Harm_PRPolGEMFarSide_hit_xout = 0;
   Harm_PRPolGEMFarSide_hit_yout = 0;
   Harm_PRPolGEMFarSide_hit_zout = 0;
   Harm_PRPolGEMFarSide_hit_txp = 0;
   Harm_PRPolGEMFarSide_hit_typ = 0;
   Harm_PRPolGEMFarSide_hit_xg = 0;
   Harm_PRPolGEMFarSide_hit_yg = 0;
   Harm_PRPolGEMFarSide_hit_zg = 0;
   Harm_PRPolGEMFarSide_hit_trid = 0;
   Harm_PRPolGEMFarSide_hit_mid = 0;
   Harm_PRPolGEMFarSide_hit_pid = 0;
   Harm_PRPolGEMFarSide_hit_vx = 0;
   Harm_PRPolGEMFarSide_hit_vy = 0;
   Harm_PRPolGEMFarSide_hit_vz = 0;
   Harm_PRPolGEMFarSide_hit_p = 0;
   Harm_PRPolGEMFarSide_hit_edep = 0;
   Harm_PRPolGEMFarSide_hit_beta = 0;
   Harm_PRPolGEMFarSide_hit_otridx = 0;
   Harm_PRPolGEMFarSide_hit_ptridx = 0;
   Harm_PRPolGEMFarSide_hit_sdtridx = 0;
   Harm_PRPolGEMFarSide_Track_TID = 0;
   Harm_PRPolGEMFarSide_Track_PID = 0;
   Harm_PRPolGEMFarSide_Track_MID = 0;
   Harm_PRPolGEMFarSide_Track_NumHits = 0;
   Harm_PRPolGEMFarSide_Track_NumPlanes = 0;
   Harm_PRPolGEMFarSide_Track_NDF = 0;
   Harm_PRPolGEMFarSide_Track_Chi2fit = 0;
   Harm_PRPolGEMFarSide_Track_Chi2true = 0;
   Harm_PRPolGEMFarSide_Track_X = 0;
   Harm_PRPolGEMFarSide_Track_Y = 0;
   Harm_PRPolGEMFarSide_Track_Xp = 0;
   Harm_PRPolGEMFarSide_Track_Yp = 0;
   Harm_PRPolGEMFarSide_Track_T = 0;
   Harm_PRPolGEMFarSide_Track_P = 0;
   Harm_PRPolGEMFarSide_Track_Sx = 0;
   Harm_PRPolGEMFarSide_Track_Sy = 0;
   Harm_PRPolGEMFarSide_Track_Sz = 0;
   Harm_PRPolGEMFarSide_Track_Xfit = 0;
   Harm_PRPolGEMFarSide_Track_Yfit = 0;
   Harm_PRPolGEMFarSide_Track_Xpfit = 0;
   Harm_PRPolGEMFarSide_Track_Ypfit = 0;
   Harm_PRPolGEMFarSide_Track_otridx = 0;
   Harm_PRPolGEMFarSide_Track_ptridx = 0;
   Harm_PRPolGEMFarSide_Track_sdtridx = 0;
   Harm_PRPolScintFarSide_hit_row = 0;
   Harm_PRPolScintFarSide_hit_col = 0;
   Harm_PRPolScintFarSide_hit_cell = 0;
   Harm_PRPolScintFarSide_hit_plane = 0;
   Harm_PRPolScintFarSide_hit_wire = 0;
   Harm_PRPolScintFarSide_hit_xcell = 0;
   Harm_PRPolScintFarSide_hit_ycell = 0;
   Harm_PRPolScintFarSide_hit_zcell = 0;
   Harm_PRPolScintFarSide_hit_xcellg = 0;
   Harm_PRPolScintFarSide_hit_ycellg = 0;
   Harm_PRPolScintFarSide_hit_zcellg = 0;
   Harm_PRPolScintFarSide_hit_xhit = 0;
   Harm_PRPolScintFarSide_hit_yhit = 0;
   Harm_PRPolScintFarSide_hit_zhit = 0;
   Harm_PRPolScintFarSide_hit_xhitg = 0;
   Harm_PRPolScintFarSide_hit_yhitg = 0;
   Harm_PRPolScintFarSide_hit_zhitg = 0;
   Harm_PRPolScintFarSide_hit_sumedep = 0;
   Harm_PRPolScintFarSide_hit_tavg = 0;
   Harm_PRPolScintFarSide_hit_trms = 0;
   Harm_PRPolScintFarSide_hit_tmin = 0;
   Harm_PRPolScintFarSide_hit_tmax = 0;
   Harm_PRPolScintFarSide_hit_otridx = 0;
   Harm_PRPolScintFarSide_hit_ptridx = 0;
   Harm_PRPolScintFarSide_hit_sdtridx = 0;
   Harm_ActAn_dighit_chan = 0;
   Harm_ActAn_dighit_adc = 0;
   Harm_ActAn_dighit_samp = 0;
   Harm_ActAn_dighit_tdc = 0;
   Harm_PRPolScintFarSide_dighit_chan = 0;
   Harm_PRPolScintFarSide_dighit_adc = 0;
   Harm_PRPolScintFarSide_dighit_tdc_l = 0;
   Harm_PRPolScintFarSide_dighit_tdc_t = 0;
   Harm_CEPolFront_dighit_module = 0;
   Harm_CEPolFront_dighit_strip = 0;
   Harm_CEPolFront_dighit_adc = 0;
   Harm_CEPolFront_dighit_samp = 0;
   Harm_CEPolRear_dighit_module = 0;
   Harm_CEPolRear_dighit_strip = 0;
   Harm_CEPolRear_dighit_adc = 0;
   Harm_CEPolRear_dighit_samp = 0;
   Harm_PRPolGEMFarSide_dighit_module = 0;
   Harm_PRPolGEMFarSide_dighit_strip = 0;
   Harm_PRPolGEMFarSide_dighit_adc = 0;
   Harm_PRPolGEMFarSide_dighit_samp = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("Harm.ActAnScint.det.esum", &Harm_ActAnScint_det_esum, &b_Harm_ActAnScint_det_esum);
   fChain->SetBranchAddress("Harm.ActAnScint.hit.nhits", &Harm_ActAnScint_hit_nhits, &b_Harm_ActAnScint_hit_nhits);
   fChain->SetBranchAddress("Harm.ActAnScint.hit.row", &Harm_ActAnScint_hit_row, &b_Harm_ActAnScint_hit_row);
   fChain->SetBranchAddress("Harm.ActAnScint.hit.col", &Harm_ActAnScint_hit_col, &b_Harm_ActAnScint_hit_col);
   fChain->SetBranchAddress("Harm.ActAnScint.hit.cell", &Harm_ActAnScint_hit_cell, &b_Harm_ActAnScint_hit_cell);
   fChain->SetBranchAddress("Harm.ActAnScint.hit.plane", &Harm_ActAnScint_hit_plane, &b_Harm_ActAnScint_hit_plane);
   fChain->SetBranchAddress("Harm.ActAnScint.hit.wire", &Harm_ActAnScint_hit_wire, &b_Harm_ActAnScint_hit_wire);
   fChain->SetBranchAddress("Harm.ActAnScint.hit.xcell", &Harm_ActAnScint_hit_xcell, &b_Harm_ActAnScint_hit_xcell);
   fChain->SetBranchAddress("Harm.ActAnScint.hit.ycell", &Harm_ActAnScint_hit_ycell, &b_Harm_ActAnScint_hit_ycell);
   fChain->SetBranchAddress("Harm.ActAnScint.hit.zcell", &Harm_ActAnScint_hit_zcell, &b_Harm_ActAnScint_hit_zcell);
   fChain->SetBranchAddress("Harm.ActAnScint.hit.xcellg", &Harm_ActAnScint_hit_xcellg, &b_Harm_ActAnScint_hit_xcellg);
   fChain->SetBranchAddress("Harm.ActAnScint.hit.ycellg", &Harm_ActAnScint_hit_ycellg, &b_Harm_ActAnScint_hit_ycellg);
   fChain->SetBranchAddress("Harm.ActAnScint.hit.zcellg", &Harm_ActAnScint_hit_zcellg, &b_Harm_ActAnScint_hit_zcellg);
   fChain->SetBranchAddress("Harm.ActAnScint.hit.xhit", &Harm_ActAnScint_hit_xhit, &b_Harm_ActAnScint_hit_xhit);
   fChain->SetBranchAddress("Harm.ActAnScint.hit.yhit", &Harm_ActAnScint_hit_yhit, &b_Harm_ActAnScint_hit_yhit);
   fChain->SetBranchAddress("Harm.ActAnScint.hit.zhit", &Harm_ActAnScint_hit_zhit, &b_Harm_ActAnScint_hit_zhit);
   fChain->SetBranchAddress("Harm.ActAnScint.hit.xhitg", &Harm_ActAnScint_hit_xhitg, &b_Harm_ActAnScint_hit_xhitg);
   fChain->SetBranchAddress("Harm.ActAnScint.hit.yhitg", &Harm_ActAnScint_hit_yhitg, &b_Harm_ActAnScint_hit_yhitg);
   fChain->SetBranchAddress("Harm.ActAnScint.hit.zhitg", &Harm_ActAnScint_hit_zhitg, &b_Harm_ActAnScint_hit_zhitg);
   fChain->SetBranchAddress("Harm.ActAnScint.hit.sumedep", &Harm_ActAnScint_hit_sumedep, &b_Harm_ActAnScint_hit_sumedep);
   fChain->SetBranchAddress("Harm.ActAnScint.hit.tavg", &Harm_ActAnScint_hit_tavg, &b_Harm_ActAnScint_hit_tavg);
   fChain->SetBranchAddress("Harm.ActAnScint.hit.trms", &Harm_ActAnScint_hit_trms, &b_Harm_ActAnScint_hit_trms);
   fChain->SetBranchAddress("Harm.ActAnScint.hit.tmin", &Harm_ActAnScint_hit_tmin, &b_Harm_ActAnScint_hit_tmin);
   fChain->SetBranchAddress("Harm.ActAnScint.hit.tmax", &Harm_ActAnScint_hit_tmax, &b_Harm_ActAnScint_hit_tmax);
   fChain->SetBranchAddress("Harm.ActAnScint.hit.otridx", &Harm_ActAnScint_hit_otridx, &b_Harm_ActAnScint_hit_otridx);
   fChain->SetBranchAddress("Harm.ActAnScint.hit.ptridx", &Harm_ActAnScint_hit_ptridx, &b_Harm_ActAnScint_hit_ptridx);
   fChain->SetBranchAddress("Harm.ActAnScint.hit.sdtridx", &Harm_ActAnScint_hit_sdtridx, &b_Harm_ActAnScint_hit_sdtridx);
   fChain->SetBranchAddress("Harm.CEPolFront.hit.nhits", &Harm_CEPolFront_hit_nhits, &b_Harm_CEPolFront_hit_nhits);
   fChain->SetBranchAddress("Harm.CEPolFront.hit.plane", &Harm_CEPolFront_hit_plane, &b_Harm_CEPolFront_hit_plane);
   fChain->SetBranchAddress("Harm.CEPolFront.hit.strip", &Harm_CEPolFront_hit_strip, &b_Harm_CEPolFront_hit_strip);
   fChain->SetBranchAddress("Harm.CEPolFront.hit.x", &Harm_CEPolFront_hit_x, &b_Harm_CEPolFront_hit_x);
   fChain->SetBranchAddress("Harm.CEPolFront.hit.y", &Harm_CEPolFront_hit_y, &b_Harm_CEPolFront_hit_y);
   fChain->SetBranchAddress("Harm.CEPolFront.hit.z", &Harm_CEPolFront_hit_z, &b_Harm_CEPolFront_hit_z);
   fChain->SetBranchAddress("Harm.CEPolFront.hit.polx", &Harm_CEPolFront_hit_polx, &b_Harm_CEPolFront_hit_polx);
   fChain->SetBranchAddress("Harm.CEPolFront.hit.poly", &Harm_CEPolFront_hit_poly, &b_Harm_CEPolFront_hit_poly);
   fChain->SetBranchAddress("Harm.CEPolFront.hit.polz", &Harm_CEPolFront_hit_polz, &b_Harm_CEPolFront_hit_polz);
   fChain->SetBranchAddress("Harm.CEPolFront.hit.t", &Harm_CEPolFront_hit_t, &b_Harm_CEPolFront_hit_t);
   fChain->SetBranchAddress("Harm.CEPolFront.hit.trms", &Harm_CEPolFront_hit_trms, &b_Harm_CEPolFront_hit_trms);
   fChain->SetBranchAddress("Harm.CEPolFront.hit.tmin", &Harm_CEPolFront_hit_tmin, &b_Harm_CEPolFront_hit_tmin);
   fChain->SetBranchAddress("Harm.CEPolFront.hit.tmax", &Harm_CEPolFront_hit_tmax, &b_Harm_CEPolFront_hit_tmax);
   fChain->SetBranchAddress("Harm.CEPolFront.hit.tx", &Harm_CEPolFront_hit_tx, &b_Harm_CEPolFront_hit_tx);
   fChain->SetBranchAddress("Harm.CEPolFront.hit.ty", &Harm_CEPolFront_hit_ty, &b_Harm_CEPolFront_hit_ty);
   fChain->SetBranchAddress("Harm.CEPolFront.hit.xin", &Harm_CEPolFront_hit_xin, &b_Harm_CEPolFront_hit_xin);
   fChain->SetBranchAddress("Harm.CEPolFront.hit.yin", &Harm_CEPolFront_hit_yin, &b_Harm_CEPolFront_hit_yin);
   fChain->SetBranchAddress("Harm.CEPolFront.hit.zin", &Harm_CEPolFront_hit_zin, &b_Harm_CEPolFront_hit_zin);
   fChain->SetBranchAddress("Harm.CEPolFront.hit.xout", &Harm_CEPolFront_hit_xout, &b_Harm_CEPolFront_hit_xout);
   fChain->SetBranchAddress("Harm.CEPolFront.hit.yout", &Harm_CEPolFront_hit_yout, &b_Harm_CEPolFront_hit_yout);
   fChain->SetBranchAddress("Harm.CEPolFront.hit.zout", &Harm_CEPolFront_hit_zout, &b_Harm_CEPolFront_hit_zout);
   fChain->SetBranchAddress("Harm.CEPolFront.hit.txp", &Harm_CEPolFront_hit_txp, &b_Harm_CEPolFront_hit_txp);
   fChain->SetBranchAddress("Harm.CEPolFront.hit.typ", &Harm_CEPolFront_hit_typ, &b_Harm_CEPolFront_hit_typ);
   fChain->SetBranchAddress("Harm.CEPolFront.hit.xg", &Harm_CEPolFront_hit_xg, &b_Harm_CEPolFront_hit_xg);
   fChain->SetBranchAddress("Harm.CEPolFront.hit.yg", &Harm_CEPolFront_hit_yg, &b_Harm_CEPolFront_hit_yg);
   fChain->SetBranchAddress("Harm.CEPolFront.hit.zg", &Harm_CEPolFront_hit_zg, &b_Harm_CEPolFront_hit_zg);
   fChain->SetBranchAddress("Harm.CEPolFront.hit.trid", &Harm_CEPolFront_hit_trid, &b_Harm_CEPolFront_hit_trid);
   fChain->SetBranchAddress("Harm.CEPolFront.hit.mid", &Harm_CEPolFront_hit_mid, &b_Harm_CEPolFront_hit_mid);
   fChain->SetBranchAddress("Harm.CEPolFront.hit.pid", &Harm_CEPolFront_hit_pid, &b_Harm_CEPolFront_hit_pid);
   fChain->SetBranchAddress("Harm.CEPolFront.hit.vx", &Harm_CEPolFront_hit_vx, &b_Harm_CEPolFront_hit_vx);
   fChain->SetBranchAddress("Harm.CEPolFront.hit.vy", &Harm_CEPolFront_hit_vy, &b_Harm_CEPolFront_hit_vy);
   fChain->SetBranchAddress("Harm.CEPolFront.hit.vz", &Harm_CEPolFront_hit_vz, &b_Harm_CEPolFront_hit_vz);
   fChain->SetBranchAddress("Harm.CEPolFront.hit.p", &Harm_CEPolFront_hit_p, &b_Harm_CEPolFront_hit_p);
   fChain->SetBranchAddress("Harm.CEPolFront.hit.edep", &Harm_CEPolFront_hit_edep, &b_Harm_CEPolFront_hit_edep);
   fChain->SetBranchAddress("Harm.CEPolFront.hit.beta", &Harm_CEPolFront_hit_beta, &b_Harm_CEPolFront_hit_beta);
   fChain->SetBranchAddress("Harm.CEPolFront.hit.otridx", &Harm_CEPolFront_hit_otridx, &b_Harm_CEPolFront_hit_otridx);
   fChain->SetBranchAddress("Harm.CEPolFront.hit.ptridx", &Harm_CEPolFront_hit_ptridx, &b_Harm_CEPolFront_hit_ptridx);
   fChain->SetBranchAddress("Harm.CEPolFront.hit.sdtridx", &Harm_CEPolFront_hit_sdtridx, &b_Harm_CEPolFront_hit_sdtridx);
   fChain->SetBranchAddress("Harm.CEPolFront.Track.ntracks", &Harm_CEPolFront_Track_ntracks, &b_Harm_CEPolFront_Track_ntracks);
   fChain->SetBranchAddress("Harm.CEPolFront.Track.TID", &Harm_CEPolFront_Track_TID, &b_Harm_CEPolFront_Track_TID);
   fChain->SetBranchAddress("Harm.CEPolFront.Track.PID", &Harm_CEPolFront_Track_PID, &b_Harm_CEPolFront_Track_PID);
   fChain->SetBranchAddress("Harm.CEPolFront.Track.MID", &Harm_CEPolFront_Track_MID, &b_Harm_CEPolFront_Track_MID);
   fChain->SetBranchAddress("Harm.CEPolFront.Track.NumHits", &Harm_CEPolFront_Track_NumHits, &b_Harm_CEPolFront_Track_NumHits);
   fChain->SetBranchAddress("Harm.CEPolFront.Track.NumPlanes", &Harm_CEPolFront_Track_NumPlanes, &b_Harm_CEPolFront_Track_NumPlanes);
   fChain->SetBranchAddress("Harm.CEPolFront.Track.NDF", &Harm_CEPolFront_Track_NDF, &b_Harm_CEPolFront_Track_NDF);
   fChain->SetBranchAddress("Harm.CEPolFront.Track.Chi2fit", &Harm_CEPolFront_Track_Chi2fit, &b_Harm_CEPolFront_Track_Chi2fit);
   fChain->SetBranchAddress("Harm.CEPolFront.Track.Chi2true", &Harm_CEPolFront_Track_Chi2true, &b_Harm_CEPolFront_Track_Chi2true);
   fChain->SetBranchAddress("Harm.CEPolFront.Track.X", &Harm_CEPolFront_Track_X, &b_Harm_CEPolFront_Track_X);
   fChain->SetBranchAddress("Harm.CEPolFront.Track.Y", &Harm_CEPolFront_Track_Y, &b_Harm_CEPolFront_Track_Y);
   fChain->SetBranchAddress("Harm.CEPolFront.Track.Xp", &Harm_CEPolFront_Track_Xp, &b_Harm_CEPolFront_Track_Xp);
   fChain->SetBranchAddress("Harm.CEPolFront.Track.Yp", &Harm_CEPolFront_Track_Yp, &b_Harm_CEPolFront_Track_Yp);
   fChain->SetBranchAddress("Harm.CEPolFront.Track.T", &Harm_CEPolFront_Track_T, &b_Harm_CEPolFront_Track_T);
   fChain->SetBranchAddress("Harm.CEPolFront.Track.P", &Harm_CEPolFront_Track_P, &b_Harm_CEPolFront_Track_P);
   fChain->SetBranchAddress("Harm.CEPolFront.Track.Sx", &Harm_CEPolFront_Track_Sx, &b_Harm_CEPolFront_Track_Sx);
   fChain->SetBranchAddress("Harm.CEPolFront.Track.Sy", &Harm_CEPolFront_Track_Sy, &b_Harm_CEPolFront_Track_Sy);
   fChain->SetBranchAddress("Harm.CEPolFront.Track.Sz", &Harm_CEPolFront_Track_Sz, &b_Harm_CEPolFront_Track_Sz);
   fChain->SetBranchAddress("Harm.CEPolFront.Track.Xfit", &Harm_CEPolFront_Track_Xfit, &b_Harm_CEPolFront_Track_Xfit);
   fChain->SetBranchAddress("Harm.CEPolFront.Track.Yfit", &Harm_CEPolFront_Track_Yfit, &b_Harm_CEPolFront_Track_Yfit);
   fChain->SetBranchAddress("Harm.CEPolFront.Track.Xpfit", &Harm_CEPolFront_Track_Xpfit, &b_Harm_CEPolFront_Track_Xpfit);
   fChain->SetBranchAddress("Harm.CEPolFront.Track.Ypfit", &Harm_CEPolFront_Track_Ypfit, &b_Harm_CEPolFront_Track_Ypfit);
   fChain->SetBranchAddress("Harm.CEPolFront.Track.otridx", &Harm_CEPolFront_Track_otridx, &b_Harm_CEPolFront_Track_otridx);
   fChain->SetBranchAddress("Harm.CEPolFront.Track.ptridx", &Harm_CEPolFront_Track_ptridx, &b_Harm_CEPolFront_Track_ptridx);
   fChain->SetBranchAddress("Harm.CEPolFront.Track.sdtridx", &Harm_CEPolFront_Track_sdtridx, &b_Harm_CEPolFront_Track_sdtridx);
   fChain->SetBranchAddress("Harm.CEPolRear.hit.nhits", &Harm_CEPolRear_hit_nhits, &b_Harm_CEPolRear_hit_nhits);
   fChain->SetBranchAddress("Harm.CEPolRear.hit.plane", &Harm_CEPolRear_hit_plane, &b_Harm_CEPolRear_hit_plane);
   fChain->SetBranchAddress("Harm.CEPolRear.hit.strip", &Harm_CEPolRear_hit_strip, &b_Harm_CEPolRear_hit_strip);
   fChain->SetBranchAddress("Harm.CEPolRear.hit.x", &Harm_CEPolRear_hit_x, &b_Harm_CEPolRear_hit_x);
   fChain->SetBranchAddress("Harm.CEPolRear.hit.y", &Harm_CEPolRear_hit_y, &b_Harm_CEPolRear_hit_y);
   fChain->SetBranchAddress("Harm.CEPolRear.hit.z", &Harm_CEPolRear_hit_z, &b_Harm_CEPolRear_hit_z);
   fChain->SetBranchAddress("Harm.CEPolRear.hit.polx", &Harm_CEPolRear_hit_polx, &b_Harm_CEPolRear_hit_polx);
   fChain->SetBranchAddress("Harm.CEPolRear.hit.poly", &Harm_CEPolRear_hit_poly, &b_Harm_CEPolRear_hit_poly);
   fChain->SetBranchAddress("Harm.CEPolRear.hit.polz", &Harm_CEPolRear_hit_polz, &b_Harm_CEPolRear_hit_polz);
   fChain->SetBranchAddress("Harm.CEPolRear.hit.t", &Harm_CEPolRear_hit_t, &b_Harm_CEPolRear_hit_t);
   fChain->SetBranchAddress("Harm.CEPolRear.hit.trms", &Harm_CEPolRear_hit_trms, &b_Harm_CEPolRear_hit_trms);
   fChain->SetBranchAddress("Harm.CEPolRear.hit.tmin", &Harm_CEPolRear_hit_tmin, &b_Harm_CEPolRear_hit_tmin);
   fChain->SetBranchAddress("Harm.CEPolRear.hit.tmax", &Harm_CEPolRear_hit_tmax, &b_Harm_CEPolRear_hit_tmax);
   fChain->SetBranchAddress("Harm.CEPolRear.hit.tx", &Harm_CEPolRear_hit_tx, &b_Harm_CEPolRear_hit_tx);
   fChain->SetBranchAddress("Harm.CEPolRear.hit.ty", &Harm_CEPolRear_hit_ty, &b_Harm_CEPolRear_hit_ty);
   fChain->SetBranchAddress("Harm.CEPolRear.hit.xin", &Harm_CEPolRear_hit_xin, &b_Harm_CEPolRear_hit_xin);
   fChain->SetBranchAddress("Harm.CEPolRear.hit.yin", &Harm_CEPolRear_hit_yin, &b_Harm_CEPolRear_hit_yin);
   fChain->SetBranchAddress("Harm.CEPolRear.hit.zin", &Harm_CEPolRear_hit_zin, &b_Harm_CEPolRear_hit_zin);
   fChain->SetBranchAddress("Harm.CEPolRear.hit.xout", &Harm_CEPolRear_hit_xout, &b_Harm_CEPolRear_hit_xout);
   fChain->SetBranchAddress("Harm.CEPolRear.hit.yout", &Harm_CEPolRear_hit_yout, &b_Harm_CEPolRear_hit_yout);
   fChain->SetBranchAddress("Harm.CEPolRear.hit.zout", &Harm_CEPolRear_hit_zout, &b_Harm_CEPolRear_hit_zout);
   fChain->SetBranchAddress("Harm.CEPolRear.hit.txp", &Harm_CEPolRear_hit_txp, &b_Harm_CEPolRear_hit_txp);
   fChain->SetBranchAddress("Harm.CEPolRear.hit.typ", &Harm_CEPolRear_hit_typ, &b_Harm_CEPolRear_hit_typ);
   fChain->SetBranchAddress("Harm.CEPolRear.hit.xg", &Harm_CEPolRear_hit_xg, &b_Harm_CEPolRear_hit_xg);
   fChain->SetBranchAddress("Harm.CEPolRear.hit.yg", &Harm_CEPolRear_hit_yg, &b_Harm_CEPolRear_hit_yg);
   fChain->SetBranchAddress("Harm.CEPolRear.hit.zg", &Harm_CEPolRear_hit_zg, &b_Harm_CEPolRear_hit_zg);
   fChain->SetBranchAddress("Harm.CEPolRear.hit.trid", &Harm_CEPolRear_hit_trid, &b_Harm_CEPolRear_hit_trid);
   fChain->SetBranchAddress("Harm.CEPolRear.hit.mid", &Harm_CEPolRear_hit_mid, &b_Harm_CEPolRear_hit_mid);
   fChain->SetBranchAddress("Harm.CEPolRear.hit.pid", &Harm_CEPolRear_hit_pid, &b_Harm_CEPolRear_hit_pid);
   fChain->SetBranchAddress("Harm.CEPolRear.hit.vx", &Harm_CEPolRear_hit_vx, &b_Harm_CEPolRear_hit_vx);
   fChain->SetBranchAddress("Harm.CEPolRear.hit.vy", &Harm_CEPolRear_hit_vy, &b_Harm_CEPolRear_hit_vy);
   fChain->SetBranchAddress("Harm.CEPolRear.hit.vz", &Harm_CEPolRear_hit_vz, &b_Harm_CEPolRear_hit_vz);
   fChain->SetBranchAddress("Harm.CEPolRear.hit.p", &Harm_CEPolRear_hit_p, &b_Harm_CEPolRear_hit_p);
   fChain->SetBranchAddress("Harm.CEPolRear.hit.edep", &Harm_CEPolRear_hit_edep, &b_Harm_CEPolRear_hit_edep);
   fChain->SetBranchAddress("Harm.CEPolRear.hit.beta", &Harm_CEPolRear_hit_beta, &b_Harm_CEPolRear_hit_beta);
   fChain->SetBranchAddress("Harm.CEPolRear.hit.otridx", &Harm_CEPolRear_hit_otridx, &b_Harm_CEPolRear_hit_otridx);
   fChain->SetBranchAddress("Harm.CEPolRear.hit.ptridx", &Harm_CEPolRear_hit_ptridx, &b_Harm_CEPolRear_hit_ptridx);
   fChain->SetBranchAddress("Harm.CEPolRear.hit.sdtridx", &Harm_CEPolRear_hit_sdtridx, &b_Harm_CEPolRear_hit_sdtridx);
   fChain->SetBranchAddress("Harm.CEPolRear.Track.ntracks", &Harm_CEPolRear_Track_ntracks, &b_Harm_CEPolRear_Track_ntracks);
   fChain->SetBranchAddress("Harm.CEPolRear.Track.TID", &Harm_CEPolRear_Track_TID, &b_Harm_CEPolRear_Track_TID);
   fChain->SetBranchAddress("Harm.CEPolRear.Track.PID", &Harm_CEPolRear_Track_PID, &b_Harm_CEPolRear_Track_PID);
   fChain->SetBranchAddress("Harm.CEPolRear.Track.MID", &Harm_CEPolRear_Track_MID, &b_Harm_CEPolRear_Track_MID);
   fChain->SetBranchAddress("Harm.CEPolRear.Track.NumHits", &Harm_CEPolRear_Track_NumHits, &b_Harm_CEPolRear_Track_NumHits);
   fChain->SetBranchAddress("Harm.CEPolRear.Track.NumPlanes", &Harm_CEPolRear_Track_NumPlanes, &b_Harm_CEPolRear_Track_NumPlanes);
   fChain->SetBranchAddress("Harm.CEPolRear.Track.NDF", &Harm_CEPolRear_Track_NDF, &b_Harm_CEPolRear_Track_NDF);
   fChain->SetBranchAddress("Harm.CEPolRear.Track.Chi2fit", &Harm_CEPolRear_Track_Chi2fit, &b_Harm_CEPolRear_Track_Chi2fit);
   fChain->SetBranchAddress("Harm.CEPolRear.Track.Chi2true", &Harm_CEPolRear_Track_Chi2true, &b_Harm_CEPolRear_Track_Chi2true);
   fChain->SetBranchAddress("Harm.CEPolRear.Track.X", &Harm_CEPolRear_Track_X, &b_Harm_CEPolRear_Track_X);
   fChain->SetBranchAddress("Harm.CEPolRear.Track.Y", &Harm_CEPolRear_Track_Y, &b_Harm_CEPolRear_Track_Y);
   fChain->SetBranchAddress("Harm.CEPolRear.Track.Xp", &Harm_CEPolRear_Track_Xp, &b_Harm_CEPolRear_Track_Xp);
   fChain->SetBranchAddress("Harm.CEPolRear.Track.Yp", &Harm_CEPolRear_Track_Yp, &b_Harm_CEPolRear_Track_Yp);
   fChain->SetBranchAddress("Harm.CEPolRear.Track.T", &Harm_CEPolRear_Track_T, &b_Harm_CEPolRear_Track_T);
   fChain->SetBranchAddress("Harm.CEPolRear.Track.P", &Harm_CEPolRear_Track_P, &b_Harm_CEPolRear_Track_P);
   fChain->SetBranchAddress("Harm.CEPolRear.Track.Sx", &Harm_CEPolRear_Track_Sx, &b_Harm_CEPolRear_Track_Sx);
   fChain->SetBranchAddress("Harm.CEPolRear.Track.Sy", &Harm_CEPolRear_Track_Sy, &b_Harm_CEPolRear_Track_Sy);
   fChain->SetBranchAddress("Harm.CEPolRear.Track.Sz", &Harm_CEPolRear_Track_Sz, &b_Harm_CEPolRear_Track_Sz);
   fChain->SetBranchAddress("Harm.CEPolRear.Track.Xfit", &Harm_CEPolRear_Track_Xfit, &b_Harm_CEPolRear_Track_Xfit);
   fChain->SetBranchAddress("Harm.CEPolRear.Track.Yfit", &Harm_CEPolRear_Track_Yfit, &b_Harm_CEPolRear_Track_Yfit);
   fChain->SetBranchAddress("Harm.CEPolRear.Track.Xpfit", &Harm_CEPolRear_Track_Xpfit, &b_Harm_CEPolRear_Track_Xpfit);
   fChain->SetBranchAddress("Harm.CEPolRear.Track.Ypfit", &Harm_CEPolRear_Track_Ypfit, &b_Harm_CEPolRear_Track_Ypfit);
   fChain->SetBranchAddress("Harm.CEPolRear.Track.otridx", &Harm_CEPolRear_Track_otridx, &b_Harm_CEPolRear_Track_otridx);
   fChain->SetBranchAddress("Harm.CEPolRear.Track.ptridx", &Harm_CEPolRear_Track_ptridx, &b_Harm_CEPolRear_Track_ptridx);
   fChain->SetBranchAddress("Harm.CEPolRear.Track.sdtridx", &Harm_CEPolRear_Track_sdtridx, &b_Harm_CEPolRear_Track_sdtridx);
   fChain->SetBranchAddress("Harm.PRPolGEMFarSide.hit.nhits", &Harm_PRPolGEMFarSide_hit_nhits, &b_Harm_PRPolGEMFarSide_hit_nhits);
   fChain->SetBranchAddress("Harm.PRPolGEMFarSide.hit.plane", &Harm_PRPolGEMFarSide_hit_plane, &b_Harm_PRPolGEMFarSide_hit_plane);
   fChain->SetBranchAddress("Harm.PRPolGEMFarSide.hit.strip", &Harm_PRPolGEMFarSide_hit_strip, &b_Harm_PRPolGEMFarSide_hit_strip);
   fChain->SetBranchAddress("Harm.PRPolGEMFarSide.hit.x", &Harm_PRPolGEMFarSide_hit_x, &b_Harm_PRPolGEMFarSide_hit_x);
   fChain->SetBranchAddress("Harm.PRPolGEMFarSide.hit.y", &Harm_PRPolGEMFarSide_hit_y, &b_Harm_PRPolGEMFarSide_hit_y);
   fChain->SetBranchAddress("Harm.PRPolGEMFarSide.hit.z", &Harm_PRPolGEMFarSide_hit_z, &b_Harm_PRPolGEMFarSide_hit_z);
   fChain->SetBranchAddress("Harm.PRPolGEMFarSide.hit.polx", &Harm_PRPolGEMFarSide_hit_polx, &b_Harm_PRPolGEMFarSide_hit_polx);
   fChain->SetBranchAddress("Harm.PRPolGEMFarSide.hit.poly", &Harm_PRPolGEMFarSide_hit_poly, &b_Harm_PRPolGEMFarSide_hit_poly);
   fChain->SetBranchAddress("Harm.PRPolGEMFarSide.hit.polz", &Harm_PRPolGEMFarSide_hit_polz, &b_Harm_PRPolGEMFarSide_hit_polz);
   fChain->SetBranchAddress("Harm.PRPolGEMFarSide.hit.t", &Harm_PRPolGEMFarSide_hit_t, &b_Harm_PRPolGEMFarSide_hit_t);
   fChain->SetBranchAddress("Harm.PRPolGEMFarSide.hit.trms", &Harm_PRPolGEMFarSide_hit_trms, &b_Harm_PRPolGEMFarSide_hit_trms);
   fChain->SetBranchAddress("Harm.PRPolGEMFarSide.hit.tmin", &Harm_PRPolGEMFarSide_hit_tmin, &b_Harm_PRPolGEMFarSide_hit_tmin);
   fChain->SetBranchAddress("Harm.PRPolGEMFarSide.hit.tmax", &Harm_PRPolGEMFarSide_hit_tmax, &b_Harm_PRPolGEMFarSide_hit_tmax);
   fChain->SetBranchAddress("Harm.PRPolGEMFarSide.hit.tx", &Harm_PRPolGEMFarSide_hit_tx, &b_Harm_PRPolGEMFarSide_hit_tx);
   fChain->SetBranchAddress("Harm.PRPolGEMFarSide.hit.ty", &Harm_PRPolGEMFarSide_hit_ty, &b_Harm_PRPolGEMFarSide_hit_ty);
   fChain->SetBranchAddress("Harm.PRPolGEMFarSide.hit.xin", &Harm_PRPolGEMFarSide_hit_xin, &b_Harm_PRPolGEMFarSide_hit_xin);
   fChain->SetBranchAddress("Harm.PRPolGEMFarSide.hit.yin", &Harm_PRPolGEMFarSide_hit_yin, &b_Harm_PRPolGEMFarSide_hit_yin);
   fChain->SetBranchAddress("Harm.PRPolGEMFarSide.hit.zin", &Harm_PRPolGEMFarSide_hit_zin, &b_Harm_PRPolGEMFarSide_hit_zin);
   fChain->SetBranchAddress("Harm.PRPolGEMFarSide.hit.xout", &Harm_PRPolGEMFarSide_hit_xout, &b_Harm_PRPolGEMFarSide_hit_xout);
   fChain->SetBranchAddress("Harm.PRPolGEMFarSide.hit.yout", &Harm_PRPolGEMFarSide_hit_yout, &b_Harm_PRPolGEMFarSide_hit_yout);
   fChain->SetBranchAddress("Harm.PRPolGEMFarSide.hit.zout", &Harm_PRPolGEMFarSide_hit_zout, &b_Harm_PRPolGEMFarSide_hit_zout);
   fChain->SetBranchAddress("Harm.PRPolGEMFarSide.hit.txp", &Harm_PRPolGEMFarSide_hit_txp, &b_Harm_PRPolGEMFarSide_hit_txp);
   fChain->SetBranchAddress("Harm.PRPolGEMFarSide.hit.typ", &Harm_PRPolGEMFarSide_hit_typ, &b_Harm_PRPolGEMFarSide_hit_typ);
   fChain->SetBranchAddress("Harm.PRPolGEMFarSide.hit.xg", &Harm_PRPolGEMFarSide_hit_xg, &b_Harm_PRPolGEMFarSide_hit_xg);
   fChain->SetBranchAddress("Harm.PRPolGEMFarSide.hit.yg", &Harm_PRPolGEMFarSide_hit_yg, &b_Harm_PRPolGEMFarSide_hit_yg);
   fChain->SetBranchAddress("Harm.PRPolGEMFarSide.hit.zg", &Harm_PRPolGEMFarSide_hit_zg, &b_Harm_PRPolGEMFarSide_hit_zg);
   fChain->SetBranchAddress("Harm.PRPolGEMFarSide.hit.trid", &Harm_PRPolGEMFarSide_hit_trid, &b_Harm_PRPolGEMFarSide_hit_trid);
   fChain->SetBranchAddress("Harm.PRPolGEMFarSide.hit.mid", &Harm_PRPolGEMFarSide_hit_mid, &b_Harm_PRPolGEMFarSide_hit_mid);
   fChain->SetBranchAddress("Harm.PRPolGEMFarSide.hit.pid", &Harm_PRPolGEMFarSide_hit_pid, &b_Harm_PRPolGEMFarSide_hit_pid);
   fChain->SetBranchAddress("Harm.PRPolGEMFarSide.hit.vx", &Harm_PRPolGEMFarSide_hit_vx, &b_Harm_PRPolGEMFarSide_hit_vx);
   fChain->SetBranchAddress("Harm.PRPolGEMFarSide.hit.vy", &Harm_PRPolGEMFarSide_hit_vy, &b_Harm_PRPolGEMFarSide_hit_vy);
   fChain->SetBranchAddress("Harm.PRPolGEMFarSide.hit.vz", &Harm_PRPolGEMFarSide_hit_vz, &b_Harm_PRPolGEMFarSide_hit_vz);
   fChain->SetBranchAddress("Harm.PRPolGEMFarSide.hit.p", &Harm_PRPolGEMFarSide_hit_p, &b_Harm_PRPolGEMFarSide_hit_p);
   fChain->SetBranchAddress("Harm.PRPolGEMFarSide.hit.edep", &Harm_PRPolGEMFarSide_hit_edep, &b_Harm_PRPolGEMFarSide_hit_edep);
   fChain->SetBranchAddress("Harm.PRPolGEMFarSide.hit.beta", &Harm_PRPolGEMFarSide_hit_beta, &b_Harm_PRPolGEMFarSide_hit_beta);
   fChain->SetBranchAddress("Harm.PRPolGEMFarSide.hit.otridx", &Harm_PRPolGEMFarSide_hit_otridx, &b_Harm_PRPolGEMFarSide_hit_otridx);
   fChain->SetBranchAddress("Harm.PRPolGEMFarSide.hit.ptridx", &Harm_PRPolGEMFarSide_hit_ptridx, &b_Harm_PRPolGEMFarSide_hit_ptridx);
   fChain->SetBranchAddress("Harm.PRPolGEMFarSide.hit.sdtridx", &Harm_PRPolGEMFarSide_hit_sdtridx, &b_Harm_PRPolGEMFarSide_hit_sdtridx);
   fChain->SetBranchAddress("Harm.PRPolGEMFarSide.Track.ntracks", &Harm_PRPolGEMFarSide_Track_ntracks, &b_Harm_PRPolGEMFarSide_Track_ntracks);
   fChain->SetBranchAddress("Harm.PRPolGEMFarSide.Track.TID", &Harm_PRPolGEMFarSide_Track_TID, &b_Harm_PRPolGEMFarSide_Track_TID);
   fChain->SetBranchAddress("Harm.PRPolGEMFarSide.Track.PID", &Harm_PRPolGEMFarSide_Track_PID, &b_Harm_PRPolGEMFarSide_Track_PID);
   fChain->SetBranchAddress("Harm.PRPolGEMFarSide.Track.MID", &Harm_PRPolGEMFarSide_Track_MID, &b_Harm_PRPolGEMFarSide_Track_MID);
   fChain->SetBranchAddress("Harm.PRPolGEMFarSide.Track.NumHits", &Harm_PRPolGEMFarSide_Track_NumHits, &b_Harm_PRPolGEMFarSide_Track_NumHits);
   fChain->SetBranchAddress("Harm.PRPolGEMFarSide.Track.NumPlanes", &Harm_PRPolGEMFarSide_Track_NumPlanes, &b_Harm_PRPolGEMFarSide_Track_NumPlanes);
   fChain->SetBranchAddress("Harm.PRPolGEMFarSide.Track.NDF", &Harm_PRPolGEMFarSide_Track_NDF, &b_Harm_PRPolGEMFarSide_Track_NDF);
   fChain->SetBranchAddress("Harm.PRPolGEMFarSide.Track.Chi2fit", &Harm_PRPolGEMFarSide_Track_Chi2fit, &b_Harm_PRPolGEMFarSide_Track_Chi2fit);
   fChain->SetBranchAddress("Harm.PRPolGEMFarSide.Track.Chi2true", &Harm_PRPolGEMFarSide_Track_Chi2true, &b_Harm_PRPolGEMFarSide_Track_Chi2true);
   fChain->SetBranchAddress("Harm.PRPolGEMFarSide.Track.X", &Harm_PRPolGEMFarSide_Track_X, &b_Harm_PRPolGEMFarSide_Track_X);
   fChain->SetBranchAddress("Harm.PRPolGEMFarSide.Track.Y", &Harm_PRPolGEMFarSide_Track_Y, &b_Harm_PRPolGEMFarSide_Track_Y);
   fChain->SetBranchAddress("Harm.PRPolGEMFarSide.Track.Xp", &Harm_PRPolGEMFarSide_Track_Xp, &b_Harm_PRPolGEMFarSide_Track_Xp);
   fChain->SetBranchAddress("Harm.PRPolGEMFarSide.Track.Yp", &Harm_PRPolGEMFarSide_Track_Yp, &b_Harm_PRPolGEMFarSide_Track_Yp);
   fChain->SetBranchAddress("Harm.PRPolGEMFarSide.Track.T", &Harm_PRPolGEMFarSide_Track_T, &b_Harm_PRPolGEMFarSide_Track_T);
   fChain->SetBranchAddress("Harm.PRPolGEMFarSide.Track.P", &Harm_PRPolGEMFarSide_Track_P, &b_Harm_PRPolGEMFarSide_Track_P);
   fChain->SetBranchAddress("Harm.PRPolGEMFarSide.Track.Sx", &Harm_PRPolGEMFarSide_Track_Sx, &b_Harm_PRPolGEMFarSide_Track_Sx);
   fChain->SetBranchAddress("Harm.PRPolGEMFarSide.Track.Sy", &Harm_PRPolGEMFarSide_Track_Sy, &b_Harm_PRPolGEMFarSide_Track_Sy);
   fChain->SetBranchAddress("Harm.PRPolGEMFarSide.Track.Sz", &Harm_PRPolGEMFarSide_Track_Sz, &b_Harm_PRPolGEMFarSide_Track_Sz);
   fChain->SetBranchAddress("Harm.PRPolGEMFarSide.Track.Xfit", &Harm_PRPolGEMFarSide_Track_Xfit, &b_Harm_PRPolGEMFarSide_Track_Xfit);
   fChain->SetBranchAddress("Harm.PRPolGEMFarSide.Track.Yfit", &Harm_PRPolGEMFarSide_Track_Yfit, &b_Harm_PRPolGEMFarSide_Track_Yfit);
   fChain->SetBranchAddress("Harm.PRPolGEMFarSide.Track.Xpfit", &Harm_PRPolGEMFarSide_Track_Xpfit, &b_Harm_PRPolGEMFarSide_Track_Xpfit);
   fChain->SetBranchAddress("Harm.PRPolGEMFarSide.Track.Ypfit", &Harm_PRPolGEMFarSide_Track_Ypfit, &b_Harm_PRPolGEMFarSide_Track_Ypfit);
   fChain->SetBranchAddress("Harm.PRPolGEMFarSide.Track.otridx", &Harm_PRPolGEMFarSide_Track_otridx, &b_Harm_PRPolGEMFarSide_Track_otridx);
   fChain->SetBranchAddress("Harm.PRPolGEMFarSide.Track.ptridx", &Harm_PRPolGEMFarSide_Track_ptridx, &b_Harm_PRPolGEMFarSide_Track_ptridx);
   fChain->SetBranchAddress("Harm.PRPolGEMFarSide.Track.sdtridx", &Harm_PRPolGEMFarSide_Track_sdtridx, &b_Harm_PRPolGEMFarSide_Track_sdtridx);
   fChain->SetBranchAddress("Harm.PRPolScintFarSide.det.esum", &Harm_PRPolScintFarSide_det_esum, &b_Harm_PRPolScintFarSide_det_esum);
   fChain->SetBranchAddress("Harm.PRPolScintFarSide.hit.nhits", &Harm_PRPolScintFarSide_hit_nhits, &b_Harm_PRPolScintFarSide_hit_nhits);
   fChain->SetBranchAddress("Harm.PRPolScintFarSide.hit.row", &Harm_PRPolScintFarSide_hit_row, &b_Harm_PRPolScintFarSide_hit_row);
   fChain->SetBranchAddress("Harm.PRPolScintFarSide.hit.col", &Harm_PRPolScintFarSide_hit_col, &b_Harm_PRPolScintFarSide_hit_col);
   fChain->SetBranchAddress("Harm.PRPolScintFarSide.hit.cell", &Harm_PRPolScintFarSide_hit_cell, &b_Harm_PRPolScintFarSide_hit_cell);
   fChain->SetBranchAddress("Harm.PRPolScintFarSide.hit.plane", &Harm_PRPolScintFarSide_hit_plane, &b_Harm_PRPolScintFarSide_hit_plane);
   fChain->SetBranchAddress("Harm.PRPolScintFarSide.hit.wire", &Harm_PRPolScintFarSide_hit_wire, &b_Harm_PRPolScintFarSide_hit_wire);
   fChain->SetBranchAddress("Harm.PRPolScintFarSide.hit.xcell", &Harm_PRPolScintFarSide_hit_xcell, &b_Harm_PRPolScintFarSide_hit_xcell);
   fChain->SetBranchAddress("Harm.PRPolScintFarSide.hit.ycell", &Harm_PRPolScintFarSide_hit_ycell, &b_Harm_PRPolScintFarSide_hit_ycell);
   fChain->SetBranchAddress("Harm.PRPolScintFarSide.hit.zcell", &Harm_PRPolScintFarSide_hit_zcell, &b_Harm_PRPolScintFarSide_hit_zcell);
   fChain->SetBranchAddress("Harm.PRPolScintFarSide.hit.xcellg", &Harm_PRPolScintFarSide_hit_xcellg, &b_Harm_PRPolScintFarSide_hit_xcellg);
   fChain->SetBranchAddress("Harm.PRPolScintFarSide.hit.ycellg", &Harm_PRPolScintFarSide_hit_ycellg, &b_Harm_PRPolScintFarSide_hit_ycellg);
   fChain->SetBranchAddress("Harm.PRPolScintFarSide.hit.zcellg", &Harm_PRPolScintFarSide_hit_zcellg, &b_Harm_PRPolScintFarSide_hit_zcellg);
   fChain->SetBranchAddress("Harm.PRPolScintFarSide.hit.xhit", &Harm_PRPolScintFarSide_hit_xhit, &b_Harm_PRPolScintFarSide_hit_xhit);
   fChain->SetBranchAddress("Harm.PRPolScintFarSide.hit.yhit", &Harm_PRPolScintFarSide_hit_yhit, &b_Harm_PRPolScintFarSide_hit_yhit);
   fChain->SetBranchAddress("Harm.PRPolScintFarSide.hit.zhit", &Harm_PRPolScintFarSide_hit_zhit, &b_Harm_PRPolScintFarSide_hit_zhit);
   fChain->SetBranchAddress("Harm.PRPolScintFarSide.hit.xhitg", &Harm_PRPolScintFarSide_hit_xhitg, &b_Harm_PRPolScintFarSide_hit_xhitg);
   fChain->SetBranchAddress("Harm.PRPolScintFarSide.hit.yhitg", &Harm_PRPolScintFarSide_hit_yhitg, &b_Harm_PRPolScintFarSide_hit_yhitg);
   fChain->SetBranchAddress("Harm.PRPolScintFarSide.hit.zhitg", &Harm_PRPolScintFarSide_hit_zhitg, &b_Harm_PRPolScintFarSide_hit_zhitg);
   fChain->SetBranchAddress("Harm.PRPolScintFarSide.hit.sumedep", &Harm_PRPolScintFarSide_hit_sumedep, &b_Harm_PRPolScintFarSide_hit_sumedep);
   fChain->SetBranchAddress("Harm.PRPolScintFarSide.hit.tavg", &Harm_PRPolScintFarSide_hit_tavg, &b_Harm_PRPolScintFarSide_hit_tavg);
   fChain->SetBranchAddress("Harm.PRPolScintFarSide.hit.trms", &Harm_PRPolScintFarSide_hit_trms, &b_Harm_PRPolScintFarSide_hit_trms);
   fChain->SetBranchAddress("Harm.PRPolScintFarSide.hit.tmin", &Harm_PRPolScintFarSide_hit_tmin, &b_Harm_PRPolScintFarSide_hit_tmin);
   fChain->SetBranchAddress("Harm.PRPolScintFarSide.hit.tmax", &Harm_PRPolScintFarSide_hit_tmax, &b_Harm_PRPolScintFarSide_hit_tmax);
   fChain->SetBranchAddress("Harm.PRPolScintFarSide.hit.otridx", &Harm_PRPolScintFarSide_hit_otridx, &b_Harm_PRPolScintFarSide_hit_otridx);
   fChain->SetBranchAddress("Harm.PRPolScintFarSide.hit.ptridx", &Harm_PRPolScintFarSide_hit_ptridx, &b_Harm_PRPolScintFarSide_hit_ptridx);
   fChain->SetBranchAddress("Harm.PRPolScintFarSide.hit.sdtridx", &Harm_PRPolScintFarSide_hit_sdtridx, &b_Harm_PRPolScintFarSide_hit_sdtridx);
   fChain->SetBranchAddress("Harm.ActAn.dighit.nchan", &Harm_ActAn_dighit_nchan, &b_Harm_ActAn_dighit_nchan);
   fChain->SetBranchAddress("Harm.ActAn.dighit.chan", &Harm_ActAn_dighit_chan, &b_Harm_ActAn_dighit_chan);
   fChain->SetBranchAddress("Harm.ActAn.dighit.adc", &Harm_ActAn_dighit_adc, &b_Harm_ActAn_dighit_adc);
   fChain->SetBranchAddress("Harm.ActAn.dighit.samp", &Harm_ActAn_dighit_samp, &b_Harm_ActAn_dighit_samp);
   fChain->SetBranchAddress("Harm.ActAn.dighit.tdc", &Harm_ActAn_dighit_tdc, &b_Harm_ActAn_dighit_tdc);
   fChain->SetBranchAddress("Harm.PRPolScintFarSide.dighit.nchan", &Harm_PRPolScintFarSide_dighit_nchan, &b_Harm_PRPolScintFarSide_dighit_nchan);
   fChain->SetBranchAddress("Harm.PRPolScintFarSide.dighit.chan", &Harm_PRPolScintFarSide_dighit_chan, &b_Harm_PRPolScintFarSide_dighit_chan);
   fChain->SetBranchAddress("Harm.PRPolScintFarSide.dighit.adc", &Harm_PRPolScintFarSide_dighit_adc, &b_Harm_PRPolScintFarSide_dighit_adc);
   fChain->SetBranchAddress("Harm.PRPolScintFarSide.dighit.tdc_l", &Harm_PRPolScintFarSide_dighit_tdc_l, &b_Harm_PRPolScintFarSide_dighit_tdc_l);
   fChain->SetBranchAddress("Harm.PRPolScintFarSide.dighit.tdc_t", &Harm_PRPolScintFarSide_dighit_tdc_t, &b_Harm_PRPolScintFarSide_dighit_tdc_t);
   fChain->SetBranchAddress("Harm.CEPolFront.dighit.nstrips", &Harm_CEPolFront_dighit_nstrips, &b_Harm_CEPolFront_dighit_nstrips);
   fChain->SetBranchAddress("Harm.CEPolFront.dighit.module", &Harm_CEPolFront_dighit_module, &b_Harm_CEPolFront_dighit_module);
   fChain->SetBranchAddress("Harm.CEPolFront.dighit.strip", &Harm_CEPolFront_dighit_strip, &b_Harm_CEPolFront_dighit_strip);
   fChain->SetBranchAddress("Harm.CEPolFront.dighit.adc", &Harm_CEPolFront_dighit_adc, &b_Harm_CEPolFront_dighit_adc);
   fChain->SetBranchAddress("Harm.CEPolFront.dighit.samp", &Harm_CEPolFront_dighit_samp, &b_Harm_CEPolFront_dighit_samp);
   fChain->SetBranchAddress("Harm.CEPolRear.dighit.nstrips", &Harm_CEPolRear_dighit_nstrips, &b_Harm_CEPolRear_dighit_nstrips);
   fChain->SetBranchAddress("Harm.CEPolRear.dighit.module", &Harm_CEPolRear_dighit_module, &b_Harm_CEPolRear_dighit_module);
   fChain->SetBranchAddress("Harm.CEPolRear.dighit.strip", &Harm_CEPolRear_dighit_strip, &b_Harm_CEPolRear_dighit_strip);
   fChain->SetBranchAddress("Harm.CEPolRear.dighit.adc", &Harm_CEPolRear_dighit_adc, &b_Harm_CEPolRear_dighit_adc);
   fChain->SetBranchAddress("Harm.CEPolRear.dighit.samp", &Harm_CEPolRear_dighit_samp, &b_Harm_CEPolRear_dighit_samp);
   fChain->SetBranchAddress("Harm.PRPolGEMFarSide.dighit.nstrips", &Harm_PRPolGEMFarSide_dighit_nstrips, &b_Harm_PRPolGEMFarSide_dighit_nstrips);
   fChain->SetBranchAddress("Harm.PRPolGEMFarSide.dighit.module", &Harm_PRPolGEMFarSide_dighit_module, &b_Harm_PRPolGEMFarSide_dighit_module);
   fChain->SetBranchAddress("Harm.PRPolGEMFarSide.dighit.strip", &Harm_PRPolGEMFarSide_dighit_strip, &b_Harm_PRPolGEMFarSide_dighit_strip);
   fChain->SetBranchAddress("Harm.PRPolGEMFarSide.dighit.adc", &Harm_PRPolGEMFarSide_dighit_adc, &b_Harm_PRPolGEMFarSide_dighit_adc);
   fChain->SetBranchAddress("Harm.PRPolGEMFarSide.dighit.samp", &Harm_PRPolGEMFarSide_dighit_samp, &b_Harm_PRPolGEMFarSide_dighit_samp);
   Notify();
}

Bool_t genrp_tree_digitized::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void genrp_tree_digitized::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t genrp_tree_digitized::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef genrp_tree_digitized_cxx
