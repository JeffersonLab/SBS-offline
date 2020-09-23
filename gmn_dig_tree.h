//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Sep 22 22:47:14 2020 by ROOT version 6.14/04
// from TTree T/Geant4 SBS Simulation
// found on file: digitized/simdigtest_0.root
//////////////////////////////////////////////////////////

#ifndef gmn_dig_tree_h
#define gmn_dig_tree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"

#define treeName "T"

using namespace std;

class gmn_dig_tree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
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
   Int_t           Earm_BBGEM_hit_nhits;
   vector<int>     *Earm_BBGEM_hit_plane;
   vector<int>     *Earm_BBGEM_hit_strip;
   vector<double>  *Earm_BBGEM_hit_x;
   vector<double>  *Earm_BBGEM_hit_y;
   vector<double>  *Earm_BBGEM_hit_z;
   vector<double>  *Earm_BBGEM_hit_polx;
   vector<double>  *Earm_BBGEM_hit_poly;
   vector<double>  *Earm_BBGEM_hit_polz;
   vector<double>  *Earm_BBGEM_hit_t;
   vector<double>  *Earm_BBGEM_hit_trms;
   vector<double>  *Earm_BBGEM_hit_tmin;
   vector<double>  *Earm_BBGEM_hit_tmax;
   vector<double>  *Earm_BBGEM_hit_tx;
   vector<double>  *Earm_BBGEM_hit_ty;
   vector<double>  *Earm_BBGEM_hit_xin;
   vector<double>  *Earm_BBGEM_hit_yin;
   vector<double>  *Earm_BBGEM_hit_zin;
   vector<double>  *Earm_BBGEM_hit_xout;
   vector<double>  *Earm_BBGEM_hit_yout;
   vector<double>  *Earm_BBGEM_hit_zout;
   vector<double>  *Earm_BBGEM_hit_txp;
   vector<double>  *Earm_BBGEM_hit_typ;
   vector<double>  *Earm_BBGEM_hit_xg;
   vector<double>  *Earm_BBGEM_hit_yg;
   vector<double>  *Earm_BBGEM_hit_zg;
   vector<int>     *Earm_BBGEM_hit_trid;
   vector<int>     *Earm_BBGEM_hit_mid;
   vector<int>     *Earm_BBGEM_hit_pid;
   vector<double>  *Earm_BBGEM_hit_vx;
   vector<double>  *Earm_BBGEM_hit_vy;
   vector<double>  *Earm_BBGEM_hit_vz;
   vector<double>  *Earm_BBGEM_hit_p;
   vector<double>  *Earm_BBGEM_hit_edep;
   vector<double>  *Earm_BBGEM_hit_beta;
   vector<int>     *Earm_BBGEM_hit_otridx;
   vector<int>     *Earm_BBGEM_hit_ptridx;
   vector<int>     *Earm_BBGEM_hit_sdtridx;
   Int_t           Earm_BBGEM_Track_ntracks;
   vector<int>     *Earm_BBGEM_Track_TID;
   vector<int>     *Earm_BBGEM_Track_PID;
   vector<int>     *Earm_BBGEM_Track_MID;
   vector<int>     *Earm_BBGEM_Track_NumHits;
   vector<int>     *Earm_BBGEM_Track_NumPlanes;
   vector<int>     *Earm_BBGEM_Track_NDF;
   vector<double>  *Earm_BBGEM_Track_Chi2fit;
   vector<double>  *Earm_BBGEM_Track_Chi2true;
   vector<double>  *Earm_BBGEM_Track_X;
   vector<double>  *Earm_BBGEM_Track_Y;
   vector<double>  *Earm_BBGEM_Track_Xp;
   vector<double>  *Earm_BBGEM_Track_Yp;
   vector<double>  *Earm_BBGEM_Track_T;
   vector<double>  *Earm_BBGEM_Track_P;
   vector<double>  *Earm_BBGEM_Track_Sx;
   vector<double>  *Earm_BBGEM_Track_Sy;
   vector<double>  *Earm_BBGEM_Track_Sz;
   vector<double>  *Earm_BBGEM_Track_Xfit;
   vector<double>  *Earm_BBGEM_Track_Yfit;
   vector<double>  *Earm_BBGEM_Track_Xpfit;
   vector<double>  *Earm_BBGEM_Track_Ypfit;
   vector<int>     *Earm_BBGEM_Track_otridx;
   vector<int>     *Earm_BBGEM_Track_ptridx;
   vector<int>     *Earm_BBGEM_Track_sdtridx;
   Double_t        Earm_BBHodoScint_det_esum;
   Int_t           Earm_BBHodoScint_hit_nhits;
   vector<int>     *Earm_BBHodoScint_hit_row;
   vector<int>     *Earm_BBHodoScint_hit_col;
   vector<int>     *Earm_BBHodoScint_hit_cell;
   vector<int>     *Earm_BBHodoScint_hit_plane;
   vector<int>     *Earm_BBHodoScint_hit_wire;
   vector<double>  *Earm_BBHodoScint_hit_xcell;
   vector<double>  *Earm_BBHodoScint_hit_ycell;
   vector<double>  *Earm_BBHodoScint_hit_zcell;
   vector<double>  *Earm_BBHodoScint_hit_xcellg;
   vector<double>  *Earm_BBHodoScint_hit_ycellg;
   vector<double>  *Earm_BBHodoScint_hit_zcellg;
   vector<double>  *Earm_BBHodoScint_hit_xhit;
   vector<double>  *Earm_BBHodoScint_hit_yhit;
   vector<double>  *Earm_BBHodoScint_hit_zhit;
   vector<double>  *Earm_BBHodoScint_hit_xhitg;
   vector<double>  *Earm_BBHodoScint_hit_yhitg;
   vector<double>  *Earm_BBHodoScint_hit_zhitg;
   vector<double>  *Earm_BBHodoScint_hit_sumedep;
   vector<double>  *Earm_BBHodoScint_hit_tavg;
   vector<double>  *Earm_BBHodoScint_hit_trms;
   vector<double>  *Earm_BBHodoScint_hit_tmin;
   vector<double>  *Earm_BBHodoScint_hit_tmax;
   vector<int>     *Earm_BBHodoScint_hit_otridx;
   vector<int>     *Earm_BBHodoScint_hit_ptridx;
   vector<int>     *Earm_BBHodoScint_hit_sdtridx;
   Int_t           Earm_BBPS_hit_nhits;
   vector<int>     *Earm_BBPS_hit_PMT;
   vector<int>     *Earm_BBPS_hit_row;
   vector<int>     *Earm_BBPS_hit_col;
   vector<int>     *Earm_BBPS_hit_plane;
   vector<double>  *Earm_BBPS_hit_xcell;
   vector<double>  *Earm_BBPS_hit_ycell;
   vector<double>  *Earm_BBPS_hit_zcell;
   vector<double>  *Earm_BBPS_hit_xgcell;
   vector<double>  *Earm_BBPS_hit_ygcell;
   vector<double>  *Earm_BBPS_hit_zgcell;
   vector<int>     *Earm_BBPS_hit_NumPhotoelectrons;
   vector<double>  *Earm_BBPS_hit_Time_avg;
   vector<double>  *Earm_BBPS_hit_Time_rms;
   vector<double>  *Earm_BBPS_hit_Time_min;
   vector<double>  *Earm_BBPS_hit_Time_max;
   vector<int>     *Earm_BBPS_hit_otridx;
   vector<int>     *Earm_BBPS_hit_ptridx;
   vector<int>     *Earm_BBPS_hit_sdtridx;
   Double_t        Earm_BBPSTF1_det_esum;
   Int_t           Earm_BBPSTF1_hit_nhits;
   vector<int>     *Earm_BBPSTF1_hit_row;
   vector<int>     *Earm_BBPSTF1_hit_col;
   vector<int>     *Earm_BBPSTF1_hit_cell;
   vector<int>     *Earm_BBPSTF1_hit_plane;
   vector<int>     *Earm_BBPSTF1_hit_wire;
   vector<double>  *Earm_BBPSTF1_hit_xcell;
   vector<double>  *Earm_BBPSTF1_hit_ycell;
   vector<double>  *Earm_BBPSTF1_hit_zcell;
   vector<double>  *Earm_BBPSTF1_hit_xcellg;
   vector<double>  *Earm_BBPSTF1_hit_ycellg;
   vector<double>  *Earm_BBPSTF1_hit_zcellg;
   vector<double>  *Earm_BBPSTF1_hit_xhit;
   vector<double>  *Earm_BBPSTF1_hit_yhit;
   vector<double>  *Earm_BBPSTF1_hit_zhit;
   vector<double>  *Earm_BBPSTF1_hit_xhitg;
   vector<double>  *Earm_BBPSTF1_hit_yhitg;
   vector<double>  *Earm_BBPSTF1_hit_zhitg;
   vector<double>  *Earm_BBPSTF1_hit_sumedep;
   vector<double>  *Earm_BBPSTF1_hit_tavg;
   vector<double>  *Earm_BBPSTF1_hit_trms;
   vector<double>  *Earm_BBPSTF1_hit_tmin;
   vector<double>  *Earm_BBPSTF1_hit_tmax;
   vector<int>     *Earm_BBPSTF1_hit_otridx;
   vector<int>     *Earm_BBPSTF1_hit_ptridx;
   vector<int>     *Earm_BBPSTF1_hit_sdtridx;
   Int_t           Earm_BBSH_hit_nhits;
   vector<int>     *Earm_BBSH_hit_PMT;
   vector<int>     *Earm_BBSH_hit_row;
   vector<int>     *Earm_BBSH_hit_col;
   vector<int>     *Earm_BBSH_hit_plane;
   vector<double>  *Earm_BBSH_hit_xcell;
   vector<double>  *Earm_BBSH_hit_ycell;
   vector<double>  *Earm_BBSH_hit_zcell;
   vector<double>  *Earm_BBSH_hit_xgcell;
   vector<double>  *Earm_BBSH_hit_ygcell;
   vector<double>  *Earm_BBSH_hit_zgcell;
   vector<int>     *Earm_BBSH_hit_NumPhotoelectrons;
   vector<double>  *Earm_BBSH_hit_Time_avg;
   vector<double>  *Earm_BBSH_hit_Time_rms;
   vector<double>  *Earm_BBSH_hit_Time_min;
   vector<double>  *Earm_BBSH_hit_Time_max;
   vector<int>     *Earm_BBSH_hit_otridx;
   vector<int>     *Earm_BBSH_hit_ptridx;
   vector<int>     *Earm_BBSH_hit_sdtridx;
   Double_t        Earm_BBSHTF1_det_esum;
   Int_t           Earm_BBSHTF1_hit_nhits;
   vector<int>     *Earm_BBSHTF1_hit_row;
   vector<int>     *Earm_BBSHTF1_hit_col;
   vector<int>     *Earm_BBSHTF1_hit_cell;
   vector<int>     *Earm_BBSHTF1_hit_plane;
   vector<int>     *Earm_BBSHTF1_hit_wire;
   vector<double>  *Earm_BBSHTF1_hit_xcell;
   vector<double>  *Earm_BBSHTF1_hit_ycell;
   vector<double>  *Earm_BBSHTF1_hit_zcell;
   vector<double>  *Earm_BBSHTF1_hit_xcellg;
   vector<double>  *Earm_BBSHTF1_hit_ycellg;
   vector<double>  *Earm_BBSHTF1_hit_zcellg;
   vector<double>  *Earm_BBSHTF1_hit_xhit;
   vector<double>  *Earm_BBSHTF1_hit_yhit;
   vector<double>  *Earm_BBSHTF1_hit_zhit;
   vector<double>  *Earm_BBSHTF1_hit_xhitg;
   vector<double>  *Earm_BBSHTF1_hit_yhitg;
   vector<double>  *Earm_BBSHTF1_hit_zhitg;
   vector<double>  *Earm_BBSHTF1_hit_sumedep;
   vector<double>  *Earm_BBSHTF1_hit_tavg;
   vector<double>  *Earm_BBSHTF1_hit_trms;
   vector<double>  *Earm_BBSHTF1_hit_tmin;
   vector<double>  *Earm_BBSHTF1_hit_tmax;
   vector<int>     *Earm_BBSHTF1_hit_otridx;
   vector<int>     *Earm_BBSHTF1_hit_ptridx;
   vector<int>     *Earm_BBSHTF1_hit_sdtridx;
   Int_t           Earm_GRINCH_hit_nhits;
   vector<int>     *Earm_GRINCH_hit_PMT;
   vector<int>     *Earm_GRINCH_hit_row;
   vector<int>     *Earm_GRINCH_hit_col;
   vector<double>  *Earm_GRINCH_hit_xpmt;
   vector<double>  *Earm_GRINCH_hit_ypmt;
   vector<double>  *Earm_GRINCH_hit_zpmt;
   vector<double>  *Earm_GRINCH_hit_xgpmt;
   vector<double>  *Earm_GRINCH_hit_ygpmt;
   vector<double>  *Earm_GRINCH_hit_zgpmt;
   vector<int>     *Earm_GRINCH_hit_NumPhotoelectrons;
   vector<double>  *Earm_GRINCH_hit_Time_avg;
   vector<double>  *Earm_GRINCH_hit_Time_rms;
   vector<double>  *Earm_GRINCH_hit_Time_min;
   vector<double>  *Earm_GRINCH_hit_Time_max;
   vector<int>     *Earm_GRINCH_hit_mTrackNo;
   vector<double>  *Earm_GRINCH_hit_xhit;
   vector<double>  *Earm_GRINCH_hit_yhit;
   vector<double>  *Earm_GRINCH_hit_zhit;
   vector<double>  *Earm_GRINCH_hit_pxhit;
   vector<double>  *Earm_GRINCH_hit_pyhit;
   vector<double>  *Earm_GRINCH_hit_pzhit;
   vector<double>  *Earm_GRINCH_hit_pvx;
   vector<double>  *Earm_GRINCH_hit_pvy;
   vector<double>  *Earm_GRINCH_hit_pvz;
   vector<double>  *Earm_GRINCH_hit_ppx;
   vector<double>  *Earm_GRINCH_hit_ppy;
   vector<double>  *Earm_GRINCH_hit_ppz;
   vector<int>     *Earm_GRINCH_hit_volume_flag;
   vector<int>     *Earm_GRINCH_hit_otridx;
   vector<int>     *Earm_GRINCH_hit_ptridx;
   vector<int>     *Earm_GRINCH_hit_sdtridx;
   Int_t           Harm_HCal_hit_nhits;
   vector<int>     *Harm_HCal_hit_PMT;
   vector<int>     *Harm_HCal_hit_row;
   vector<int>     *Harm_HCal_hit_col;
   vector<int>     *Harm_HCal_hit_plane;
   vector<double>  *Harm_HCal_hit_xcell;
   vector<double>  *Harm_HCal_hit_ycell;
   vector<double>  *Harm_HCal_hit_zcell;
   vector<double>  *Harm_HCal_hit_xgcell;
   vector<double>  *Harm_HCal_hit_ygcell;
   vector<double>  *Harm_HCal_hit_zgcell;
   vector<int>     *Harm_HCal_hit_NumPhotoelectrons;
   vector<double>  *Harm_HCal_hit_Time_avg;
   vector<double>  *Harm_HCal_hit_Time_rms;
   vector<double>  *Harm_HCal_hit_Time_min;
   vector<double>  *Harm_HCal_hit_Time_max;
   vector<int>     *Harm_HCal_hit_otridx;
   vector<int>     *Harm_HCal_hit_ptridx;
   vector<int>     *Harm_HCal_hit_sdtridx;
   Double_t        Harm_HCalScint_det_esum;
   Int_t           Harm_HCalScint_hit_nhits;
   vector<int>     *Harm_HCalScint_hit_row;
   vector<int>     *Harm_HCalScint_hit_col;
   vector<int>     *Harm_HCalScint_hit_cell;
   vector<int>     *Harm_HCalScint_hit_plane;
   vector<int>     *Harm_HCalScint_hit_wire;
   vector<double>  *Harm_HCalScint_hit_xcell;
   vector<double>  *Harm_HCalScint_hit_ycell;
   vector<double>  *Harm_HCalScint_hit_zcell;
   vector<double>  *Harm_HCalScint_hit_xcellg;
   vector<double>  *Harm_HCalScint_hit_ycellg;
   vector<double>  *Harm_HCalScint_hit_zcellg;
   vector<double>  *Harm_HCalScint_hit_xhit;
   vector<double>  *Harm_HCalScint_hit_yhit;
   vector<double>  *Harm_HCalScint_hit_zhit;
   vector<double>  *Harm_HCalScint_hit_xhitg;
   vector<double>  *Harm_HCalScint_hit_yhitg;
   vector<double>  *Harm_HCalScint_hit_zhitg;
   vector<double>  *Harm_HCalScint_hit_sumedep;
   vector<double>  *Harm_HCalScint_hit_tavg;
   vector<double>  *Harm_HCalScint_hit_trms;
   vector<double>  *Harm_HCalScint_hit_tmin;
   vector<double>  *Harm_HCalScint_hit_tmax;
   vector<int>     *Harm_HCalScint_hit_otridx;
   vector<int>     *Harm_HCalScint_hit_ptridx;
   vector<int>     *Harm_HCalScint_hit_sdtridx;
   Int_t           OTrack_ntracks;
   vector<int>     *OTrack_TID;
   vector<int>     *OTrack_MID;
   vector<int>     *OTrack_PID;
   vector<double>  *OTrack_posx;
   vector<double>  *OTrack_posy;
   vector<double>  *OTrack_posz;
   vector<double>  *OTrack_momx;
   vector<double>  *OTrack_momy;
   vector<double>  *OTrack_momz;
   vector<double>  *OTrack_polx;
   vector<double>  *OTrack_poly;
   vector<double>  *OTrack_polz;
   vector<double>  *OTrack_Etot;
   vector<double>  *OTrack_T;
   Int_t           PTrack_ntracks;
   vector<int>     *PTrack_TID;
   vector<int>     *PTrack_PID;
   vector<double>  *PTrack_posx;
   vector<double>  *PTrack_posy;
   vector<double>  *PTrack_posz;
   vector<double>  *PTrack_momx;
   vector<double>  *PTrack_momy;
   vector<double>  *PTrack_momz;
   vector<double>  *PTrack_polx;
   vector<double>  *PTrack_poly;
   vector<double>  *PTrack_polz;
   vector<double>  *PTrack_Etot;
   vector<double>  *PTrack_T;
   Int_t           SDTrack_ntracks;
   vector<int>     *SDTrack_TID;
   vector<int>     *SDTrack_MID;
   vector<int>     *SDTrack_PID;
   vector<double>  *SDTrack_posx;
   vector<double>  *SDTrack_posy;
   vector<double>  *SDTrack_posz;
   vector<double>  *SDTrack_momx;
   vector<double>  *SDTrack_momy;
   vector<double>  *SDTrack_momz;
   vector<double>  *SDTrack_polx;
   vector<double>  *SDTrack_poly;
   vector<double>  *SDTrack_polz;
   vector<double>  *SDTrack_Etot;
   vector<double>  *SDTrack_T;
   Int_t           Earm_BBPS_dighit_nchan;
   vector<int>     *Earm_BBPS_dighit_chan;
   vector<int>     *Earm_BBPS_dighit_adc;
   Int_t           Earm_BBSH_dighit_nchan;
   vector<int>     *Earm_BBSH_dighit_chan;
   vector<int>     *Earm_BBSH_dighit_adc;
   Int_t           Earm_BBHodo_dighit_nchan;
   vector<int>     *Earm_BBHodo_dighit_chan;
   vector<int>     *Earm_BBHodo_dighit_adc;
   vector<int>     *Earm_BBHodo_dighit_tdc_l;
   vector<int>     *Earm_BBHodo_dighit_tdc_t;
   Int_t           Earm_GRINCH_dighit_nchan;
   vector<int>     *Earm_GRINCH_dighit_chan;
   vector<int>     *Earm_GRINCH_dighit_adc;
   vector<int>     *Earm_GRINCH_dighit_tdc_l;
   vector<int>     *Earm_GRINCH_dighit_tdc_t;
   Int_t           Earm_BBGEM_1x_dighit_nstrips;
   vector<int>     *Earm_BBGEM_1x_dighit_strip;
   vector<int>     *Earm_BBGEM_1x_dighit_adc_0;
   vector<int>     *Earm_BBGEM_1x_dighit_adc_1;
   vector<int>     *Earm_BBGEM_1x_dighit_adc_2;
   vector<int>     *Earm_BBGEM_1x_dighit_adc_3;
   vector<int>     *Earm_BBGEM_1x_dighit_adc_4;
   vector<int>     *Earm_BBGEM_1x_dighit_adc_5;
   Int_t           Earm_BBGEM_1y_dighit_nstrips;
   vector<int>     *Earm_BBGEM_1y_dighit_strip;
   vector<int>     *Earm_BBGEM_1y_dighit_adc_0;
   vector<int>     *Earm_BBGEM_1y_dighit_adc_1;
   vector<int>     *Earm_BBGEM_1y_dighit_adc_2;
   vector<int>     *Earm_BBGEM_1y_dighit_adc_3;
   vector<int>     *Earm_BBGEM_1y_dighit_adc_4;
   vector<int>     *Earm_BBGEM_1y_dighit_adc_5;
   Int_t           Earm_BBGEM_2x_dighit_nstrips;
   vector<int>     *Earm_BBGEM_2x_dighit_strip;
   vector<int>     *Earm_BBGEM_2x_dighit_adc_0;
   vector<int>     *Earm_BBGEM_2x_dighit_adc_1;
   vector<int>     *Earm_BBGEM_2x_dighit_adc_2;
   vector<int>     *Earm_BBGEM_2x_dighit_adc_3;
   vector<int>     *Earm_BBGEM_2x_dighit_adc_4;
   vector<int>     *Earm_BBGEM_2x_dighit_adc_5;
   Int_t           Earm_BBGEM_2y_dighit_nstrips;
   vector<int>     *Earm_BBGEM_2y_dighit_strip;
   vector<int>     *Earm_BBGEM_2y_dighit_adc_0;
   vector<int>     *Earm_BBGEM_2y_dighit_adc_1;
   vector<int>     *Earm_BBGEM_2y_dighit_adc_2;
   vector<int>     *Earm_BBGEM_2y_dighit_adc_3;
   vector<int>     *Earm_BBGEM_2y_dighit_adc_4;
   vector<int>     *Earm_BBGEM_2y_dighit_adc_5;
   Int_t           Earm_BBGEM_3x_dighit_nstrips;
   vector<int>     *Earm_BBGEM_3x_dighit_strip;
   vector<int>     *Earm_BBGEM_3x_dighit_adc_0;
   vector<int>     *Earm_BBGEM_3x_dighit_adc_1;
   vector<int>     *Earm_BBGEM_3x_dighit_adc_2;
   vector<int>     *Earm_BBGEM_3x_dighit_adc_3;
   vector<int>     *Earm_BBGEM_3x_dighit_adc_4;
   vector<int>     *Earm_BBGEM_3x_dighit_adc_5;
   Int_t           Earm_BBGEM_3y_dighit_nstrips;
   vector<int>     *Earm_BBGEM_3y_dighit_strip;
   vector<int>     *Earm_BBGEM_3y_dighit_adc_0;
   vector<int>     *Earm_BBGEM_3y_dighit_adc_1;
   vector<int>     *Earm_BBGEM_3y_dighit_adc_2;
   vector<int>     *Earm_BBGEM_3y_dighit_adc_3;
   vector<int>     *Earm_BBGEM_3y_dighit_adc_4;
   vector<int>     *Earm_BBGEM_3y_dighit_adc_5;
   Int_t           Earm_BBGEM_4x_dighit_nstrips;
   vector<int>     *Earm_BBGEM_4x_dighit_strip;
   vector<int>     *Earm_BBGEM_4x_dighit_adc_0;
   vector<int>     *Earm_BBGEM_4x_dighit_adc_1;
   vector<int>     *Earm_BBGEM_4x_dighit_adc_2;
   vector<int>     *Earm_BBGEM_4x_dighit_adc_3;
   vector<int>     *Earm_BBGEM_4x_dighit_adc_4;
   vector<int>     *Earm_BBGEM_4x_dighit_adc_5;
   Int_t           Earm_BBGEM_4y_dighit_nstrips;
   vector<int>     *Earm_BBGEM_4y_dighit_strip;
   vector<int>     *Earm_BBGEM_4y_dighit_adc_0;
   vector<int>     *Earm_BBGEM_4y_dighit_adc_1;
   vector<int>     *Earm_BBGEM_4y_dighit_adc_2;
   vector<int>     *Earm_BBGEM_4y_dighit_adc_3;
   vector<int>     *Earm_BBGEM_4y_dighit_adc_4;
   vector<int>     *Earm_BBGEM_4y_dighit_adc_5;
   Int_t           Earm_BBGEM_5x_dighit_nstrips;
   vector<int>     *Earm_BBGEM_5x_dighit_strip;
   vector<int>     *Earm_BBGEM_5x_dighit_adc_0;
   vector<int>     *Earm_BBGEM_5x_dighit_adc_1;
   vector<int>     *Earm_BBGEM_5x_dighit_adc_2;
   vector<int>     *Earm_BBGEM_5x_dighit_adc_3;
   vector<int>     *Earm_BBGEM_5x_dighit_adc_4;
   vector<int>     *Earm_BBGEM_5x_dighit_adc_5;
   Int_t           Earm_BBGEM_5y_dighit_nstrips;
   vector<int>     *Earm_BBGEM_5y_dighit_strip;
   vector<int>     *Earm_BBGEM_5y_dighit_adc_0;
   vector<int>     *Earm_BBGEM_5y_dighit_adc_1;
   vector<int>     *Earm_BBGEM_5y_dighit_adc_2;
   vector<int>     *Earm_BBGEM_5y_dighit_adc_3;
   vector<int>     *Earm_BBGEM_5y_dighit_adc_4;
   vector<int>     *Earm_BBGEM_5y_dighit_adc_5;
   Int_t           Earm_BBGEM_dighit_nchan;
   vector<int>     *Harm_HCal_dighit_chan;
   vector<int>     *Harm_HCal_dighit_adc_0;
   vector<int>     *Harm_HCal_dighit_adc_1;
   vector<int>     *Harm_HCal_dighit_adc_2;
   vector<int>     *Harm_HCal_dighit_adc_3;
   vector<int>     *Harm_HCal_dighit_adc_4;
   vector<int>     *Harm_HCal_dighit_adc_5;
   vector<int>     *Harm_HCal_dighit_adc_6;
   vector<int>     *Harm_HCal_dighit_adc_7;
   vector<int>     *Harm_HCal_dighit_adc_8;
   vector<int>     *Harm_HCal_dighit_adc_9;
   vector<int>     *Harm_HCal_dighit_adc_10;
   vector<int>     *Harm_HCal_dighit_adc_11;
   vector<int>     *Harm_HCal_dighit_adc_12;
   vector<int>     *Harm_HCal_dighit_adc_13;
   vector<int>     *Harm_HCal_dighit_adc_14;
   vector<int>     *Harm_HCal_dighit_adc_15;
   vector<int>     *Harm_HCal_dighit_adc_16;
   vector<int>     *Harm_HCal_dighit_adc_17;
   vector<int>     *Harm_HCal_dighit_adc_18;
   vector<int>     *Harm_HCal_dighit_adc_19;
   vector<int>     *Harm_HCal_dighit_tdc;

   // List of branches
   TBranch        *b_ev;   //!
   TBranch        *b_Earm_BBGEM_hit_nhits;   //!
   TBranch        *b_Earm_BBGEM_hit_plane;   //!
   TBranch        *b_Earm_BBGEM_hit_strip;   //!
   TBranch        *b_Earm_BBGEM_hit_x;   //!
   TBranch        *b_Earm_BBGEM_hit_y;   //!
   TBranch        *b_Earm_BBGEM_hit_z;   //!
   TBranch        *b_Earm_BBGEM_hit_polx;   //!
   TBranch        *b_Earm_BBGEM_hit_poly;   //!
   TBranch        *b_Earm_BBGEM_hit_polz;   //!
   TBranch        *b_Earm_BBGEM_hit_t;   //!
   TBranch        *b_Earm_BBGEM_hit_trms;   //!
   TBranch        *b_Earm_BBGEM_hit_tmin;   //!
   TBranch        *b_Earm_BBGEM_hit_tmax;   //!
   TBranch        *b_Earm_BBGEM_hit_tx;   //!
   TBranch        *b_Earm_BBGEM_hit_ty;   //!
   TBranch        *b_Earm_BBGEM_hit_xin;   //!
   TBranch        *b_Earm_BBGEM_hit_yin;   //!
   TBranch        *b_Earm_BBGEM_hit_zin;   //!
   TBranch        *b_Earm_BBGEM_hit_xout;   //!
   TBranch        *b_Earm_BBGEM_hit_yout;   //!
   TBranch        *b_Earm_BBGEM_hit_zout;   //!
   TBranch        *b_Earm_BBGEM_hit_txp;   //!
   TBranch        *b_Earm_BBGEM_hit_typ;   //!
   TBranch        *b_Earm_BBGEM_hit_xg;   //!
   TBranch        *b_Earm_BBGEM_hit_yg;   //!
   TBranch        *b_Earm_BBGEM_hit_zg;   //!
   TBranch        *b_Earm_BBGEM_hit_trid;   //!
   TBranch        *b_Earm_BBGEM_hit_mid;   //!
   TBranch        *b_Earm_BBGEM_hit_pid;   //!
   TBranch        *b_Earm_BBGEM_hit_vx;   //!
   TBranch        *b_Earm_BBGEM_hit_vy;   //!
   TBranch        *b_Earm_BBGEM_hit_vz;   //!
   TBranch        *b_Earm_BBGEM_hit_p;   //!
   TBranch        *b_Earm_BBGEM_hit_edep;   //!
   TBranch        *b_Earm_BBGEM_hit_beta;   //!
   TBranch        *b_Earm_BBGEM_hit_otridx;   //!
   TBranch        *b_Earm_BBGEM_hit_ptridx;   //!
   TBranch        *b_Earm_BBGEM_hit_sdtridx;   //!
   TBranch        *b_Earm_BBGEM_Track_ntracks;   //!
   TBranch        *b_Earm_BBGEM_Track_TID;   //!
   TBranch        *b_Earm_BBGEM_Track_PID;   //!
   TBranch        *b_Earm_BBGEM_Track_MID;   //!
   TBranch        *b_Earm_BBGEM_Track_NumHits;   //!
   TBranch        *b_Earm_BBGEM_Track_NumPlanes;   //!
   TBranch        *b_Earm_BBGEM_Track_NDF;   //!
   TBranch        *b_Earm_BBGEM_Track_Chi2fit;   //!
   TBranch        *b_Earm_BBGEM_Track_Chi2true;   //!
   TBranch        *b_Earm_BBGEM_Track_X;   //!
   TBranch        *b_Earm_BBGEM_Track_Y;   //!
   TBranch        *b_Earm_BBGEM_Track_Xp;   //!
   TBranch        *b_Earm_BBGEM_Track_Yp;   //!
   TBranch        *b_Earm_BBGEM_Track_T;   //!
   TBranch        *b_Earm_BBGEM_Track_P;   //!
   TBranch        *b_Earm_BBGEM_Track_Sx;   //!
   TBranch        *b_Earm_BBGEM_Track_Sy;   //!
   TBranch        *b_Earm_BBGEM_Track_Sz;   //!
   TBranch        *b_Earm_BBGEM_Track_Xfit;   //!
   TBranch        *b_Earm_BBGEM_Track_Yfit;   //!
   TBranch        *b_Earm_BBGEM_Track_Xpfit;   //!
   TBranch        *b_Earm_BBGEM_Track_Ypfit;   //!
   TBranch        *b_Earm_BBGEM_Track_otridx;   //!
   TBranch        *b_Earm_BBGEM_Track_ptridx;   //!
   TBranch        *b_Earm_BBGEM_Track_sdtridx;   //!
   TBranch        *b_Earm_BBHodoScint_det_esum;   //!
   TBranch        *b_Earm_BBHodoScint_hit_nhits;   //!
   TBranch        *b_Earm_BBHodoScint_hit_row;   //!
   TBranch        *b_Earm_BBHodoScint_hit_col;   //!
   TBranch        *b_Earm_BBHodoScint_hit_cell;   //!
   TBranch        *b_Earm_BBHodoScint_hit_plane;   //!
   TBranch        *b_Earm_BBHodoScint_hit_wire;   //!
   TBranch        *b_Earm_BBHodoScint_hit_xcell;   //!
   TBranch        *b_Earm_BBHodoScint_hit_ycell;   //!
   TBranch        *b_Earm_BBHodoScint_hit_zcell;   //!
   TBranch        *b_Earm_BBHodoScint_hit_xcellg;   //!
   TBranch        *b_Earm_BBHodoScint_hit_ycellg;   //!
   TBranch        *b_Earm_BBHodoScint_hit_zcellg;   //!
   TBranch        *b_Earm_BBHodoScint_hit_xhit;   //!
   TBranch        *b_Earm_BBHodoScint_hit_yhit;   //!
   TBranch        *b_Earm_BBHodoScint_hit_zhit;   //!
   TBranch        *b_Earm_BBHodoScint_hit_xhitg;   //!
   TBranch        *b_Earm_BBHodoScint_hit_yhitg;   //!
   TBranch        *b_Earm_BBHodoScint_hit_zhitg;   //!
   TBranch        *b_Earm_BBHodoScint_hit_sumedep;   //!
   TBranch        *b_Earm_BBHodoScint_hit_tavg;   //!
   TBranch        *b_Earm_BBHodoScint_hit_trms;   //!
   TBranch        *b_Earm_BBHodoScint_hit_tmin;   //!
   TBranch        *b_Earm_BBHodoScint_hit_tmax;   //!
   TBranch        *b_Earm_BBHodoScint_hit_otridx;   //!
   TBranch        *b_Earm_BBHodoScint_hit_ptridx;   //!
   TBranch        *b_Earm_BBHodoScint_hit_sdtridx;   //!
   TBranch        *b_Earm_BBPS_hit_nhits;   //!
   TBranch        *b_Earm_BBPS_hit_PMT;   //!
   TBranch        *b_Earm_BBPS_hit_row;   //!
   TBranch        *b_Earm_BBPS_hit_col;   //!
   TBranch        *b_Earm_BBPS_hit_plane;   //!
   TBranch        *b_Earm_BBPS_hit_xcell;   //!
   TBranch        *b_Earm_BBPS_hit_ycell;   //!
   TBranch        *b_Earm_BBPS_hit_zcell;   //!
   TBranch        *b_Earm_BBPS_hit_xgcell;   //!
   TBranch        *b_Earm_BBPS_hit_ygcell;   //!
   TBranch        *b_Earm_BBPS_hit_zgcell;   //!
   TBranch        *b_Earm_BBPS_hit_NumPhotoelectrons;   //!
   TBranch        *b_Earm_BBPS_hit_Time_avg;   //!
   TBranch        *b_Earm_BBPS_hit_Time_rms;   //!
   TBranch        *b_Earm_BBPS_hit_Time_min;   //!
   TBranch        *b_Earm_BBPS_hit_Time_max;   //!
   TBranch        *b_Earm_BBPS_hit_otridx;   //!
   TBranch        *b_Earm_BBPS_hit_ptridx;   //!
   TBranch        *b_Earm_BBPS_hit_sdtridx;   //!
   TBranch        *b_Earm_BBPSTF1_det_esum;   //!
   TBranch        *b_Earm_BBPSTF1_hit_nhits;   //!
   TBranch        *b_Earm_BBPSTF1_hit_row;   //!
   TBranch        *b_Earm_BBPSTF1_hit_col;   //!
   TBranch        *b_Earm_BBPSTF1_hit_cell;   //!
   TBranch        *b_Earm_BBPSTF1_hit_plane;   //!
   TBranch        *b_Earm_BBPSTF1_hit_wire;   //!
   TBranch        *b_Earm_BBPSTF1_hit_xcell;   //!
   TBranch        *b_Earm_BBPSTF1_hit_ycell;   //!
   TBranch        *b_Earm_BBPSTF1_hit_zcell;   //!
   TBranch        *b_Earm_BBPSTF1_hit_xcellg;   //!
   TBranch        *b_Earm_BBPSTF1_hit_ycellg;   //!
   TBranch        *b_Earm_BBPSTF1_hit_zcellg;   //!
   TBranch        *b_Earm_BBPSTF1_hit_xhit;   //!
   TBranch        *b_Earm_BBPSTF1_hit_yhit;   //!
   TBranch        *b_Earm_BBPSTF1_hit_zhit;   //!
   TBranch        *b_Earm_BBPSTF1_hit_xhitg;   //!
   TBranch        *b_Earm_BBPSTF1_hit_yhitg;   //!
   TBranch        *b_Earm_BBPSTF1_hit_zhitg;   //!
   TBranch        *b_Earm_BBPSTF1_hit_sumedep;   //!
   TBranch        *b_Earm_BBPSTF1_hit_tavg;   //!
   TBranch        *b_Earm_BBPSTF1_hit_trms;   //!
   TBranch        *b_Earm_BBPSTF1_hit_tmin;   //!
   TBranch        *b_Earm_BBPSTF1_hit_tmax;   //!
   TBranch        *b_Earm_BBPSTF1_hit_otridx;   //!
   TBranch        *b_Earm_BBPSTF1_hit_ptridx;   //!
   TBranch        *b_Earm_BBPSTF1_hit_sdtridx;   //!
   TBranch        *b_Earm_BBSH_hit_nhits;   //!
   TBranch        *b_Earm_BBSH_hit_PMT;   //!
   TBranch        *b_Earm_BBSH_hit_row;   //!
   TBranch        *b_Earm_BBSH_hit_col;   //!
   TBranch        *b_Earm_BBSH_hit_plane;   //!
   TBranch        *b_Earm_BBSH_hit_xcell;   //!
   TBranch        *b_Earm_BBSH_hit_ycell;   //!
   TBranch        *b_Earm_BBSH_hit_zcell;   //!
   TBranch        *b_Earm_BBSH_hit_xgcell;   //!
   TBranch        *b_Earm_BBSH_hit_ygcell;   //!
   TBranch        *b_Earm_BBSH_hit_zgcell;   //!
   TBranch        *b_Earm_BBSH_hit_NumPhotoelectrons;   //!
   TBranch        *b_Earm_BBSH_hit_Time_avg;   //!
   TBranch        *b_Earm_BBSH_hit_Time_rms;   //!
   TBranch        *b_Earm_BBSH_hit_Time_min;   //!
   TBranch        *b_Earm_BBSH_hit_Time_max;   //!
   TBranch        *b_Earm_BBSH_hit_otridx;   //!
   TBranch        *b_Earm_BBSH_hit_ptridx;   //!
   TBranch        *b_Earm_BBSH_hit_sdtridx;   //!
   TBranch        *b_Earm_BBSHTF1_det_esum;   //!
   TBranch        *b_Earm_BBSHTF1_hit_nhits;   //!
   TBranch        *b_Earm_BBSHTF1_hit_row;   //!
   TBranch        *b_Earm_BBSHTF1_hit_col;   //!
   TBranch        *b_Earm_BBSHTF1_hit_cell;   //!
   TBranch        *b_Earm_BBSHTF1_hit_plane;   //!
   TBranch        *b_Earm_BBSHTF1_hit_wire;   //!
   TBranch        *b_Earm_BBSHTF1_hit_xcell;   //!
   TBranch        *b_Earm_BBSHTF1_hit_ycell;   //!
   TBranch        *b_Earm_BBSHTF1_hit_zcell;   //!
   TBranch        *b_Earm_BBSHTF1_hit_xcellg;   //!
   TBranch        *b_Earm_BBSHTF1_hit_ycellg;   //!
   TBranch        *b_Earm_BBSHTF1_hit_zcellg;   //!
   TBranch        *b_Earm_BBSHTF1_hit_xhit;   //!
   TBranch        *b_Earm_BBSHTF1_hit_yhit;   //!
   TBranch        *b_Earm_BBSHTF1_hit_zhit;   //!
   TBranch        *b_Earm_BBSHTF1_hit_xhitg;   //!
   TBranch        *b_Earm_BBSHTF1_hit_yhitg;   //!
   TBranch        *b_Earm_BBSHTF1_hit_zhitg;   //!
   TBranch        *b_Earm_BBSHTF1_hit_sumedep;   //!
   TBranch        *b_Earm_BBSHTF1_hit_tavg;   //!
   TBranch        *b_Earm_BBSHTF1_hit_trms;   //!
   TBranch        *b_Earm_BBSHTF1_hit_tmin;   //!
   TBranch        *b_Earm_BBSHTF1_hit_tmax;   //!
   TBranch        *b_Earm_BBSHTF1_hit_otridx;   //!
   TBranch        *b_Earm_BBSHTF1_hit_ptridx;   //!
   TBranch        *b_Earm_BBSHTF1_hit_sdtridx;   //!
   TBranch        *b_Earm_GRINCH_hit_nhits;   //!
   TBranch        *b_Earm_GRINCH_hit_PMT;   //!
   TBranch        *b_Earm_GRINCH_hit_row;   //!
   TBranch        *b_Earm_GRINCH_hit_col;   //!
   TBranch        *b_Earm_GRINCH_hit_xpmt;   //!
   TBranch        *b_Earm_GRINCH_hit_ypmt;   //!
   TBranch        *b_Earm_GRINCH_hit_zpmt;   //!
   TBranch        *b_Earm_GRINCH_hit_xgpmt;   //!
   TBranch        *b_Earm_GRINCH_hit_ygpmt;   //!
   TBranch        *b_Earm_GRINCH_hit_zgpmt;   //!
   TBranch        *b_Earm_GRINCH_hit_NumPhotoelectrons;   //!
   TBranch        *b_Earm_GRINCH_hit_Time_avg;   //!
   TBranch        *b_Earm_GRINCH_hit_Time_rms;   //!
   TBranch        *b_Earm_GRINCH_hit_Time_min;   //!
   TBranch        *b_Earm_GRINCH_hit_Time_max;   //!
   TBranch        *b_Earm_GRINCH_hit_mTrackNo;   //!
   TBranch        *b_Earm_GRINCH_hit_xhit;   //!
   TBranch        *b_Earm_GRINCH_hit_yhit;   //!
   TBranch        *b_Earm_GRINCH_hit_zhit;   //!
   TBranch        *b_Earm_GRINCH_hit_pxhit;   //!
   TBranch        *b_Earm_GRINCH_hit_pyhit;   //!
   TBranch        *b_Earm_GRINCH_hit_pzhit;   //!
   TBranch        *b_Earm_GRINCH_hit_pvx;   //!
   TBranch        *b_Earm_GRINCH_hit_pvy;   //!
   TBranch        *b_Earm_GRINCH_hit_pvz;   //!
   TBranch        *b_Earm_GRINCH_hit_ppx;   //!
   TBranch        *b_Earm_GRINCH_hit_ppy;   //!
   TBranch        *b_Earm_GRINCH_hit_ppz;   //!
   TBranch        *b_Earm_GRINCH_hit_volume_flag;   //!
   TBranch        *b_Earm_GRINCH_hit_otridx;   //!
   TBranch        *b_Earm_GRINCH_hit_ptridx;   //!
   TBranch        *b_Earm_GRINCH_hit_sdtridx;   //!
   TBranch        *b_Harm_HCal_hit_nhits;   //!
   TBranch        *b_Harm_HCal_hit_PMT;   //!
   TBranch        *b_Harm_HCal_hit_row;   //!
   TBranch        *b_Harm_HCal_hit_col;   //!
   TBranch        *b_Harm_HCal_hit_plane;   //!
   TBranch        *b_Harm_HCal_hit_xcell;   //!
   TBranch        *b_Harm_HCal_hit_ycell;   //!
   TBranch        *b_Harm_HCal_hit_zcell;   //!
   TBranch        *b_Harm_HCal_hit_xgcell;   //!
   TBranch        *b_Harm_HCal_hit_ygcell;   //!
   TBranch        *b_Harm_HCal_hit_zgcell;   //!
   TBranch        *b_Harm_HCal_hit_NumPhotoelectrons;   //!
   TBranch        *b_Harm_HCal_hit_Time_avg;   //!
   TBranch        *b_Harm_HCal_hit_Time_rms;   //!
   TBranch        *b_Harm_HCal_hit_Time_min;   //!
   TBranch        *b_Harm_HCal_hit_Time_max;   //!
   TBranch        *b_Harm_HCal_hit_otridx;   //!
   TBranch        *b_Harm_HCal_hit_ptridx;   //!
   TBranch        *b_Harm_HCal_hit_sdtridx;   //!
   TBranch        *b_Harm_HCalScint_det_esum;   //!
   TBranch        *b_Harm_HCalScint_hit_nhits;   //!
   TBranch        *b_Harm_HCalScint_hit_row;   //!
   TBranch        *b_Harm_HCalScint_hit_col;   //!
   TBranch        *b_Harm_HCalScint_hit_cell;   //!
   TBranch        *b_Harm_HCalScint_hit_plane;   //!
   TBranch        *b_Harm_HCalScint_hit_wire;   //!
   TBranch        *b_Harm_HCalScint_hit_xcell;   //!
   TBranch        *b_Harm_HCalScint_hit_ycell;   //!
   TBranch        *b_Harm_HCalScint_hit_zcell;   //!
   TBranch        *b_Harm_HCalScint_hit_xcellg;   //!
   TBranch        *b_Harm_HCalScint_hit_ycellg;   //!
   TBranch        *b_Harm_HCalScint_hit_zcellg;   //!
   TBranch        *b_Harm_HCalScint_hit_xhit;   //!
   TBranch        *b_Harm_HCalScint_hit_yhit;   //!
   TBranch        *b_Harm_HCalScint_hit_zhit;   //!
   TBranch        *b_Harm_HCalScint_hit_xhitg;   //!
   TBranch        *b_Harm_HCalScint_hit_yhitg;   //!
   TBranch        *b_Harm_HCalScint_hit_zhitg;   //!
   TBranch        *b_Harm_HCalScint_hit_sumedep;   //!
   TBranch        *b_Harm_HCalScint_hit_tavg;   //!
   TBranch        *b_Harm_HCalScint_hit_trms;   //!
   TBranch        *b_Harm_HCalScint_hit_tmin;   //!
   TBranch        *b_Harm_HCalScint_hit_tmax;   //!
   TBranch        *b_Harm_HCalScint_hit_otridx;   //!
   TBranch        *b_Harm_HCalScint_hit_ptridx;   //!
   TBranch        *b_Harm_HCalScint_hit_sdtridx;   //!
   TBranch        *b_OTrack_ntracks;   //!
   TBranch        *b_OTrack_TID;   //!
   TBranch        *b_OTrack_MID;   //!
   TBranch        *b_OTrack_PID;   //!
   TBranch        *b_OTrack_posx;   //!
   TBranch        *b_OTrack_posy;   //!
   TBranch        *b_OTrack_posz;   //!
   TBranch        *b_OTrack_momx;   //!
   TBranch        *b_OTrack_momy;   //!
   TBranch        *b_OTrack_momz;   //!
   TBranch        *b_OTrack_polx;   //!
   TBranch        *b_OTrack_poly;   //!
   TBranch        *b_OTrack_polz;   //!
   TBranch        *b_OTrack_Etot;   //!
   TBranch        *b_OTrack_T;   //!
   TBranch        *b_PTrack_ntracks;   //!
   TBranch        *b_PTrack_TID;   //!
   TBranch        *b_PTrack_PID;   //!
   TBranch        *b_PTrack_posx;   //!
   TBranch        *b_PTrack_posy;   //!
   TBranch        *b_PTrack_posz;   //!
   TBranch        *b_PTrack_momx;   //!
   TBranch        *b_PTrack_momy;   //!
   TBranch        *b_PTrack_momz;   //!
   TBranch        *b_PTrack_polx;   //!
   TBranch        *b_PTrack_poly;   //!
   TBranch        *b_PTrack_polz;   //!
   TBranch        *b_PTrack_Etot;   //!
   TBranch        *b_PTrack_T;   //!
   TBranch        *b_SDTrack_ntracks;   //!
   TBranch        *b_SDTrack_TID;   //!
   TBranch        *b_SDTrack_MID;   //!
   TBranch        *b_SDTrack_PID;   //!
   TBranch        *b_SDTrack_posx;   //!
   TBranch        *b_SDTrack_posy;   //!
   TBranch        *b_SDTrack_posz;   //!
   TBranch        *b_SDTrack_momx;   //!
   TBranch        *b_SDTrack_momy;   //!
   TBranch        *b_SDTrack_momz;   //!
   TBranch        *b_SDTrack_polx;   //!
   TBranch        *b_SDTrack_poly;   //!
   TBranch        *b_SDTrack_polz;   //!
   TBranch        *b_SDTrack_Etot;   //!
   TBranch        *b_SDTrack_T;   //!
   TBranch        *b_Earm_BBPS_dighit_nchan;   //!
   TBranch        *b_Earm_BBPS_dighit_chan;   //!
   TBranch        *b_Earm_BBPS_dighit_adc;   //!
   TBranch        *b_Earm_BBSH_dighit_nchan;   //!
   TBranch        *b_Earm_BBSH_dighit_chan;   //!
   TBranch        *b_Earm_BBSH_dighit_adc;   //!
   TBranch        *b_Earm_BBHodo_dighit_nchan;   //!
   TBranch        *b_Earm_BBHodo_dighit_chan;   //!
   TBranch        *b_Earm_BBHodo_dighit_adc;   //!
   TBranch        *b_Earm_BBHodo_dighit_tdc_l;   //!
   TBranch        *b_Earm_BBHodo_dighit_tdc_t;   //!
   TBranch        *b_Earm_GRINCH_dighit_nchan;   //!
   TBranch        *b_Earm_GRINCH_dighit_chan;   //!
   TBranch        *b_Earm_GRINCH_dighit_adc;   //!
   TBranch        *b_Earm_GRINCH_dighit_tdc_l;   //!
   TBranch        *b_Earm_GRINCH_dighit_tdc_t;   //!
   TBranch        *b_Earm_BBGEM_1x_dighit_nstrips;   //!
   TBranch        *b_Earm_BBGEM_1x_dighit_strip;   //!
   TBranch        *b_Earm_BBGEM_1x_dighit_adc_0;   //!
   TBranch        *b_Earm_BBGEM_1x_dighit_adc_1;   //!
   TBranch        *b_Earm_BBGEM_1x_dighit_adc_2;   //!
   TBranch        *b_Earm_BBGEM_1x_dighit_adc_3;   //!
   TBranch        *b_Earm_BBGEM_1x_dighit_adc_4;   //!
   TBranch        *b_Earm_BBGEM_1x_dighit_adc_5;   //!
   TBranch        *b_Earm_BBGEM_1y_dighit_nstrips;   //!
   TBranch        *b_Earm_BBGEM_1y_dighit_strip;   //!
   TBranch        *b_Earm_BBGEM_1y_dighit_adc_0;   //!
   TBranch        *b_Earm_BBGEM_1y_dighit_adc_1;   //!
   TBranch        *b_Earm_BBGEM_1y_dighit_adc_2;   //!
   TBranch        *b_Earm_BBGEM_1y_dighit_adc_3;   //!
   TBranch        *b_Earm_BBGEM_1y_dighit_adc_4;   //!
   TBranch        *b_Earm_BBGEM_1y_dighit_adc_5;   //!
   TBranch        *b_Earm_BBGEM_2x_dighit_nstrips;   //!
   TBranch        *b_Earm_BBGEM_2x_dighit_strip;   //!
   TBranch        *b_Earm_BBGEM_2x_dighit_adc_0;   //!
   TBranch        *b_Earm_BBGEM_2x_dighit_adc_1;   //!
   TBranch        *b_Earm_BBGEM_2x_dighit_adc_2;   //!
   TBranch        *b_Earm_BBGEM_2x_dighit_adc_3;   //!
   TBranch        *b_Earm_BBGEM_2x_dighit_adc_4;   //!
   TBranch        *b_Earm_BBGEM_2x_dighit_adc_5;   //!
   TBranch        *b_Earm_BBGEM_2y_dighit_nstrips;   //!
   TBranch        *b_Earm_BBGEM_2y_dighit_strip;   //!
   TBranch        *b_Earm_BBGEM_2y_dighit_adc_0;   //!
   TBranch        *b_Earm_BBGEM_2y_dighit_adc_1;   //!
   TBranch        *b_Earm_BBGEM_2y_dighit_adc_2;   //!
   TBranch        *b_Earm_BBGEM_2y_dighit_adc_3;   //!
   TBranch        *b_Earm_BBGEM_2y_dighit_adc_4;   //!
   TBranch        *b_Earm_BBGEM_2y_dighit_adc_5;   //!
   TBranch        *b_Earm_BBGEM_3x_dighit_nstrips;   //!
   TBranch        *b_Earm_BBGEM_3x_dighit_strip;   //!
   TBranch        *b_Earm_BBGEM_3x_dighit_adc_0;   //!
   TBranch        *b_Earm_BBGEM_3x_dighit_adc_1;   //!
   TBranch        *b_Earm_BBGEM_3x_dighit_adc_2;   //!
   TBranch        *b_Earm_BBGEM_3x_dighit_adc_3;   //!
   TBranch        *b_Earm_BBGEM_3x_dighit_adc_4;   //!
   TBranch        *b_Earm_BBGEM_3x_dighit_adc_5;   //!
   TBranch        *b_Earm_BBGEM_3y_dighit_nstrips;   //!
   TBranch        *b_Earm_BBGEM_3y_dighit_strip;   //!
   TBranch        *b_Earm_BBGEM_3y_dighit_adc_0;   //!
   TBranch        *b_Earm_BBGEM_3y_dighit_adc_1;   //!
   TBranch        *b_Earm_BBGEM_3y_dighit_adc_2;   //!
   TBranch        *b_Earm_BBGEM_3y_dighit_adc_3;   //!
   TBranch        *b_Earm_BBGEM_3y_dighit_adc_4;   //!
   TBranch        *b_Earm_BBGEM_3y_dighit_adc_5;   //!
   TBranch        *b_Earm_BBGEM_4x_dighit_nstrips;   //!
   TBranch        *b_Earm_BBGEM_4x_dighit_strip;   //!
   TBranch        *b_Earm_BBGEM_4x_dighit_adc_0;   //!
   TBranch        *b_Earm_BBGEM_4x_dighit_adc_1;   //!
   TBranch        *b_Earm_BBGEM_4x_dighit_adc_2;   //!
   TBranch        *b_Earm_BBGEM_4x_dighit_adc_3;   //!
   TBranch        *b_Earm_BBGEM_4x_dighit_adc_4;   //!
   TBranch        *b_Earm_BBGEM_4x_dighit_adc_5;   //!
   TBranch        *b_Earm_BBGEM_4y_dighit_nstrips;   //!
   TBranch        *b_Earm_BBGEM_4y_dighit_strip;   //!
   TBranch        *b_Earm_BBGEM_4y_dighit_adc_0;   //!
   TBranch        *b_Earm_BBGEM_4y_dighit_adc_1;   //!
   TBranch        *b_Earm_BBGEM_4y_dighit_adc_2;   //!
   TBranch        *b_Earm_BBGEM_4y_dighit_adc_3;   //!
   TBranch        *b_Earm_BBGEM_4y_dighit_adc_4;   //!
   TBranch        *b_Earm_BBGEM_4y_dighit_adc_5;   //!
   TBranch        *b_Earm_BBGEM_5x_dighit_nstrips;   //!
   TBranch        *b_Earm_BBGEM_5x_dighit_strip;   //!
   TBranch        *b_Earm_BBGEM_5x_dighit_adc_0;   //!
   TBranch        *b_Earm_BBGEM_5x_dighit_adc_1;   //!
   TBranch        *b_Earm_BBGEM_5x_dighit_adc_2;   //!
   TBranch        *b_Earm_BBGEM_5x_dighit_adc_3;   //!
   TBranch        *b_Earm_BBGEM_5x_dighit_adc_4;   //!
   TBranch        *b_Earm_BBGEM_5x_dighit_adc_5;   //!
   TBranch        *b_Earm_BBGEM_5y_dighit_nstrips;   //!
   TBranch        *b_Earm_BBGEM_5y_dighit_strip;   //!
   TBranch        *b_Earm_BBGEM_5y_dighit_adc_0;   //!
   TBranch        *b_Earm_BBGEM_5y_dighit_adc_1;   //!
   TBranch        *b_Earm_BBGEM_5y_dighit_adc_2;   //!
   TBranch        *b_Earm_BBGEM_5y_dighit_adc_3;   //!
   TBranch        *b_Earm_BBGEM_5y_dighit_adc_4;   //!
   TBranch        *b_Earm_BBGEM_5y_dighit_adc_5;   //!
   TBranch        *b_Earm_BBGEM_dighit_nchan;   //!
   TBranch        *b_Harm_HCal_dighit_chan;   //!
   TBranch        *b_Harm_HCal_dighit_adc_0;   //!
   TBranch        *b_Harm_HCal_dighit_adc_1;   //!
   TBranch        *b_Harm_HCal_dighit_adc_2;   //!
   TBranch        *b_Harm_HCal_dighit_adc_3;   //!
   TBranch        *b_Harm_HCal_dighit_adc_4;   //!
   TBranch        *b_Harm_HCal_dighit_adc_5;   //!
   TBranch        *b_Harm_HCal_dighit_adc_6;   //!
   TBranch        *b_Harm_HCal_dighit_adc_7;   //!
   TBranch        *b_Harm_HCal_dighit_adc_8;   //!
   TBranch        *b_Harm_HCal_dighit_adc_9;   //!
   TBranch        *b_Harm_HCal_dighit_adc_10;   //!
   TBranch        *b_Harm_HCal_dighit_adc_11;   //!
   TBranch        *b_Harm_HCal_dighit_adc_12;   //!
   TBranch        *b_Harm_HCal_dighit_adc_13;   //!
   TBranch        *b_Harm_HCal_dighit_adc_14;   //!
   TBranch        *b_Harm_HCal_dighit_adc_15;   //!
   TBranch        *b_Harm_HCal_dighit_adc_16;   //!
   TBranch        *b_Harm_HCal_dighit_adc_17;   //!
   TBranch        *b_Harm_HCal_dighit_adc_18;   //!
   TBranch        *b_Harm_HCal_dighit_adc_19;   //!
   TBranch        *b_Harm_HCal_dighit_tdc;   //!

   gmn_dig_tree(TTree *tree=0);
   virtual ~gmn_dig_tree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

//#ifdef gmn_dig_tree_cxx

//#endif // #ifdef gmn_dig_tree_cxx
