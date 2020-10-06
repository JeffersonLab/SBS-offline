#include "g4sbs_data.h"
#include <TTree.h>
#include <TBranch.h>
#include <iostream>


namespace TSBSGeant4 {
  template<typename T>
  int VDetData_t::SetupBranch(TTree *tree, const char* prefix,
      const char* varname, T &var)
  {
    TString branchname = TString::Format("%s.%s",prefix,varname);
    if(!tree)
      return 1;
    var = 0;
    int ret = tree->SetBranchAddress(branchname.Data(),&var);
    if( ret != 0 ) {
      std::cerr << "Unable to set branch '" << branchname
        << "' failed with error code: " << ret << std::endl;
      return 1;
    }

    return 0;
  }

  bool CalData_t::SetupBranches(TTree *tree, const char* prefix)
  {
    int ret = 0;
    ret += SetupBranch(tree,prefix,"nhits", nhits);
    ret += SetupBranch(tree,prefix,"row", row);
    ret += SetupBranch(tree,prefix,"col", col);
    ret += SetupBranch(tree,prefix,"cell", cell);
    ret += SetupBranch(tree,prefix,"plane", plane);
    ret += SetupBranch(tree,prefix,"xcell", xcell);
    ret += SetupBranch(tree,prefix,"ycell", ycell);
    ret += SetupBranch(tree,prefix,"zcell", zcell);
    ret += SetupBranch(tree,prefix,"xcellg", xcellg);
    ret += SetupBranch(tree,prefix,"ycellg", ycellg);
    ret += SetupBranch(tree,prefix,"zcellg", zcellg);
    ret += SetupBranch(tree,prefix,"xhit", xhit);
    ret += SetupBranch(tree,prefix,"yhit", yhit);
    ret += SetupBranch(tree,prefix,"zhit", zhit);
    ret += SetupBranch(tree,prefix,"xhitg", xhitg);
    ret += SetupBranch(tree,prefix,"yhitg", yhitg);
    ret += SetupBranch(tree,prefix,"zhitg", zhitg);
    ret += SetupBranch(tree,prefix,"sumedep", sumedep);
    ret += SetupBranch(tree,prefix,"tavg", tavg);
    ret += SetupBranch(tree,prefix,"trms", trms);
    ret += SetupBranch(tree,prefix,"tmin", tmin);
    ret += SetupBranch(tree,prefix,"tmax", tmax);
    return (ret ==0);
  }
  
  /*
  bool ECalData_t::SetupBranches(TTree *tree, const char* prefix)
  {
    int ret = 0;
    ret += SetupBranch(tree,prefix,"nhits", nhits);
    ret += SetupBranch(tree,prefix,"PMT", PMT);
    ret += SetupBranch(tree,prefix,"row", row);
    ret += SetupBranch(tree,prefix,"col", col);
    ret += SetupBranch(tree,prefix,"plane", plane);
    ret += SetupBranch(tree,prefix,"xcell", xcell);
    ret += SetupBranch(tree,prefix,"ycell", ycell);
    ret += SetupBranch(tree,prefix,"zcell", zcell);
    ret += SetupBranch(tree,prefix,"xgcell", xgcell);
    ret += SetupBranch(tree,prefix,"ygcell", ygcell);
    ret += SetupBranch(tree,prefix,"zgcell", zgcell);
    ret += SetupBranch(tree,prefix,"NumPhotoelectrons", NumPhotoelectrons);
    ret += SetupBranch(tree,prefix,"Time_avg", Time_avg);
    ret += SetupBranch(tree,prefix,"Time_rms", Time_rms);
    ret += SetupBranch(tree,prefix,"Time_min", Time_min);
    ret += SetupBranch(tree,prefix,"Time_max", Time_max);
    return (ret==0);
  }

  bool ECalPartData_t::SetupBranches(TTree* tree, const char* prefix)
  {
    int ret = 0;
    ret += SetupBranch(tree,prefix,"npart_ECAL", npart_ECAL);
    ret += SetupBranch(tree,prefix,"E", E);
    ret += SetupBranch(tree,prefix,"t", t);
    ret += SetupBranch(tree,prefix,"part_PMT", part_PMT);
    ret += SetupBranch(tree,prefix,"detected", detected);
    return (ret==0);
    //return ECalData_t::SetupBranches(tree,prefix);
  }
  */
  
  bool RICHData_t::SetupBranches(TTree* tree, const char* prefix)
  {
    int ret = 0;
    ret += SetupBranch(tree,prefix,"nhits",nhits);
    ret += SetupBranch(tree,prefix,"PMT",PMT);
    ret += SetupBranch(tree,prefix,"row",row);
    ret += SetupBranch(tree,prefix,"col",col);
    ret += SetupBranch(tree,prefix,"xpmt",xpmt);
    ret += SetupBranch(tree,prefix,"ypmt",ypmt);
    ret += SetupBranch(tree,prefix,"zpmt",zpmt);
    ret += SetupBranch(tree,prefix,"xgpmt",xgpmt);
    ret += SetupBranch(tree,prefix,"ygpmt",ygpmt);
    ret += SetupBranch(tree,prefix,"zgpmt",zgpmt);
    ret += SetupBranch(tree,prefix,"NumPhotoelectrons",NumPhotoelectrons);
    ret += SetupBranch(tree,prefix,"Time_avg",Time_avg);
    ret += SetupBranch(tree,prefix,"Time_rms",Time_rms);
    ret += SetupBranch(tree,prefix,"Time_min",Time_min);
    ret += SetupBranch(tree,prefix,"Time_max",Time_max);
    ret += SetupBranch(tree,prefix,"mTrackNo",mTrackNo);
    ret += SetupBranch(tree,prefix,"xhit",xhit);
    ret += SetupBranch(tree,prefix,"yhit",yhit);
    ret += SetupBranch(tree,prefix,"zhit",zhit);
    ret += SetupBranch(tree,prefix,"pxhit",pxhit);
    ret += SetupBranch(tree,prefix,"pyhit",pyhit);
    ret += SetupBranch(tree,prefix,"pzhit",pzhit);
    ret += SetupBranch(tree,prefix,"pvx",pvx);
    ret += SetupBranch(tree,prefix,"pvy",pvy);
    ret += SetupBranch(tree,prefix,"pvz",pvz);
    ret += SetupBranch(tree,prefix,"ppx",ppx);
    ret += SetupBranch(tree,prefix,"ppy",ppy);
    ret += SetupBranch(tree,prefix,"ppz",ppz);
    ret += SetupBranch(tree,prefix,"volume_flag",volume_flag);
    return (ret==0);
  }

  bool GEMData_t::SetupBranches(TTree* tree, const char *prefix)
  {
    int ret = 0;
    ret += SetupBranch(tree,prefix,"nhits",nhits);
    ret += SetupBranch(tree,prefix,"plane",plane);
    ret += SetupBranch(tree,prefix,"strip",strip);
    ret += SetupBranch(tree,prefix,"x",x);
    ret += SetupBranch(tree,prefix,"y",y);
    ret += SetupBranch(tree,prefix,"z",z);
    ret += SetupBranch(tree,prefix,"polx",polx);
    ret += SetupBranch(tree,prefix,"poly",poly);
    ret += SetupBranch(tree,prefix,"polz",polz);
    ret += SetupBranch(tree,prefix,"t",t);
    ret += SetupBranch(tree,prefix,"trms",trms);
    ret += SetupBranch(tree,prefix,"tmin",tmin);
    ret += SetupBranch(tree,prefix,"tmax",tmax);
    ret += SetupBranch(tree,prefix,"tx",tx);
    ret += SetupBranch(tree,prefix,"ty",ty);
    ret += SetupBranch(tree,prefix,"txp",txp);
    ret += SetupBranch(tree,prefix,"typ",typ);
    ret += SetupBranch(tree,prefix,"xg",xg);
    ret += SetupBranch(tree,prefix,"yg",yg);
    ret += SetupBranch(tree,prefix,"zg",zg);
    ret += SetupBranch(tree,prefix,"trid",trid);
    ret += SetupBranch(tree,prefix,"mid",mid);
    ret += SetupBranch(tree,prefix,"pid",pid);
    ret += SetupBranch(tree,prefix,"vx",vx);
    ret += SetupBranch(tree,prefix,"vy",vy);
    ret += SetupBranch(tree,prefix,"vz",vz);
    ret += SetupBranch(tree,prefix,"p",p);
    ret += SetupBranch(tree,prefix,"edep",edep);
    ret += SetupBranch(tree,prefix,"beta",beta);
    ret += SetupBranch(tree,prefix,"xin",xin);
    ret += SetupBranch(tree,prefix,"yin",yin);
    ret += SetupBranch(tree,prefix,"zin",zin);
    ret += SetupBranch(tree,prefix,"xout",xout);
    ret += SetupBranch(tree,prefix,"yout",yout);
    ret += SetupBranch(tree,prefix,"zout",zout);
    return (ret==0);
  }
  
  bool TrackerData_t::SetupBranches(TTree* tree, const char *prefix)
  {
    int ret = 0;
    ret += SetupBranch(tree,prefix,"ntracks",ntracks);
    ret += SetupBranch(tree,prefix,"TID",TID);
    ret += SetupBranch(tree,prefix,"PID",PID);
    ret += SetupBranch(tree,prefix,"MID",MID);
    ret += SetupBranch(tree,prefix,"NumHits",NumHits);
    ret += SetupBranch(tree,prefix,"NumPlanes",NumPlanes);
    ret += SetupBranch(tree,prefix,"NDF",NDF);
    ret += SetupBranch(tree,prefix,"Chi2fit",Chi2fit);
    ret += SetupBranch(tree,prefix,"Chi2true",Chi2true);
    ret += SetupBranch(tree,prefix,"X",X);
    ret += SetupBranch(tree,prefix,"Y",Y);
    ret += SetupBranch(tree,prefix,"Xp",Xp);
    ret += SetupBranch(tree,prefix,"Yp",Yp);
    ret += SetupBranch(tree,prefix,"T",T);
    ret += SetupBranch(tree,prefix,"P",P);
    ret += SetupBranch(tree,prefix,"Sx",Sx);
    ret += SetupBranch(tree,prefix,"Sy",Sy);
    ret += SetupBranch(tree,prefix,"Sz",Sz);
    ret += SetupBranch(tree,prefix,"Xfit",Xfit);
    ret += SetupBranch(tree,prefix,"Yfit",Yfit);
    ret += SetupBranch(tree,prefix,"Xpfit",Xpfit);
    ret += SetupBranch(tree,prefix,"Ypfit",Ypfit);
    return (ret==0);
  }
  
  bool DigCalData_t::SetupBranches(TTree* tree, const char *prefix)
  {
    int ret = 0;
    ret += SetupBranch(tree,prefix,"nchan",nchan);
    ret += SetupBranch(tree,prefix,"chan",chan);
    ret += SetupBranch(tree,prefix,"adc",adc);
    return (ret==0);
   /*
    if(!tree)return(false);
    chan = new std::vector<int>;
    adc = new std::vector<int>;
    b_nchan = tree->Branch(Form("%s.nchan", prefix), &nchan);
    b_chan = tree->Branch(Form("%s.chan", prefix), &chan);
    b_adc = tree->Branch(Form("%s.adc", prefix), &adc);
    return true;
    */
  }
  /*
  void DigCalData_t::ClearBranches()
  {
    if(chan){//if one var is defined they all are
      nchan = 0;
      chan->clear();
      adc->clear();
    }
  }
  
  void DigCalData_t::FillBranches()
  {
    if(b_nchan){//if one branch is defined they all are
      b_nchan->Fill();
      b_chan->Fill();
      b_adc->Fill();
    }
  }
  */
  
  bool DigTimingData_t::SetupBranches(TTree* tree, const char *prefix)
  {
    int ret = 0;
    ret += SetupBranch(tree,prefix,"nchan",nchan);
    ret += SetupBranch(tree,prefix,"chan",chan);
    ret += SetupBranch(tree,prefix,"adc",adc);
    ret += SetupBranch(tree,prefix,"tdc_l",tdc_l);
    ret += SetupBranch(tree,prefix,"tdc_t",tdc_t);
    return (ret==0);
    /*
    if(!tree)return(false);
    chan = new std::vector<int>;
    adc = new std::vector<int>;
    tdc_l = new std::vector<int>;
    tdc_t = new std::vector<int>;
    
    b_nchan = tree->Branch(Form("%s.nchan", prefix), &nchan);
    b_chan = tree->Branch(Form("%s.chan", prefix), &chan);
    b_adc = tree->Branch(Form("%s.adc", prefix), &adc);
    b_tdc_l = tree->Branch(Form("%s.tdc_l", prefix), &tdc_l);
    b_tdc_t = tree->Branch(Form("%s.tdc_t", prefix), &tdc_t);
    return true;
    */
  }
  /*
  void DigTimingData_t::ClearBranches()
  {
    if(chan){
      nchan = 0;
      chan->clear();
      adc->clear();
      tdc_l->clear();
      tdc_t->clear();
    }
  }
  
  void DigTimingData_t::FillBranches()
  {
    if(b_nchan){
      b_nchan->Fill();
      b_chan->Fill();
      b_adc->Fill();
      b_tdc_l->Fill();
      b_tdc_t->Fill();
    }
  }
  */
  bool DigSampCalData_t::SetupBranches(TTree* tree, const char *prefix)
  {
    int ret = 0;
    ret += SetupBranch(tree,prefix,"nchan",nchan);
    ret += SetupBranch(tree,prefix,"chan",chan);
    ret += SetupBranch(tree,prefix,"adc_0",adc_0);
    ret += SetupBranch(tree,prefix,"adc_1",adc_1);
    ret += SetupBranch(tree,prefix,"adc_2",adc_2);
    ret += SetupBranch(tree,prefix,"adc_3",adc_3);
    ret += SetupBranch(tree,prefix,"adc_4",adc_4);
    ret += SetupBranch(tree,prefix,"adc_5",adc_5);
    ret += SetupBranch(tree,prefix,"adc_6",adc_6);
    ret += SetupBranch(tree,prefix,"adc_7",adc_7);
    ret += SetupBranch(tree,prefix,"adc_8",adc_8);
    ret += SetupBranch(tree,prefix,"adc_9",adc_9);
    ret += SetupBranch(tree,prefix,"adc_10",adc_10);
    ret += SetupBranch(tree,prefix,"adc_11",adc_11);
    ret += SetupBranch(tree,prefix,"adc_12",adc_12);
    ret += SetupBranch(tree,prefix,"adc_13",adc_13);
    ret += SetupBranch(tree,prefix,"adc_14",adc_14);
    ret += SetupBranch(tree,prefix,"adc_15",adc_15);
    ret += SetupBranch(tree,prefix,"adc_16",adc_16);
    ret += SetupBranch(tree,prefix,"adc_17",adc_17);
    ret += SetupBranch(tree,prefix,"adc_18",adc_18);
    ret += SetupBranch(tree,prefix,"adc_19",adc_19);
    ret += SetupBranch(tree,prefix,"tdc",tdc);
    return (ret==0);
    
    /*
    if(!tree)return(false);
    chan = new std::vector<int>;
    adc_0 = new std::vector<int>;
    adc_1 = new std::vector<int>;
    adc_2 = new std::vector<int>;
    adc_3 = new std::vector<int>;
    adc_4 = new std::vector<int>;
    adc_5 = new std::vector<int>;
    adc_6 = new std::vector<int>;
    adc_7 = new std::vector<int>;
    adc_8 = new std::vector<int>;
    adc_9 = new std::vector<int>;
    adc_10 = new std::vector<int>;
    adc_11 = new std::vector<int>;
    adc_12 = new std::vector<int>;
    adc_13 = new std::vector<int>;
    adc_14 = new std::vector<int>;
    adc_15 = new std::vector<int>;
    adc_16 = new std::vector<int>;
    adc_17 = new std::vector<int>;
    adc_18 = new std::vector<int>;
    adc_19 = new std::vector<int>;
    tdc = new std::vector<int>;
    
    b_nchan = tree->Branch(Form("%s.nchan", prefix), &nchan);
    b_chan = tree->Branch(Form("%s.chan", prefix), &chan);
    b_adc_0 = tree->Branch(Form("%s.adc_0", prefix), &adc_0);
    b_adc_1 = tree->Branch(Form("%s.adc_1", prefix), &adc_1);
    b_adc_2 = tree->Branch(Form("%s.adc_2", prefix), &adc_2);
    b_adc_3 = tree->Branch(Form("%s.adc_3", prefix), &adc_3);
    b_adc_4 = tree->Branch(Form("%s.adc_4", prefix), &adc_4);
    b_adc_5 = tree->Branch(Form("%s.adc_5", prefix), &adc_5);
    b_adc_6 = tree->Branch(Form("%s.adc_6", prefix), &adc_6);
    b_adc_7 = tree->Branch(Form("%s.adc_7", prefix), &adc_7);
    b_adc_8 = tree->Branch(Form("%s.adc_8", prefix), &adc_8);
    b_adc_9 = tree->Branch(Form("%s.adc_9", prefix), &adc_9);
    b_adc_10 = tree->Branch(Form("%s.adc_10", prefix), &adc_10);
    b_adc_11 = tree->Branch(Form("%s.adc_11", prefix), &adc_11);
    b_adc_12 = tree->Branch(Form("%s.adc_12", prefix), &adc_12);
    b_adc_13 = tree->Branch(Form("%s.adc_13", prefix), &adc_13);
    b_adc_14 = tree->Branch(Form("%s.adc_14", prefix), &adc_14);
    b_adc_15 = tree->Branch(Form("%s.adc_15", prefix), &adc_15);
    b_adc_16 = tree->Branch(Form("%s.adc_16", prefix), &adc_16);
    b_adc_17 = tree->Branch(Form("%s.adc_17", prefix), &adc_17);
    b_adc_18 = tree->Branch(Form("%s.adc_18", prefix), &adc_18);
    b_adc_19 = tree->Branch(Form("%s.adc_19", prefix), &adc_19);
    b_tdc = tree->Branch(Form("%s.tdc", prefix), &tdc);
    return true;
    */
  }
  /*
  void DigSampCalData_t::ClearBranches()
  {
    if(chan){
      nchan = 0;
      chan->clear();
      adc_0->clear();
      adc_1->clear();
      adc_2->clear();
      adc_3->clear();
      adc_4->clear();
      adc_5->clear();
      adc_6->clear();
      adc_7->clear();
      adc_8->clear();
      adc_9->clear();
      adc_10->clear();
      adc_11->clear();
      adc_12->clear();
      adc_13->clear();
      adc_14->clear();
      adc_15->clear();
      adc_16->clear();
      adc_17->clear();
      adc_18->clear();
      adc_19->clear();
      tdc->clear();
    }
  }
  
  void DigSampCalData_t::FillBranches()
  {
    if(b_nchan){
      b_nchan->Fill();
      b_chan->Fill();
      b_adc_0->Fill();
      b_adc_1->Fill();
      b_adc_2->Fill();
      b_adc_3->Fill();
      b_adc_4->Fill();
      b_adc_5->Fill();
      b_adc_6->Fill();
      b_adc_7->Fill();
      b_adc_8->Fill();
      b_adc_9->Fill();
      b_adc_10->Fill();
      b_adc_11->Fill();
      b_adc_12->Fill();
      b_adc_13->Fill();
      b_adc_14->Fill();
      b_adc_15->Fill();
      b_adc_16->Fill();
      b_adc_17->Fill();
      b_adc_18->Fill();
      b_adc_19->Fill();
      b_tdc->Fill();
    }
  }
  */
  bool DigGEMData_t::SetupBranches(TTree* tree, const char *prefix)
  {
    int ret = 0;
    ret += SetupBranch(tree,prefix,"nstrips",nstrips);
    ret += SetupBranch(tree,prefix,"module",module);
    ret += SetupBranch(tree,prefix,"strip",strip);
    ret += SetupBranch(tree,prefix,"adc_0",adc_0);
    ret += SetupBranch(tree,prefix,"adc_1",adc_1);
    ret += SetupBranch(tree,prefix,"adc_2",adc_2);
    ret += SetupBranch(tree,prefix,"adc_3",adc_3);
    ret += SetupBranch(tree,prefix,"adc_4",adc_4);
    ret += SetupBranch(tree,prefix,"adc_5",adc_5);
    return (ret==0);
   
    /*
    if(!tree)return(false);
    b_nstrips = tree->Branch(Form("%s.nstrips", prefix), &nstrips);
    b_module = tree->Branch(Form("%s.module", prefix), &module);
    b_strip = tree->Branch(Form("%s.strip", prefix), &strip);
    b_adc_0 = tree->Branch(Form("%s.adc_0", prefix), &adc_0);
    b_adc_1 = tree->Branch(Form("%s.adc_1", prefix), &adc_1);
    b_adc_2 = tree->Branch(Form("%s.adc_2", prefix), &adc_2);
    b_adc_3 = tree->Branch(Form("%s.adc_3", prefix), &adc_3);
    b_adc_4 = tree->Branch(Form("%s.adc_4", prefix), &adc_4);
    b_adc_5 = tree->Branch(Form("%s.adc_5", prefix), &adc_5);
    return true;
    */
  }
  /*
  void DigGEMData_t::ClearBranches()
  {
    if(strip){
      nstrips = 0;
      strip->clear();
      module->clear();
      adc_0->clear();
      adc_1->clear();
      adc_2->clear();
      adc_3->clear();
      adc_4->clear();
      adc_5->clear();
    }
  }
  
  void DigGEMData_t::FillBranches()
  {
    if(b_nstrips){
      b_nstrips->Fill();
      b_module->Fill();
      b_strip->Fill();
      b_adc_0->Fill();
      b_adc_1->Fill();
      b_adc_2->Fill();
      b_adc_3->Fill();
      b_adc_4->Fill();
      b_adc_5->Fill();
    }
  } 
  */
}

