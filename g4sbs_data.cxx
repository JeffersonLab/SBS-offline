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
  }

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
  }

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
  }
}

