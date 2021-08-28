R__ADD_INCLUDE_PATH($SBSDEVEL/install/include)
R__ADD_LIBRARY_PATH($SBSDEVEL/install/lib64)
R__ADD_LIBRARY_PATH($SBSDEVEL/install/lib)
R__LOAD_LIBRARY(libsbs.so)

#if !defined(__CLING__) || defined(__ROOTCLING__)
#include <iostream>

#include "TSystem.h"
#include "TList.h"
#include "THaRun.h"
#include "THaEvent.h"
#include "THaAnalyzer.h"
#include "THaApparatus.h"

//#include "SBSGEMStand.h"
//#include "SBSBigBite.h"
#include "SBSEArm.h"
#include "SBSRPBeamSideHodo.h"
#include "SBSRPFarSideHodo.h"
#include "SBSCHAnalyzer.h"
#endif

void replay_rpgentest(Int_t runnum = 288, Int_t lastEvent = -1){
  
 gSystem->Load("libsbs.so");
 SBSRPBeamSideHodo *rpbeamsidehodo = new SBSRPBeamSideHodo("RPBeamSideHodo","RPBEAMSIDEHODO");
 SBSRPFarSideHodo *rpfarsidehodo = new SBSRPFarSideHodo("RPFarSideHodo","RPFARSIDEHODO");
 SBSCHAnalyzer *chanalyzer = new SBSCHAnalyzer("CHAnalyzer","CHANALYZER");
 SBSHCal *hcal = new SBSHCal("hcal","HCAL");
 
 SBSEArm *harm = new SBSEArm("sbs","Hadron Arm with RPBeamSideHodo,RPFarSideHodo and CHAnalyzer");
 harm->AddDetector(rpbeamsidehodo);
 harm->AddDetector(rpfarsidehodo);
 harm->AddDetector(chanalyzer);
 harm->AddDetector(hcal);

 THaAnalyzer* analyzer = new THaAnalyzer;
 gHaApps->Add(harm);

 THaEvent* event = new THaEvent;

 THaRun* run = new THaRun(TString::Format("data/fadc_%d.dat.0",runnum) );
  run->SetDate(TDatime());

  analyzer->SetVerbosity(0);
  analyzer->SetOdefFile("output_rpgentest.def");
  analyzer->SetEvent(event);
  analyzer->SetOutFile(TString::Format("rootfiles/rpgen_%d.root",runnum));
  analyzer->SetSummaryFile("sbs_rpgentest.log");
  
  analyzer->Process(run);
 
}
