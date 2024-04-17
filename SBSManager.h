#ifndef SBSMANAGER_H
#define SBSMANAGER_H

///////////////////////////////////////////////////////////////////////////////
//  - 2021-02-25 Juan Carlos Cornejo <cornejo@jlab.org>
//  * Create this global (singleton) class

#include "TString.h"
#include "TObject.h"
#include "THaCrateMap.h"

class SBSManager : public TObject {
public:
  SBSManager();
  static SBSManager *GetInstance();
  void SetDefaultCrateMapName(const char* name);

  Decoder::THaCrateMap* GetCrateMap( Long64_t tloc );

private:
  std::unique_ptr<Decoder::THaCrateMap> fCrateMap;
  TString fCrateMapName;
  //FIXME: With analyzer >= 1.8, use the time stored in THaCrateMap
  Long64_t fCrateMapInitTime;

  static std::unique_ptr<SBSManager> fManager;

  ClassDefNV(SBSManager,0) ///< Base class for SBS Manage
};

#endif//SBSMANAGER_H
