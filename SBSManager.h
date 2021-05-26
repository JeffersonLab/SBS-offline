#ifndef SBSMANAGER_H
#define SBSMANAGER_H

///////////////////////////////////////////////////////////////////////////////
//  - 2021-02-25 Juan Carlos Cornejo <cornejo@jlab.org>
//  * Create this global (instanton) class

#include <TString.h>
#include <TObject.h>


namespace Decoder {
  class THaCrateMap;
}

class SBSManager : public TObject {
public:
  SBSManager();
  virtual ~SBSManager();
  static SBSManager *GetInstance();
  static void SetDefaultCrateMapName(const char* name);

  Decoder::THaCrateMap* GetCrateMap();

private:
  Decoder::THaCrateMap *fCrateMap;
  static TString fCrateMapName;
  static SBSManager *fManager;

  ClassDef(SBSManager,0) ///< Base class for SBS Manage
};

#endif//SBSMANAGER_H
