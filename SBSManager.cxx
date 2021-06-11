#include "SBSManager.h"
#include "THaCrateMap.h"
#include <ctime>


SBSManager* SBSManager::fManager = 0;
TString SBSManager::fCrateMapName = "cratemap";

///////////////////////////////////////////////////////////////////////////////
// Get an instance of the singleton
SBSManager* SBSManager::GetInstance()
{
  if(!fManager) {
    fManager = new SBSManager();
  }
  return fManager;
}
  
///////////////////////////////////////////////////////////////////////////////
// Get the instance of the global crate map
Decoder::THaCrateMap* SBSManager::GetCrateMap()
{
  if(!fCrateMap) {
    fCrateMap = new Decoder::THaCrateMap(fCrateMapName);
    // Initialize to the default run time (which is right now)
    fCrateMap->init(time(0));
  }

  return fCrateMap;
}

///////////////////////////////////////////////////////////////////////////////
void SBSManager::SetDefaultCrateMapName(const char* name)
{
  fCrateMapName = name;
}

///////////////////////////////////////////////////////////////////////////////
// Basic constructor
SBSManager::SBSManager() : fCrateMap(0) {
}

///////////////////////////////////////////////////////////////////////////////
// Destructor
SBSManager::~SBSManager() {
  if(fCrateMap)
    delete fCrateMap;
};
