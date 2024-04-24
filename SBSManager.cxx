#include "SBSManager.h"

std::unique_ptr<SBSManager> SBSManager::fManager = nullptr;

///////////////////////////////////////////////////////////////////////////////
// Get an instance of the singleton
SBSManager* SBSManager::GetInstance()
{
  if(!fManager) {
    fManager.reset(new SBSManager);
  }
  return fManager.get();
}
  
///////////////////////////////////////////////////////////////////////////////
// Get the instance of the global crate map
// Initialize the map for the given Unix time 'tloc'. If the time changed since
// the previous initialization, then re-initialize the map.
Decoder::THaCrateMap* SBSManager::GetCrateMap( Long64_t tloc )
{
  if(!fCrateMap) {
    fCrateMap.reset(new Decoder::THaCrateMap(fCrateMapName));
    fCrateMapInitTime = -1;
  }
  if(fCrateMapInitTime != tloc) {
    fCrateMap->init(tloc);
    fCrateMapInitTime = tloc;
  }
  return fCrateMap.get();
}

///////////////////////////////////////////////////////////////////////////////
void SBSManager::SetDefaultCrateMapName(const char* name)
{
  if( fCrateMapName != name ) {
    // Force re-init
    fCrateMap.reset();
  }
  fCrateMapName = name;
}

///////////////////////////////////////////////////////////////////////////////
// Basic constructor
SBSManager::SBSManager()
  : fCrateMap{}
  , fCrateMapName{"cratemap"}
  , fCrateMapInitTime{-1}
{}
