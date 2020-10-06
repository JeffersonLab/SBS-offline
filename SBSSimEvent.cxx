#include <iostream>
#include "SBSSimEvent.h"
#include "TTree.h"

// ----------------------------------------------- 
// class SBSSimEvent: encapsulation of g4sbs_tree 
//
/*
//_____________________________________________________________________________
SBSSimEvent::SBSSimEvent() : g4sbs_tree()
{
  std::cout << "Initializing SBSSimEvent" << std::endl;
  RunID = EvtID = 0;
  //Weight = 1;
  Clear();
}
*/
//_____________________________________________________________________________
SBSSimEvent::SBSSimEvent(TTree* tree, std::vector<TString> det_list) : g4sbs_tree(tree, det_list)
{
  std::cout << "Initializing SBSSimEvent" << std::endl;
  RunID = EvtID = 0;
  //Weight = 1;
  Clear();
}

//_____________________________________________________________________________
void SBSSimEvent::Clear( const Option_t* opt )
{
  // do nothing...
}

//_____________________________________________________________________________
void SBSSimEvent::Print( const Option_t* opt ) const
{
  std::cout << RunID << " " << EvtID << " " << ev_sigma*ev_solang << std::endl;
}

//_____________________________________________________________________________
Int_t SBSSimEvent::GetEntry( Long64_t entry )
{
  EvtID = entry;
  // Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}

//-----------------------------------------------------------------------------
ClassImp(SBSSimEvent)

