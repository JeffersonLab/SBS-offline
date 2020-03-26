#include <iostream>
#include "SBSSimEvent.h"
#include "TTree.h"

// ----------------------------------------------- 
// class SBSSimEvent: encapsulation of digsim_tree 
//
//_____________________________________________________________________________
SBSSimEvent::SBSSimEvent() : digsim_tree()
{
  std::cout << "Initializing TSBSSimEvent" << std::endl;
  RunID = EvtID = 0;
  Weight = 1;
  Clear();
}

//_____________________________________________________________________________
SBSSimEvent::SBSSimEvent(TTree* tree) : digsim_tree(tree)
{
  std::cout << "Initializing SBSSimEvent" << std::endl;
  RunID = EvtID = 0;
  Weight = 1;
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
  std::cout << RunID << " " << RunID << " " << Weight << std::endl;
}

//-----------------------------------------------------------------------------
ClassImp(SBSSimEvent)

