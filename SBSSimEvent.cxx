#include <iostream>
#include "SBSSimEvent.h"
#include "TTree.h"

// ----------------------------------------------- 
// class SBSSimEvent: encapsulation of digsim_tree 
//
//_____________________________________________________________________________
SBSSimEvent::SBSSimEvent() : gmn_dig_tree()
{
  std::cout << "Initializing SBSSimEvent" << std::endl;
  RunID = EvtID = 0;
  //Weight = 1;
  Clear();
}

//_____________________________________________________________________________
SBSSimEvent::SBSSimEvent(TTree* tree) : gmn_dig_tree(tree)
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

//-----------------------------------------------------------------------------
ClassImp(SBSSimEvent)

