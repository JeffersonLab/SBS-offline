#ifndef __SBSSimEvent_h
#define __SBSSimEvent_h

#include "digsim_tree.h"

class TTree;

class SBSSimEvent : public digsim_tree {
public:
  SBSSimEvent();                 // Default constructor, for ROOT I/O
  SBSSimEvent(TTree* tree);
  
  virtual ~SBSSimEvent(){};
  
  virtual void Clear( const Option_t* opt="" );
  virtual void Print( const Option_t* opt="" ) const;

  ClassDef(SBSSimEvent, 1) // Simulated data for one event
};

#endif