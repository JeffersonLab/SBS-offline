#ifndef __SBSSimEvent_h
#define __SBSSimEvent_h

//#include "gmn_dig_tree.h"
#include "g4sbs_tree.h"

class TTree;

class SBSSimEvent : public g4sbs_tree {
 public:
  SBSSimEvent();                 // Default constructor, for ROOT I/O
  SBSSimEvent(TTree* tree);//, std::vector<TString> det_list);
  
  virtual ~SBSSimEvent(){};
  
  virtual Int_t    GetEntry(Long64_t entry);
  virtual void Clear( const Option_t* opt="" );
  virtual void Print( const Option_t* opt="" ) const;

  ULong64_t RunID;
  ULong64_t EvtID;

  ClassDef(SBSSimEvent, 1) // Simulated data for one event
};

#endif
