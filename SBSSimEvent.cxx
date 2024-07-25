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
//SBSSimEvent::SBSSimEvent(TTree* tree, TString experiment) {
SBSSimEvent::SBSSimEvent(TTree* tree, Exp_t experiment) {
  std::cout << "Initializing SBSSimEvent" << std::endl;
  RunID = EvtID = 0;

  fExperiment = experiment;

  //Now we need to initialize the appropriate Tree structure based on experiment:
  //We should probably use an enum or something simple to make this less clunky than doing a string comparison each time we open the file
  //or load the event:
  switch( fExperiment ){
  case kGEnRP://"genrp":
    //do nothing for now; eventually we will allocate the genrp_tree and store the pointer in the data member of this class:
    Tgmn = new gmn_tree_digitized(tree);
    Tgenrp = new genrp_tree_digitized(tree);
    
    if(Tgmn==0 || Tgenrp==0){
      std::cout << " SBSSimEvent::SBSSimEvent(): Digitized tree Tgmn / Tgenrp can't be found! Stopping the program! " << std::endl;
      exit(-1);
    }
    
    break;
  case kGEp://"gep":
    Tgep = new gep_tree_digitized(tree);
    
    if(Tgep==0){
      std::cout << " SBSSimEvent::SBSSimEvent(): Digitized tree Tgep can't be found! Stopping the program! " << std::endl;
      exit(-1);
    }
    break;
  case kSIDIS://"sidis":
    //Tsidis = new sidis_tree_digitized(tree);
    break;
  case kGMN://"gmn":
    Tgmn = new gmn_tree_digitized(tree);

    if(Tgmn==0){
      std::cout << " SBSSimEvent::SBSSimEvent(): Digitized tree Tgmn is empty! Stopping the program! " << std::endl;
      exit(-1);
    }
    
    //case //"gen":
    break;
  default:
    Tgmn = new gmn_tree_digitized(tree);
    if(Tgmn==0){
      std::cout << " SBSSimEvent::SBSSimEvent(): Digitized tree Tgmn is empty! Stopping the program! " << std::endl;
      exit(-1);
    }
    
    //Tgenrp = new genrp_tree_digitized(tree);
    break;
  }
  //Weight = 1;
  Clear();
}
/*
//_____________________________________________________________________________
SBSSimEvent::SBSSimEvent(TTree* tree, std::vector<TString> det_list) : g4sbs_tree(tree, det_list)
{
  std::cout << "Initializing SBSSimEvent" << std::endl;
  cout << det_list.size() << endl;
  RunID = EvtID = 0;
  //Weight = 1;
  Clear();
}
*/
//_____________________________________________________________________________
void SBSSimEvent::Clear( const Option_t* opt )
{
  // do nothing...
}

//_____________________________________________________________________________
void SBSSimEvent::Print( const Option_t* opt ) const
{
  //std::cout << RunID << " " << EvtID << " " << ev_sigma*ev_solang << std::endl;
}

//_____________________________________________________________________________
Int_t SBSSimEvent::GetEntry( Long64_t entry )
{
  EvtID = entry;

  //std::cout << "SBSSimEvent::GetEntry(" << entry << "): " << std::endl;
  // Read contents of entry.
  //if (!fChain) return 0;

  //fChain->Print();
  
  //int ret = fChain->GetEntry(entry);

  int ret=-1;
  //switch is not actually capable of taking strings... use a "enum"
  switch( fExperiment ){
  case kGEnRP://"genrp":
    //do nothing for now; eventually we will invoke the "GetEntry" methods of the various classes:
    ret = Tgmn->GetEntry(entry);
    ret = Tgenrp->GetEntry(entry);
    break;
  case kGEp://"gep":
    ret = Tgep->GetEntry(entry);
    break;
  case kSIDIS://"sidis":
    //ret = Tsidis->GetEntry(entry);
    break;
  case kGMN://"gmn":
    //case "gen":
    ret = Tgmn->GetEntry(entry);
    break;
  default:
    ret = Tgmn->GetEntry(entry);
    //ret = Tgenrp->GetEntry(entry);
    break;
  }
  return ret;
}

//-----------------------------------------------------------------------------
ClassImp(SBSSimEvent)

