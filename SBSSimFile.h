/////////////////////////////////////////////////////////////////////
//
//   SBSSimFile
//
//   Interface to an input file with simulated SoLID data
//
//   Ole Hansen (ole@jlab.org)  December 2011
//
/////////////////////////////////////////////////////////////////////

#ifndef __SBSSimFile_h
#define __SBSSimFile_h

#include "THaRunBase.h"
#include "TString.h"
#include "digsim_tree.h"
#include "ha_compiledata.h"

class TFile;
class TTree;
class TBranch;
//class digsim_tree;
//class SBSSimEvent;

//Class SBSSimEvent encapsulates digsim_tree
class SBSSimEvent : public digsim_tree {
public:
  SBSSimEvent();                 // Default constructor, for ROOT I/O
  SBSSimEvent(TTree* tree);
  
  virtual ~SBSSimEvent(){};
  
  virtual void Clear( const Option_t* opt="" );
  virtual void Print( const Option_t* opt="" ) const;

  ClassDef(SBSSimEvent, 1) // Simulated data for one event
};

class SBSSimFile : public THaRunBase {
 public:
  SBSSimFile(const char* filename, const char* description = "");
  SBSSimFile(const SBSSimFile &run);
  virtual ~SBSSimFile();
  virtual SBSSimFile &operator=(const THaRunBase &rhs);
  // for ROOT RTTI
  SBSSimFile() : fROOTFile(0), fTree(0), fEvent(0), fEntry(0) {}

  virtual void  Print( Option_t* opt="" ) const;

  Int_t         Close();
  virtual Int_t Compare( const TObject* obj ) const;
#if ANALYZER_VERSION_CODE >= 67072  // ANALYZER_VERSION(1,6,0)
  const UInt_t* GetEvBuffer() const;
#else
  const  Int_t* GetEvBuffer() const;
#endif
  Int_t         Init();
  const char*   GetFileName() const { return fROOTFileName.Data(); }
  Int_t         Open();
  Int_t         ReadEvent();
  void          SetFileName( const char* name ) { fROOTFileName = name; }

 protected:
  virtual Int_t ReadDatabase() {return 0;}

  TString fROOTFileName;  //  Name of input file
  TFile* fROOTFile;       //! Input ROOT file
  TTree* fTree;           //! Input Tree with simulation data
  SBSSimEvent* fEvent;   //! Current event

  ULong64_t fNEntries;    //! Number of entries in tree
  ULong64_t fEntry;       //! Current entry number

  ClassDef(SBSSimFile,1) // Interface to input file with simulated SoLID data
};

#endif
