//*-- Author :    Ole Hansen (ole@jlab.org)    9-Dec-2011

/////////////////////////////////////////////////////////////////////
//
//   SBSSimFile
//
//   Interface to an input file with simulated SoLID spectrometer data
//
//   Takes raw digitized simulation data from ROOT input file and
//   uses them to fill a SBSSimEvent object. A pointer to the event
//   object is available via GetEvBuffer() for use by the decoder.
//
/////////////////////////////////////////////////////////////////////

#include "THaGlobals.h"
#include "THaApparatus.h"
#include "SBSSimFile.h"
#include "SBSSimEvent.h"
#include "SBSBBShower.h"
#include "SBSBBTotalShower.h"
//#include "evio.h"     // for S_SUCCESS
// We really only need S_SUCCESS from evio.h, so why not just define
// it ourselves so we don't have to pull the whole header file.
#ifndef S_SUCCESS
#define S_SUCCESS 0
#define S_FAILURE -1
#endif

#include "TFile.h"
#include "TTree.h"
#include "TError.h"
#include "TClonesArray.h"
#include "TString.h"
#include "TMath.h"

#include <cstring>
#include <libgen.h>    // for POSIX basename()
#include <cstdlib>
#include <iostream>

using namespace std;

//-----------------------------------------------------------------------------
SBSSimFile::SBSSimFile(const char* filename, const char *experiment, const char* description) :
  THaRunBase(description), fROOTFileName(filename), //fExperiment(experiment), 
  fROOTFile(0), fTree(0), 
  fEvent(0), fNEntries(0), fEntry(0), fVerbose(0)
{
  // Constructor

  // Use default if no file name given
  if( fROOTFileName.IsNull() ) {
    Info( __FUNCTION__, "Using default input file MCdata.root" );
    fROOTFileName = "MCdata.root";
  }
  /*
  fDetList.clear();
  TIter aiter(gHaApps);
  THaApparatus* app = 0;
  while( (app=(THaApparatus*)aiter()) ){
    TList* listdet = app->GetDetectors();
    TIter diter(listdet);
    TObject* det = 0;
    while( (det=(TObject*)diter()) ){
      fDetList.push_back(det->GetName() );
      if(strcmp(app->GetDetector(det->GetName())->GetClassName(),"SBSBBTotalShower")==0){
	SBSBBTotalShower* TS = (SBSBBTotalShower*)app->GetDetector(det->GetName());
	fDetList.push_back(TS->GetShower()->GetName());
	fDetList.push_back(TS->GetPreShower()->GetName());
	
      }
    }
  }
  */
  GetExperiment(experiment);
  
  // fValidExperiments.insert("gmn");
  // fValidExperiments.insert("genrp");
  // fValidExperiments.insert("gep");
  // fValidExperiments.insert("sidis");
  
  fValidExperiments.insert(kGMN);
  fValidExperiments.insert(kGEnRP);
  fValidExperiments.insert(kGEp);
  fValidExperiments.insert(kSIDIS);

  if( fValidExperiments.find( fExperiment ) == fValidExperiments.end() ){ //This is not a valid experiment. Default to gmn and print a warning:
    TString fWarn;
    fWarn.Form( "SBSSimFile(%s, %s, %s)", filename, experiment, description );
    
    Warning(Here(fWarn.Data()), "Invalid simulated experiment choice... defaulting to gmn");

    fExperiment = kGMN;//"gmn";
  }
  
}

/*
//-----------------------------------------------------------------------------
SBSSimFile::SBSSimFile(const char* filename, const char* description, std::vector<TString> det_list) :
  THaRunBase(description), fROOTFileName(filename), fROOTFile(0), fTree(0), 
  fEvent(0), fNEntries(0), fEntry(0), fVerbose(0)
{
  // Constructor

  // Use default if no file name given
  if( fROOTFileName.IsNull() ) {
    Info( __FUNCTION__, "Using default input file MCdata.root" );
    fROOTFileName = "MCdata.root";
  }
  for(size_t i = 0; i<det_list.size(); i++)fDetList.push_back(det_list[i]);
}
*/
//-----------------------------------------------------------------------------
SBSSimFile::SBSSimFile(const SBSSimFile &run)
  : THaRunBase(run), fROOTFileName(run.fROOTFileName), 
    fROOTFile(0), fTree(0), fEvent(0), fNEntries(0), fEntry(0), fVerbose(0)
{
}

//-----------------------------------------------------------------------------
SBSSimFile& SBSSimFile::operator=(const THaRunBase& rhs)
{
  if (this != &rhs) {
    THaRunBase::operator=(rhs);
    if( rhs.InheritsFrom("SBSSimFile") )
      fROOTFileName = static_cast<const SBSSimFile&>(rhs).fROOTFileName;
    fROOTFile = 0;
    fTree = 0;
    fEvent = 0;
    fNEntries = fEntry = 0;
  }
  return *this;
}

//_____________________________________________________________________________
Int_t SBSSimFile::Compare( const TObject* obj ) const
{
  // Compare two SBSSimFiles. They are different if either their number
  // or their input file name differs.

  if (this == obj) return 0;
  const THaRunBase* rhs = dynamic_cast<const THaRunBase*>(obj);
  if( !rhs ) return -1;
  // operator< compares fNumber
  if( *this < *rhs )       return -1;
  else if( *rhs < *this )  return  1;
  const SBSSimFile* rhsr = dynamic_cast<const SBSSimFile*>(rhs);
  if( !rhsr ) return 0;
  if( fROOTFileName < rhsr->fROOTFileName ) return -1;
  else if( rhsr->fROOTFileName < fROOTFileName ) return 1;
  return 0;
}

//-----------------------------------------------------------------------------
Int_t SBSSimFile::Init()
{
  // Initialize the run. Sets run date, reads run database etc.

  // TODO: get date from MC production file?
  //fDate.Set(2012,1,1,0,0,0);
  fAssumeDate = kTRUE;
  fDataSet |= kDate;

  Int_t ret = THaRunBase::Init();
  if( !ret ) {
    char* s = strdup(fROOTFileName);
    fName = basename(s);
    free(s);
    fNumber = 1;
  }
  return ret;
}

//-----------------------------------------------------------------------------
Int_t SBSSimFile::ReadDatabase()
{
  
  static const char* const here = "ReadDatabase";
  
  THaRunBase::ReadDatabase();
  
  FILE* f = Podd::OpenDBFile("run", fDate, "SBSSimFile::ReadDatabase", "r", 1);
  if( !f )  return -1;
  TString expt;
  DBRequest request[] = {
    { "experiment",  &expt, kTString },
    { nullptr }
  };
  Int_t err = THaAnalysisObject::LoadDB( f, fDate, request, "");
  if(fVerbose>0)std::cout << " expt " << expt.Data() << " DBreq status: " << err << std::endl; 
  
  fclose(f);
  if( err ){
    std::cout << "Warning: can't read experiment flag in MC run database! " << std::endl 
	      << "Defaulting to GMN - if you're not running GMN simulation, it will likely cause issues." << std::endl 
	      << "Fix database in this case. " << std::endl;
    //default!!!
    fExperiment = kGMN;
  }else{
    GetExperiment(expt.Data());
  }
  //return err;
  
  fDBRead = true;
  return READ_OK;
}

//-----------------------------------------------------------------------------
void SBSSimFile::GetExperiment(const char *experiment)
{
  if(fVerbose>1)cout << "using experiment configuration: " << experiment << endl;
  
  if(strcmp(experiment,"gmn")==0 || strcmp(experiment,"gen")==0){
    //(experiment=="gmn" || experiment=="gen")//{
    fExperiment = kGMN;
  }
  if(strcmp(experiment,"genrp")==0){
    //(experiment=="genrp")//{
    fExperiment = kGEnRP;
  }
  if(strcmp(experiment,"gep")==0){
    //(experiment=="gep")//{
    fExperiment = kGEp;
  }
  if(strcmp(experiment,"sidis")==0){
    //(experiment=="sidis")//{
    fExperiment = kSIDIS;
  }
  
}

//-----------------------------------------------------------------------------
Int_t SBSSimFile::Open()
{
  ReadDatabase();
  // Open ROOT input file
  if(fVerbose>0)std::cout << "SBSSimFile::Open(): initialize file with experiment: " << fExperiment << std::endl;
  
  fROOTFile = new TFile(fROOTFileName, "READ", "SBS MC data");
  if( !fROOTFile || fROOTFile->IsZombie() ) {
    Error( __FUNCTION__, "Cannot open input file %s", fROOTFileName.Data() );
    Close();
    return -1;
  }

  fTree = static_cast<TTree*>( fROOTFile->Get(treeName) );
  if( !fTree ) {
    Error( __FUNCTION__, "Tree %s does not exist in the input file.",
	   treeName );
    Close();
    return -2;
  }

  //  fTree->SetBranchStatus("*", kFALSE);

  /*
  // Set up reading of the event branch
  delete fEvent; fEvent = 0;

  // UInt_t found = 0;
  // fTree->SetBranchStatus( eventBranchName, kTRUE, &found );
  // if( found > 0 ) {
  //   fTree->SetBranchAddress( eventBranchName, &fEvent );
  // }
  TBranch* br = fTree->GetBranch(treeName);
  if( br ) {
    br->SetAddress(&fEvent);
  }
  else {
    Error( __FUNCTION__, "No event branch \"%s\" in the input tree.",
	   treeName );
    Close();
    return -3;
  }

  */
  
  fNEntries = fTree->GetEntries();
  fEntry = 0;
  
  //fEvent is actually the tree!!!
  // and if it works it turns out to make the thing actually much simpler.
  delete fEvent; fEvent = 0;// really needed anymore ?
  fEvent = new SBSSimEvent(fTree, fExperiment);
  //cout << fDetList.size() << endl;
  //fTree->Print();
  //fEvent = new SBSSimEvent(fTree, fDetList);
  
  fOpened = kTRUE;
  return 0;
}

//-----------------------------------------------------------------------------
Int_t SBSSimFile::Close()
{
  if(fVerbose>0)std::cout << " SBSSimFile::Close(): closing file and destroy previously configured SBSSimEvent " << std::endl;
  
  delete fTree; fTree = 0;
  if (fROOTFile) {
    fROOTFile->Close();
    delete fROOTFile; fROOTFile = 0;
  }
  delete fEvent; fEvent = 0;
  fOpened = kFALSE;
  return 0;
}

//-----------------------------------------------------------------------------
Int_t SBSSimFile::ReadEvent()
{
  // Read next event from ROOT file
  if(fVerbose>1 || fEntry%(fNEntries/100)==0 )
    std::cout << " SBSSimFile::ReadEvent() -> fEntry = " 
	      << fEntry << " / " << fNEntries
	      << std::endl;
  if( fEntry >= fNEntries )
    return EOF;

  Int_t ret;
  if( !IsOpen() ) {
    std::cout << " Warning: file not open when ReadEvent() called... unintended consequences?" << std::endl;
    ret = Open();
    if( ret ) return ret;
  }

  // Clear the event to get rid of anything still hanging around
  if( fEvent ) fEvent->Clear();

  // Read input file
  //ret = fTree->GetEntry(fEntry++);

  //std::cout << "trying to load event " << fEntry << std::endl;
  
  ret = fEvent->GetEntry(fEntry++);
  if( ret == 0 )
    return EOF;
  if( ret < 0 )
    return -128;  // CODA_ERR
  return S_SUCCESS;
}

//-----------------------------------------------------------------------------
#if ANALYZER_VERSION_CODE >= ANALYZER_VERSION(1,6,0)
const UInt_t *SBSSimFile::GetEvBuffer() const
#else
const  Int_t *SBSSimFile::GetEvBuffer() const
#endif
{
  if( !IsOpen() ) return 0;
  // EPAF: this is the "reinterpret_cast" that is essential.
  // It transforms the "tree" into a (stl?) vector of integers
#if ANALYZER_VERSION_CODE >= ANALYZER_VERSION(1,6,0)
  return reinterpret_cast<UInt_t*>(fEvent);
#else
  return reinterpret_cast<Int_t*>(fEvent);
#endif
}

//-----------------------------------------------------------------------------
SBSSimFile::~SBSSimFile()
{
  if( IsOpen() )
    Close();
}

//_____________________________________________________________________________
void SBSSimFile::Print( Option_t* opt ) const
{
  // Print run info and status

  TString sopt(opt);
  bool do_header = sopt.Contains("start", TString::kIgnoreCase);
  if( sopt.IsNull() || do_header ) {
    THaRunBase::Print("STARTINFO");
  }

  if( sopt.IsNull() || sopt.Contains("status", TString::kIgnoreCase) ) {
    cout << "Analyzed events:       " << fNumAnalyzed << endl;
    cout << "Initialized:           " << fIsInit << endl;
    cout << "Opened:                " << fOpened << endl;
  }

  if( fVerbose > 3 && fTree ) fTree->Print();

}

//-----------------------------------------------------------------------------
ClassImp(SBSSimFile)
