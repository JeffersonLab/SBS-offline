#include "SBSAnalyzer.h"


// Pointer to single instance of this object
SBSAnalyzer* SBSAnalyzer::fgAnalyzer = nullptr;

//_____________________________________________________________________________
SBSAnalyzer::SBSAnalyzer() : THaAnalyzer()
{
}

//_____________________________________________________________________________
SBSAnalyzer::~SBSAnalyzer()
{
  // Destructor.

  SBSAnalyzer::Close();
  DeleteContainer(fPostProcess);
  DeleteContainer(fEvtHandlers);
  DeleteContainer(fInterStage);
  delete fExtra; fExtra = nullptr;
  delete fBench;
  if( fgAnalyzer == this )
    fgAnalyzer = nullptr;
}

//_____________________________________________________________________________
void SBSAnalyzer::PrintCutSummary() const
{
  // Print summary of cuts etc.
  // Only print to screen if fVerbose>1, but always print to
  // the summary file if a summary file is requested.

  if( gHaCuts->GetSize() > 0 ) {
    cout << "Cut summary:" << endl;
    if( fVerbose>1 )
      gHaCuts->Print("STATS");
    if( fSummaryFileName.Length() > 0 ) {
      ofstream ostr(fSummaryFileName);
      if( ostr ) {
	// Write to file via cout
	streambuf* cout_buf = cout.rdbuf();
	cout.rdbuf(ostr.rdbuf());
	TDatime now;
	cout << "Cut Summary for run " << fRun->GetNumber()
	     << " completed " << now.AsString()
	     << endl << endl;
	gHaCuts->Print("STATS");
	cout.rdbuf(cout_buf);
	ostr.close();
      }
    }
  }
}

//_____________________________________________________________________________

ClassImp(THaAnalyzer)
