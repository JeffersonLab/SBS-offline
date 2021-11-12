#ifndef SBSAnalyzer_h_
#define SBSAnalyzer_h_

#include "THaAnalyzer.h"

class SBSAnalyzer : public THaAnalyzer {
 public:
  SBSAnalyzer();
  virtual ~SBSAnalyzer();
  
 protected:
  //functions we want to override
  virtual void   PrintCutSummary() const;
  
 private:
  THaAnalyzer( const SBSAnalyzer& );
  THaAnalyzer& operator=( const SBSAnalyzer& );

  ClassDef(SBSAnalyzer, 0)
}
