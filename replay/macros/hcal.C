#include <TH2F.h>
#include <TChain.h>
#include <TCanvas.h>
#include <vector>
#include <iostream>
#include <TSystem.h>


const Int_t numModules = 192;
const Int_t MAX_FADC_SAMPLES = 250;


TH2F* MakeHisto(Int_t module, Int_t bins)
{
  return new TH2F(TString::Format("h%d",module),"",bins,0,bins,100,0,4095);
}

Int_t hcal(Int_t run = 931)
{
  TChain *T = new TChain("T");
  T->Add(TString::Format("rootfiles/fadc_%d.root",run));


  Double_t nhit = 0;
  Int_t numSamples[numModules] = {0};
  Int_t foundModules = 0;
  Double_t a[numModules][MAX_FADC_SAMPLES];
  T->SetBranchStatus("*",0);
  T->SetBranchStatus("Ndata.sbs.hcal.a.m*",1);
  T->SetBranchStatus("sbs.hcal.a.m*",1);
  //T->SetBranchStatus("sbs.hcal.nhit",1);
  for(Int_t m = 0; m < numModules; m++) {
    //a[m] = 0;
    T->SetBranchAddress(TString::Format("sbs.hcal.a.m%d",m),a[m]);
    T->SetBranchAddress(TString::Format("Ndata.sbs.hcal.a.m%d",m),
        &(numSamples[m]));
  }
  // Get the first entry so we can initialize the histograms
  //T->GetEntry(0);

  std::vector<TH2F*> h;
  // Create histograms
  for(Int_t m = 0; m < numModules; m++) {
  //  h.push_back(new TH2F(TString::Format("h%d",m),
  //      "",numSamples[m],0,numSamples[m],100,0,4095));
    h.push_back(0);
  }

  // Now loop through the tree and fill in the histograms
  Int_t num_entries = T->GetEntries();
  for(Int_t entry = 0; entry < num_entries; entry++) {
    // Get next entry
    T->GetEntry(entry);

    // Fill histograms
    for(Int_t m = 0; m < numModules; m++) {
      //std::cout << "numSamples[" << m << "]: " << numSamples[m] << std::endl;
      if(h[m]) {
        for(Int_t s = 0;  s < numSamples[m]; s++) {
          h[m]->Fill(s,a[m][s]);
          //if(m>0 && a[s][m] !=0 ) {
          //  std::cout << "i=" << entry << ", m: " << m << ", s: " << s << ", a: "
          //    << a[s][m] << std::endl;
          //}
        }
      } else if ( numSamples[m] > 0) {
        h[m] = MakeHisto(m,numSamples[m] );
        foundModules++;
        std::cout << "Entry: " << entry << " found channel: " <<
          m << std::endl;
      }
    }
  }

  if(foundModules == 0) {
    std::cerr << "Found no modules with samples!" << std::endl;
    return -5;
  }

  Int_t numCols = 4;
  Int_t numRows = 3;
  Int_t numPerCanvas = numCols*numRows;

  std::vector<TCanvas*> canvVector;
  Int_t countDrawn = 0;
  TCanvas *canv = 0;
  Int_t canvID = -1;
  for(Int_t m = 0; m < numModules; m++) {
    if(h[m]) {
      if(canv == 0 || countDrawn == numPerCanvas) {
        canvID++;
        canv = new TCanvas(TString::Format("canvas%d",canvID),
            TString::Format("Canvas%3d For Run %d",canvID,run),
            1200,1000);
        canv->Divide(numCols,numRows);
        canvVector.push_back(canv);
        countDrawn = 0;
      }
      countDrawn++;
      canv->cd(countDrawn);
      std::cout << "Drawing histogram for module: " << m << std::endl;
      h[m]->Draw();
    }
  }

  gSystem->mkdir("images/",kTRUE);
  for(size_t i = 0; i < canvVector.size(); i++) {
    canvVector[i]->SaveAs(TString::Format("images/canvas_%d.png",int(i)));
    canvVector[i]->SaveAs(TString::Format("images/canvas_%d.C",int(i)));
  }

  return 0;
}
