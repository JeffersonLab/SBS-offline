#include <TH2F.h>
#include <TChain.h>
#include <TCanvas.h>
#include <vector>
#include <iostream>
#include <TSystem.h>


const Int_t numSamples = 10;

Int_t hcal(Int_t run = 931)
{
  TChain *T = new TChain("T");
  T->Add(TString::Format("rootfiles/fadc_%d.root",run));

  Double_t nhit = 0;
  Int_t numModule[numSamples] = {10};
  Double_t a[numSamples][300];
  T->SetBranchStatus("*",0);
  T->SetBranchStatus("Ndata.sbs.hcal.a*",1);
  T->SetBranchStatus("sbs.hcal.a*",1);
  T->SetBranchStatus("sbs.hcal.nhit",1);
  for(Int_t s = 0; s < numSamples; s++) {
    T->SetBranchAddress("sbs.hcal.nhit",&nhit);
    T->SetBranchAddress(TString::Format("Ndata.sbs.hcal.a%d",s),
        &(numModule[s]));
    T->SetBranchAddress(TString::Format("sbs.hcal.a%d",s),a[s]);
  }
  T->GetEntry(0);
  std::cout << "numModule: " << numModule[0] << std::endl;

  // Create histograms
  std::vector<TH2F*> h;
  for(Int_t m = 0; m < numModule[0]; m++) {
    h.push_back(new TH2F(TString::Format("h%d",m),
        "",numSamples,0,numSamples,20,0,4095));
  }

  Int_t num_entries = T->GetEntries()*0.05;
  for(Int_t entry = 0; entry < num_entries; entry++) {
    // Get next entry
    T->GetEntry(entry);

    // Fill histograms
    for(Int_t s = 0;  s < numSamples; s++) {
      for(Int_t m = 0; m < numModule[s]; m++) {
        h[m]->Fill(s,a[s][m]);
        //std::cout << "i=" << entry << ", m: " << m << ", s: " << s << ", a: "
        //  << a[s][m] << std::endl;
      }
    }
  }


  Int_t numCols = 4;
  Int_t numRows = 3;
  Int_t numPerCanvas = numCols*numRows;

  Int_t numCanvas = numModule[0]/numPerCanvas;
  if ( numModule[0] % numPerCanvas != 0)
    numCanvas++;
  std::cout << "numCanvas: " << numCanvas << std::endl;
  std::cout << "numCols: " << numCols << std::endl;
  std::cout << "numRows: " << numRows << std::endl;
  std::cout << "numModule: " << numModule[0] << std::endl;

  std::vector<TCanvas*> canvas;
  for(Int_t c = 0; c < numCanvas; c++) {
    TCanvas *canv = new TCanvas(TString::Format("canvas%d",c),
          TString::Format("canvas%d",c),
          1200,1000);
    canv->Divide(numCols,numRows);
    canvas.push_back(canv);
  }

  Int_t cc = 0; // Current canvas number
  for(Int_t m = 0; m < numModule[0]; m++) {
    Int_t cCanv = m/(numPerCanvas);
    Int_t cPad  = m%(numPerCanvas);
    std::cerr << "Canvas: " << cCanv << ", Pad: " << cPad << std::endl;
    canvas[cCanv]->cd(cPad+1);
    h[m]->Draw();
  }

  gSystem->mkdir("images/",kTRUE);
  for(size_t i = 0; i < canvas.size(); i++) {
    canvas[i]->SaveAs(TString::Format("images/canvas_%d.png",int(i)));
  }

  return 0;
}
