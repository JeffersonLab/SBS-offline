#ifndef HCAL_H
#define HCAL_H

const Int_t MAX_FADC_SAMPLES = 250;
const Int_t MAX_HCAL_MODULES = 288;
namespace hcalt {
  Double_t samps[MAX_HCAL_MODULES*MAX_FADC_SAMPLES+1000] = {0};
  Double_t nsamps[MAX_HCAL_MODULES+1000] = {0};
  Double_t row[MAX_HCAL_MODULES+1000] = {0};
  Double_t col[MAX_HCAL_MODULES+1000] = {0};
  Double_t samps_idx[MAX_HCAL_MODULES+1000] = {0};
  Double_t a[MAX_HCAL_MODULES+1000] = {0};
  Int_t ndata = 0;
  Double_t ledbit = -1;
  Double_t ledcount = 0;
};

#endif // HCAL_H
