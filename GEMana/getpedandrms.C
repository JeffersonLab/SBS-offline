#include "TH1F.h"

void getpedandrms(){
    TChain *T = new TChain("T");
    T->Add("Afile.root");
    const char *prefix = "sbs.gems.y1";

    FILE *output = fopen("pedrms.dat", "w");

    int minadc = 0;
    int maxadc = 1000;
    int nadc = 100;

    // Get strip numbers - we'll assume contiguous

    Double_t nch;
    Double_t strip[50000];
    Double_t adcval[50000];


    int i,j;

    int minstrip = 1e9;
    int maxstrip = -1e9;

    T->SetBranchAddress(Form("%s.nch", prefix), &nch);
    T->SetBranchAddress(Form("%s.strip", prefix), strip);
    T->SetBranchAddress(Form("%s.adc0", prefix), adcval);

    for( i = 0; i < T->GetEntries(); i++ ){
        T->GetEntry(i);
        for( j = 0; j < (int) nch; j++ ){
            if( strip[j] < minstrip ){
                minstrip = (int) strip[j];
            }
            if( strip[j] > maxstrip ){
                maxstrip = (int) strip[j];
            }
        }

    }

    printf("min/max/strip %d %d\n", minstrip, maxstrip);

    TH1F **h;

    printf("Allocating %d\n", (maxstrip-minstrip)+1);
    h = new TH1F *[(maxstrip-minstrip)+1];
    for( i = 0; i < (maxstrip-minstrip)+1; i++ ){
        h[i] = new TH1F(Form("strip_%d", i+minstrip),Form("strip_%d", i+minstrip), nadc, minadc, maxadc);
    }

    for( i = 0; i < T->GetEntries(); i++ ){
        T->GetEntry(i);
        for( j = 0; j < (int) nch; j++ ){
            int hidx = ((int) strip[j])-minstrip;
            h[hidx]->Fill(adcval[j]);
        }
    }

    fprintf(output, "%s.ped =", prefix);
    for( i = 0; i < (maxstrip-minstrip)+1; i++ ){
        fprintf(output, "%d %f ", i+minstrip, h[i]->GetMean());
//        fprintf(output, "%d %f ", i+minstrip, 0.0);
    }
    fprintf(output, "\n\n%s.rms =", prefix);
    for( i = 0; i < (maxstrip-minstrip)+1; i++ ){
        fprintf(output, "%d %f ", i+minstrip, h[i]->GetRMS());
//        fprintf(output, "%d %f ", i+minstrip, 0.0);
    }
    fclose(output);

}
