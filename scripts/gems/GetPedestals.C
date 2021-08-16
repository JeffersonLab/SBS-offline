#include "TH2F.h"
#include "TFile.h"
#include "TH1F.h"
#include <iostream>
#include <fstream>
#include "TString.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TPad.h"
#include "TStyle.h"
#include "TClonesArray.h"

void GetPedestals( const char *rootfilename, int nmodules=20, const char *prefix="sbs.uvagem.m", const char *outfilename="pedtemp.root", int plotflag=0 ){

  gStyle->SetOptFit();
  
  TCanvas *c1 = new TCanvas("c1","c1",1600,1200);

  c1->Divide(8,8,0.001,0.001);
  
  TFile *fin = new TFile( rootfilename, "READ" );

  //For "results" plots:
  TFile *fout = new TFile( outfilename, "RECREATE" );

  fout->cd();
  
  TH2F *hPedRMSU_distribution_vs_module = new TH2F("hPedRMSU_distribution_vs_module", "", nmodules, -0.5, nmodules-0.5, 200, 0, 200 );
  TH2F *hPedSigmaU_distribution_vs_module = new TH2F("hPedSigmaU_distribution_vs_module", "", nmodules, -0.5, nmodules-0.5, 200, 0, 200 );

  TH2F *hPedRMSV_distribution_vs_module = new TH2F("hPedRMSV_distribution_vs_module", "", nmodules, -0.5, nmodules-0.5, 200, 0, 200 );
  TH2F *hPedSigmaV_distribution_vs_module = new TH2F("hPedSigmaV_distribution_vs_module", "", nmodules, -0.5, nmodules-0.5, 200, 0, 200 );

  TH2F *hPedMeanU_distribution_vs_module = new TH2F("hPedMeanU_distribution_vs_module", "", nmodules, -0.5, nmodules-0.5, 250, -250, 250 );
  TH2F *hPedMeanV_distribution_vs_module = new TH2F("hPedMeanV_distribution_vs_module", "", nmodules, -0.5, nmodules-0.5, 250, -250, 250 );
  TH2F *hPedFitMeanU_distribution_vs_module = new TH2F("hPedFitMeanU_distribution_vs_module", "", nmodules, -0.5, nmodules-0.5, 250, -250, 250 );
  TH2F *hPedFitMeanV_distribution_vs_module = new TH2F("hPedFitMeanV_distribution_vs_module", "", nmodules, -0.5, nmodules-0.5, 250, -250, 250 );
  //TH2D *hPedestal_distribution_vs_module = new TH2D("hPedNoise_distribution_vs_module", "", nmodules -0.5, nmodules-0.5, 500,-250,250);

  TString dbfilename = outfilename;
  dbfilename.ReplaceAll(".root",".dat");

  ofstream outfile(dbfilename.Data());

  //TString detname = prefix;
  
  for( int i=0; i<nmodules; i++ ){
    TString histname;
    TString header;
    
    TH2F *hpedu, *hpedv; //common-mode corrected only
    TH2F *hcmmu, *hcmmv; //Let's treat the "common-mode mean" on an APV-card by APV-card basis, but this would of course require us to define more database parameters. The main question is whether to add the "common-mode mean" by APV card to the individual strip pedestals
    TH2D *h
    
    fin->GetObject( histname.Format( "hpedestalU_%s%d", prefix, i ), hpedu );
    fin->GetObject( histname.Format( "hpedestalV_%s%d", prefix, i ), hpedv );
    fin->GetObject( histname.Format( "hrawADCpedsubU_%s%d", prefix, i ), hcmmu );
    fin->GetObject( histname.Format( "hrawADCpedsubV_%s%d", prefix, i ), hcmmv );
    
    //Do U strips first:

    int nstripsU = hpedu->GetNbinsX();
    int nstripsV = hpedv->GetNbinsX();

    int nAPVsU = nstripsU/128;
    int nAPVsV = nstripsV/128;

    
    
    
    TH1D *htemp;
 
    //analyze U strips:
    vector<double> pedmeanU(nstripsU), pedrmsU(nstripsU);

    bool firstpage=true;
    
    cout << "getting pedestal mean and rms, module " << i << " U strips" << endl;

    TClonesArray *hPedU_1Dhistos = new TClonesArray("TH1D", nstripsU);
    
    for( int istrip=1; istrip <= nstripsU; istrip++ ){
      c1->cd(1+(istrip-1)%64);
      
      htemp = hpedu->ProjectionY( histname.Format( "hPedU_m%d_strip%d", i, istrip-1 ), istrip, istrip );

      new( (*hPedU_1Dhistos)[istrip-1] ) TH1D( *htemp );

      ( (TH1D*) (*hPedU_1Dhistos)[istrip-1] )->SetTitle(histname.Format("Ped. dist. module %d, U strip %d", i, istrip-1) );
      
      double histmean = ( (TH1D*) (*hPedU_1Dhistos)[istrip-1] )->GetMean();
      double histrms = ( (TH1D*) (*hPedU_1Dhistos)[istrip-1] )->GetRMS();
      
      if( plotflag != 0 ) ( (TH1D*) (*hPedU_1Dhistos)[istrip-1] )->Draw();
      ( (TH1D*) (*hPedU_1Dhistos)[istrip-1] )->Fit("gaus","SQ","",histmean-2.5*histrms,histmean+2.5*histrms);
      if( plotflag != 0 ) ( (TH1D*) (*hPedU_1Dhistos)[istrip-1] )->GetXaxis()->SetRangeUser(histmean-5.0*histrms,histmean+5.0*histrms);
      double fitmean = ( (TF1*) ( (TH1D*) (*hPedU_1Dhistos)[istrip-1] )->GetListOfFunctions()->FindObject("gaus") )->GetParameter(1);
      double fitsigma = ( (TF1*) ( (TH1D*) (*hPedU_1Dhistos)[istrip-1] )->GetListOfFunctions()->FindObject("gaus") )->GetParameter(2);

      pedmeanU[istrip-1] = fitmean;
      pedrmsU[istrip-1] = fitsigma;

      hPedRMSU_distribution_vs_module->Fill( i, histrms );
      hPedSigmaU_distribution_vs_module->Fill( i, fitsigma );
      hPedMeanU_distribution_vs_module->Fill( i, histmean );
      hPedFitMeanU_distribution_vs_module->Fill( i, fitmean );

      //gPad->Modified();
      //c1->Update(); 
      //gSystem->ProcessEvents();

      
      
      if( plotflag != 0 && istrip%64 == 0 ) {
	TString pdffilename;
	pdffilename.Form("pedplots_module%d.pdf", i);
	if( firstpage ){
	  pdffilename += "(";
	  firstpage=false;
	}
	
	c1->Print(pdffilename.Data());
	c1->Update();
	gSystem->ProcessEvents();
	
	c1->Clear();
	
	c1->Divide(8,8,.001,.001);

	
      }

      
      
      //fout->cd();
      //htemp->Write();
      //htemp->Delete();

      htemp->Delete();
    }

    hPedU_1Dhistos->Delete();

    //output pedestal mean, U strips:
    header.Form("%s%d.pedu = ", prefix, i );

    outfile << header << endl;

    for( int istrip=0; istrip<nstripsU; istrip++ ){
      TString entry;
      entry.Form("  %15.5g ", pedmeanU[istrip] );
      outfile << entry;
      if( (istrip+1)%16 == 0 ) outfile << endl;
    }
    outfile << endl;
    
    //output pedestal rms, U strips: 
    header.Form("%s%d.rmsu = ", prefix, i );

    outfile << header << endl;

    for( int istrip=0; istrip<nstripsU; istrip++ ){
      TString entry;
      entry.Form("  %15.5g ", pedrmsU[istrip] );
      outfile << entry;
      if( (istrip+1)%16 == 0 ) outfile << endl;
    }

    outfile << endl;

    //V strips:

    TClonesArray *hPedV_1Dhistos = new TClonesArray("TH1D", nstripsV);

    vector<double> pedmeanV(nstripsV), pedrmsV(nstripsV);

    cout << "getting pedestal mean and rms, module " << i << " V strips" << endl;
    for( int istrip=1; istrip <= nstripsV; istrip++ ){
      
      c1->cd(1+(istrip-1)%64);
      
      htemp = hpedv->ProjectionY( histname.Format( "hPedV_m%d_strip%d", i, istrip-1 ), istrip, istrip );

      new( (*hPedV_1Dhistos)[istrip-1] ) TH1D( *htemp );
      
      ( (TH1D*) (*hPedV_1Dhistos)[istrip-1] )->SetTitle(histname.Format("Ped. dist. module %d, V strip %d", i, istrip-1) );
      
      double histmean = ( (TH1D*) (*hPedV_1Dhistos)[istrip-1] )->GetMean();
      double histrms = ( (TH1D*) (*hPedV_1Dhistos)[istrip-1] )->GetRMS();

      if( plotflag != 0 ) ( (TH1D*) (*hPedV_1Dhistos)[istrip-1] )->Draw();
      ( (TH1D*) (*hPedV_1Dhistos)[istrip-1] )->Fit("gaus","SQ","",histmean-2.5*histrms,histmean+2.5*histrms);
      if( plotflag != 0 ) ( (TH1D*) (*hPedV_1Dhistos)[istrip-1] )->GetXaxis()->SetRangeUser(histmean-5.0*histrms,histmean+5.0*histrms);
      
      double fitmean = ( (TF1*) ( (TH1D*) (*hPedV_1Dhistos)[istrip-1] )->GetListOfFunctions()->FindObject("gaus") )->GetParameter(1);
      double fitsigma = ( (TF1*) ( (TH1D*) (*hPedV_1Dhistos)[istrip-1] )->GetListOfFunctions()->FindObject("gaus") )->GetParameter(2);

      pedmeanV[istrip-1] = fitmean;
      pedrmsV[istrip-1] = fitsigma;

      hPedRMSV_distribution_vs_module->Fill( i, histrms );
      hPedSigmaV_distribution_vs_module->Fill( i, fitsigma );

      hPedMeanV_distribution_vs_module->Fill( i, histmean );
      hPedFitMeanV_distribution_vs_module->Fill( i, fitmean );

      
      
      if( plotflag != 0 && istrip%64 == 0 ) {
	TString pdffilename;
	pdffilename.Form("pedplots_module%d.pdf", i);
       
	if( istrip == nstripsV ){
	  pdffilename += ")";
	}
	
	c1->Print(pdffilename.Data());
	c1->Update();
	gSystem->ProcessEvents();
	
	c1->Clear();
	
	c1->Divide(8,8,.001,.001);

	
      }

      
      //htemp->Delete();

      
      
      //fout->cd();
      //htemp->Write();
      //htemp->Delete();

      htemp->Delete();
      
    }

    //delete all histos after the fact
    hPedV_1Dhistos->Delete();

    //output pedestal mean, U strips:
    header.Form("%s%d.pedv = ", prefix, i );

    outfile << header << endl;

    for( int istrip=0; istrip<nstripsV; istrip++ ){
      TString entry;
      entry.Form("  %15.5g ", pedmeanV[istrip] );
      outfile << entry;
      if( (istrip+1)%16 == 0 ) outfile << endl;
    }

    outfile << endl;
    //output pedestal rms, U strips: 
    header.Form("%s%d.rmsv = ", prefix, i );

    outfile << header << endl;

    for( int istrip=0; istrip<nstripsV; istrip++ ){
      TString entry;
      entry.Form("  %15.5g ", pedrmsV[istrip] );
      outfile << entry;
      if( (istrip+1)%16 == 0 ) outfile << endl;
    }

    outfile << endl;
    
  }


  hPedRMSU_distribution_vs_module->Write();
  hPedSigmaU_distribution_vs_module->Write();
  hPedRMSV_distribution_vs_module->Write();
  hPedSigmaV_distribution_vs_module->Write();

  hPedMeanU_distribution_vs_module->Write();
  hPedFitMeanU_distribution_vs_module->Write();

  hPedMeanV_distribution_vs_module->Write();
  hPedFitMeanV_distribution_vs_module->Write();
  //fout->Write(); don't write out everything
}
