#include "Riostream.h"
#include <TFile.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TGraphErrors.h>
#include <TTree.h>
#include <TROOT.h>
#include <TLegend.h>
#include <string>
#include <array>
#include <vector>

void CalibratePedestal(const char* filename = "bbcal_179.root"){
	//const Int_t entrySize = rows * cols;
	//Open the file.
	TFile *events = TFile::Open(filename);

	//Get the Tree.
	TTree* tree = 0;
	events->GetObject("T",tree);

	//Set the variable to hold the values.
	Double_t data[189];
	tree->SetBranchAddress("bb.sh.a_p",&data);

	//Get the number of events.
	Int_t nEvents =(Int_t) tree->GetEntries();

	//c1->SetLogy();
	vector<Double_t> rmsValues;
	vector<Double_t> meanValues;

	for(Int_t cell = 0; cell < 189; cell++){
		//Zoom out
		auto outDisplay = new TH1D("outDisplay","Cell",600,-1000,5000);

		//Put in the values here.
		vector<Double_t> values;

		//Read in data
		for(int event = 0; event < nEvents; event++){
			//Read in that data.
			tree->GetEntry(event);
			//Fill the bins with content.
			values.push_back(data[cell]);
			event++;
		}

		//Put data on the zoomed-out histogram.
		for(int event = 0; event < values.size(); event++){
			outDisplay->Fill(values[event]);
		}
		
		//Find max so we know where to zoom in.
		auto zoomBin = outDisplay->GetMaximumBin();		
		auto zoomMaxLow = outDisplay->GetBinLowEdge(zoomBin-1);
		auto zoomMaxHigh = outDisplay->GetBinLowEdge(zoomBin+1) + outDisplay->GetBinWidth(zoomBin);

		//Zoom in.
		auto inDisplay = new TH1D("inDisplay","Cell",30,zoomMaxLow,zoomMaxHigh);
		
		delete gROOT->FindObject("outDisplay");
		
		//Fill in this histogram.
		for(int event = 0; event < values.size(); event++){
			inDisplay->Fill(values.at(event));
		}
			
		//Get the sum of all values within the histogram		
		Double_t sum = 0;
		Double_t squareSum = 0;
		Double_t entries = 0;
		for(int i = 1; i < 30; i++){
			entries+=inDisplay->GetBinContent(i);
			Double_t area = inDisplay->GetBinCenter(i)*inDisplay->GetBinContent(i);	
			sum += area;
		}
		
		//Calculate the mean of the histogram
		meanValues.push_back(sum/entries);

		//Calculate the RMS of the histogram
		for(int i = 1; i < 30; i++){
			squareSum+=pow(inDisplay->GetBinCenter(i)-meanValues[cell],2)*inDisplay->GetBinContent(i);
		}
		rmsValues.push_back(sqrt(squareSum/entries));
		
		//Memory Clean up
		delete gROOT->FindObject("inDisplay");
	}

	//Do the same for the pre-shower.
	//Set the variable to hold the values.
	Double_t data1[54];
	tree->SetBranchAddress("bb.ps.a_p",&data1);
	vector<Double_t> rmsValues1;
	vector<Double_t> meanValues1;

	for(Int_t cell = 0; cell < 54; cell++){
		//Zoom out
		auto outDisplay = new TH1D("outDisplay","Cell",600,-1000,5000);

		//Put in the values here.
		vector<Double_t> values;

		//Read in data
		for(int event = 0; event < nEvents; event++){
			//Read in that data.
			tree->GetEntry(event);
			//Fill the bins with content.
			values.push_back(data1[cell]);
			event++;
		}

		//Put data on the zoomed-out histogram.
		for(int event = 0; event < values.size(); event++){
			outDisplay->Fill(values[event]);
		}
		
		//Find max so we know where to zoom in.
		auto zoomBin = outDisplay->GetMaximumBin();		
		auto zoomMaxLow = outDisplay->GetBinLowEdge(zoomBin-1);
		auto zoomMaxHigh = outDisplay->GetBinLowEdge(zoomBin+1) + outDisplay->GetBinWidth(zoomBin);

		//Zoom in.
		auto inDisplay = new TH1D("inDisplay","Cell",30,zoomMaxLow,zoomMaxHigh);
		
		delete gROOT->FindObject("outDisplay");
		
		//Fill in this histogram.
		for(int event = 0; event < values.size(); event++){
			inDisplay->Fill(values.at(event));
		}
			
		//Get the sum of all values within the histogram		
		Double_t sum = 0;
		Double_t squareSum = 0;
		Double_t entries = 0;
		for(int i = 1; i < 30; i++){
			entries+=inDisplay->GetBinContent(i);
			Double_t area = inDisplay->GetBinCenter(i)*inDisplay->GetBinContent(i);	
			sum += area;
		}
		
		//Calculate the mean of the histogram
		meanValues1.push_back(sum/entries);

		//Calculate the RMS of the histogram
		for(int i = 1; i < 30; i++){
			squareSum+=pow(inDisplay->GetBinCenter(i)-meanValues1[cell],2)*inDisplay->GetBinContent(i);
		}
		rmsValues1.push_back(sqrt(squareSum/entries));
		
		//Memory Clean up
		delete gROOT->FindObject("inDisplay");
	}	

	//Write this to a file.
	TFile outputFile ("pedestalcalibrated.root","RECREATE");
	outputFile.WriteObjectAny(&meanValues, "std::vector<Double_t>", "Mean");
	outputFile.WriteObjectAny(&rmsValues, "std::vector<Double_t>", "RMS");
	outputFile.WriteObjectAny(&meanValues1, "std::vector<Double_t>", "MeanPS");
	outputFile.WriteObjectAny(&rmsValues1, "std::vector<Double_t>", "RMSPS");
	outputFile.Close();
}
