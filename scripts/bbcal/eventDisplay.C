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
#include <TCanvas.h>
#include <TStyle.h>

//Display an array, set a minimum value to display.
void Draw(Double_t data[], vector<Double_t>* means, Double_t min, const char* title, Int_t rows = 27, Int_t cols = 7){

	//Create a 2D histogram to store the data.
	auto display = new TH2I("display",title,cols,0,cols,rows,0,rows);



	//Fill the bins with data.
	for(int i = 0; i < cols; i++){
		for(int j = 0; j < rows; j++){

			//Only display if value(channel) > min.
			if(data[(j*cols)+i] - means->at((j*cols)+i) > min){

				//Fill the bins
				display->SetBinContent(i+1,j+1,(Int_t)(data[i+(j*cols)] - means->at(i+(j*cols))));
			} else {
				//Else display 0.
				display->SetBinContent(i+1,j+1,0);
			}
		}
	}
	display->GetYaxis()->SetNdivisions(rows);
	display->GetXaxis()->SetNdivisions(cols);

	//Have to completely replace axes to position axis labels where they are	
	display->GetYaxis()->SetLabelSize(0);
	display->GetXaxis()->SetLabelSize(0);	

	//Let's make a function to set the axis values
	TF1 *xfunc = new TF1("xfunc","x",1,cols+1);
	TF1 *yfunc = new TF1("yfunc","x",1,rows+1);

	//Set colour, box text size and remove the stats box.
	display->SetFillColor(kRed-9);
	display->SetStats(0);
	display->SetMarkerSize(1.2+((7-cols)/1.7));

	//Have to make the old axis labels invisible (unable to edit some axis properties directly).
	display->GetYaxis()->SetLabelSize(0);
	display->GetXaxis()->SetLabelSize(0);
	display->SetTitleSize(1/cols);

	//Display the event.
	display->Draw("box,text");	

	//Create completely new axes to get label in the middle of the divisions)
	TGaxis *x = new TGaxis(0,0,cols,0,"xfunc",cols+1,"M");
	x->SetLabelSize(0.25/cols);
	x->SetLabelOffset(0.015*(cols-7));
	x->Draw();

	TGaxis *y = new TGaxis(0,0,0,rows,"yfunc",rows+1,"M");
	y->SetLabelSize(0.25/cols);
	y->Draw();

	
	//Vertical lines.
	for (int i = 1; i < cols; i++){
		TLine *line = new TLine(i,0,i,rows);
		line->SetLineStyle(kDotted);
		line->Draw();	
	}

	//Horizontal lines.
	//Vertical lines.
	for (int i = 1; i < rows; i++){
		TLine *line = new TLine(0,i,cols,i);
		line->SetLineStyle(kDotted);
		line->Draw();	
	}

	//Memory clean up.
	//delete gROOT->FindObject("display");
}

void eventDisplay(){

	//Get the mean and RMS values.
	vector<Double_t>* rmsValues;
	vector<Double_t>* meanValues;
	vector<Double_t>* rmsValuesPS;
	vector<Double_t>* meanValuesPS;
	TFile *calibration = TFile::Open("pedestalcalibrated.root");
	meanValues = (vector<Double_t>*)calibration->Get("Mean");
	rmsValues = (vector<Double_t>*)calibration->Get("RMS");
	meanValuesPS = (vector<Double_t>*)calibration->Get("MeanPS");
	rmsValuesPS = (vector<Double_t>*)calibration->Get("RMSPS");
	calibration->Close();
	
	//Create a Canvas
	auto c1 = new TCanvas("c1","Event Display",500,1200);
	TPad *shower = new TPad("shower","Shower",0.01,0.01,0.7,0.99);
	shower->Draw();
	TPad *preshower = new TPad("preshower","Pre-Shower",.8,.01,.99,.99);
	preshower->SetLeftMargin(0.15);
	preshower->Draw();
	
	//Open the file.
	TFile *events = TFile::Open("bbcal_179.root");

	//Get the Tree.
	TTree* tree = 0;
	events->GetObject("T",tree);

	//Set the variable to hold the values.
	Double_t data[189];
	Double_t dataPS[54];
	tree->SetBranchAddress("bb.sh.a_p",&data);
	tree->SetBranchAddress("bb.ps.a_p",&dataPS);

	//Get the number of events.
	Int_t nEvents =(Int_t) tree->GetEntries();
	
	Int_t event = 0;
	Int_t cell = 0;

	for(Int_t event = 0; event < nEvents; event++){

		//Clear the data array
//		for(int i = 0; i < 189;i++){
//			data[i] = 0.0;
//			if(i < 54){
//				dataPS[i] = 0;			
//			}
//		}
		//Read in that data.
		tree->GetEntry(event);		

		//Get the total calibrated ADC value
		Double_t sum = 0;
		for(int i = 0; i < 189; i++){
			sum += data[i]-meanValues->at(i);	
		}

		//Don't display events with a total ADC value < 50
		if(sum < 50){
			continue;
		}

		//Create the histogram to draw this event.
		std::string title = "Event ";
		title += std::to_string(event);

		//Display the event
		//c1->cd(1);
		shower->cd();
		Draw(data,meanValues, 10, title.c_str(),27,7);
		//c1->cd(2);
		preshower->cd();
		Draw(dataPS,meanValuesPS, 10, "" ,27,2);
		gPad -> WaitPrimitive();

		//Memory clean up
		delete gROOT->FindObject("display");
	}
}


