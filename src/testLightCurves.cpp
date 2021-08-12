/*
    Test Light Curve Class
*/
// Standard
#include <iostream>
#include <vector>
#include <string>
#include <stdlib.h>

// ROOT
#include "TCanvas.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TTree.h"
#include "TFile.h"
#include "TLine.h"
#include "TGraphErrors.h"

// Custom
#include "TK95.h"
#include "EMP13.h"
#include "PSDTools.h"
#include "LightCurve.h"
#include "DCF.h"

using namespace std;

int main()
{
    string infile = "/mnt/Storage/VERITAS/SourceAnalysis/OJ287/2017_Flare/Data/LCData/SwiftXRT_Soft_LC.dat";
    /* 
    This code is specific to the my own system...
    ToDo: Replace with a publically available dataset
    */
    // string infile = "/Users/obriens/OJ287_Update/Data/LCData/SwiftXRT_Soft_LC.dat";

    // Read in data from TTree
    TTree *tLC = new TTree();
    tLC->ReadFile(infile.c_str(), "MJD/D:Flux/D:FluxErrL/D:FluxErrH/D");

    double imjd = 0;
    double iflux = 0;
    double iflux_errl = 0;
    double iflux_erru = 0;

    tLC->SetBranchAddress("MJD", &imjd);
    tLC->SetBranchAddress("Flux", &iflux);
    tLC->SetBranchAddress("FluxErrL", &iflux_errl);
    tLC->SetBranchAddress("FluxErrH", &iflux_erru);


    vector <double> mjd(tLC->GetEntries());
    vector <double> flux(tLC->GetEntries());
    vector <double> flux_err(tLC->GetEntries());

    // Reading data
    for (int i = 0; i < tLC->GetEntries(); i ++)
    {
        tLC->GetEntry(i);
        mjd[i] = imjd;
        flux[i] = iflux;
        flux_err[i] = 0.5*(iflux_errl + iflux_erru);
    }

    LightCurve *LC1 = new LightCurve( mjd, flux, flux_err);
    LC1->CalculateProperties();


    LC1->PrintDetails();


    // Look at the DCF
    DCF *myDCF = new DCF(LC1, LC1);
    myDCF->SetTimeDetails(1, -100, 100);

    TGraphErrors *gDCF = myDCF->CalculateDCF();
    TCanvas *c1 = new TCanvas();
    gDCF->SetLineStyle(1);
    gDCF->SetMarkerStyle(8);
    gDCF->SetTitle("DCF;Tau[#tau]; DCF Value");
    gDCF->Draw("APL");

    c1->SetGrid();
    c1->Print("DCF.png");
    

    return 0;
}