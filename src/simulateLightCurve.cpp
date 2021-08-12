/*

    Check the PSD distribution from the EMP13 light curves
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

// Custom
#include "TK95.h"
#include "EMP13.h"
#include "PSDTools.h"



using namespace std;

void help_message()
{
    cout << "This is a helpful message" << endl;
}


int main(int argc, char *argv[])
{


    // inital default values
    // string fFilename = "/mnt/Storage/VERITAS/SourceAnalysis/OJ287/2017_Flare/Data/LCData/SwiftXRT_Soft_LC.dat";;
    string fFilename = "/Users/obriens/OJ287_Update/Data/LCData/SwiftXRT_Soft_LC.dat";
    string fOutfile = "test.root";
    bool bBetaForced = false;
    double fBeta = 0;
    int fNSim = 100;
    double fDelT = 1; // Time spacing between observations (days)

    string arg; 

    if (argc == 1)
    {
        help_message();
        return 0;
    }

    // read the CMD options
    for (int i = 1; i < argc; i++)
    {
        arg = argv[i];
        // help message
        if ((arg == "-h") || (arg == "--help"))
        {
            help_message();
            return 0;
        }
        // Overwrite default file name
        else if ((arg == "-i") || (arg == "--infile"))
        {
            fFilename = argv[i +1];
            i++;
        }
        // Force an assumed Beta for PSD
        else if ((arg == "-b") || (arg == "--beta"))
        {
            bBetaForced = true;
            fBeta = atof(argv[i+1]);
            i++;
        }
        // Overwrite file name
        else if ((arg == "-o") || (arg == "--outfile"))
        {
            fOutfile = argv[i+1];
            i++;
        }
        // Overwrite number of sims
        else if ((arg == "-n") || (arg == "--nsim"))
        {
            fNSim = atoi(argv[i+1]);
            i++;
        }
        // Set time bin
        else if ((arg == "-t") || (arg == "--time"))
        {
            fDelT = atof(argv[i+1]);
            i++;
        }
        else
        {
            help_message();
            return 0;
        }
    }


    // Read in the data file using a ROOT TTree
    TTree *tLC = new TTree();
    // Assuming the XRT format
    tLC->ReadFile(fFilename.c_str(), "MJD/D:Flux/D:FluxErrL/D:FluxErrH/D");

    double imjd = 0;
    double iflux = 0;
    double iflux_errl = 0;
    double iflux_erru = 0;

    tLC->SetBranchAddress("MJD", &imjd);
    tLC->SetBranchAddress("Flux", &iflux);
    tLC->SetBranchAddress("FluxErrL", &iflux_errl);
    tLC->SetBranchAddress("FluxErrH", &iflux_erru);


    int fNPoints = tLC->GetEntries();

    vector <double> mjd(fNPoints);
    vector <double> flux(fNPoints);
    vector <double> flux_err(fNPoints);

    // Read the data
    cout << "Reading Data" << endl;
    for (int i = 0; i < fNPoints; i ++)
    {
        tLC->GetEntry(i);
        mjd[i] = imjd;
        flux[i] = iflux;
        flux_err[i] = 0.5*(iflux_errl + iflux_erru); // take the mean error
    }

    /*
        Assume that the passed data is approximately equi-sampled and well behaved on short time scales
        Resample the data into equi-spaced bins.
        Using the TGraph::Eval method to interpolate between points
    */
    // Create TGraph of input data
    TGraph *gInput = new TGraph (fNPoints, &(mjd[0]), &(flux[0]));
    
    // Define the time of interest
    int fMJD_min = int(mjd[0]);
    int fMJD_max = int(mjd[fNPoints-1]) +1;
    int fNData = (fMJD_max - fMJD_min) / fDelT;

    // Get interpolated light curves
    vector <double> inter_mjd(fNData);
    vector <double> inter_flux(fNData);

    for (int i = 0 ; i < fNData; i ++ )
    {
        inter_mjd[i] = fMJD_min + i*fDelT;
        inter_flux[i] = gInput->Eval(inter_mjd[i]);
    }
    
    // Create EMP13 object
    EMP13 *lc_emp = new EMP13();
    lc_emp->SetLightCurve(fNData, fDelT, &(inter_mjd[0]), &(inter_flux[0]));
    lc_emp->SetModel(0);  // power law model

    // Create a PSD object
    PSDTools *psd = new PSDTools();
    psd->SetModel(0);

    // Get the best-fit PSD if required
    if (!bBetaForced)
    {
        psd->SetLightCurve( fNData, fDelT, inter_mjd, inter_flux);
        double *inParms =  psd->FitPSD();
        fBeta = inParms[1];
        cout << "Best Fit Beta: " << fBeta << endl;
    }

    // Set PSD index to be simulated
    lc_emp->SetModelParameter(0, fBeta);



    // For debugging/quality testing
    // Best fit beta
    TH1D *hBeta = new TH1D("hBeta", "hBeta", 100, fBeta-0.5, fBeta+0.5);
    // Number of itterations until convergence
    TH1D *hConv = new TH1D("hConv", "hConv", 100, 0, 1000);
    
    vector <double> iSimTime;
    vector <double> iSimFlux;


    // Open TFile for later storage
    TFile *fFile = new TFile(fOutfile.c_str(), "RECREATE");
    TTree *fTree = new TTree("SimulationTree", "SimulationTree");

    // Store resampled random light curves
    vector <double> iResampleMJD;
    vector <double> iResampleFlux;


    TGraph *gSimLC = 0;
    fTree->Branch("fNPoints", &fNPoints, "fNpoints/I");
    fTree->Branch("Beta", &fBeta, "fBeta/d");
    fTree->Branch("mjd", &iResampleMJD);
    fTree->Branch("flux", &iResampleFlux);

    

    int nitter = 0;
    TGraph *gResampler = 0;
    // Main simulation loop
    while (nitter < fNSim)
    {
        iResampleMJD.assign(fNPoints, 0);
        iResampleFlux.assign(fNPoints, 0);
        // if ((10*nitter / fNSim ) % 10 == 0)
        // {
        //     cout << nitter << endl;
        // }

        // if (gSimLC) {delete gSimLC;}ls
        iSimTime.clear();
        iSimFlux.clear();
        // Simulate light curve
        int conv = lc_emp->GetRandomLightCurveEMP13(iSimTime, iSimFlux, fDelT, 1, fNData);
        
        psd->SetLightCurve( fNData, fDelT, iSimTime, iSimFlux);
        double *fitParms =  psd->FitPSD();
        // Only take converged fits
        if (psd->GetFitStatus() != 0)
        {continue;}

        // Book keeping
        hBeta->Fill(fitParms[1]);
        hConv->Fill(conv);
        nitter++;

        // cout << "Resampling" << endl;
        if (gResampler) {delete gResampler;}
        gResampler = new TGraph(fNData, &(inter_mjd[0]), &(iSimFlux[0]));
        for (int i = 0; i < fNPoints; i ++)
        {
            iResampleMJD[i] = mjd[i];
            iResampleFlux[i] = gResampler->Eval(mjd[i]);
        }

        gSimLC = new TGraph(fNData, &(iResampleMJD[0]), &(iResampleFlux[0]) );
        fTree->Fill();
    }


    fTree->Write();
    TCanvas *cBeta = new TCanvas();
    hBeta->Draw();
    cBeta->Print("Beta.png");
    cBeta->Write();
    TCanvas *cConv = new TCanvas();
    hConv->Draw();
    cConv->SetLogx();
    cConv->Print("Conv.png");
    cConv->Write();


    fFile->Close();
    return 1;
}