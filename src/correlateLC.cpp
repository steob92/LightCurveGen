/*

    Correlate Light Curves against each other
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
#include "TH2D.h"
#include "TTree.h"
#include "TFile.h"
#include "TLine.h"

// Custom
#include "TK95.h"
#include "EMP13.h"
#include "PSDTools.h"
#include "LightCurve.h"
#include "DCF.h"


using namespace std;

int main()
{

    vector <string> fDataFiles(3);
    // fDataFiles[0] = "/mnt/Storage/VERITAS/SourceAnalysis/OJ287/2017_Flare/Data/LCData/SwiftXRT_Soft_LC.dat";
    // fDataFiles[1] = "/mnt/Storage/VERITAS/SourceAnalysis/OJ287/2017_Flare/Data/LCData/SwiftXRT_Mod_LC.dat";
    // fDataFiles[2] = "/mnt/Storage/VERITAS/SourceAnalysis/OJ287/2017_Flare/Data/LCData/SwiftXRT_Hard_LC.dat";

    fDataFiles[0] = "/Users/obriens/OJ287_Update/Data/LCData/SwiftXRT_Soft_LC.dat";
    fDataFiles[1] = "/Users/obriens/OJ287_Update/Data/LCData/SwiftXRT_Mod_LC.dat";
    fDataFiles[2] = "/Users/obriens/OJ287_Update/Data/LCData/SwiftXRT_Hard_LC.dat";

    vector <string> fSimFiles(3);
    fSimFiles[0] = "XRT_Soft.root";
    fSimFiles[1] = "XRT_Mod.root";
    fSimFiles[2] = "XRT_Hard.root";



    vector <string> fLabels(3);
    fLabels[0] = "Soft";
    fLabels[1] = "Mod";
    fLabels[2] = "Hard";

    vector <vector <double>> fMJD;
    vector <vector <double>> fFlux;
    vector <vector <double>> fFluxErr;

    for (int i = 0 ; i < fLabels.size(); i ++)
    {
        cout << fLabels[i] << endl;
        TTree *tLC = new TTree();
        tLC->ReadFile(fDataFiles[i].c_str(), "MJD/D:Flux/D:FluxErrL/D:FluxErrH/D");

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
        for (int j = 0; j < tLC->GetEntries(); j ++)
        {
            tLC->GetEntry(j);
            mjd[j] = imjd;
            flux[j] = iflux;
            flux_err[j] =  0.5*(iflux_errl + iflux_erru) ;
        }
        fMJD.push_back(mjd);
        fFlux.push_back(flux);
        fFluxErr.push_back(flux_err);

    }

    // for (int i = 0; i < )
    TH2D* hSimResults = 0;

    LightCurve *simLC1 = 0;
    LightCurve *simLC2 = 0;

    for (int i = 0; i < fLabels.size(); i ++)
    {
        TCanvas *c1 = new TCanvas();
        c1->Divide(2,2);



        LightCurve *LC1 = new LightCurve( fMJD[i], fFlux[i], fFluxErr[i]);
        TFile *fSim1 = new TFile(fSimFiles[i].c_str());
        int sim_fNpoints1;
        // double sim_mjd1[10000];
        // double sim_flux1[10000];
        // TGraph *gSimDataGraph1 = 0;
        // TGraph *gSimDataGraph2 = 0;

        vector <double> *simMJD1 = 0;
        vector <double> *simMJD2 = 0;
        vector <double> *simFlux1 = 0;
        vector <double> *simFlux2 = 0;




        
        TTree *SimData1 = (TTree*)fSim1->Get("SimulationTree");
        TTree *SimData2 = 0;

        SimData1->SetBranchAddress("fNPoints", &sim_fNpoints1);
        // SimData1->SetBranchAddress("gSimLC", &gSimDataGraph1);
        SimData1->SetBranchAddress("mjd", &simMJD1);
        SimData1->SetBranchAddress("flux", &simFlux1);


        for (int j = 0; j < fLabels.size(); j++)
        {

            LightCurve *LC2 = new LightCurve( fMJD[j], fFlux[j], fFluxErr[j]);

            DCF myDCF(LC1, LC2);
            myDCF.SetTimeDetails(1, -50, 50);

            TGraphErrors *gDCF = myDCF.CalculateDCF(false);

            hSimResults = new TH2D("hSimResults", "hSimResults", gDCF->GetN(), -50-0.5,50 + 0.5, 200, -2, 2 );

            TFile *fSim2 = new TFile(fSimFiles[j].c_str());


            int sim_fNpoints2;
            SimData2 = (TTree*)fSim2->Get("SimulationTree");


            SimData2->SetBranchAddress("fNPoints", &sim_fNpoints2);
            SimData2->SetBranchAddress("mjd", &simMJD2);
            SimData2->SetBranchAddress("flux", &simFlux2);

            for (int isim = 0; isim < 10; isim++)
            {
                // cout << isim << endl;
                SimData1->GetEntry(isim);

                // cout << "Pointing" << endl;
                vector <double> &vmjd1 = *simMJD1 ;
                vector <double> &vflux1 = *simFlux1 ;
                // cout << "Creating" << endl;

                // cout << vmjd1.size() << endl;
                simLC1 = new LightCurve( vmjd1, vflux1);
                // SimData1->GetEntry(isim);
                // cout << "Looping" << endl;

                for (int jsim = 0; jsim < 10; jsim++)
                {
                    SimData2->GetEntry(jsim);
                    vector <double> &vmjd2 = *simMJD2 ;
                    vector <double> &vflux2 = *simFlux2 ;

                    // vmjd2 = *simMJD2;
                    // vflux2 = *simFlux2;

                    simLC2 = new LightCurve( vmjd2, vflux2);

                    DCF myDCFSim(simLC1, simLC2);
                    myDCFSim.SetTimeDetails(1, -50, 50);

                    TGraphErrors *gDCFSim = myDCF.CalculateDCF(false);

                    for (int k = 0; k < gDCFSim->GetN(); k++)
                    {
                        hSimResults->Fill(gDCFSim->GetX()[k], gDCFSim->GetY()[k]);
                    }

                    // delete myDCFSim;
                    delete simLC2;
                    // delete gDCFSim;
                    

                }
                delete simLC1;
                // simLC1 = 0;
                
            }

            c1->cd(j+1);
            gDCF->SetLineStyle(1);
            gDCF->SetMarkerStyle(8);
            gDCF->SetTitle("DCF; Tau[#tau]; DCF Value");
            gDCF->Draw("APL");
            // hSimResults->Draw("COLZ");
            // gDCF->Draw("SAME,PL");
            
            gPad->SetGrid();

            delete hSimResults;
            delete SimData2;
            fSim2->Close();

        }
        fSim1->Close();
    
        c1->Print((fLabels[i] + ".png").c_str());

    }
    

    return 0;

}
