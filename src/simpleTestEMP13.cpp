/*
    
*/
// Standard
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

// ROOT
#include "TCanvas.h"
#include "TGraph.h"
#include "TTree.h"

// Custom
#include "TK95.h"
#include "EMP13.h"
#include "PSDTools.h"

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


  double *mjd = new double[tLC->GetEntries()];
  double *flux = new double[tLC->GetEntries()];
  double *flux_err = new double[tLC->GetEntries()];

  // Reading data
  for (int i = 0; i < tLC->GetEntries(); i ++)
  {
    tLC->GetEntry(i);
    mjd[i] = imjd;
    flux[i] = iflux;
    flux_err[i] = 0.5*(iflux_errl + iflux_erru);
  }


  // To a TGraph
  TGraph *gFlux = new TGraph(tLC->GetEntries(), &(mjd[0]), &(flux[0]));
  gFlux->SetTitle("Input LC;Time[MJD];Flux");
  // Get number of days to sim
  int nPoints = int(mjd[tLC->GetEntries()-1]) + 1 - int(mjd[0]);

  double *inter_time = new double[nPoints];
  double *inter_flux = new double[nPoints];
  for (int i = 0; i < nPoints; i++)
  {
    inter_time[i] = int(mjd[0]) + i;
    // inter_flux[i] = gFlux->Eval(inter_time[i], 0 , "S");
    inter_flux[i] = gFlux->Eval(inter_time[i]);
  }

  TGraph *gFlux_inter = new TGraph(nPoints, &(inter_time[0]), &(inter_flux[0]));
  gFlux_inter->SetMarkerStyle(8);
  gFlux_inter->SetMarkerColor(2);
  gFlux_inter->SetLineColor(2);





  // Simulate a light curve
  EMP13 *lc_emp = new EMP13();
  lc_emp->SetLightCurve(nPoints, 1, &(inter_time[0]), &(inter_flux[0]));
  lc_emp->SetModel(0);
  lc_emp->SetModelParameter(0, 1.8);

  vector <double> iSimTime;
  vector <double> iSimFlux;

  lc_emp->GetRandomLightCurveEMP13(iSimTime, iSimFlux, 1, 1, nPoints);
  // // double *iSimTime = new double[100*nPoints];
  // // double *iSimFlux = new double[100*nPoints];
  // // lc_emp->GetRandomLightCurve2(iSimFlux, iSimTime, 100, 100*nPoints);
  // cout << "Done" << endl;
  // // for (int i = 0 ; i < nPoints; i++)
  // // {
  // //   cout << i << " " << inter_time[i] << " " << inter_flux[i] << " " << iSimFlux[i] << endl;
  // // }
  PSDTools *psd = new PSDTools();

  psd->SetModel(0);
  psd->SetLightCurve( nPoints, 1.0, inter_time, inter_flux);
  double *inParms =  psd->FitPSD();
  psd->SetLightCurve( nPoints, 1.0, inter_time, &(iSimFlux[0]));
  double *outParms =  psd->FitPSD();

  TGraph *gSim = new TGraph(nPoints, &(inter_time[0]), &(iSimFlux[0]));
  gSim->SetLineColor(3);

  TCanvas *c1 = new TCanvas();
  gFlux->Draw("APL");
  gFlux_inter->Draw("SAME,L");
  gSim->Draw("SAME,L");

  c1->Print("InterLC.png");



  return 0;
}
