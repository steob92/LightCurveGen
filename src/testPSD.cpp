/*

    Check the PSD distribution from the EMP13 light curves
*/
// Standard
#include <iostream>
#include <vector>

// ROOT
#include "TCanvas.h"
#include "TGraph.h"
#include "TH1D.h"

// Custom
#include "TK95.h"
#include "PSDTools.h"


using namespace std;

int main()
{
    
  // Creating new TK95 light curve as a starting point


  int nPoints = 1000; \\ number of light curve points
  double delT = 1.;   \\ spacing of the light curve 
  double beta = 2.0;  \\ PSD index to be simed

  // Gen and analysis objects
  TK95 *tk95 = new TK95();
  PSDTools *psd = new PSDTools();

  // Power Law Model
  psd->SetModel(0);
  tk95->SetModel(0);
  tk95->SetModelParameter(0, beta);
  
  
  // double arrays to hold the light curves
  double *time = new double[nPoints];
  double *flux = new double[nPoints];

  // Fit Parameters
  double *fitParms = 0;
  
  TH1D *hBeta = new TH2D("hBeta", "hBeta;Simulated Index;", 100, 0,5);
  
  int nsimmed = 0;
  
 
  
  // While used here so we can require accurate PSD estimate
  while (nsimmed < 1000)
  {
    
    
    tk95->GetRandomLightCurve(nPoints, delT , time, flux);
    lc_graphs[i] = new TGraph(nPoints, &(time[0]), &(flux[0]));
    lc_graphs[i]->SetLineColor(i+1);

    psd->SetLightCurve( nPoints, delT, time, flux);
    fitParms = psd->FitPSD();
    
    cout << psd->fFitStatus << endl;
    
    nsimmed++;
  }

  return 0;
}
