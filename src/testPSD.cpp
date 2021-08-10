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
#include "EMP13.h"
#include "PSDTools.h"


using namespace std;

int main()
{
    
  // Creating new TK95 light curve as a starting point


  int nPoints = 1000; // number of light curve points
  double delT = 1.;   // spacing of the light curve 
  double beta = 2.0;  // PSD index to be simed

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
  double *scaled_flux = new double[nPoints];

  // Get TK95 and rescale
  tk95->GetRandomLightCurve(nPoints, delT , time, flux);

  double i_min = 9999;
  double i_max = -9999;

  for (int i = 0; i < nPoints; i ++)
    {
      if (flux[i] < i_min){i_min = flux[i];}
      if (flux[i] > i_max){i_max = flux[i];}
    }

  for (int i = 0; i < nPoints; i ++)
    {
      scaled_flux[i] = 100 * (i_max - flux[i]) / (i_max - i_min);
      
    }
  

  EMP13 *lc_emp = new EMP13();
  lc_emp->SetLightCurve(nPoints, 1, &(time[0]), &(scaled_flux[0]));
  lc_emp->SetModel(0);
  lc_emp->SetModelParameter(0, 2.0);

  
  // Fit Parameters
  double *fitParms = 0;
  
  TH1D *hBeta = new TH1D("hBeta", "hBeta;Simulated Index;", 100, 1.5,2.5);
  TH1D *hConv = new TH1D("hConv", "hConv;Number of itterations;", 100, 0,100);
  
  int nsimmed = 0;
  
 
  
  // While used here so we can require accurate PSD estimate
  int nitter = 0;
  while (nsimmed < 1000)
  {

    vector <double> iSimTime;
    vector <double> iSimFlux;

    nitter = lc_emp->GetRandomLightCurveEMP13(iSimTime, iSimFlux, delT, 1, nPoints);


    psd->SetLightCurve( nPoints, delT, time, &(iSimFlux[0]));
			
    fitParms = psd->FitPSD();

    // Only take converged fits
    if (psd->GetFitStatus() != 0)
      {continue;}

    hBeta->Fill(fitParms[1]);
    hConv->Fill(nitter);
    nsimmed++;
  }

  TCanvas *c1 = new TCanvas();
  hBeta->Draw();
  c1->Print("PSD_index.png");
  hConv->Draw();
  c1->Print("PSD_indexConv.png");
  return 0;
}
