/*


*/
// Standard
#include <iostream>
#include <vector>

// ROOT
#include "TCanvas.h"
#include "TGraph.h"

// Custom
#include "TK95.h"
#include "PSDTools.h"

using namespace std;

int main()
{

  // Gen and analysis objects
  TK95 *tk95 = new TK95();
  PSDTools *psd = new PSDTools();

  // Power Law Model
  psd->SetModel(0);
  tk95->SetModel(0);
  tk95->SetModelParameter(0, 1);


  // indicies to sim
  int nSims = 5;
  int nPoints = 1000;
  double *indx = new double[nSims];

  // Vector to hold graphs
  vector <TGraph*> lc_graphs(nSims);
  vector <TGraph*> psd_graphs(nSims);
  vector <TF1*> psd_fits(nSims);


  // double arrays to hold the light curves
  double *time = new double[nPoints];
  double *flux = new double[nPoints];

  // Fit Parameters
  double *fitParms = 0;


  for (int i = 0; i < nSims; i ++)
  {
    // Itterate in 0.5
    indx[i] = 1.0 + 0.5*i;
    tk95->SetModelParameter(0, indx[i]);

    tk95->GetRandomLightCurve(nPoints, 1.0 , time, flux);
    lc_graphs[i] = new TGraph(nPoints, &(time[0]), &(flux[0]));
    lc_graphs[i]->SetLineColor(i+1);

    psd->SetLightCurve( nPoints, 1.0, time, flux);
    psd_graphs[i] = psd->GetPSD();
    psd_graphs[i]->SetLineColor(i+1);

    fitParms = psd->FitPSD();

    psd_fits[i] = (TF1*)psd->fPSDModel->Clone();
    psd_fits[i]->SetLineColor(i+1);
    cout << "Simulated: " << indx[i] << " " << fitParms[1] << endl;
  }


  // Plot Light Curves
  TCanvas *cLC = new TCanvas();
  lc_graphs[0]->SetTitle("Light Curve;MJD;Flux");

  lc_graphs[0]->Draw("AL");

  for (int i = 1 ; i < nSims; i++)
  {
    lc_graphs[i]->Draw("SAME,L");
  }
  cLC->SetGrid();
  cLC->Print("LightCurves.png");

  // PSD
  TCanvas *cPSD = new TCanvas();
  psd_graphs[0]->SetTitle("PSD; Freq ; PSD(Freq)");
  psd_graphs[0]->Draw("AL");
  psd_fits[0]->Draw("SAME");

  for (int i = 1 ; i < nSims; i++)
  {
    psd_graphs[i]->Draw("SAME,L");
    psd_fits[i]->Draw("SAME,L");
  }

  cPSD->SetLogy();
  cPSD->SetLogx();
  cPSD->SetGrid();
  cPSD->Print("PSD.png");

  return 0;

}
