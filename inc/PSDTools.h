// ROOT
#include "TVirtualFFT.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TF1.h"
#include "TGraph.h"
#include "TRandom3.h"

// Minimizer stuff
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

// Custom
#include "TK95.h"


// Standard
#include <fstream>
#include <algorithm>
#include <cmath>
#include <vector>

/*
  PSD Tools for analyzing light curves

  Taking the likelihood equations from
  https://ui.adsabs.harvard.edu/abs/2013MNRAS.433..907E/abstract
*/


#ifndef __PSDTOOLS_H_
#define __PSDTOOLS_H_

class PSDTools {
private:

  TK95 *LCGen;

  // For the light curve
  int fNpoints;
  double fdT;
  double *fTime;
  double *fFlux;
  double fMeanFlux;

  // FFT tools
  TVirtualFFT *fFFT;    // Using ROOT's FFT
  double* fRe;          // Real Component of FFT
  double* fIm;          // Complex Component of FFT
  double* fAmp;         // Amplitude of FFT
  double* fPhi;         // Phase of FFT
  double* fOmega;       // Frequency of FFT
  double* fPSD;         // PSD
  double fMeanTimeDiff;

  bool fEven;
  int fNlim;
  int fFitStatus;

  // PSD Models
  TGraph *fPDF;
  int fModelID;
  int fNParms;

  double BrokenPowerLaw(double *x, double* parms);

  // Random Number Generator
  TRandom3 *fRand;

  void ClearObjs();


public:

  PSDTools ();
  ~PSDTools ();

  double fOmegaMax;
  TF1 *fPSDModel;



  // PSD Model
  void SetModel(int modelID = 0);
  void SetModelParameters(int n, double* parms);
  void SetModelParameter(int i, double parm);

  void SetLightCurve(int npoints, double dT, double *time, double *flux);

  void CalculatePSD();
  // void CalculatePSDUnEqual();

  double Get2LogL(const double* parms);

  TGraph* GetPSD();
  double* FitPSD();

  int GetFitStatus();

  // void GetAmpPhi(double *iAmp, double *iPhi);
  void GetAmpPhi(std::vector <double> &iAmp, std::vector <double> &iPhi);


};
#endif
