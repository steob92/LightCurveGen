/*
  TK 95 https://ui.adsabs.harvard.edu/abs/1995A%26A...300..707T/abstract

  Generating Light Curves base on an input PSD
  Simulated light curves will have a Gaussian PDF

  Rely on ROOT's packages for FFT functionality
*/


#include "TVirtualFFT.h"
#include "TRandom3.h"
#include "TF1.h"

#ifndef __TK95_H_
#define __TK95_H_

class TK95 {
private:
  /* data */
  // Assumed PSD Model
  TF1 *fPSDModel;

  // FFT Items
  std::vector <double> fFreq;
  std::vector <double> fRe;
  std::vector <double> fIm;

  double fPI;

  TRandom3 *fRand;

public:
  TK95 ();
  virtual ~TK95 ();

  // Set PSD Model
  void SetModel(int modelID = 0);
  void SetModel(TF1* model);

  // Setting the PSD Parameters
  void SetModelParameters(int n, double* parms);
  void SetModelParameter(int i, double parm);

  // Get the random Light curve
  void GetRandomLightCurve(int n, double dt,  std::vector <double> &tLC, std::vector <double> &fLC);
};
#endif
