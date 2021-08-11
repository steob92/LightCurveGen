// Standard
#include <fstream>
#include <algorithm>
#include <numeric>
#include <vector>

// ROOT
#include "TVirtualFFT.h"
#include "TRandom3.h"
#include "TF1.h"
#include "TGraph.h"

// Custom
#include "TK95.h"
#include "PSDTools.h"

#ifndef __EMP13_H_
#define __EMP13_H_

class EMP13: public TK95 {
private:
  /* data */

  // FFT Items
  double *fFreq;
  double *fRe;
  double *fIm;
  TVirtualFFT *fFFT;

  // Periodogram things
  int fNlim;

  bool bEven;
  double fMeanFlux;
  double fdT;

  TGraph *fPDF;
  double *fAmp;
  double *fPhi;
  double *fOmega;
  double *fPSD;

  double fPI;

  TRandom3 *fRand;


  // Light curve book keeping
  int fNpoints;
  double *fTime;
  double *fFlux;
  PSDTools *fPSDInput;

public:
  EMP13 ();
  ~EMP13 ();


  void SetLightCurve(int npoints, double dT, double *time, double *flux);
  void CalculatePSD();
  void SetPDF();
  TGraph *GetPDF(){return (TGraph*)fPDF->Clone();}
  void GetRandomFlux(double *iflux, int npoints = -1);
  int GetRandomLightCurveEMP13(std::vector <double>&iTime, std::vector <double>&iFlux, double dt, int nRed, int npoints = -1);

  int GetNPoints(){return fNpoints;}

};


// Shamelessly grabbed from Stackoverflow
class sort_indices
{
   private:
     double* mparr;
   public:
     sort_indices(double* parr) : mparr(parr) {}
     bool operator()(int i, int j) const { return mparr[i]>mparr[j]; }
};

template <typename T>
std::vector<size_t> sort_indexes(const std::vector<T> &v) {

  // initialize original index locations
  std::vector<size_t> idx(v.size());
  iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing values in v
  // using std::stable_sort instead of std::sort
  // to avoid unnecessary index re-orderings
  // when v contains elements of equal values
  stable_sort(idx.begin(), idx.end(),
       [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

  return idx;
}
#endif
