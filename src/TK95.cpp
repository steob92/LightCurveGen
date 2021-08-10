#include "TK95.h"


using namespace std;


TK95::TK95 ()
{

  fPSDModel = 0;

  fRand = new TRandom3();
  fRand->SetSeed(0);

  fPI = TMath::Pi();
}

TK95::~TK95 ()
{
  delete fPSDModel;
  delete fFreq;
  delete fRe;
  delete fIm;
  delete fRand;
}


/*
  Set the Model for the PSD to generate from:
  0 - Power Law
  Put in more later

*/
void TK95::SetModel(int modelID)
{

  if (fPSDModel)
  {
    delete fPSDModel;
  }

  // Power Law Model
  if (modelID == 0)
  {
    fPSDModel = new TF1("fPSDModel", "TMath::Power(x, -[0])", 0, 1e5);
  }
  else
  {
    std::cout << "Model not implemented, defaulting to Power Law" << std::endl;
    SetModel(0);
  }
}



/*
  Set the Model for the PSD to generate from:
  Custom input model

*/
void TK95::SetModel(TF1* model)
{
  fPSDModel = (TF1*)model->Clone();
}


// Setting the PSD Parameters
void TK95::SetModelParameters(int n, double* parms)
{
  for (unsigned int i = 0; i < n; i++)
  {
    SetModelParameter(i, parms[i]);
  }
}

void TK95::SetModelParameter(int i, double parm)
{
  fPSDModel->SetParameter(i, parm);
}


/*
  Getting random light curve based on TK95 method
  n - number of points
  dt - time spacing
  fLC - flux points to be returned
  tLC - time points to be returned
*/
void TK95::GetRandomLightCurve(int n, double dt, double *tLC, double *fLC)
{

  double dOmega = (2 * fPI / dt / n);
  double omega = 0;

  fRe = new double[n];
  fIm = new double[n];

  bool bEven = n % 2 == 0;
  int nlim = bEven ? n/2 -1  : (n-1) / 2;

  for (int i = 1; i <= nlim ; i++)
  {
    omega = (i) / dt / n;
    // Gaussian random number to sample PSD
    fRe[i] = fRand->Gaus(0,1) * TMath::Sqrt(0.5 * fPSDModel->Eval(omega));
    fIm[i] = fRand->Gaus(0,1) * TMath::Sqrt(0.5 * fPSDModel->Eval(omega));

    // Symetry
    fRe[nlim + i] = fRe[i];
    fIm[nlim + i] = -fIm[i];

  }

  fRe[0] = 0;
  fIm[0] = 0;

  if (bEven)
  {
    fIm[nlim] = 0;
  }

  TVirtualFFT *fft = TVirtualFFT::FFT(1, &n, "C2R");
  fft->SetPointsComplex(fRe,fIm);
  fft->Transform();

  fft->GetPoints(fLC);
  double fmin = *std::min_element(fLC,fLC+n);
  double fmax = *std::max_element(fLC,fLC+n);

  // min/max normalization
  for (unsigned int i = 0; i < n; i++)
  {
    tLC[i] = i * dt;
    fLC[i] = (fLC[i] - fmin) /  (fmax - fmin);
  }

  delete []fRe;
  delete []fIm;
  delete fft;
}
