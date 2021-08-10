#include "PSDTools.h"


PSDTools::PSDTools ()
{
  LCGen = new TK95();

  //
  fPSDModel = 0;
  fPDF = 0;
  fFitStatus = 0;

  // Light Curve Bits
  fNpoints = 0;
  fTime = 0;
  fFlux = 0;
  fdT = 0;
  fMeanFlux = 0;

  // FFT Components
  fFFT = 0;
  fRe = 0;
  fIm = 0;
  fAmp = 0;
  fPhi = 0;
  fPSD = 0;
  fOmega = 0;
  // fOmegaMax = 0;
  fOmegaMax = 100;
  fEven = false;
  fNlim = 0;

  fRand = new TRandom3(0);


}

PSDTools::~PSDTools ()
{

  ClearObjs();
  delete []fTime;
  delete []fFlux;


  delete fPSDModel;
  delete fPDF;
  delete fRand;

}



/*
  Set the Model for the PSD to generate from:
  0 - Power Law
  1 - Log Parabola
  2 - Broken Power Law

*/
void PSDTools::SetModel(int modelID)
{

  if (fPSDModel)
  {
    delete fPSDModel;
  }

  // Power Law Model
  if (modelID == 0)
  {
    fModelID = 0;
    fNParms = 2;
    fPSDModel = new TF1("fPSDModel", "[0] * TMath::Power(x / 1e-2, -1.*[1])",  0, 1e5);
    fPSDModel->SetParameter(0,1.e-3);
    fPSDModel->SetParameter(1,1);

  }
  // Log parabola Model
  else if (modelID == 1)
  {
    fModelID = 1;
    fNParms = 3;
    fPSDModel = new TF1("fPSDModel", "[0] * TMath::Power(x / 1e-2, -1.*[1] + [2] * TMath::Log(x / 1e-2))",  0, 1e5);
    fPSDModel->SetParameter(0,1.e-3);
    fPSDModel->SetParameter(1,1);
    fPSDModel->SetParameter(2,0.1);

  }
  // Broken Power Law Model
  else if (modelID == 2)
  {
    fModelID = 2;
    fNParms = 3;
    fPSDModel = new TF1("fPSDModel", this, &PSDTools::BrokenPowerLaw,  1e-6, 1e5, 3);
    fPSDModel->SetParameter(0,1.e-3);
    fPSDModel->SetParameter(1,1);
    fPSDModel->SetParameter(2,0.);
    // fPSDModel->SetParameter(3,1e-2);


  }
  else
  {
    std::cout << "Model not implemented, defaulting to Power Law" << std::endl;
    SetModel(0);
  }
}


double PSDTools::BrokenPowerLaw(double* x, double* parms)
{
  // parms[0] = norm
  // parms[1] = gamma1
  // parms[2] = gamma2
  // parms[3] = EBreak

  if (x[0] <= fOmegaMax)
  {
    return parms[0] * TMath::Power(x[0] / fOmegaMax, -1. * parms[1]);
  }
  else
  {
    return parms[0] * TMath::Power(x[0] / fOmegaMax, -1. * parms[2]);
  }
}

// Setting the PSD Parameters
void PSDTools::SetModelParameters(int n, double* parms)
{
  for (unsigned int i = 0; i < n; i++)
  {
    SetModelParameter(i, parms[i]);
  }
}

// Set individual parameters
void PSDTools::SetModelParameter(int i, double parm)
{
  fPSDModel->SetParameter(i, parm);
}


// Assign the light curve
void PSDTools::SetLightCurve(int npoints, double dT, double *time, double *flux)
{
  fNpoints = npoints;

  // Check if we need to delete
  if (fTime){delete []fTime;}
  if (fFlux){delete []fFlux;}

  // Create new arrays of appropriate size
  fTime = new double[fNpoints];
  fFlux = new double[fNpoints];
  fdT = dT;

  double ifluxmean = 0;


  for (int i = 0; i < fNpoints; i++)
  {
    fTime[i] = time[i];
    fFlux[i] = flux[i];
    ifluxmean += flux[i];
  }

  ifluxmean /= fNpoints;
  if (fMeanFlux == 0)
  {
    fMeanFlux = ifluxmean;
  }


  CalculatePSD();
}


void PSDTools::CalculatePSD()
{
  // Clear objects just incase
  ClearObjs();

  // Create FFT object
  fFFT = TVirtualFFT::FFT(1, &fNpoints, "R2C M K");
  fFFT->SetPoints(fFlux);
  fFFT->Transform();

  // Get FFT
  fRe = new double[fNpoints];
  fIm = new double[fNpoints];



  fEven = fNpoints % 2 == 0;
  fNlim = fEven ? fNpoints/2  : (fNpoints-1) / 2;

  fAmp = new double[fNpoints];
  fPhi = new double[fNpoints];
  fOmega = new double[fNpoints];
  fPSD = new double[fNpoints];



  // Calculate the PSD, Amplitude and Phase of the FFT
  for (int i = 0; i < fNpoints; i++)
  {
    fFFT->GetPointComplex(i, fRe[i], fIm[i]);
    fAmp[i] = TMath::Sqrt(fRe[i] * fRe[i] + fIm[i] * fIm[i] ) / fNpoints;
    fPhi[i] = TMath::ATan2(fRe[i], fIm[i]);
    fOmega[i] = (i) / fdT / fNpoints;
    fPSD[i] = 2 * fdT * (fRe[i] * fRe[i] + fIm[i] * fIm[i] );
    fPSD[i] /= fMeanFlux * fMeanFlux * fNpoints;
  }
}


/*
  Equation A11 of Emmanoulopoulos et al 2013
*/
double PSDTools::Get2LogL(const double* parms)
{

  for (int i = 0; i < fNParms; i++)
  {
    SetModelParameter(i, parms[i]);
  }

  float psd = 0;
  float model = 0;
  float c = 0;

  int itot = fEven ? fNlim/2 -1: (fNlim-1)/2;


  for (int i = 1; i <= itot; i++)
  {

    // Add small offset so not 0
    model = fPSDModel->Eval(fOmega[i]) + 1e-9;
    psd = fPSD[i] + 1e-9;
    c += TMath::Log(model) + psd / model;
  }
  c *=2;
  if (fEven)
  {
    model = fPSDModel->Eval(fOmega[fNlim/2]) + 1e-9;
    psd = fPSD[fNlim/2] + 1e-9;
    c += 2*(TMath::Log(TMath::Pi() * psd * model) + psd / model );
  }

  return c;
}


// Return the PSD as a TGraph
TGraph* PSDTools::GetPSD()
{
  TGraph *gPSD = 0;
  if (fEven) { gPSD = new TGraph(fNlim/2, &(fOmega[1]), &(fPSD[1])); }
  else { gPSD = new TGraph(int(fNlim/2), &(fOmega[1]), &(fPSD[1])); }
  return gPSD;
}


double* PSDTools::FitPSD()
{
  ROOT::Math::Minimizer* minimum =
  ROOT::Math::Factory::CreateMinimizer("Minuit2", "Minos");

  // set tolerance , etc...
  minimum->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2
  minimum->SetMaxIterations(1000000);  // for GSL
  minimum->SetTolerance(10.);
  minimum->SetPrintLevel(1);
  minimum->SetErrorDef(1.);

  // create function wrapper for minimizer
  // a IMultiGenFunction type
  // std::cout << fModelID << " " << fNParms << std::endl;
  ROOT::Math::Functor f(this ,&PSDTools::Get2LogL, fNParms);

  double *step = new double[fNParms];
  double *variable = new double[fNParms];

  minimum->SetFunction(f);



  if (fModelID == 0)
  {
    variable[0] = 1;
    step[0] = 1e-3;
    variable[1] = 1.0;
    step[1] = 0.1;

    // Is limited variable needed?
    minimum->SetLimitedVariable(0, "Norm", variable[0], step[0], 0, 1e5);
    // minimum->SetVariable(1, "Beta", variable[1], step[1]);
    minimum->SetLimitedVariable(1, "Beta", variable[0], step[0], 0, 5);
  }
  else if (fModelID == 2)
  {
    variable[0] = fPSD[int(fNpoints/6)];
    step[0] = 1.e-2 * fPSD[int(fNpoints/6)];
    variable[1] = 1.;
    step[1] = 1e-3;
    variable[2] = 0.;
    step[2] = 1e-3;
    // variable[3] = 1e-2;
    // step[3] = 1e-3;

    minimum->SetLimitedVariable(0, "Norm", variable[0], step[0], 0, 1e1);
    minimum->SetVariable(1, "gamma1", variable[1], step[1]);
    // Expect gamma2 to be close to 0
    // minimum->SetVariable(2, "gamma2", variable[2], step[2]);
    minimum->SetLimitedVariable(2, "gamma2", variable[2], step[2], -1e-3, 1e-3);
    // We know that we aren't seeing less than 1 day sampling
    // minimum->SetLimitedVariable(3, "FBreak", variable[3], step[3], 1e-3, 1e1);
  }


  // do the minimization
  minimum->Minimize();

  const double *xs = minimum->X();
  fFitStatus = minimum->Status();

  double *ret = new double[fNParms];
  for (int i =0; i < fNParms; i++)
  {
    ret[i] = xs[i];
  }
  delete minimum;
  delete [] step;
  delete [] variable;
  // delete xs;
  return ret;

}


int PSDTools::GetFitStatus()
{
  return fFitStatus;
}

// Get a copy of the amplitude and phase
// void PSDTools::GetAmpPhi(double *iAmp, double *iPhi)
void PSDTools::GetAmpPhi(std::vector <double> &iAmp, std::vector <double> &iPhi)
{
  // Make sure we're working with the same length
  // if (iAmp){delete []iAmp;}
  // if (iPhi){delete []iPhi;}
  // iAmp = new double[fNlim];
  // iPhi = new double[fNlim];
  iAmp.clear();
  iAmp.assign(fNlim, 0);
  iPhi.clear();
  iPhi.assign(fNlim, 0);

  for (int i = 0; i < fNlim; i++)
  {
    iAmp[i] = fAmp[i];
    iPhi[i] = fPhi[i];
    // std::cout << i  << " " << fAmp[i] << " " << fPhi[i] << std::endl;
  }

}




// Clear regularly used pointers
void PSDTools::ClearObjs()
{
  if(fFFT) {delete fFFT;}
  if (fRe) {delete [] fRe;}
  if (fIm) {delete [] fIm;}
  if (fAmp) {delete [] fAmp;}
  if (fPhi) {delete [] fPhi;}
  if (fOmega) {delete [] fOmega;}
  if (fPSD) {delete [] fPSD;}
}
