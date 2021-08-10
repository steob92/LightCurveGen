#include "EMP13.h"

EMP13::EMP13() : TK95()
{
  fNpoints = 0;
  fTime = 0;
  fFlux = 0;
  fFFT = 0;
  fNlim = 0;
  fMeanFlux = 0;

  bEven = false;
  fPDF = 0;
  fAmp = 0;
  fPhi = 0;
  fOmega = 0;
  // fPSD = 0;

  fPSDInput = 0;

  fRand = new TRandom3(0);
}

// Assign the light curve
void EMP13::SetLightCurve(int npoints, double dT, double *time, double *flux)
{
  fNpoints = npoints;

  // Check if we need to delete
  if (fTime){delete []fTime;}
  if (fFlux){delete []fFlux;}

  // Create new arrays of appropriate size
  fTime = new double[fNpoints];
  fFlux = new double[fNpoints];
  fdT = dT;

  fMeanFlux = 0;
  for (int i = 0; i < fNpoints; i++)
  {
    fTime[i] = time[i];
    fFlux[i] = flux[i];
    fMeanFlux += flux[i];
  }
  fMeanFlux /= fNpoints;

  // Calculate the PDF

  fPSDInput = new PSDTools();
  fPSDInput->SetModel(0);
  fPSDInput->SetLightCurve( fNpoints, dT, fTime, fFlux);

  SetPDF();

}

// Set the PDF using the passed light curve
// Defining the PDF using the CDF...
// Yes this is counter intuative...
void EMP13::SetPDF()
{

  // Create a new array so we don't overwrite the original data
  double *iflux = new double[fNpoints];
  double *iprob = new double[fNpoints];

  // Create a copy and min/max sort it
  std::copy( fFlux, fFlux + fNpoints, iflux);
  std::sort(iflux, iflux + fNpoints);

  // Creating range in [0,1] seperated by delta
  double delta = 1. / (fNpoints - 1);

  // Obtain the CDF probability
  for (int i = 0; i < fNpoints;  i++) {iprob[i] = delta * i;}

  // Write this to a TGraph
  if (fPDF){delete fPDF;}
  fPDF = new TGraph(fNpoints, iprob, iflux);

  delete []iflux;
  delete []iprob;
}

// Get Get array of random numbers described by the PDF
void EMP13::GetRandomFlux(double *iflux, int npoints)
{
  int nsim = 0;
  if (npoints < 0){nsim = fNpoints;}
  else{nsim = npoints;}
  // Get Random Array [0,1]
  double *irand = new double[nsim];
  fRand->RndmArray(nsim, irand);
  // Evaluate the CDF
  for (int i = 0; i < nsim; i ++ ){iflux[i] = fPDF->Eval(irand[i]);}

  delete []irand;
}



int EMP13::GetRandomLightCurveEMP13(std::vector <double> &iTime, std::vector <double> &iFlux, double dt, int nRed, int npoints)
{
  // Defaulting to passed LC length
  if (npoints == -1){npoints = fNpoints;}

  // Red noise burnin
  int i_nsim = nRed * npoints;

  std::vector <double> i_tk95_time(i_nsim);
  std::vector <double> i_tk95_flux(i_nsim);

  // Step 1.
  // Get a TK95 LC matching the PSD
  // std::cout << "Step 1. \n\t Random TK95" << std::endl;
  GetRandomLightCurve(i_nsim, dt, &(i_tk95_time[0]), &(i_tk95_flux[0]));

  // Step 2.
  // Get random flux values
  // std::cout << "Step 2. \n\t Random Flux"<< std::endl;
  std::vector <double> i_flux_sim(i_nsim);
  std::vector <double> i_flux_sim_previous(i_nsim); // Previous step

  GetRandomFlux(&(i_flux_sim[0]), i_nsim);



  PSDTools *i_flux_psd = new PSDTools();
  std::vector <double> iAmp_flux_sim;
  std::vector <double> iPhi_flux_sim;
  // std::cout << "\tDone" << std::endl;

  // Step 3. FFTs
  // std::cout << "Step 3.\n\t FFTs" << std::endl;
  // std::cout << "\t Set PSD" << std::endl;

  PSDTools *i_tk95_psd = new PSDTools();
  i_tk95_psd->SetLightCurve( i_nsim, dt, &(i_tk95_time[0]), &(i_tk95_flux[0]));
  std::vector <double> iAmp_tk95;
  std::vector <double> iPhi_tk95;



  bool iConverged = false;
  int i_itterations = 0;
  std::cout << "\t FFTs" << std::endl;

  TVirtualFFT* iFFT_TK95 = TVirtualFFT::FFT(1, &i_nsim, "R2C M K");
  TVirtualFFT* iFFT_Flux = TVirtualFFT::FFT(1, &i_nsim, "R2C M K");
  TVirtualFFT* iFFTInverse = TVirtualFFT::FFT(1, &i_nsim, "C2R M K");

  iFFT_TK95->SetPoints(&(i_tk95_flux[0]));
  iFFT_TK95->Transform();


  std::vector <double> i_flux_sim_adj(i_nsim); // adjusted flux


  // Array of indices (used for ranking later)
  std::vector <long unsigned int> indices_sim;
  std::vector <long unsigned int> indices_adj;

  // std::cout << "\tDone" << std::endl;
  //
  // std::cout << "Step 3.5 (looping) FFT Combine" << std::endl;
  while (!iConverged)
  {
    if (i_itterations > 0)
    {
      double delta = 0;
      for (int i = 0; i < i_nsim; i++)
      {
        delta += TMath::Abs(i_flux_sim_previous[i] - i_flux_sim[i]);
        // std::cout << i_flux_sim_previous[i] << " " << i_flux_sim[i] << std::endl;
      }
      delta /= i_nsim;

      std::cout << i_itterations << " " << delta << std::endl;

      if ((delta < 1e-19) || (i_itterations > 1000))
      {
        std::cout << "Converged!!" << std::endl;
        iConverged = true;
        break;
      }
    }

    // Copy Previous
    for (int i = 0 ; i < i_nsim; i ++)
    {
      i_flux_sim_previous[i] = i_flux_sim[i];
    }


    // i_flux_psd->SetLightCurve(i_nsim, dt, &(i_tk95_time[0]), &(i_flux_sim[0]));
    // i_flux_psd->GetAmpPhi(iAmp_flux_sim, iPhi_flux_sim);

    iFFT_Flux->SetPoints(&(i_flux_sim[0]));
    iFFT_Flux->Transform();


    // Get the phase
    for (int i = 0; i < i_nsim; i++)
    {
      double a,b;
      iFFT_TK95->GetPointComplex(i, a, b);

      double x,y;
      iFFT_Flux->GetPointComplex(i, x, y);

      // Step 3. Combine and FFT
      double re, im;
      re = TMath::Sqrt(a*a + b*b) * TMath::Cos(TMath::ATan2(y, x));
      im = TMath::Sqrt(a*a + b*b) * TMath::Sin(TMath::ATan2(y, x));
      // std::cout << i << " " << re << " " << im << std::endl;
      // re = TMath::Sqrt(x*x + y*y) * TMath::Cos(TMath::ATan2(b, a));
      // im = TMath::Sqrt(x*x + y*y) * TMath::Sin(TMath::ATan2(b, a));
      iFFTInverse->SetPoint(i, re, im );
    }

    iFFTInverse->Transform();
    iFFTInverse->GetPoints(&(i_flux_sim_adj[0]));

    indices_sim = sort_indexes(i_flux_sim);
    indices_adj = sort_indexes(i_flux_sim_adj);


    for (int i = 0; i < i_nsim; i ++ )
    {
      // std::cout << i  << " " << i_flux_sim[i] << " " << i_flux_sim_adj[i] << std::endl;
      i_flux_sim_adj[indices_adj[i]] = i_flux_sim[indices_sim[i]];
    }
    for (int i = 0; i< i_nsim; i++)
    {
      i_flux_sim[i] = i_flux_sim_adj[i];
    }

    i_itterations++;

  }

  iTime.clear();
  iTime.resize(i_nsim);
  iFlux.clear();
  iFlux.resize(i_nsim);
  for (int i = 0; i < i_nsim; i ++ )
  {
    iTime[i] = i*dt;
    iFlux[i] = i_flux_sim[i];
  }

  // std::cout << "Exiting..." << std::endl;
  return i_itterations;
}



EMP13::~EMP13 ()
{

  delete []fFreq;
  delete []fRe;
  delete []fIm;
  delete []fPDF;
  delete []fAmp;
  delete []fPhi;
  delete []fOmega;
  delete []fPSD;
  delete []fTime;
  delete []fFlux;

  delete fFFT;
  delete fRand;
  delete fPSDInput;
}
