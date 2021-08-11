/*
  Discreet corelation function

*/
#include "LightCurve.h"

#ifndef __DCF_H_
#define __DCF_H_
class DCF {
private:
    /* data */

    int fNTimeBin; // Number of time bins
    double fTimeBin; // Width of the time bins
    double* fTimeBinning; // Actual Binning
    double fTimeMin;
    double fTimeMax;

    // DCF Value
    double* fDCF;
    double* fDCFErr;

    // Light Curve Objects
    LightCurve *fLC1;
    LightCurve *fLC2;

public:

    DCF();
    DCF(LightCurve *iLC1, LightCurve *iLC2);
    DCF(int iNPoints1, double *iMJD1, double *iFlux1, int iNPoints2, double *iMJD2, double *iFlux2);
    DCF(int iNPoints1, double *iMJD1, double *iFlux1, double *iFluxErr1, int iNPoints2, double *iMJD2, double *iFlux2, double *iFluxErr2 );
    ~DCF();


    // Passing light curves on the fly
    // iLC is the light curve number (1 or 2)
    void SetLightCurve(int iNLC, LightCurve *iLC);
    void SetLightCurve(int iNLC, int iNPoints, double *iMJD, double *iFlux);
    void SetLightCurve(int iNLC, int iNPoints, double *iMJD, double *iFlux, double *iFluxErr);

    // Setting time details
    void SetTimeDetails(double iDelT, double iTimeMin, double iTimeMax);
    void GetSubLC(double tau, double delT, LightCurve *&iSub1,  LightCurve *&iSub2);


    TGraphErrors *CalculateDCF();
    

};
#endif
