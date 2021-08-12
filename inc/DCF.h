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
    std::vector <double> fTimeBinning; // Actual Binning
    double fTimeMin;
    double fTimeMax;

    // DCF Value
    std::vector <double> fDCF;
    std::vector <double> fDCFErr;

    // Light Curve Objects
    LightCurve *fLC1;
    LightCurve *fLC2;

public:

    DCF();
    DCF(LightCurve *iLC1, LightCurve *iLC2);
    DCF(int iNPoints1, double *iMJD1, double *iFlux1, int iNPoints2, double *iMJD2, double *iFlux2);
    DCF(int iNPoints1, double *iMJD1, double *iFlux1, double *iFluxErr1, int iNPoints2, double *iMJD2, double *iFlux2, double *iFluxErr2 );
    
    ~DCF()
    {
        // delete[] fDCF;
        // delete []fDCFErr;
        delete fLC1;
        delete fLC2;

    }

    // Passing light curves on the fly
    // iLC is the light curve number (1 or 2)
    void SetLightCurve(int iNLC, LightCurve *iLC);
    void SetLightCurve(int iNLC, int iNPoints, double *iMJD, double *iFlux);
    void SetLightCurve(int iNLC, int iNPoints, double *iMJD, double *iFlux, double *iFluxErr);

    // Setting time details
    void SetTimeDetails(double iDelT, double iTimeMin, double iTimeMax);
    void GetSubLC(double tau, double delT, LightCurve *&iSub1,  LightCurve *&iSub2);


    TGraphErrors *CalculateDCF(bool bPlotErrors = true);
    

};
#endif
