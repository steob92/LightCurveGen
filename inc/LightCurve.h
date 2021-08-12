/*
    General Light curve class
        Holds data and some statistical properties
        ToDos:
            * Plotting
            * Link with PSDTools?
            * Interpolating
*/

// Standard
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>

// ROOT
#include "TGraphErrors.h"
#include "TCanvas.h"


// Custom
#include "EMP13.h"

#ifndef __LIGHTCURVE_H_
#define __LIGHTCURVE_H_
class LightCurve {
private:
    /* data */
    double *fMJD;
    double *fFlux;
    double *fFluxErr;
    

    // Some Properties of the light curves
    double fMJDMin;
    double fMJDMax;
    double fFluxMin;
    double fFluxMax;
    double fFluxMean;
    double fFluxSTD;
    double fFluxErrMean;
    

public:
    LightCurve ( std::vector <double> iMJD, std::vector <double> iFlux );
    LightCurve ( std::vector <double> iMJD, std::vector <double> iFlux, std::vector <double> iFluxErr );
    LightCurve ( int iNPoints, double *iMJD, double *iFlux );
    LightCurve ( int iNPoints, double *iMJD, double *iFlux, double* iFluxErr );
    
    ~LightCurve ()
    {
        // delete []fMJD;
        // delete []fFlux;
        // delete []fFluxErr;
    }

    int fNPoints;


    // Calculate light curve properties
    void CalculateProperties();

    double GetMJDMin(){return fMJDMin;}
    double GetMJDMax(){return fMJDMax;}

    double GetFluxMin(){return fFluxMin;}
    double GetFluxMax(){return fFluxMax;}

    double GetFluxMean(){return fFluxMean;}
    double GetFluxSTD(){return fFluxSTD;}
    double GetFluxErrMean(){return fFluxErrMean;}
    
    std::vector <double> GetMJD()
    {
        std::vector <double> vMJD(fNPoints);
        for (int i = 0; i < fNPoints; i++){vMJD[i] = fMJD[i];}
        return vMJD;
    }
    std::vector <double> GetFlux()
    {
        std::vector <double> vFlux(fNPoints);
        for (int i = 0; i < fNPoints; i++){vFlux[i] = fFlux[i];}
        return vFlux;
    }
    std::vector <double> GetFluxErr()
    {
        std::vector <double> vFluxErr(fNPoints);
        for (int i = 0; i < fNPoints; i++){vFluxErr[i] = fFluxErr[i];}
        return vFluxErr;
    }

    int PrintDetails();


};
#endif