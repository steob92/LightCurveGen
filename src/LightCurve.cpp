#include "LightCurve.h"



// Main constructor that all other instantces will call
LightCurve::LightCurve ( std::vector <double> iMJD, std::vector <double> iFlux, std::vector <double> iFluxErr )
{

    fMJD = 0;
    fFlux = 0;
    fFluxErr = 0;
    
    fNPoints = iMJD.size();
    fMJD = new double[fNPoints];
    fFlux = new double[fNPoints];
    fFluxErr = new double[fNPoints];

    for (int i = 0; i <fNPoints ; i++ )
    {
        fMJD[i] = iMJD[i];
        fFlux[i] = iFlux[i];
        fFluxErr[i] = iFluxErr[i];
    }


    fMJDMin = 99999;
    fMJDMax = -99999;
    fFluxMin = 99999;
    fFluxMax = -99999;
    fFluxMean = 0;
    fFluxSTD = 0;
    fFluxErrMean = 0;

    CalculateProperties();
}


// LightCurve::LightCurve ( std::vector <double> iMJD, std::vector <double> iFlux )
// {
//     std::vector <double> iFluxErr(iMJD.size(), 0 );
//     LightCurve(iMJD, iFlux, iFluxErr);
// }

// LightCurve::LightCurve ( int iNPoints, double *iMJD, double *iFlux )
// {
//     std::vector <double> ivMJD(iNPoints);
//     std::vector <double> ivFlux(iNPoints);

//     for (int i = 0; i < iNPoints; i++)
//     {
//         ivMJD[i] = iMJD[i];
//         ivFlux[i] = iFlux[i];
//     }
//     LightCurve(ivMJD, ivFlux);
// }


// LightCurve::LightCurve ( int iNPoints, double *iMJD, double *iFlux, double* iFluxErr )
// {
//     std::vector <double> ivMJD(iNPoints);
//     std::vector <double> ivFlux(iNPoints);
//     std::vector <double> ivFluxErr(iNPoints);

//     for (int i = 0; i < iNPoints; i++)
//     {
//         ivMJD[i] = iMJD[i];
//         ivFlux[i] = iFlux[i];
//         ivFluxErr[i] = iFluxErr[i];
//     }
//     LightCurve(ivMJD, ivFlux, ivFluxErr);
// }



LightCurve::~LightCurve ()
{
    delete []fMJD;
    delete []fFlux;
    delete []fFluxErr;
}



// Get Some Basic Statistical Properties
void LightCurve::CalculateProperties()
{

    fFluxMean = 0;
    fFluxSTD = 0;
    fFluxErrMean = 0;

    // Making sure the vectors are correctly sorted
    // std::vector <long unsigned int> indx = sort_indexes(fMJD);
    // std::vector <long unsigned int> indx = sort_indices(fMJD);

    std::vector<int> indx(fNPoints);
    std::iota(indx.begin(),indx.end(),0); //Initializing
    std::stable_sort( indx.begin(),indx.end(), [&](int i,int j){return fMJD[i]<fMJD[j];} );

    double* iMJD = new double[indx.size()];
    double* iFlux = new double[indx.size()];
    double* iFluxErr = new double[indx.size()];
    
    // Assigning to temp Vectors
    for (int i = 0; i < indx.size(); i ++)
    {
        iMJD[i] = fMJD[indx[i]];
        iFlux[i] = fFlux[indx[i]];
        iFluxErr[i] = fFluxErr[indx[i]];

        // Calculate the means
        fFluxMean += iFlux[i];
        fFluxErrMean += iFluxErr[i];
    }

    fFluxMean /= fNPoints;
    fFluxErrMean /= fNPoints;


    for (int i = 0; i < fNPoints; i ++)
    {   
        fMJD[i] = iMJD[i];
        fFlux[i] = iFlux[i];
        fFluxErr[i] = iFluxErr[i];

    }

    // Getting Min/Max properties
    // fMJDMin = std::min_element(fMJD.begin(), fMJD.end()));
    // double iMJDMax = *std::max_element(fMJD.begin(), fMJD.end());
    // double iFluxMin = *std::min_element(fFlux.begin(), fFlux.end());
    // double iFluxMax = *std::max_element(fFlux.begin(), fFlux.end());   
    fMJDMin = *std::min_element(fMJD, fMJD + fNPoints);
    fMJDMax = *std::max_element(fMJD, fMJD + fNPoints);
    fFluxMin = *std::min_element(fFlux, fFlux + fNPoints);
    fFluxMax = *std::max_element(fFlux, fFlux + fNPoints);   

    // // fMJDMin = iMJDMin;
    // fMJDMax = iMJDMax;
    // fFluxMin = iFluxMin;
    // fFluxMax = iFluxMax;


    // Get the Standard Deviation
    for (int i = 0; i < fNPoints; i ++){fFluxSTD += (fFlux[i] - fFluxMean) * (fFlux[i] - fFluxMean);}
    fFluxSTD = std::sqrt(fFluxSTD / fNPoints);


}


int LightCurve::PrintDetails()
{

    std::cout << "MJD Min: " <<  fMJDMin << std::endl;
    std::cout << "MJD Max: " <<  fMJDMax << std::endl;
    std::cout << "Flux Min: " <<  fFluxMin << std::endl;
    std::cout << "Flux Max: " <<  fFluxMax << std::endl;
    std::cout << "Flux Mean: " <<  fFluxMean << std::endl;
    std::cout << "Flux STD: " <<  fFluxSTD << std::endl;

    return 0;
}