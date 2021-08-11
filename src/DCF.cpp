#include "DCF.h"


DCF::DCF()
{
    fNTimeBin = 0;
    fTimeBin = 0;
    fTimeBinning = 0;
    fTimeMin = 0;
    fTimeMax = 0;

    fDCF = 0;
    fDCFErr = 0;
    fLC1 = 0;
    fLC2 = 0;
}


DCF::DCF(LightCurve *iLC1, LightCurve *iLC2)
{
    DCF();
    
    fLC1 = iLC1;
    fLC2 = iLC2;

}


DCF::DCF(int iNPoints1, double *iMJD1, double *iFlux1, int iNPoints2, double *iMJD2, double *iFlux2)
{
    double *iFluxErr1 = new double [iNPoints1];
    double *iFluxErr2 = new double [iNPoints2];
    DCF(iNPoints1, iMJD1, iFlux1, iFluxErr1, iNPoints2, iMJD2, iFlux2, iFluxErr2 );

    delete []iFluxErr1;
    delete []iFluxErr2;
}

DCF::DCF(int iNPoints1, double *iMJD1, double *iFlux1, double *iFluxErr1, int iNPoints2, double *iMJD2, double *iFlux2, double *iFluxErr2 )
{
    // Create light curve objects
    std::vector<double> vMJD1(iNPoints1);
    std::vector<double> vFlux1(iNPoints1);
    std::vector<double> vFluxErr1(iNPoints1);
    LightCurve *iLC1 = new LightCurve( vMJD1, vFlux1, vFluxErr1);

    std::vector<double> vMJD2(iNPoints2);
    std::vector<double> vFlux2(iNPoints2);
    std::vector<double> vFluxErr2(iNPoints2);
    LightCurve *iLC2 = new LightCurve( vMJD2, vFlux2, vFluxErr2);

    DCF(iLC1, iLC2);
}



void DCF::SetLightCurve(int iNLC, LightCurve* iLC)
{
    // ToDo vector<LightCurve>...
    if (iNLC == 1)
    {
        if (fLC1){delete fLC1;}
        fLC1 = iLC;
    }

    if (iNLC == 2)
    {
        if (fLC2){delete fLC2;}
        fLC2 = iLC;
    }
}

void DCF::SetLightCurve(int iNLC, int iNPoints,  double *iMJD, double *iFlux)
{
    double *iFluxErr = new double[iNPoints];
    SetLightCurve( iNLC, iNPoints, iMJD, iFlux, iFluxErr);
    delete []iFluxErr;

}

void DCF::SetLightCurve(int iNLC, int iNPoints, double *iMJD, double *iFlux, double *iFluxErr)
{
    std::vector <double> vMJD(iNPoints);
    std::vector <double> vFlux(iNPoints);
    std::vector <double> vFluxErr(iNPoints);
    for (int i =0; i < iNPoints; i++)
    {
        vMJD[i] = iMJD[i];
        vFlux[i] = iFlux[i];
        vFluxErr[i] = iFluxErr[i];
    }
    LightCurve *iLC = new LightCurve( vMJD, vFlux, vFluxErr);
    SetLightCurve(iNLC, iLC);

}


/*
    Return a subset LC matching time constraints
    There is definitly a more elegant solution
*/
void DCF::GetSubLC(double tau, double delT, LightCurve *&iSub1,  LightCurve *&iSub2)
{

    // Grab a copy of vectors
    std::vector <double> iMJD1_org = fLC1->GetMJD();
    std::vector <double> iFlux1_org = fLC1->GetFlux();
    std::vector <double> iFluxErr1_org = fLC1->GetFluxErr();
    std::vector <double> iMJD2_org = fLC2->GetMJD();
    std::vector <double> iFlux2_org = fLC2->GetFlux();
    std::vector <double> iFluxErr2_org = fLC2->GetFluxErr();
    

    // Create new vectors to work with
    std::vector <double> iMJD1;
    std::vector <double> iFlux1;
    std::vector <double> iFluxErr1;
    
    std::vector <double> iMJD2;
    std::vector <double> iFlux2;
    std::vector <double> iFluxErr2;
    
    // std::cout << "Setting LCs: " << fLC1->fNPoints << " " << fLC2->fNPoints << std::endl;

    // std::cout << "\t\t" << iMJD1.size() << " " << iMJD2.size() << std::endl;
    for (int i = 0; i < iMJD1_org.size(); i++)
    {
        for (int j = 0; j < iMJD2_org.size(); j++)
        {
            // Within the defined time range
            // std::cout << iMJD1_org[i] <<  " " <<  iMJD2_org[j] << tau << abs(iMJD1_org[i] - iMJD2_org[j] - tau)  << std::endl;
            if (abs(iMJD1_org[i] - iMJD2_org[j] - tau) < delT)
            {
                iMJD1.push_back(iMJD1_org[i]);
                iFlux1.push_back(iFlux1_org[i]);
                iFluxErr1.push_back(iFluxErr1_org[i]);

                iMJD2.push_back(iMJD2_org[j]);
                iFlux2.push_back(iFlux2_org[j]);
                iFluxErr2.push_back(iFluxErr2_org[j]);
            }

        }
    }

    // Create new light curves with the new subsets
    // if (iSub1) {delete iSub1;}
    iSub1 = new LightCurve(iMJD1, iFlux1, iFluxErr1);
    // if (iSub2) {delete iSub2;}
    iSub2 = new LightCurve(iMJD2, iFlux2, iFluxErr2);


    // std::cout << iSub1 << " " << iSub2 << std::endl;
    
    
}


// Set up the timing arrays
void DCF::SetTimeDetails(double iDelT, double iTimeMin, double iTimeMax)
{
    fTimeBin = iDelT;
    fTimeMin = iTimeMin;
    fTimeMax = iTimeMax;

    // Round up...
    fNTimeBin = std::ceil((fTimeMax - fTimeMin) / fTimeBin);

    // if (fTimeBinning){delete fTimeBinning;}
    fTimeBinning = new double[fNTimeBin];
    
    for (int i = 0; i < fNTimeBin; i++)
    {
        // std::cout << i << " " << fTimeMin + i*fTimeBin << std::endl;

        fTimeBinning[i] = fTimeMin + i*fTimeBin;
    }

}


TGraphErrors *DCF::CalculateDCF()
{
    // Delete any existing arrays
    // std::cout << "Deleting... " << std::endl;
    // if (fDCF) {delete []fDCF;}
    // if (fDCFErr) {delete []fDCFErr;}
    // std::cout << "Done... " << std::endl;

    fDCF = new double[fNTimeBin];
    fDCFErr = new double[fNTimeBin];

    LightCurve *iSub1 = 0;
    LightCurve *iSub2 = 0;

    double* iUDCF = 0;

    for (int i = 0; i < fNTimeBin; i++)
    {
        fDCF[i] = 0;
        fDCFErr[i] = 0;

        // Get the half bin
        GetSubLC(fTimeBinning[i], 0.5*fTimeBin, iSub1, iSub2);

        std::vector <double> iFlux1 = iSub1->GetFlux();
        std::vector <double> iFlux2 = iSub2->GetFlux();
        std::vector <double> iFluxErr1 = iSub1->GetFluxErr();
        std::vector <double> iFluxErr2 = iSub2->GetFluxErr();


        double iMeanFlux1 = iSub1->GetFluxMean();
        double iMeanFlux2 = iSub2->GetFluxMean();
        double iStdFlux1 = iSub1->GetFluxSTD();
        double iStdFlux2 = iSub2->GetFluxSTD();
        double iMeanFluxErr1 = iSub1->GetFluxErrMean();
        double iMeanFluxErr2 = iSub2->GetFluxErrMean();


        // std::cout << i << " " << fTimeBinning[i] << " " <<  0.5*fTimeBin << " " << iFlux1.size() << std::endl;

        if (iFlux1.size() < 5){continue;}

        // Calculate the DCF
        // Binned DCF values
        // if (iUDCF) {delete []iUDCF;}
        iUDCF = new double[iFlux1.size()];

        for (int j = 0; j < iFlux1.size(); j++)
        {
            iUDCF[j] = (iFlux1[j] - iMeanFlux1) * (iFlux2[j] - iMeanFlux2) \
                        / sqrt( (pow(iStdFlux1,2) - pow(iMeanFluxErr1,2)) * (pow(iStdFlux2,2) - pow(iMeanFluxErr2,2)) );
                        // / sqrt( pow(iStdFlux1,2) * pow(iStdFlux2,2) );

            fDCF[i] += iUDCF[j];
            // std::cout << "\t\t" << iUDCF[j] << " " << (iFlux1[j] - iMeanFlux1) * (iFlux2[j] - iMeanFlux2)
            //         //   << " " << (pow(iStdFlux1,2) - pow(iMeanFluxErr1,2)) * (pow(iStdFlux2,2) - pow(iMeanFluxErr2,2))
            //           << " " << sqrt( pow(iStdFlux1,2) * pow(iStdFlux2,2) )
            //           << " " << iFlux1[j] << " " << iFlux2[j] 
            //           << std::endl;
        }
        fDCF[i] /= iFlux1.size();

        // Calculate the error on DCF
        for (int j = 0; j < iFlux1.size(); j++ )
        {
            fDCFErr[i] += pow((iUDCF[j] - fDCF[i]),2);
        }
        fDCFErr[i] = sqrt(fDCFErr[i]) / (iFlux1.size() -1);

        // std::cout << i << " " << fTimeBinning[i] << " " << fDCF[i] << " " << fDCFErr[i] 
        //           << " " << iSub1->GetFluxMean() << " " << iSub2->GetFluxMean()
        //           << " " << iFlux1.size() << " " << iFlux2.size()
        //           << std::endl;

    }


    TGraphErrors *gDCF = new TGraphErrors(fNTimeBin, fTimeBinning, fDCF, 0, fDCFErr);
    return gDCF;
}



DCF::~DCF()
{
    delete []fTimeBinning;
    delete []fTimeBinning;
    delete []fDCF;
    delete []fTimeBinning;
    delete []fDCFErr;
    delete fLC1;
    delete fLC2;
}
