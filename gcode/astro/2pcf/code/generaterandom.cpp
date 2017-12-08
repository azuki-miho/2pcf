#include<iostream>
#include<cmath>
#include<ctime>
#include<cstdlib>
#include "tools.h"

using namespace std;

int main()
{
    long n=5000;
    long galaxynumber;
    double alpha=-1.117, log_L_star=10.095, phi_star = 0.03167;
    double min_log_L = 7.8, max_log_L = 11.0;
    double delta_log_L = (max_log_L-min_log_L)/n;
    double normalfactor;
    double x_min=0, x_max=500, y_min = 0, y_max = 500, z_min = 0, z_max = 500;
    double volumn = (x_max-x_min)*(y_max-y_min)*(z_max-z_min);
//    long randomaccuracy=10000;
    double probabilityarray[n+1];               //This shows the intergral probability for logL
    double normalprobabilityarray[n+1];         //This is the normal probability array
    galaxy *galaxyarray;
    for (int i = 0; i < n+1; i++)
    {
        probabilityarray[i]  = 0;
    }
    //has multiply all the things but not normalization
    initprobabilityarrayv1(probabilityarray,min_log_L,delta_log_L,n,alpha,log_L_star);
    galaxynumber = floor(probabilityarray[n]*phi_star*volumn);
    galaxynumber = 1000;
    //galaxynumber = floor(probabilityarray[n]*log(10)*phi_star*volumn); //suitable for initprobabilityarrayv2
//    cout << galaxynumber << endl;
    normalfactor = probabilityarray[n];
    for (int i = 0; i < n+1; i++)
    {
        normalprobabilityarray[i] = probabilityarray[i]/normalfactor;
    }
    galaxyarray = (galaxy*)malloc(sizeof(galaxy)*galaxynumber);
    srand(time(NULL));
    for (int i = 0; i < galaxynumber; i++)
    {
        double xrandom, yrandom, zrandom, log_Lrandom;
        double log_L;
//        xrandom = (double)(rand()%(randomaccuracy+1))/randomaccuracy;
//        yrandom = (double)(rand()%(randomaccuracy+1))/randomaccuracy;
//        zrandom = (double)(rand()%(randomaccuracy+1))/randomaccuracy;
//        log_Lrandom = (double)(rand()%(randomaccuracy+1))/randomaccuracy;
        xrandom = (double)(rand())/(RAND_MAX);
        yrandom = (double)(rand())/(RAND_MAX);
        zrandom = (double)(rand())/(RAND_MAX);
        log_Lrandom = (double)(rand())/(RAND_MAX);
//        cout << log_Lrandom << endl;
        galaxyarray[i].xyz[0] = x_min + (x_max-x_min)*xrandom;
        galaxyarray[i].xyz[1] = y_min + (y_max-y_min)*yrandom;
        galaxyarray[i].xyz[2] = z_min + (z_max-z_min)*zrandom;
        log_L = min_log_L + findinlog_Ltable(normalprobabilityarray,log_Lrandom,n)*delta_log_L;
        galaxyarray[i].luminosity = pow(10,log_L);
    }
    cout << RAND_MAX << endl;
    free(galaxyarray);
    return 0;
}
