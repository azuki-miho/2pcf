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
    double x_min=0, x_max=1000, y_min = 0, y_max = 1000, z_min = 0, z_max = 1000;
    long randomaccuracy=10000;
    double probabilityarray[n+1];               //This shows the intergral probability for logL
    double normalprobabilityarray[n+1];         //This is the normal probability array
    galaxy *galaxyarray;
    for (int i = 0; i < n+1; i++)
    {
        probabilityarray[i]  = 0;
    }
    //has multiply all the things but not normalization
    initprobabilityarray(probabilityarray,min_log_L,delta_log_L,n,alpha,log_L_star);
    galaxynumber = floor(probabilityarray[n]);
/*
    normalfactor = probabilityarray[n];
    for (int i = 0; i < n+1; i++)
    {
        normalprobabilityarray[i] = probabilityarray[i]/normalfactor;
    }
    galaxyarray = (galaxy*)malloc(sizeof(galaxy)*galaxynumber);
    for (int i = 0; i < galaxynumber; i++)
    {
        double xrandom, yrandom, zrandom, log_Lrandom;
        double log_L;
        srand(time(0));
        xrandom = (double)(rand()%randomaccuracy)*randomaccuracy;
        srand(time(0));
        yrandom = (double)(rand()%randomaccuracy)*randomaccuracy;
        srand(time(0));
        zrandom = (double)(rand()%randomaccuracy)*randomaccuracy;
        srand(time(0));
        log_Lrandom = (double)(rand()%randomaccuracy)*randomaccuracy;
        galaxyarray[i].xyz[0] = x_min + (x_max-x_min)*xrandom;
        galaxyarray[i].xyz[1] = y_min + (y_max-y_min)*yrandom;
        galaxyarray[i].xyz[2] = z_min + (z_max-z_min)*zrandom;
        log_L = findinlog_Ltable(normalprobabilityarray,log_Lrandom);
	galaxyarray[i].luminosity = pow(10,log_L);
    }
*/
    free(galaxyarray);
    return 0;
}
