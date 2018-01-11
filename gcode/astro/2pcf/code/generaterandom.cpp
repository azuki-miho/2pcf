#include<iostream>
#include<fstream>
#include<cmath>
#include<ctime>
#include<cstdlib>
#include "tools.h"

#define H_0 100
#define c 299792.458

using namespace std;

int main()
{
//define the variable
    time_t timebegin, timeend;
    timebegin = time(NULL);
    long n=5000;
    long galaxynumber;
    long redshiftn = 4000;
    double alpha=-1.117, log_L_star=10.095, phi_star = 0.03167;
    double min_log_L = 7.8, max_log_L = 11.0;
    double delta_log_L = (max_log_L-min_log_L)/n;
    double normalfactor;
    double x_min=0, x_max=800, y_min = 0, y_max = 800, z_min = 0, z_max = 800;
    double diagnal = pow((pow(x_max,2)+pow(y_max,2)+pow(z_max,2)),0.5);
    double volumn = (x_max-x_min)*(y_max-y_min)*(z_max-z_min);
    double bg = 0.008, ed = 0.13;
    double obs_x = x_max/2, obs_y = y_max/2, obs_z = z_max/2;
    double xinobs, yinobs, zinobs, rinobs, redinobs, apparentm;
//    long randomaccuracy=10000;
    double probabilityarray[n+1];               //This shows the intergral probability for logL
    double normalprobabilityarray[n+1];         //This is the normal probability array
    galaxy *galaxyarray;
    redtor *redtorarray;
    rtored *rtoredarray;
    ofstream outfilemyrandom, outfilemyrandomreadme;
    string filemyrandom = "/home/gongjingyu/gcode/astro/2pcf/MY_RANDOM/myrandom";
    string filemyrandomreadme = "/home/gongjingyu/gcode/astro/2pcf/MY_RANDOM/readme";
//initiate
    redtorarray = (redtor*)malloc(sizeof(redtor)*(redshiftn+1));
    rtoredarray = (rtored*)malloc(sizeof(rtored)*(redshiftn+1));
    initredtortable(bg,ed,redshiftn,redtorarray,H_0,c);
    invredtortable(redtorarray,rtoredarray,redshiftn,(ed-bg)/redshiftn);
    outfilemyrandom.open(filemyrandom.data());
    outfilemyrandomreadme.open(filemyrandomreadme.data());
    for (long i = 0; i < n+1; i++)
    {
        probabilityarray[i]  = 0;
    }
    //has multiply all the things but not normalization
    initprobabilityarrayv1(probabilityarray,min_log_L,delta_log_L,n,alpha,log_L_star);
    galaxynumber = floor(probabilityarray[n]*phi_star*volumn);
    cout << galaxynumber << endl;
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
    long visiblen = 0;
    for (long i = 0; i < galaxynumber; i++)
    {
        xinobs = galaxyarray[i].xyz[0] - obs_x;
        yinobs = galaxyarray[i].xyz[1] - obs_y;
        zinobs = galaxyarray[i].xyz[2] - obs_z;
        rinobs = pow((pow(xinobs,2)+pow(yinobs,2)+pow(zinobs,2)),0.5);
        if (rinobs > rtoredarray[0].r && rinobs < rtoredarray[redshiftn].r)
        {
            redinobs = findinrtoredtable(rtoredarray,redshiftn,rinobs);
            apparentm = calculateapparentm(galaxyarray[i].luminosity,redinobs,rinobs);
            if (apparentm <= 17.77)
            {
                visiblen += 1;
                outfilemyrandom << visiblen << "  " << xinobs << "  " << yinobs << "  " << zinobs << "  " << apparentm << endl;
            }
        }
    }
    cout << visiblen << endl;
    outfilemyrandomreadme << "galaxy ID:" << endl << "xinobs:" << endl << "yinobs:" << endl << "zinobs:" << "apparent magnitude:" << endl << "coordinate in h^{-1}Mpc";
    outfilemyrandom.close();
    outfilemyrandomreadme.close();
    free(redtorarray);
    free(rtoredarray);
    free(galaxyarray);
    timeend = time(NULL);
    cout << timeend - timebegin << endl;
    return 0;
}
