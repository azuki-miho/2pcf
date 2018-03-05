#ifndef TOOLS_H
#define TOOLS_H

#include<iostream>
#include<cstdlib>
#include<cmath>
#include<vector>

#define MAX_CAP_NUMBER 20

using namespace std;

struct galaxy1d
{
    double loc;
    long galid;
};

struct galaxy
{
    double xyz[3];
    double radecr[3];
    long locorder[3];
    double luminosity;
};

struct galaxyv2
{
    double xyz[3];
    double radecr[3];
    galaxyv2 *next;
    double luminosity;
};

struct redtor
{
    double red;
    double r;
};

struct rtored
{
    double r;
    double red;
};

void addarraytotpcf(double **tpcf, double rpirange, double rprange, int rpin, int rpn, double x1, double y1, double z1,galaxy *gala2, vector<long> &gals);

void addtotpcfv1(double **tpcf, double rpirange, double rprange, int rpin, int rpn, double rpi, double rp);

void addtotpcfv2(double **tpcf, double rpirange, double rprange, int rpin, int rpn, double rpi, double rp);

double calculateapparentm(double luminosity, double red, double r);

void calculatetpcf(double **tpcf, double rpirange,double rprange, int rpin, int rpn, galaxy *gala1, long n1, galaxy *gala2, galaxy1d *xa2, galaxy1d *ya2, galaxy1d *za2, long n2);

void calculatetpcfv2(double **tpcf, double rpirange,double  rprange,int rpin, int rpn, galaxyv2 *galaxyarray, long totalnumber, galaxyv2 ****box, double xmin, double ymin, double zmin, double xstep, double ystep, double zstep, long xnum, long ynum, long znum);

void filltpcf(double **tpcf, int rpin, int rpn);

void findgalsv1(double x, double y, double z, double rpirange, double rprange, galaxy1d *xa2, galaxy1d *ya2, galaxy1d *za2, long n2,vector<long> &gals);

void findgalsv2(double x, double y, double z, double rpirange, double rprange, galaxy1d *xa2, galaxy1d *ya2, galaxy1d *za2, long n2, galaxy *gala2, vector<long> &gals);

void findgals1d(double x, double radius, galaxy1d *xa, long n, vector<long> &xgals);

long findgals1dlftedge(double x, galaxy1d *xa, long n);

long findgals1drtedge(double x, galaxy1d *xa, long);

double findinlog_Ltable(double *normalprobabilityarray, double log_Lrandom, long n);

double findinredtortable(redtor *redtorarray, long n, double red);

double findinrtoredtable(rtored *rtoredarray, long n, double r);

void galintersect(vector<long> &gals1, vector<long> &gals2, vector<long> &galsi);

void galaxysphtocar(double * raarray, double * decarray, double * rarray, galaxy *galarray, long n);

void galaxysphtocar(double * raarray, double * decarray, double * rarray, galaxyv2 *galarray, long n);

void init1darray(galaxy1d *g1d, galaxy *galarray, long n, int xyzp);

void initlinklist(galaxyv2 ****box, double xmin, double ymin, double zmin, double xstep, double ystep, double zstep, galaxyv2 *galaxyarray, long totalnumber);

void initprobabilityarrayv1(double *probabilityarray, double min_log_L, double delta_log_L, long n, double alpha, double log_L_star);

void initprobabilityarrayv2(double *probabilityarray, double min_log_L, double delta_log_L, long n, double alpha, double log_L_star);

void initredtortable(double bg, double ed, int n, redtor *redtorarray, double H_0, double c);

void inittpcf(double **tpcf, int rpin, int rpn);

void invredtortable(redtor *redtorarray, rtored *rtoredarray, long n, double redstep);

double linearinterpolate(double *xvalue, double *yvalue, long n, double x_0);

void quicksortgalaxy1d(galaxy1d *galo, galaxy *galarray,long n, long befnum, int xyz);

void quicksortgals(vector<long>::iterator bg, long n);

void radecrtoxyz(galaxy *galarr, long n);

double redshift(double z);

double trapequadrature(double bottom, double top, int n, double (*f)(double));

void xyztoradecr(galaxy* galarr, long n);
#endif
