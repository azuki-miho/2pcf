#ifndef TOOLS_H
#define TOOLS_H

#include<iostream>
#include<cstdlib>
#include<cmath>
#include<vector>

using namespace std;

struct galaxy1d
{
    double loc;
    long galid;
};

struct galaxy
{
    double xyz[3];
    long locorder[3];
    double luminosity;
};

struct redtor
{
    double red;
    double r;
};

void addarraytotpcf(double **tpcf, double rpirange, double rprange, int rpin, int rpn, double x1, double y1, double z1,galaxy *gala2, vector<long> &gals);

void addtotpcfv1(double **tpcf, double rpirange, double rprange, int rpin, int rpn, double rpi, double rp);

void addtotpcfv2(double **tpcf, double rpirange, double rprange, int rpin, int rpn, double rpi, double rp);

void calculatetpcf(double **tpcf, double rpirange,double rprange, int rpin, int rpn, galaxy *gala1, long n1, galaxy *gala2, galaxy1d *xa2, galaxy1d *ya2, galaxy1d *za2, long n2);

void filltpcf(double **tpcf, int rpin, int rpn);

void findgalsv1(double x, double y, double z, double rpirange, double rprange, galaxy1d *xa2, galaxy1d *ya2, galaxy1d *za2, long n2,vector<long> &gals);

void findgalsv2(double x, double y, double z, double rpirange, double rprange, galaxy1d *xa2, galaxy1d *ya2, galaxy1d *za2, long n2, galaxy *gala2, vector<long> &gals);

void findgals1d(double x, double radius, galaxy1d *xa, long n, vector<long> &xgals);

long findgals1dlftedge(double x, galaxy1d *xa, long n);

long findgals1drtedge(double x, galaxy1d *xa, long);

double findinredtortable(redtor * redtorarray, long n, double z);

void galintersect(vector<long> &gals1, vector<long> &gals2, vector<long> &galsi);

void galaxysphtocar(double * raarray, double * decarray, double * rarray, galaxy *galarray, long n);

void init1darray(galaxy1d *g1d, galaxy *galarray, long n, int xyzp);

void initredtortable(double bg, double ed, int n, redtor *redtorarray, double H_0, double c);

void inittpcf(double **tpcf, int rpin, int rpn);

void quicksortgalaxy1d(galaxy1d *galo, galaxy *galarray,long n, long befnum, int xyz);

void quicksortgals(vector<long>::iterator bg, long n);

double redshift(double z);

double trapequadrature(double bottom, double top, int n, double (*f)(double));
#endif
