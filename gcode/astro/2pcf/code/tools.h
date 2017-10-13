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
};

struct redtor
{
    double red;
    double r;
};

void addarraytotpcf(double **tpcf, int rpn, int rpin, double x1, double y1, double z1,galaxy *gala2, vector<long> &gals);

void addtotpcf(double **tpcf, int rpn, int rpin, double rp, double rpi);

void calculatetpcf(double **tpcf,int rpn, int rpin, galaxy *gala1, long n1, galaxy *gala2, galaxy1d *xa2, galaxy1d *ya2, galaxy1d *za2, long n2);

double findinredtortable(redtor * redtorarray, long n, double z);

void galaxysphtocar(double * raarray, double * decarray, double * rarray, galaxy *galarray, long n);

void init1darray(galaxy1d *g1d, galaxy *galarray, long n, int xyzp);

void initredtortable(double bg, double ed, int n, redtor *redtorarray, double H_0, double c);

void quicksortgalaxy1d(galaxy1d *galo, galaxy *galarray,long n, long befnum, int xyz);

double redshift(double z);

double trapequadrature(double bottom, double top, int n, double (*f)(double));
#endif
