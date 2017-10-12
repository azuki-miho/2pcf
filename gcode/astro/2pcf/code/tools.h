#include<iostream>
#include<cstdlib>
#include<cmath>

#ifndef TOOLS_H
#define TOOLS_H
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

double findinredtortable(redtor * redtorarray, long n, double z);

void galaxysphtocar(double * raarray, double * decarray, double * rarray, galaxy *galarray, long n);

void init1darray(galaxy1d *g1d, galaxy *galarray, long n, int xyzp);

void initredtortable(double bg, double ed, int n, redtor *redtorarray, double H_0, double c);

void quicksortgalaxy1d(galaxy1d *galo, galaxy *galarray,long n, long befnum, int xyz);

double redshift(double z);

double trapequadrature(double bottom, double top, int n, double (*f)(double));
#endif
