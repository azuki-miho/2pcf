#include<iostream>
#include<cstdlib>
#include<cmath>

#ifndef TOOLS_H
#define TOOLS_H
using namespace std;

struct galaxy
{
    galaxy() : x(0),y(0),z(0) {  }
    galaxy(double x0, double y0, double z0) : x(x0),y(y0),z(x0){  }
    double x,y,z;
};

/*inline galaxy::galaxy()
{
    x = 0; y = 0; z = 0;
}*/

/*inline galaxy::galaxy(double x0, double y0, double z0)
{
    x = x0; y = y0; z = z0;
}*/

double redshift(double z);

double trapequadrature(double bottom, double top, int n, double (*f)(double));
#endif
