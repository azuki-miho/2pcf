#include<iostream>
#include<cstdlib>
#include<cmath>

using namespace std;

double redshift(double z);

double trapequadrature(double bottom, double top, int n, double (*f)(double));
