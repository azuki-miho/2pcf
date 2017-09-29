#include<iostream>
#include<cmath>
#include<cstdlib>
#include "tools.h"

using namespace std;

double trapequadrature(double bottom, double top, int n, double (*f)(double))
{
    double step = (top - bottom)/n;
    double sum = 0;
    for (int i = 1; i < n; i++)
    {
        sum += f(bottom+i*step)*step;
    }
    sum += (f(bottom)+f(top))*step/2;
    return sum;
}

double redshift(double z)
{
    double OmegaLambda = 0.718, Omegam = 0.282;
    return 1/pow((OmegaLambda+Omegam*pow(1+z,3)),0.5);
}

double test(double x)
{
    return 1/x;
}

//int main()
//{
//    double integral;
//    integral = trapequadrature(1,10,800,test);
//    cout << integral << endl;
//    return 0;
//}
