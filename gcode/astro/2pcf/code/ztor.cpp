#include<iostream>
#include<cstdlib>
#include<cstdio>
#include "tools.h"

#define H_0 100
#define c 299792.458

using namespace std;

int main(int argc,const char *argv[])
{
    if (argc > 0)
    {
        double z = atof(argv[1]);
        double r = trapequadrature(0,z,1000,redshift)/H_0*c;
        cout << r << endl;
    }
    return 0;
}
