#ifndef poly_H
#define poly_H

#include<iostream>
#include<cstdlib>
#include<cmatch>
#include<vector>

using namespace std;


struct galv2
{
    double xyz[3];
    double radecr[3];
    galv2 *next;
    double luminosity;
};

struct poly
{
    double weight;
    long max_cap_number = MAX_CAP_NUMBER;
    long cap_number;
    double caps_xyz[MAX_CAP_NUMER][3];
    double cm[MAX_CAP_NUMBER];
}



#endif
