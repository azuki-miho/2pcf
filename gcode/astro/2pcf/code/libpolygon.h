#ifndef poly_H
#define poly_H

#include<iostream>
#include<cstdlib>
#include<cmath>
#include<vector>

#define MAX_CAP_NUMBER 20

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
    double str;
    long max_cap_number = MAX_CAP_NUMBER;
    long cap_number;
    double caps_xyz[MAX_CAP_NUMBER][3];
    double cm[MAX_CAP_NUMBER];
};

double angdis8(double galra, double galdec, double polyra, double polydec, long para);

double dalp(double stm2, double galdec, long jq1, double hc1, double min_dec, int para);

void initgal_match_poly(long num_gal, double *galra, double *galdec, long num_poly, double *poly_racen, double *poly_deccen,double radius, long match_max, long *Npoly_match, long **gal_match_poly, double **gal_match_poly_dis, long **head_of_chain, long nc1, long nc2, long *linklist);

void inithead_of_chain(long **head_of_chain,long nc1,long nc2,long *linklist,long num_poly,double *poly_racen,double *poly_deccen,double min_ra, double min_dec,double hc1, double hc2);

int in_polygon(poly *poly_array, long polyindx, double x, double y, double z);

long match_one_poly(long Nmatch, long *matchindx, galv2* gal_array, long galindx, poly *poly_array);

void quicksortdist(long nm, double *tdist, long *tindx);


#endif
