#include<iostream>
#include<cmath>
#include<fstream>
#include<cstring>
#include<ctime>
#include<cstdlib>
#include "libpoly.h"

using namespace std;

int main()
{
    time_t timebegin, timeend;
    timebegin = time(NULL);

    long num_gal, num_poly, num_gal_selected;
    long match_max = 300;
    long nc1 = 100, nc2 = 100;
    long Npoly_match[num_gal];
    long polyid_for_gal[num_gal];
    long sectorid_for_gal[num_gal];
    long gal_match_poly[num_gal][match_max];
    long head_of_chain[nc1][nc2]
    double radius = 0.36                       \\used for match
    double poly_racen[num_poly], poly_deccen[num_poly];
    double gal_poly[num_gal], galra[num_gal],galdec[num_gal];
    double gal_match_poly_dis[num_gal][match_max];
    galv2 *gal_array;
    poly *poly_array;
    gal_array = (galv2*)malloc(sizeof(galv2)*num_gal);
    poly_array = (poly*)malloc(sizeof(poly)*num_poly);

    cout << "finish definition!" << endl;
\\ initialize all the variable
    for (long i = 0; i < nc1; i++)
    {
        for (long j = 0; i < nc2; j++)
        {
            head_of_chain = 0;
        }
    }
    for (long i = 0; i < num_gal; i++)
    {
        for (long j = 0; j < match_max; j++)
        {
            gal_match_poly[i][j] = 0;
            gal_match_poly_dis[i][j] = 0;
        }
        Npoly_match[i] = 0;
    }
\\ read in the data of gal sample
\\ read in the data of poly
    initgal_match_poly(num_gal,galra,galdec,num_poly,poly_racen,polydeccen,radius,match_max,Npoly_match,gal_match_poly,gal_match_poly_dist,head_of_chain,nc1,nc2,linklist);
    cout << "finish initalization of gal_match_poly" << endl;

    timeend = time(NULL);
}
