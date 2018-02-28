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
    long isec[num_poly];
    long gal_match_poly[num_gal][match_max];
    long head_of_chain[nc1][nc2]
    double radius = 0.36                       \\used for match
    double poly_racen[num_poly], poly_deccen[num_poly];
    double maglim[num_poly];
    double fgot[num_poly];
    double galra[num_gal],galdec[num_gal];
    double completeness[num_gal],apm_limit[num_gal];
    double gal_match_poly_dis[num_gal][match_max];
    galv2 *gal_array;
    poly *poly_array;
    gal_array = (galv2*)malloc(sizeof(galv2)*num_gal);
    poly_array = (poly*)malloc(sizeof(poly)*num_poly);

    cout << "finish definition!" << endl;
\\ initialize all the variable
    for (long i = 0; i < num_gal; i++)
    {
        Npoly_match[i] = 0;
        polyid_for_gal[i] = 0;
        sectorid_for_gal[i] = 0;
        Npoly_match[i] = 0;
        completeness[i] = 0;
        apm_limit[i] = 0;
        for (long j = 0; j < match_max)
        {
            gal_match_poly[i][j] = 0;
            gal_match_poly_dis[i][j] = 0;
        }
        
    }
    for (long i = 0; i < nc1; i++)
    {
        for (long j = 0; i < nc2; j++)
        {
            head_of_chain = 0;
        }
    }
    
    initgal_match_poly(num_gal,galra,galdec,num_poly,poly_racen,polydeccen,radius,match_max,Npoly_match,gal_match_poly,gal_match_poly_dist,head_of_chain,nc1,nc2,linklist);
    cout << "finish initalization of gal_match_poly" << endl;
\\ read in the data of gal sample
\\ read in the data of poly
\\ find the exact polygon ID for every galaxy
    for (long i = 0; i < num_gal; i++)
    {
        if (Npoly_match[i] = 0)
        {
            continue;
        }
        long tempid = match_one_poly();
        if (tempid > 0)
        {
            polyid_for_gal[i] = tempid;
            sectorid_for_gal[i] = isec[tempid];
            apm_limit[i] = maglim[tempid];
            completeness[i] = fgot[tempid];
        }
    }
\\ select whether one galaxy sample can be observed
    timeend = time(NULL);
}
