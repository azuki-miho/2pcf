#include<iostream>
#include<cmath>
#include<cstdlib>
#include "libpolygon.h"

using namespace std;

void initgal_match_poly(long num_gal, double *galra, double *galdec, long num_poly, double *poly_racen, double *poly_deccen,double radius, long match_max, long *Npoly_match, long **gal_match_poly, double **gal_match_poly_dis, long **head_of_chain, long nc1, long nc2, long *linklist);
{
    double min_dec = galdec[0];
    double max_dec = galdec[0];
    for (long i = 0; i < num_gal; i++)
    {
        if (min_dec > galdec[i])
        {
            min_dec = galdec[i];
        }
        if (max_dec < galdec[i])
        {
            max_dec = galdec[i];
        }
    }
    min_dec -= 0.1;
    max_dec += 0.1;

    double hc1, hc2;
    hc1 = (max_dec-min_dec)/((double)nc1);
    hc2 = (max_ra-min_ra)/((double)nc2);
    inithead_of_chain(head_of_chain,nc1,nc2,linklist,num_poly,poly_racen,poly_deccen);

    long iq1, iq2;
    double stm2 = sin(radius/2./180.*M_PI);
    for (int i = 0; i < num_gal; i++)
    {
        iq1 = floor((galdec[i]-min_dec)/hc1);
        iq2 = floor((galra[i]-min_ra)/hc2);
        long jq1m = floor(radius/hc1)+1;
        long jq2m;
        for (long jq1 = iq1-jq1m; jq1 <= iq1+jq1m; jq1++)
        {
            if ((jq1 >= nc1)or(jq1 < 0))
            {
                continue;
            }
            if (jq1 = iq1)
            {
                dltra = dalp(stm2,galdec[i],jq1,hc1,min_dec,1);
            }
            else
            {
                dltra = dalp(stm2,galdec[i],jq1,hc1,min_dec,0);
            }
            jq2m = floor(dltra/hc2)+1;
            jq2max = iq2+jq2m;
            jq2min = iq2-jq2m;
            if (jq2max-jq2min+1>nc2)
            {
                jq2max = jq2min-1+nc2
            }
            for (long jq2 = jq2min; jq2 <= jq2max; jq2++)
            {
                jq2t = jq2
                for (;jq2t >= nc2;)
                {
                    jq2t -= nc2;
                }
                for (;jq2t < 0;)
                {
                    jq2t += nc2;
                }
                j = head_of_chain[jq1][jq2t];
                for (;j != 0;)
                {
                    adist = angdis8(gal_ra[i],gal_dec[i],poly_racen[j],poly_deccen[j],1);
                    if (adist < radius)
                    {
                        if (Npoly_match[i] < match_max)
                        {
                            Npoly_match[i] += 1;
                            gal_match_poly[i][Npoly_match[i]-1] = j;
                            gal_match_poly_dist[i][Npoly_match[i]-1] = adist*3600.;
                        }
                        else
                        {
                            cout << "match_max has been reached";
                        }
                    }
                    j = linklist[j];
                }
            }
        }
    }
    long nm;
    double tdist[match_max];
    long isort[match_max];
    long tindx[match_max];
    for (long i = 0; i < num_gal; i++)
    {
        nm = Npoly_match[i];
        for (;nm > 1;)
        {
            for (long j = 0; j < nm ; j++)
            {
                tdist[j] = gal_match_poly_dist[i][j];
                tindx[j] = gal_match_poly[i][j];
            }
            sort_by_dist(nm,tdist,isort);    \\sort_by_dist need to be written
            for (long j = 0; j < nm ; j++)
            {
                gal_match_poly[i][j] = tindx[isort[j]];
                gal_match_poly_dist[i][j] = tdist[isort[j]];
            }
        }
    }
    return;
}
