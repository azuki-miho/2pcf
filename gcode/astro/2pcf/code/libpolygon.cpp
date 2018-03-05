#include<iostream>
#include<cmath>
#include<cstdlib>
#include "libpolygon.h"

using namespace std;

double angdis8(double galra, double galdec, double polyra, double polydec, long para)
{
\\if para is 1, it mean the unit is degree, if it is 0 the unit is rad
    if (para == 1)
    {
        galra = galra/180.*M_PI;
        galdec = galdec/180.*M_PI;
        polyra = polyra/180.*M_PI;
        polydec = polydec/180.*M_PI;
    }
    double theta1 = galdec + 0.5*M_PI;
    double thet12 = polydec + 0.5*M_PI;
    double cosgamma = sin(theta1)*sin(theta2)*cos(galra-polyra)+cos(theta1)*cos(theta2); \\dot product of the two unit vector
    double angdis = 0;
    if (cosgamma < 1)
    {
        angdis = arccos(cosgamma);
    }
    if (para == 1)
    {
        angdis = angdis*180/M_PI;
    }
    return angdis;
}

double dalp(double stm2, double galdec, long jq1, double hc1, double min_dec, int para)
{
    double cdeccl = cos((min_dec + jq1*hc1)/180*M_PI);
    double cdeccu = cos((min_dec + (jq1+1)*hc1)/180*M_PI);
    double cdeccmin;
    double sdalp2;
    double sdalp22;
    double atmp, btmp;
    if (cdeccl < cdeccu)
    {
        cdeccmin = cdeccl;
    }
    else
    {
        cdeccmin = cdeccu;
    }
    if (cdeccmin == 0)
    {
        return 180.;
    }
    if (para == 1)
    {
        sdalp2 = stm2/pow(cos(galdec/180.*M_PI)*cdeccmin);
    }
    else
    {
        atmp = sin((galdec-(min_dec+jq1*hc1)*0.5)/180.*M_PI);
        btmp = sin((galdec-(min_dec+(jq1+1)*hc1)*0.5)/180.*M_PI);
        atmp = atmp*atmp;
        btmp = btmp*btmp;
        if (btmp < atmp)
        {
            atmp = btmp;
        }
        sdalp22 = (stm2*stm2-atmp)/(cos(galdec/180.*M_PI)*cdeccmin);
        if (sdalp22 < 0)
        {
            sdalp22 = 0.;
        }
        sdalp2 = pow(sdalp22,0.5);
    }
    if (sdalp2 > 1)
    {
        return 180.;
    }
    else
    {
        return arcsin(sdalp2)*2.;
    }

}

void initgal_match_poly(long num_gal, double *galra, double *galdec, long num_poly, double *poly_racen, double *poly_deccen,double radius, long match_max, long *Npoly_match, long **gal_match_poly, double **gal_match_poly_dis, long **head_of_chain, long nc1, long nc2, long *linklist)
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
    inithead_of_chain(head_of_chain,nc1,nc2,linklist,num_poly,poly_racen,poly_deccen,min_ra,min_dec,hc1,hc2);

    long iq1, iq2;
    double stm2 = sin(radius/2./180.*M_PI);
    double dltra;
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
            quicksortdist(nm,tdist,tindx);    \\sort_by_dist need to be written
            for (long j = 0; j < nm ; j++)
            {
                gal_match_poly[i][j] = tindx[j];
                gal_match_poly_dist[i][j] = tdist[j];
            }
        }
    }
    return;
}

void inithead_of_chain(long **head_of_chain,long nc1,long nc2,long *linklist,long num_poly,double *poly_racen,double *poly_deccen,double min_ra, double min_dec, ,double hc1, double hc2)
{
    long q1, q2;
    for (long i = 0; i < num_poly; i++)
    {
        q1 = floor((poly_deccen[i]-min_dec)/hc1);
        q2 = floor((poly_racen[i]-min_ra)/hc2);
        if ((q1 >= nc1) || (q1 < 0))
        {
            continue;
        }
        if (q2 >= nc2)
        {
            q2 -= nc2;
        }
        else if(q2 < 0)
        {
            q2 += nc2;
        }
        linklist[i] = head_of_chain[q1][q2];
        head_of_chain[q1][q2] = i;
    }
}

bool in_polygon(polygon *poly_array, long, polyindx, double x, double y, double z)
{
    Ncap = poly_array[polyindx].cap_number;
    for (long i = 0; i < Ncap; i++)
    {
        if (poly_array[polyindx].cm[i] == 0 || poly_array[polyindx].cm[i] < -2.0)
        {
            return false;
        }

    }
    double cm_temp;
    double cmdist_temp;
    double xdist;
    double ydist;
    double zdist;
    for (long i = 0; i < Ncap; i++)
    {
        if (poly_array[polyindx].cm[i] > 2)
        {
            continue;
        }
        cm_temp = fabs(poly_array[polyindx].cm[i]);
        xdist = x - poly_array[polyindx].caps_xyz[i][0];
        ydist = y - poly_array[polyindx].caps_xyz[i][1];
        zdist = z - poly_array[polyindx].caps_xyz[i][2];
        cmdist_temp = (xdist*xdist+ydist*ydist+zdist*zdist)/2.0;
        if (poly_array[polyindx].cm[i] > 0)
        {
            if (cmdist_temp > cmtemp)
            {
                return false;
            }
        }
        else
        {
            if (cmdist_temp < cmtemp)
            {
                return false;
            }
        }
    }
    return true;
}

long match_one_poly(long Nmatch, long *matchidx, galaxyv2* gal_array, long galindx, polygon *poly_array, long polyindx)
{
    for (int i = 0; i < Nmatch ; i++)
    {
        long polyindx = matchindx[i];
        if (in_polygon(poly_array,polyindx,gal_array[galindx].xyz[0],gal_array[galindx].xyz[1],gal_array[galindx].xyz[2]));
        {
            return polyindx;
        }
    }
    return 0;
}

void quicksortdist(long nm, double *tdist, long *tindx)
{
    if (nm == 1)
    {
        return;
    }
    if (nm <= 2)
    {
        if (tdist[0] > tdist[1])
        {
            double tempdist = tdist[0];
            long temdindx = tindx[0];
            tdist[0] = tdist[1];
            tindx[0] = tdist[1];
        }
    }
    else
    {
        long item_location = 0;
        long compare_location = 1;
        for (; compare_location < nm;)
        {
            if (tdist[item_location] < tdist[compare_location])
            {
                compare_location += 1;
            }
            else
            {
                double tempdist = tdist[item_location];
                long tempindx = tindx[item_location];
                tdist[item_location] = tdist[compare_location];
                tindx[item_location] = tindx[compare_location];
                tdist[compare_location] = tdist[item_location+1];
                tindx[compare_location] = tindx[item_location+1];
                tdist[item_location+1] = tempdist;
                tindx[item_location+1] = tempindx;
                item_location += 1;
                compare_location += 1;
            }
        }
        quicksortdist(item_location,tdist,tindx);
        quicksortdist(nm-item_location-1,tdist+item_location+1,tindx+item_location+1);
    }
}
