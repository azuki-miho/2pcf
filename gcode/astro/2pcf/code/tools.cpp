#include<iostream>
#include<cmath>
#include<cstdlib>
#include<cstdio>
#include "tools.h"

using namespace std;

double findinredtortable(redtor * redtorarray, long n, double z)
{
    long left, right;
    double r;
    left = floor((z-redtorarray[0].red)/(redtorarray[n+1].red-redtorarray[0].red)*n);
    right = ceil((z-redtorarray[0].red)/(redtorarray[n+1].red-redtorarray[0].red)*n);
    if (left == right)
    {
        r = retorarray[left].r;
    }
    else
    {
    }
}

void galaxysphtocar(double * raarray,double * decarray, double * rarray, galaxy *galarray,long n)
{
    for (long i = 0; i < n; i++)
    {
        galarray[i].xyz[0] = rarray[i]*cos(decarray[i]*M_PI/180.)*sin(raarray[i]*M_PI/180.);
        galarray[i].xyz[1] = rarray[i]*cos(decarray[i]*M_PI/180.)*cos(raarray[i]*M_PI/180.);
        galarray[i].xyz[2] = rarray[i]*sin(decarray[i]*M_PI/180.);
    }
    return;
}

void init1darray(galaxy1d *g1d, galaxy *galarray, long n, int xyzp)
{
    for (int i = 0; i < n; i++)
    {
        g1d[i].loc = galarray[i].xyz[xyzp];
        g1d[i].galid = i;
        galarray[i].locorder[xyzp] = i;
    }
    return;
}

void initredtortable(double bg, double ed, int n, redtor *,double H_0)
{
    double interval, step;
    interval = ed - bg;
    step = interval/n;
    long intergraln;
    intergraln = ceil(bg/step);
    for (int i = 0; i < n+1; i++)
    {
        redtor[i].red = bg + i * step;
        if (i == 0)
        {
            redtor[i].r = trapequadrature(0,bg,intergraln,redshift)/H_0;

        else
        {
            redtor[i].r = redtor[i-1].r + step * redshift(bg+i*step)/H_0;
        }
    }
    return;
}

void quicksortgalaxy1d(galaxy1d *galo, galaxy *galarray, long n, long befnum, int xyzp)
{
    if (n >= 2)
    {
        int d = 0;
        for (int i = 1; i < n; i++)
        {
            if (galo[i].loc < galo[0].loc)
            {
                if (d != 0)
                {
                    galarray[galo[d].galid].locorder[xyzp] = befnum + i;
                    galarray[galo[i].galid].locorder[xyzp] = befnum + d;
                    galaxy1d temp;
                    temp.loc = galo[i].loc;
                    temp.galid = galo[i].galid;
                    galo[i].loc = galo[d].loc;
                    galo[i].galid = galo[d].galid;
                    galo[d].loc = temp.loc;
                    galo[d].galid = temp.galid;
                    d++;
                }
            }
            else
            {
                if (d == 0)
                {
                    d = i;
                }
            }
        }
        if (d==0)
        {
            d = n;
        }
//        cout << 2 << endl;
        galarray[galo[0].galid].locorder[xyzp] = befnum+d-1;
        galarray[galo[d-1].galid].locorder[xyzp] = befnum;
        galaxy1d temp;
//        cout << 3 << endl;
        temp.loc = galo[d-1].loc;
        temp.galid = galo[d-1].galid;
        galo[d-1].loc = galo[0].loc;
        galo[d-1].galid = galo[0].galid;
        galo[0].loc = temp.loc;
        galo[0].galid = temp.galid;

        quicksortgalaxy1d(galo,galarray,d-1,befnum,xyzp);
        quicksortgalaxy1d(galo+d,galarray,n-d,befnum+d,xyzp);
    }
    return;
}

double redshift(double z)
{
    double OmegaLambda = 0.718, Omegam = 0.282;
    return 1/pow((OmegaLambda+Omegam*pow(1+z,3)),0.5);
}


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

/*
double test(double x)
{
    return 1/x;
}

int main()
{
    double integral;
    integral = trapequadrature(1,10,800,test);
    cout << integral << endl;
    printf("%.8f",M_PI);
    return 0;
}*/
