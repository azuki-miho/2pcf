#include<iostream>
#include<cmath>
#include<cstdlib>
#include<cstdio>
#include<vector>
#include "tools.h"

using namespace std;

void addarraytotpcf(double **tpcf, double rprange, double rpirange, int rpn, int rpin, double x1, double y1, double z1, galaxy *gala2, vector<long> &gals)
{
    for (int i = 0; i < gals.size(); i++)
    {
        double rp, rpi;
        double x2, y2, z2;
        x2 = gala2[gals[i]].xyz[0]; y2 = gala2[gals[i]].xyz[1]; z2 = gala2[gals[i]].xyz[2];
        double lx,ly,lz,sx,sy,sz;
        lx = x2+x1; ly = y2+y1; lz = z2+z1;
        sx = x2-x1; sy = y2-y1; sz = z2-z1;
        rpi = (lx*sx+ly*sy+lz*sz)/sqrt(lx*lx+ly*ly+lz*lz);
        rp = sqrt(sx*sx+sy*sy+sz*sz-rpi*rpi);
        if (rp*rp + rpi*rpi > 0)
        {
            addtotpcf(tpcf,rprange,rpirange,rpn,rpn,rp,rpi);
        }
    }
    return;
}

void addtotpcf(double **tpcf, double rprange, double rpirange, int rpn, int rpin, double rp, double rpi)
{
    if (abs(rp)<rprange && abs(rpi)<rpirange)
    {
        int rpo, rpio;
        rpio = floor(rpi/rprange*rpin)+rpin;
        if (rp < rprange*pow(2,-rpn+1))
        {
            rpo = rpn;
        }
        else
        {
            rpo = floor(log(rp/rprange)/log(2))+2*rpn;
        }
        tpcf[rpio][rpo] += 1;
    }
    return;
}

void calculatetpcf(double **tpcf, double rprange, double rpirange, int rpn, int rpin, galaxy * gala1, long n1, galaxy *gala2, galaxy1d *xa2, galaxy1d *ya2, galaxy1d *za2, long n2)
{
    for (int i = 0; i < n1; i++)
    {
        vector<long> gals;
        findgals(gala1[i].xyz[0],gala1[i].xyz[1],gala1[i].xyz[2],rprange,rpirange,xa2,ya2,za2,n2,gals);
        addarraytotpcf(tpcf,rprange,rpirange,rpn,rpin,gala1[i].xyz[0],gala1[i].xyz[1],gala1[i].xyz[2],gala2,gals);
    }
    return;
}

void findgals(double x, double y, double z, double rprange, double rpirange, galaxy1d *xa2, galaxy1d *ya2, galaxy1d *za2, long n2, vector<long> &gals)
{
    double radius = sqrt(rprange*rprange+rpirange*rpirange);
    vector<long> xgals, ygals, zgals, xygals;
    findgals1d(x,radius,xa2,n2,xgals);
    findgals1d(y,radius,ya2,n2,ygals);
    findgals1d(z,radius,za2,n2,ygals);
    quicksortgals(xgals.begin(),xgals.size());
    quicksortgals(ygals.begin(),ygals.size());
    quicksortgals(zgals.begin(),zgals.size());
    galintersect(xgals,ygals,xygals);
    galintersect(xygals,zgals,gals);
    return;
}

void findgals1d(double x, double radius, galaxy1d *xa, long n, vector<long> &xgals)
{
    long lft, rt;
    lft = findgals1dlftedge(x-radius,xa,n);
    rt = findgals1drtedge(x+radius,xa,n);
    for (long i = lft; i < rt; i++)
    {
        xgals.push_back(xa[i].galid);
    }
    return;
}

long findgals1dlftedge(double x, galaxy1d *xa, long n)
{
    long mid,lft,rt;
    lft = 0; rt = n;
    mid = floor((lft+rt)/2.);
    for (;lft != rt;)
    {
        if (xa[mid].loc >= x)
        {
            rt = mid;
            mid = floor((lft+rt)/2.);
        }
        else
        {
            lft = mid + 1;
            mid = floor((lft+rt)/2.);
        }
    }
    return mid;
}

long findgals1drtedge(double x, galaxy1d *xa, long n)
{
    long mid,lft,rt;
    lft = 0; rt = n;
    mid = floor((lft+rt)/2.);
    for (;lft != rt;)
    {
        if (xa[mid].loc > x)
        {
            rt = mid;
            mid = floor((lft+rt)/2.);
        }
        else
        {
            lft = mid + 1;
            mid = floor((lft+rt)/2.);
        }
    }
    return mid;
}

double findinredtortable(redtor * redtorarray, long n, double z)
{
    long left, right;
    double r;
    left = floor((z-redtorarray[0].red)/(redtorarray[n].red-redtorarray[0].red)*n);
    right = ceil((z-redtorarray[0].red)/(redtorarray[n].red-redtorarray[0].red)*n);
    if (left == right)
    {
        r = redtorarray[left].r;
    }
    else
    {
        r = (z-redtorarray[left].red)/(redtorarray[right].red-redtorarray[left].red)*(redtorarray[right].r-redtorarray[left].r)+redtorarray[left].r;
    }
    return r;
}

void galaxysphtocar(double * raarray,double * decarray, double * rarray, galaxy *galarray,long n)
{
    for (long i = 0; i < n; i++)
    {
        galarray[i].xyz[0] = rarray[i]*cos(decarray[i]*M_PI/180.)*cos(raarray[i]*M_PI/180.);
        galarray[i].xyz[1] = rarray[i]*cos(decarray[i]*M_PI/180.)*sin(raarray[i]*M_PI/180.);
        galarray[i].xyz[2] = rarray[i]*sin(decarray[i]*M_PI/180.);
    }
    return;
}

void galintersect(vector<long> &gals1, vector<long> &gals2, vector<long> &galsi)
{
    long n1, n2;
    n1 = gals1.size(); n2 = gals2.size();
    for (int i1=0,i2=0; i1<n1 && i2<n2; )
    {
        if (gals1[i1] == gals2[i2])
        {
            galsi.push_back(gals1[i1]);
            i1++;i2++;
        }
        else if (gals1[i1] < gals2[i2])
        {
            i1++;
        }
        else
        {
            i2++;
        }
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

void initredtortable(double bg, double ed, int n, redtor *redtorarray,double H_0,double c)
{
    double interval, step;
    interval = ed - bg;
    step = interval/n;
    long intergraln;
    intergraln = ceil(bg/step);
    for (int i = 0; i < n+1; i++)
    {
        redtorarray[i].red = bg + i * step;
        if (i == 0)
        {
            redtorarray[i].r = trapequadrature(0,bg,intergraln,redshift)/H_0*c;
        }
        else
        {
            redtorarray[i].r = redtorarray[i-1].r + step * redshift(bg+i*step)/H_0*c;
        }
    }
    return;
}

void quicksortgalaxy1d(galaxy1d *galo, galaxy *galarray, long n, long befnum, int xyzp)
{
    if (n >= 2)
    {
        long d = 0;
        for (long i = 1; i < n; i++)
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

void quicksortgals(vector<long>::iterator bg, long n)
{
    if (n >= 2)
    {
        long d = 0;
        for (long i = 1; i < n; i++)
        {
            if (bg[i] < bg[0])
            {
                if (d != 0)
                {
                    long temp;
                    temp = bg[i];
                    bg[i] = bg[d];
                    bg[d] = temp;
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
        if (d == 0)
        {
            d = n;
        }
        long temp;
        temp = bg[d-1];
        bg[d-1] = bg[0];
        bg[0] = temp;

        quicksortgals(bg,d-1);
        quicksortgals(bg+d,n-d);
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
