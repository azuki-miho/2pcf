#include<iostream>
#include<fstream>
#include<cstring>
#include<sstream>
#include<cstdlib>
#include<ctime>
#include<vector>
#include "tools.h"

#define H_0 100
#define c 299792.458

using namespace std;

int main()
{
    time_t timebegin, timeend;            //not necessary
    timebegin = time(NULL);               //not necessary
    ifstream infiledata, infilerandom;
    ofstream outfiletpcfdata, outfiletpcf;
    string filedata="/home/gongjingyu/gcode/astro/2pcf/SDSS7_REAL/SDSS7_real";
    string filerandom = "/home/gongjingyu/gcode/astro/2pcf/SDSS7_RANDOMSELECT/randomSelect";
    string filetpcfdata="/home/gongjingyu/gcode/astro/2pcf/outcome/20180111/20180111_tpcfv2data";
    string filetpcf = "/home/gongjingyu/gcode/astro/2pcf/outcome/20180111/20180111_tpcfv2";
    infiledata.open(filedata.data());
    infilerandom.open(filerandom.data());
    outfiletpcfdata.open(filetpcfdata.data());
    outfiletpcf.open(filetpcf.data());
    long totalnumberdata = 396068;       //acutally 396068
    long totalnumberrandom = 585639;    //acutally 585639
    //long totalnumberrandom = 80000;     acutally 672238
    long numberdata  = 0;
    long numberrandom = 0;
    double bg=0.008, ed=0.13;
    long redshiftn = 4000;
    double *raarraydata, *decarraydata, *rarraydata;
    double xmin_data, xmax_data, ymin_data, ymax_data, zmin_data, zmax_data;
    double xmin_random, xmax_random, ymin_random, ymax_random, zmin_random, zmax_random;
    double tempx, tempy, tempz;
    double box_xstep=10, box_ystep=10, box_zstep=10;
    long box_xnum_data, box_ynum_data, box_znum_data;
    long box_xnum_random, box_ynum_random, box_znum_random;
    galaxyv2 *galaxyarraydata;
    galaxyv2 *galaxyarrayrandom;
    redtor * redtorarray;
    raarraydata = (double*)malloc(sizeof(double)*totalnumberdata);
    decarraydata = (double*)malloc(sizeof(double)*totalnumberdata);
    rarraydata = (double*)malloc(sizeof(double)*totalnumberdata);
    galaxyarraydata = (galaxyv2*)malloc(sizeof(galaxyv2)*totalnumberdata);
    galaxyarrayrandom = (galaxyv2*)malloc(sizeof(galaxyv2)*totalnumberrandom);
    redtorarray = (redtor*)malloc(sizeof(redtor)*(redshiftn+1));
    initredtortable(bg,ed,redshiftn,redtorarray,H_0,c);
//    cout << redtorarray[redshiftn-1].r << endl;
//    cout << redtorarray[redshiftn].r << endl;
//    cout << redtorarray[0].r << endl;
//    cout << trapequadrature(0,0.01,200,redshift) << endl;
//    raarrayd = new double [totalnumber];
//    decarrayd = new double [totalnumber];
//    rarrayd = new double [totalnumber];
    while (/*!infiledata.eof()*/numberdata<totalnumberdata)
    {
        string s;
        getline(infiledata,s);
        istringstream is(s);
        string str1,str2,str3,str4,str5,str6,str7,str8,str9,str10,str11,str12,str13,str14,str15,str16,str17;
        const char  *str1c, *str2c, *str3c, *str4c, *str5c, *str6c, *str7c, *str8c, *str9c, *str10c, *str11c, *str12c, *str13c, *str14c, *str15c, *str16c, *str17c;
        is >> str1 >> str2 >> str3 >> str4 >> str5 >> str6 >> str7 >> str8 >> str9 >> str10 >> str11 >> str12 >> str13 >> str14 >> str15 >> str16 >> str17;
        str1c = str1.c_str();
        str4c = str4.c_str();
        str5c = str5.c_str();
        str17c = str17.c_str();
        raarraydata[numberdata] = atof(str4c);
        decarraydata[numberdata] = atof(str5c);
//        rarrayd[number] = trapequadrature(0,atof(str17c),200,redshift)/H_0;
        rarraydata[numberdata] = findinredtortable(redtorarray,redshiftn,atof(str17c));
        galaxyarraydata[numberdata].next = NULL;
/*        if (number%1000 == 0)
        {
            cout << rarrayd[number] << endl;
        }*/
        numberdata += 1;
    }
    galaxysphtocar(raarraydata,decarraydata,rarraydata,galaxyarraydata,totalnumberdata);
    while (/*!infilerandom.eof()*/numberrandom<totalnumberrandom)
    {
        string s;
        getline(infilerandom,s);
        istringstream is(s);
        string str1,str2,str3,str4;
        const char  *str1c, *str2c, *str3c, *str4c;
        is >> str1 >> str2 >> str3 >> str4;
        str2c = str2.c_str();
        str3c = str3.c_str();
        str4c = str4.c_str();
        galaxyarrayrandom[numberrandom].xyz[0] = atof(str2c);
        galaxyarrayrandom[numberrandom].xyz[1] = atof(str3c);
        galaxyarrayrandom[numberrandom].xyz[2] = atof(str4c);
        galaxyarrayrandom[numberrandom].next = NULL;
        numberrandom += 1;
    }
//decide the range of this box
    xmin_data = xmax_data = galaxyarraydata[0].xyz[0];
    ymin_data = ymax_data = galaxyarraydata[0].xyz[1];
    zmin_data = zmax_data = galaxyarraydata[0].xyz[2];
    for (long i = 1; i < totalnumberdata; i++)
    {
        tempx = galaxyarraydata[i].xyz[0];
        tempy = galaxyarraydata[i].xyz[1];
        tempz = galaxyarraydata[i].xyz[2];
        xmin_data = (tempx < xmin_data) ? tempx : xmin_data;
        xmax_data = (tempx > xmax_data) ? tempx : xmax_data;
        ymin_data = (tempy < ymin_data) ? tempy : ymin_data;
        ymax_data = (tempy > ymax_data) ? tempy : ymax_data;
        zmin_data = (tempz < zmin_data) ? tempz : zmin_data;
        zmax_data = (tempz > zmax_data) ? tempz : zmax_data;
    }
    xmin_data = xmin_data - 1; xmax_data = xmax_data + 1;
    ymin_data = ymin_data - 1; ymax_data = ymax_data + 1;
    zmin_data = zmin_data - 1; zmax_data = zmax_data + 1;
    box_xnum_data = ceil((xmax_data-xmin_data)/box_xstep);
    box_ynum_data = ceil((ymax_data-ymin_data)/box_ystep);
    box_znum_data = ceil((zmax_data-zmin_data)/box_zstep);

    xmin_random = xmax_random = galaxyarrayrandom[0].xyz[0];
    ymin_random = ymax_random = galaxyarrayrandom[0].xyz[1];
    zmin_random = zmax_random = galaxyarrayrandom[0].xyz[2];
    for (long i = 1; i < totalnumberrandom; i++)
    {
        tempx = galaxyarrayrandom[i].xyz[0];
        tempy = galaxyarrayrandom[i].xyz[1];
        tempz = galaxyarrayrandom[i].xyz[2];
        xmin_random = (tempx < xmin_random) ? tempx : xmin_random;
        xmax_random = (tempx > xmax_random) ? tempx : xmax_random;
        ymin_random = (tempy < ymin_random) ? tempy : ymin_random;
        ymax_random = (tempy > ymax_random) ? tempy : ymax_random;
        zmin_random = (tempz < zmin_random) ? tempz : zmin_random;
        zmax_random = (tempz > zmax_random) ? tempz : zmax_random;
    }
    xmin_random = xmin_random - 1; xmax_random = xmax_random + 1;
    ymin_random = ymin_random - 1; ymax_random = ymax_random + 1;
    zmin_random = zmin_random - 1; zmax_random = zmax_random + 1;
    box_xnum_random = ceil((xmax_random-xmin_random)/box_xstep);
    box_ynum_random = ceil((ymax_random-ymin_random)/box_ystep);
    box_znum_random = ceil((zmax_random-zmin_random)/box_zstep);
// generate the linklist
    galaxyv2 ****box_data;
    box_data = (galaxyv2****)malloc(sizeof(galaxyv2***)*box_xnum_data);
    for (long i = 0; i < box_xnum_data; i++)
    {
        box_data[i] = (galaxyv2***)malloc(sizeof(galaxyv2**)*box_ynum_data);
        for (long j = 0; j < box_ynum_data; j++)
        {
            box_data[i][j] = (galaxyv2**)malloc(sizeof(galaxyv2*)*box_znum_data);
            for (long k = 0; k < box_znum_data; k++)
            {
                box_data[i][j][k] = NULL;
            }
        }
    }
    initlinklist(box_data, xmin_data, ymin_data, zmin_data, box_xstep, box_ystep, box_zstep, galaxyarraydata,totalnumberdata);

    galaxyv2 ****box_random;
    box_random = (galaxyv2****)malloc(sizeof(galaxyv2***)*box_xnum_random);
    for (long i = 0; i < box_xnum_random; i++)
    {
        box_random[i] = (galaxyv2***)malloc(sizeof(galaxyv2**)*box_ynum_random);
        for (long j = 0; j < box_ynum_random; j++)
        {
            box_random[i][j] = (galaxyv2**)malloc(sizeof(galaxyv2*)*box_znum_random);
            for (long k = 0; k < box_znum_random; k++)
            {
                box_random[i][j][k] = NULL;
            }
        }
    }

    initlinklist(box_random, xmin_random, ymin_random, zmin_random, box_xstep, box_ystep, box_zstep, galaxyarrayrandom,totalnumberrandom);
    double rprange = 20, rpirange = 20;
    int rpn  = 40, rpin = 40;
    double **tpcfdd;
    double **tpcfrr;
    double **tpcfdr;
    double **tpcf;
    tpcfdd = new double *[rpin*2];
    tpcfrr = new double *[rpin*2];
    tpcfdr = new double *[rpin*2];
    tpcf = new double *[rpin*2];
    for (int i = 0; i < rpin*2; i++)
    {
        tpcfdd[i] = new double [rpn*2];
        tpcfrr[i] = new double [rpn*2];
        tpcfdr[i] = new double [rpn*2];
        tpcf[i] = new double [rpn*2];
    }
    inittpcf(tpcfdd,rpin,rpn);
    inittpcf(tpcfrr,rpin,rpn);
    inittpcf(tpcfdr,rpin,rpn);
    inittpcf(tpcf,rpin,rpn);
// calculate the 2pcf
/*    int check = 0;
    for (long i = 0; i < box_xnum_data; i++)
    {
        for (long j = 0; j < box_ynum_data; j++)
        {
            for (long k = 0; k < box_znum_data; k++)
            {
                if (box_data[i][j][k] != NULL)
                {
                    check = 1;
                }
            }
        }
    }
    cout << check << endl;
*/
    calculatetpcfv2(tpcfdd,rpirange,rprange,rpin,rpn,galaxyarraydata,totalnumberdata,box_data,xmin_data,ymin_data,zmin_data,box_xstep,box_ystep,box_zstep,box_xnum_data,box_ynum_data,box_znum_data);
    calculatetpcfv2(tpcfrr,rpirange,rprange,rpin,rpn,galaxyarrayrandom,totalnumberrandom,box_random,xmin_random,ymin_random,zmin_random,box_xstep,box_ystep,box_zstep,box_xnum_random,box_ynum_random,box_znum_random);
    calculatetpcfv2(tpcfdr,rpirange,rprange,rpin,rpn,galaxyarraydata,totalnumberdata,box_random,xmin_random,ymin_random,zmin_random,box_xstep,box_ystep,box_zstep,box_xnum_random,box_ynum_random,box_znum_random);

    filltpcf(tpcfdd,rpin,rpn);
    filltpcf(tpcfrr,rpin,rpn);
    filltpcf(tpcfdr,rpin,rpn);
    for (int i = 0; i < rpin*2; i++)
    {
        for (int j = 0; j < rpn*2; j++)
        {
            tpcf[i][j] = ((double)tpcfdd[i][j])*((double)tpcfrr[i][j])/(((double)tpcfdr[i][j])*((double)tpcfdr[i][j]))-1;
        }
    }
    /*
    long testn = 0;
    for (int i = 0; i < rpin*2; i++)
    {
        for (int j = 0; j < rpn*2; j++)
        {
            testn += tpcfdd[i][j];
        }
    }
    cout << testn << endl;
    */
    for (int i = 0; i < rpin*2; i++)
    {
        for (int j = 0; j < rpn*2; j++)
        {
            outfiletpcfdata << tpcfdd[i][j];
            if (j!=(rpn*2-1))
            {
                outfiletpcfdata << " ";
            }
        }
        if (i!=(rpin*2-1))
        {
            outfiletpcfdata << endl;
        }
    }
    for (int i = 0; i < rpin*2; i++)
    {
        for (int j = 0; j < rpn*2; j++)
        {
            outfiletpcf << tpcf[i][j];
            if (j!=(rpn*2-1))
            {
                outfiletpcf << " ";
            }
        }
        if (i!=(rpin*2-1))
        {
            outfiletpcf << endl;
        }
    }
    for (int i = 0; i < rpin*2; i++)
    {
        delete [] tpcfdd[i];
        delete [] tpcfrr[i];
        delete [] tpcfdr[i];
        delete [] tpcf[i];
    }
    delete [] tpcfdd;
    delete [] tpcfrr;
    delete [] tpcfdr;
    delete [] tpcf;

    free(redtorarray);
    free(raarraydata);free(decarraydata);free(rarraydata);
    free(galaxyarraydata);
    free(galaxyarrayrandom);
    outfiletpcfdata.close();
    outfiletpcf.close();
    infiledata.close();                             //not necessaryi
    infilerandom.close();
    timeend = time(NULL);
    cout << timeend-timebegin << endl;          //not necessary
    return 0;
}
