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
    string filetpcfdata="/home/gongjingyu/gcode/astro/2pcf/outcome/20180228/20180228_tpcfv1data";
    string filetpcf = "/home/gongjingyu/gcode/astro/2pcf/outcome/20180228/20180228_tpcfv1";
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
    galaxy *galaxyarraydata;
    galaxy *galaxyarrayrandom;
    redtor * redtorarray;
    raarraydata = (double*)malloc(sizeof(double)*totalnumberdata);
    decarraydata = (double*)malloc(sizeof(double)*totalnumberdata);
    rarraydata = (double*)malloc(sizeof(double)*totalnumberdata);
    galaxyarraydata = (galaxy*)malloc(sizeof(galaxy)*totalnumberdata);
    galaxyarrayrandom = (galaxy*)malloc(sizeof(galaxy)*totalnumberrandom);
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
        numberrandom += 1;
    }

    galaxy1d *xarraydata,*yarraydata,*zarraydata;
    galaxy1d *xarrayrandom, *yarrayrandom, *zarrayrandom;

    xarraydata = (galaxy1d*)malloc(sizeof(galaxy1d)*totalnumberdata);
    yarraydata = (galaxy1d*)malloc(sizeof(galaxy1d)*totalnumberdata);
    zarraydata = (galaxy1d*)malloc(sizeof(galaxy1d)*totalnumberdata);
    xarrayrandom = (galaxy1d*)malloc(sizeof(galaxy1d)*totalnumberrandom);
    yarrayrandom = (galaxy1d*)malloc(sizeof(galaxy1d)*totalnumberrandom);
    zarrayrandom = (galaxy1d*)malloc(sizeof(galaxy1d)*totalnumberrandom);
    init1darray(xarraydata,galaxyarraydata,totalnumberdata,0);
    init1darray(yarraydata,galaxyarraydata,totalnumberdata,1);
    init1darray(zarraydata,galaxyarraydata,totalnumberdata,2);
    init1darray(xarrayrandom,galaxyarrayrandom,totalnumberrandom,0);
    init1darray(yarrayrandom,galaxyarrayrandom,totalnumberrandom,1);
    init1darray(zarrayrandom,galaxyarrayrandom,totalnumberrandom,2);
/*    for (int i = 0; i < 10; i++)
    {
        cout << xarraydata[i].loc << endl;
    }
    for (int i = 0; i < totalnumber; i++)
    {
        if (xarrayd[i].loc < -10000)
        {
            cout << "wroing" << endl;
        }
    }*/
    quicksortgalaxy1d(xarraydata,galaxyarraydata,totalnumberdata,0,0);
    quicksortgalaxy1d(yarraydata,galaxyarraydata,totalnumberdata,0,1);
    quicksortgalaxy1d(zarraydata,galaxyarraydata,totalnumberdata,0,2);
    quicksortgalaxy1d(xarrayrandom,galaxyarrayrandom,totalnumberrandom,0,0);
    quicksortgalaxy1d(yarrayrandom,galaxyarrayrandom,totalnumberrandom,0,1);
    quicksortgalaxy1d(zarrayrandom,galaxyarrayrandom,totalnumberrandom,0,2);
    /*
    cout << xarraydata[0].loc << " " << xarraydata[totalnumberdata-1].loc << endl;
    cout << yarraydata[0].loc << " " << yarraydata[totalnumberdata-1].loc << endl;
    cout << zarraydata[0].loc << " " << zarraydata[totalnumberdata-1].loc << endl;
    cout << xarrayrandom[0].loc << " " << xarrayrandom[totalnumberrandom-1].loc << endl;
    cout << yarrayrandom[0].loc << " " << yarrayrandom[totalnumberrandom-1].loc << endl;
    cout << zarrayrandom[0].loc << " " << zarrayrandom[totalnumberrandom-1].loc << endl;
    */
    /*
    for (int i = 0; i < 20; i++)
    {
        cout << xarraydata[i].loc << endl;
    }
    for (int i = 0; i < 10; i++)
    {
        cout << galaxyarrayd[i].locorder[0] << endl;
    }
    */
    double rprange = 5, rpirange = 5;
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
    calculatetpcf(tpcfdd,rpirange,rprange,rpin,rpn,galaxyarraydata,totalnumberdata,galaxyarraydata,xarraydata,yarraydata,zarraydata,totalnumberdata);
    calculatetpcf(tpcfrr,rpirange,rprange,rpin,rpn,galaxyarrayrandom,totalnumberrandom,galaxyarrayrandom,xarrayrandom,yarrayrandom,zarrayrandom,totalnumberrandom);
    calculatetpcf(tpcfdr,rpirange,rprange,rpin,rpn,galaxyarraydata,totalnumberdata,galaxyarrayrandom,xarrayrandom,yarrayrandom,zarrayrandom,totalnumberrandom);

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
    cout << tpcfdd[0][0] << endl;
    cout << tpcfrr[0][0] << endl;
    cout << tpcfdr[0][0] << endl;
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
    free(xarraydata);free(yarraydata);free(zarraydata);
    free(xarrayrandom);free(yarrayrandom);free(zarrayrandom);
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
