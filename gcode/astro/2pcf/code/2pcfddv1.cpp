#include<iostream>
#include<fstream>
#include<cstring>
#include<sstream>
#include<cstdlib>
#include<ctime>
#include<vector>
#include "tools.h"

#define H_0 69.7
#define c 299792.458

using namespace std;

int main()
{
    time_t timebegin, timeend;            //not necessary
    timebegin = time(NULL);               //not necessary
    ifstream infiledata;
    ofstream outfiletpcfdata;
    string filedata="/home/gongjingyu/gcode/astro/2pcf/SDSS7_REAL/SDSS7_real";
    string filetpcfdata="/home/gongjingyu/gcode/astro/2pcf/outcome/20171110/20171110_tpcfv2data";
    infiledata.open(filedata.data());
    outfiletpcfdata.open(filetpcfdata.data());
    long totalnumberdata = 80000;       //acutally 396068
    long numberdata  = 0;
    double bg=0.008, ed=0.13;
    long redshiftn = 4000;
    long *numberarraydata;
    double *raarraydata, *decarraydata, *rarraydata;
    galaxy *galaxyarraydata;
    redtor * redtorarray;
    numberarraydata = (long*)malloc(sizeof(long)*totalnumberdata);
    raarraydata = (double*)malloc(sizeof(double)*totalnumberdata);
    decarraydata = (double*)malloc(sizeof(double)*totalnumberdata);
    rarraydata = (double*)malloc(sizeof(double)*totalnumberdata);
    galaxyarraydata = (galaxy*)malloc(sizeof(galaxy)*totalnumberdata);
    redtorarray = (redtor*)malloc(sizeof(redtor)*(redshiftn+1));
    initredtortable(bg,ed,redshiftn,redtorarray,H_0,c);
//    cout << redtorarray[redshiftn-1].r << endl;
//    cout << redtorarray[redshiftn].r << endl;
//    cout << redtorarray[0].r << endl;
//    cout << trapequadrature(0,0.01,200,redshift) << endl;
//    numberarrayd = new long [totalnumber];
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
        numberarraydata[numberdata] = atol(str1c);
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

    galaxy1d *xarraydata,*yarraydata,*zarraydata;
    xarraydata = (galaxy1d*)malloc(sizeof(galaxy1d)*totalnumberdata);
    yarraydata = (galaxy1d*)malloc(sizeof(galaxy1d)*totalnumberdata);
    zarraydata = (galaxy1d*)malloc(sizeof(galaxy1d)*totalnumberdata);
    init1darray(xarraydata,galaxyarraydata,totalnumberdata,0);
    init1darray(yarraydata,galaxyarraydata,totalnumberdata,1);
    init1darray(zarraydata,galaxyarraydata,totalnumberdata,2);
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
/*    for (int i = 0; i < 20; i++)
    {
        cout << xarraydata[i].loc << endl;
    }
    for (int i = 0; i < 10; i++)
    {
        cout << galaxyarrayd[i].locorder[0] << endl;
    }*/
    double rprange = 40, rpirange = 40;
    int rpn  = 40, rpin = 40;
    double **tpcfdd;
    tpcfdd = new double *[rpin*2];
    for (int i = 0; i < rpin*2; i++)
    {
        tpcfdd[i] = new double [rpn*2];
    }
    inittpcf(tpcfdd,rpin,rpn);
    calculatetpcf(tpcfdd,rpirange,rprange,rpin,rpn,galaxyarraydata,totalnumberdata,galaxyarraydata,xarraydata,yarraydata,zarraydata,totalnumberdata);

    filltpcf(tpcfdd,rpin,rpn);
    long testn = 0;
    for (int i = 0; i < rpin*2; i++)
    {
        for (int j = 0; j < rpn*2; j++)
        {
            testn += tpcfdd[i][j];
        }
    }
    cout << testn << endl;
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
        delete [] tpcfdd[i];
    }
    delete [] tpcfdd;

    free(redtorarray);
    free(xarraydata);free(yarraydata);free(zarraydata);
    free(numberarraydata);
    free(raarraydata);free(decarraydata);free(rarraydata);
    free(galaxyarraydata);
    outfiletpcfdata.close();
    infiledata.close();                             //not necessary
    timeend = time(NULL);
    cout << timeend-timebegin << endl;          //not necessary
    return 0;
}
