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
    ifstream infile;
    string file="/home/gongjingyu/gcode/astro/2pcf/SDSS7_REAL/SDSS7_real";
    infile.open(file.data());
    long totalnumber = 396068;
    long number  = 0;
    double bg=0.008, ed=0.13;
    long redshiftn = 4000;
    long *numberarray;
    double *raarray, *decarray, *rarray;
    galaxy *galaxyarray;
    redtor * redtorarray;
    numberarray = (long*)malloc(sizeof(long)*totalnumber);
    raarray = (double*)malloc(sizeof(double)*totalnumber);
    decarray = (double*)malloc(sizeof(double)*totalnumber);
    rarray = (double*)malloc(sizeof(double)*totalnumber);
    galaxyarray = (galaxy*)malloc(sizeof(galaxy)*totalnumber);
    redtorarray = (redtor*)malloc(sizeof(redtor)*(redshiftn+1));
    initredtortable(bg,ed,redshiftn,redtorarray,H_0,c);
//    cout << redtorarray[redshiftn-1].r << endl;
//    cout << redtorarray[redshiftn].r << endl;
//    cout << redtorarray[0].r << endl;
//    cout << trapequadrature(0,0.01,200,redshift) << endl;
//    numberarray = new long [totalnumber];
//    raarray = new double [totalnumber];
//    decarray = new double [totalnumber];
//    rarray = new double [totalnumber];
    while (!infile.eof())
    {
        string s;
        getline(infile,s);
        istringstream is(s);
        string str1,str2,str3,str4,str5,str6,str7,str8,str9,str10,str11,str12,str13,str14,str15,str16,str17;
        const char  *str1c, *str2c, *str3c, *str4c, *str5c, *str6c, *str7c, *str8c, *str9c, *str10c, *str11c, *str12c, *str13c, *str14c, *str15c, *str16c, *str17c;
        is >> str1 >> str2 >> str3 >> str4 >> str5 >> str6 >> str7 >> str8 >> str9 >> str10 >> str11 >> str12 >> str13 >> str14 >> str15 >> str16 >> str17;
        str1c = str1.c_str();
        str4c = str4.c_str();
        str5c = str5.c_str();
        str17c = str17.c_str();
        numberarray[number] = atol(str1c);
        raarray[number] = atof(str4c);
        decarray[number] = atof(str5c);
//        rarray[number] = trapequadrature(0,atof(str17c),200,redshift)/H_0;
        rarray[number] = findinredtortable(redtorarray,redshiftn,atof(str17c));
/*        if (number%1000 == 0)
        {
            cout << rarray[number] << endl;
        }*/
        number += 1;
    }
    galaxysphtocar(raarray,decarray,rarray,galaxyarray,totalnumber);

    galaxy1d *xarray,*yarray,*zarray;
    xarray = (galaxy1d*)malloc(sizeof(galaxy1d)*totalnumber);
    yarray = (galaxy1d*)malloc(sizeof(galaxy1d)*totalnumber);
    zarray = (galaxy1d*)malloc(sizeof(galaxy1d)*totalnumber);
    init1darray(xarray,galaxyarray,totalnumber,0);
    init1darray(yarray,galaxyarray,totalnumber,1);
    init1darray(zarray,galaxyarray,totalnumber,2);
    for (int i = 0; i < 10; i++)
    {
        cout << xarray[i].loc << endl;
    }
/*    for (int i = 0; i < totalnumber; i++)
    {
        if (xarray[i].loc < -10000)
        {
            cout << "wroing" << endl;
        }
    }*/
    quicksortgalaxy1d(xarray,galaxyarray,totalnumber,0,0);
    quicksortgalaxy1d(yarray,galaxyarray,totalnumber,0,1);
    quicksortgalaxy1d(zarray,galaxyarray,totalnumber,0,2);
/*    for (int i = 0; i < totalnumber; i++)
    {
        cout << zarray[i].loc << endl;
    }*/
/*    for (int i = 0; i < 10; i++)
    {
        cout << galaxyarray[i].locorder[0] << endl;
    }*/
/*    double rprange = 40, rpirange = 40;
    int rpn  = 40, rpin = 40;
    double **tpcfdd;
    tpcfdd = new double *[rpin*2];
    for (int i = 0; i < rpin*2; i++)
    {
        tpcfdd[i] = new double [rpn*2];
    }
    for (int i = 0; i < rpin*2; i++)
    {
        delete [] tpcfdd[i];
    }
    delete [] tpcfdd;*/
    free(redtorarray);
    free(xarray);free(yarray);free(zarray);
    free(numberarray);
    free(raarray);free(decarray);free(rarray);
    free(galaxyarray);
    infile.close();                             //not necessary
    timeend = time(NULL);
    cout << timeend-timebegin << endl;          //not necessary
    return 0;
}
