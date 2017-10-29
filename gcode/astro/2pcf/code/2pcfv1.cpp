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
    ifstream infiled;
    string filed="/home/gongjingyu/gcode/astro/2pcf/SDSS7_REAL/SDSS7_real";
    infiled.open(filed.data());
    long totalnumberd = 396068;
    long number  = 0;
    double bg=0.008, ed=0.13;
    long redshiftn = 4000;
    long *numberarrayd;
    double *raarrayd, *decarrayd, *rarrayd;
    galaxy *galaxyarrayd;
    redtor * redtorarray;
    numberarrayd = (long*)malloc(sizeof(long)*totalnumberd);
    raarrayd = (double*)malloc(sizeof(double)*totalnumberd);
    decarrayd = (double*)malloc(sizeof(double)*totalnumberd);
    rarrayd = (double*)malloc(sizeof(double)*totalnumberd);
    galaxyarrayd = (galaxy*)malloc(sizeof(galaxy)*totalnumberd);
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
    while (!infiled.eof())
    {
        string s;
        getline(infiled,s);
        istringstream is(s);
        string str1,str2,str3,str4,str5,str6,str7,str8,str9,str10,str11,str12,str13,str14,str15,str16,str17;
        const char  *str1c, *str2c, *str3c, *str4c, *str5c, *str6c, *str7c, *str8c, *str9c, *str10c, *str11c, *str12c, *str13c, *str14c, *str15c, *str16c, *str17c;
        is >> str1 >> str2 >> str3 >> str4 >> str5 >> str6 >> str7 >> str8 >> str9 >> str10 >> str11 >> str12 >> str13 >> str14 >> str15 >> str16 >> str17;
        str1c = str1.c_str();
        str4c = str4.c_str();
        str5c = str5.c_str();
        str17c = str17.c_str();
        numberarrayd[number] = atol(str1c);
        raarrayd[number] = atof(str4c);
        decarrayd[number] = atof(str5c);
//        rarrayd[number] = trapequadrature(0,atof(str17c),200,redshift)/H_0;
        rarrayd[number] = findinredtortable(redtorarray,redshiftn,atof(str17c));
/*        if (number%1000 == 0)
        {
            cout << rarrayd[number] << endl;
        }*/
        number += 1;
    }
    galaxysphtocar(raarrayd,decarrayd,rarrayd,galaxyarrayd,totalnumber);

    galaxy1d *xarrayd,*yarrayd,*zarrayd;
    xarrayd = (galaxy1d*)malloc(sizeof(galaxy1d)*totalnumber);
    yarrayd = (galaxy1d*)malloc(sizeof(galaxy1d)*totalnumber);
    zarrayd = (galaxy1d*)malloc(sizeof(galaxy1d)*totalnumber);
    init1darray(xarrayd,galaxyarrayd,totalnumber,0);
    init1darray(yarrayd,galaxyarrayd,totalnumber,1);
    init1darray(zarrayd,galaxyarrayd,totalnumber,2);
    for (int i = 0; i < 10; i++)
    {
        cout << xarrayd[i].loc << endl;
    }
/*    for (int i = 0; i < totalnumber; i++)
    {
        if (xarrayd[i].loc < -10000)
        {
            cout << "wroing" << endl;
        }
    }*/
    quicksortgalaxy1d(xarrayd,galaxyarrayd,totalnumber,0,0);
    quicksortgalaxy1d(yarrayd,galaxyarrayd,totalnumber,0,1);
    quicksortgalaxy1d(zarrayd,galaxyarrayd,totalnumber,0,2);
    for (int i = 0; i < 20; i++)
    {
        cout << xarrayd[i].loc << endl;
    }
/*    for (int i = 0; i < 10; i++)
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
    calculatetpcf(tpcfdd,rpirange,rprange,rpin,rpn,galaxyarrayd,totalnumberd,galaxyarrayd,xarrayd,yarrayd,zarrayd,totalnumberd);

    filltpcf(tpcfdd,rpin,rpn);
    for (int i = 0; i < rpin*2; i++)
    {
        delete [] tpcfdd[i];
    }
    delete [] tpcfdd;

    free(redtorarray);
    free(xarrayd);free(yarrayd);free(zarrayd);
    free(numberarrayd);
    free(raarrayd);free(decarrayd);free(rarrayd);
    free(galaxyarrayd);
    infile.close();                             //not necessary
    timeend = time(NULL);
    cout << timeend-timebegin << endl;          //not necessary
    return 0;
}
