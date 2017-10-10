#include<iostream>
#include<fstream>
#include<cstring>
#include<sstream>
#include<cstdlib>
#include<ctime>
#include "tools.h"

#define H_0 69.7

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
    double bg=0.01, ed=0.12;
    long redshiftn = 2000;
    long *numberarray;
    double *raarray, *decarray, *rarray;
    galaxy *galaxyarray;
    redtor * redtorarray;
    numberarray = (long*)malloc(sizeof(long)*totalnumber);
    raarray = (double*)malloc(sizeof(double)*totalnumber);
    decarray = (double*)malloc(sizeof(double)*totalnumber);
    rarray = (double*)malloc(sizeof(double)*totalnumber);
    galaxyarray = (galaxy*)malloc(sizeof(galaxy)*totalnumber);
    redtorarray = (redtor*)malloc(sizeof(redtor)*totalnumber);
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
        rarray[number] = findinredtortable(redtorarray,redshiftn,atof(str17c));
/*        if (number%1000 == 0)
        {
            cout << rarray[number] << endl;
        }*/
        number += 1;
    }
    galaxysphtocar(raarray,decarray,rarray,galaxyarray,totalnumber);

    galaxy1d *x1d,*y1d,*z1d;
    x1d = (galaxy1d*)malloc(sizeof(galaxy1d)*totalnumber);
    y1d = (galaxy1d*)malloc(sizeof(galaxy1d)*totalnumber);
    z1d = (galaxy1d*)malloc(sizeof(galaxy1d)*totalnumber);
    init1darray(x1d,galaxyarray,totalnumber,0);
    init1darray(y1d,galaxyarray,totalnumber,1);
    init1darray(z1d,galaxyarray,totalnumber,2);
/*    for (int i = 0; i < 10; i++)
    {
        cout << x1d[i].loc << endl;
    }*/
    quicksortgalaxy1d(x1d,galaxyarray,totalnumber,0,0);
    quicksortgalaxy1d(y1d,galaxyarray,totalnumber,0,1);
    quicksortgalaxy1d(z1d,galaxyarray,totalnumber,0,2);
/*    for (int i = 0; i < totalnumber; i++)
    {
        cout << x1d[i].loc << endl;
    }*/
/*    for (int i = 0; i < 10; i++)
    {
        cout << galaxyarray[i].locorder[0] << endl;
    }*/
    free(x1d);free(y1d);free(z1d);
    free(numberarray);
    free(raarray);free(decarray);free(rarray);
    free(galaxyarray);
    infile.close();                             //not necessary
    timeend = time(NULL);
    cout << timeend-timebegin << endl;          //not necessary
    return 0;
}
