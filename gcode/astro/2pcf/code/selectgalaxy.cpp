#include<iostream>
#include<fstream>
#include<cstring>
#include<sstream>
#include<cstdlib>
#include<vector>
#include<cmath>
#include "tools.h"

#define H_0 69.7
#define c 299892.458

int main()
{
    ifstream infiledata, infilerandom;
    ofstream outfileselect, outfilereadme;
    string filedata="/home/gongjingyu/gcode/astro/2pcf/SDSS7_REAL/SDSS7_real";
    string filerandom = "/home/gongjingyu/gcode/astro/2pcf/SDSS7_RANDOM/randomSamp";
    string fileselect = "/home/gongjingyu/gcode/astro/2pcf/SDSS7_RANDOMSELECT/randomSelect";
    string filereadme = "/home/gongjingyu/gcode/astro/2pcf/SDSS7_RANDOMSELECT/readme";
    infiledata.open(filedata.data());
    infilerandom.open(filerandom.data());
    outfileselect.open(fileselect.data());
    outfilereadme.open((filereadme.data()));
    long totalnumberdata = 396068;       //acutally 396068
    long totalnumberrandom = 672238;     //acutally 672238
    long numberdata  = 0;
    long numberrandom = 0;
    long numberselect = 0;
    double bg=0.008, ed=0.13;
    double ramin=0,ramax=0,decmin=0,decmax=0,rmin=0,rmax=0;
    long redshiftn = 4000;
    redtor * redtorarray;
    redtorarray = (redtor*)malloc(sizeof(redtor)*(redshiftn+1));
    initredtortable(bg,ed,redshiftn,redtorarray,H_0,c);
    while (numberdata<totalnumberdata)
    {
        string s;
        double ra,dec,r;
        getline(infiledata,s);
        istringstream is(s);
        string str1,str2,str3,str4,str5,str6,str7,str8,str9,str10,str11,str12,str13,str14,str15,str16,str17;
        const char  *str1c, *str2c, *str3c, *str4c, *str5c, *str6c, *str7c, *str8c, *str9c, *str10c, *str11c, *str12c, *str13c, *str14c, *str15c, *str16c, *str17c;
        is >> str1 >> str2 >> str3 >> str4 >> str5 >> str6 >> str7 >> str8 >> str9 >> str10 >> str11 >> str12 >> str13 >> str14 >> str15 >> str16 >> str17;
        str1c = str1.c_str();
        str4c = str4.c_str();
        str5c = str5.c_str();
        str17c = str17.c_str();
        ra = atof(str4c);
        dec = atof(str5c);
//        rarrayd[number] = trapequadrature(0,atof(str17c),200,redshift)/H_0;
        r = findinredtortable(redtorarray,redshiftn,atof(str17c));
        if (ra < ramin)
        {
            ramin = ra;
        }
        if (ra > ramax)
        {
            ramax = ra;
        }
        if (dec < decmin)
        {
            decmin = dec;
        }
        if (dec > decmax)
        {
            decmax = dec;
        }
        if (r < rmin)
        {
            rmin = r;
        }
        if (r > rmax)
        {
            rmax = r;
        }
/*        if (number%1000 == 0)
        {
            cout << rarrayd[number] << endl;
        }*/
        numberdata += 1;
    }
    while (/*!infilerandom.eof()*/numberrandom<totalnumberrandom)
    {
        string s;
        double x, y, z;
        double ra, dec, r;
        getline(infilerandom,s);
        istringstream is(s);
        string str1,str2,str3,str4,str5,str6,str7,str8,str9,str10,str11,str12,str13,str14,str15,str16,str17;
        const char  *str1c, *str2c, *str3c, *str4c, *str5c, *str6c, *str7c, *str8c, *str9c, *str10c, *str11c, *str12c, *str13c, *str14c, *str15c, *str16c, *str17c;
        is >> str1 >> str2 >> str3 >> str4 >> str5 >> str6 >> str7 >> str8 >> str9 >> str10 >> str11 >> str12 >> str13 >> str14 >> str15 >> str16 >> str17;
        str6c = str6.c_str();
        str7c = str7.c_str();
        str8c = str8.c_str();
        x = atof(str6c);
        y = atof(str7c);
        z = atof(str8c);
        r = pow(x*x+y*y+z*z,0.5);
        dec = acos(pow(x*x+y*y,0.5)/r)*180./M_PI;
        if (z <= 0)
        {
            dec = -dec;
        }
        ra = asin(abs(x)/pow(x*x+y*y,0.5))*180./M_PI;
        if (x >= 0)
        {
            if (y >= 0)
            {
                ra = ra;
            }
            else
            {
                ra = 360-ra;
            }
        }
        else
        {
            if (y >= 0)
            {
                ra = 180 - ra;
            }
            else
            {
                ra = 180+ra;
            }
        }
        if (ramin < ra && ra < ramax & decmin < dec && dec < decmax && rmin < r && r < rmax)
        {
            numberselect += 1;
            outfileselect  << numberselect << "  " << x << "  " << y << "  " << z << endl;
        }
        numberrandom += 1;
    }
    outfilereadme << "galaxy ID:" << endl << "x:" << endl << "y:" << endl << "z:" << endl << "totalnumber:" << numberselect << endl;
    outfileselect.close();
    outfilereadme.close();
    infiledata.close();
    infilerandom.close();
    free(redtorarray);
    return 0;
}
