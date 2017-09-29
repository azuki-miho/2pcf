#include<iostream>
#include<fstream>
#include<cstring>
#include<sstream>
#include<cstdlib>


using namespace std;

int main()
{
    ifstream infile;
    string file="/home/gongjingyu/gcode/astro/2pcf/SDSS7_REAL/SDSS7_real";
    string s;
    infile.open(file.data());
    long totalnumber = 396068;
    long number  = 0;
    long *numberarray;
    double *raarray, *decarray, *rrrarray;
    numberarray = (long*)malloc(sizeof(long)*totalnumber);
    raarray = (double*)malloc(sizeof(double)*totalnumber);
    decarray = (double*)malloc(sizeof(double)*totalnumber);
    rrrarray = (double*)malloc(sizeof(double)*totalnumber);
    numberarray = new long [totalnumber];
    raarray = new double [totalnumber];
    decarray = new double [totalnumber];
    rrrarray = new double [totalnumber];
    while (!infile.eof())
    {
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
        rrrarray[number] = atof(str17c);
        if (number%1000 == 0)
        {
            cout << rrrarray[number] << endl;
        }
        number += 1;
    }
    cout << number << endl;
    free(numberarray);
    free(raarray);free(decarray);free(rrrarray);
    infile.close();
    return 0;
}
