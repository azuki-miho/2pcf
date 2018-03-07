#include<iostream>
#include<cmath>
#include<fstream>
#include<cstring>
#include<ctime>
#include<cstdlib>
#include<vector>
#include<sstream>
#include<unistd.h>
#include "libpolygon.h"

using namespace std;

int main()
{
    time_t timebegin, timeend;
    timebegin = time(NULL);

    long num_gal = 1086583;
    long num_poly = 653897;
    long num_gal_selected = 0;
    long match_max = 300;
    long max_cap_number = 20;
    long nc1 = 100, nc2 = 100;
    long idx_sel_gal = 0;
    long N_1 = 0, N_2 = 0, N_3 = 0;
    long cn1 = 0, cn2 = 0, cn3 = 0, cn4 = 0;
    long *Npoly_match;
    Npoly_match = (long*)malloc(sizeof(long)*num_gal);
    long *polyid_for_gal;
    polyid_for_gal = (long*)malloc(sizeof(long)*num_gal);
    long *sectorid_for_gal;
    sectorid_for_gal = (long*)malloc(sizeof(long)*num_gal);

    long *isec;
    isec = (long*)malloc(sizeof(long)*num_poly);
    long *linklist;
    linklist = (long*)malloc(sizeof(long)*num_poly);
    long **head_of_chain;
    head_of_chain = (long**)malloc(sizeof(long*)*nc1);
    for (long i = 0; i < nc1; i++)
    {
        head_of_chain[i] = (long*)malloc(sizeof(long)*nc2);
    }
//nc1 for dec and nc2 for ra
    int *cn4_visible;
    cn4_visible = (int*)malloc(sizeof(int)*num_gal);
    long **gal_match_poly;
    double radius = 0.36;
//used for match
    double *poly_racen, *poly_deccen;
    poly_racen = (double*)malloc(sizeof(double)*num_poly);
    poly_deccen = (double*)malloc(sizeof(double)*num_poly);
    double *poly_maglim;
    poly_maglim = (double*)malloc(sizeof(double)*num_poly);
    double *fgot;
    fgot = (double*)malloc(sizeof(double)*num_poly);
    double *galra,*galdec;
    galra = (double*)malloc(sizeof(double)*num_gal);
    galdec = (double*)malloc(sizeof(double)*num_gal);
    double *completeness,*apm_limit;
    completeness = (double*)malloc(sizeof(double)*num_gal);
    apm_limit = (double*)malloc(sizeof(double)*num_gal);
    double *gal_apm;
    gal_apm = (double*)malloc(sizeof(double)*num_gal);
    double **gal_match_poly_dis;
    galv2 *gal_array;
    poly *poly_array;
    gal_match_poly = (long**)malloc(sizeof(long*)*num_gal);
    gal_match_poly_dis = (double**)malloc(sizeof(double*)*num_gal);
    gal_array = (galv2*)malloc(sizeof(galv2)*num_gal);
    poly_array = (poly*)malloc(sizeof(poly)*num_poly);
    for (long i = 0; i < num_gal; i++)
    {
        gal_match_poly[i] = (long*)malloc(sizeof(long)*match_max);
        gal_match_poly_dis[i] = (double*)malloc(sizeof(double)*match_max);
    }
    cout << "finish definition!" << endl;
// initialize all the variable
    for (long i = 0; i < num_gal; i++)
    {
        Npoly_match[i] = 0;
        polyid_for_gal[i] = -1;
        sectorid_for_gal[i] = -1;
        Npoly_match[i] = 0;
        cn4_visible[i] = 0;
        completeness[i] = 0;
        apm_limit[i] = 0;
        for (long j = 0; j < match_max; j++)
        {
            gal_match_poly[i][j] = 0;
            gal_match_poly_dis[i][j] = 0;
        }
    }
    for (long i = 0; i < num_poly; i++)
    {
        linklist[i] = 0;
        poly_array[i].max_cap_number = max_cap_number;
    }
    for (long i = 0; i < nc1; i++)
    {
        for (long j = 0; j < nc2; j++)
        {
            head_of_chain[i][j] = -1;
        }
    }
    cout << "make all the element 0" << endl;
// read in the data of gal sample
    ifstream galfile;
    string galfilestr = "/home/gongjingyu/gcode/astro/2pcf/MY_RANDOM/myrandomSelect";
    galfile.open(galfilestr.data());
    for (long i; i < num_gal; i++)
    {
        string s;
        double x,y,z;
        double ra, dec, r;
        getline(galfile,s);
        istringstream is(s);
        string str1, str2, str3, str4, str5;
        const char *str1c, *str2c, *str3c, *str4c, *str5c;
        is >> str1 >> str2 >> str3 >> str4 >> str5;
        str2c = str2.c_str();
        str3c = str3.c_str();
        str4c = str4.c_str();
        str5c = str5.c_str();
        x = gal_array[i].xyz[0] = atof(str2c);
        y = gal_array[i].xyz[1] = atof(str3c);
        z = gal_array[i].xyz[2] = atof(str4c);
        gal_apm[i] = atof(str5c);
        r = pow(x*x+y*y+z*z,0.5);
        dec = asin(z/r)*180./M_PI;
        ra = asin(abs(y)/pow(x*x+y*y,0.5))*180./M_PI;
        if (x >= 0)
        {
            if (y >= 0)
            {
                ra = ra;
            }
            else
            {
                ra = 360 - ra;
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
                ra = 180 + ra;
            }
        }
        galra[i] = ra;
        galdec[i] = dec;

    }
    cout << "finish read the galfile" << endl;
// read in the data of poly
    ifstream polyfile;
    string polyfilestr = "/home/gongjingyu/gcode/astro/2pcf/polygonfile/lss_combmask.dr72.ply";
    polyfile.open(polyfilestr.data());
    for (long i = 0; i < 1; i++)
    {
        string s;
        getline(polyfile,s);
        istringstream is(s);
        string str1;
        const char *str1c;
        is >> str1;
        str1c = str1.c_str();
        //num_poly = atof(str1c);
    }
    for (long i = 0; i < num_poly; i++)
    {
        string s;
        getline(polyfile,s);
        istringstream is(s);
        string str1, str2, str3, str4, str5, str6, str7, str8;
        const char *str1c, *str2c, *str3c, *str4c, *str5c, *str6c, *str7c, *str8c;
        is >> str1 >> str2 >> str3 >> str4 >> str5 >> str6 >> str7 >> str8;
        str4c = str4.c_str();
        str6c = str6.c_str();
        str8c = str8.c_str();
        poly_array[i].cap_number = atoi(str4c);
        poly_array[i].weight = atof(str6c);
        poly_array[i].str = atof(str8c);
        if (poly_array[i].cap_number > poly_array[i].max_cap_number)
        {
            cout << "larger max_cap_number needed" << endl;
            return 0;
        }
        for (long j = 0; j < poly_array[i].cap_number; j++)
        {
            string st;
            getline(polyfile,st);
            istringstream ist(st);
            string str1t, str2t, str3t, str4t;
            const char *str1ct, *str2ct, *str3ct, *str4ct;
            ist >> str1t >> str2t >> str3t >> str4t;
            str1ct = str1t.c_str();
            str2ct = str2t.c_str();
            str3ct = str3t.c_str();
            str4ct = str4t.c_str();
            //cout << str1ct << "," << str2ct << "," << str3ct << "," << str4ct << endl;
            poly_array[i].caps_xyz[j][0] = atof(str1ct);
            poly_array[i].caps_xyz[j][1] = atof(str2ct);
            poly_array[i].caps_xyz[j][2] = atof(str3ct);
            poly_array[i].cm[j] = atof(str4ct);
        }
    }
    cout << "finish read in the polygonfile" << endl;
// read in the ra dec apm and fgot
    ifstream polydatafile;
    string polydatafilestr = "/home/gongjingyu/gcode/astro/2pcf/polygonfile/lss_combmask.dr72.dat";
    polydatafile.open(polydatafilestr.data());
    for (long i = 0; i < num_poly; i++)
    {
        string s;
        getline(polydatafile,s);
        istringstream is(s);
        string str1, str2, str3, str4, str5;
        const char *str1c, *str2c, *str3c, *str4c, *str5c;
        is >> str1 >> str2 >> str3 >> str4 >> str5;
        str1c = str1.c_str();
        str2c = str2.c_str();
        str3c = str3.c_str();
        str4c = str4.c_str();
        str5c = str5.c_str();
        poly_racen[i] = atof(str1c);
        poly_deccen[i] = atof(str2c);
        isec[i] = atof(str3c);
        poly_maglim[i] = atof(str4c);
        fgot[i] = atof(str5c);
    }
    cout << "finish read in the polygon data file" << endl;
    initgal_match_poly(num_gal,galra,galdec,num_poly,poly_racen,poly_deccen,radius,match_max,Npoly_match,gal_match_poly,gal_match_poly_dis,head_of_chain,nc1,nc2,linklist);
    cout << "finish initalization of gal_match_poly" << endl;
// find the exact polygon ID for every galaxy
    for (long i = 0; i < num_gal; i++)
    {
        if (Npoly_match[i] == 0)
        {
            continue;
        }
        long tempid = match_one_poly(Npoly_match[i],gal_match_poly[i],gal_array,i,poly_array);
        if (tempid >= 0)
        {
            polyid_for_gal[i] = tempid;
            sectorid_for_gal[i] = isec[tempid];
            apm_limit[i] = poly_maglim[tempid];
            completeness[i] = fgot[tempid];
        }
    }
// select whether one galaxy sample can be observed
    ofstream galpolymatchfile;
    string galpolymatchfilestr = "/home/gongjingyu/gcode/astro/2pcf/MY_RANDOM/galpolymatch";
    galpolymatchfile.open(galpolymatchfilestr.data());
    for (long i = 0; i < num_gal; i++)
    {
        if (polyid_for_gal[i] < 0)
        {
            cn1 += 1;
            continue;
        }
        if (gal_apm[i] > apm_limit[i])
        {
            cn2 += 1;
            continue;
        }
        if (completeness[i] < 0.7)
        {
            cn3 += 1;
            continue;
        }
        if (((double)(rand()))/(RAND_MAX) > completeness[i])
        {
            cn4 += 4;
            cn4_visible[i] = -1;
        }
//output the xyz or radecr, cn4visible, apparent magnitude, polygon ID, sector ID, completeness
        num_gal_selected += 1;
        galpolymatchfile << num_gal_selected << "  " << gal_array[i].xyz[0] << "  " << gal_array[i].xyz[1] << "  " << gal_array[i].xyz[2] << "  " << gal_apm[i] << "  " << cn4_visible[i] << "  " << polyid_for_gal[i] << "  " << sectorid_for_gal[i] << "  " << completeness[i] << endl;
    }
    ofstream readmefile;
    string readmefilestr = "/home/gongjingyu/gcode/astro/2pcf/MY_RANDOM/readme_galpolymatch";
    readmefile.open(readmefilestr.data());
    readmefile << "index:" << endl << "x:" << endl << "y:" << endl << "z:" << endl << "apparent magnitude:" << endl << "cn4_visible:" << endl << "polygonID:" << endl << "sectorID:" << endl << "completeness:" << endl <<"total number:" << num_gal_selected << endl;
    galpolymatchfile.close();
    readmefile.close();
    timeend = time(NULL);
    cout << (timeend-timebegin) << endl;
}
