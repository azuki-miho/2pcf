#include<iostream>
#include<cmath>
#include<fstream>
#include<cstring>
#include<ctime>
#include<cstdlib>
#include "libpoly.h"

using namespace std;

int main()
{
    time_t timebegin, timeend;
    timebegin = time(NULL);

    long num_gal, num_poly, num_gal_selected;
    long match_max = 300;
    long nc1 = 100, nc2 = 100;
    long idx_sel_gal = 0;
    long N_1 = 0, N_2 = 0, N_3 = 0;
    long cn1 = 0, cn2 = 0, cn3 = 0, cn4 = 0;
    long Npoly_match[num_gal];
    long polyid_for_gal[num_gal];
    long sectorid_for_gal[num_gal];

    long isec[num_poly];
    long linklist[num_poly];
    long gal_match_poly[num_gal][match_max];
    long head_of_chain[nc1][nc2]; \\nc1 for dec and nc2 for ra
    int cn4_visible[num_gal];
    double radius = 0.36;                       \\used for match
    double poly_racen[num_poly], poly_deccen[num_poly];
    double poly_maglim[num_poly];
    double fgot[num_poly];
    double galra[num_gal],galdec[num_gal];
    double completeness[num_gal],apm_limit[num_gal];
    double gal_apm[num_gal];
    double gal_match_poly_dis[num_gal][match_max];
    galv2 *gal_array;
    poly *poly_array;
    gal_array = (galv2*)malloc(sizeof(galv2)*num_gal);
    poly_array = (poly*)malloc(sizeof(poly)*num_poly);

    cout << "finish definition!" << endl;
\\ initialize all the variable
    for (long i = 0; i < num_gal; i++)
    {
        Npoly_match[i] = 0;
        polyid_for_gal[i] = 0;
        sectorid_for_gal[i] = 0;
        Npoly_match[i] = 0;
        cn4_visible[i] = 0;
        completeness[i] = 0;
        apm_limit[i] = 0;
        for (long j = 0; j < match_max)
        {
            gal_match_poly[i][j] = 0;
            gal_match_poly_dis[i][j] = 0;
        }
    }
    for (long i = 0; i < num_poly; i++)
    {
        linklist[i] = 0;
    }
    for (long i = 0; i < nc1; i++)
    {
        for (long j = 0; i < nc2; j++)
        {
            head_of_chain = 0;
        }
    }

    initgal_match_poly(num_gal,galra,galdec,num_poly,poly_racen,polydeccen,radius,match_max,Npoly_match,gal_match_poly,gal_match_poly_dist,head_of_chain,nc1,nc2,linklist);
    cout << "finish initalization of gal_match_poly" << endl;
\\ read in the data of gal sample\
    ifstream galfile;
    string galfilestr = "/home/gongjingyu/gcode/astro/2pcf/MY_RANDOM/myrandomSelect";
    galfile.open(galfilestr.data());
    for (long i; i < num_gal; i++)
    {
        string s;
        getline(galfile,s);
        istringstream is(s);
        string str1, str2, str3, str4, str5;
        const char *str1c, *str2c, *str3c, *str4c, *str5c;
        is >> str1 >> str2 >> str3 >> str4 >> str5;
        str2c = str2.c_str();
        str3c = str3.c_str();
        str4c = str4.c_str();
        str5c = str5.c_str();
        gal_array[i].xyz[0] = atof(str2c);
        gal_array[i].xyz[1] = atof(str3c);
        gal_array[i].xyz[2] = atof(str4c);
        gal_apm[i] = atof(str5c);
    }
\\ read in the data of poly
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
        num_poly = atof(str1c);
    }
    for long (i = 0; i < polygons; i++)
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
        for (j = 0; j < poly_array[i].cap_number; j++)
        {
            string st;
            getline(polyfile,st);
            istringstream ist(st);
            string str1t, str2t, str3t, str4t;
            const char *str1ct, *str2ct, *str3ct, *str4ct;
            is >> str1t >> str2t >> str3t >> str4t;
            str1ct = str1t.c_str();
            str2ct = str2t.c_str();
            str3ct = str3t.c_str();
            str4ct = str4t.c_str();
            poly_array[i].caps_xyz[j][0] = atof(str1ct);
            poly_array[i].caps_xyz[j][1] = atof(str2ct);
            poly_array[i].caps_xyz[j][2] = atof(str3ct);
            poly_array[i].cm[j] = atof(str4ct);
        }
    }
\\ read in the ra dec apm and fgot
    for (long i = 0; i < num_poly; i++)
    {
        string s;
        getline(polyfile,s);
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
\\ find the exact polygon ID for every galaxy
    for (long i = 0; i < num_gal; i++)
    {
        if (Npoly_match[i] = 0)
        {
            continue;
        }
        long tempid = match_one_poly();
        if (tempid > 0)
        {
            polyid_for_gal[i] = tempid;
            sectorid_for_gal[i] = isec[tempid];
            apm_limit[i] = poly_maglim[tempid];
            completeness[i] = fgot[tempid];
        }
    }
\\ select whether one galaxy sample can be observed
    long selnum = 0
    ofstream galpolymatchfile;
    string galpolymatchfilestr = "/home/gongjingyu/gcode/astro/MY_RANDOM/galpolymatch";
    galpolymatchfile.open(galpolymatchfilestr.data());
    for (long i = 0; i < num_gal; i++)
    {
        if (polyid_for_gal[i] <= 0)
        {
            cn1 += 1;
            continue;
        }
        if (gal_apm[i] > apm_limit[i])
        {
            cn2 += 1;
            continue;
        }
        if (completeness[i] < 7)
        {
            cn3 += 1;
            continue;
        }
        if ( > completeness[i])
        {
            cn4 += 4;
            cn4_visible[i] = -1;
        }
\\output the xyz or radecr, cn4visible, apparent magnitude, polygon ID, sector ID, completeness
        selnum += 1;
        galpolymatchfile << selnum << "  " << gal_array[i].xyz[0] << "  " << gal_array[i].xyz[1] << "  " << gal_array[i].xyz[2] << "  " << gal_apm[i] << "  " << cn4_visible[i] << "  " << polyid_for_gal[i] << "  " << sectorid_for_gal[i] << "  " << completeness[i] << endl;
    }
    ofstream readmefile;
    string readmefilestr = "/home/gongjingyu/gcode/astro/2pcf/MY_RANDOM/readme_galpolymatch";
    readmefile.open(readmefilestr.data());
    readmefile << "index:" << endl << "x:" << endl << "y:" << endl << "z:" << endl << "apparent magnitude:" << endl << "cn4_visible:" << endl << "polygonID:" << endl << "sectorID:" << endl << "completeness:" << endl <<"total number:" << selnum;
    timeend = time(NULL);
}
