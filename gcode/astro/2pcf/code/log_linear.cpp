#include<iostream>
#include<fstream>
#include<cstring>
#include<sstream>
#include<ctime>
#include<vector>
#include "tools.h"

int main()
{
    ifstream tpcflogfile;
    ofstream tpcflinearfile;
    string tpcflogfilename = "/home/gongjingyu/gcode/astro/2pcf/outcome/20180227/20180227_tpcfv2";
    string tpcflinearfilename = "/home/gongjingyu/gcode/astro/2pcf/outcome/20180227/20180227_tpcfv2tolinear";
    tpcflogfile.open(tpcflogfilename.data());
    tpcflinearfile.open(tpcflinearfilename.data());
    long rpn = 40, rpin = 40;
    double **tpcflog, **tpcflinear;
    tpcflog = new double *[rpin*2];
    tpcflinear = new double *[rpin*2];
    for (int i = 0; i < rpin*2; i++)
    {
        tpcflog[i] = new double [rpn*2];
        tpcflinear[i] = new double [rpn*2];
    }
    for (int i = 0; i < rpin*2; i++)
    {
        string s;
        getline(tpcflogfile,s);
        istringstream is(s);
        for (int j = 0; j < rpn*2; j++)
        {
            string str;
            const char *strc;
            is >> str;
            strc = str.c_str();
            tpcflog[i][j] = atof(strc);
        }
    }
    double tpcflogtmpvalue[rpn];
    double tpcflogtmpcoordinate[rpn];
    double factor;
    factor = pow(200.,0.025);
    for (int i = rpn-1; i > -1; i--)
    {
        if (i == rpn-1)
        {
            tpcflogtmpcoordinate[i] = 40;
        }
        else
        {
            tpcflogtmpcoordinate[i] = tpcflogtmpcoordinate[i+1]/factor;
        }
    }
    for (int i = 0; i < rpin*2; i++)
    {
        for (int j = 0; j < rpn; j++)
        {
            tpcflogtmpvalue[j] = tpcflog[i][j+rpn];
        }
        for (int j = 0; j < rpn; j++)
        {
            tpcflinear[i][j+rpn] = linearinterpolate(tpcflogtmpcoordinate,tpcflogtmpvalue,rpn,j+1);
        }
    }
    filltpcf(tpcflinear,rpin,rpn);
    for (int i = 0; i < rpin*2; i++)
    {
        for (int j = 0; j < rpn*2; j++)
        {
            tpcflinearfile << tpcflinear[i][j];
            if (j != (rpn*2-1))
            {
                tpcflinearfile << " ";
            }
        }
        if (i != (rpin*2-1))
        {
            tpcflinearfile << endl;
        }
    }
    tpcflogfile.close();
    tpcflinearfile.close();
    for (int i = 0; i < rpin*2; i++)
    {
        delete [] tpcflog[i];
        delete [] tpcflinear[i];
    }
    delete [] tpcflog;
    delete [] tpcflinear;
    return 0;
}
