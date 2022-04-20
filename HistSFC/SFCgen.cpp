/***Create text file with SFC keys***/

#include "CoordTransform.h"
#include "SFCConversion.h"
#include <math.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <string>
#include <numeric>


int main()
{
        ifstream filename;
        string x;

    filename.open("/nas2/usrdata/haicheng/sigmod/files.txt");
    if (!filename) {
        cout << "Unable to open file";
        exit(1); // terminate with error
    }

    //while (filename >> x)
        //while (1)
        {
        //FILE *input_file = fopen(("/nas2/usrdata/haicheng/sigmod/" + x).c_str(), "r");
        FILE *input_file = fopen("/nas2/usrdata/haicheng/sigmod/test.csv", "r");
        char buf[1024];
        char * pch, *lastpos;
        char ele[64];

        int j;
        int dimnum = 4;

        CoordTransform<double, int> cotrans(dimnum);
        double delta[4] = { 13427.64, 359007.3, -8.79, 12000 };
        double scale[4] = { 100,100,1000,-100 };
        cotrans.SetTransform(delta, scale);
        SFCConversion sfc(dimnum, 30);
        SFCConversion sfc2(dimnum-1, 30);

        Point<double> inPt(dimnum);
        Point<int> pt(dimnum);
        Point<int> pt2(dimnum-1);
        //ofstream output_file(("/nas2/usrdata/haicheng/sigmod/sfc_" + x).c_str());
          
        ofstream output_file("/nas2/usrdata/haicheng/ahn3_1_sfc.csv");
        while (1) //always true
        {
                memset(buf, 0, 1024);
                fgets(buf, 1024, input_file);

                if (strlen(buf) == 0) break; // no more data

                j = 0;
                lastpos = buf;
                pch = strchr(buf, ',');
                while (pch != 0)
                {
                        memset(ele, 0, 64);
                        strncpy(ele, lastpos, pch - lastpos);
                        if (strlen(ele) != 0)
                        {
                                inPt[j] = atof(ele);
                                j++;
                        }

                        lastpos = pch + 1;
                        pch = strchr(lastpos, ',');
                }

                if (strlen(lastpos) != 0 && strcmp(lastpos, "\n") != 0)//final part
                {
                        inPt[j] = atof(lastpos);
                        j++;
                }
                pt = cotrans.Transform(inPt);
                pt2[0] = pt[0];
                pt2[1] = pt[1];
                pt2[2] = pt[3];

                sfc_bigint Okey3D = sfc2.MortonEncode(pt2);
                sfc_bigint Okey4D = sfc.MortonEncode(pt);

                cout << Okey4D << ", " << Okey3D << ", " << fixed << setprecision(2)<< inPt[2] << "\n";
        }
        }
}

