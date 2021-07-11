/*Generating cLoD values,refer to http://nd-pc.org/documents/vario-nD-PC-v7.pdf*/

#pragma once
#include "SFCConversion.h"
#include "Point.h"
#include <math.h>
#include <random>
#include <iostream>
#include <fstream>
#include <numeric>
#include <chrono>
#include <iomanip>

using namespace std;

void LoDgen(string infile, string outfile)
{
	FILE *input_file;
	input_file = fopen(infile.c_str(), "r");
	if (!input_file)
	{
		throw "No input!";
	}
	else
	{
		NDPoint<double> inPt;
		char buf[1024];
		char * pch, *lastpos;
		char ele[64];
		int j;
		ofstream output_file(outfile);
		unsigned seed = chrono::system_clock::now().time_since_epoch().count();
		mt19937_64 generator(seed);
		uniform_real_distribution<double> dist(0.0, 1.0);
		auto real_rand = bind(dist, generator);
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

			/*setting of LoD range, depending on span of other dimensions*/
			double lod = 900;
			while (lod > 899.995)
			{
				double U = real_rand();
				lod = log2(262143 * U + 1) / 2 * 100;
			}

			output_file << fixed << setprecision(2) << inPt[0] << "," << inPt[1] << "," << inPt[2] << "," << lod << '\n';
		}
	}

	fclose(input_file);
	cout << "Time cost: " << clock() / CLOCKS_PER_SEC << "s" << endl;
}
