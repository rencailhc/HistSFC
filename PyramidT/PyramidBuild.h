#pragma once
#include <iostream>
#include <fstream>
#include <iomanip>
#include "Window.h"

#define V_MAX 1<<12

using namespace std;

void Pvalue(string InFile, string OutFile, int tot_dim)
{
	char buf[1024];

	ofstream output_file(OutFile);
	ifstream data(InFile);
	if (!data.is_open())
	{
		cout << "Error opening sfc file";
		exit(1);
	}

	char * pch, *lastpos;
	char ele[64];

	int i, j;
	i = 0;
	NDPoint<int> Pt_O;
	NDPoint<double> inPt;

	while (!data.eof())
	{
		data.getline(buf, 1024);

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
				Pt_O[j] = atof(ele);
				inPt[j] = atof(ele) / (V_MAX);
				j++;
			}

			lastpos = pch + 1;
			pch = strchr(lastpos, ',');
		}

		if (strlen(lastpos) != 0 && strcmp(lastpos, "\n") != 0)//final part
		{
			Pt_O[j] = atof(lastpos);
			inPt[j] = atof(lastpos) / (V_MAX);
			j++;
		}

		for (int i = 2; i <= tot_dim; i += 2)
		{
			double dvalue = 0;
			double d_max = 0;
			double h_v = abs(0.5 - inPt[0]);
			for (int k = 1; k < i; k++)
			{
				if (h_v < abs(0.5 - inPt[k]))
				{
					d_max = k;
					h_v = abs(0.5 - inPt[k]);
				}
			}
			if (inPt[d_max] < 0.5) dvalue = d_max;
			else dvalue = i + d_max;

			double p_v = dvalue + h_v;
			//cout << dvalue << "," <<fixed << setprecision(10) << h_v << endl;
			output_file << fixed << setprecision(8) << p_v << ", ";
		}

		for (int i = 0; i < tot_dim - 1; i++)
			output_file << Pt_O[i] << ",";

		output_file << Pt_O[tot_dim - 1] << endl;

	}

}

void Pvalue_shift(string InFile, string OutFile, int* medians, int tot_dim)
{
	char buf[1024];

	ofstream output_file(OutFile);
	ifstream data(InFile);
	if (!data.is_open())
	{
		cout << "Error opening sfc file";
		exit(1);
	}

	char * pch, *lastpos;
	char ele[64];

	int i, j;
	i = 0;
	NDPoint<int> Pt_O;
	NDPoint<double> inPt;

	int* datamedian = medians;

	while (!data.eof())
	{
		data.getline(buf, 1024);

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
				Pt_O[j] = atof(ele);
				inPt[j] = pow(atof(ele) / (V_MAX), -1 / (log2(datamedian[j]) - log2(V_MAX)));
				j++;
			}

			lastpos = pch + 1;
			pch = strchr(lastpos, ',');
		}

		if (strlen(lastpos) != 0 && strcmp(lastpos, "\n") != 0)//final part
		{
			Pt_O[j] = atof(ele);
			inPt[j] = pow(atof(ele) / (V_MAX), -1 / (log2(datamedian[j]) - log2(V_MAX)));
			j++;
		}

		for (int i = 2; i <= tot_dim; i += 2)
		{
			double dvalue = 0;
			double d_max = 0;
			double h_v = abs(0.5 - inPt[0]);
			for (int k = 1; k < i; k++)
			{
				if (h_v < abs(0.5 - inPt[k]))
				{
					d_max = k;
					h_v = abs(0.5 - inPt[k]);
				}
			}
			if (inPt[d_max] < 0.5) dvalue = d_max;
			else dvalue = i + d_max;

			double p_v = dvalue + h_v;
			//cout << dvalue << "," <<fixed << setprecision(10) << h_v << endl;
			output_file << fixed << setprecision(8) << p_v << ", ";
		}

		for (int i = 0; i < tot_dim - 1; i++)
			output_file << Pt_O[i] << ",";

		output_file << Pt_O[tot_dim - 1] << endl;

	}

}