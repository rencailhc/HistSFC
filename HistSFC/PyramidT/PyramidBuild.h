/*
Building the Pyramid-Technique solution using Oracle IOT. For more details, please refer to 
Berchtold, Stefan et al., "The pyramid-technique: Towards breaking the curse of dimensionality." SIGMOD 1998.
Shi, Qingxiu and Bradford Nickerson. "Decreasing radius k-nearest neighbor search using mapping-based indexing schemes." University of New Brunswick, 2006.
*/


#pragma once
#include <iostream>
#include <fstream>
#include <iomanip>
#include "Window.h"

using namespace std;

/*points uniformly distributed*/
template <typename T>
void Pvalue_uni(string InFile, string OutFile, const CoordTrans& trans, T* dimmax)
{
	char buf[1024];
	short dimnum = trans.dimnum;
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
	NDPoint<T> Pt_O;
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
				inPt[j] = (atof(ele) - trans._delta[j]) / (dimmax[j] - trans._delta[j]);
				j++;
			}

			lastpos = pch + 1;
			pch = strchr(lastpos, ',');
		}

		if (strlen(lastpos) != 0 && strcmp(lastpos, "\n") != 0)//final part
		{
			Pt_O[j] = atof(lastpos);
			inPt[j] = (atof(lastpos) - trans._delta[j]) / (dimmax[j] - trans._delta[j]);
			j++;
		}
		/*
		//setting organizing dimensions
		double dim2 = inPt[2];
		inPt[2] = inPt[3];
		short odims = dimnum - 1;
		
		double dvalue = 0;
		double d_max = 0;
		double h_v = abs(0.5 - inPt[0]);
		for (int k = 1; k < odims; k++)
		{
			if (h_v < abs(0.5 - inPt[k]))
			{
				d_max = k;
				h_v = abs(0.5 - inPt[k]);
			}
		}
		if (inPt[d_max] < 0.5) dvalue = d_max;
		else dvalue = odims + d_max;

		double p_v = dvalue + h_v;
		//cout << dvalue << "," <<fixed << setprecision(10) << h_v << endl;
		output_file << fixed << setprecision(8) << p_v << ", ";	//the precision also depends on total number of points	

		for (int i = 0; i < dimnum - 1; i++)
			output_file << setprecision(2) << Pt_O[i] << ",";

		output_file << Pt_O[dimnum - 1] << endl;
		*/

		for (int i = 2; i <= dimnum; i ++)
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
			output_file << setprecision(10) << p_v << ", ";
		}

		for (int i = 0; i < dimnum - 1; i++)
			output_file << Pt_O[i] << ",";

		output_file << Pt_O[dimnum - 1] << endl;

	}

}


/*points unevenly distributed*/
template <typename T>
void Pvalue_shift(string InFile, string OutFile, const CoordTrans& trans, T* dimmax, T* medians)
{
	//dimmax and medians can be acquired in the metadata table
	char buf[1024];
	short dimnum = trans.dimnum;
	ofstream output_file(OutFile);
	ifstream data(InFile);
	if (!data.is_open())
	{
		cout << "Error opening sfc file";
		exit(1);
	}

	char * pch, *lastpos;
	char ele[64];

	int j;
	NDPoint<T> Pt_O;
	NDPoint<double> inPt;

	T* datamedian = medians;

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
				inPt[j] = pow((atof(ele) - trans._delta[j]) / (dimmax[j]- trans._delta[j]), -1 / (log2(datamedian[j] - trans._delta[j]) - log2(dimmax[j] - trans._delta[j])));
				j++;
			}

			lastpos = pch + 1;
			pch = strchr(lastpos, ',');
		}

		if (strlen(lastpos) != 0 && strcmp(lastpos, "\n") != 0)//final part
		{
			Pt_O[j] = atof(lastpos);
			inPt[j] = pow((atof(lastpos) - trans._delta[j]) / (dimmax[j] - trans._delta[j]), -1 / (log2(datamedian[j] - trans._delta[j]) - log2(dimmax[j] - trans._delta[j])));
			j++;
		}

		/*
		//setting organizing dimensions
		double dim2 = inPt[2];
		inPt[2] = inPt[3];  
		short odims = dimnum - 1;

		double dvalue = 0;
		double d_max = 0;
		double h_v = abs(0.5 - inPt[0]);
		for (int k = 1; k < odims; k++)
		{
			if (h_v < abs(0.5 - inPt[k]))
			{
				d_max = k;
				h_v = abs(0.5 - inPt[k]);
			}
		}
		if (inPt[d_max] < 0.5) dvalue = d_max;
		else dvalue = odims + d_max;

		double p_v = dvalue + h_v;
		//cout << dvalue << "," <<fixed << setprecision(10) << h_v << endl;
		output_file << fixed << setprecision(PREC) << p_v << ", " << setprecision(2);	//the precision also depends on total number of points

		for (int i = 0; i < dimnum - 1; i++)
			output_file << Pt_O[i] << ",";

		output_file << Pt_O[dimnum - 1] << endl;
		*/
		for (int i = 2; i <= dimnum; i ++)
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
			output_file << setprecision(10) << p_v << ", ";
		}

		for (int i = 0; i < dimnum - 1; i++)
			output_file << Pt_O[i] << ",";

		output_file << Pt_O[dimnum - 1] << endl;

	}

}


void IOTbuild();
