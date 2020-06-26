// PyramidT.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <string>
#include "PyramidBuild.h"
#include "PyramidQuery.h"

int main()
{
	Pvalue("E:/sigmod/simres/pt/sim1_data.csv", "E:/sigmod/simres/pt/pt_sim1_uni_data.csv",16);
	/*
	///////////
	int datamedian[] = { 1757,2186,2606,2253,2074,2106,1383,1753,1436, 2910, 2452, 1217, 2685,2046, 2180, 3089 };
	///////////

	for (int Dims = 4; Dims <= 16; Dims += 2)
	{
		string output = "E:/sigmod/simres/pt/" + to_string(Dims) + "D" + to_string(Dims) + "D_cube_compare_sim1.csv";
		ofstream output_file(output);
		long long trueres = 0;
		long long plainres = 0;
		long long histres = 0;

		//cout << "------------" << Dims << "D--------------" << endl;
		//string filename = "E:/sigmod/simres/" + to_string(i) + "Dwindow.txt";
		string filename = "E:/sigmod/simres/sim1/" + to_string(Dims) + "D" + to_string(Dims) + "D_hyperwindow.txt";
		string IOTTab = "pv_sim1_" + to_string(Dims) + "D_iot";
		FILE *input_file = fopen(filename.c_str(), "r");
		char buf[1024];
		char * pch, *lastpos;
		char ele[64];
		int j;

		int i = Dims;

		NDPoint<double> QL(Dims);
		NDPoint<double> QH(Dims);
		for (int d = 0; d < Dims; d++)
		{
			QL[d] = 0;
			QH[d] = 4096 / (V_MAX);
		}

		int linenum = 1;

		while (linenum) //always true
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
				if (i == 1)
				{
					if (j == 0)
					{
						strncpy(ele, lastpos + 1, pch - lastpos - 2);
						//QL[j] = atof(ele)/(V_MAX);
						QL[j] = pow(atof(ele) / (V_MAX), -1 / (log2(datamedian[j]) - log2(V_MAX)));
						j++;
					}
					else
					{
						strncpy(ele, lastpos + 2, pch - lastpos - 3);
						//QH[j - i] = atof(ele) / (V_MAX);
						QH[j - i] = pow(atof(ele) / (V_MAX), -1 / (log2(datamedian[j - i]) - log2(V_MAX)));
						j++;
					}
				}
				else
				{
					if (j <= i - 1)
					{
						if (j == 0)
							strncpy(ele, lastpos + 1, pch - lastpos - 1);
						else if (j == i - 1)
							strncpy(ele, lastpos + 1, pch - lastpos - 2);
						else
							strncpy(ele, lastpos + 1, pch - lastpos - 1);
						if (strlen(ele) != 0)
						{
							//QL[j] = atof(ele) / (V_MAX);
							QL[j] = pow(atof(ele) / (V_MAX), -1 / (log2(datamedian[j]) - log2(V_MAX)));
							j++;
						}
					}
					else
					{
						if (j == i)
							strncpy(ele, lastpos + 2, pch - lastpos - 2);
						else if (j == 2 * i - 1)
							strncpy(ele, lastpos + 1, pch - lastpos - 2);
						else
							strncpy(ele, lastpos + 1, pch - lastpos - 1);
						if (strlen(ele) != 0)
						{
							//QH[j - i] = atof(ele) / (V_MAX);
							QH[j - i] = pow(atof(ele) / (V_MAX), -1 / (log2(datamedian[j - i]) - log2(V_MAX)));
							j++;
						}
					}
				}

				lastpos = pch + 1;
				pch = strchr(lastpos, ',');

			}

			if (strlen(lastpos) != 0 && strcmp(lastpos, "\n") != 0)//final part
				trueres = (long long)atof(lastpos);


			//cout << "########" << to_string(linenum) << "th window" << "########" << endl;
			NDWindow<double> qrec(QL, QH);
			//map<double, double> ranges = QueryRange<double>(qrec);

			//for (auto it = ranges.begin(); it != ranges.end(); it++)
				//cout << fixed<<setprecision(8) << it->first << ", " << it->second << endl;

			//cout << linenum<< ", another" << endl;

			//plainres = QueryIOT<double>(IOTTab, qrec);
			output_file << trueres << ", " << plainres << ", " << (plainres - trueres)*1.0f / trueres << endl;

			linenum++;
		}
	}*/
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
