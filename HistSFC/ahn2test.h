#pragma once
#include <initializer_list>
#include <array>
#include <iostream>
#include <chrono>
#include <random>
#include <math.h>
#include "Window.h"

using namespace std;

template <typename T>
class WindowSim {
public:
	int nDims;
	T* DimLow;  //assuming the order is X, Y, Z/LoD, Time...
	T* DimHigh;

public:

	WindowSim()
	{
		nDims = 0;
		DimLow = nullptr;
		DimHigh = nullptr;
	}

	WindowSim(int len) : nDims(len), DimLow(new T[nDims]), DimHigh(new T[nDims])
	{

	}

	WindowSim(initializer_list<T> l1, initializer_list<T> l2)
	{
		nDims = (int)l1.size();
		DimLow = new T[nDims];
		DimHigh = new T[nDims];
		uninitialized_copy(l1.begin(), l1.end(), DimLow);
		uninitialized_copy(l2.begin(), l2.end(), DimHigh);
	}

	~WindowSim()
	{
		delete[] DimLow;
		delete[] DimHigh;
		nDims = 0;
	}

	vector<NDWindow<T>> WindowGen(int num)  //generate certain number of query windows
	{
		//used for computing edge length
			//array<int, 12> level = { 0,1,2,3,4,5,6,7,8,9,10,11 };
		array<double, 12> prob = { 0.001908156,0.056718847,0.112310723,0.17019676,0.23482707,0.312265428,0.400281378,0.576441766,0.713853012,0.84425464,0.94305549,1 };
		double baseunit = 125.0;  //the edge length at level 11, in meter
		unsigned seed = chrono::system_clock::now().time_since_epoch().count();
		mt19937 gen(seed);
		uniform_real_distribution<> disedge(0, 1);
		uniform_real_distribution<> disX(DimLow[0], DimHigh[0]);
		uniform_real_distribution<> disY(DimLow[1], DimHigh[1]);
		double randedge = 0;   //edge length of the XY box
		double randX = 0;
		double randY = 0;
		double randLoD = 0;
		double randZ = 0;
		double randprob = 0; //used for generating other dimensions
		vector<NDWindow<T>> windowList;

		for (int i = 0; i < num; i++)
		{
			randprob = disedge(gen);
			for (int j = 0; j < prob.size(); j++)
			{
				if (randprob <= prob[j] and randprob >= prob[0])
				{
					double randlevel = j - (prob[j] - randprob) / (prob[j] - prob[j - 1]);
					randedge = baseunit * pow(2, prob.size() - 1 - randlevel);   //unit is meter
					break;
				}
			}
			randX = disX(gen);
			randY = disY(gen);

			switch (nDims) {
			case 2:	//XY domain, according to Table 3 in (van Oosterom et al., 2017) 
			{
				NDWindow<T> windowRes(NDPoint<T>({ randX, randY }), NDPoint<T>({ randX + randedge, randY + randedge }));
				windowList.push_back(windowRes);
				break;
			}

			case 3:	//XYLoD domain according to logs, or synthetic XYZ domain
			{
				/*XYLoD box*/
				//LoD distribution from log
				array<double, 14> LoDprob = { 0.007455621,0.029167063,0.080426088,0.18326968,0.312803326,0.376729529,0.434892926,0.494352701,0.558317004,0.627281623,0.688470867,0.75517647,0.839536063,1 };
				uniform_real_distribution<> disLoD(0, 1);
				randprob = disLoD(gen);
				for (int j = 0; j < LoDprob.size(); j++)
				{
					if (randprob <= LoDprob[j] and randprob >= LoDprob[0])
					{
						randLoD = j - (LoDprob[j] - randprob) / (LoDprob[j] - LoDprob[j - 1]);
						randLoD = randLoD * DimHigh[2]/14;   // scale up
						break;
					}
				}
				NDWindow<T> windowRes(NDPoint<T>({ randX, randY, 0 }), NDPoint<T>({ randX + randedge, randY + randedge, randLoD }));
				windowList.push_back(windowRes);
				break;

				/*XYZ box*/
				/*
				array<double, 3> Z = { 0,30,100 };
				array<double, 3> Zprob = { 0.5,0.8,1 };  //corresponding to 0, 30, 100 of elevation
				uniform_real_distribution<> disZ(0, 1);
				randprob = disZ(gen);
				for (int j = 0; j < Zprob.size(); j++)
				{
					if (randprob <= Zprob[j])
					{
						randZ = Z[j];
						break;
					}
				}
				NDWindow<T> windowRes(NDPoint<T>({ randX, randY, randZ }), NDPoint<T>({ randX + randedge, randY + randedge,999 }));
				windowList.push_back(windowRes);
				break;*/
			}

			case 4:	//XYZLoD domain, synthetic
			{
				array<double, 14> LoDprob = { 0.007455621,0.029167063,0.080426088,0.18326968,0.312803326,0.376729529,0.434892926,0.494352701,0.558317004,0.627281623,0.688470867,0.75517647,0.839536063,1 };
				uniform_real_distribution<> disLoD(0, 1);
				randprob = disLoD(gen);
				for (int j = 0; j < LoDprob.size(); j++)
				{
					if (randprob <= LoDprob[j] and randprob >= LoDprob[0])
					{
						randLoD = j - (LoDprob[j] - randprob) / (LoDprob[j] - LoDprob[j - 1]);
						break;
					}
				}

				array<double, 3> Z = { -3, 30, 65 };
				array<double, 3> Zprob = { 0.2, 0.8, 1 };  //corresponding to 0, 30, 100 of elevation
				uniform_real_distribution<> disZ(0, 1);
				randprob = disZ(gen);
				for (int j = 0; j < Zprob.size(); j++)
				{
					if (randprob < Zprob[j])
					{
						randZ = Z[j] * disZ(gen);
						break;
					}
				}
				NDWindow<T> windowRes(NDPoint<T>({ randX, randY, randZ, 0 }), NDPoint<T>({ randX + randedge, randY + randedge, DimHigh[2],randLoD }));
				windowList.push_back(windowRes);
				break;
			}
			default:
				cout << "currently not implemented!" << endl;
			}
		}

		return windowList;
	}

};


void ahn2test();