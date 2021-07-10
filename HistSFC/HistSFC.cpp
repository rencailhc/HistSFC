// HistSFC.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include "ahn2test.h"
#include "LoDgen.h"
#include "HistBuild.h"
#include "PyramidBuild.h"
#include "Simulation.h"
#include "testHS.h"
#include "omp.h"


int main()
{
	/*
	auto tri = simplex(2,1);
	for (int i = 0; i < tri.getSize(); i++)
	{
		tri.faces[i].print_eq();
	}
	*/
	//disttest();
	//ahn2test();
	//LoDgen("E:/91.txt", "E:/91_lod.csv");
	//CoordTrans<double> trans({ 19000,369000,0 }, { 100,100,100 });
	//HistML("E:/91_lod.csv", trans);
	/*
	for (int i = 6; i <= 6; i++)
	{
		ExTree(HistIOT("simccv_real_"+to_string(i)+"d_5_iot", 1000, i), "hist_simccv_real_" + to_string(i) + "d_5");
	}*/
	//ExTree(HistML("E:/91_lod.csv", trans), "HistTree_91_xyl_100_3");
	//for(int i=2;i<=6;i++)
		//ExTree(HistHeight("simccv_ideal_6d_6_iot",8,6), "hist_simccv_ideal_6d_6");
	//IOTbuild();
	//idealsim_ccv();
	//realsim_ccv();
	//simgen();
	//realsim();
	//idealsim_skew();
	//test_uni();
	//test_uni_batch();
	testView();
	//viewgen();
	//idealsim_uni();
	/*
	int simnum = 100;
	for 2-nD query test
	unsigned seed = chrono::system_clock::now().time_since_epoch().count();
	mt19937 gen(seed);
	uniform_real_distribution<> dis(0.01*(1 << maxBits_sim2), 0.05*(1 << maxBits_sim2));	//for the 2-nD query
	double delta = dis(gen);

	RandWindow<long long> windowSim(2, 0, 1 << maxBits_sim2, delta, "sim2_ind_data");
	auto windowList = windowSim.Gen(simnum);
	int mDims = 2;
	PointCloudDB<long long, sfc_bigint> PCDB_ind("sim2_ind_3D_iot", 3);
	//PCDB_ind.HistTab = "hist_sim2_ind_3D";
	auto PCDBhist_ind = PCDB_ind;
	PCDBhist_ind.HIST = true;
	PCDBhist_ind.HistTab = "hist_sim2_ind_3D";

	Query<long long, sfc_bigint> testPlain_ind(PCDB_ind);
	Query<long long, sfc_bigint> testHist_ind(PCDBhist_ind);

	string fnSFC_ind = "E:/plainsfc_" + to_string(mDims) + "D_3D_cube_compare_ind.csv";
	string fnSFChist_ind = "E:/histsfc_" + to_string(mDims) + "D_3D_cube_compare_ind.csv";

	NDPoint<long long> P1(3);
	NDPoint<long long> P2(3);
	for (int i = 0; i < 3; i++)
	{
		P1[i] = 0;
		P2[i] = 1 << maxBits_sim2;
	}
	NDWindow<long long> windowFull(P1, P2);


	for (int i = 0; i < simnum; i++)
	{
		for (int j = 0; j < mDims; j++)
		{
			windowFull.minPoint[j] = windowList[i].minPoint[j];
			windowFull.maxPoint[j] = windowList[i].maxPoint[j];
		}

		testPlain_ind.QueryIOT(windowFull);
		testPlain_ind.ExMeasurement_batch(fnSFC_ind);
		testHist_ind.QueryIOT(windowFull);
		testHist_ind.ExMeasurement_batch(fnSFChist_ind);
	}
	*/

	
}

