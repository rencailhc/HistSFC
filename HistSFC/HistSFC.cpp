// HistSFC.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include "ahn2test.h"
#include "LoDgen.h"
#include "HistBuild.h"
#include "PyramidBuild.h"
#include "Simulation.h"

int main()
{
	//ahn2test();
	//LoDgen("E:/91.txt", "E:/91_lod.csv");
	//CoordTrans<double> trans({ 19000,369000,0 }, { 100,100,100 });
	//HistML("E:/91_lod.csv", trans);
	/*
	for (int i = 2; i <= 6; i++)
	{
		ExTree(HistIOT("sim2_ind_"+to_string(i)+"d_iot_4", 100, i), "hist_sim2_ind_" + to_string(i) + "d_4");
	}*/
	//ExTree(HistML("E:/91_lod.csv", trans), "HistTree_91_xyl_100_3");
	//ExTree(HistIOT("sim2_ind_2d_iot_3",100,2), "hist_sim2_ind_2d_3");
	//IOTbuild();
	//simgen();
	realsim();
	//idealsim_skew();
}

