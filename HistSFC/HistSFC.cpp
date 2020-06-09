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
	//ExTree(HistML("E:/91_lod.csv", trans), "HistTree_91_xyl_100_3");
	//ExTree(HistIOT("ahn2iot_91_xyl_plainsfc"), "HistTree_91_xyl_100");
	//IOTbuild();
	//simgen();
	simtest();
	
}

