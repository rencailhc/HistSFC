#include <iostream>
#include <ctime>
#include "ahn2test.h"
#include "Query.h"
#include "PyramidQuery.h"

void ahn2test()
{
	int simnum = 1;
	WindowSim<double> ndw( {19000,369000,0}, {20000, 370000,900} );
	vector<NDWindow<double>> windowList = ndw.WindowGen(simnum);
	CoordTrans<double> trans( {19000,369000,0},{100,100,100} );
	CoordTrans<double> transPyramid({ 19000,369000,0 }, { 1.0/1000,1.0/1000,1.0/900 });	//scale computed as 1/(dimmax-delta)
	PointCloudDB<double, sfc_bigint> PCDB("ahn2iot_91_xyl_plainsfc", 3, 28992, trans, true, "histtree_91_xyl_100");
	PyramidDB<double> PCDBPyramid("ahn2iot_91_xyl_pyramid", 3, 28992, transPyramid, true, { 19491.54, 369555.01, 849.98 });
	Query<double, sfc_bigint> testHist(PCDB);
	PCDB.HIST = false;
	Query<double, sfc_bigint> testPlain(PCDB);
	PyramidQuery<double> testPyramid(PCDBPyramid);

	time_t t = time(0);   // get time now
	struct tm * now = localtime(&t);
	char buffer[80];
	strftime(buffer, 80, "%Y-%m-%d %Hh%Mm", now);
	string timetag(buffer);
	string filename = "E:/AHN2TEST " + timetag +".txt";
	for (int i = 0; i < simnum; i++)
	{
		//testHist.QueryIOT(windowList[i]);
		//testHist.ExMeasurement(filename);
		//testPlain.QueryIOT(windowList[i]);
		//testPlain.ExMeasurement(filename);
		testPyramid.QueryIOT(windowList[i]);
		testPyramid.ExMeasurement(filename);
	}


}