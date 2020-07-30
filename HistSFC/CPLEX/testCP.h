#pragma once
#include "QueryCP.h"
#include "../PointCloud.h"
#include "../Simulation.h"
#include "../Geom.h"

void test_uni()
{
	short nDims = 4;
	PointCloudDB<long long, sfc_bigint> PCDB("simuni_" + to_string(nDims) + "D_iot", nDims);
	QueryExtend<long long, sfc_bigint> testPlain(PCDB);

	string fnSFC = "E:/nDalgorithm/uni/plainsfc_" + to_string(nDims) + "D_geom vs window.csv";

	/*
	halfspace h1({ 1,1,1,1 }, 4000, Sign::le);
	halfspace h2({ -1,1,0,0 }, 0, Sign::le);
	halfspace h3({ 0,1,0,1 }, 400, Sign::ge);
	*/
	halfspace h1({ 1,0,0,0 }, 410, Sign::ge);
	halfspace h2({ 1,0,0,0 }, 1410, Sign::le);
	halfspace h3({ 0,1,0,0 }, 720, Sign::ge);
	halfspace h4({ 0,1,0,0 }, 1720, Sign::le);
	halfspace h5({ 0,0,1,0 }, 410, Sign::ge);
	halfspace h6({ 0,0,1,0 }, 1410, Sign::le);
	halfspace h7({ 0,0,0,1 }, 720, Sign::ge);
	halfspace h8({ 0,0,0,1 }, 1720, Sign::le);

	NDGeom geom(nDims);
	geom.add(h1);
	geom.add(h2);
	geom.add(h3);
	geom.add(h4);
	geom.add(h5);
	geom.add(h6);
	geom.add(h7);
	geom.add(h8);
	testPlain.QueryIOT(geom);
	testPlain.ExMeasurement(fnSFC);

	Query<long long, sfc_bigint> testPlain2(PCDB);
	NDWindow<long long> wind({ 410, 720, 410, 720 }, { 1410, 1720, 1410, 1720 });
	testPlain2.QueryIOT(wind);
	testPlain2.ExMeasurement(fnSFC);
}
