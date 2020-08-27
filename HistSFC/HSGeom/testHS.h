#pragma once
#include "QueryCP.h"
#include "QuerySweep.h"
#include "../PointCloud.h"
#include "../Simulation.h"
#include "../Geom.h"

void test_uni()
{
	short nDims = 2;
	PointCloudDB<long long, sfc_bigint> PCDB("simuni_" + to_string(nDims) + "D_iot", nDims);
	QueryExtendCP<long long, sfc_bigint> testPlainCP(PCDB);
	QueryExtendSweep<long long, sfc_bigint> testPlainSweep(PCDB);

	string fnSFC = "E:/nDalgorithm/uni/HSsfc_" + to_string(nDims) + "D_geom.csv";

	halfspace h1({ 1,1 }, 2000, Sign::le);
	halfspace h2({ -1,1 }, 0, Sign::le);
	halfspace h3({ 0,1 }, 600, Sign::ge);
	
	/*
	halfspace h1({ 1,0,0,0 }, 410, Sign::ge);
	halfspace h2({ 1,0,0,0 }, 1410, Sign::le);
	halfspace h3({ 0,1,0,0 }, 720, Sign::ge);
	halfspace h4({ 0,1,0,0 }, 1720, Sign::le);
	halfspace h5({ 0,0,1,0 }, 410, Sign::ge);
	halfspace h6({ 0,0,1,0 }, 1410, Sign::le);
	halfspace h7({ 0,0,0,1 }, 720, Sign::ge);
	halfspace h8({ 0,0,0,1 }, 1720, Sign::le);
	*/
	NDGeom geom(nDims);
	geom.add(h1);
	geom.add(h2);
	geom.add(h3);
	//geom.add(h4);
	//geom.add(h5);
	//geom.add(h6);
	//geom.add(h7);
	//geom.add(h8);
	testPlainCP.QueryIOT(geom);
	testPlainCP.ExMeasurement(fnSFC);

	testPlainSweep.QueryIOT(geom);
	testPlainSweep.ExMeasurement(fnSFC);

	/*
	Query<long long, sfc_bigint> testPlain2(PCDB);
	NDWindow<long long> wind({ 410, 720, 410, 720 }, { 1410, 1720, 1410, 1720 });
	testPlain2.QueryIOT(wind);
	testPlain2.ExMeasurement(fnSFC);
	*/
}

void test_uni_batch()
{
	/*dimension effect*/
	/*for (int ndim = 2; ndim <= 16; ndim += 2)
	{
		PointCloudDB<long long, sfc_bigint> PCDB("simuni_" + to_string(ndim) + "D_iot", ndim);
		QueryExtendCP<long long, sfc_bigint> testPlainCP(PCDB);
		QueryExtendSweep<long long, sfc_bigint> testPlainSweep(PCDB);

		string fnSFC = "E:/nDalgorithm/uni/HSsfc_dimension_effect.csv";

		halfspace h1(ndim);
		halfspace h2(ndim);
		halfspace h3(ndim);
		for (int i = 0; i < ndim; i++)
		{
			h1.w[i] = 1;
			h2.w[i] = 0;
			h3.w[i] = 0;
		}
		h1.b = 1000 * ndim;
		h1.s = Sign::le;

		h2.w[0] = -1;
		h2.w[1] = 1;
		h2.b = 0;
		h2.s = Sign::le;

		h3.w[ndim-1] = 1;
		h3.b = 600-35*ndim;
		h3.s = Sign::ge;

		NDGeom geom(ndim);
		geom.add(h1);
		geom.add(h2);
		geom.add(h3);

		testPlainCP.QueryIOT(geom);
		testPlainCP.ExMeasurement(fnSFC);

		testPlainSweep.QueryIOT(geom);
		testPlainSweep.ExMeasurement(fnSFC);
	}*/


	/*geom complexity effect*/
	short nDims = 6;
	PointCloudDB<long long, sfc_bigint> PCDB("simuni_" + to_string(nDims) + "D_iot", nDims);
	QueryExtendCP<long long, sfc_bigint> testPlainCP(PCDB);
	QueryExtendSweep<long long, sfc_bigint> testPlainSweep(PCDB);

	string fnSFC = "E:/nDalgorithm/uni/HSsfc_" + to_string(nDims) + "D_geom_complexity.csv";

	halfspace h2(nDims);
	for (int k = 0; k < nDims; k++)
		h2.w[k] = 0;
	h2.w[0] = -1;
	h2.w[1] = 1;
	h2.b = 0;
	h2.s = Sign::le;

	NDGeom geom(nDims);
	geom.add(h2);

	int face_num = 1;
	for (int i = 0; i < nDims; i++)
	{
		for (int j = i+1; j < 4; j++)
		{
			face_num++;
			halfspace h1(nDims);
			for (int k = 0; k < nDims; k++)
				h1.w[k] = 0;
			h1.w[i] = 1;
			h1.w[j] = 1;
			//h1.b = face_num*400;
			h1.s = Sign::le;
			geom.add(h1);
			for (int k = 1; k < face_num; k++)
				geom.faces[k].b = log2(face_num) * 500;

			ofstream output(fnSFC, ios::app | ios::out);
			output << "Number of geom faces: " << geom.faces.size() << endl;

			testPlainCP.QueryIOT(geom);
			testPlainCP.ExMeasurement(fnSFC);

			testPlainSweep.QueryIOT(geom);
			testPlainSweep.ExMeasurement(fnSFC);

		}

	}
	
}
