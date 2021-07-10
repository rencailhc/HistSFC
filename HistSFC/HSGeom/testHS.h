#pragma once
#define _USE_MATH_DEFINES
#include <math.h>
#include "QueryCP.h"
#include "QuerySweep.h"
#include "../PointCloud.h"
#include "../Simulation.h"
#include "../Geom.h"
#include "../ViewStruct.h"
#include "../ViewQuery.h"

#define D 1000  //largest distance in view


NDGeom simplex(int dimnum, int scale) {
	NDGeom simp(dimnum);
	for (int i = 0; i < dimnum; i++)
	{
		halfspace hs(dimnum);
		for (int j = 0; j <= dimnum; j++)
		{
			if (j == i)
				hs.w[j] = 1;
			else
				hs.w[j] = 0;
		}
		hs.b = 193;
		hs.s = Sign::ge;

		simp.add(hs);
	}

	halfspace hs_n(dimnum);
	for (int j = 0; j <= dimnum; j++)
	{
		hs_n.w[j] = 1;
	}
	hs_n.b = scale + 193* (dimnum-2);
	hs_n.s = Sign::le;

	simp.add(hs_n);

	return simp;
}

int factorial(int n)
{
	if (n > 1)
		return n * factorial(n - 1);
	else
		return 1;
}

/*regular simplex with starting vertex at the origin*/
NDGeom simplexR(int dimnum, int scale) {
	NDGeom simp(dimnum);
	double v10 = 0.0010091; //ratio of simplex volume at 10d, with residual ratio 0.1
	int n = dimnum;
	double r = (n*sqrt(2) / (sqrt(n + 1) + 1.0) - pow(v10*sqrt(pow(2, n) / (n + 1.0))*factorial(n), 1.0 / n))*sqrt(n + 1) / (n*sqrt(2));
	double m = (1 + (1 + sqrt(n + 1)) / n) / (n + 1);

	for (int i = 0; i < dimnum; i++)
	{
		halfspace hs(dimnum);
		for (int j = 0; j <= dimnum; j++)
		{
			if (j == i)
				hs.w[j] = 1-m;
			else
				hs.w[j] = -m;
		}
		hs.b = 0;
		hs.s = Sign::le;

		simp.add(hs);
	}

	halfspace hs_n(dimnum);
	for (int j = 0; j < dimnum; j++)
	{
		hs_n.w[j] = 1;
	}
	hs_n.b = scale * (1 / m - r * n);
	hs_n.s = Sign::le;

	simp.add(hs_n);

	return simp;
}

NDGeom octagonXD(int dimnum, int scale, int f) {
	//f defines number of half-spaces which is 2f
	NDGeom oct(dimnum);

	double selectivity = 0.001;

	for (int i = -f + 1; i <= f; i++)
	{
		double theta = M_PI * i / f;
		halfspace hs(dimnum);
		hs.w[0] = cos(theta);
		hs.w[1] = sin(theta);
		hs.b = sqrt(selectivity / M_PI)*scale;
		hs.b += (hs.w[0] + hs.w[1]) * scale / 2;	//move to the center of the data region
		hs.s = Sign::le;

		oct.add(hs);
	}

	//boudary confined by other dimensions are not included as they are embeded in the data region
	
	return oct;
}

void test_uni()
{
	short nDims = 6;
	PointCloudDB<long long, sfc_bigint> PCDB("simccv_ideal_" + to_string(nDims) + "D_1_iot", nDims);
	//PCDB.HIST = true;
	//PCDB.HistTab = "hist_simccv_ideal_" + to_string(nDims) + "d_1";
	QueryExtendCP<long long, sfc_bigint> testPlainCP(PCDB);
	QueryExtendSweep<long long, sfc_bigint> testPlainSweep(PCDB);
	QueryExtendSweep<long long, sfc_bigint> testPlainSphere(PCDB);
	QueryExtendSweep<long long, sfc_bigint> testPlainVertex(PCDB);

	string fnSFC = "E:/nDalgorithm/uni/HSsfc_" + to_string(nDims) + "D_geom.csv";
	/*
	halfspace h1({ 1,1,0,0,0,0 }, 2000, Sign::le);
	halfspace h2({ -1,1,0,0,0,0 }, 0, Sign::le);
	halfspace h3({ 0,1,0,0,0,0 }, 600, Sign::ge);
	
	
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
	testPlainCP.QueryIOT(geom);
	testPlainCP.ExMeasurement(fnSFC);
	*/

	//auto geom = simplexR(nDims,1<< maxBits_sim1);
	for (int i = 1; i <= 4; i++)
	{
		auto geom = octagonXD(nDims, 1 << maxBits_sim1, 8*i);
		testPlainSweep.QueryIOT(geom, ALG::SWEEP);
		testPlainSweep.ExMeasurement_batch(fnSFC);

		//testPlainSphere.QueryIOT(geom, ALG::SPHERE);
		//testPlainSphere.ExMeasurement_batch(fnSFC);

		//testPlainVertex.QueryIOT(geom, ALG::VERTEX);
		//testPlainVertex.ExMeasurement_batch(fnSFC);
	}
	

	/*
	testPlainSphere.QueryIOT(geom, ALG::SPHERE);
	testPlainSphere.ExMeasurement(fnSFC);

	testPlainVertex.QueryIOT(geom, ALG::VERTEX);
	testPlainVertex.ExMeasurement(fnSFC);

	
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

			testPlainSweep.QueryIOT(geom, ALG::SWEEP);
			testPlainSweep.ExMeasurement(fnSFC);

		}

	}
	
}


void testView()
{
	viewlib::viewPos testv = { {1000,1000,1000},{0.57735,0.57735,0.57735},D };
	auto vpoly = viewlib::view2poly(testv, 1<< maxBits_sim1);
	auto geom = viewlib::view2geom(testv, 1 << maxBits_sim1);
	PointCloudDB<int, sfc_bigint> PCDB("simccv_ideal_4D_1_iot", 4);
	QueryExtendSweep<int, sfc_bigint> testPlainSweep(PCDB);
	QueryExtendCP<int, sfc_bigint> testPlainCP(PCDB);
	ViewQuery<int, sfc_bigint> testview(PCDB);

	//testview.QueryIOT(testv, 1 << maxBits_sim1);
	testPlainSweep.QueryView(vpoly, geom, ALG::SWEEP);
	//testPlainCP.QueryView(vpoly, geom);
}