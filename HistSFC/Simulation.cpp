#include "Simulation.h"
#include "Query.h"
#include "PyramidQuery.h"

#define N (long long) 1000	//number of points to generate

void simgen() {
	PartUSeries sim1_data(16);
	sim1_data.Gen(N);
	sim1_data.Exfile("E:/sim1_data.csv");

	RandSeries sim2_data(6);
	sim2_data.Gen(N);
	sim2_data.Exfile("E:/sim2_cor_data.csv");
	sim2_data.statisPlot(0);
	sim2_data.statisPlot(1);
	sim2_data.statisPlot(2);
	sim2_data.statisPlot(3);
	sim2_data.statisPlot(4);
	sim2_data.statisPlot(5);
}


void simtest() {
	int simnum = 100;
	RandWindow windowSim(2, 0, 1 << 12, 200, "simuni_data");
	auto windowList = windowSim.Gen(simnum);

	for (int nDims = 2; nDims <= 16; nDims += 2)
	{
		int mDims = nDims;
		string output = "E:/sigmod/simres/simuni/" + to_string(mDims) + "D" + to_string(nDims) + "D_cube_compare.csv";
		ofstream output_file(output);

		PointCloudDB<long long, sfc_bigint> PCDB("simuni_" + to_string(nDims) + "D_iot", nDims);
		auto PCDBhist = PCDB;
		PCDBhist.HIST = true;
		PCDBhist.HistTab = "hist_sim2_ind_" + to_string(nDims) + "D";
		PyramidDB<double> PCDBPyramid("pv_simuni_" + to_string(nDims) + "D_iot", nDims);
		for (int i = 0; i < nDims; i++)
		{
			PCDBPyramid.trans._scale[i] = 1.0f / (1 << maxBits_sim1);
		}
		//PyramidDB<double> PCDBPyramidEx("pv_extend_simuni_" + to_string(nDims) + "D_iot", nDims);
		
		Query<long long, sfc_bigint> testPlain(PCDB);
		Query<long long, sfc_bigint> testHist(PCDBhist);
		PyramidQuery<double> testPyramid(PCDBPyramid);
		//PyramidQuery<double> testPyramidEx(PCDBPyramidEx);

		time_t t = time(0);   // get time now
		struct tm * now = localtime(&t);
		char buffer[80];
		strftime(buffer, 80, "%Y-%m-%d %Hh%Mm", now);
		string timetag(buffer);
		string filename = "E:/AHN2TEST " + timetag + ".txt";
		for (int i = 0; i < simnum; i++)
		{
			testHist.QueryIOT(windowList[i]);
			testHist.ExMeasurement(filename);
			//testPlain.QueryIOT(windowList[i]);
			//testPlain.ExMeasurement(filename);
			testPyramid.QueryIOT(windowList[i]);
			testPyramid.ExMeasurement(filename);
			testPyramid_uni.QueryIOT(windowList[i]);
			testPyramid_uni.ExMeasurement(filename);
		}
		


	
	}
}