#include "Simulation.h"
#include "Query.h"
#include "PyramidQuery.h"

#define N (long long) 1000000	//number of points to generate

void simgen() {
	/*
	PartUSeries sim1_data(16);
	sim1_data.Gen(N);
	sim1_data.Exfile("E:/sim1_data.csv");
	*/

	RandSeries sim2_data(6);
	sim2_data.Gen(N);
	sim2_data.Exfile("E:/sigmod/simres/sim2/real_cor4.csv");
	/*
	sim2_data.statisPlot(0);
	sim2_data.statisPlot(1);
	sim2_data.statisPlot(2);
	sim2_data.statisPlot(3);
	sim2_data.statisPlot(4);
	sim2_data.statisPlot(5);
	*/
}


void idealsim_uni()
{
	int simnum = 100;
	//RandWindow<long long> windowSim(2, 0, 1 << maxBits_sim1, 410, "simuni_data");
	//auto windowList = windowSim.Gen(simnum);
	for (int nDims = 2; nDims <= 16; nDims += 2)
	{
		int mDims = nDims;
		RandWindow<long long> windowSim(mDims, 0, 1 << maxBits_sim1, 2000, "simuni_data");
		auto windowList = windowSim.Gen(simnum);

		PointCloudDB<long long, sfc_bigint> PCDB("simuni_" + to_string(nDims) + "D_iot", nDims);
		PyramidDB<long long, double> PCDBPyramid("pv_simuni_" + to_string(nDims) + "D_iot", nDims);
		for (int i = 0; i < nDims; i++)
		{
			PCDBPyramid.trans._scale[i] = 1.0f / (1 << maxBits_sim1);
		}

		Query<long long, sfc_bigint> testPlain(PCDB);
		PyramidQuery<long long, double> testPyramid(PCDBPyramid);

		string fnSFC = "E:/sigmod/simres/simuni/plainsfc_" + to_string(mDims) + "D_" + to_string(nDims) + "D_cube_compare_nn.csv";
		string fnPyramid = "E:/sigmod/simres/simuni/pyramid_" + to_string(mDims) + "D_" + to_string(nDims) + "D_cube_compare_nn.csv";
		
		NDPoint<long long> P1(nDims);
		NDPoint<long long> P2(nDims);
		for (int i = 0; i < nDims; i++)
		{
			P1[i] = 0;
			P2[i] = 1 << maxBits_sim1;
		}
		NDWindow<long long> windowFull(P1, P2);

		for (int i = 0; i < simnum; i++)
		{
			for (int j = 0; j < mDims; j++)
			{
				windowFull.minPoint[j] = windowList[i].minPoint[j];
				windowFull.maxPoint[j] = windowList[i].maxPoint[j];
			}
			testPlain.QueryIOT(windowFull);
			testPlain.ExMeasurement_batch(fnSFC);
			testPyramid.QueryIOT(windowFull);
			testPyramid.ExMeasurement_batch(fnPyramid);
		}

	}
}


void idealsim_skew() {
	int simnum = 100;
	//RandWindow<long long> windowSim(2, 0, 1 << maxBits_sim1, 410, "sim1_data");
	//auto windowList = windowSim.Gen(simnum);
	for (int nDims = 2; nDims <= 16; nDims += 2)
	{
		int mDims = nDims;
		RandWindow<long long> windowSim(mDims, 0, 1 << 12, 2000, "sim1_data");
		auto windowList = windowSim.Gen(simnum);

		PointCloudDB<long long, sfc_bigint> PCDB("sim1_" + to_string(nDims) + "D_iot", nDims);
		auto PCDBhist = PCDB;
		PCDBhist.HIST = true;
		PCDBhist.HistTab = "hist_sim1_" + to_string(nDims) + "D";
		PyramidDB<long long, double> PCDBPyramid("pv_sim1_uni_" + to_string(nDims) + "D_iot", nDims);
		PyramidDB<long long, double> PCDBPyramidEx("pv_sim1_shift_" + to_string(nDims) + "D_iot", nDims);
		for (int i = 0; i < nDims; i++)
		{
			PCDBPyramid.trans._scale[i] = 1.0f / (1 << maxBits_sim1);
			PCDBPyramidEx.trans._scale[i] = 1.0f / (1 << maxBits_sim1);
		}
		PCDBPyramidEx.Extend = true;
		double medians[] = { 1757,2186,2606,2253,2074,2106,1383,1753,1436, 2910, 2452, 1217, 2685,2046, 2180, 3089 };  //idealsim-skewed
		PCDBPyramidEx._medians = medians;

		Query<long long, sfc_bigint> testPlain(PCDB);
		Query<long long, sfc_bigint> testHist(PCDBhist);
		PyramidQuery<long long, double> testPyramid(PCDBPyramid);
		PyramidQuery<long long, double> testPyramidEx(PCDBPyramidEx);


		string fnSFC = "E:/sigmod/simres/sim1/plainsfc_" + to_string(mDims) + "D_" + to_string(nDims) + "D_cube_compare_nn.csv";
		string fnSFChist = "E:/sigmod/simres/sim1/histsfc_" + to_string(mDims) + "D_" + to_string(nDims) + "D_cube_compare_nn.csv";
		string fnPyramid = "E:/sigmod/simres/sim1/pyramid_" + to_string(mDims) + "D_" + to_string(nDims) + "D_cube_compare_nn.csv";
		string fnPyramidEx = "E:/sigmod/simres/sim1/pyramidEx_" + to_string(mDims) + "D_" + to_string(nDims) + "D_cube_compare_nn.csv";

		NDPoint<long long> P1(nDims);
		NDPoint<long long> P2(nDims);
		for (int i = 0; i < nDims; i++)
		{
			P1[i] = 0;
			P2[i] = 1 << maxBits_sim1;
		}
		NDWindow<long long> windowFull(P1, P2);

		for (int i = 0; i < simnum; i++)
		{
			for (int j = 0; j < mDims; j++)
			{
				windowFull.minPoint[j] = windowList[i].minPoint[j];
				windowFull.maxPoint[j] = windowList[i].maxPoint[j];
			}
			testPlain.QueryIOT(windowFull);
			testPlain.ExMeasurement_batch(fnSFC);
			testHist.QueryIOT(windowFull);
			testHist.ExMeasurement_batch(fnSFChist);
			testPyramid.QueryIOT(windowFull);
			testPyramid.ExMeasurement_batch(fnPyramid);
			testPyramidEx.QueryIOT(windowFull);
			testPyramidEx.ExMeasurement_batch(fnPyramidEx);
		}

	}
}

void realsim() {
	int simnum = 100;
	/*for 2-nD query test*/
	unsigned seed = chrono::system_clock::now().time_since_epoch().count();
	mt19937 gen(seed);
	uniform_real_distribution<> dis(0.01*(1 << maxBits_sim2), 0.05*(1 << maxBits_sim2));	//for the 2-nD query
	double delta = dis(gen);

	RandWindow<long long> windowSim(2, 0, 1 << maxBits_sim2, delta, "sim2_ind_data");
	auto windowList = windowSim.Gen(simnum);
	/*************************************/

	for (int nDims = 2; nDims <= 6; nDims++)
	{
		int mDims = 2;	//querying dimension mDims, organizing dimension nDims

		//RandWindow<long long> windowSim(mDims, 0, 1 << maxBits_sim2, 0, "pv_sim2_uni_cor");
		//auto windowList = windowSim.Gen(simnum);
		//windowList[0] = NDWindow<long long>({ 84913, 87709, 815598 }, { 130298, 938691, 1037835 });
		
		PointCloudDB<long long, sfc_bigint> PCDB_ind("sim2_ind_" + to_string(nDims) + "D_iot", nDims);
		auto PCDBhist_ind = PCDB_ind;
		PCDBhist_ind.HIST = true;
		PCDBhist_ind.HistTab = "hist_sim2_ind_" + to_string(nDims) + "D";

		PointCloudDB<long long, sfc_bigint> PCDB_cor("sim2_cor_" + to_string(nDims) + "D_iot", nDims);
		auto PCDBhist_cor = PCDB_cor;
		PCDBhist_cor.HIST = true;
		PCDBhist_cor.HistTab = "hist_sim2_cor_" + to_string(nDims) + "D";

		PyramidDB<long long, double> PCDBPyramid_ind("pv_sim2_uni_ind_" + to_string(nDims) + "D_iot", nDims);
		PyramidDB<long long, double> PCDBPyramid_cor("pv_sim2_uni_cor_" + to_string(nDims) + "D_iot", nDims);
		for (int i = 0; i < nDims; i++)
		{
			PCDBPyramid_ind.trans._scale[i] = 1.0f / (1 << maxBits_sim2);
			PCDBPyramid_cor.trans._scale[i] = 1.0f / (1 << maxBits_sim2);
		}

		PyramidDB<long long, double> PCDBPyramidEx_ind("pv_sim2_shift_ind_" + to_string(nDims) + "D_iot", nDims);
		PyramidDB<long long, double> PCDBPyramidEx_cor("pv_sim2_shift_cor_" + to_string(nDims) + "D_iot", nDims);
		for (int i = 0; i < nDims; i++)
		{
			PCDBPyramidEx_ind.trans._scale[i] = 1.0f / (1 << maxBits_sim2);
			PCDBPyramidEx_cor.trans._scale[i] = 1.0f / (1 << maxBits_sim2);
		}
		PCDBPyramidEx_ind.Extend = true;
		PCDBPyramidEx_cor.Extend = true;
		double sim2midians_ind[] = { 524555,524328.5,176956.5,272217,165051.5,879803 };
		double sim2midians_cor[] = { 524301,527249.5,176467,269097,164903,880241 };
		PCDBPyramidEx_ind._medians = sim2midians_ind;
		PCDBPyramidEx_cor._medians = sim2midians_cor;
		
		Query<long long, sfc_bigint> testPlain_ind(PCDB_ind);
		Query<long long, sfc_bigint> testHist_ind(PCDBhist_ind);
		PyramidQuery<long long, double> testPyramid_ind(PCDBPyramid_ind);
		PyramidQuery<long long, double> testPyramidEx_ind(PCDBPyramidEx_ind);
		Query<long long, sfc_bigint> testPlain_cor(PCDB_cor);
		Query<long long, sfc_bigint> testHist_cor(PCDBhist_cor);
		PyramidQuery<long long, double> testPyramid_cor(PCDBPyramid_cor);
		PyramidQuery<long long, double> testPyramidEx_cor(PCDBPyramidEx_cor);

		/*
		time_t t = time(0);   // get time now
		struct tm * now = localtime(&t);
		char buffer[80];
		strftime(buffer, 80, "%Y-%m-%d %Hh%Mm", now);
		string timetag(buffer);
		string filename = "E:/SIMTEST " + timetag + ".txt";
		*/
		
		string fnSFC_ind = "E:/sigmod/simres/sim2/plainsfc_" + to_string(mDims) + "D_" + to_string(nDims) + "D_cube_compare_ind.csv";
		string fnSFChist_ind = "E:/sigmod/simres/sim2/histsfc_" + to_string(mDims) + "D_" + to_string(nDims) + "D_cube_compare_ind.csv";
		string fnPyramid_ind = "E:/sigmod/simres/sim2/pyramid_" + to_string(mDims) + "D_" + to_string(nDims) + "D_cube_compare_ind.csv";
		string fnPyramidEx_ind = "E:/sigmod/simres/sim2/pyramidEx_" + to_string(mDims) + "D_" + to_string(nDims) + "D_cube_compare_ind.csv";
		string fnSFC_cor = "E:/sigmod/simres/sim2/plainsfc_" + to_string(mDims) + "D_" + to_string(nDims) + "D_cube_compare_cor.csv";
		string fnSFChist_cor = "E:/sigmod/simres/sim2/histsfc_" + to_string(mDims) + "D_" + to_string(nDims) + "D_cube_compare_cor.csv";
		string fnPyramid_cor = "E:/sigmod/simres/sim2/pyramid_" + to_string(mDims) + "D_" + to_string(nDims) + "D_cube_compare_cor.csv";
		string fnPyramidEx_cor = "E:/sigmod/simres/sim2/pyramidEx_" + to_string(mDims) + "D_" + to_string(nDims) + "D_cube_compare_cor.csv";
		


		NDPoint<long long> P1(nDims);
		NDPoint<long long> P2(nDims);
		for (int i = 0; i < nDims; i++)
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
			testPyramid_ind.QueryIOT(windowFull);
			testPyramid_ind.ExMeasurement_batch(fnPyramid_ind);
			testPyramidEx_ind.QueryIOT(windowFull); 
			testPyramidEx_ind.ExMeasurement_batch(fnPyramidEx_ind);
			
		/*
			testPlain_cor.QueryIOT(windowFull);
			testPlain_cor.ExMeasurement_batch(fnSFC_cor);
			testHist_cor.QueryIOT(windowFull);
			testHist_cor.ExMeasurement_batch(fnSFChist_cor);
			testPyramid_cor.QueryIOT(windowFull);
			testPyramid_cor.ExMeasurement_batch(fnPyramid_cor);
			testPyramidEx_cor.QueryIOT(windowFull);
			testPyramidEx_cor.ExMeasurement_batch(fnPyramidEx_cor);
			*/
		}
		


	
	}
}