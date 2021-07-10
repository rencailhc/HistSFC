#include "Simulation.h"
#include "Query.h"
#include "PyramidQuery.h"

#define N (long long) 10000000	//number of points to generate

void simgen() {
	
	PartUSeries sim1_data(6);
	sim1_data.Gen(N);
	sim1_data.Exfile("E:/ISPRS/idealsim/simccv_ideal_6d_6.csv");
	
	/*
	RandSeries sim2_data(6);
	sim2_data.Gen(N);
	sim2_data.Exfile("E:/ISPRS/realsim/simccv_real_6d_1.csv");
	
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
	int simnum = 1;
	//RandWindow<long long> windowSim(2, 0, 1 << maxBits_sim1, 410, "simuni_data");
	//auto windowList = windowSim.Gen(simnum);
	for (int nDims = 8; nDims <= 16; nDims += 2)
	{
		int mDims = nDims;
		RandWindow<long long> windowSim(mDims, 0, 1 << maxBits_sim1, 3500, "simuni_data");
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

void disttest() {
	int simnum = 10;
	/*for 2-nD query test*/
	unsigned seed = chrono::system_clock::now().time_since_epoch().count();
	mt19937 gen(seed);
	uniform_real_distribution<> dis(0.001*(1 << maxBits_sim2), 0.1*(1 << maxBits_sim2));	//for the 2-nD query
	double delta = dis(gen);

	/*************************************/

	for (int nDims = 6; nDims <= 6; nDims++)
	{
			
		RandWindow<long long> windowSim(nDims, 0, 1 << maxBits_sim2, delta, "sim2_sharp_data");
		auto windowList = windowSim.Gen(simnum);
		//auto windowsim = NDWindow<long long>({ 322076, 99033, 71730, 95148, 84557, 91191 }, { 640455, 506425, 230118, 296388, 860998, 237002 });
		/*mixed skew dimensions*/
		PointCloudDB<long long, sfc_bigint> PCDB_ind("sim2_ind_" + to_string(nDims) + "D_iot", nDims);
		auto PCDBhist_ind = PCDB_ind;
		PCDBhist_ind.HIST = true;
		PCDBhist_ind.HistTab = "hist_sim2_ind_" + to_string(nDims) + "D";
		Query<long long, sfc_bigint> testPlain_ind(PCDB_ind);
		Query<long long, sfc_bigint> testHist_ind(PCDBhist_ind);
		string fnSFC_ind = "E:/ISPRS/disttest/plainsfc_" + to_string(nDims) + "D_" + to_string(nDims) + "D_cube_compare_ind.csv";
		string fnSFChist_ind = "E:/ISPRS/disttest/histsfc_" + to_string(nDims) + "D_" + to_string(nDims) + "D_cube_compare_ind.csv";

		/*uniform dimensions*/
		PointCloudDB<long long, sfc_bigint> PCDB_uni("sim2_uni_" + to_string(nDims) + "D_iot", nDims);
		auto PCDBhist_uni = PCDB_uni;
		PCDBhist_uni.HIST = true;
		PCDBhist_uni.HistTab = "hist_sim2_uni_" + to_string(nDims) + "D";
		Query<long long, sfc_bigint> testPlain_uni(PCDB_uni);
		Query<long long, sfc_bigint> testHist_uni(PCDBhist_uni);
		string fnSFC_uni = "E:/ISPRS/disttest/plainsfc_" + to_string(nDims) + "D_" + to_string(nDims) + "D_cube_compare_uni.csv";
		string fnSFChist_uni = "E:/ISPRS/disttest/histsfc_" + to_string(nDims) + "D_" + to_string(nDims) + "D_cube_compare_uni.csv";

		/*sharp dimensions*/
		PointCloudDB<long long, sfc_bigint> PCDB_sharp("sim2_sharp_" + to_string(nDims) + "D_iot", nDims);
		auto PCDBhist_sharp = PCDB_sharp;
		PCDBhist_sharp.HIST = true;
		PCDBhist_sharp.HistTab = "hist_sim2_sharp_" + to_string(nDims) + "D";
		Query<long long, sfc_bigint> testPlain_sharp(PCDB_sharp);
		Query<long long, sfc_bigint> testHist_sharp(PCDBhist_sharp);
		string fnSFC_sharp = "E:/ISPRS/disttest/plainsfc_" + to_string(nDims) + "D_" + to_string(nDims) + "D_cube_compare_sharp.csv";
		string fnSFChist_sharp = "E:/ISPRS/disttest/histsfc_" + to_string(nDims) + "D_" + to_string(nDims) + "D_cube_compare_sharp.csv";

		/*chess board dimensions*/
		/*PointCloudDB<long long, sfc_bigint> PCDB_board("sim2_board_" + to_string(nDims) + "D_iot", nDims);
		auto PCDBhist_board = PCDB_board;
		PCDBhist_board.HIST = true;
		PCDBhist_board.HistTab = "hist_sim2_board_" + to_string(nDims) + "D";
		Query<long long, sfc_bigint> testPlain_board(PCDB_board);
		Query<long long, sfc_bigint> testHist_board(PCDBhist_board);
		string fnSFC_board = "E:/ISPRS/disttest/plainsfc_" + to_string(nDims) + "D_" + to_string(nDims) + "D_cube_compare_board.csv";
		string fnSFChist_board = "E:/ISPRS/disttest/histsfc_" + to_string(nDims) + "D_" + to_string(nDims) + "D_cube_compare_board.csv";
		*/
		for (int i = 0; i < simnum; i++)
		{
			
			testPlain_uni.QueryIOT(windowList[i]);
			testPlain_uni.ExMeasurement_batch(fnSFC_uni);
			testHist_uni.QueryIOT(windowList[i]);
			testHist_uni.ExMeasurement_batch(fnSFChist_uni);
			/*
			testPlain_board.QueryIOT(windowList[i]);
			testPlain_board.ExMeasurement_batch(fnSFC_board);
			testHist_board.QueryIOT(windowList[i]);
			testHist_board.ExMeasurement_batch(fnSFChist_board);
			*/
			testPlain_sharp.QueryIOT(windowList[i]);
			testPlain_sharp.ExMeasurement_batch(fnSFC_sharp);
			testHist_sharp.QueryIOT(windowList[i]);
			testHist_sharp.ExMeasurement_batch(fnSFChist_sharp);
			
			testPlain_ind.QueryIOT(windowList[i]);
			testPlain_ind.ExMeasurement_batch(fnSFC_ind);
			testHist_ind.QueryIOT(windowList[i]);
			testHist_ind.ExMeasurement_batch(fnSFChist_ind);
			
		}
	}
}

void idealsim_ccv()
{
	int simnum = 500;
	int nDims = 6;
	int edge_len = 410;
	RandWindow<long long> windowSim(nDims, 0, 1 << maxBits_sim1, edge_len, "simccv_ideal_6d_6_flat");
	auto windowList = windowSim.Gen_ideal(simnum);

	PointCloudDB<long long, sfc_bigint> PCDB("simccv_ideal_" + to_string(nDims) + "d_6_iot", nDims);
	auto PCDBhist = PCDB;
	PCDBhist.HIST = true;
	PCDBhist.HistTab = "hist_simccv_ideal_" + to_string(nDims) + "d_6";
	Query<long long, sfc_bigint> testPlain(PCDB);
	Query<long long, sfc_bigint> testHist(PCDBhist);

	string fnSFC = "E:/ISPRS/idealsim/plainsfc_" + to_string(nDims) + "D_" + to_string(nDims) + "D_cube_compare_nn_6.csv";
	string fnSFChist = "E:/ISPRS/idealsim/histsfc_" + to_string(nDims) + "D_" + to_string(nDims) + "D_cube_compare_nn_6.csv";

	//double cell_len = (1 << maxBits_sim1)/pow(N, 1.0/nDims);

	for (int i = 0; i < simnum; i++)
	{
		testPlain.QueryIOT(windowList[i]);
		testPlain.ExMeasurement_batch(fnSFC);
		testHist.QueryIOT(windowList[i]);
		testHist.ExMeasurement_batch(fnSFChist);
	}
}

void realsim_ccv()
{
	int simnum = 500;
	int nDims = 5;
	
	RandWindow<long long> windowSim(nDims, 0, 1 << maxBits_sim2, 0, "simccv_real_6d_4");
	auto windowList = windowSim.Gen(simnum);

	PointCloudDB<long long, sfc_bigint> PCDB("simccv_real_" + to_string(nDims) + "d_4_iot", nDims);
	auto PCDBhist = PCDB;
	PCDBhist.HIST = true;
	PCDBhist.HistTab = "hist_simccv_real_" + to_string(nDims) + "d_4";
	Query<long long, sfc_bigint> testPlain(PCDB);
	Query<long long, sfc_bigint> testHist(PCDBhist);

	string fnSFC = "E:/ISPRS/realsim/" + to_string(nDims) + "D/plainsfc_" + to_string(nDims) + "D_" + to_string(nDims) + "D_compare_nn_4.csv";
	string fnSFChist = "E:/ISPRS/realsim/" + to_string(nDims) + "D/histsfc_" + to_string(nDims) + "D_" + to_string(nDims) + "D_compare_nn_4.csv";

	//double cell_len = (1 << maxBits_sim1)/pow(N, 1.0/nDims);

	for (int i = 0; i < simnum; i++)
	{
		testPlain.QueryIOT(windowList[i]);
		testPlain.ExMeasurement_batch(fnSFC);
		testHist.QueryIOT(windowList[i]);
		testHist.ExMeasurement_batch(fnSFChist);
	}
		
	
}