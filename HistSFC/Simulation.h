#pragma once
#include <iostream>
#include <random>
#include <chrono>
#include <fstream>
#include <map>
#include "typedef.h"
#include "Point.h"
#include "Window.h"
#include "math.h"
#include "SFCConversion.h"
#include "Oracle.h"

#define minR (short) 0	//min value of a dimension
#define maxBits_sim1 (short) 12	//max bits of a dimension
#define maxBits_sim2 (short) 20	//max bits of a dimension

using namespace std;

/*Uniform or partial uniform simulation*/

class PartUSeries
{
public:

	PartUSeries()
	{
		nDims = 1;
	}

	PartUSeries(int dimn)
	{
		nDims = dimn;
	}

	void Gen(long long p_num)
	{
		long long* windowL = new long long[nDims];
		long long* windowU = new long long[nDims];
		unsigned seed = chrono::system_clock::now().time_since_epoch().count();
		mt19937 gen(seed);

		uniform_real_distribution<> dis(minR, 1<<maxBits_sim1);
		uniform_real_distribution<> dis2(minR, 1 << (maxBits_sim1-1));

		SFCConversion<long long> sfc(nDims, maxBits_sim1);
		for (int i = 0; i < nDims; i++)
		{
			double D1 = dis2(gen);
			double D2 = D1 + (1 << (maxBits_sim1 - 1));
			/*when partial uniformly distributed*/
			windowL[i] = min(D1, D2);
			windowU[i] = max(D1, D2);
			/*when uniformly distributed*/
			
		}
		

		while (S.size() < p_num)
		{
			long long * arr = new long long[nDims];
			for (int i = 0; i < nDims; i++)
			{
				uniform_real_distribution<> distribution(windowL[i], windowU[i]);
				arr[i] = distribution(gen);
			}

			NDPoint<long long> node(arr, nDims);	//set a dimensionality m <= nDims, to avoid redundancy of the key
			sfc_bigint key = sfc.MortonEncode(node);
			if (S.find(key) == S.end())
				S.insert(make_pair(key, arr));
		}

	}

	void Exfile(string filename)
	{
		ofstream output_file(filename);
		vector<SFCConversion<long long>> sfclst;	//SFC conversions for different dimensionality
		/*for (int i = 2; i <= nDims; i += 2)	//Interval of dimensionality is 2, can be modified to any depending on the simulation
		{
			SFCConversion<long long> sfc(i, maxBits_sim1);
			sfclst.push_back(sfc);
		}*/

		SFCConversion<long long> sfc(nDims, maxBits_sim1);

		for (auto it = S.begin(); it != S.end(); it++)
		{
			long long * arr = it->second;
			/*for (int i = 2; i <= nDims; i += 2)	//Interval of dimensionality is 2, can be modified to any depending on the simulation
			{
				NDPoint<long long> node(arr, i);
				output_file << sfclst[i/2-1].MortonEncode(node) << ", ";
			}*/

			NDPoint<long long> node(arr, nDims);
			output_file << sfc.MortonEncode(node) << ", ";

			for (int j = 0; j < nDims; j++)
			{
				if (j != nDims - 1)
					output_file << it->second[j] << ", ";
				else
					output_file << it->second[j] << '\n';
			}
		}
	}

private:
	int nDims;
	map <sfc_bigint, long long *> S;
};


/*Random simulation*/

enum class Dist {
	Uniform, Gamma1, Gamma2, Gamma3, Normal1, Normal2, Normal3, NormalDelta
};

template <Dist D>
inline long long ValGen2()
{
	unsigned seed = chrono::system_clock::now().time_since_epoch().count();
	mt19937 gen(seed);
	uniform_real_distribution<> distribution(0, 1<<maxBits_sim2);
	long long val = (long long)distribution(gen);
	return val;
};

template<>
inline long long ValGen2<Dist::Gamma1>()
{
	unsigned seed = chrono::system_clock::now().time_since_epoch().count();
	mt19937 generator(seed);
	gamma_distribution<double> distribution(1, 2);
	//gamma_distribution<double> distribution(.05, 1);
	double tmp = distribution(generator) * pow(2, 17);
	while (tmp<minR or tmp>1 << maxBits_sim2)
		tmp = distribution(generator) * pow(2, 17);
	long long val = (long long)tmp;
	return val;
};


template<>
inline long long ValGen2<Dist::Gamma2>()
{
	unsigned seed = chrono::system_clock::now().time_since_epoch().count();
	mt19937 generator(seed);
	gamma_distribution<double> distribution(2, 3);
	//gamma_distribution<double> distribution(10, 0.1);
	double tmp = distribution(generator)* pow(2, 15);
	while (tmp < minR or tmp>1 << maxBits_sim2)
		tmp = distribution(generator)* pow(2, 15);
	long long val = (long long)tmp;
	return val;
};

template<>
inline long long ValGen2<Dist::Gamma3>()
{
	unsigned seed = chrono::system_clock::now().time_since_epoch().count();
	mt19937 generator(seed);
	gamma_distribution<double> distribution(10, 2);
	//gamma_distribution<double> distribution(820, 0.02);
	double tmp = distribution(generator)* pow(2, 16);
	while (tmp < minR or tmp>1 << maxBits_sim2)
		tmp = distribution(generator)* pow(2, 16);
	long long val = (long long)tmp;
	return val;
};


template<>
inline long long ValGen2<Dist::Normal1>()
{
	unsigned seed = chrono::system_clock::now().time_since_epoch().count();
	mt19937 generator(seed);
	normal_distribution<double> distribution(pow(2, 19), pow(3, 0.5)*pow(2, 17));
	//normal_distribution<double> distribution(pow(2, 19), pow(2, 18));
	double tmp = distribution(generator);
	while (tmp < minR or tmp>1 << maxBits_sim2)
		tmp = distribution(generator);
	long long val = (long long)tmp;
	return val;
};

template<>
inline long long ValGen2<Dist::Normal2>()
{
	unsigned seed = chrono::system_clock::now().time_since_epoch().count();
	mt19937 generator(seed);
	normal_distribution<double> distribution(pow(2, 19), pow(2, 18));
	double tmp = distribution(generator);
	while (tmp < minR or tmp>1 << maxBits_sim2)
		tmp = distribution(generator);
	long long val = (long long)tmp;
	return val;
};

template<>
inline long long ValGen2<Dist::Normal3>()
{
	unsigned seed = chrono::system_clock::now().time_since_epoch().count();
	mt19937 generator(seed);
	normal_distribution<double> distribution(pow(2, 19), pow(2, 17));
	double tmp = distribution(generator);
	while (tmp < minR or tmp>1 << maxBits_sim2)
		tmp = distribution(generator);
	auto val = (long long)tmp;
	return val;
};


template<>
inline long long ValGen2<Dist::NormalDelta>()
{
	unsigned seed = chrono::system_clock::now().time_since_epoch().count();
	mt19937 generator(seed);
	normal_distribution<double> distribution(0, pow(2, 17));
	double tmp = distribution(generator);
	long long val = (long long)tmp;
	return val;
};


class RandSeries
{
public:

	RandSeries()
	{
		nDims = 1;
	}

	RandSeries(int dimn)
	{
		nDims = dimn;
	}

	void Gen(long long p_num)
	{
		SFCConversion<long long> sfc(nDims, maxBits_sim2);
		while (S.size() < p_num)
		{
			long long * arr = new long long[nDims];
			arr[0] = ValGen2<Dist::Uniform>();
			arr[1] = ValGen2<Dist::Uniform>();
			arr[2] = ValGen2<Dist::Gamma1>();
			arr[3] = ValGen2<Dist::Normal3>();
			arr[4] = ValGen2<Dist::Gamma2>();
			arr[5] = ValGen2<Dist::Gamma3>();

			//create correlation
			/*arr[1] = arr[0] + ValGen2<Dist::NormalDelta>();
			arr[3] = 0.8 * arr[1];
			while (arr[1] < 0 or arr[3] < 0)
			{
				arr[1] = arr[0] + ValGen2<Dist::NormalDelta>();
				arr[3] = 0.8 * arr[1];
			}
			*/
			//for (int i = 0; i < nDims; i++)
				//arr[i] = ValGen2<Dist::Normal1>();

			NDPoint<long long> node(arr, 3);
			sfc_bigint key = sfc.MortonEncode(node);
			if (S.find(key) == S.end())
				S.insert(make_pair(key, arr));

		}

	}

	int getDim()
	{
		return nDims;
	}

	void print()
	{
		for (auto it = S.begin(); it != S.end(); it++)
		{
			for (int i = 0; i < nDims; i++)
			{
				if (i != nDims - 1)
					cout << it->second[i] << ", ";
				else
					cout << it->second[i] << endl;
			}
		}
	}

	void Exfile(string filename)
	{
		ofstream output_file(filename);
		vector<SFCConversion<long long>> sfclst;	//SFC conversions for different dimensionality
		for (int i = 2; i <= nDims; i ++)	
		{
			SFCConversion<long long> sfc(i, maxBits_sim2);
			sfclst.push_back(sfc);
		}
		for (auto it = S.begin(); it != S.end(); it++)
		{
			long long * arr = it->second;
			for (int i = 2; i <= nDims; i++)
			{
				NDPoint<long long> node(arr, i);
				output_file << sfclst[i-2].MortonEncode(node) << ", ";
			}

			for (int j = 0; j < nDims; j++)
			{
				if (j != nDims - 1)
					output_file << it->second[j] << ", ";
				else
					output_file << it->second[j] << '\n';
			}
		}
	}

	void statisPlot(short tDim)
	{
		const int nbins = 100;

		int p[nbins] = {};

		for (auto it = S.begin(); it != S.end(); it++)
		{
			if ((it->second[tDim] >= minR) && (it->second[tDim] < 1<<maxBits_sim2)) ++p[(int)((it->second[tDim] - minR)*1.0f / ((1 << maxBits_sim2) - minR) * 100)];
		}

		for (int i = 0; i < nbins; ++i) {
			cout << i << "-" << (i + 1) << ": ";
			for (int j = 0; j < p[i] * 100 / S.size(); ++j)
				cout << '*';
			cout << endl;
		}

	}


private:
	int nDims;
	map <sfc_bigint, long long *> S;
};

/*Gernate nD window randomly*/
template<typename T>
class RandWindow {
public:
	short nDims;
	T DimLow;
	T DimHigh;
	T delta;	//length of window box
	string tab;	//table for querying

public:

	RandWindow()
	{
		nDims = 0;
		DimLow = 0;
		DimHigh = 0;
		delta = 0;
		tab = "";
	}

	RandWindow(short dims) : nDims(dims), DimLow(0), DimHigh(0), delta(0), tab("")
	{

	}

	RandWindow(short dims, T dimmin, T dimmax, T delta, string inTab)
	{
		nDims = dims;
		DimLow = dimmin;
		DimHigh = dimmax;
		this->delta = delta;
		tab = inTab;
	}

	vector<NDWindow<T>> Gen(int num)
	{
		Environment *env = Environment::createEnvironment(Environment::DEFAULT);
		Connection  *con = env->createConnection(orclconn().User, orclconn().Password, orclconn().Database);
		Statement *stmt = NULL;
		stmt = con->createStatement();
		stmt->setPrefetchRowCount(FETCH_SIZE);
		ResultSet *rs = NULL;

		string dimnum = to_string(nDims);

		string sql = "select count(*) from " + tab;
		rs = stmt->executeQuery(sql);
		rs->next();
		//long long p_num = stoll(rs->getString(1));
		//stmt->closeResultSet(rs);
		long long p_num = 10000000;

		auto windowL = new T[nDims];
		auto windowU = new T[nDims];
		unsigned seed = chrono::system_clock::now().time_since_epoch().count();
		mt19937 gen(seed);
		vector<NDWindow<T>> windowList;

		auto regionL = new long long[nDims];
		auto regionU = new long long[nDims];
		for (int i = 0; i < 5; i++)
		{
			regionL[i] = 0;
			regionU[i] = 1<< maxBits_sim2;
		}
		regionL[2] = 0;
		regionU[2] = 0.1*pow(2, 19);
		regionL[3] = 0.4*pow(2, 19);
		regionU[3] = 0.6*pow(2, 19);
		regionL[4] = 0;
		regionU[4] = 0.1*pow(2, 19);
		
		NDWindow<T> region(NDPoint<T>(regionL, nDims), NDPoint<T>(regionU, nDims));

		int success = 0;
		while (success < num)
		{
			for (int i = 0; i < nDims; i++)
			{
				uniform_real_distribution<> dis(DimLow, DimHigh);
				T D1 = dis(gen);
				T D2 = dis(gen);
				windowL[i] = min(D1, D2);
				windowU[i] = max(D1, D2);

				////////////
				//windowU[i] = windowL[i] + delta;
				//windowU[i] = max(D1, D2);
				//if (windowU[i] > DimHigh)	windowU[i] = DimHigh;
				///////////
			}
			
			short ncmp = 1;
			for (int i = 0; i < nDims; i++)
			{
				ncmp &= windowL[i] < region.maxPoint[i] && windowU[i] > region.minPoint[i];
			}

			if (ncmp)
			{
				sql = "select count(*) from " + tab + " WHERE ";
				for (int i = 0; i < nDims; i++)
				{
					//column name: a1, a2, a3,...
					if (i != nDims - 1)
						sql += "a" + to_string(i + 1) + ">= " + to_string(windowL[i]) + " and a" + to_string(i + 1) + "<= "
						+ to_string(windowU[i]) + " and ";
					else
						sql += "a" + to_string(i + 1) + ">= " + to_string(windowL[i]) + " and a" + to_string(i + 1) + "<= "
						+ to_string(windowU[i]);
				}
				rs = stmt->executeQuery(sql);
				rs->next();
				long long pnum = stoll(rs->getString(1));
				//long long pnum = 400;

				if (pnum >= 0.00001*p_num and pnum <= 0.001*p_num)
				{
					NDWindow<T> windowRes(NDPoint<T>(windowL, nDims), NDPoint<T>(windowU, nDims));
					windowList.push_back(windowRes);
					success++;
					cout << windowList.size() << "th window, " << pnum << endl;
				}
				stmt->closeResultSet(rs);
			}
		}


		con->terminateStatement(stmt);
		env->terminateConnection(con);
		Environment::terminateEnvironment(env);

		return windowList;
	}


	vector<NDWindow<T>> Gen_ideal(int num)
	{
		auto regionL = new long long[nDims];
		auto regionU = new long long[nDims];

		Environment *env = Environment::createEnvironment(Environment::DEFAULT);
		Connection  *con = env->createConnection(orclconn().User, orclconn().Password, orclconn().Database);
		Statement *stmt = NULL;
		stmt = con->createStatement();
		stmt->setPrefetchRowCount(FETCH_SIZE);
		ResultSet *rs = NULL;

		string dimnum = to_string(nDims);

		string sql = "select ";
		for (int i = 0; i < nDims; i++)
		{
			//column name: a1, a2, a3,...
			if (i != nDims - 1)
				sql += "min(a" + to_string(i + 1) + "), max(a" + to_string(i + 1) + "), ";
			else
				sql += "min(a" + to_string(i + 1) + "), max(a" + to_string(i + 1) + ") ";
		}
		sql += "from " + tab;

		rs = stmt->executeQuery(sql);
		rs->next();
		for (int i = 0; i < nDims; i++)
		{
			regionL[i] = stoll(rs->getString(2*i + 1));
			regionU[i] = stoll(rs->getString(2*i + 2));
		}
		stmt->closeResultSet(rs);

		NDWindow<T> region(NDPoint<T>(regionL, nDims), NDPoint<T>(regionU, nDims));

		auto windowL = new T[nDims];
		auto windowU = new T[nDims];
		unsigned seed = chrono::system_clock::now().time_since_epoch().count();
		mt19937 gen(seed);
		vector<NDWindow<T>> windowList;

		int success = 0;
		while (success < num)
		{
			for (int i = 0; i < nDims; i++)
			{
				uniform_real_distribution<> dis(DimLow, DimHigh);
				T D1 = dis(gen);
				windowL[i] = D1;
				windowU[i] = windowL[i] + delta;
				if (windowU[i] > DimHigh)	windowU[i] = DimHigh;
			}

			short ncmp = 1;
			for (int i = 0; i < nDims; i++)
			{
				ncmp &= windowL[i] < region.maxPoint[i] && windowU[i] > region.minPoint[i];
			}
			
			if (ncmp)
			{
				NDWindow<T> windowRes(NDPoint<T>(windowL, nDims), NDPoint<T>(windowU, nDims));
				windowList.push_back(windowRes);
				success++;
			}

		}


		con->terminateStatement(stmt);
		env->terminateConnection(con);
		Environment::terminateEnvironment(env);

		return windowList;
	}


};


void simgen();
void idealsim_uni();
void idealsim_skew();
void realsim();
void disttest();
void idealsim_ccv();
void realsim_ccv();