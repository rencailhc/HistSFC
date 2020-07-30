#pragma once
#include <map>
#include <string> 
#include <fstream>
#include <iomanip>
#include <chrono>
#include <cctype>
#include "Window.h"
#include "PointCloud.h"
#include "Oracle.h"
#include "BaseStruct.h"

template <typename T, typename U>	
class PyramidQuery {
protected:
	Measurement measure;	//collecting performance info
	NDWindow<T> windowquery;	//query window
public:
	PyramidDB<T, U> PCDBP;

private:
	U tmin(U a, U b) {
		if (a <= 0 and b >= 0)
			return 0;
		else
			return min<U>(abs(a), abs(b));
	}

	map <U, U> QueryRange(NDWindow<double> window)
	{
		measure.histLoad = 0;
		map <U, U> ranges;
		int ndims = window.nDims;
		U* L = new U[ndims];
		U* H = new U[ndims];
		U QImin = 0;
		U QImax = 0;

		U h_l;
		U h_h;

		for (int i = 0; i < ndims; i++)
		{
			L[i] = window.minPoint[i] - 0.5;
			H[i] = window.maxPoint[i] - 0.5;
			if (H[i] > 0.5) H[i] = 0.5; //in case the query box exeeds the data region
		}

		for (int i = 0; i < 2 * ndims; i++)
		{
			int intersect = 0;
			if (i < ndims)
			{
				for (int j = 0; j < ndims; j++)
				{
					if (L[i] <= -tmin(L[j], H[j]) and j != i) intersect++;
				}
				QImin = L[i];
				QImax = min<U>(H[i], 0);
			}
			else
			{
				for (int j = 0; j < ndims; j++)
				{
					if (H[i - ndims] >= tmin(L[j], H[j]) and j != i - ndims) intersect++;
				}
				QImin = max<U>(L[i - ndims], 0);
				QImax = H[i - ndims];
			}

			if (intersect == ndims - 1)
			{
				int m = 0;
				for (int j = 0; j < ndims; j++)
				{
					if (L[j] <= 0 and H[j] >= 0) m++;
				}
				//cout << m << ",";
				if (m == ndims)
				{
					h_l = 0;
				}
				else
				{
					double Qh = 0;
					double Ql = 0;
					for (int j = 0; j < ndims; j++)
					{
						if (j != (i%ndims))
						{
							Qh = max<U>(tmin(QImin, QImax), tmin(L[j], H[j]));
							if (Ql < Qh) Ql = Qh;
						}
					}
					h_l = Ql;
					//cout << Ql<<", ";
				}
				h_h = max<U>(abs(QImin), abs(QImax));
				ranges.insert(make_pair(i + h_l, i + h_h));
			}
		}

		measure.rangeNum = ranges.size();

		delete[] L;
		delete[] H;

		return ranges;

	}

protected:
	template<typename A, typename B>
	bool Inside(NDPoint<A> &pt, const NDWindow<B> &wd)	//used for second filter
	{
		if (pt.returnSize() != wd.nDims)
		{
			throw("Dimensionality does not match!");
		}
		else
		{
			short ncmp = 1;
			for (int i = 0; i < wd.nDims; i++)
			{
				ncmp &= pt[i] >= wd.minPoint[i] && pt[i] <= wd.maxPoint[i];
			}
			if (ncmp) return 1;
		}
		return 0;
	}

public:
	PyramidQuery(const PyramidDB<T, U>& PC)
	{
		measure = {};
		windowquery = {};
		PCDBP = PC;
	}

	void QueryIOT(const NDWindow<T> & window)
	{
		windowquery = window;
		short dimnum = PCDBP.nDims;
		NDWindow<double> windowQ = window.Transform<double>(PCDBP.trans);
		for (int i = 0; i < windowQ.nDims; i++)
		{
			if (windowQ.minPoint[i] < 0) windowQ.minPoint[i] = 0;
			if (windowQ.maxPoint[i] > 1) windowQ.maxPoint[i] = 1;
		}

		if (PCDBP.Extend)
		{
			NDPoint<double> medianP(PCDBP._medians,PCDBP.nDims);
			auto MPshift = medianP.Transform<double>(PCDBP.trans);
			for (int i = 0; i < window.nDims; i++)
			{
				windowQ.minPoint[i] = pow(windowQ.minPoint[i], -1 / log2(MPshift[i]));
				windowQ.maxPoint[i] = pow(windowQ.maxPoint[i], -1 / log2(MPshift[i]));
			}
		}

		auto start = chrono::high_resolution_clock::now();
		map <U, U> ranges;
		ranges = QueryRange(windowQ);
		auto end = chrono::high_resolution_clock::now();
		measure.rangeComp = chrono::duration_cast<chrono::milliseconds>(end - start).count();

		Environment *env = Environment::createEnvironment(Environment::DEFAULT);
		Connection  *con = env->createConnection(orclconn().User, orclconn().Password, orclconn().Database);

		auto start1 = chrono::high_resolution_clock::now();

		Statement *f_stmt = NULL;
		f_stmt = con->createStatement();
		string sql = "create table range_packs (lower number, upper number)";
		f_stmt->executeUpdate(sql);
		con->terminateStatement(f_stmt);
		Statement *stmt = con->createStatement();
		stmt->setSQL("insert into range_packs values(round(:l," + to_string(PREC) + "), round(:u," + to_string(PREC) +"))");
		stmt->setMaxIterations(ranges.size() + 1);
		stmt->setMaxParamSize(1, 100);
		stmt->setMaxParamSize(2, 100);
		for (auto it = ranges.begin(); it != ranges.end(); it++)
		{
			double pv_l = it->first;
			stmt->setDouble(1, pv_l);
			double pv_h = it->second;
			stmt->setDouble(2, pv_h);
			stmt->addIteration();
		}
		stmt->executeUpdate();

		sql = "select /*+ use_nl (t r) */ * from " + PCDBP.Table + " t, range_packs r where (t.pv >= r.lower and t.pv <= r.upper)";
		stmt->setSQL(sql);
		stmt->setPrefetchRowCount(FETCH_SIZE);
		ResultSet *rs = stmt->executeQuery();

		auto end1 = chrono::high_resolution_clock::now();
		measure.firstCost = measure.rangeComp + chrono::duration_cast<chrono::milliseconds>(end1 - start1).count();

		//ofstream output_file("E:/query_res_pyramid.csv");
		string outp = "E:/query_res_pyramid_uni.csv";
		if(PCDBP.Extend)
			outp = "E:/query_res_pyramid_extend.csv";
		ofstream output_file(outp);

		long long apnum = 0;
		long long spnum = 0;
		NDPoint<T> pt(dimnum);

		while (rs->next())
		{
			apnum++;
			for (int i = 0; i < dimnum; i++) pt[i] = rs->getDouble(i+2);
			//pt[dimnum-1] = rs->getDouble(dimnum + 2);	//select different dimensions, could be improved later
			//cout << pt[0] << ", " << pt[1] << ", " << pt[2] << ", " << pt[3] << "\n";
			//output_file << pt[0] << ", " << pt[1] << ", " << pt[2] << "\n";
			if (this->Inside<T, T>(pt, window) == 1)
			{
				spnum++;
				//output_file << pt[0] << ", " << pt[1] << ", " << pt[2] << "\n";
			}
		}
		auto end2 = chrono::high_resolution_clock::now();
		measure.secondCost = chrono::duration_cast<chrono::milliseconds>(end2 - end1).count();
		measure.appPNum = apnum;
		measure.accPNum = spnum;
		measure.FPR = (apnum - spnum)*1.0f / spnum;

		stmt->closeResultSet(rs);
		stmt->executeUpdate("drop table range_packs");
		con->commit();
		con->terminateStatement(stmt);
		env->terminateConnection(con);
		Environment::terminateEnvironment(env);

	}

	void ExMeasurement(string filename) 
	{
		ofstream output(filename, ios::app | ios::out);
		output << "Query table: " << PCDBP.Table << ", query geometry: " << "[" << fixed << setprecision(2);
		for (int i = 0; i < windowquery.nDims; i++)
		{
			output << windowquery.minPoint[i] << " ";
		}
		output << ", ";
		for (int i = 0; i < windowquery.nDims; i++)
		{
			output << windowquery.maxPoint[i] << " ";
		}
		output << "]\n";
		if (PCDBP.Extend) output << "Extended Pyramid-T\n";
		else output << "Pyramid-T\n";

		output << "rangeNum, appPNum, accPNum, FPR, rangeComp, firstCost, secondCost\n";
		output << measure.rangeNum << ", " << measure.appPNum << ", " << measure.accPNum << ", " << measure.FPR << ", "
			<< measure.rangeComp << ", " << measure.firstCost << ", " << measure.secondCost << "\n";
		output << "\n";
	}

	void ExMeasurement_batch(string filename) 
	{
		ofstream output(filename, ios::app | ios::out);
		output << measure.rangeNum << ", " << measure.appPNum << ", " << measure.accPNum << ", " << measure.FPR << ", "
			<< measure.rangeComp << ", " << measure.firstCost << ", " << measure.secondCost << "\n";
		//output << measure.FPR << ",";
	}

};