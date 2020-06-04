#pragma once
#include "Query.h"

template <typename T>	//real number types
class PyramidQuery:public Query<T, T> {
	using Query<T, T>::windowquery;
	using Query<T, T>::measure;

public:
	PyramidDB<T> PCDBP;

private:
	T tmin(T a, T b) {
		if (a <= 0 and b >= 0)
			return 0;
		else
			return min<T>(abs(a), abs(b));
	}

	map <T, T> QueryRange(NDWindow<T> window)
	{
		measure.histLoad = 0;
		map <T, T> ranges;
		int ndims = window.nDims;
		T* L = new T[ndims];
		T* H = new T[ndims];
		T QImin = 0;
		T QImax = 0;

		T h_l;
		T h_h;

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
				QImax = min<T>(H[i], 0);
			}
			else
			{
				for (int j = 0; j < ndims; j++)
				{
					if (H[i - ndims] >= tmin(L[j], H[j]) and j != i - ndims) intersect++;
				}
				QImin = max<T>(L[i - ndims], 0);
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
					double Ql = 0.5;
					for (int j = 0; j < ndims; j++)
					{
						if (j != (i%ndims))
						{
							Qh = max<T>(tmin(QImin, QImax), tmin(L[j], H[j]));
							if (Ql > Qh) Ql = Qh;
						}
					}
					h_l = Ql;
					//cout << Ql<<", ";
				}
				h_h = max<T>(abs(QImin), abs(QImax));
				ranges.insert(make_pair(i + h_l, i + h_h));
			}
		}

		measure.rangeNum = ranges.size();
		cout << ranges.size() << endl;

		delete[] L;
		delete[] H;

		return ranges;

	}

public:
	PyramidQuery(const PyramidDB<T>& PC)
	{
		measure = {};
		windowquery = {};
		PCDBP = PC;
	}

	void QueryIOT(const NDWindow<T> & window) override
	{
		windowquery = window;
		short dimnum = PCDBP.nDims;
		NDWindow<T> windowQ = window.Transform(PCDBP.trans);
		if (PCDBP.Extend)
		{
			NDPoint<T> medianP(PCDBP._medians,PCDBP.nDims);
			NDPoint<T> MPshift = medianP.Transform(PCDBP.trans);
			for (int i = 0; i < window.nDims; i++)
			{
				cout << fixed << setprecision(2) << windowQ.minPoint[i] << ", " << windowQ.maxPoint[i] << endl;
				windowQ.minPoint[i] = pow(windowQ.minPoint[i], -1 / log2(MPshift[i]));
				windowQ.maxPoint[i] = pow(windowQ.maxPoint[i], -1 / log2(MPshift[i]));
				cout << fixed << setprecision(2) << windowQ.minPoint[i] << ", " << windowQ.maxPoint[i] << endl;
			}
		}

		auto start = chrono::high_resolution_clock::now();
		map <T, T> ranges;
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
		stmt->setSQL("insert into range_packs values(:l, :u)");
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

		ofstream output_file("E:/query_res_pyramid.csv");
		long long apnum = 0;
		long long spnum = 0;
		NDPoint<T> pt(dimnum);

		while (rs->next())
		{
			apnum++;
			for (int i = 0; i < dimnum; i++) pt[i] = stod(rs->getString(i+2));
			pt[dimnum-1] = stod(rs->getString(dimnum + 2));	//select different dimensions, could be improved later
			//cout << pt[0] << ", " << pt[1] << ", " << pt[2] << ", " << pt[3] << "\n";
			if (this->Inside<T, T>(pt, window) == 1)
			{
				spnum++;
				output_file << fixed << setprecision(2) << pt[0] << ", " << pt[1] << ", " << pt[2] << "\n";
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

};