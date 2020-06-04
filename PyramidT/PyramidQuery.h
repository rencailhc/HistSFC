#pragma once
#include <iostream>
#include <fstream>
#include <iomanip>
#include <map>
#include "oracle.h"
#include "window.h"
#include "PointCloud.h"

template <typename T>
class PyramidQuery {
public:
	PointCloudDB PCDB;

private:
	T tmin(T a, T b) {
		if (a <= 0 and b >= 0)
			return 0;
		else
			return min<T>(abs(a), abs(b));
	}

	map <T, T> QueryRange(NDWindow<T> window)
	{
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

		delete[] L;
		delete[] H;

		return ranges;

	}

public:
	long long QueryIOT(NDWindow<T> window)
	{
		int dimnum = window.nDims;
		map <T, T> ranges;

		ranges = QueryRange<T>(window);

		Environment *env = Environment::createEnvironment(Environment::DEFAULT);
		Connection  *con = env->createConnection(User, Password, Database);

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

		sql = "select /*+ use_nl (t r) */ t.pv from " + iot_table + " t, range_packs r where (t.pv >= r.lower and t.pv <= r.upper)";
		stmt->setSQL(sql);
		stmt->setPrefetchRowCount(FETCH_SIZE);
		ResultSet *rs = stmt->executeQuery();

		long long apnum = 0;
		while (rs->next()) apnum++;

		stmt->closeResultSet(rs);
		stmt->executeUpdate("drop table range_packs");
		con->commit();
		con->terminateStatement(stmt);
		env->terminateConnection(con);
		Environment::terminateEnvironment(env);

		return apnum;
	}

};