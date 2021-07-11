/*4D window querying based on the 4D PostGIS solution*/

#pragma once
#include <string> 
#include <fstream>
#include <chrono>
#include <cctype>
#include <stdio.h>
#include <stdlib.h>
#include <libpq-fe.h>
#include "../Window.h"
#include "../PointCloud.h"
#include "../BaseStruct.h"

template <typename T, typename U>
class PG4DQuery {
protected:
	Measurement measure;	//collecting performance info
	NDWindow<T> windowquery;	//query window
public:
	PointCloudDB<T, U> pgPCDB;

private:
	bool Inside(double *pt, const NDWindow<T> &wd)	//used for second filter
	{
		short ncmp = 1;
		for (int i = 0; i < wd.nDims; i++)
		{
			ncmp &= pt[i] >= wd.minPoint[i] && pt[i] <= wd.maxPoint[i];
		}
		if (ncmp) return 1;

		return 0;
	}

	vector<double *> inside(string mpoint, NDWindow<T> window) {
		/*filter out all points inside*/
		vector<double *> insider;
		string delimiter1 = "(";
		string delimiter2 = ",";
		string delimiter3 = " ";
		string delimiter4 = ")";

		size_t pos1 = 0;
		size_t pos2 = 0;
		string token1, token2;
		int i = 0;
		if ((pos1 = mpoint.find(delimiter1)) != string::npos) {
			mpoint.erase(0, pos1 + delimiter1.length());
			while ((pos1 = mpoint.find(delimiter2)) != string::npos) {
				token1 = mpoint.substr(0, pos1);
				i = 0;
				double *pt = new double[4];
				while ((pos2 = token1.find(delimiter3)) != string::npos) {
					token2 = token1.substr(0, pos2);
					pt[i] = stod(token2);
					i++;
					token1.erase(0, pos2 + delimiter3.length());
				}
				pt[i] = stod(token1);
				measure.appPNum++;
				if (Inside(pt, window))
					insider.push_back(pt);
					
				mpoint.erase(0, pos1 + delimiter2.length());
			}

			if ((pos1 = mpoint.find(delimiter4)) != string::npos){
				token1 = mpoint.substr(0, pos1);
				i = 0;
				double *pt = new double[4];
				while ((pos2 = token1.find(delimiter3)) != string::npos) {
					token2 = token1.substr(0, pos2);
					pt[i] = stod(token2);
					i++;
					token1.erase(0, pos2 + delimiter3.length());
				}
				pt[i] = stod(token1);
				measure.appPNum++;
				if (Inside(pt, window))
					insider.push_back(pt);

			}
			
		}

		return insider;
	};

public:
	PG4DQuery() {
		measure = {};
		windowquery = {};
		pgPCDB = {};
	}

	PG4DQuery(const PointCloudDB<T, U>& PC)
	{
		measure = {};
		windowquery = {};
		pgPCDB = PC;
	}

	void QueryBlock(const NDWindow<T> & window)
	{
		measure.rangeNum = 0;
		measure.histLoad = 0;
		measure.rangeComp = 0;
		measure.appPNum = 0;
		measure.accPNum = 0;
		windowquery = window;

		PGconn *conn = PQconnectdb("user=haicheng dbname=haicheng");

		if (PQstatus(conn) == CONNECTION_BAD) {

			fprintf(stderr, "Connection to database failed: %s\n",
				PQerrorMessage(conn));
			PQfinish(conn);
			exit(1);
		}

		auto start = chrono::high_resolution_clock::now();
		string windowRep = "LINESTRING(" + to_string(window.minPoint[0]) + " " + to_string(window.minPoint[1]) + " " + to_string(window.minPoint[2]) + " " + to_string(window.minPoint[3]) + ", "
			+ to_string(window.maxPoint[0]) + " " + to_string(window.maxPoint[1]) + " " + to_string(window.maxPoint[2]) + " " + to_string(window.maxPoint[3]) + ")";
		PGresult *res = PQexec(conn, ("select st_asewkt(mpoint) from " + pgPCDB.Table + " where mpoint &&& '" + windowRep + "'").c_str()); 
		auto end1 = chrono::high_resolution_clock::now();
		measure.firstCost = chrono::duration_cast<chrono::milliseconds>(end1 - start).count();

		if (PQresultStatus(res) != PGRES_TUPLES_OK) {

			printf("No data retrieved\n");
			PQclear(res);
			PQfinish(conn);
			exit(1);
		}

		int rows = PQntuples(res);

		for (int i = 0; i < rows; i++) {
			string mpoint = PQgetvalue(res, i, 0);
			auto insider = inside(mpoint, window);
			for (auto it = insider.begin(); it != insider.end(); it++)
			{
				measure.accPNum++;
				//cout << (*it)[0] << ", " << (*it)[1] << ", " << (*it)[2] << ", " << (*it)[3] << endl;
			}

		}
		auto end2 = chrono::high_resolution_clock::now();
		measure.secondCost = chrono::duration_cast<chrono::milliseconds>(end2 - end1).count();
		measure.FPR= (measure.appPNum - measure.accPNum)*1.0f / measure.accPNum;

		PQclear(res);

		res = PQexec(conn, ("select count(*), sum(ST_NPoints(mpoint)) from " + pgPCDB.Table + " where mpoint &&& '" + windowRep + "'").c_str());
		measure.rangeNum = atoi(PQgetvalue(res, 0, 0));
		measure.appPNum = atol(PQgetvalue(res, 0, 1));

		PQclear(res);
		PQfinish(conn);
	}

	void ExMeasurement(string filename)
	{
		ofstream output(filename, ios::app | ios::out);
		output << "Query table: " << pgPCDB.Table << ", query geometry: " << "[" << fixed << setprecision(2);
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
		output << "PostGIS-Rtree\n";

		output << "rangeNum, appPNum, accPNum, FPR, rangeComp, histLoad, firstCost, secondCost\n";
		output << measure.rangeNum << ", " << measure.appPNum << ", " << measure.accPNum << ", " << measure.FPR << ", "
			<< measure.rangeComp << ", " << measure.histLoad << ", " << measure.firstCost << ", " << measure.secondCost << "\n";
		output << "\n";
	}

	void ExMeasurement_batch(string filename)
	{
		ofstream output(filename, ios::app | ios::out);
		output << measure.rangeNum << ", " << measure.appPNum << ", " << measure.accPNum << ", " << measure.FPR << ", "
			<< measure.rangeComp << ", " << measure.histLoad << ", " << measure.firstCost << ", " << measure.secondCost << "\n";
	}

};
