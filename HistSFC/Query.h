#pragma once
#include <map>
#include <fstream>
#include <iomanip>
#include <chrono>
#include "Window.h"
#include "PointCloud.h"
#include "Oracle.h"
#include "BaseStruct.h"
#include "SFCConversion.h"

#define MAX_CYCLE (unsigned int) 1000000

template <typename T, typename U>   //type for window and database key, respectively
class Query {
protected:
	Measurement measure;	//collecting performance info
	NDWindow<T> windowquery;	//query window
public:
	PointCloudDB<T, U> PCDB;

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

private:
	template<typename A, typename B>
	short Intersect(const NDWindow<A> & queryrect, const NDWindow<B> & node)
	{
		short ncmp = 1;
		short dimnum = queryrect.nDims;
		for (int i = 0; i < dimnum; i++)
		{
			ncmp &= node.minPoint[i] >= queryrect.minPoint[i] && node.maxPoint[i] <= queryrect.maxPoint[i];
		}
		if (ncmp) return 1;   //contain

		/*
		intersect:
		//http://stackoverflow.com/questions/306316/determine-if-two-rectangles-overlap-each-other
		RectA.Left < RectB.Right && RectA.Right > RectB.Left && RectA.Top > RectB.Bottom && RectA.Bottom < RectB.Top
		this can be extended more dimensions
		//http://stackoverflow.com/questions/5009526/overlapping-cubes
		if (nrt.x0 < qrt.x1 && nrt.x1 > qrt.x0 &&
		nrt.y0 < qrt.y1 && nrt.y1 > qrt.y0)
		return 2;
		*/
		ncmp = 1;
		for (int i = 0; i < dimnum; i++)
		{
			ncmp &= node.minPoint[i] < queryrect.maxPoint[i] && node.maxPoint[i] > queryrect.minPoint[i];
		}
		if (ncmp)
		{
			for (int i = 0; i < dimnum; i++)
			{
				ncmp += node.minPoint[i] >= queryrect.minPoint[i] && node.maxPoint[i] <= queryrect.maxPoint[i];
			}

			return 1 + ncmp;  //overlap
		}

		//not overlap
		return 0;
	}

	HistNodeND *HistLoad(string HistTab)
	{
		HistNodeND *HistRoot = (HistNodeND*)malloc(sizeof(HistNodeND));
		map <long long, HistNodeND *> HistNodePool;

		try
		{
			Environment *env = Environment::createEnvironment(Environment::DEFAULT);

			Connection  *con = env->createConnection(orclconn().User, orclconn().Password, orclconn().Database);
			Statement *stmt = NULL;
			ResultSet *rs = NULL;
			string sql = "select id from " + HistTab + " where sfc = 0";
			stmt = con->createStatement();
			stmt->setPrefetchRowCount(FETCH_SIZE);
			rs = stmt->executeQuery(sql);
			rs->next();

			long long id = stoll(rs->getString(1));
			HistNodePool.insert(make_pair(id, HistRoot));

			sql = "select * from " + HistTab;
			rs = stmt->executeQuery(sql);

			while (rs->next())
			{
				HistNodeND *curnode;
				id = stoll(rs->getString(1));
				auto it = HistNodePool.find(id);
				if (it != HistNodePool.end())
				{
					curnode = it->second;
				}
				else
				{
					curnode = (HistNodeND*)malloc(sizeof(HistNodeND));
					HistNodePool.insert(make_pair(id, curnode));
				}
				curnode->key = (sfc_bigint)rs->getString(2);
				curnode->pnum = stoll(rs->getString(3));
				curnode->cnum = rs->getInt(4);
				curnode->height = rs->getInt(5);

				long long child_id = stoll(rs->getString(6));
				it = HistNodePool.find(child_id);
				if (it != HistNodePool.end()) curnode->child = it->second;
				else
				{
					HistNodeND *child = (HistNodeND*)malloc(sizeof(HistNodeND));
					curnode->child = child;
					HistNodePool.insert(make_pair(child_id, child));
				}

				long long neighbor_id = stoll(rs->getString(7));
				it = HistNodePool.find(neighbor_id);
				if (it != HistNodePool.end()) curnode->neighbor = it->second;
				else
				{
					HistNodeND *neighbor = (HistNodeND*)malloc(sizeof(HistNodeND));
					curnode->neighbor = neighbor;
					HistNodePool.insert(make_pair(neighbor_id, neighbor));
				}

			}

			HistNodePool.clear();

			stmt->closeResultSet(rs);
			con->terminateStatement(stmt);
			env->terminateConnection(con);
			Environment::terminateEnvironment(env);
		}

		catch (std::exception &ex)
		{
			cout << ex.what() << endl;
		}

		return HistRoot;
	}


	map <U, U> PlainWindowRange(NDWindow<T> const & window, short maxdepth, short dimbits = 20)
	{
		measure.histLoad = 0;
		map <U, U> ranges;
		int dimnum = PCDB.nDims;
		vector <NodeND> SearchT;
		short order = dimbits;
		NDPoint<int> NodeL(dimnum); //lower bound of child node
		NDPoint<int> NodeH(dimnum); //upper bound of child node
		SFCConversion<int> sfc(dimnum, order + 1); //with 1 more bit for the sake of rounding the double type

		if (!maxdepth) maxdepth = order + 1;

		NodeND rootnode = {0, order};
		SearchT.push_back(rootnode);
		NDWindow<int> cell;
		NodeND node;
		sfc_bigint rangeL;   //lower boundary of a range 
		sfc_bigint rangeH;   //upper boundary of a range 
		int i = 0;

		unsigned int cycle_time = MAX_CYCLE;  //threshold for the refinement cycles
		unsigned int cycle = SearchT.size();
		vector <NodeND> SearchTC;
		while (!SearchT.empty())
		{
			while (!SearchT.empty())
			{
				i++;
				node = SearchT.back();
				//free(SearchT.back());   //remove node built from malloc
				SearchT.pop_back();

				rangeL = node.key << (node.height*dimnum);
				rangeH = ((node.key + 1) << (node.height*dimnum)) - 1;
				NodeL = sfc.MortonDecode(rangeL);
				NodeH = sfc.MortonDecode(rangeH);
				cell.SetMinPoint(NodeL);
				cell.SetMaxPoint(NodeH);

				if (Intersect<T, int>(window, cell) == 1)
				{
					ranges.insert(make_pair(rangeL, rangeH));
				}

				else if (Intersect<T, int>(window, cell) > 1)
				{
					if (order - node.height < maxdepth and cycle <= cycle_time)
					{
						for (int i = 0; i < 1 << dimnum; i++)
						{
							NodeND child = { (node.key << dimnum) + i, node.height-1 };
							SearchTC.push_back(child);
							cycle++;
						}
					}
					else ranges.insert(make_pair(rangeL, rangeH));

				}

				else continue;

			}
			SearchT.swap(SearchTC);
			SearchTC.clear();

		}
		measure.rangeNum = cycle;

		map <U, U> ranges_merged;
		U sfc_s = ranges.begin()->first;
		U sfc_e = ranges.begin()->second;
		for (auto it = ranges.begin(); it != ranges.end(); ++it)
		{
			if (it->first - sfc_e < 2) sfc_e = it->second;
			else
			{
				ranges_merged.insert(make_pair(sfc_s, sfc_e));
				sfc_s = it->first;
				sfc_e = it->second;
			}
		}
		ranges_merged.insert(make_pair(sfc_s, sfc_e));
		return ranges;
	}

	map <U, U> HistWindowRange(NDWindow<T> const & window, short dimbits = 20)
	{
		//Directly search the histogram tree
		auto start = chrono::high_resolution_clock::now();
		map <U, U> ranges;
		int dimnum = PCDB.nDims;
		vector <HistNodeND *> bNodes;  //nodes located on the boundary
		vector <HistNodeND *> SearchT;
		NDPoint<int> NodeL(dimnum);
		NDPoint<int> NodeH(dimnum);
		sfc_bigint rangeL;
		sfc_bigint rangeH;
		SFCConversion<int> sfc(dimnum, dimbits);

		NDWindow<int> cell;
		HistNodeND *histroot = HistLoad(PCDB.HistTab);

		auto end = chrono::high_resolution_clock::now();
		measure.histLoad = chrono::duration_cast<chrono::milliseconds>(end - start).count();

		SearchT.push_back(histroot);
		HistNodeND *node;
		long long sumBP = 0;  //points covered by boundary nodes
		long long sumIP = 0;  //points covered by inner nodes

		while (!SearchT.empty())
		{
			node = SearchT.back();
			SearchT.pop_back();
			rangeL = node->key << (node->height*dimnum);
			rangeH = (((node->key + 1) << (node->height*dimnum)) - 1);
			NodeL = sfc.MortonDecode(rangeL);
			NodeH = sfc.MortonDecode(rangeH);
			cell.SetMinPoint(NodeL);
			cell.SetMaxPoint(NodeH);

			if (Intersect<T, int>(window, cell) == 1)
			{
				ranges.insert(make_pair(rangeL, rangeH));
				sumIP += node->pnum;
			}

			else if (Intersect<T, int>(window, cell) > 1)
			{
				if (node->cnum != 0)
				{
					HistNodeND *curnode = node->child;
					for (int i = 0; i < node->cnum; i++)
					{
						SearchT.push_back(curnode);
						curnode = curnode->neighbor;
					}
				}
				else
				{
					//ranges.insert(make_pair(keyL, keyH));
					node->cnum = dimnum - Intersect(window, cell) + 2;  //store the number of plains intersects temporarily
					bNodes.push_back(node);
					sumBP += node->pnum / pow(2, node->cnum);
					//cout << "Space: " << NodeLR[0] << "," << NodeLR[1] << "," << NodeLR[2] << "," << NodeLR[3] << ";" << NodeHR[0] << "," << NodeHR[1] << "," << NodeHR[2] << "," << NodeHR[3] << ", number of points: " << node->pnum << ", ratio: "<< 1.0*node->pnum/(double)(keyH - keyL) << endl;
				}

			}

			else continue;
		}

		auto end1 = chrono::high_resolution_clock::now();
		cout << "Hist search costs: " << chrono::duration_cast<chrono::milliseconds>(end1 - end).count() << "ms" << endl;
		cout << "Inner points: " << sumIP << ", boundary points: " << sumBP << ", boundary nodes: " << bNodes.size() << endl;
		long long accuracy = trunc((sumIP + sumBP) * 0.5);  //accuracy threshold
		double avgthres = accuracy / bNodes.size();
		cout << "accuracy: " << accuracy << endl;
		for (int i = 0; i < bNodes.size(); i++)
		{
			//cout << "number of points: "<<bNodes[i]->pnum << ", interdim: " << bNodes[i]->cnum << ", ";
			bNodes[i]->cnum = bNodes[i]->height - trunc(log2(bNodes[i]->pnum / avgthres) / bNodes[i]->cnum) - 1;  //store the original height in cnum temporarily																							  //sumADnodes += pow(16, bNodes[i]->height - bNodes[i]->cnum);
		}

		unsigned int cycle_time = MAX_CYCLE;  //threshold for the refinement cycles
		unsigned int cycle = bNodes.size();
		if (sumBP > accuracy)     //10^7 is a threshold for performance
		{
			vector <HistNodeND *> bNodesC;   //children pool
			while (!bNodes.empty())
			{
				while (!bNodes.empty())
				{
					node = bNodes.back();
					bNodes.pop_back();
					rangeL = node->key << (node->height*dimnum);
					rangeH = (((node->key + 1) << (node->height*dimnum)) - 1);
					NodeL = sfc.MortonDecode(rangeL);
					NodeH = sfc.MortonDecode(rangeH);
					cell.SetMinPoint(NodeL);
					cell.SetMaxPoint(NodeH);

					if (Intersect<T, int>(window, cell) == 1)
					{
						ranges.insert(make_pair(rangeL, rangeH));
					}
					else if (Intersect<T, int>(window, cell) > 1)
					{
						if (node->height > node->cnum and cycle <= cycle_time)
						{
							//cout << height << ", " << node->cnum << endl;
							for (int i = 0; i < 1 << dimnum; i++)
							{
								HistNodeND *child = (HistNodeND*)malloc(sizeof(HistNodeND));
								child->key = (node->key << dimnum) + i;
								child->cnum = node->cnum;
								child->height = node->height - 1;
								child->pnum = 0;
								bNodesC.push_back(child);
								cycle++;
							}
						}
						else
						{
							ranges.insert(make_pair(rangeL, rangeH));
						}
					}
					node->cnum = 0;

				}
				bNodes.swap(bNodesC);
				bNodesC.clear();
			}

		}
		else
		{
			for (int i = 0; i < bNodes.size(); i++)
			{
				ranges.insert(make_pair(bNodes[i]->key << bNodes[i]->height * dimnum, ((bNodes[i]->key + 1) << bNodes[i]->height * dimnum) - 1));
				//bNodes[i]->cnum = 0;
			}
		}
		
		auto end2 = chrono::high_resolution_clock::now();
		cout << "Adaptive decomposition costs: " << chrono::duration_cast<chrono::milliseconds>(end2 - end1).count() << "ms" << endl;
		measure.rangeNum = cycle;

		map <U, U> ranges_merged;
		U sfc_s = ranges.begin()->first;
		U sfc_e = ranges.begin()->second;
		for (auto it = ranges.begin(); it != ranges.end(); ++it)
		{
			if (it->first - sfc_e < 2) sfc_e = it->second;
			else
			{
				ranges_merged.insert(make_pair(sfc_s, sfc_e));
				sfc_s = it->first;
				sfc_e = it->second;
			}
		}
		ranges_merged.insert(make_pair(sfc_s, sfc_e));

		return ranges;
	}


public:
	Query(){
		measure = {};
		windowquery = {};
		PCDB = {};
	}

	Query(const PointCloudDB<T, U>& PC)
	{
		measure = {};
		windowquery = {};
		PCDB = PC;
	}

	virtual void QueryIOT(const NDWindow<T> & window)
	{
		windowquery = window;
		int dimnum = PCDB.nDims;
		const NDWindow<T> windowQ = window.Transform(PCDB.trans);   //transform original to cater to storage, to improve efficiency
		map <U, U> ranges;
		short dimbits = 20;  //maximum number of bits for a dimension retrieved from database

		if (PCDB.HIST)
		{
			auto start = chrono::high_resolution_clock::now();
			ranges = HistWindowRange(windowQ, dimbits);
			auto end = chrono::high_resolution_clock::now();
			measure.rangeComp = chrono::duration_cast<chrono::milliseconds>(end - start).count() - measure.histLoad;
		}
		else
		{
			auto start = chrono::high_resolution_clock::now();
			ranges = PlainWindowRange(windowQ, 11, dimbits);
			auto end = chrono::high_resolution_clock::now();
			measure.rangeComp = chrono::duration_cast<chrono::milliseconds>(end - start).count();
		}

		Environment *env = Environment::createEnvironment(Environment::DEFAULT);
		Connection  *con = env->createConnection(orclconn().User, orclconn().Password, orclconn().Database);

		auto start = chrono::high_resolution_clock::now();

		Statement *f_stmt = NULL;
		f_stmt = con->createStatement();
		string sql = "create table range_packs (lower number, upper number)";
		f_stmt->executeUpdate(sql);
		con->terminateStatement(f_stmt);

		Statement *stmt = con->createStatement();
		stmt->setSQL("insert into range_packs values(:l, :u)");
		stmt->setMaxIterations(ranges.size() + 1);
		stmt->setMaxParamSize(1, 30);
		stmt->setMaxParamSize(2, 30);
		for (auto it = ranges.begin(); it != ranges.end(); it++)
		{
			stringstream keystr;
			keystr << it->first;
			string key_l = keystr.str();
			stmt->setString(1, key_l);
			keystr.clear();
			keystr.str("");
			keystr << it->second;
			string key_h = keystr.str();
			stmt->setString(2, key_h);
			stmt->addIteration();
		}
		stmt->executeUpdate();

		sql = "select /*+ use_nl (t r)*/ t.sfc from " + PCDB.Table + " t, range_packs r where (t.sfc between r.lower and r.upper)";
		stmt->setSQL(sql);
		stmt->setPrefetchRowCount(FETCH_SIZE);
		ResultSet *rs = stmt->executeQuery();

		auto end1 = chrono::high_resolution_clock::now();
		measure.firstCost = measure.rangeComp + chrono::duration_cast<chrono::milliseconds>(end1 - start).count();

		ofstream output_file("E:/query_res_plain.csv");
		SFCConversion<int> sfc(dimnum, dimbits);
		unsigned long long apnum = 0;  //approximate number of point
		unsigned long long spnum = 0;  //accurate number of point
		NDPoint<int> pt;
		NDPoint<double> ptreal;
		while (rs->next())
		{
			apnum++;

			sfc_bigint key = (sfc_bigint)rs->getString(1);
			pt = sfc.MortonDecode(key);
			//output_file << pt[0] << ", " << pt[1] << ", " << pt[2] << ", " << pt[3] << "\n";
			if (Inside<int, T>(pt, windowQ) == 1)
			{
				spnum++;
				ptreal = pt.InverseTransform<double>(PCDB.trans);
				//output_file << setprecision(2) << ptreal[0] << ", " << ptreal[1] << ", " << ptreal[2] << "\n";
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

		cout << endl;
	}

	virtual void ExMeasurement(string filename)
	{
		ofstream output(filename, ios::app|ios::out);
		output << "Query table: " << PCDB.Table << ", query geometry: " << "[" << fixed << setprecision(2);
		for (int i = 0; i < windowquery.nDims; i++)
		{
			output << windowquery.minPoint[i]<<" ";
		}
		output << ", ";
		for (int i = 0; i < windowquery.nDims; i++)
		{
			output << windowquery.maxPoint[i] << " ";
		}
		output << "]\n";
		if (PCDB.HIST) output << "HistSFC\n";
		else output << "PlainSFC\n";

		output << "rangeNum, appPNum, accPNum, FPR, rangeComp, histLoad, firstCost, secondCost\n";
		output << measure.rangeNum << ", " << measure.appPNum << ", " << measure.accPNum << ", " << measure.FPR << ", "
			<< measure.rangeComp << ", " << measure.histLoad << ", " << measure.firstCost << ", " << measure.secondCost << "\n";
		output << "\n";
	}


};