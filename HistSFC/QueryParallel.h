/*Balanced parallel decoding*/
#pragma once
#include <map>
#include <deque>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <cctype>
#include <omp.h>
#include "Window.h"
#include "PointCloud.h"
#include "Oracle.h"
#include "BaseStruct.h"
#include "SFCConversion.h"

#define MAX_CYCLE (unsigned int) 1000000

template <typename T, typename U>   //type for coordinate and database key, respectively
class Query {
private:
	NDWindow<T> windowquery;        //query window
	long long miss = 0;
protected:
	HistNodeND *histroot;   //root of HistTree
	Measurement measure;    //collecting performance info
public:
	PointCloudDB<T, U> PCDB;

protected:
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

private:
	template<typename A, typename B>
	bool Inside(NDPoint<A> &pt, const NDWindow<B> &wd)      //used for second filter
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

	template<typename A, typename B>
	short Intersect(const NDWindow<A> & queryrect, const NDWindow<B> & node)
	{
		short ncmp = 1;
		short dimnum = queryrect.nDims;
		for (int i = 0; i < dimnum; i++)
		{
			ncmp &= node.minPoint[i] <= queryrect.maxPoint[i] && node.maxPoint[i] >= queryrect.minPoint[i];
		}
		if (ncmp)
		{
			for (int i = 0; i < dimnum; i++)
			{
				ncmp &= node.minPoint[i] >= queryrect.minPoint[i] && node.maxPoint[i] <= queryrect.maxPoint[i];
			}
			if (ncmp)
				return 1;   //contain
			else return 2; //overlap
		}

		//not overlap
		miss++;
		return 0;
	}

	template<typename A, typename B>
	double IntersectR(const NDWindow<A> & queryrect, const NDWindow<B> & node)
	{
		double ratio = 1;
		for (int i = 0; i < node.nDims; i++)
		{
			ratio *= (min<double>(queryrect.maxPoint[i], node.maxPoint[i]) - max<double>(queryrect.minPoint[i], node.minPoint[i])) / (node.maxPoint[i] - node.minPoint[i]);
		}
		return 1 - ratio;
	}

	multimap <U, U> PlainWindowRange(const NDWindow<double> & window, short maxdepth, unsigned short dimbits = 20)
	{
		measure.histLoad = 0;
		map <U, U> ranges;
		int dimnum = PCDB.nDims;
		multimap <unsigned short, sfc_bigint> SearchT;
		unsigned short order = dimbits;
		NDPoint<int> NodeL(dimnum); //lower bound of child node
		NDPoint<int> NodeH(dimnum); //upper bound of child node
		SFCConversion<int> sfc(dimnum, order + 1); //with 1 more bit for the sake of rounding the double type

		if (!maxdepth) maxdepth = order + 1;

		SearchT.insert(make_pair(order, 0));  //insert root node
		NDWindow<int> cell;
		NodeND node;
		sfc_bigint rangeL;   //lower boundary of a range
		sfc_bigint rangeH;   //upper boundary of a range

		unsigned int cycle_time = MAX_CYCLE;  //threshold for the refinement cycles
		int cycle = SearchT.size();
		while (!SearchT.empty())
		{
			if (cycle <= cycle_time)
			{
				auto it = SearchT.rbegin();
				node = { it->second, it->first };
				//free(SearchT.back());   //remove node built from malloc
				SearchT.erase(--it.base());
				for (int i = 0; i < 1 << dimnum; i++)
				{
					//for children
					sfc_bigint key = (node.key << dimnum) + i;
					unsigned short height = node.height - 1;
					rangeL = key << (height*dimnum);
					rangeH = ((key + 1) << (height*dimnum)) - 1;
					NodeL = sfc.MortonDecode(rangeL);
					NodeH = sfc.MortonDecode(rangeH);
					cell.SetMinPoint(NodeL);
					cell.SetMaxPoint(NodeH);
					if (Intersect<double, int>(window, cell) == 1)
					{
						ranges.insert(make_pair(rangeL, rangeH));
						cycle++;
					}

					else if (Intersect<double, int>(window, cell) > 1)
					{
						cycle++;
						if (order - height < maxdepth)
						{
							SearchT.insert(make_pair(height, key));
						}
						else
						{
							ranges.insert(make_pair(rangeL, rangeH));
						}

					}

				}
				cycle--;
			}

			else break;

		}

		for (auto it = SearchT.begin(); it != SearchT.end(); it++)
		{
			ranges.insert(make_pair(it->second << (it->first*dimnum), ((it->second + 1) << (it->first*dimnum)) - 1));
		}

		measure.rangeNum = ranges.size();
		cout << "range number: " << ranges.size() << ", miss: " << miss << endl;

		multimap <U, U> ranges_merged;
		if (ranges.size())
		{
			U sfc_s = ranges.begin()->first;
			U sfc_e = ranges.begin()->second;
			for (auto it = next(ranges.begin(), 1); it != ranges.end(); ++it)
			{
				if (it->first - sfc_e < 2) sfc_e = it->second;
				else
				{
					ranges_merged.insert(make_pair(sfc_e - sfc_s, sfc_s));
					sfc_s = it->first;
					sfc_e = it->second;
				}
			}
			ranges_merged.insert(make_pair(sfc_e - sfc_s, sfc_s));
		}
		return ranges_merged;
	}

	multimap <U, U> HistWindowRange(const NDWindow<double>& window, short dimbits = 20)
	{
		//Directly search the histogram tree
		auto start = chrono::high_resolution_clock::now();
		map <U, U> ranges;
		int dimnum = PCDB.nDims;
		deque <HistNodeND *> bNodes;  //nodes located on the boundary
		vector <HistNodeND *> SearchT;
		NDPoint<int> NodeL(dimnum);
		NDPoint<int> NodeH(dimnum);
		sfc_bigint rangeL;
		sfc_bigint rangeH;
		SFCConversion<int> sfc(dimnum, dimbits);
		NDWindow<int> cell;
		SearchT.push_back(histroot);
		HistNodeND *node;
		int childnum = 1 << dimnum;
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

			if (Intersect<double, int>(window, cell) == 1)
			{
				ranges.insert(make_pair(rangeL, rangeH));
				sumIP++;
			}

			else if (Intersect<double, int>(window, cell) > 1)
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
					if (node->pnum <= childnum)
					{
						ranges.insert(make_pair(rangeL, rangeH));
						sumIP++;
					}
					else
					{
						bNodes.push_back(node);
						sumBP++;
					}
				}

			}

			else continue;
		}

		auto end1 = chrono::high_resolution_clock::now();

		unsigned int cycle_time = MAX_CYCLE;  //threshold for the refinement cycles
		while (!bNodes.empty())
		{
			if (sumIP + sumBP <= cycle_time)
			{
				node = bNodes[0];
				bNodes.pop_front();
				for (int i = 0; i < childnum; i++)
				{
					//for children
					HistNodeND *child = (HistNodeND*)malloc(sizeof(HistNodeND));
					child->key = (node->key << dimnum) + i;
					child->height = node->height - 1;
					rangeL = child->key << (child->height*dimnum);
					rangeH = ((child->key + 1) << (child->height*dimnum)) - 1;
					NodeL = sfc.MortonDecode(rangeL);
					NodeH = sfc.MortonDecode(rangeH);
					cell.SetMinPoint(NodeL);
					cell.SetMaxPoint(NodeH);
					if (Intersect<double, int>(window, cell) == 1)
					{
						ranges.insert(make_pair(rangeL, rangeH));
						sumIP++;
					}

					else if (Intersect<double, int>(window, cell) > 1)
					{
						sumBP++;
						bNodes.push_back(child);
					}

				}
				sumBP--;
			}

			else break;

		}

		for (int i = 0; i < bNodes.size(); i++)
		{
			ranges.insert(make_pair(bNodes[i]->key << (bNodes[i]->height*dimnum), ((bNodes[i]->key + 1) << (bNodes[i]->height*dimnum)) - 1));
		}

		auto end2 = chrono::high_resolution_clock::now();
		cout << "Adaptive decomposition costs: " << chrono::duration_cast<chrono::milliseconds>(end2 - end1).count() << "ms" << endl;

		measure.rangeNum = ranges.size();

		multimap <U, U> ranges_merged;
		if (ranges.size())
		{
			U sfc_s = ranges.begin()->first;
			U sfc_e = ranges.begin()->second;
			for (auto it = next(ranges.begin(), 1); it != ranges.end(); ++it)
			{
				if (it->first - sfc_e < 2) sfc_e = it->second;
				else
				{
					ranges_merged.insert(make_pair(sfc_e - sfc_s, sfc_s));
					sfc_s = it->first;
					sfc_e = it->second;
				}
			}
			ranges_merged.insert(make_pair(sfc_e - sfc_s, sfc_s));
		}

		return ranges_merged;
	}

	multimap <U, U> HistWindowRange2(const NDWindow<double>& window, short dimbits = 20)
	{
		//Directly search the histogram tree
		auto start = chrono::high_resolution_clock::now();
		map <U, U> ranges;
		int dimnum = PCDB.nDims;
		multimap <double, HistNodeND *> bNodes;  //intersect ration times pnum, nodes located on the boundary
		vector <HistNodeND *> SearchT;
		NDPoint<int> NodeL(dimnum);
		NDPoint<int> NodeH(dimnum);
		sfc_bigint rangeL;
		sfc_bigint rangeH;
		SFCConversion<int> sfc(dimnum, dimbits);
		NDWindow<int> cell;
		SearchT.push_back(histroot);
		HistNodeND* node;
		int childnum = 1 << dimnum;
		long long sumBP = 0;  //number of boundary nodes
		long long sumIP = 0;  //number of inner nodes

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

			if (Intersect<double, int>(window, cell) == 1)
			{
				ranges.insert(make_pair(rangeL, rangeH));
				sumIP++;
			}

			else if (Intersect<double, int>(window, cell) > 1)
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
					double ratio = IntersectR<double, int>(window, cell);
					bNodes.insert(make_pair(ratio*node->pnum, node));
					sumBP++;
				}

			}

			else continue;

		}

		auto end1 = chrono::high_resolution_clock::now();
		//cout << "Hist search costs: " << chrono::duration_cast<chrono::milliseconds>(end1 - start).count() << "ms" << endl;
		cout << "Hist leaf nodes: " << bNodes.size() << endl;

		unsigned int cycle_time = MAX_CYCLE;  //threshold for the refinement cycle
		while (!bNodes.empty())
		{
			if (sumIP + sumBP <= cycle_time)
			{
				auto it = bNodes.rbegin();
				node = it->second;
				bNodes.erase(--it.base());
				for (int i = 0; i < childnum; i++)
				{
					HistNodeND *child = (HistNodeND*)malloc(sizeof(HistNodeND));
					child->key = (node->key << dimnum) + i;
					child->cnum = childnum;
					child->height = node->height - 1;
					child->pnum = node->pnum / childnum;
					rangeL = child->key << (child->height*dimnum);
					rangeH = (((child->key + 1) << (child->height*dimnum)) - 1);
					NodeL = sfc.MortonDecode(rangeL);
					NodeH = sfc.MortonDecode(rangeH);
					cell.SetMinPoint(NodeL);
					cell.SetMaxPoint(NodeH);
					if (Intersect<double, int>(window, cell) == 1)
					{
						ranges.insert(make_pair(rangeL, rangeH));
						sumIP++;
					}
					else if (Intersect<double, int>(window, cell) > 1)
					{
						double ratio = IntersectR<double, int>(window, cell);
						bNodes.insert(make_pair(ratio*child->pnum, child));
						sumBP++;
					}
				}
				sumBP--;
			}
			else break;

		}

		for (auto it = bNodes.begin(); it != bNodes.end(); it++)
		{
			ranges.insert(make_pair(it->second->key << it->second->height * dimnum, ((it->second->key + 1) << it->second->height * dimnum) - 1));
		}

		auto end2 = chrono::high_resolution_clock::now();
		measure.rangeNum = ranges.size();
		cout << "range number: " << ranges.size() << endl;

		multimap <U, U> ranges_merged;
		if (ranges.size())
		{
			U sfc_s = ranges.begin()->first;
			U sfc_e = ranges.begin()->second;
			for (auto it = next(ranges.begin(), 1); it != ranges.end(); ++it)
			{
				if (it->first - sfc_e < 2) sfc_e = it->second;
				else
				{
					ranges_merged.insert(make_pair(sfc_e - sfc_s, sfc_s));
					sfc_s = it->first;
					sfc_e = it->second;
				}
			}
			ranges_merged.insert(make_pair(sfc_e - sfc_s, sfc_s));
		}

		return ranges_merged;
	}



public:
	Query() {
		measure = {};
		windowquery = {};
		PCDB = {};
	}

	Query(const PointCloudDB<T, U>& PC)
	{
		measure = {};
		windowquery = {};
		PCDB = PC;
		histroot = nullptr;
		if (PCDB.HIST)
		{
			auto start = chrono::high_resolution_clock::now();
			histroot = HistLoad(PCDB.HistTab);
			auto end = chrono::high_resolution_clock::now();
			measure.histLoad = chrono::duration_cast<chrono::milliseconds>(end - start).count();
		}
	}

	virtual void QueryIOT(const NDWindow<T> & window)
	{
		windowquery = window;
		int dimnum = PCDB.nDims;
		NDWindow<double> windowQ = window.template Transform<double>(PCDB.trans);   //transform original to cater to storage, to improve efficiency
		multimap <U, U> ranges;
		short dimbits = 30;  //maximum number of bits for a dimension retrieved from database

		if (PCDB.HIST)
		{
			if (!histroot)
			{
				auto start = chrono::high_resolution_clock::now();
				histroot = HistLoad(PCDB.HistTab);
				auto end = chrono::high_resolution_clock::now();
				measure.histLoad = chrono::duration_cast<chrono::milliseconds>(end - start).count();
			}
			auto start = chrono::high_resolution_clock::now();
			ranges = HistWindowRange(windowQ, dimbits);
			auto end = chrono::high_resolution_clock::now();
			measure.rangeComp = chrono::duration_cast<chrono::milliseconds>(end - start).count();
		}
		else
		{
			auto start = chrono::high_resolution_clock::now();
			ranges = PlainWindowRange(windowQ, 25, dimbits);       //search depth can be modified depending on accuracy requirement
			auto end = chrono::high_resolution_clock::now();
			measure.rangeComp = chrono::duration_cast<chrono::milliseconds>(end - start).count();

		}

		if (ranges.size()) {
			Environment *env = Environment::createEnvironment(Environment::DEFAULT);
			Connection  *con = env->createConnection(orclconn().User, orclconn().Password, orclconn().Database);

			auto start = chrono::high_resolution_clock::now();

			bool varchar_mode = false;      //open this mode when key is larger than the number type of Oracle
			if ((double)ranges.rbegin()->second > pow(10, 38)) varchar_mode = true;

			Statement *f_stmt = NULL;
			f_stmt = con->createStatement();
			string sql = "create table range_packs (id integer, lower number, upper number)";
			if (varchar_mode) sql = "create table range_packs (id integer, lower VARCHAR2(60), upper VARCHAR2(60))";
			f_stmt->executeUpdate(sql);
			con->terminateStatement(f_stmt);

			Statement *stmt = con->createStatement();
			stmt->setSQL("insert into range_packs values(:i, :l, :u)");
			if (varchar_mode) stmt->setSQL("insert into range_packs values(:i, Lpad(:l, 60, ' '), Lpad(:u, 60, ' '))");
			stmt->setMaxIterations(ranges.size() + 1);
			stmt->setMaxParamSize(2, 100);
			stmt->setMaxParamSize(3, 100);
			int num_threads = omp_get_max_threads();

			for (auto it = ranges.begin(); it != ranges.end(); it++)
			{
				stmt->setInt(1, rand() % num_threads);
				stringstream keystr;
				keystr << it->second;
				string key_l = keystr.str();
				stmt->setString(2, key_l);
				keystr.clear();
				keystr.str("");
				keystr << it->second + it->first;
				string key_h = keystr.str();
				stmt->setString(3, key_h);
				stmt->addIteration();
			}
			stmt->executeUpdate();
			con->commit();

			long long apnum = 0;
			long long spnum = 0;
			long long ap, sp;

			#pragma omp parallel firstprivate(ap,sp)
			{
				int pid = omp_get_thread_num();
				Environment *envp = Environment::createEnvironment(Environment::THREADED_UNMUTEXED);
				Connection  *conp = envp->createConnection(orclconn().User, orclconn().Password, orclconn().Database);
				Statement *pstmt = conp->createStatement();
				pstmt->setPrefetchRowCount(FETCH_SIZE);
				string sql_p = "select /*+ use_nl (t r)*/ t.sfc from " + PCDB.Table + " t, range_packs r where r.id =" + to_string(pid) + " and (t.sfc between r.lower and r.upper)";
				pstmt->setSQL(sql_p);
				ResultSet *rs = pstmt->executeQuery();
				auto start2 = chrono::high_resolution_clock::now();

				//ofstream output_file("/home/haicheng/pcdata/query_res_app.csv");
				SFCConversion<int> sfc(dimnum, dimbits);
				while (rs->next())
				{
					ap++;

					string s = rs->getString(1);
					if (varchar_mode) s.erase(std::remove_if(s.begin(), s.end(), ::isspace), s.end());
					sfc_bigint key = (sfc_bigint)s;
					NDPoint<int> pt = sfc.MortonDecode(key);
					//output_file << pt[0] << ", " << pt[1] << ", " << pt[2] << ", " << pt[3] << "\n";
					if (Inside<int>(pt, windowQ))
					{
						sp++;
						NDPoint<double> ptreal = pt.InverseTransform<double>(PCDB.trans);
						//output_file << fixed << setprecision(2) << ptreal[1] + x_d<< ", " << ptreal[2] + y_d << ", " << ptreal[3] << "\n";
					}
				}

				auto end2 = chrono::high_resolution_clock::now();
				//cout << "pid: "<<pid<<" , points: " <<ap<<endl;
				pstmt->closeResultSet(rs);
				conp->terminateStatement(pstmt);
				envp->terminateConnection(conp);
				Environment::terminateEnvironment(envp);
				#pragma omp critical
				{
					apnum += ap;
					spnum += sp;
					measure.firstCost = measure.rangeComp + chrono::duration_cast<chrono::milliseconds>(start2 - start).count();
					measure.secondCost = chrono::duration_cast<chrono::milliseconds>(end2 - start2).count();
				}
			}

			measure.appPNum = apnum;
			measure.accPNum = spnum;
			measure.FPR = (apnum - spnum)*1.0f / spnum;

			stmt->executeUpdate("drop table range_packs");
			con->commit();
			con->terminateStatement(stmt);
			env->terminateConnection(con);
			Environment::terminateEnvironment(env);


		}
	}

	virtual void ExMeasurement(string filename)
	{
		if (measure.rangeNum) {
			ofstream output(filename, ios::app | ios::out);
			output << "Query table: " << PCDB.Table << ", query geometry: " << "[" << fixed << setprecision(2);
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
			if (PCDB.HIST) output << "HistSFC\n";
			else output << "PlainSFC\n";

			output << "rangeNum, appPNum, accPNum, FPR, rangeComp, histLoad, firstCost, secondCost\n";
			output << measure.rangeNum << ", " << measure.appPNum << ", " << measure.accPNum << ", " << measure.FPR << ", "
				<< measure.rangeComp << ", " << measure.histLoad << ", " << measure.firstCost << ", " << measure.secondCost << "\n";
			output << "\n";
		}
	}

	virtual void ExMeasurement_batch(string filename)
	{

		ofstream output(filename, ios::app | ios::out);
		output << measure.rangeNum << ", " << measure.appPNum << ", " << measure.accPNum << ", " << measure.FPR << ", "
			<< measure.rangeComp << ", " << measure.histLoad << ", " << measure.firstCost << ", " << measure.secondCost << "\n";
		//output << measure.FPR << ",";

	}
};
