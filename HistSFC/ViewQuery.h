/*4D perspective view (x, y, z, cLoD) selection based on 4D geometrical computation*/

#pragma once
#include "ViewStruct.h"
#include "Query.h"

using namespace viewlib;

template <typename T, typename U>       //data type of coordinates and SFC key
class ViewQuery : public Query <T, U> {
	using Query<T, U>::measure;
	using Query<T, U>::histroot;
	using Query<T, U>::PCDB;

private:
	map <U, U> PlainViewRange(const Frustum3D & vf, Cone4D & cn, short maxdepth, short dimbits = 20)
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
					auto cellreal = cell.InverseTransform<T>(PCDB.trans);
					if (intersectFV3D<T>(vf, cellreal) == 2 and intersectCN4D<T>(cn, cellreal) == 2)
					{
						ranges.insert(make_pair(rangeL, rangeH));
						cycle++;
					}

					else if (intersectFV3D<T>(vf, cellreal) and intersectCN4D<T>(cn, cellreal))
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
		//cout << "range number: " << ranges.size() << ", miss: " << miss << endl;
		map <U, U> ranges_merged;
		U sfc_s = ranges.begin()->first;
		U sfc_e = ranges.begin()->second;
		for (auto it = next(ranges.begin(), 1); it != ranges.end(); ++it)
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
		return ranges_merged;
	}

	map <U, U> HistViewRange(const Frustum3D & vf, Cone4D & cn, short dimbits = 20)
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
			auto cellreal = cell.InverseTransform<T>(PCDB.trans);
			if (intersectFV3D<T>(vf, cellreal) == 2 and intersectCN4D<T>(cn, cellreal) == 2)
			{
				ranges.insert(make_pair(rangeL, rangeH));
				sumIP++;
			}

			else if (intersectFV3D<T>(vf, cellreal) and intersectCN4D<T>(cn, cellreal))
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

		unsigned int cycle_time = 2;  //threshold for the refinement cycles
		unsigned int cycle = bNodes.size();
		while (!bNodes.empty())
		{
			if (cycle <= cycle_time)
			{
				node = bNodes[0];
				bNodes.erase(bNodes.begin());
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
					auto cellreal = cell.InverseTransform<T>(PCDB.trans);
					if (intersectFV3D<T>(vf, cellreal) == 2 and intersectCN4D<T>(cn, cellreal) == 2)
					{
						ranges.insert(make_pair(rangeL, rangeH));
						cycle++;
					}

					else if (intersectFV3D<T>(vf, cellreal) and intersectCN4D<T>(cn, cellreal))
					{
						cycle++;
						bNodes.push_back(child);
					}

				}
				cycle--;
			}

			else break;

		}

		for (int i = 0; i < bNodes.size(); i++)
		{
			ranges.insert(make_pair(bNodes[i]->key << (bNodes[i]->height*dimnum), ((bNodes[i]->key + 1) << (bNodes[i]->height*dimnum)) - 1));
		}

		auto end2 = chrono::high_resolution_clock::now();
		//cout << "Adaptive decomposition costs: " << chrono::duration_cast<chrono::milliseconds>(end2 - end1).count() << "ms" << endl;
		measure.rangeNum = ranges.size();

		map <U, U> ranges_merged;
		U sfc_s = ranges.begin()->first;
		U sfc_e = ranges.begin()->second;
		for (auto it = next(ranges.begin(), 1); it != ranges.end(); ++it)
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

		return ranges_merged;
	}

public:
	ViewQuery() {
		measure = {};
		histroot = nullptr;
		PCDB = {};
	}

	ViewQuery(const PointCloudDB<T, U>& PC)
	{
		measure = {};
		PCDB = PC;
		histroot = nullptr;
		if (PCDB.HIST)
		{
			auto start = chrono::high_resolution_clock::now();
			histroot = this->HistLoad(PCDB.HistTab);
			auto end = chrono::high_resolution_clock::now();
			measure.histLoad = chrono::duration_cast<chrono::milliseconds>(end - start).count();
		}
	}

	void QueryIOT(const viewPos& qbox, T mlod)
	{
		int dimnum = PCDB.nDims;
		map <U, U> ranges;
		short dimbits = 12;

		auto v3d = FrustumBuild(qbox);
		Cone4D LoDview = { {qbox.P[0],qbox.P[1],qbox.P[2]}, qbox.distance, qbox.distance / mlod, mlod };

		if (PCDB.HIST)
		{
			if (!histroot)
			{
				auto start = chrono::high_resolution_clock::now();
				histroot = this->HistLoad(PCDB.HistTab);
				auto end = chrono::high_resolution_clock::now();
				measure.histLoad = chrono::duration_cast<chrono::milliseconds>(end - start).count();
			}
			auto start = chrono::high_resolution_clock::now();
			ranges = HistViewRange(v3d, LoDview, dimbits);
			auto end = chrono::high_resolution_clock::now();
			measure.rangeComp = chrono::duration_cast<chrono::milliseconds>(end - start).count();
		}
		else
		{
			auto start = chrono::high_resolution_clock::now();
			ranges = PlainViewRange(v3d, LoDview, 12, dimbits);	//search depth can be modified depending on accuracy requirement
			auto end = chrono::high_resolution_clock::now();
			measure.rangeComp = chrono::duration_cast<chrono::milliseconds>(end - start).count();
		}

		Environment *env = Environment::createEnvironment(Environment::DEFAULT);
		Connection  *con = env->createConnection(orclconn().User, orclconn().Password, orclconn().Database);

		auto start = chrono::high_resolution_clock::now();

		bool varchar_mode = false;	//open this mode when key is larger than the number type of Oracle
		if ((double)ranges.rbegin()->second > pow(10, 38)) varchar_mode = true;

		Statement *f_stmt = NULL;
		f_stmt = con->createStatement();
		string sql = "create table range_packs (lower number, upper number)";
		if (varchar_mode) sql = "create table range_packs (lower VARCHAR2(60), upper VARCHAR2(60))";
		f_stmt->executeUpdate(sql);
		con->terminateStatement(f_stmt);

		Statement *stmt = con->createStatement();
		stmt->setSQL("insert into range_packs values(:l, :u)");
		if (varchar_mode) stmt->setSQL("insert into range_packs values(Lpad(:l, 60, ' '), Lpad(:u, 60, ' '))");
		stmt->setMaxIterations(ranges.size() + 1);
		stmt->setMaxParamSize(1, 100);
		stmt->setMaxParamSize(2, 100);
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

		ofstream output_file("E:/query_4dgeom_apx.csv");
		SFCConversion<int> sfc(dimnum, dimbits);
		unsigned long long apnum = 0;  //approximate number of point
		unsigned long long spnum = 0;  //accurate number of point
		NDPoint<int> pt;
		NDPoint<T> ptreal;
		while (rs->next())
		{
			apnum++;

			string s = rs->getString(1);
			if (varchar_mode) s.erase(std::remove_if(s.begin(), s.end(), ::isspace), s.end());
			sfc_bigint key = (sfc_bigint)s;
			pt = sfc.MortonDecode(key);
			ptreal = pt.InverseTransform<T>(PCDB.trans);
			output_file << pt[0] << ", " << pt[1] << ", " << pt[2] << ", " << pt[3] << "\n";
			if (intersectFPP3D<T>(v3d, ptreal) and intersectCPP4D<T>(LoDview, ptreal))
			{
				spnum++;
				//output_file << ptreal[0] << ", " << ptreal[1] << ", " << ptreal[2] << ", " << ptreal[3] << "\n";
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
