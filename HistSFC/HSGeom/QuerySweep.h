#pragma once
#include "../Window.h"
#include "../Geom.h"
#include "../Query.h"
#include <vector>
#include <cmath> 

#define PREC (unsigned int)(1000)	//precison used for distance comparison

enum class ALG {
	SWEEP, SPHERE, VERTEX
};


typedef struct NodeNDQ
{
	NodeND * node;
	vector<int> *intersects;	//to record intersection planes
} NodeNDQ;

typedef struct HistNodeNDQ
{
	HistNodeND *histnode;
	vector<int> *intersects;
} HistNodeNDQ;


template <typename T, typename U>	//data type of coordinates and SFC key
class QueryExtendSweep : public Query <T, U> {
	using Query<T, U>::measure;
	using Query<T, U>::histroot;
	using Query<T, U>::PCDB;
private:
	NDGeom QueryGeom;
	ALG alg;
	unsigned int cycle_counter;

private:
	template<typename A>
	NDPoint<A> sweep_box_with_plane_enter(const halfspace& h, const NDWindow<A> & node)
	{
		NDPoint<A> enter(node.nDims);
		for (int i = 0; i < h.dimnum; i++)
		{
			if (h.w[i] < 0)
				enter[i] = node.maxPoint[i];
			else
				enter[i] = node.minPoint[i];
		}

		return enter;
	}

	template<typename A>
	NDPoint<A> sweep_box_with_plane_exit(const halfspace& h, const NDWindow<A> & node)
	{
		NDPoint<A> exit(node.nDims);
		for (int i = 0; i < h.dimnum; i++)
		{
			if (h.w[i] < 0)
				exit[i] = node.minPoint[i];
			else
				exit[i] = node.maxPoint[i];
		}

		return exit;
	}

	template<typename A>
	double signed_distance(const NDPoint<A> & enter, const halfspace& h)
	{
		double dist = 0;
		/*
		cout << "node: ";
		for (int i = 0; i < 4; i++)
		{
			cout << enter[i] << ", ";
		}
		cout << endl;

		for (int i = 0; i < h.dimnum; i++)
		{
			cout << h.w[i] << ", ";
		}
		cout << endl;
		*/
		for (int i = 0; i < h.dimnum; i++)
		{
			dist += enter[i] * h.w[i];
		}
		

		dist += -h.b;
		
		return dist;
	}

	template<typename A>
	short Intersect(const NDGeom & geom, const NDWindow<A> & node, vector<int> *vec_o, ALG al)
	{
		if (node.nDims != geom.dimnum)
		{
			throw("Dimensionality does not match!");
		}

		short res;
		double side_dist = node.maxPoint[0] - node.minPoint[0];
		double diagonal_dist = side_dist * pow(node.nDims, 0.5);
		int intersect_count = vec_o->size();
		vector<int> vec = *vec_o;

		switch (al)
		{
		case ALG::SPHERE:
		{
			res = 2; //first assume contain

			NDPoint<A> center(node.nDims);
			for (int i = 0; i < node.nDims; i++)
				center[i] = 0.5*(node.minPoint[i] + node.maxPoint[i]);

			for (auto it = vec.begin(); it != vec.end();)
			{
				cycle_counter++;	//cycles of iteration

				if (signed_distance(center, geom.faces[*it]) + diagonal_dist*0.5 < 0)
				{
					res = 0;
					break;
				}

				else {
					if (abs(signed_distance(center, geom.faces[*it])) < diagonal_dist*0.5)
					{
						res = 1;
						++it;
					}
					else
					{
						vec.erase(it);
						intersect_count--;	//contain
					}
				}
				
			}

			if (res == 1 and intersect_count!= vec_o->size())
			{
				vector<int> vec_res(intersect_count);
				for (int i = 0; i < vec_res.size(); i++)
					vec_res[i] = vec[i];
				*vec_o = vec_res;
			}
		}
			break;
		
		case ALG::SWEEP:
		{
			res = 1;	//intersect

			for (auto it = vec.begin(); it != vec.end();)
			{
				cycle_counter++;

				NDPoint<A> enter = sweep_box_with_plane_enter(geom.faces[*it], node);
				double enter_dist = signed_distance(enter, geom.faces[*it]);
				if (floor(enter_dist*PREC) >= 0)
				{
					vec.erase(it);
					intersect_count--;	//contain
				}
				else
				{
					double abs_enter_dist = abs(enter_dist);
					if (ceil(abs_enter_dist*PREC) > ceil(diagonal_dist*PREC))
					{
						res = 0;	//no-overlap
						break;
					}
					else
					{
						if (ceil(abs_enter_dist*PREC) > ceil(side_dist*PREC))
						{
							NDPoint<A> exit = sweep_box_with_plane_exit(geom.faces[*it], node);
							double exit_dist = signed_distance(exit, geom.faces[*it]);
							if (ceil(exit_dist*PREC) < 0)
							{
								res = 0;
								break;
							}
						}

					}

					++it;
				}

			}

			if (intersect_count)
			{
				if (res and intersect_count != vec_o->size())
				{
					vector<int> vec_res(intersect_count);
					for (int i = 0; i < vec_res.size(); i++)
						vec_res[i] = vec[i];
					*vec_o = vec_res;
				}
				
			}
			else
				res = 2;	//contain
		}
			break;

		case ALG::VERTEX:
		{
			res = 2;
			NDPoint<A> nodeVertex(node.nDims);
			for (auto it = vec.begin(); it != vec.end();)
			{
				int num_intersects = 0;   //how many vertices intersect a half-plane
				for (int i = 0; i < 1 << node.nDims; i++)
				{
					cycle_counter++;

					int bit_idx = i;
					for (int j = 0; j < node.nDims; j++)
					{
						if (bit_idx % 2)
							nodeVertex[j] = node.maxPoint[j];
						else
							nodeVertex[j] = node.minPoint[j];
						
						bit_idx = bit_idx >> 1;
					}
					
					double vdist = signed_distance(nodeVertex, geom.faces[*it]);
					num_intersects += (vdist >= 0);
				}

				if (num_intersects == 0)
				{
					res = 0;
					break;
				}	
				else {
					if (num_intersects == (1 << node.nDims))
					{
						vec.erase(it);
						intersect_count--;	//contain
					}
					else
					{
						res = 1;
						++it;
					}
						
				}

				
			}

			if (res == 1 and intersect_count != vec_o->size())
			{
				vector<int> vec_res(intersect_count);
				for (int i = 0; i < vec_res.size(); i++)
					vec_res[i] = vec[i];
				*vec_o = vec_res;
			}
		}

			break;
		}

		return res;

	}

	template<typename A>
	bool Inside(NDPoint<A> &pt, const NDGeom &geom)	//used for second filter
	{
		if (pt.returnSize() != geom.dimnum)
		{
			throw("Dimensionality does not match!");
		}
		else
		{
			short ncmp = 1;
			double expr = 0;
			for (auto it = geom.faces.begin(); it != geom.faces.end(); it++)
			{
				for (int i = 0; i < geom.dimnum; i++)
				{
					expr += it->w[i] * pt[i];
				}

				switch (it->s)
				{
				case Sign::le:
					ncmp &= expr <= it->b;
					break;
				case Sign::ge:
					ncmp &= expr >= it->b;
					break;
				case Sign::eq:
					ncmp &= expr == it->b;
					break;
				default:
					break;
				}
				expr = 0;
			}
			if (ncmp) return 1;
		}
		return 0;
	}


	map <U, U> PlainGeomRange(const NDGeom & geom, ALG al, short maxdepth, short dimbits = 20)
	{
		measure.histLoad = 0;
		map <U, U> ranges;
		int dimnum = PCDB.nDims;
		multimap <unsigned short, NodeNDQ> SearchT;
		short order = dimbits;
		NDPoint<int> NodeL(dimnum); //lower bound of child node
		NDPoint<int> NodeH(dimnum); //upper bound of child node
		SFCConversion<int> sfc(dimnum, order + 1); //with 1 more bit for the sake of rounding the double type

		if (!maxdepth) maxdepth = order + 1;

		NodeND rootnode = { 0, order };
		vector<int> init(geom.faces.size());
		for (int i = 0; i < init.size(); i++)
			init[i] = i;
		NodeNDQ rootnodeq = { &rootnode, &init };
		SearchT.insert(make_pair(order, rootnodeq));
		NDWindow<int> cell;
		NodeNDQ node;
		sfc_bigint rangeL;   //lower boundary of a range 
		sfc_bigint rangeH;   //upper boundary of a range 

		unsigned int cycle_time = MAX_CYCLE;  //threshold for the refinement cycles
		int cycle = SearchT.size();
		while (!SearchT.empty())
		{
			if (cycle <= cycle_time)
			{
				auto it = SearchT.rbegin();
				node = it->second;
				//free(SearchT.back());   //remove node built from malloc
				SearchT.erase(--it.base());
				for (int i = 0; i < 1 << dimnum; i++)
				{
					//for children
					sfc_bigint key = (node.node->key << dimnum) + i;
					unsigned short height = node.node->height - 1;
					NodeND* childnode = new NodeND{ key, height };
					vector<int> *intersectinfo = new vector<int>{ *(node.intersects) };
					NodeNDQ child = { childnode, intersectinfo };
					rangeL = key << (height*dimnum);
					rangeH = ((key + 1) << (height*dimnum)) - 1;
					NodeL = sfc.MortonDecode(rangeL);
					NodeH = sfc.MortonDecode(rangeH);
					cell.SetMinPoint(NodeL);
					cell.SetMaxPoint(NodeH);

					short mark = Intersect<int>(geom, cell, child.intersects, al);
					if (mark == 2)
					{
						ranges.insert(make_pair(rangeL, rangeH));
						cycle++;
					}


					else if (mark == 1)
					{
						cycle++;
						if (order - height < maxdepth)
						{
							SearchT.insert(make_pair(height, child));
						}

						else ranges.insert(make_pair(rangeL, rangeH));

					}

				}

				cycle--;
			}
			else break;

		}

		for (auto it = SearchT.begin(); it != SearchT.end(); it++)
		{
			node = it->second;
			ranges.insert(make_pair(node.node->key << (it->first*dimnum), ((node.node->key + 1) << (it->first*dimnum)) - 1));
		}

		measure.rangeNum = cycle;
		//cout << ranges.size() << endl;
		map <U, U> ranges_merged;
		U sfc_s = ranges.begin()->first;
		U sfc_e = ranges.begin()->second;
		for (auto it = next(ranges.begin(), 1); it != ranges.end(); ++it)
		{
			//cout << it->second - it->first << endl;
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


	map <U, U> HistGeomRange(const NDGeom & geom, ALG al, short dimbits = 20)
	{
		//Directly search the histogram tree
		auto start = chrono::high_resolution_clock::now();
		map <U, U> ranges;
		int dimnum = PCDB.nDims;
		vector <HistNodeNDQ *> bNodes;  //nodes located on the boundary
		vector <HistNodeNDQ *> SearchT;
		NDPoint<int> NodeL(dimnum);
		NDPoint<int> NodeH(dimnum);
		sfc_bigint rangeL;
		sfc_bigint rangeH;
		SFCConversion<int> sfc(dimnum, dimbits);
		NDWindow<int> cell;
		vector<int> init(geom.faces.size());
		for (int i = 0; i < init.size(); i++)
			init[i] = i;
		HistNodeNDQ rootQ = { histroot, &init };
		SearchT.push_back(&rootQ);
		HistNodeNDQ *node;
		int childnum = 1 << dimnum;
		long long sumBP = 0;  //points covered by boundary nodes
		long long sumIP = 0;  //points covered by inner nodes

		while (!SearchT.empty())
		{
			node = SearchT.back();
			SearchT.pop_back();
			rangeL = node->histnode->key << (node->histnode->height*dimnum);
			rangeH = (((node->histnode->key + 1) << (node->histnode->height*dimnum)) - 1);
			NodeL = sfc.MortonDecode(rangeL);
			NodeH = sfc.MortonDecode(rangeH);
			cell.SetMinPoint(NodeL);
			cell.SetMaxPoint(NodeH);

			short mark = Intersect<int>(geom, cell, node->intersects, al);
			if (mark == 2)
			{
				ranges.insert(make_pair(rangeL, rangeH));
				sumIP ++;
			}

			else if (mark == 1)
			{
				if (node->histnode->cnum != 0)
				{
					HistNodeND *curnode = node->histnode->child;
					for (int i = 0; i < node->histnode->cnum; i++)
					{
						vector<int> *intersectinfo = new vector<int>{ *(node->intersects) };
						HistNodeNDQ *child = new HistNodeNDQ{ curnode, intersectinfo };
						SearchT.push_back(child);
						curnode = curnode->neighbor;
					}
				}
				else
				{
					if (node->histnode->pnum <= childnum)
					{
						ranges.insert(make_pair(rangeL, rangeH));
						sumIP++;
					}
					else
					{
						bNodes.push_back(node);
						sumBP++;
					}
					//cout << "Space: " << NodeLR[0] << "," << NodeLR[1] << "," << NodeLR[2] << "," << NodeLR[3] << ";" << NodeHR[0] << "," << NodeHR[1] << "," << NodeHR[2] << "," << NodeHR[3] << ", number of points: " << node->pnum << ", ratio: "<< 1.0*node->pnum/(double)(keyH - keyL) << endl;
				}

			}

			else continue;
		}

		auto end1 = chrono::high_resolution_clock::now();
		cout << "Hist search costs: " << chrono::duration_cast<chrono::milliseconds>(end1 - start).count() << "ms" << endl;
		cout << "Inner points: " << sumIP << ", boundary points: " << sumBP << ", boundary nodes: " << bNodes.size() << endl;

		unsigned int cycle_time = MAX_CYCLE;  //threshold for the refinement cycles
		while (!bNodes.empty())
		{
			if (sumIP + sumBP <= cycle_time)
			{
				auto it = bNodes.begin();
				node = *it;
				bNodes.erase(it);
				for (int i = 0; i < childnum; i++)
				{
					HistNodeND *child = (HistNodeND*)malloc(sizeof(HistNodeND));
					child->key = (node->histnode->key << dimnum) + i;
					child->cnum = node->histnode->cnum;
					child->height = node->histnode->height - 1;
					child->pnum = node->histnode->pnum / childnum;
					rangeL = child->key << (child->height*dimnum);
					rangeH = (((child->key + 1) << (child->height*dimnum)) - 1);
					NodeL = sfc.MortonDecode(rangeL);
					NodeH = sfc.MortonDecode(rangeH);
					cell.SetMinPoint(NodeL);
					cell.SetMaxPoint(NodeH);
				
					short mark = Intersect<int>(geom, cell, node->intersects, al);
					if (mark == 2)
					{
						ranges.insert(make_pair(rangeL, rangeH));
						sumIP++;
					}
					else if (mark == 1)
					{
						if (child->pnum <= childnum)
						{
							ranges.insert(make_pair(rangeL, rangeH));
							sumIP++;
						}
						else
						{
							vector<int> *intersectinfo = new vector<int>{ *(node->intersects) };
							HistNodeNDQ *childnode = new HistNodeNDQ{ child, intersectinfo };
							bNodes.push_back(childnode);
							sumBP++;
						}
						
					}
				}

				sumBP--;
			}

			else break;
		}

		for (auto it = bNodes.begin(); it != bNodes.end(); it++)
		{
			ranges.insert(make_pair((*it)->histnode->key << (*it)->histnode->height * dimnum, (((*it)->histnode->key + 1) << (*it)->histnode->height * dimnum) - 1));
		}

		auto end2 = chrono::high_resolution_clock::now();
		cout << "Adaptive decomposition costs: " << chrono::duration_cast<chrono::milliseconds>(end2 - end1).count() << "ms" << endl;
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
	QueryExtendSweep() {
		measure = {};
		QueryGeom = {};
		PCDB = {};
	}

	QueryExtendSweep(const PointCloudDB<T, U>& PC)
	{
		measure = {};
		QueryGeom = {};
		PCDB = PC;
		cycle_counter = 0;

		if (PCDB.HIST)
		{
			auto start = chrono::high_resolution_clock::now();
			histroot = this->HistLoad(PCDB.HistTab);
			auto end = chrono::high_resolution_clock::now();
			measure.histLoad = chrono::duration_cast<chrono::milliseconds>(end - start).count();
		}
	}
	
	void QueryIOT(const NDGeom & geom, ALG al)
	{
		alg = al;
		cycle_counter = 0;
		QueryGeom = geom;
		int dimnum = PCDB.nDims;
		NDGeom geomtrans = geom.Transform(PCDB.trans);
		NDGeom geomtrans2 = geom.Transform(PCDB.trans).Normalize();

		map <U, U> ranges;
		short dimbits = 12;  //maximum number of bits for a dimension retrieved from database
		
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
			ranges = HistGeomRange(geomtrans2, al, dimbits);
			auto end = chrono::high_resolution_clock::now();
			measure.rangeComp = chrono::duration_cast<chrono::milliseconds>(end - start).count();
		}
		else
		{
			auto start = chrono::high_resolution_clock::now();
			ranges = PlainGeomRange(geomtrans2, al, 10, dimbits);	//search depth can be modified depending on accuracy requirement
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

		ofstream output_file("E:/querysweep_res_plain.csv");
		SFCConversion<int> sfc(dimnum, dimbits);
		unsigned long long apnum = 0;  //approximate number of point
		unsigned long long spnum = 0;  //accurate number of point
		NDPoint<int> pt;
		NDPoint<double> ptreal;
		while (rs->next())
		{
			apnum++;

			string s = rs->getString(1);
			if (varchar_mode) s.erase(std::remove_if(s.begin(), s.end(), ::isspace), s.end());
			sfc_bigint key = (sfc_bigint)s;
			pt = sfc.MortonDecode(key);
			//output_file << pt[0] << ", " << pt[1] << ", " << pt[2] << ", " << pt[3] << "\n";
			if (Inside<int>(pt, geomtrans))
			{
				spnum++;
				//output_file << pt[0] << ", " << pt[1] << endl;
				//ptreal = pt.InverseTransform<double>(PCDB.trans);
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
	}

	void ExMeasurement(string filename) override
	{
		ofstream output(filename, ios::app | ios::out);
		output << "Query table: " << PCDB.Table << ", query geometry: " << "[halfspaces]\n";
		output << "Method: ";
		switch (alg)
		{
		case ALG::SWEEP:
			output << "SWEEP";
			break;
		case ALG::SPHERE:
			output << "SPHERE";
			break;
		case ALG::VERTEX:
			output << "VERTEX";
			break;
		default:
			break;
		}
		output << endl;
		if (PCDB.HIST) output << "HistSFC\n";
		else output << "PlainSFC\n";

		output << "rangeNum, appPNum, accPNum, FPR, cycleCount, rangeComp, histLoad, firstCost, secondCost\n";
		output << measure.rangeNum << ", " << measure.appPNum << ", " << measure.accPNum << ", " << measure.FPR << ", " << cycle_counter
			<< ", " << measure.rangeComp << ", " << measure.histLoad << ", " << measure.firstCost << ", " << measure.secondCost << "\n";
		output << "\n";
	}


};
