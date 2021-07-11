/*nD-polytope querying based on CPLEX, i.e., the linear optimization method*/ 

#pragma once
#include "../Window.h"
#include "../Geom.h"
#include "../Query.h"
#include <ilcplex/ilocplex.h>
ILOSTLBEGIN

double cptime = 0;

//convert sfc cell to variable ranges in CPLEX
template< typename T>
IloNumVarArray cell2var(IloEnv env, const NDWindow<T> &cell)
{
	IloNumVarArray x(env);
	for (int i = 0; i < cell.nDims; i++)
		x.add(IloNumVar(env, cell.minPoint[i], cell.maxPoint[i]));

	return x;
}

//convert halfspace to range in cplex
IloRange halfspace2con(IloEnv env, const halfspace &face, IloNumVarArray x)
{
	if (face.dimnum != x.getSize()) throw "size does not match!";

	IloRange con;
	switch (face.s)
	{
	case Sign::le:
		con = IloRange(env, -IloInfinity, face.b);
		break;
	case Sign::ge:
		con = IloRange(env, face.b, IloInfinity);
		break;
	case Sign::eq:
		con = IloRange(env, face.b, face.b);
		break;
	default:
		break;
	}

	for (int i = 0; i < face.dimnum; i++)
		con.setLinearCoef(x[i], face.w[i]);

	return con;
}

//convert to range constraint in cplex
IloRangeArray & geom2cons(IloEnv env, const NDGeom &geom, IloNumVarArray x)
{
	if (geom.dimnum != x.getSize()) throw "size does not match!";

	IloRange con;
	IloRangeArray cons(env);
	for (auto it = geom.faces.begin(); it != geom.faces.end(); it++)
	{
		con = halfspace2con(env, *it, x);
		cons.add(con);
	}

	/*specifically for view frustum*/
	for (auto it = geom.curves.begin(); it != geom.curves.end(); it++)
	{
		//cons.add(-1* x[0] * x[0] <= -100);
		cons.add(it->w[0] * x[0] * x[0] + it->w[1] * x[0] + it->w[2] * x[1] * x[1] + it->w[3] * x[1] + it->w[4] * x[2] * x[2] + it->w[5] * x[2] + 244.14 * x[3] <= it->b);
	}

	return cons;

}

template <typename T, typename U>       //data type of coordinates and SFC key
class QueryExtendCP : public Query <T, U> {
	using Query<T, U>::measure;
	using Query<T, U>::histroot;
	using Query<T, U>::PCDB;
private:
	NDGeom QueryGeom;

private:
	template<typename A>
	bool Inside(NDPoint<A> &pt, const NDGeom &geom) //used for second filter
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
				if (ncmp) {
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
				else return 0;
			}

			for (auto it = geom.curves.begin(); it != geom.curves.end(); it++)
			{
				if (ncmp) {
					for (int i = 0; i < geom.dimnum; i++)
					{
						expr += it->w[2 * i] * pt[i] * pt[i] + it->w[2 * i + 1] * pt[i];
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
				else return 0;
			}
			if (ncmp) return 1;
		}
		return 0;
	}

	template<typename A>
	short Intersect(IloCplex & cpx, IloRangeArray & cons, IloNumVarArray & x, const NDGeom & geom, const NDWindow<A> & node)        //efficient implementation
	{
		IloEnv env = cpx.getEnv();
		for (int i = 0; i < node.nDims; i++)
		{
			x[i].setBounds(node.minPoint[i], node.maxPoint[i]);
		}

		if (cpx.solve())
		{
			return 1;
			double lb = 0;
			double ub = 0;
			for (int i = 0; i < geom.faces.size(); i++)
			{
				lb = cons[i].getLb();
				ub = cons[i].getUb();
				cons[i].setBounds(geom.faces[i].b, geom.faces[i].b);
				if (cpx.solve())
				{
					cons[i].setBounds(lb, ub);
					return 1;       //intersect
				}
				cons[i].setBounds(lb, ub);
			}
			return 2;   //contain
		}

		return 0;
	}


	void RangeTest(const NDGeom & geom, string infile)
	{
		auto dimnum = geom.dimnum;
		SFCConversion<int> sfc(dimnum, 12);
		IloEnv env;     //initialize cplex environment
		IloModel model(env);
		IloObjective obj(env, 0);
		model.add(obj);
		IloNumVarArray x(env, dimnum, 0, 0);
		IloRangeArray cons = geom2cons(env, geom, x);
		model.add(cons);
		IloCplex cplex(model);
		cplex.setOut(env.getNullStream());

		char buf[1024];
		char * pch, *lastpos;
		char ele[64];
		FILE *input_file = fopen(infile.c_str(), "r");
		if (!input_file)
		{
			throw "No input!";
		}
		else
		{
			NDPoint<int> NodeL(dimnum); //lower bound of child node
			NDPoint<int> NodeH(dimnum); //upper bound of child node
			NDWindow<int> cell;
			NodeND node;
			sfc_bigint rangeL;   //lower boundary of a range
			sfc_bigint rangeH;   //upper boundary of a range

			int j;
			int height;
			sfc_bigint key;
			int posNum = 0;
			int negNum = 0;
			while (1)
			{
				memset(buf, 0, 1024);
				fgets(buf, 1024, input_file);

				if (strlen(buf) == 0) break; // no more data

				j = 0;
				lastpos = buf;
				pch = strchr(buf, ',');
				while (pch != 0)
				{
					memset(ele, 0, 64);
					strncpy(ele, lastpos, pch - lastpos);
					if (strlen(ele) != 0)
					{
						height = atoi(ele);
						j++;
					}

					lastpos = pch + 1;
					pch = strchr(lastpos, ',');
				}

				if (strlen(lastpos) != 0 && strcmp(lastpos, "\n") != 0)//final part
				{
					key = (sfc_bigint)lastpos;
					j++;
				}

				rangeL = key << (height*dimnum);
				rangeH = ((key + 1) << (height*dimnum)) - 1;
				NodeL = sfc.MortonDecode(rangeL);
				NodeH = sfc.MortonDecode(rangeH);
				cell.SetMinPoint(NodeL);
				cell.SetMaxPoint(NodeH);
				short mark = Intersect<int>(cplex, cons, x, geom, cell);
				if (mark)
				{
					posNum++;
				}
				else
				{
					negNum++;
				}

			}

			cout << "Positive nodes: " << posNum << ", " << "Negative nodes: " << negNum << endl;
		}
	}

	map <U, U> PlainGeomRange(const NDGeom & geom, short maxdepth, short dimbits = 20)
	{
		measure.histLoad = 0;
		map <U, U> ranges;
		int dimnum = PCDB.nDims;
		multimap <unsigned short, sfc_bigint> SearchT;
		short order = dimbits;
		NDPoint<int> NodeL(dimnum); //lower bound of child node
		NDPoint<int> NodeH(dimnum); //upper bound of child node
		SFCConversion<int> sfc(dimnum, order + 1); //with 1 more bit for the sake of rounding the double type

		IloEnv env;     //initialize cplex environment
		IloModel model(env);
		IloObjective obj(env, 0);
		model.add(obj);
		IloNumVarArray x(env, dimnum, 0, 0);
		IloRangeArray cons = geom2cons(env, geom, x);
		model.add(cons);
		IloCplex cplex(model);
		cplex.setOut(env.getNullStream());

		if (!maxdepth) maxdepth = order + 1;

		SearchT.insert(make_pair(order, 0));
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
				node.key = it->second;
				node.height = it->first;
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
					short mark = Intersect<int>(cplex, cons, x, geom, cell);
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
		env.end();

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

	map <U, U> HistGeomRange(const NDGeom & geom, short dimbits = 20)
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
		long long sumBP = 0;  //points covered by boundary nodes
		long long sumIP = 0;  //points covered by inner nodes

		IloEnv env;     //initialize cplex environment
		IloModel model(env);
		IloObjective obj(env, 0);
		model.add(obj);
		IloNumVarArray x(env, dimnum, 0, 0);
		IloRangeArray cons = geom2cons(env, geom, x);
		model.add(cons);
		IloCplex cplex(model);
		cplex.setOut(env.getNullStream());

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
			short mark = Intersect<int>(cplex, cons, x, geom, cell);
			if (mark == 2)
			{
				ranges.insert(make_pair(rangeL, rangeH));
				sumIP++;
			}

			else if (mark == 1)
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
					bNodes.push_back(node);
					sumBP++;
					//cout << "Space: " << NodeLR[0] << "," << NodeLR[1] << "," << NodeLR[2] << "," << NodeLR[3] << ";" << NodeHR[0] << "," << NodeHR[1] << "," << NodeHR[2] << "," << NodeHR[3] << ", number of points: " << node->pnum << ", ratio: "<< 1.0*node->pnum/(double)(keyH - keyL) << endl;
				}

			}

			else continue;
		}

		auto end1 = chrono::high_resolution_clock::now();
		cout << "Hist search costs: " << chrono::duration_cast<chrono::milliseconds>(end1 - start).count() << "ms" << endl;
		cout << "Inner points: " << sumIP << ", boundary points: " << sumBP << ", boundary nodes: " << bNodes.size() << endl;

		unsigned int cycle_time = MAX_CYCLE;  //threshold for the refinement cycles
		int childnum = 1 << dimnum;
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
					short mark = Intersect<int>(cplex, cons, x, geom, cell);
					if (mark == 2)
					{
						ranges.insert(make_pair(rangeL, rangeH));
						sumIP++;
					}
					else if (mark == 1)
					{
						bNodes.push_back(child);
						sumBP++;
					}
				}
				sumBP--;
			}
			else break;

		}

		for (auto it = bNodes.begin(); it != bNodes.end(); it++)
		{
			ranges.insert(make_pair((*it)->key << (*it)->height * dimnum, (((*it)->key + 1) << (*it)->height * dimnum) - 1));
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
		env.end();
		return ranges_merged;
	}

public:
	QueryExtendCP() {
		measure = {};
		QueryGeom = {};
		PCDB = {};
	}

	QueryExtendCP(const PointCloudDB<T, U>& PC)
	{
		measure = {};
		QueryGeom = {};
		PCDB = PC;
		if (PCDB.HIST)
		{
			auto start = chrono::high_resolution_clock::now();
			histroot = this->HistLoad(PCDB.HistTab);
			auto end = chrono::high_resolution_clock::now();
			measure.histLoad = chrono::duration_cast<chrono::milliseconds>(end - start).count();
		}
	}

	/*specifically for frustum view selection*/
	void QueryView(const NDGeom & poly, const NDGeom & geom)
	{
		QueryGeom = geom;
		int dimnum = PCDB.nDims;
		NDGeom polytrans = poly.Transform(PCDB.trans);
		NDGeom geomtrans = geom.Transform(PCDB.trans);

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
			ranges = HistGeomRange(polytrans, dimbits);
			auto end = chrono::high_resolution_clock::now();
			measure.rangeComp = chrono::duration_cast<chrono::milliseconds>(end - start).count();
		}
		else
		{
			auto start = chrono::high_resolution_clock::now();
			ranges = PlainGeomRange(polytrans, 12, dimbits); //search depth can be modified depending on accuracy requirement
			auto end = chrono::high_resolution_clock::now();
			measure.rangeComp = chrono::duration_cast<chrono::milliseconds>(end - start).count();
		}

		Environment *env = Environment::createEnvironment(Environment::DEFAULT);
		Connection  *con = env->createConnection(orclconn().User, orclconn().Password, orclconn().Database);

		auto start = chrono::high_resolution_clock::now();

		bool varchar_mode = false;      //open this mode when key is larger than the number type of Oracle
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

		ofstream output_file("E:/queryCP_res_acc.csv");
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
				ptreal = pt.InverseTransform<double>(PCDB.trans);
				output_file << ptreal[0] << ", " << ptreal[1] << ", " << ptreal[2] << ", " << ptreal[3] << "\n";
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


	void QueryIOT(const NDGeom & geom)
	{
		QueryGeom = geom;
		int dimnum = PCDB.nDims;
		NDGeom geomtrans = geom.Transform(PCDB.trans);

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
			ranges = HistGeomRange(geomtrans, dimbits);
			auto end = chrono::high_resolution_clock::now();
			measure.rangeComp = chrono::duration_cast<chrono::milliseconds>(end - start).count();
		}
		else
		{
			auto start = chrono::high_resolution_clock::now();
			ranges = PlainGeomRange(geomtrans, 12, dimbits); //search depth can be modified depending on accuracy requirement
			auto end = chrono::high_resolution_clock::now();
			measure.rangeComp = chrono::duration_cast<chrono::milliseconds>(end - start).count();
		}

		Environment *env = Environment::createEnvironment(Environment::DEFAULT);
		Connection  *con = env->createConnection(orclconn().User, orclconn().Password, orclconn().Database);

		auto start = chrono::high_resolution_clock::now();

		bool varchar_mode = false;      //open this mode when key is larger than the number type of Oracle
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

		ofstream output_file("E:/queryCP_res_plain.csv");
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
			output_file << pt[0] << ", " << pt[1] << ", " << pt[2] << ", " << pt[3] << "\n";
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
		output << "Method: CPLEX" << endl;
		if (PCDB.HIST) output << "HistSFC\n";
		else output << "PlainSFC\n";

		output << "rangeNum, appPNum, accPNum, FPR, rangeComp, histLoad, firstCost, secondCost\n";
		output << measure.rangeNum << ", " << measure.appPNum << ", " << measure.accPNum << ", " << measure.FPR << ", "
			<< measure.rangeComp << ", " << measure.histLoad << ", " << measure.firstCost << ", " << measure.secondCost << "\n";
		output << "\n";
	}
};
