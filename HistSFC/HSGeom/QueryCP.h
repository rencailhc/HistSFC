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
	return cons;

}

template <typename T, typename U>	//data type of coordinates and SFC key
class QueryExtendCP : public Query <T, U> {
	using Query<T, U>::measure;
	using Query<T, U>::histroot;
	using Query<T, U>::PCDB;
private:
	NDGeom QueryGeom;

private:
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
			for (auto it=geom.faces.begin(); it!=geom.faces.end(); it++)
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

	template<typename A>
	short Intersect(IloCplex & cpx, IloRangeArray & cons, IloNumVarArray & x, const NDGeom & geom, const NDWindow<A> & node)	//efficient implementation
	{
		IloEnv env = cpx.getEnv();
		for (int i = 0; i < node.nDims; i++)
		{
			x[i].setBounds(node.minPoint[i], node.maxPoint[i]);
		}
		
		if (cpx.solve())
		{
			double lb = 0;
			double ub = 0;
			for (int i=0; i<geom.faces.size();i++)
			{
				lb = cons[i].getLb();
				ub = cons[i].getUb();
				cons[i].setBounds(geom.faces[i].b, geom.faces[i].b);
				if (cpx.solve())
				{
					cons[i].setBounds(lb, ub);
					return 1;	//intersect
				}
				cons[i].setBounds(lb, ub);
			}
			return 2;   //contain
		}

		return 0;
	}

	template<typename A>
	short Intersect(IloEnv env, const NDGeom & geom, const NDWindow<A> & node)	//slow implementation
	{
		if (node.nDims != geom.dimnum)
		{
			throw("Dimensionality does not match!");
		}

		IloModel model(env);
		IloObjective obj(env, 0);
		model.add(obj);
		IloNumVarArray x = cell2var<A>(env, node);
		IloRangeArray cons = geom2cons(env, geom, x);
		model.add(cons);
		IloCplex cplex(model);
		cplex.setOut(env.getNullStream());
		
		if (cplex.solve())
		{
			model.remove(cons);
			IloRange con;
			for (auto it = geom.faces.begin(); it != geom.faces.end(); it++)
			{
				halfspace h = *it;
				h.s = Sign::eq;
				con = halfspace2con(env, h, x);
				model.add(con);
				if (cplex.solve()) return 1;	//intersect
				model.remove(con);
			}	
			return 2;   //contain
		}

		return 0;
	}

	map <U, U> PlainGeomRange(const NDGeom & geom, short maxdepth, short dimbits = 20)
	{
		measure.histLoad = 0;
		map <U, U> ranges;
		int dimnum = PCDB.nDims;
		vector <NodeND> SearchT;
		short order = dimbits;
		NDPoint<int> NodeL(dimnum); //lower bound of child node
		NDPoint<int> NodeH(dimnum); //upper bound of child node
		SFCConversion<int> sfc(dimnum, order + 1); //with 1 more bit for the sake of rounding the double type
		
		IloEnv env;	//initialize cplex environment
		IloModel model(env);
		IloObjective obj(env, 0);
		model.add(obj);
		IloNumVarArray x(env, dimnum, 0, 0);
		IloRangeArray cons = geom2cons(env, geom, x);
		model.add(cons);
		IloCplex cplex(model);
		cplex.setOut(env.getNullStream());

		if (!maxdepth) maxdepth = order + 1;

		NodeND rootnode = { 0, order };
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
				short mark = Intersect<int>(cplex, cons, x, geom, cell);
				if (mark == 2)
				{
					ranges.insert(make_pair(rangeL, rangeH));
				}

				else if (mark == 1)
				{
					if (order - node.height < maxdepth and cycle <= cycle_time)
					{
						for (int i = 0; i < 1 << dimnum; i++)
						{
							NodeND child = { (node.key << dimnum) + i, node.height - 1 };
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
		env.end();

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

		IloEnv env;	//initialize cplex environment
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
				sumIP += node->pnum;
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
					sumBP += node->pnum / pow(2, node->cnum);
					//cout << "Space: " << NodeLR[0] << "," << NodeLR[1] << "," << NodeLR[2] << "," << NodeLR[3] << ";" << NodeHR[0] << "," << NodeHR[1] << "," << NodeHR[2] << "," << NodeHR[3] << ", number of points: " << node->pnum << ", ratio: "<< 1.0*node->pnum/(double)(keyH - keyL) << endl;
				}

			}

			else continue;
		}

		auto end1 = chrono::high_resolution_clock::now();
		cout << "Hist search costs: " << chrono::duration_cast<chrono::milliseconds>(end1 - start).count() << "ms" << endl;
		cout << "Inner points: " << sumIP << ", boundary points: " << sumBP << ", boundary nodes: " << bNodes.size() << endl;

		unsigned int cycle_time = MAX_CYCLE;  //threshold for the refinement cycles
		unsigned int cycle = bNodes.size();
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
				short mark = Intersect<int>(cplex, cons, x, geom, cell);
				if (mark == 2)
				{
					ranges.insert(make_pair(rangeL, rangeH));
				}
				else if (mark == 1)
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
		env.end();
		return ranges;
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
			ranges = PlainGeomRange(geomtrans, 7, dimbits);	//search depth can be modified depending on accuracy requirement
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
		output << "Method: CPLEX" << endl;
		if (PCDB.HIST) output << "HistSFC\n";
		else output << "PlainSFC\n";

		output << "rangeNum, appPNum, accPNum, FPR, rangeComp, histLoad, firstCost, secondCost\n";
		output << measure.rangeNum << ", " << measure.appPNum << ", " << measure.accPNum << ", " << measure.FPR << ", "
			<< measure.rangeComp << ", " << measure.histLoad << ", " << measure.firstCost << ", " << measure.secondCost << "\n";
		output << "\n";
	}
};