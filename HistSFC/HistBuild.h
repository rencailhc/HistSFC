#pragma once
#include <map>
#include <fstream>
#include <string>
#include "Point.h"
#include "BaseStruct.h"
#include "SFCConversion.h"
#include "Oracle.h"

using namespace std;
#define FETCH_SIZE (unsigned int)(1 << 20)

/*Two methods for computing HistTree, HistML and HistIOT*/

HistNodeND *HistML(string infile, const CoordTrans<double>& trans, int mincap = 100, short dimnum = 3, short maxbits=30)
{
	//HistTree built during loading, based on original xyz file, return root node of the histogram tree
	//maxbits: maximum number of bits per dimension
	map<sfc_bigint, HistNodeND *> L;
	SFCConversion<int> sfc(dimnum, maxbits);
	int Theight = 10;   //can be derived by testing on partial data

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
		int i, j;
		i = 0;

		NDPoint<double> inPt;
		NDPoint<double> pt(dimnum);
		NDPoint<int> convertpt(dimnum);
		ofstream output_file("E:/91_sfc_xyl.csv");
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
					inPt[j] = atof(ele);
					j++;
				}

				lastpos = pch + 1;
				pch = strchr(lastpos, ',');
			}

			if (strlen(lastpos) != 0 && strcmp(lastpos, "\n") != 0)//final part
			{
				inPt[j] = atof(lastpos);
				j++;
			}
			//extract organizing dimensions
			pt[0] = inPt[0];
			pt[1] = inPt[1];
			pt[2] = inPt[3];

			inPt = pt.Transform<double>(trans);
			convertpt = { (int)inPt[0], (int)inPt[1],(int)inPt[2] };
			sfc_bigint Okey = sfc.MortonEncode(convertpt);
			output_file << Okey << "\n";
			
			sfc_bigint key = Okey >> dimnum * Theight;
			auto it = L.find(key);
			if (it != L.end())
			{
				it->second->pnum++;
			}
			else
			{
				//nodenum++;
				HistNodeND *node = (HistNodeND*)malloc(sizeof(HistNodeND));
				node->child = 0;
				node->neighbor = 0;
				node->key = key;
				node->pnum = 1;
				node->cnum = 0;
				node->height = Theight;
				L.insert(make_pair(key, node));
			}
		}
		//cout << "First batch node number: " << nodenum << endl;
	}
	fclose(input_file);


	map <sfc_bigint, HistNodeND *> FixPool;   //Fix the HistTree below the base level
	for (auto it = L.begin(); it != L.end(); it++)
	{
		if (it->second->pnum >= mincap)
		{
			FixPool.insert(make_pair(it->second->key, it->second));
		}
	}
	//cout << FixPool.size() << endl;
	//tbb::memory_pool<std::allocator<HistNodeND> > ExtPool;

	ifstream input_sfc("E:/91_sfc_xyl.csv");
	if (!input_sfc.is_open())
	{
		cout << "Error opening sfc file";
		exit(1);
	}

	//attach the children of overloaded nodes
	if (!FixPool.empty())
	{
		while (!input_sfc.eof())
		{
			input_sfc.getline(buf, 1024);
			if (strlen(buf) == 0) break;
			sfc_bigint key = (sfc_bigint)buf;
			sfc_bigint key_pre = key >> dimnum * Theight;
			auto it = FixPool.find(key_pre);
			if (it != FixPool.end())
			{
				HistNodeND *node = (HistNodeND*)malloc(sizeof(HistNodeND));
				if (!it->second->child)
				{
					it->second->child = node;
					node->neighbor = 0;
				}
				else
				{
					node->neighbor = it->second->child;
					it->second->child = node;
				}
				node->child = 0;
				node->key = key;
				node->pnum = 1;
				node->cnum = 0;
				node->height = 0;
			}

		}
	}

	//Fix nodes below middle level
	while (!FixPool.empty())
	{
		auto it = FixPool.begin();
		HistNodeND *PriorNode = 0;  //used for passing neighbor pointer
		int cnum = 0;

		for (int i = 0; i < 1 << dimnum; i++)
		{
			sfc_bigint key = (it->first << dimnum) + i;
			HistNodeND *node = it->second->child;
			long pnum = 0;
			while (node)
			{
				if (node->key >> ((it->second->height - 1)*dimnum) == key)
				{
					pnum++;
				}
				node = node->neighbor;
			}
			if (pnum == 0) continue;
			else
			{
				cnum++;
				HistNodeND *cnode = (HistNodeND*)malloc(sizeof(HistNodeND));
				cnode->neighbor = PriorNode;
				cnode->key = key;
				cnode->pnum = pnum;
				cnode->cnum = 0;
				if (pnum < mincap)
				{
					cnode->child = 0;
				}
				else
				{
					cnode->child = it->second->child;   //attach all node at bottom level to it
					FixPool.insert(make_pair(key, cnode));
				}
				cnode->height = it->second->height - 1;
				PriorNode = cnode;
			}

		}
		it->second->child = PriorNode;
		it->second->cnum = cnum;
		FixPool.erase(it);
	}
	//ExtPool.recycle();   //free all nodes at bottom level

	//build nodes above middle level
	sfc_bigint key_prefix = L.rbegin()->first;
	while (key_prefix != 0)   //if it is root node
	{
		map<sfc_bigint, HistNodeND *> PL;   //parent node list

		//cout << L.size() << endl;
		for (auto it = L.begin(); it != L.end(); it++)
		{
			sfc_bigint p_key = it->second->key >> dimnum;
			auto it2 = PL.find(p_key);
			if (it2 != PL.end())
			{
				it2->second->pnum += it->second->pnum;
				it2->second->cnum++;
				it->second->neighbor = it2->second->child;
				it2->second->child = it->second;
			}
			else
			{
				HistNodeND *pnode = (HistNodeND*)malloc(sizeof(HistNodeND));
				pnode->child = it->second;
				pnode->neighbor = 0;
				pnode->key = p_key;
				pnode->pnum = it->second->pnum;
				pnode->cnum = 1;
				pnode->height = it->second->height + 1;
				PL.insert(make_pair(p_key, pnode));
			}
			//if(it->second->height<=log2(t)/3) free(it->second);  //free original node to save memory

		}
		//cout << "height," << PL.begin()->second->height << endl;
		for (auto it2 = PL.begin(); it2 != PL.end(); it2++)
		{
			if (it2->second->pnum < mincap)
			{
				HistNodeND *curnode = it2->second->child;
				HistNodeND *nextnode = 0;
				for (int i = 0; i < it2->second->cnum; i++)
				{
					nextnode = curnode->neighbor;
					free(curnode);   //free original node to save memory
					curnode = nextnode;
				}

				it2->second->child = 0;
				it2->second->cnum = 0;
			}
		}
		L.swap(PL);
		PL.clear();
		key_prefix = L.rbegin()->first;  //the last element has the largest key

	}

	return L.begin()->second;

}


HistNodeND *HistIOT(string iot_table, int mincap = 100, short dimnum = 3)
{
	//Histogram built from the IOT table, linear time cost, but can keep memory usage small
	//mincap: minimum capacity of leaf node
	HistNodeND *HistRoot = (HistNodeND*)malloc(sizeof(HistNodeND));
	HistRoot->key = 0;
	HistRoot->child = 0;
	HistRoot->neighbor = 0;
	map<sfc_bigint, HistNodeND *> L;

	Environment *env = Environment::createEnvironment(Environment::DEFAULT);
	Connection  *con = env->createConnection(orclconn().User, orclconn().Password, orclconn().Database);
	Statement *stmt = NULL;
	ResultSet *rs = NULL;
	string sql = "select sfc from " + iot_table;
	stmt = con->createStatement(sql);
	stmt->setPrefetchRowCount(FETCH_SIZE);
	rs = stmt->executeQuery();
	rs->next();
	sfc_bigint key_start = (sfc_bigint)rs->getString(1);
	sfc_bigint key_prefix = key_start >> dimnum;

	/////////Fix pool size//////////
	//long long nodesize = 0;    
	/////////////////////////////

	long long pnum = 1;
	short height = 1;
	int thres = mincap / (1 << dimnum);
	SFCcell *FixP = (SFCcell*)malloc(sizeof(SFCcell));
	FixP->key = 0;
	FixP->next = 0;
	SFCcell *fp = (SFCcell*)malloc(sizeof(SFCcell));
	fp->key = key_start;
	fp->next = 0;
	FixP->next = fp;
	SFCcell * gpointer = fp; //global pointer

	sfc_bigint key;
	while (rs->next()) {
		key = (sfc_bigint)rs->getString(1);
		sfc_bigint key_c = (key >> dimnum * height);
		if (pnum <= thres or key_c == key_prefix)
		{
			pnum++;
			int cyc = 0;
			while (key_c != key_prefix and ++cyc)
			{
				if (key_c >> dimnum * cyc == key_prefix >> dimnum * cyc) break;
			}
			//if (key_c >> dimnum * cyc == 51615649175) cout << HistRoot->neighbor->key << ", "<<abs(height + cyc - HistRoot->neighbor->height) <<endl;
			if (HistRoot->neighbor and key_c >> dimnum * cyc == HistRoot->neighbor->key >> dimnum * abs(height + cyc - HistRoot->neighbor->height))
			{
				//judge whether later SFC block built contains previous block
				HistNodeND *node = (HistNodeND*)malloc(sizeof(HistNodeND));
				node->child = 0;
				node->neighbor = HistRoot->neighbor;
				HistRoot->neighbor = node;
				node->key = key_prefix;
				node->pnum = pnum - 1;
				node->cnum = 0;
				node->height = height;
				pnum = 1;
				height = 1;
				key_prefix = key >> dimnum;
				gpointer = FixP;
			}
			else
			{
				key_prefix = key_prefix >> dimnum * cyc;
				height += cyc;
			}
			//cout << key_prefix << ", " << height << ", " << pnum << endl;
		}
		else
		{
			//cout << key_prefix <<", " << pnum << endl;
			HistNodeND *node = (HistNodeND*)malloc(sizeof(HistNodeND));
			node->child = 0;
			node->neighbor = HistRoot->neighbor;
			HistRoot->neighbor = node;
			node->key = key_prefix;
			node->pnum = pnum;
			node->cnum = 0;
			node->height = height;
			pnum = 1;
			height = 1;
			key_prefix = key >> dimnum;
			if (node->pnum >= mincap)
			{
				map <HistNodeND *, SFCcell *> FixN;
				FixN.insert(make_pair(node, FixP->next));
				SFCcell *tpointer;  //turning pointer
				while (!FixN.empty())
				{
					auto it = FixN.begin();
					HistNodeND *PriorNode = 0;  //used for passing neighbor pointer
					int cnum = 0;
					SFCcell *fpointer = it->second;
					//output_file << "new node: " << it->first->key << ", pnum: " << it->first->pnum << "\n";
					//cout << fpointer->key << ", " << fpointer->next->key << ", " << fpointer->next->next->key << endl;
					for (int i = 0; i < (1 << dimnum); i++)
					{
						tpointer = fpointer;
						sfc_bigint key_fix = (it->first->key << dimnum) + i;
						int pnum = 0;
						while (fpointer != gpointer->next)
						{
							//cout << (PPool.back() >> ((it->first->height - 1)*dimnum)) << ", " << key_fix << endl;
							if ((fpointer->key >> ((it->first->height - 1)*dimnum)) == key_fix)
							{
								pnum++;
								fpointer = fpointer->next;
							}
							else break;
						}
						if (pnum > 0)
						{
							//output_file << "child: " << key_fix << ", pnum: " << pnum << "\n";
							cnum++;
							HistNodeND *cnode = (HistNodeND*)malloc(sizeof(HistNodeND));
							cnode->neighbor = PriorNode;
							cnode->child = 0;
							cnode->cnum = 0;
							cnode->key = key_fix;
							cnode->pnum = pnum;
							cnode->height = it->first->height - 1;
							PriorNode = cnode;
							if (pnum >= mincap)
							{
								FixN.insert(make_pair(cnode, tpointer));
							}
						}

					}
					it->first->child = PriorNode;
					it->first->cnum = cnum;
					FixN.erase(it);
				}
			}
			gpointer = FixP;

		}

		if (!gpointer->next)
		{
			//nodesize++;
			fp = (SFCcell*)malloc(sizeof(SFCcell));
			fp->next = gpointer->next;
			gpointer->next = fp;
		}
		else fp = gpointer->next;
		fp->key = key;
		gpointer = fp;

	}

	HistNodeND *node = (HistNodeND*)malloc(sizeof(HistNodeND)); //The last node built
	node->child = 0;
	node->neighbor = HistRoot->neighbor;
	HistRoot->neighbor = node;
	node->key = key_prefix;
	node->pnum = pnum;

	node->cnum = 0;
	node->height = height;
	if (node->pnum >= mincap)
	{
		map <HistNodeND *, SFCcell *> FixN;
		FixN.insert(make_pair(node, FixP->next));
		SFCcell *tpointer;  //turning pointer
		while (!FixN.empty())
		{
			auto it = FixN.begin();
			SFCcell *fpointer = it->second;
			HistNodeND *PriorNode = 0;  //used for passing neighbor pointer
			int cnum = 0;

			for (int i = 0; i < (1 << dimnum); i++)
			{
				tpointer = fpointer;
				sfc_bigint key_fix = (it->first->key << dimnum) + i;
				int pnum = 0;
				while (fpointer != gpointer->next)
				{
					//cout << (PPool.back() >> ((it->first->height - 1)*dimnum)) << ", " << key_fix << endl;
					if ((fpointer->key >> ((it->first->height - 1)*dimnum)) == key_fix)
					{
						pnum++;
						fpointer = fpointer->next;
					}
					else break;
				}
				if (pnum > 0)
				{
					cnum++;
					HistNodeND *cnode = (HistNodeND*)malloc(sizeof(HistNodeND));
					cnode->neighbor = PriorNode;
					cnode->child = 0;
					cnode->cnum = 0;
					cnode->key = key_fix;
					cnode->pnum = pnum;
					cnode->height = it->first->height - 1;
					PriorNode = cnode;
					if (pnum >= mincap)
					{
						FixN.insert(make_pair(cnode, tpointer));
					}
				}
			}
			it->first->child = PriorNode;
			it->first->cnum = cnum;
			FixN.erase(it);
		}
	}
	gpointer = FixP;

	/*build nodes at upper levels*/
	node = HistRoot->neighbor;
	short minh = node->height;
	short maxh = node->height;
	while (node)
	{
		if (node->height < minh) minh = node->height;
		if (node->height > maxh) maxh = node->height;
		node = node->neighbor;
	}

	//ofstream output_fil2("E:/middle_layer.csv");
	node = HistRoot->neighbor;
	HistNodeND * prev = HistRoot;
	while (node)
	{
		//output_fil2 << node->key << "\n";
		if (node->height == minh)
		{
			prev->neighbor = node->neighbor;
			L.insert(make_pair(node->key, node));
		}
		else prev = node;
		node = node->neighbor;
	}

	key_prefix = L.rbegin()->first;
	while (key_prefix != 0)   //if it is root node
	{
		map<sfc_bigint, HistNodeND *> PL;   //parent node list

		if (L.begin()->second->height + 1 <= maxh)
		{
			node = HistRoot->neighbor;
			prev = HistRoot;
			while (node)
			{
				if (node->height == L.begin()->second->height + 1)
				{
					PL.insert(make_pair(node->key, node));
					prev->neighbor = node->neighbor;
				}
				else prev = node;
				node = node->neighbor;
			}
		}

		for (auto it = L.begin(); it != L.end(); it++)
		{
			sfc_bigint p_key = it->first >> dimnum;
			auto it2 = PL.find(p_key);
			if (it2 != PL.end())
			{
				it2->second->pnum += it->second->pnum;
				it2->second->cnum++;
				it->second->neighbor = it2->second->child;
				it2->second->child = it->second;
			}
			else
			{
				HistNodeND *pnode = (HistNodeND*)malloc(sizeof(HistNodeND));
				pnode->child = it->second;
				pnode->neighbor = 0;
				pnode->key = p_key;
				pnode->pnum = it->second->pnum;
				pnode->cnum = 1;
				pnode->height = it->second->height + 1;
				PL.insert(make_pair(p_key, pnode));
			}
			//if(it->second->height<=log2(t)/3) free(it->second);  //free original node to save memory

		}

		for (auto it2 = PL.begin(); it2 != PL.end(); it2++)
		{
			if (it2->second->pnum < mincap)
			{
				HistNodeND *curnode = it2->second->child;
				HistNodeND *nextnode = 0;
				for (int i = 0; i < it2->second->cnum; i++)
				{
					nextnode = curnode->neighbor;
					free(curnode);   //free original node to save memory
					curnode = nextnode;
				}

				it2->second->child = 0;
				it2->second->cnum = 0;
			}
		}

		L.swap(PL);
		PL.clear();
		key_prefix = L.rbegin()->first;  //the last element has the largest key

	}

	HistRoot = L.begin()->second;

	stmt->closeResultSet(rs);
	con->terminateStatement(stmt);
	env->terminateConnection(con);

	//cout << "Refinement pool size: " << nodesize << endl;

	return HistRoot;
}



void ExTree(HistNodeND *histroot, string table_name)
{
	//export HistTree to a flat table, table_name: HistTree table name
	vector <HistNodeND *> SearchT;
	SearchT.push_back(histroot);
	HistNodeND *node;

	try
	{
		Environment *env = Environment::createEnvironment(Environment::DEFAULT);
		Connection  *con = env->createConnection(orclconn().User, orclconn().Password, orclconn().Database);

		string sql = "create table " + table_name +
			" (id number, sfc number, pnum number, cnum number, height number, child number, neighbor number)";
		Statement *f_stmt = con->createStatement();
		f_stmt->executeUpdate(sql);
		con->terminateStatement(f_stmt);

		sql = "insert into " + table_name +
			" values ( :val_id, :val_key, :val_pnum, :val_cnum," +
			" :val_height, :val_child, :val_neighbor)";

		Statement *stmt = con->createStatement();
		stmt->setSQL(sql);
		stmt->setMaxParamSize(1, 30);
		stmt->setMaxParamSize(2, 30);
		stmt->setMaxParamSize(3, 30);
		stmt->setMaxParamSize(6, 30);
		stmt->setMaxParamSize(7, 30);
		while (!SearchT.empty())
		{
			node = SearchT.back();
			SearchT.pop_back();

			if (node->cnum != 0)
			{
				HistNodeND *curnode = node->child;
				for (int i = 0; i < node->cnum; i++)
				{
					SearchT.push_back(curnode);
					curnode = curnode->neighbor;
				}
			}

			ostringstream keystr;
			keystr << (long long)node;
			string id = keystr.str();
			stmt->setString(1, id);
			keystr.clear();
			keystr.str("");

			keystr << node->key;
			string key = keystr.str();
			stmt->setString(2, key);
			keystr.clear();
			keystr.str("");

			keystr << node->pnum;
			string pnum = keystr.str();
			stmt->setString(3, pnum);
			keystr.clear();
			keystr.str("");

			stmt->setInt(4, node->cnum);
			stmt->setInt(5, node->height);

			keystr << (long long)node->child;
			string child = keystr.str();
			stmt->setString(6, child);
			keystr.clear();
			keystr.str("");

			keystr << (long long)node->neighbor;
			string neighbor = keystr.str();
			stmt->setString(7, neighbor);
			keystr.clear();
			keystr.str("");

			stmt->executeUpdate();
		}

		con->commit();
		con->terminateStatement(stmt);
		env->terminateConnection(con);
		Environment::terminateEnvironment(env);
	}

	catch (std::exception &ex)
	{
		cout << "Problem: " << ex.what() << endl;
	}

}

