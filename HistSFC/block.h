#pragma once
#include <map>
#include <fstream>
#include <string>
#include "Point.h"
#include "BaseStruct.h"
#include "SFCConversion.h"
#include "Oracle.h"

HistNodeND *BlockLoad(string HistTab)
{
	HistNodeND *HistRoot = (HistNodeND*)malloc(sizeof(HistNodeND));
	map <long long, HistNodeND *> HistNodePool;

	try
	{
		Environment *env = Environment::createEnvironment(Environment::DEFAULT);

		Connection  *con = env->createConnection(orclconn().User, orclconn().Password, orclconn().Database);
		Statement *stmt = NULL;
		ResultSet *rs = NULL;
		string sql = "select id from " + HistTab + " where height = (select max(height) from " + HistTab + " )";
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
			curnode->pnum = id;
			curnode->cnum = rs->getInt(4);
			curnode->height = rs->getInt(5);
			curnode->child = 0;
			curnode->neighbor = 0;

			long long child_id = stoll(rs->getString(6));
			if (child_id)
			{
				it = HistNodePool.find(child_id);
				if (it != HistNodePool.end()) curnode->child = it->second;
				else
				{
					HistNodeND *child = (HistNodeND*)malloc(sizeof(HistNodeND));
					curnode->child = child;
					HistNodePool.insert(make_pair(child_id, child));
				}
			}
			
			long long neighbor_id = stoll(rs->getString(7));
			if (neighbor_id)
			{
				it = HistNodePool.find(neighbor_id);
				if (it != HistNodePool.end()) curnode->neighbor = it->second;
				else
				{
					HistNodeND *neighbor = (HistNodeND*)malloc(sizeof(HistNodeND));
					curnode->neighbor = neighbor;
					HistNodePool.insert(make_pair(neighbor_id, neighbor));
				}
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

void BlockCalc(string infile, string outfile, const CoordTrans& trans, HistNodeND *blocktree)
{

	auto dimnum = trans.dimnum;
	SFCConversion<int> sfc(dimnum, blocktree->height);
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
		NDPoint<double> inPt(dimnum);
		NDPoint<int> convertpt(dimnum);
		ofstream output_file(outfile);
		while (1)
		{
			memset(buf, 0, 1024);
			fgets(buf, 1024, input_file);

			if (strlen(buf) == 0) break; // no more data

			int j = 0;
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

			convertpt = inPt.Transform<int>(trans);
			sfc_bigint Okey = sfc.MortonEncode(convertpt);
			//output_file << Okey << "\n";

			HistNodeND *iter = blocktree;
			HistNodeND *block;
			while (iter)
			{
				block = iter;
				if (Okey >> block->height*dimnum == block->key)
				{
					iter = block->child;
				}
				else
				{
					iter = block->neighbor;
				}
			}

			if (Okey >> block->height*dimnum != block->key)
			{
				cout << "wrong!" << endl;
				break;
			}

			auto id = block->pnum;  //pnum is replaced by id in blocktree

			for (int i = 0; i < dimnum; i++)
			{
				output_file << inPt[i] << ", ";
			}

			output_file << id << endl;
		}
	}
	fclose(input_file);
}

void BlockWrite(string infile, string outfile, string HistTab)
{
	CoordTrans trans({ 13427.64, 359007.3, -8.79, 12000 }, { 100,100,1000,-100 });
	BlockCalc(infile, outfile, trans, BlockLoad(HistTab));
}