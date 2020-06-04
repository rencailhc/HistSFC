#pragma once
#include "typedef.h"
#include <initializer_list>


typedef struct HistNodeND   //HistTree node
{
	HistNodeND* child;     // pointer to the first child
	HistNodeND* neighbor;  // pointer to a sibling
	sfc_bigint key; //sfc prefix
	unsigned long long pnum;  // number of points covered by the node
	unsigned short cnum;  // number of children
	unsigned short height;  //height of the node in the tree
} HistNodeND;


typedef struct NodeND  //used for PlainSFC recursion
{
	sfc_bigint key;			//sfc key
	unsigned short height;     //height of the node during searching
	//dimensionality is not stored, as it is computed during searching
} NodeND;


typedef struct SFCcell    //used for building HistTree
{
	sfc_bigint key;
	SFCcell *next;
} SFCcell;


typedef struct ConnParam   //used for connecting database
{
	const char* Database;
	const char* User;
	const char* Password;
} ConnParam;


typedef struct Measurement    //used for measuring performance
{
	//all time cost is in ms
	unsigned int rangeNum;	//number of ranges generated
	unsigned long long appPNum;	//number of point of first filter
	unsigned long long accPNum;	//exact answer
	double FPR;	//false positive rate
	double rangeComp;	//time to compute the ranges
	double histLoad;	//time to load HistTree
	double firstCost;	//time cost of the first filter
	double secondCost;	//time cost of the second filter
} Measurement;


template <typename T>    //used for transforming coordinates
struct CoordTrans
{
	const short dimnum;
	T* _delta;
	T* _scale;

	CoordTrans():dimnum(0)
	{
		_delta = nullptr;
		_scale = nullptr;
	}

	CoordTrans(short dims):dimnum(dims)
	{
		_delta = new T[dimnum];
		_scale = new T[dimnum];
	}

	CoordTrans(std::initializer_list<T> l1, std::initializer_list<T> l2) :dimnum((int)l1.size())
	{
		_delta = new T[dimnum];
		_scale = new T[dimnum];
		std::uninitialized_copy(l1.begin(), l1.end(), _delta);
		std::uninitialized_copy(l2.begin(), l2.end(), _scale);
	}

	CoordTrans& operator=(CoordTrans&& other)
	{
		//assignment by move
		(short &)dimnum = other.dimnum;
		delete[] _delta;
		delete[] _scale;
		_delta = other._delta;
		_scale = other._scale;
		return *this;
	}

	CoordTrans& operator=(const CoordTrans& other)
	{
		//assignment by copy
		(short &)dimnum = other.dimnum;
		delete[] _delta;
		delete[] _scale;
		_delta = new T[dimnum];
		_scale = new T[dimnum];
		for (int i = 0; i < dimnum; i++)
		{
			_delta[i] = other._delta[i];
			_scale[i] = other._scale[i];
		}
		
		return *this;
	}
};