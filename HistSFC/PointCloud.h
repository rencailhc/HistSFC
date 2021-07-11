/*PointCloud database class*/

#pragma once
#include "Point.h"

template <typename T, typename U>  //the type of original record and database key, respectively
class PointCloudDB
{
public:
	string Table;
	short nDims;
	int SRID;
	CoordTrans trans;
	bool HIST;
	string HistTab;

public:
	PointCloudDB() :nDims(0)
	{
		Table = "";
		SRID = 0;
		trans = {};
		HIST = false;
		HistTab = "";
	}

	PointCloudDB(string Tab, const int dim) :nDims(dim)
	{
		Table = Tab;
		SRID = 0;
		trans.dimnum = nDims;
		trans._delta = new double[nDims];
		trans._scale = new double[nDims];
		for (int i = 0; i < nDims; i++)
		{
			trans._delta[i] = 0;
			trans._scale[i] = 1;
		}
		HIST = false;
		HistTab = "";
	}

	PointCloudDB(string Tab, const int dim, int srs, const CoordTrans& tr) :nDims(dim)
	{
		Table = Tab;
		SRID = srs;
		trans = tr;
		HIST = false;
		HistTab = "";
	}

	PointCloudDB(string Tab, const int dim, int srs, const CoordTrans& tr, bool hist, string histtab) :nDims(dim)
	{
		Table = Tab;
		SRID = srs;
		trans = tr;
		HIST = hist;
		HistTab = histtab;
	}

	PointCloudDB(const PointCloudDB& other)
	{
		Table = other.Table;
		nDims = other.nDims;
		SRID = other.SRID;
		trans = other.trans;
		HIST = other.HIST;
		HistTab = other.HistTab;
	}


	PointCloudDB& operator=(const PointCloudDB& other)
	{
		Table = other.Table;
		nDims = other.nDims;
		SRID = other.SRID;
		trans = other.trans;
		HIST = other.HIST;
		HistTab = other.HistTab;
		return *this;
	}

};

template <typename T, typename U>
class PyramidDB :public PointCloudDB<T, U> {
public:
	bool Extend;
	double* _medians;

public:
	PyramidDB()
	{
		this->Table = "";
		this->nDims = 0;
		this->SRID = 0;
		this->trans = {};
		Extend = false;
		_medians = nullptr;
	}

	PyramidDB(string Tab, const int dim)
	{
		this->Table = Tab;
		this->nDims = dim;
		this->SRID = 0;
		this->trans.dimnum = dim;
		this->trans._delta = new double[dim];
		this->trans._scale = new double[dim];
		for (int i = 0; i < dim; i++)
		{
			this->trans._delta[i] = 0;
			this->trans._scale[i] = 1;
		}
		Extend = 0;
		_medians = nullptr;
	}

	PyramidDB(string Tab, const int dim, int srs, const CoordTrans& tr)
	{
		this->Table = Tab;
		this->nDims = dim;
		this->SRID = srs;
		this->trans = tr;
		Extend = 0;
		_medians = nullptr;
	}

	PyramidDB(string Tab, const int dim, int srs, const CoordTrans& tr, bool ex, T* medians)
	{
		this->Table = Tab;
		this->nDims = dim;
		this->SRID = srs;
		this->trans = tr;
		Extend = ex;
		_medians = medians;
	}

	PyramidDB(string Tab, const int dim, int srs, const CoordTrans& tr, bool ex, std::initializer_list<T> list)
	{
		this->Table = Tab;
		this->nDims = dim;
		this->SRID = srs;
		this->trans = tr;
		Extend = ex;
		_medians = new double[dim];
		uninitialized_copy(list.begin(), list.end(), _medians);
	}

	PyramidDB& operator=(const PyramidDB& other)
	{
		this->Table = other.Table;
		(int &)this->nDims = other.nDims;
		this->SRID = other.SRID;
		this->trans = other.trans;
		Extend = other.Extend;
		_medians = other._medians;
		return *this;
	}

};
