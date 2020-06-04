#pragma once
#include "Point.h"

template< typename T>
class NDWindow
{
public:
	NDPoint<T> minPoint;
	NDPoint<T> maxPoint;
	int  nDims;

public:
	NDWindow() :nDims(0)
	{

	}

	NDWindow(const NDPoint<T> & P1, const NDPoint<T> & P2) //, nDims, nDims
	{
		this->nDims = P1.returnSize();
		this->minPoint = P1;
		this->maxPoint = P2;

	}
	

	int GetDimensions()
	{
		return this->nDims;
	}

	NDPoint<T> GetMinPoint() const //, nDims
	{
		return this->minPoint;
	}

	NDPoint<T> GetMaxPoint() const //, nDims
	{
		return this->maxPoint;
	}

	void SetMinPoint(const NDPoint<T> & minpt) //, nDims
	{
		this->minPoint = minpt;
		this->nDims = minpt.returnSize();
	}

	void SetMaxPoint(const NDPoint<T> & maxpt)
	{
		this->maxPoint = maxpt;
		this->nDims = maxpt.returnSize();
	}

	void SetPoints(const NDPoint<T>& minPoint, const NDPoint<T>& maxPoint)
	{
		STATIC_ASSERT(minPoint.returnSize() == maxPoint.returnSize());

		this->minPoint = minPoint;
		this->maxPoint = maxPoint;

		this->nDims = minPoint.returnSize();
	}

	T GetDimWidth(int idx)
	{
		return this->maxPoint[idx] - this->minPoint[idx];
	}

	NDWindow<T> Transform(const CoordTrans<T> & trans) const
	{
		NDWindow<T> windowTrans(minPoint.Transform(trans), maxPoint.Transform(trans));
		return windowTrans;
	}

};
