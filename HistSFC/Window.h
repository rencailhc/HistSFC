/*Defining nD query window*/

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

	template <typename U>
	NDWindow<U> Transform(const CoordTrans & trans) const
	{
		NDWindow<U> windowTrans(minPoint.Transform<U>(trans), maxPoint.Transform<U>(trans));
		for (int i = 0; i < nDims; i++)
		{
			if (windowTrans.minPoint[i] > windowTrans.maxPoint[i])
			{
				U vtmp = windowTrans.minPoint[i];
				windowTrans.minPoint[i] = windowTrans.maxPoint[i];
				windowTrans.maxPoint[i] = vtmp;
			}
		}
		return windowTrans;
	}

	template <typename U>
	NDWindow<U> InverseTransform(const CoordTrans& trans) const
	{
		NDWindow<U> windowInvTrans(minPoint.InverseTransform<U>(trans), maxPoint.InverseTransform<U>(trans));
		for (int i = 0; i < nDims; i++)
		{
			if (windowInvTrans.minPoint[i] > windowInvTrans.maxPoint[i])
			{
				U vtmp = windowInvTrans.minPoint[i];
				windowInvTrans.minPoint[i] = windowInvTrans.maxPoint[i];
				windowInvTrans.maxPoint[i] = vtmp;
			}
		}
		return windowInvTrans;
	}

	void print() const
	{
		for (int i = 0; i < nDims-1; i++)
		{
			std::cout << std::fixed<< std::setprecision(3)<<minPoint[i] << ",";
		}
		std::cout << minPoint[nDims - 1] << "; ";
		for (int i = 0; i < nDims-1; i++)
		{
			std::cout << maxPoint[i] << ",";
		}
		std::cout << maxPoint[nDims - 1] << std::endl;
	}

};
