#pragma once

#include "boost/multiprecision/cpp_int.hpp"
#include "BaseStruct.h"
#include <array>
#include <initializer_list>

using namespace boost::multiprecision;

#define DIM_MAX 40
#define PREC (short) 8	//the precision of the pyramid key, depending on how many points there are

template<typename T>
class NDPoint
{
private:
	std::array< T, DIM_MAX > elements_; 
	int nDims;

public:
	
	NDPoint() :nDims(0)
	{

	}

	NDPoint(std::size_t dims)
	{
		nDims = dims;
		for (int i = 0; i < nDims; ++i)
		{
			elements_[i] = 0;
		}
	}

	NDPoint(T *coordinates, int dims)
	{
		nDims = dims;
		for (int i = 0; i < this->nDims; i++)
		{
			elements_[i] = *(coordinates + i);
		}
	}

	NDPoint(std::initializer_list<T> list)
	{
		nDims = (int)list.size();
		uninitialized_copy(list.begin(), list.end(), elements_.begin());
	}

	int returnSize() const
	{
		//return this->elements_.size();
		return nDims;
	}


	T& operator[](int const i)
	{
		return elements_[i];
	}

	T const& operator[](int const i) const
	{
		return elements_[i];
	}

	void operator+=(NDPoint<T> const& other)
	{
		if (nDims != other.returnSize()) throw "Size does not match!";

		for (int i = 0; i < nDims; ++i)
		{
			elements_[i] += other.elements_[i];
		}
	}


	void operator=(const NDPoint<T> & other)
	{
		nDims = other.returnSize();

		for (int i = 0; i < nDims; ++i)
		{
			elements_[i] = other.elements_[i];
		}
	}

	void operator-=(NDPoint<T> const& other)
	{
		if (nDims != other.returnSize()) throw "Size does not match!";

		for (int i = 0; i < nDims; ++i)
		{
			elements_[i] -= other.elements_[i];
		}
	}

	template <typename U>
	NDPoint<U> Transform(const CoordTrans<U>& trans) const 
	{
		NDPoint<U> outPt(nDims);
		for (int i = 0; i < nDims; i++)
		{
			outPt[i] = (elements_[i] - trans._delta[i])*trans._scale[i];
		}
		return outPt;
	}

	template <typename U>
	NDPoint<U> InverseTransform(const CoordTrans<U>& trans) const
	{
		NDPoint<U> outPt(nDims);
		for (int i = 0; i < nDims; i++)
		{
			//outPt[i] = lround((inPt[i] - _delta[i])*_scale[i]);
			outPt[i] = ((double)elements_[i]) / trans._scale[i] + trans._delta[i]; //decoding
		}
		return outPt;
	}
};