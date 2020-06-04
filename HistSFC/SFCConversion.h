#pragma once
#include <array>            // std::array
#include "Point.h"

template<typename T>
class SFCConversion
{
private:
	std::vector<unsigned long> g_mask;//nDims
	int nDims;  
	int mBits;

public:
	//Point<long, nDims> ptCoord; //n*m
	//Point<long,  mBits> ptBits; //m*n

	SFCConversion(int dims, int bits) : g_mask(dims, 0), nDims(dims), mBits(bits)
	{
		for (int i = 0; i < nDims; i++)
		{
			g_mask[i] = ((unsigned long)1) << (nDims - i - 1);
		}
	}

	template<typename J = T>
	typename std::enable_if<!std::is_integral<J>::value, sfc_bigint>::type
		MortonEncode(NDPoint<J>& ptCoord)
	{
		throw("Point is not integral type!");
		return 0;
	}

	template<typename J = T>
	typename std::enable_if<std::is_integral<J>::value, sfc_bigint>::type
		MortonEncode(NDPoint<J>& ptCoord)
	{
		NDPoint<J> ptBits(mBits);

		for (int i = 0; i < mBits; i++)//m
		{
			ptBits[i] = 0;
			J mask = ((J)1) << (mBits - i - 1); //move to the ith bit

			for (int j = 0; j < nDims; j++) //get one bit from each nDims
			{
				if (ptCoord[j] & mask) // both 1
					ptBits[i] |= (J)1 << (nDims - j - 1);// push this bit to dim position(zyx...) nDims -----(nDims - j) , << j just for N curve
			}//
		}//m group

		sfc_bigint idx = BitSequence2Value(ptBits);
		return idx;
	}// from n*m coords to m*n bitsequence //, nDims

	template<typename J = T>
	typename std::enable_if<!std::is_integral<J>::value, NDPoint<J>&>::type
		 MortonDecode(sfc_bigint idx)
	{
		throw("Point is not integral type!");
		return 0;
	}


	template<typename J = T>
	typename std::enable_if<std::is_integral<J>::value, NDPoint<J>&>::type
		MortonDecode(sfc_bigint idx)// from m*n bitsequence to n*m coords //, nDims
	{
		NDPoint<J> ptCoord(nDims); //n*m
		NDPoint<J> ptBits = Value2BitSequence(idx); //m*n, mBits

		for (int i = 0; i < nDims; i++)//m n-bits
		{
			ptCoord[i] = 0;
			int mask = ((unsigned int)1) << (nDims - i - 1);

			for (int j = 0; j < mBits; j++)
			{
				if (ptBits[j] & mask) //both 1 
					ptCoord[i] |= (J)1 << (mBits - j - 1); //get the i-th bit from  j-th bits
			}//
		}//n nDims

		return ptCoord;
	}


private:
	sfc_bigint BitSequence2Value(NDPoint<T>& ptBits)
	{
		//if (mBits * nDims >= 64) return 0;
		sfc_bigint  result = 0;
		for (int i = 0; i < mBits; i++)
		{
			if (ptBits[i])
			{
				sfc_bigint a = (((sfc_bigint)ptBits[i]) << (mBits - i - 1)*nDims);
				result |= a;
			}

		}
		return result;
	}


	sfc_bigint BitSequence2Value2(NDPoint<T>& ptBits)
	{
		//if (mBits * nDims >= 64) return 0;
		sfc_bigint  result = 0;
		for (int i = 0; i < mBits; i++)
		{
			if (ptBits[i])
			{
				sfc_bigint a = (((sfc_bigint)ptBits[i]) << i * 8);
				result |= a;
			}

		}
		return result;
	}

	NDPoint<T>& Value2BitSequence(sfc_bigint value)
	{
		NDPoint<T> ptOutput(mBits);
		//if (mBits * nDims >= 64) return ptOutput;

		T mask = (((T)1 << nDims) - 1);
		for (int i = 0; i < mBits; i++)
		{
			ptOutput[mBits - i - 1] = (T)((value >> (i*nDims)) & mask);
		}

		return ptOutput;
	}

	NDPoint<T>& Value2BitSequence2(sfc_bigint value)
	{
		NDPoint<T> ptOutput(mBits);
		//if (mBits * nDims >= 64) return ptOutput;

		long mask = (((T)1 << 8) - 1);
		for (int i = 0; i < mBits; i++)
		{
			ptOutput[i] = (T)((value >> (i * 8)) & mask);
		}

		return ptOutput;
	}
};
