﻿#pragma once
#include "typedef.h"
#include <initializer_list>

enum class Sign {
	//<=, >=,  =
	le, ge, eq
};

class halfspace
{
	//Any hyperplane can be written as the set of points x satisfying
	//ω ⋅ x - β <= or >= or = 0,
	//where ω is the(not necessarily normalized) normal vector to the hyperplane.

	//w   list of float(ω)
	//b   float(β)
public:
	short dimnum;
	double* w;
	double b;
	Sign s;
public:
	halfspace() :dimnum(0)
	{
		w = nullptr;
		b = 0;
		s = Sign::eq;
	}

	halfspace(short dims) :dimnum(dims)
	{
		w = new double[dimnum];
		b = 0;
		s = Sign::eq;
	}

	halfspace(std::initializer_list<double> l1, double B, Sign S) :dimnum((int)l1.size())
	{
		w = new double[dimnum];
		std::uninitialized_copy(l1.begin(), l1.end(), w);
		b = B;
		s = S;
	}

	halfspace & Normalize() {
		double base = 0;
		for (int i = 0; i < dimnum; i++)
		{
			base += w[i] * w[i];
		}
		base = pow(base, 0.5);

		for (int i = 0; i < dimnum; i++)
		{
			w[i] = w[i] / base;
		}
		b = b / base;

		if (s == Sign::le) {
			for (int i = 0; i < dimnum; i++)
			{
				w[i] = -w[i];
			}
			b = -b;
			s = Sign::ge;
		}

		return *this;
	}

	double * Normal() {
		halfspace copy = *this;
		copy.Normalize();
		return copy.w;
	}

	halfspace& operator=(const halfspace& other)
	{
		//assignment by copy
		(short &)dimnum = other.dimnum;
		delete[] w;
		b = other.b;
		s = other.s;
		w = new double[dimnum];
		for (int i = 0; i < dimnum; i++)
		{
			w[i] = other.w[i];
		}

		return *this;
	}

};

class NDGeom {
public:
	short dimnum;
	std::vector <halfspace> faces;
public:
	NDGeom()
	{
		dimnum = 0;
		faces = {};
	}

	NDGeom(short dims):dimnum(dims)
	{
		faces = {};
	}

	void add(halfspace h) {
		if (!dimnum) dimnum = h.dimnum;
		if (dimnum != h.dimnum) throw "Dimensionality does not match!";
		faces.push_back(h);
	}

	int getSize() const
	{
		return faces.size();
	}

	NDGeom Transform(const CoordTrans & trans) const
	{
		NDGeom geomtrans(dimnum);
		for (auto it = faces.begin(); it != faces.end(); it++)
		{
			halfspace h(dimnum);
			h = *it;
			for (int i = 0; i < dimnum; i++)
			{
				h.w[i] = it->w[i] / trans._scale[i];
				h.b -= it->w[i] * trans._delta[i];
			}
			geomtrans.add(h);
		}

		return geomtrans;
	}

	NDGeom & Normalize()
	{
		for (auto it = faces.begin(); it != faces.end(); it++)
		{
			double base = 0;
			for (int i = 0; i < dimnum; i++)
			{
				base += it->w[i] * it->w[i];
			}
			base = pow(base, 0.5);

			for (int i = 0; i < dimnum; i++)
			{
				it->w[i] = it->w[i] / base;
			}
			it->b = it->b / base;

			if (it->s == Sign::le) {
				for (int i = 0; i < dimnum; i++)
				{
					it->w[i] = -it->w[i];
				}
				it->b = -it->b;
				it->s = Sign::ge;
			}
		}
		return *this;
	}

	halfspace& operator[](int const i)
	{
		return faces[i];
	}

	NDGeom& operator=(const NDGeom& other)
	{
		//assignment by copy
		(short &)dimnum = other.dimnum;
		faces.clear();
		for (auto it=other.faces.begin();it!=other.faces.end();it++)
		{
			faces.push_back(*it);
		}

		return *this;
	}
};