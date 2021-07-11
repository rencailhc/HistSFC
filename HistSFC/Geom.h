/*Half-space, quadratic-space and nD-polytope model*/ 

#pragma once
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
		for (int i = 0; i < dimnum; i++)
			w[i] = 0;
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

	void print_eq() {
		cout << "Equation: ";
		for (int i = 0; i < dimnum-1; i++) {
			cout << to_string(w[i]) << "x" << to_string(i) << " + ";
		}
		cout << to_string(w[dimnum - 1]) << "x" << to_string(dimnum - 1);
		switch (s)
		{
		case Sign::le:
			cout << " <= ";
			break;
		case Sign::ge:
			cout << " >= ";
			break;
		case Sign::eq:
			cout << " = ";
			break;
		default:
			break;
		}

		cout << b << endl;
	}

};


class quadspace
{
	//Quadratic space
	//ω ⋅ x^2 + ω ⋅ x - β <= or >= or = 0,
	//w   list of float(ω)
	//b   float(β)
public:
	short dimnum;
	double* w;
	double b;
	Sign s;
public:
	quadspace() :dimnum(0)
	{
		w = nullptr;
		b = 0;
		s = Sign::eq;
	}

	quadspace(short dims) :dimnum(dims)
	{
		w = new double[dimnum*2]; //also include quadratic terms
		for (int i = 0; i < dimnum*2; i++)
			w[i] = 0;
		b = 0;
		s = Sign::eq;
	}

	quadspace(std::initializer_list<double> l1, double B, Sign S) :dimnum((int)(l1.size()/2))
	{
		w = new double[dimnum*2];
		std::uninitialized_copy(l1.begin(), l1.end(), w);
		b = B;
		s = S;
	}

	quadspace& operator=(const quadspace& other)
	{
		//assignment by copy
		(short &)dimnum = other.dimnum;
		delete[] w;
		b = other.b;
		s = other.s;
		w = new double[dimnum*2];
		for (int i = 0; i < dimnum*2; i++)
		{
			w[i] = other.w[i];
		}

		return *this;
	}

	void print_eq() {
		cout << "Equation: ";
		for (int i = 0; i < dimnum - 1; i++) {
			cout << to_string(w[i * 2]) << "x" << to_string(i) << "^2" << " + " << to_string(w[i * 2 + 1]) << "x" << to_string(i) << " + ";
		}
		cout << to_string(w[2 * dimnum - 2]) << "x" << to_string(dimnum - 1) << "^2" << " + " << to_string(w[2 * dimnum - 1]) << "x" << to_string(dimnum - 1);
		switch (s)
		{
		case Sign::le:
			cout << " <= ";
			break;
		case Sign::ge:
			cout << " >= ";
			break;
		case Sign::eq:
			cout << " = ";
			break;
		default:
			break;
		}

		cout << b << endl;
	}

};

class NDGeom {
public:
	short dimnum;
	std::vector <halfspace> faces;
	std::vector <quadspace> curves;
public:
	NDGeom()
	{
		dimnum = 0;
		faces = {};
		curves = {};
	}

	NDGeom(short dims):dimnum(dims)
	{
		faces = {};
		curves = {};
	}

	void add(halfspace h) {
		if (!dimnum) dimnum = h.dimnum;
		if (dimnum != h.dimnum) throw "Dimensionality does not match!";
		faces.push_back(h);
	}

	void add(quadspace c) {
		if (!dimnum) dimnum = c.dimnum;
		if (dimnum != c.dimnum) throw "Dimensionality does not match!";
		curves.push_back(c);
	}

	int getSize() const
	{
		return faces.size() + curves.size();
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

		for (auto it = curves.begin(); it != curves.end(); it++)
		{
			quadspace c(dimnum);
			c = *it;
			for (int i = 0; i < dimnum; i++)
			{
				c.w[i * 2] = it->w[i * 2] / (trans._scale[i] * trans._scale[i]);
				c.w[i * 2 + 1] = it->w[i * 2 + 1] / trans._scale[i] + 2 * it->w[i * 2] * trans._delta[i] / trans._scale[i];
				c.b -= it->w[2 * i + 1] * trans._delta[i] + it->w[2 * i] * trans._delta[i] * trans._delta[i];
			}
			geomtrans.add(c);
		}

		return geomtrans;
	}

	NDGeom & Normalize()
	{
		for (auto it = faces.begin(); it != faces.end(); it++)
		{
			(*it).Normalize();
		}
		return *this;
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

		curves.clear();
		for (auto it = other.curves.begin(); it != other.curves.end(); it++)
		{
			curves.push_back(*it);
		}

		return *this;
	}

	void print_eq() {
		for (int i = 0; i < faces.size(); i++)
		{
			faces[i].print_eq();
		}

		for (int i = 0; i < curves.size(); i++)
		{
			curves[i].print_eq();
		}
	}
};
