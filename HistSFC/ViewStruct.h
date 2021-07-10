#pragma once
#include "BaseStruct.h"
#include "Geom.h"


/*data stuctures and operators for 4D perspective view selection*/
namespace viewlib
{
	struct viewPos
	{
		double P[3];
		double dir[3];
		double distance;
	};

	struct Cone4D
	{
		//expression: (x-u)^2+(y-v)^2+(z-w)^2 < (a-b*LoD)^2
		double P[3];
		double a;
		double b;
		double maxl;  //max lod value, min cLoD is 0
	};

	typedef struct Cube3D
	{
		double minP[3];
		double maxP[3];
	} Cube3D;

	struct Sphere3D
	{
		double o[3];
		double radius;
	};

	struct Cone3D
	{
		//equation: (x-u)^2 + (y-v)^2 + c <= (a-b*LoD)^2
		double P[2];
		double c;
		double a;
		double b;
		double maxl; //maxLoD
	};

	struct Plane3D
	{
		/*
			A 3d plane
			  v
			  ^----------+
			  |          |
			  |          |
			  +----------> u
			  P
		*/
		double P[3];
		double u[3];
		double v[3];
	};

	struct Line3D
	{
		double P[3];
		double dir[3];
	};


	struct Box3D
	{
		Plane3D faces[6];
	};


	struct Frustum3D
	{
		double P[3];    //vertex
		double CP[4][3];    //conner vertex, ordered counter-clockwise
		double e[6][3];   //6 edges out of 8, except 2 unnecessary at the bottome face
		double normal[5][3];   //5 normals
	};


	Frustum3D FrustumBuild(viewPos params)
	{
		Frustum3D F;
		for (int i = 0; i < 3; i++) F.P[i] = params.P[i];

		double right[3];
		if (params.dir[1])
		{
			double s_x = sqrt(1 / (params.dir[0] * params.dir[0] / (params.dir[1] * params.dir[1]) + 1));
			double s_y = -params.dir[0] * s_x / params.dir[1];
			double side1[] = { s_x, s_y, 0 };
			double side2[] = { -s_x, -s_y, 0 };

			if (s_x*params.dir[1] - s_y * params.dir[0] > 0)
			{
				for (int i = 0; i < 3; i++) right[i] = side1[i];
			}
			else
			{
				for (int i = 0; i < 3; i++) right[i] = side2[i];
			}
		}
		else
		{
			right[0] = 0;
			right[2] = 0;
			if (params.dir[0] > 0) 	right[1] = -1;
			else right[1] = 1;
		}

		double up[] = { right[1] * params.dir[2], -right[0] * params.dir[2], right[0] * params.dir[1] - right[1] * params.dir[0] };
		double fc[] = { params.P[0] + params.dir[0] * params.distance, params.P[1] + params.dir[1] * params.distance,
			params.P[2] + params.dir[2] * params.distance };
		for (int i = 0; i < 3; i++)
		{
			F.P[i] = params.P[i];
			F.CP[0][i] = fc[i] + 2 * up[i] * params.distance - 2 * right[i] * params.distance;
			F.CP[1][i] = fc[i] + 2 * up[i] * params.distance + 2 * right[i] * params.distance;
			F.CP[2][i] = fc[i] - 2 * up[i] * params.distance + 2 * right[i] * params.distance;
			F.CP[3][i] = fc[i] - 2 * up[i] * params.distance - 2 * right[i] * params.distance;
			F.e[0][i] = F.CP[0][i] - F.P[i];
			F.e[1][i] = F.CP[1][i] - F.P[i];
			F.e[2][i] = F.CP[2][i] - F.P[i];
			F.e[3][i] = F.CP[3][i] - F.P[i];
			F.e[4][i] = F.e[0][i] - F.e[1][i];
			F.e[5][i] = F.e[2][i] - F.e[1][i];
		}

		F.normal[0][0] = F.e[0][1] * F.e[1][2] - F.e[0][2] * F.e[1][1];
		F.normal[0][1] = F.e[0][2] * F.e[1][0] - F.e[0][0] * F.e[1][2];
		F.normal[0][2] = F.e[0][0] * F.e[1][1] - F.e[0][1] * F.e[1][0];
		F.normal[1][0] = F.e[1][1] * F.e[2][2] - F.e[1][2] * F.e[2][1];
		F.normal[1][1] = F.e[1][2] * F.e[2][0] - F.e[1][0] * F.e[2][2];
		F.normal[1][2] = F.e[1][0] * F.e[2][1] - F.e[1][1] * F.e[2][0];
		F.normal[2][0] = F.e[2][1] * F.e[3][2] - F.e[2][2] * F.e[3][1];
		F.normal[2][1] = F.e[2][2] * F.e[3][0] - F.e[2][0] * F.e[3][2];
		F.normal[2][2] = F.e[2][0] * F.e[3][1] - F.e[2][1] * F.e[3][0];
		F.normal[3][0] = F.e[3][1] * F.e[0][2] - F.e[3][2] * F.e[0][1];
		F.normal[3][1] = F.e[3][2] * F.e[0][0] - F.e[3][0] * F.e[0][2];
		F.normal[3][2] = F.e[3][0] * F.e[0][1] - F.e[3][1] * F.e[0][0];
		F.normal[4][0] = F.e[4][1] * F.e[5][2] - F.e[4][2] * F.e[5][1];
		F.normal[4][1] = F.e[4][2] * F.e[5][0] - F.e[4][0] * F.e[5][2];
		F.normal[4][2] = F.e[4][0] * F.e[5][1] - F.e[4][1] * F.e[5][0];

		return F;
	}


	Box3D Cube2Box(Cube3D R)
	{
		double xdimwidth = R.maxP[0] - R.minP[0];
		double ydimwidth = R.maxP[1] - R.minP[1];
		double loddimwidth = R.maxP[2] - R.minP[2];
		Plane3D face0 = { {R.minP[0], R.minP[1], R.minP[2]}, {xdimwidth,0,0}, {0,ydimwidth,0} };
		Plane3D face1 = { {R.minP[0], R.minP[1], R.minP[2]}, {0,0,loddimwidth}, {0,ydimwidth,0} };
		Plane3D face2 = { {R.minP[0], R.minP[1], R.minP[2]}, {xdimwidth,0,0}, {0,0,loddimwidth} };
		Plane3D face3 = { {R.maxP[0], R.maxP[1], R.maxP[2]}, {-xdimwidth,0,0}, {0,-ydimwidth,0} };
		Plane3D face4 = { {R.maxP[0], R.maxP[1], R.maxP[2]}, {0,0,-loddimwidth}, {0,-ydimwidth,0} };
		Plane3D face5 = { {R.maxP[0], R.maxP[1], R.maxP[2]}, {-xdimwidth,0,0}, {0,0,-loddimwidth} };
		Box3D B = { face0, face1, face2, face3, face4, face5 };
		return B;
	}


	/*Detect intersection between a 3D line segment and a 3D rectangle*/
	bool intersectLF3D(Line3D s, Plane3D face)
	{
		double a = face.u[0];
		double b = face.v[0];
		double c = -s.dir[0];
		double d = s.P[0] - face.P[0];
		double e = face.u[1];
		double f = face.v[1];
		double g = -s.dir[1];
		double h = s.P[1] - face.P[1];
		double i = face.u[2];
		double j = face.v[2];
		double k = -s.dir[2];
		double l = s.P[2] - face.P[2];

		double delta = (a*f*k) + (b*g*i) + (c*e*j) - (c*f*i) - (a*g*j) - (b*e*k);

		if (delta)
		{
			double v1 = ((d*f*k) + (b*g*l) + (c*h*j) - (c*f*l) - (d*g*j) - (b*h*k)) / delta;
			double v2 = ((a*h*k) + (d*g*i) + (c*e*l) - (c*h*i) - (a*g*l) - (d*e*k)) / delta;
			double v3 = ((a*f*l) + (b*h*i) + (d*e*j) - (d*f*i) - (a*h*j) - (b*e*l)) / delta;
			return v1 >= 0 && v1 <= 1 && v2 >= 0 && v2 <= 1 && v3 >= 0 && v3 <= 1;
		}
		else return false;

	}

	/*Detect intersection between a 3D line segment and a 3D triangle*/
	bool intersectLT3D(Line3D s, Plane3D triangle)
	{
		double a = triangle.u[0];
		double b = triangle.v[0];
		double c = -s.dir[0];
		double d = s.P[0] - triangle.P[0];
		double e = triangle.u[1];
		double f = triangle.v[1];
		double g = -s.dir[1];
		double h = s.P[1] - triangle.P[1];
		double i = triangle.u[2];
		double j = triangle.v[2];
		double k = -s.dir[2];
		double l = s.P[2] - triangle.P[2];

		double delta = (a*f*k) + (b*g*i) + (c*e*j) - (c*f*i) - (a*g*j) - (b*e*k);

		if (delta != 0)
		{
			double v1 = ((d*f*k) + (b*g*l) + (c*h*j) - (c*f*l) - (d*g*j) - (b*h*k)) / delta;
			double v2 = ((a*h*k) + (d*g*i) + (c*e*l) - (c*h*i) - (a*g*l) - (d*e*k)) / delta;
			double v3 = ((a*f*l) + (b*h*i) + (d*e*j) - (d*f*i) - (a*h*j) - (b*e*l)) / delta;
			return v1 >= 0 && v1 <= 1 && v2 >= 0 && v2 <= 1 - v1 && v3 >= 0 && v3 <= 1;
		}
		else return false;

	}


	/*Detect intersection between a 3D triangle and a 3D face*/
	bool intersectTF3D(Plane3D triangle, Plane3D face)
	{
		Line3D e1, e2, e3, e4;
		for (int i = 0; i < 3; i++)
		{
			e1.P[i] = triangle.P[i];
			e1.dir[i] = triangle.u[i];
			e2.P[i] = triangle.P[i];
			e2.dir[i] = triangle.v[i];
			e3.P[i] = triangle.P[i] + triangle.u[i];
			e3.dir[i] = triangle.v[i] - triangle.u[i];
		}

		bool intersect1 = intersectLF3D(e1, face) || intersectLF3D(e2, face) || intersectLF3D(e3, face);

		//reverse check
		for (int i = 0; i < 3; i++)
		{
			e1.P[i] = face.P[i];
			e1.dir[i] = face.u[i];
			e2.P[i] = face.P[i];
			e2.dir[i] = face.v[i];
			e3.P[i] = face.P[i] + face.u[i] + face.v[i];
			e3.dir[i] = -face.u[i];
			e4.P[i] = face.P[i] + face.u[i] + face.v[i];
			e4.dir[i] = -face.v[i];
		}

		bool intersect2 = intersectLT3D(e1, triangle) || intersectLT3D(e2, triangle) || intersectLT3D(e3, triangle) || intersectLT3D(e4, triangle);

		return intersect1 || intersect2;
	}

	/*Detect intersection between 2 3D rectangles*/
	bool intersectFF3D(Plane3D face1, Plane3D face2)
	{
		Line3D e1, e2, e3, e4;
		for (int i = 0; i < 3; i++)
		{
			e1.P[i] = face1.P[i];
			e1.dir[i] = face1.u[i];
			e2.P[i] = face1.P[i];
			e2.dir[i] = face1.v[i];
			e3.P[i] = face1.P[i] + face1.u[i] + face1.v[i];
			e3.dir[i] = -face1.u[i];
			e4.P[i] = face1.P[i] + face1.u[i] + face1.v[i];
			e4.dir[i] = -face1.v[i];
		}

		bool intersect1 = intersectLF3D(e1, face2) || intersectLF3D(e2, face2) || intersectLF3D(e3, face2) || intersectLF3D(e4, face2);

		//reverse check
		for (int i = 0; i < 3; i++)
		{
			e1.P[i] = face2.P[i];
			e1.dir[i] = face2.u[i];
			e2.P[i] = face2.P[i];
			e2.dir[i] = face2.v[i];
			e3.P[i] = face2.P[i] + face2.u[i] + face2.v[i];
			e3.dir[i] = -face2.u[i];
			e4.P[i] = face2.P[i] + face2.u[i] + face2.v[i];
			e4.dir[i] = -face2.v[i];
		}

		bool intersect2 = intersectLF3D(e1, face1) || intersectLF3D(e2, face1) || intersectLF3D(e3, face1) || intersectLF3D(e4, face1);

		return intersect1 || intersect2;

	}

	/*Detect intersection between a 3D sphere and a 3D cube*/
	bool intersectSV3D(Sphere3D ball, Cube3D cube)
	{
		//referencing https://stackoverflow.com/questions/4578967/cube-sphere-intersection-test
		double dist_squared = ball.radius * ball.radius;
		if (ball.o[0] < cube.minP[0]) dist_squared -= (ball.o[0] - cube.minP[0]) * (ball.o[0] - cube.minP[0]);
		else if (ball.o[0] > cube.maxP[0]) dist_squared -= (ball.o[0] - cube.maxP[0]) * (ball.o[0] - cube.maxP[0]);
		if (ball.o[1] < cube.minP[1]) dist_squared -= (ball.o[1] - cube.minP[1]) * (ball.o[1] - cube.minP[1]);
		else if (ball.o[1] > cube.maxP[1]) dist_squared -= (ball.o[1] - cube.maxP[1]) * (ball.o[1] - cube.maxP[1]);
		if (ball.o[2] < cube.minP[2]) dist_squared -= (ball.o[2] - cube.minP[2]) * (ball.o[2] - cube.minP[2]);
		else if (ball.o[2] > cube.maxP[2]) dist_squared -= (ball.o[2] - cube.maxP[2]) * (ball.o[2] - cube.maxP[2]);
		return dist_squared >= 0;
	}

	/*Detect intersection between a 3D line segment and a 3D cone*/
	bool intersectLC3D(Line3D l, Cone3D cone)
	{
		// Beware, indigest formulaes !
		double A = l.dir[0] * l.dir[0] + l.dir[1] * l.dir[1] - l.dir[2] * l.dir[2] * cone.b*cone.b;
		double B = 2 * (l.P[0] - cone.P[0]) * l.dir[0] + 2 * (l.P[1] - cone.P[1]) * l.dir[1] + 2 * (cone.a - cone.b*l.P[2]) * l.dir[2] * cone.b;
		double C = (l.P[0] - cone.P[0]) * (l.P[0] - cone.P[0]) + (l.P[1] - cone.P[1]) * (l.P[1] - cone.P[1]) + cone.c - (cone.a - cone.b*l.P[2]) * (cone.a - cone.b*l.P[2]);

		// Now, we solve the polynom At² + Bt + C = 0
		double delta = B * B - 4 * A * C;
		if (delta < 0)
			return false; // No intersection between the cone and the line
		else if (A != 0)
		{
			// Check the two solutions (there might be only one, but that does not change a lot of things)
			double t1 = (-B + sqrt(delta)) / (2 * A);
			double z1 = l.P[2] + t1 * l.dir[2];
			bool t1_intersect = (t1 >= 0 && t1 <= 1 && z1 >= 0 && z1 <= cone.maxl);

			double t2 = (-B - sqrt(delta)) / (2 * A);
			double z2 = l.P[2] + t2 * l.dir[2];
			bool t2_intersect = (t2 >= 0 && t2 <= 1 && z2 >= 0 && z2 <= cone.maxl);

			return t1_intersect || t2_intersect;
		}
		else if (B != 0)
		{
			double t = -C / B;
			double z = l.P[2] + t * l.dir[2];
			return t >= 0 && t <= 1 && z >= 0 && z <= cone.maxl;
		}
		else return C == 0;
	}


	/*Detect intersection between the 3D xyz view frustum and 3D xyz node*/
	template <typename T>
	short intersectFV3D(Frustum3D view, NDWindow<T> node)
	{
		//2 indicates inside, 1 indicates intersection
		short ncmp = 1;

		//view inside cube
		for (int i = 0; i < 3; i++)
		{
			ncmp &= view.P[i] >= node.minPoint[i] && view.P[i] <= node.maxPoint[i];
		}
		if (ncmp) return 1;

		//Voxel inside view
		ncmp = 1;

		double D;
		for (int i = 0; i < 5; i++)
		{
			if (i < 4) D = -view.P[0] * view.normal[i][0] - view.P[1] * view.normal[i][1] - view.P[2] * view.normal[i][2];
			else D = -view.CP[1][0] * view.normal[i][0] - view.CP[1][1] * view.normal[i][1] - view.CP[1][2] * view.normal[i][2];
			ncmp &= node.minPoint[0] * view.normal[i][0] + node.minPoint[1] * view.normal[i][1] + node.minPoint[2] * view.normal[i][2] + D > 0 &&
				node.minPoint[0] * view.normal[i][0] + node.maxPoint[1] * view.normal[i][1] + node.minPoint[2] * view.normal[i][2] + D > 0 &&
				node.minPoint[0] * view.normal[i][0] + node.maxPoint[1] * view.normal[i][1] + node.maxPoint[2] * view.normal[i][2] + D > 0 &&
				node.minPoint[0] * view.normal[i][0] + node.minPoint[1] * view.normal[i][1] + node.minPoint[2] * view.normal[i][2] + D > 0 &&
				node.maxPoint[0] * view.normal[i][0] + node.maxPoint[1] * view.normal[i][1] + node.minPoint[2] * view.normal[i][2] + D > 0 &&
				node.maxPoint[0] * view.normal[i][0] + node.minPoint[1] * view.normal[i][1] + node.minPoint[2] * view.normal[i][2] + D > 0 &&
				node.maxPoint[0] * view.normal[i][0] + node.minPoint[1] * view.normal[i][1] + node.maxPoint[2] * view.normal[i][2] + D > 0 &&
				node.maxPoint[0] * view.normal[i][0] + node.maxPoint[1] * view.normal[i][1] + node.maxPoint[2] * view.normal[i][2] + D > 0;
			if (!ncmp) break;
		}

		if (ncmp) return 2;

		//normal intersection
		Plane3D triangle1 = { {view.P[0],view.P[1],view.P[2]},{view.e[0][0],view.e[0][1],view.e[0][2]},{view.e[1][0],view.e[1][1],view.e[1][2]} };
		Plane3D triangle2 = { {view.P[0],view.P[1],view.P[2]},{view.e[1][0],view.e[1][1],view.e[1][2]},{view.e[2][0],view.e[2][1],view.e[2][2]} };
		Plane3D triangle3 = { {view.P[0],view.P[1],view.P[2]},{view.e[2][0],view.e[2][1],view.e[2][2]},{view.e[3][0],view.e[3][1],view.e[3][2]} };
		Plane3D triangle4 = { {view.P[0],view.P[1],view.P[2]},{view.e[3][0],view.e[3][1],view.e[3][2]},{view.e[0][0],view.e[0][1],view.e[0][2]} };
		Plane3D bottom = { {view.CP[1][0],view.CP[1][1],view.CP[1][2]},{view.e[4][0],view.e[4][1],view.e[4][2]},{view.e[5][0],view.e[5][1],view.e[5][2]} };

		double arrmin[] = { node.minPoint[0], node.minPoint[1], node.minPoint[2] };
		double arrmax[] = { node.maxPoint[0], node.maxPoint[1], node.maxPoint[2] };

		Box3D box = Cube2Box(Cube3D({ { arrmin[0], arrmin[1], arrmin[2] }, { arrmax[0], arrmax[1], arrmax[2] } }));

		for (int i = 0; i < 6; i++)
		{
			ncmp |= intersectTF3D(triangle1, box.faces[i]) ||
				intersectTF3D(triangle2, box.faces[i]) ||
				intersectTF3D(triangle3, box.faces[i]) ||
				intersectTF3D(triangle4, box.faces[i]) ||
				intersectFF3D(bottom, box.faces[i]);
			if (ncmp) return 1;
		}

		return 0;
	}


	/*Detect intersection between a 3D cone and a 3D rectangle*/
	bool intersectCF3D(Cone3D cone, Plane3D rect)
	{
		Line3D l1, l2, l3, l4;
		for (int i = 0; i < 3; i++)
		{
			l1.P[i] = rect.P[i];
			l1.dir[i] = rect.u[i];
			l2.P[i] = rect.P[i];
			l2.dir[i] = rect.v[i];
			l3.P[i] = rect.P[i] + rect.u[i];
			l3.dir[i] = rect.v[i];
			l4.P[i] = rect.P[i] + rect.v[i];
			l4.dir[i] = rect.u[i];
		}
		bool intersection = intersectLC3D(l1, cone)
			|| intersectLC3D(l2, cone)
			|| intersectLC3D(l3, cone)
			|| intersectLC3D(l4, cone);

		if (!intersection)
		{
			// It is possible that either the part of the plan lie
			// entirely in the cone, or the inverse. We need to check.

			if (rect.P[2] >= 0 && rect.P[2] <= cone.maxl)
			{
				if ((rect.P[0] - cone.P[0])*(rect.P[0] - cone.P[0]) + (rect.P[1] - cone.P[1])*(rect.P[1] - cone.P[1]) + cone.c <= (cone.a - cone.b*rect.P[2])*(cone.a - cone.b*rect.P[2]))
					return true;
			}

			// Is the cone inside the face (this one is more tricky) ?
			// It can be resolved by finding whether the axis of the cone crosses the face.
			// First, find the plane coefficient (descartes equation)
			//Equation: p = P + a*v + b*u
			double a = (cone.P[0] * rect.u[1] - cone.P[1] * rect.u[0] + rect.P[1] * rect.u[0] - rect.P[0] * rect.u[1]) / (rect.u[1] * rect.v[0] - rect.u[0] * rect.v[1]);
			double b = (cone.P[1] * rect.v[0] - cone.P[0] * rect.v[1] + rect.P[0] * rect.v[1] - rect.P[1] * rect.v[0]) / (rect.u[1] * rect.v[0] - rect.u[0] * rect.v[1]);
			double z = rect.P[2] + a * rect.v[2] + b * rect.u[2];
			intersection = (a >= 0 && a <= 1 && b >= 0 && b <= 1 && z >= 0 && z <= cone.maxl);
			if (intersection) return true;

			//the face intersects the bottom circle of the cone
			double a1 = rect.u[0] - rect.u[2] * rect.v[0] / rect.v[2];
			double a2 = rect.u[1] - rect.u[2] * rect.v[1] / rect.v[2];
			double b1 = rect.P[0] - cone.P[0] - rect.P[2] * rect.v[0] / rect.v[2];
			double b2 = rect.P[1] - cone.P[1] - rect.P[2] * rect.v[1] / rect.v[2];

			double A = a1 * a1 + a2 * a2;
			double B = 2 * a1*b1 + 2 * a2*b2;
			double C = b1 * b1 + b2 * b2 + cone.c - cone.a *cone.a;

			// Now, we solve the polynom At² + Bt + C = 0
			double delta = B * B - 4 * A * C;
			if (delta < 0)
				return false; // No intersection between the cone and the face
			else if (A != 0)
			{
				// Check the two solutions (there might be only one, but that does not change a lot of things)
				double t11 = (-B + sqrt(delta)) / (2 * A);
				double t12 = (-rect.P[2] - t11 * rect.u[2]) / rect.v[2];
				bool t1_intersect = (t11 >= 0 && t11 <= 1 && t12 >= 0 && t12 <= 1);

				double t21 = (-B - sqrt(delta)) / (2 * A);
				double t22 = (-rect.P[2] - t21 * rect.u[2]) / rect.v[2];
				bool t2_intersect = (t21 >= 0 && t21 <= 1 && t22 >= 0 && t22 <= 1);

				return t1_intersect || t2_intersect;
			}
			else
			{
				return false;
			}

		}
		return intersection;
	}

	bool intersectCV3D(Cone3D cone, Cube3D cube)
	{
		if (cube.minP[3] >= 0 && cube.maxP[3] <= cone.maxl &&
			(cube.minP[0] - cone.P[0])*(cube.minP[0] - cone.P[0]) + (cube.minP[1] - cone.P[1])*(cube.minP[1] - cone.P[1]) + cone.c <= (cone.a - cone.b*cube.minP[2])*(cone.a - cone.b*cube.minP[2]) &&
			(cube.minP[0] - cone.P[0])*(cube.minP[0] - cone.P[0]) + (cube.maxP[1] - cone.P[1])*(cube.maxP[1] - cone.P[1]) + cone.c <= (cone.a - cone.b*cube.minP[2])*(cone.a - cone.b*cube.minP[2]) &&
			(cube.minP[0] - cone.P[0])*(cube.minP[0] - cone.P[0]) + (cube.minP[1] - cone.P[1])*(cube.minP[1] - cone.P[1]) + cone.c <= (cone.a - cone.b*cube.maxP[2])*(cone.a - cone.b*cube.maxP[2]) &&
			(cube.minP[0] - cone.P[0])*(cube.minP[0] - cone.P[0]) + (cube.maxP[1] - cone.P[1])*(cube.maxP[1] - cone.P[1]) + cone.c <= (cone.a - cone.b*cube.maxP[2])*(cone.a - cone.b*cube.maxP[2]) &&
			(cube.maxP[0] - cone.P[0])*(cube.maxP[0] - cone.P[0]) + (cube.minP[1] - cone.P[1])*(cube.minP[1] - cone.P[1]) + cone.c <= (cone.a - cone.b*cube.minP[2])*(cone.a - cone.b*cube.minP[2]) &&
			(cube.maxP[0] - cone.P[0])*(cube.maxP[0] - cone.P[0]) + (cube.maxP[1] - cone.P[1])*(cube.maxP[1] - cone.P[1]) + cone.c <= (cone.a - cone.b*cube.minP[2])*(cone.a - cone.b*cube.minP[2]) &&
			(cube.maxP[0] - cone.P[0])*(cube.maxP[0] - cone.P[0]) + (cube.minP[1] - cone.P[1])*(cube.minP[1] - cone.P[1]) + cone.c <= (cone.a - cone.b*cube.maxP[2])*(cone.a - cone.b*cube.maxP[2]) &&
			(cube.maxP[0] - cone.P[0])*(cube.maxP[0] - cone.P[0]) + (cube.maxP[1] - cone.P[1])*(cube.maxP[1] - cone.P[1]) + cone.c <= (cone.a - cone.b*cube.maxP[2])*(cone.a - cone.b*cube.maxP[2])) return true;
		else
		{
			Box3D box = Cube2Box(cube);
			return intersectCF3D(cone, box.faces[0])
				|| intersectCF3D(cone, box.faces[1])
				|| intersectCF3D(cone, box.faces[2])
				|| intersectCF3D(cone, box.faces[3])
				|| intersectCF3D(cone, box.faces[4])
				|| intersectCF3D(cone, box.faces[5]);
		}

	}

	template <typename T>
	short intersectCN4D(Cone4D cone, NDWindow<T> node)
	{
		//2 indicates inside, 1 indicates intersection
		short ncmp = 1;

		//cone inside cube
		for (int i = 0; i < 4; i++)
		{
			ncmp &= cone.P[i] >= node.minPoint[i] && cone.P[i] <= node.maxPoint[i];
		}
		if (ncmp) return 1;

		//cube inside cone
		ncmp = 1;
		ncmp &= (node.minPoint[0] - cone.P[0])*(node.minPoint[0] - cone.P[0]) + (node.minPoint[1] - cone.P[1])*(node.minPoint[1] - cone.P[1]) + (node.minPoint[2] - cone.P[2])*(node.minPoint[2] - cone.P[2]) <= (cone.a - cone.b*node.minPoint[3])*(cone.a - cone.b*node.minPoint[3]) &&
			(node.minPoint[0] - cone.P[0])*(node.minPoint[0] - cone.P[0]) + (node.maxPoint[1] - cone.P[1])*(node.maxPoint[1] - cone.P[1]) + (node.minPoint[2] - cone.P[2])*(node.minPoint[2] - cone.P[2]) <= (cone.a - cone.b*node.minPoint[3])*(cone.a - cone.b*node.minPoint[3]) &&
			(node.minPoint[0] - cone.P[0])*(node.minPoint[0] - cone.P[0]) + (node.minPoint[1] - cone.P[1])*(node.minPoint[1] - cone.P[1]) + (node.maxPoint[2] - cone.P[2])*(node.maxPoint[2] - cone.P[2]) <= (cone.a - cone.b*node.minPoint[3])*(cone.a - cone.b*node.minPoint[3]) &&
			(node.minPoint[0] - cone.P[0])*(node.minPoint[0] - cone.P[0]) + (node.minPoint[1] - cone.P[1])*(node.minPoint[1] - cone.P[1]) + (node.minPoint[2] - cone.P[2])*(node.minPoint[2] - cone.P[2]) <= (cone.a - cone.b*node.maxPoint[3])*(cone.a - cone.b*node.maxPoint[3]) &&
			(node.minPoint[0] - cone.P[0])*(node.minPoint[0] - cone.P[0]) + (node.maxPoint[1] - cone.P[1])*(node.maxPoint[1] - cone.P[1]) + (node.maxPoint[2] - cone.P[2])*(node.maxPoint[2] - cone.P[2]) <= (cone.a - cone.b*node.minPoint[3])*(cone.a - cone.b*node.minPoint[3]) &&
			(node.minPoint[0] - cone.P[0])*(node.minPoint[0] - cone.P[0]) + (node.maxPoint[1] - cone.P[1])*(node.maxPoint[1] - cone.P[1]) + (node.minPoint[2] - cone.P[2])*(node.minPoint[2] - cone.P[2]) <= (cone.a - cone.b*node.maxPoint[3])*(cone.a - cone.b*node.maxPoint[3]) &&
			(node.minPoint[0] - cone.P[0])*(node.minPoint[0] - cone.P[0]) + (node.minPoint[1] - cone.P[1])*(node.minPoint[1] - cone.P[1]) + (node.maxPoint[2] - cone.P[2])*(node.maxPoint[2] - cone.P[2]) <= (cone.a - cone.b*node.maxPoint[3])*(cone.a - cone.b*node.maxPoint[3]) &&
			(node.minPoint[0] - cone.P[0])*(node.minPoint[0] - cone.P[0]) + (node.maxPoint[1] - cone.P[1])*(node.maxPoint[1] - cone.P[1]) + (node.maxPoint[2] - cone.P[2])*(node.maxPoint[2] - cone.P[2]) <= (cone.a - cone.b*node.maxPoint[3])*(cone.a - cone.b*node.maxPoint[3]) &&
			(node.maxPoint[0] - cone.P[0])*(node.maxPoint[0] - cone.P[0]) + (node.minPoint[1] - cone.P[1])*(node.minPoint[1] - cone.P[1]) + (node.minPoint[2] - cone.P[2])*(node.minPoint[2] - cone.P[2]) <= (cone.a - cone.b*node.minPoint[3])*(cone.a - cone.b*node.minPoint[3]) &&
			(node.maxPoint[0] - cone.P[0])*(node.maxPoint[0] - cone.P[0]) + (node.maxPoint[1] - cone.P[1])*(node.maxPoint[1] - cone.P[1]) + (node.minPoint[2] - cone.P[2])*(node.minPoint[2] - cone.P[2]) <= (cone.a - cone.b*node.minPoint[3])*(cone.a - cone.b*node.minPoint[3]) &&
			(node.maxPoint[0] - cone.P[0])*(node.maxPoint[0] - cone.P[0]) + (node.minPoint[1] - cone.P[1])*(node.minPoint[1] - cone.P[1]) + (node.maxPoint[2] - cone.P[2])*(node.maxPoint[2] - cone.P[2]) <= (cone.a - cone.b*node.minPoint[3])*(cone.a - cone.b*node.minPoint[3]) &&
			(node.maxPoint[0] - cone.P[0])*(node.maxPoint[0] - cone.P[0]) + (node.minPoint[1] - cone.P[1])*(node.minPoint[1] - cone.P[1]) + (node.minPoint[2] - cone.P[2])*(node.minPoint[2] - cone.P[2]) <= (cone.a - cone.b*node.maxPoint[3])*(cone.a - cone.b*node.maxPoint[3]) &&
			(node.maxPoint[0] - cone.P[0])*(node.maxPoint[0] - cone.P[0]) + (node.maxPoint[1] - cone.P[1])*(node.maxPoint[1] - cone.P[1]) + (node.maxPoint[2] - cone.P[2])*(node.maxPoint[2] - cone.P[2]) <= (cone.a - cone.b*node.minPoint[3])*(cone.a - cone.b*node.minPoint[3]) &&
			(node.maxPoint[0] - cone.P[0])*(node.maxPoint[0] - cone.P[0]) + (node.maxPoint[1] - cone.P[1])*(node.maxPoint[1] - cone.P[1]) + (node.minPoint[2] - cone.P[2])*(node.minPoint[2] - cone.P[2]) <= (cone.a - cone.b*node.maxPoint[3])*(cone.a - cone.b*node.maxPoint[3]) &&
			(node.maxPoint[0] - cone.P[0])*(node.maxPoint[0] - cone.P[0]) + (node.minPoint[1] - cone.P[1])*(node.minPoint[1] - cone.P[1]) + (node.maxPoint[2] - cone.P[2])*(node.maxPoint[2] - cone.P[2]) <= (cone.a - cone.b*node.maxPoint[3])*(cone.a - cone.b*node.maxPoint[3]) &&
			(node.maxPoint[0] - cone.P[0])*(node.maxPoint[0] - cone.P[0]) + (node.maxPoint[1] - cone.P[1])*(node.maxPoint[1] - cone.P[1]) + (node.maxPoint[2] - cone.P[2])*(node.maxPoint[2] - cone.P[2]) <= (cone.a - cone.b*node.maxPoint[3])*(cone.a - cone.b*node.maxPoint[3]);
		if (ncmp) return 2;


		//one vertex of the cube inside cone, unnecessary, but may prompte performance
		ncmp |= (node.minPoint[0] - cone.P[0])*(node.minPoint[0] - cone.P[0]) + (node.minPoint[1] - cone.P[1])*(node.minPoint[1] - cone.P[1]) + (node.minPoint[2] - cone.P[2])*(node.minPoint[2] - cone.P[2]) <= (cone.a - cone.b*node.minPoint[3])*(cone.a - cone.b*node.minPoint[3]) ||
			(node.minPoint[0] - cone.P[0])*(node.minPoint[0] - cone.P[0]) + (node.maxPoint[1] - cone.P[1])*(node.maxPoint[1] - cone.P[1]) + (node.minPoint[2] - cone.P[2])*(node.minPoint[2] - cone.P[2]) <= (cone.a - cone.b*node.minPoint[3])*(cone.a - cone.b*node.minPoint[3]) ||
			(node.minPoint[0] - cone.P[0])*(node.minPoint[0] - cone.P[0]) + (node.minPoint[1] - cone.P[1])*(node.minPoint[1] - cone.P[1]) + (node.maxPoint[2] - cone.P[2])*(node.maxPoint[2] - cone.P[2]) <= (cone.a - cone.b*node.minPoint[3])*(cone.a - cone.b*node.minPoint[3]) ||
			(node.minPoint[0] - cone.P[0])*(node.minPoint[0] - cone.P[0]) + (node.minPoint[1] - cone.P[1])*(node.minPoint[1] - cone.P[1]) + (node.minPoint[2] - cone.P[2])*(node.minPoint[2] - cone.P[2]) <= (cone.a - cone.b*node.maxPoint[3])*(cone.a - cone.b*node.maxPoint[3]) ||
			(node.minPoint[0] - cone.P[0])*(node.minPoint[0] - cone.P[0]) + (node.maxPoint[1] - cone.P[1])*(node.maxPoint[1] - cone.P[1]) + (node.maxPoint[2] - cone.P[2])*(node.maxPoint[2] - cone.P[2]) <= (cone.a - cone.b*node.minPoint[3])*(cone.a - cone.b*node.minPoint[3]) ||
			(node.minPoint[0] - cone.P[0])*(node.minPoint[0] - cone.P[0]) + (node.maxPoint[1] - cone.P[1])*(node.maxPoint[1] - cone.P[1]) + (node.minPoint[2] - cone.P[2])*(node.minPoint[2] - cone.P[2]) <= (cone.a - cone.b*node.maxPoint[3])*(cone.a - cone.b*node.maxPoint[3]) ||
			(node.minPoint[0] - cone.P[0])*(node.minPoint[0] - cone.P[0]) + (node.minPoint[1] - cone.P[1])*(node.minPoint[1] - cone.P[1]) + (node.maxPoint[2] - cone.P[2])*(node.maxPoint[2] - cone.P[2]) <= (cone.a - cone.b*node.maxPoint[3])*(cone.a - cone.b*node.maxPoint[3]) ||
			(node.minPoint[0] - cone.P[0])*(node.minPoint[0] - cone.P[0]) + (node.maxPoint[1] - cone.P[1])*(node.maxPoint[1] - cone.P[1]) + (node.maxPoint[2] - cone.P[2])*(node.maxPoint[2] - cone.P[2]) <= (cone.a - cone.b*node.maxPoint[3])*(cone.a - cone.b*node.maxPoint[3]) ||
			(node.maxPoint[0] - cone.P[0])*(node.maxPoint[0] - cone.P[0]) + (node.minPoint[1] - cone.P[1])*(node.minPoint[1] - cone.P[1]) + (node.minPoint[2] - cone.P[2])*(node.minPoint[2] - cone.P[2]) <= (cone.a - cone.b*node.minPoint[3])*(cone.a - cone.b*node.minPoint[3]) ||
			(node.maxPoint[0] - cone.P[0])*(node.maxPoint[0] - cone.P[0]) + (node.maxPoint[1] - cone.P[1])*(node.maxPoint[1] - cone.P[1]) + (node.minPoint[2] - cone.P[2])*(node.minPoint[2] - cone.P[2]) <= (cone.a - cone.b*node.minPoint[3])*(cone.a - cone.b*node.minPoint[3]) ||
			(node.maxPoint[0] - cone.P[0])*(node.maxPoint[0] - cone.P[0]) + (node.minPoint[1] - cone.P[1])*(node.minPoint[1] - cone.P[1]) + (node.maxPoint[2] - cone.P[2])*(node.maxPoint[2] - cone.P[2]) <= (cone.a - cone.b*node.minPoint[3])*(cone.a - cone.b*node.minPoint[3]) ||
			(node.maxPoint[0] - cone.P[0])*(node.maxPoint[0] - cone.P[0]) + (node.minPoint[1] - cone.P[1])*(node.minPoint[1] - cone.P[1]) + (node.minPoint[2] - cone.P[2])*(node.minPoint[2] - cone.P[2]) <= (cone.a - cone.b*node.maxPoint[3])*(cone.a - cone.b*node.maxPoint[3]) ||
			(node.maxPoint[0] - cone.P[0])*(node.maxPoint[0] - cone.P[0]) + (node.maxPoint[1] - cone.P[1])*(node.maxPoint[1] - cone.P[1]) + (node.maxPoint[2] - cone.P[2])*(node.maxPoint[2] - cone.P[2]) <= (cone.a - cone.b*node.minPoint[3])*(cone.a - cone.b*node.minPoint[3]) ||
			(node.maxPoint[0] - cone.P[0])*(node.maxPoint[0] - cone.P[0]) + (node.maxPoint[1] - cone.P[1])*(node.maxPoint[1] - cone.P[1]) + (node.minPoint[2] - cone.P[2])*(node.minPoint[2] - cone.P[2]) <= (cone.a - cone.b*node.maxPoint[3])*(cone.a - cone.b*node.maxPoint[3]) ||
			(node.maxPoint[0] - cone.P[0])*(node.maxPoint[0] - cone.P[0]) + (node.minPoint[1] - cone.P[1])*(node.minPoint[1] - cone.P[1]) + (node.maxPoint[2] - cone.P[2])*(node.maxPoint[2] - cone.P[2]) <= (cone.a - cone.b*node.maxPoint[3])*(cone.a - cone.b*node.maxPoint[3]) ||
			(node.maxPoint[0] - cone.P[0])*(node.maxPoint[0] - cone.P[0]) + (node.maxPoint[1] - cone.P[1])*(node.maxPoint[1] - cone.P[1]) + (node.maxPoint[2] - cone.P[2])*(node.maxPoint[2] - cone.P[2]) <= (cone.a - cone.b*node.maxPoint[3])*(cone.a - cone.b*node.maxPoint[3]);
		if (ncmp) return 1;

		//normal intersection, 3D cone or 3D sphere intersects 8 3D cubes
		//2 sphere boundaries
		Cube3D bcube = { {node.minPoint[0],node.minPoint[1],node.minPoint[2]},{node.maxPoint[0],node.maxPoint[1],node.maxPoint[2]} };
		bool b1, b2;
		if (node.minPoint[3] < 0 || node.minPoint[3] > cone.maxl) b1 = 0;
		else
		{
			Sphere3D ball1 = { {cone.P[0],cone.P[1],cone.P[2]}, cone.a - cone.b*node.minPoint[3] };
			b1 = intersectSV3D(ball1, bcube);
		}
		if (b1) return 1;
		if (node.maxPoint[3] < 0 || node.maxPoint[3] > cone.maxl) b2 = 0;
		else
		{
			Sphere3D ball2 = { {cone.P[0],cone.P[1],cone.P[2]}, cone.a - cone.b*node.maxPoint[3] };
			b2 = intersectSV3D(ball2, bcube);
		}
		if (b2) return 1;


		//6 cone boundaries
		bool b3, b4;
		Cube3D cubesub1 = { {node.minPoint[1], node.minPoint[2], node.minPoint[3]},{node.maxPoint[1], node.maxPoint[2], node.maxPoint[3]} };
		if ((node.minPoint[0] - cone.P[0])*(node.minPoint[0] - cone.P[0]) <= (cone.a - cone.b*cubesub1.minP[2])*(cone.a - cone.b*cubesub1.minP[2]))
		{
			Cone3D conesub1 = { {cone.P[1], cone.P[2]}, (node.minPoint[0] - cone.P[0])*(node.minPoint[0] - cone.P[0]), cone.a, cone.b, cone.maxl };
			b3 = intersectCV3D(conesub1, cubesub1);
		}
		else b3 = 0;
		if (b3) return 1;

		if ((node.maxPoint[0] - cone.P[0])*(node.maxPoint[0] - cone.P[0]) <= (cone.a - cone.b*cubesub1.minP[2])*(cone.a - cone.b*cubesub1.minP[2]))
		{
			Cone3D conesub2 = { {cone.P[1], cone.P[2]}, (node.maxPoint[0] - cone.P[0])*(node.maxPoint[0] - cone.P[0]), cone.a, cone.b, cone.maxl };
			b4 = intersectCV3D(conesub2, cubesub1);
		}
		else b4 = 0;
		if (b4) return 1;

		bool b5, b6;
		Cube3D cubesub2 = { {node.minPoint[0], node.minPoint[2], node.minPoint[3]},{node.maxPoint[0], node.maxPoint[2], node.maxPoint[3]} };
		if ((node.minPoint[1] - cone.P[1])*(node.minPoint[1] - cone.P[1]) <= (cone.a - cone.b*cubesub2.minP[2])*(cone.a - cone.b*cubesub2.minP[2]))
		{
			Cone3D conesub3 = { {cone.P[0], cone.P[2]}, (node.minPoint[1] - cone.P[1])*(node.minPoint[1] - cone.P[1]), cone.a, cone.b, cone.maxl };
			b5 = intersectCV3D(conesub3, cubesub2);
		}
		else b5 = 0;
		if (b5) return 1;

		if ((node.maxPoint[1] - cone.P[1])*(node.maxPoint[1] - cone.P[1]) <= (cone.a - cone.b*cubesub2.minP[2])*(cone.a - cone.b*cubesub2.minP[2]))
		{
			Cone3D conesub4 = { {cone.P[0], cone.P[2]}, (node.maxPoint[1] - cone.P[1])*(node.maxPoint[1] - cone.P[1]), cone.a, cone.b, cone.maxl };
			b6 = intersectCV3D(conesub4, cubesub2);
		}
		else b6 = 0;
		if (b6) return 1;

		bool b7, b8;
		Cube3D cubesub3 = { {node.minPoint[0], node.minPoint[1], node.minPoint[3]},{node.maxPoint[0], node.maxPoint[1], node.maxPoint[3]} };
		if ((node.minPoint[2] - cone.P[2])*(node.minPoint[2] - cone.P[2]) <= (cone.a - cone.b*cubesub3.minP[2])*(cone.a - cone.b*cubesub3.minP[2]))
		{
			Cone3D conesub5 = { {cone.P[0], cone.P[1]}, (node.minPoint[2] - cone.P[2])*(node.minPoint[2] - cone.P[2]), cone.a, cone.b, cone.maxl };
			b7 = intersectCV3D(conesub5, cubesub3);
		}
		else b7 = 0;
		if (b7) return 1;

		if ((node.maxPoint[2] - cone.P[2])*(node.maxPoint[2] - cone.P[2]) <= (cone.a - cone.b*cubesub3.minP[2])*(cone.a - cone.b*cubesub3.minP[2]))
		{
			Cone3D conesub6 = { {cone.P[0], cone.P[1]}, (node.maxPoint[2] - cone.P[2])*(node.maxPoint[2] - cone.P[2]), cone.a, cone.b, cone.maxl };
			b8 = intersectCV3D(conesub6, cubesub3);
		}
		else b8 = 0;
		if (b8) return 1;

		return 0;
	}


	/*Convert 4D (XYZcLoD) view frustum to an nD polytope*/
	NDGeom view2poly(viewPos params, double maxLoD) {
		NDGeom frustum(4);
		auto F = FrustumBuild(params);

		/*5 faces of the 3D XYZ view frustum*/
		frustum.add(halfspace({ F.normal[0][0], F.normal[0][1], F.normal[0][2], 0 }, 0, Sign::ge));
		frustum.add(halfspace({ F.normal[1][0], F.normal[1][1], F.normal[1][2], 0 }, 0, Sign::ge));
		frustum.add(halfspace({ F.normal[2][0], F.normal[2][1], F.normal[2][2], 0 }, 0, Sign::ge));
		frustum.add(halfspace({ F.normal[3][0], F.normal[3][1], F.normal[3][2], 0 }, 0, Sign::ge));

		auto oppoface = halfspace({ F.normal[4][0], F.normal[4][1], F.normal[4][2], 0 }, 0, Sign::eq).Normalize();
		oppoface.b = -params.distance;
		oppoface.s = Sign::ge;
		frustum.add(oppoface);

		/*the last cLoD face approximation*/
		//double disC = sqrt(params.dir[0] * params.dir[0] + params.dir[1] * params.dir[1] + params.dir[2] * params.dir[2]);
		/*another halfspace - approximation*/
		//frustum.add(halfspace({ params.dir[0] / disC, params.dir[1] / disC, params.dir[2] / disC, params.distance / maxLoD }, params.distance, Sign::le));

		double gamma = atan2f(maxLoD, params.distance);

		//rotate the half space clockwise alt radians around the y axis(viewed from y = infinity)
		double alt = atan2(params.dir[2], abs(params.dir[0]));

		//rotate the half space clockwise az radians around the z axis (viewed from z = infinity)
		double az = atan2(params.dir[1], params.dir[0]);
		
		//Up to 5 half spaces covering the cLoI getting more restricted as the distance in the XY plane from the origin increases
		for (int i = -2; i <= 2; i++)
		{
			double ang = M_PI * i / 6;

			double w[4];
			w[0] = sin(gamma) * cos(ang);
			w[1] = sin(gamma) * sin(ang);
			w[2] = 0.0;
			w[3] = cos(gamma);
			auto h = halfspace({ w[0], w[1], w[2], w[3] }, maxLoD*cos(gamma), Sign::le);
			h.w[0] = w[0] * cos(alt) - w[2] * sin(alt);
			h.w[2] = w[0] * sin(alt) + w[2] * cos(alt);
			w[0] = h.w[0];
			w[2] = h.w[2];
			h.w[0] = w[0] * cos(az) - w[1] * sin(az);
			h.w[1] = w[0] * sin(az) + w[1] * cos(az);
			frustum.add(h);
		}
		
		//XZ plane
		for (int i = -2; i <= 2; i++)
		{
			double ang = M_PI * i / 6;

			double w[4];
			w[0] = sin(gamma) * cos(ang);
			w[1] = 0.0;
			w[2] = sin(gamma) * sin(ang);
			w[3] = cos(gamma);
			auto h = halfspace({ w[0], w[1], w[2], w[3] }, maxLoD*cos(gamma), Sign::le);
			h.w[0] = w[0] * cos(alt) - w[2] * sin(alt);
			h.w[2] = w[0] * sin(alt) + w[2] * cos(alt);
			w[0] = h.w[0];
			w[2] = h.w[2];
			h.w[0] = w[0] * cos(az) - w[1] * sin(az);
			h.w[1] = w[0] * sin(az) + w[1] * cos(az);
			frustum.add(h);
		}
		
		CoordTrans trans({ -params.P[0], -params.P[1], -params.P[2], 0 }, { 1,1,1,1 });
		return frustum.Transform(trans);

	}

	NDGeom view2geom(viewPos params, double maxLoD) {
		NDGeom frustum(4);
		auto F = FrustumBuild(params);

		/*5 faces of the 3D XYZ view frustum*/
		frustum.add(halfspace({ F.normal[0][0], F.normal[0][1], F.normal[0][2], 0 }, 0, Sign::ge));
		frustum.add(halfspace({ F.normal[1][0], F.normal[1][1], F.normal[1][2], 0 }, 0, Sign::ge));
		frustum.add(halfspace({ F.normal[2][0], F.normal[2][1], F.normal[2][2], 0 }, 0, Sign::ge));
		frustum.add(halfspace({ F.normal[3][0], F.normal[3][1], F.normal[3][2], 0 }, 0, Sign::ge));

		auto oppoface = halfspace({ F.normal[4][0], F.normal[4][1], F.normal[4][2], 0 }, 0, Sign::eq).Normalize();
		oppoface.b = -params.distance;
		oppoface.s = Sign::ge;
		frustum.add(oppoface);

		/*the real quadratic equation*/
		auto dd = params.distance *params.distance;
		frustum.add(quadspace({ 1,0,1,0,1,0,-dd / (maxLoD*maxLoD),2 * dd / maxLoD }, dd, Sign::le));

		CoordTrans trans({ -params.P[0], -params.P[1], -params.P[2], 0 }, { 1,1,1,1 });
		return frustum.Transform(trans);

	}

	/*****Point operators****/
	template <typename T>
	bool intersectFPP3D(Frustum3D view, NDPoint<T> P)
	{
		short ncmp = 1;
		double D;
		for (int i = 0; i < 5; i++)
		{
			if (i < 4) D = -view.P[0] * view.normal[i][0] - view.P[1] * view.normal[i][1] - view.P[2] * view.normal[i][2];
			else D = -view.CP[1][0] * view.normal[i][0] - view.CP[1][1] * view.normal[i][1] - view.CP[1][2] * view.normal[i][2];
			ncmp &= P[0] * view.normal[i][0] + P[1] * view.normal[i][1] + P[2] * view.normal[i][2] + D >= 0;
		}

		if (ncmp) return true;
		else return false;
	}

	template <typename T>
	bool intersectCPP4D(Cone4D cone, NDPoint<T> P)
	{
		return (P[0] - cone.P[0])*(P[0] - cone.P[0]) + (P[1] - cone.P[1])*(P[1] - cone.P[1]) + (P[2] - cone.P[2])*(P[2] - cone.P[2]) <= (cone.a - cone.b*P[3])*(cone.a - cone.b*P[3]);
	}

}