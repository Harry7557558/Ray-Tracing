#pragma once

#include "Matrix.h"
#include "objdef.h"
#include "D:\Explore\Math\Graph\GraphFun\GraphFun\BitMap.h"
#include <vector>
#include <initializer_list>

#define ERR_EPSILON 1e-6	// minimum distance of intersection
#define ERR_UPSILON 1e+12
#define ERR_ZETA 1e-4	// minumum distance of SDF test, may cause threads jointing problems when too large


/* object parent class */
#define Object_Sign 0xFFFFFFFF
class object {
public:

	object() { /*cout << "\aobject::Constructor is called. This function should never be called. \n";*/ }
	object(object* a) {
		//cout << "\aobject::CpyConstructor is called. This function should never be called. \n";
	}
	~object() {}

	virtual void setcolor(const WebSafeColour &c) {}
	virtual void setcolor(const pixel &c) {}
	virtual void setcolor(const rgblight &c) {}
	virtual void setcolor(const unsigned &c) {}

	/* Intersection Test */
	virtual void meet(intersect &R, const ray &a) { cout << "\aobject::meet is called. This function should never be called. \n"; return; }

	/* Transformations, no longer useful */
	virtual inline void operator *= (matrix<double> M) { cout << "\aobject::operator*=(Matrix) is called. This function should never be called. \n"; }
	virtual inline void operator *= (const double &t) { cout << "\aobject::operator*=(double) is called. This function should never be called. \n"; }
	virtual inline void operator += (const point &a) { cout << "\aobject::operator+= is called. This function should never be called. \n"; }

	/* Rotation, first x, then origin y, then origin z */
	virtual inline void rotate(const double &rx, const double &ry, const double &rz) { cout << "\aobject::rotate is called. This function should never be called. \n"; }

	virtual inline point Max() { cout << "\aobject::Max is called. This function should never be called. \n"; return point(); }
	virtual inline point Min() { cout << "\aobject::Min is called. This function should never be called. \n"; return point(); }

	/* return the type of an object, undefined is -1 */
	virtual int telltype() const { cout << "\aobject::telltype is called. This function should never be called. \n"; return Object_Sign; }
	friend ostream& operator << (ostream& os, const object &a);
};


#ifndef _INC_object2D
#define _INC_object2D

/* Smooth opacity surface, without volumn */
class objectSF : public object {
public:
	rgblight reflect;
	objectSF() {}
	objectSF(objectSF* a) {
		reflect = a->reflect;
	}
	~objectSF() {}

	void setcolor(const WebSafeColour &c) { reflect = color(c); }
	void setcolor(const pixel &c) { reflect = c; }
	void setcolor(const rgblight &c) { reflect = c; }
	void setcolor(const unsigned &hex) { reflect = pixel(hex); }

};

/* infinite large plain */
#define Plane_Sign 0x00000000
class plane :public objectSF {
public:
	point N; double D;	// Ax+By+Cz=D

	// default, z=0
	plane() :D(0) { N.z = 1; }
	// normal N, through origin
	plane(const point &N) {
		this->N = N, D = 0;
	}
	// horizontal plain, z=D
	plane(const double &D) {
		N.z = 1, this->D = D;
	}
	// point P and normal N
	plane(const point &P, const point &N) {
		this->N = N; D = dot(P, N);
	}
	// through three points
	plane(const point &A, const point &B, const point &C) {
		N = cross(B - A, C - A); D = dot(A, N);
	}
	// one point and two direction vectors
	/*plane(const point &P, const point &D1, const point &D2) {
		N = cross(D1, D2); D = dot(P, N);
	}*/
	// x, y, and z intercepts
	plane(const double &x_int, const double &y_int, const double &z_int) {
		N.x = y_int * z_int, N.y = x_int * z_int, N.z = x_int * y_int;
		D = x_int * y_int * z_int;
	}
	plane(const plane &p) :N(p.N), D(p.D) {
		reflect = p.reflect;
	}
	plane(const plane *p) :N(p->N), D(p->D) {
		reflect = p->reflect;
	}
	~plane() {}

	inline point Max() {
		if (N.x == 0 && N.y == 0 && N.z == 0) return point(0, 0, 0);
		if (N.x == 0 && N.y == 0) return point(INFINITY, INFINITY, D / N.z + 0.01);
		if (N.x == 0 && N.z == 0) return point(INFINITY, D / N.y + 0.01, INFINITY);
		if (N.y == 0 && N.z == 0) return point(D / N.x + 0.01, INFINITY, INFINITY);
		return point(INFINITY, INFINITY, INFINITY);
	}
	inline point Min() {
		if (N.x == 0 && N.y == 0 && N.z == 0) return point(0, 0, 0);
		if (N.x == 0 && N.y == 0) return point(-INFINITY, -INFINITY, D / N.z - 0.01);
		if (N.x == 0 && N.z == 0) return point(-INFINITY, D / N.y - 0.01, -INFINITY);
		if (N.y == 0 && N.z == 0) return point(D / N.x - 0.01, -INFINITY, -INFINITY);
		return point(-INFINITY, -INFINITY, -INFINITY);
	}
	/*inline void operator *= (matrix<double> M) {
		point A, B, C;
		if (N.z != 0) {
			A = point(0, 0, D / N.z); B = point(1, 0, (D - N.x) / N.z); C = point(0, 1, (D - N.y) / N.z);
		}
		else if (N.y != 0) {
			A = point(0, D / N.y, 0); B = point(1, (D - N.x) / N.y, 0); C = point(0, (D - N.z) / N.y, 1);
		}
		else {
			A = point(D / N.x, 0, 0); B = point((D - N.y) / N.x, 1, 0); C = point((D - N.z) / N.x, 0, 1);
		}
		A *= M, B *= M, C *= M;
		N = cross(B - A, C - A); D = dot(N, A);
	}*/
	inline void operator *= (const double &t) {
		D *= t;
	}
	inline void operator += (const point &a) {
		D += dot(a, N);
	}

	inline void rotate(const double &rx, const double &ry, const double &rz) {
		matrix<double> M = matrix<double>({ {cos(rz),-sin(rz),0}, {sin(rz),cos(rz),0}, {0,0,1} })
			* matrix<double>({ {cos(ry),0,sin(ry)}, {0,1,0}, {-sin(ry),0,cos(ry)} })
			* matrix<double>({ {1,0,0}, {0,cos(rx),-sin(rx)}, {0,sin(rx),cos(rx)} });
		point A, B, C;
		if (N.z != 0) {
			A = point(0, 0, D / N.z); B = point(1, 0, (D - N.x) / N.z); C = point(0, 1, (D - N.y) / N.z);
		}
		else if (N.y != 0) {
			A = point(0, D / N.y, 0); B = point(1, (D - N.x) / N.y, 0); C = point(0, (D - N.z) / N.y, 1);
		}
		else {
			A = point(D / N.x, 0, 0); B = point((D - N.y) / N.x, 1, 0); C = point((D - N.z) / N.x, 0, 1);
		}
		A *= M, B *= M, C *= M;
		N = cross(B - A, C - A); D = dot(N, A);
	}

	void meet(intersect &R, const ray &a) {
		R.meet = 0;
		double t = (D - dot(N, a.orig)) / dot(N, a.dir);
		if (t < ERR_EPSILON || t > ERR_UPSILON) return;
		R.intrs = t * a.dir + a.orig;
		R.dist = t * a.dir.mod();
		R.reflect = 2 * (-dot(a.dir, N) / dot(N, N)) * N + a.dir;
		R.meet = 1;
		return;
	}

	friend ostream& operator << (ostream& os, const plane &a) {
		if (abs(a.N.x) > ERR_EPSILON) os << a.N.x << "*x";
		if (abs(a.N.y) > ERR_EPSILON) os << showpos << a.N.y << "*y";
		if (abs(a.N.z) > ERR_EPSILON) os << showpos << a.N.z << "*z";
		os << "=" << noshowpos << a.D;
		return os;
	}
	int telltype() const { return Plane_Sign; }
};


/* triangle class */
#define Triangle_Sign 0x00000001
class triangle :public objectSF {
public:
	point A, B, C;
	triangle() {}
	triangle(triangle* a) {
		reflect = a->reflect;
		A = a->A, B = a->B, C = a->C;
	}
	triangle(const point &A, const point &B, const point &C) {
		this->A = A, this->B = B, this->C = C;
	}
	triangle(const initializer_list<double> &A, const initializer_list<double> &B, const initializer_list<double> &C) {
		this->A = A, this->B = B, this->C = C;
	}
	void operator = (const initializer_list<initializer_list<double>> &V) {
		A = *V.begin(), B = *(V.begin() + 1), C = *(V.begin() + 2);
	}
	~triangle() {}
	void meet(intersect &R, const ray &a) {
		// Algorithm: http://www.graphics.cornell.edu/pubs/1997/MT97.pdf
		R.meet = 0;
		point E1 = B - A, E2 = C - A, T, P = cross(a.dir, E2), Q;
		double det = dot(E1, P);
		if (abs(det) < ERR_EPSILON) return;
		T = a.orig - A;
		double t, u, v;
		// t: the multiple that is extended when the vector intersects the triangle; 
		// u,v: (1-u-v)·A + u·B + v·C  is where the vector intersects the triangle; 
		// " O + t·D = (1-u-v)·A + u·B + v·C "
		// This three variables may be useful. 
		u = dot(T, P) / det;
		if (u < 0.0 || u > 1.0) return;
		Q = cross(T, E1);
		v = dot(a.dir, Q) / det;
		if (v < 0.0 || u + v > 1.0) return;
		t = dot(E2, Q) / det;
		if (t < ERR_EPSILON) return;
		R.meet = 1;
		R.dist = t * a.dir.mod();
		R.intrs = (1 - u - v)*A + u * B + v * C;
		point OB = B - R.intrs, OC = C - R.intrs;
		point ON = cross(OB, OC);
		ON *= dot(a.dir, ON) / dot(ON, ON);
		R.reflect = a.dir - 2 * ON;
		return;
	}
	inline point Max() { return point(max({ A.x, B.x, C.x }), max({ A.y, B.y, C.y }), max({ A.z, B.z, C.z })); }
	inline point Min() { return point(min({ A.x, B.x, C.x }), min({ A.y, B.y, C.y }), min({ A.z, B.z, C.z })); }
	inline friend triangle operator * (matrix<double> M, triangle a) {
		a.A = M * a.A, a.B = M * a.B, a.C = M * a.C;
		return a;
	}
	/*inline void operator *= (matrix<double> M) {
		A = M * A, B = M * B, C = M * C;
	}*/
	inline void rotate(const double &rx, const double &ry, const double &rz) {
		matrix<double> M = matrix<double>({ {cos(rz),-sin(rz),0}, {sin(rz),cos(rz),0}, {0,0,1} })
			* matrix<double>({ {cos(ry),0,sin(ry)}, {0,1,0}, {-sin(ry),0,cos(ry)} })
			* matrix<double>({ {1,0,0}, {0,cos(rx),-sin(rx)}, {0,sin(rx),cos(rx)} });
		A = M * A, B = M * B, C = M * C;
	}
	inline void operator *= (const double &t) {
		A *= t, B *= t, C *= t;
	}
	inline friend triangle operator + (triangle t, const point &a) {
		t.A += a, t.B += a, t.C += a;
		return t;
	}
	inline void operator += (const point &a) {
		A += a, B += a, C += a;
	}
	friend ostream& operator << (ostream& os, const triangle &a) {
		/*os << "(" << a.A.x << "," << a.A.y << "," << a.A.z << "), ("
			<< a.B.x << "," << a.B.y << "," << a.B.z << "), (" << a.C.x << "," << a.C.y << "," << a.C.z << ")";*/
		os << "Polyline(" << a.A << "," << a.B << "," << a.C << "," << a.A << ")";
		return os;
	}
	int telltype() const { return Triangle_Sign; }
};

/* parallelogram class */
#define Parallelogram_Sign 0x00000002
class parallelogram :public objectSF {
public:
	point O, A, B;
	parallelogram() {}
	parallelogram(parallelogram* a) {
		reflect = a->reflect;
		A = a->A, B = a->B, O = a->O;
	}
	parallelogram(const point &O, const point &A, const point &B) {
		this->O = O, this->A = O + A, this->B = O + B;
	}
	parallelogram(const point &O, const point &A, const point &B, bool absolute) {
		this->O = O, this->A = A, this->B = B;
		if (!absolute) this->A += O, this->B += O;
	}
	parallelogram(const initializer_list<double> &O, const initializer_list<double> &A, const initializer_list<double> &B) {
		this->O = O, this->A = A, this->B = B;
		this->A += O, this->B += O;
	}
	void operator = (const initializer_list<initializer_list<double>> &V) {
		O = *V.begin(), A = *(V.begin() + 1), B = *(V.begin() + 2);
	}
	~parallelogram() {}
	void meet(intersect &R, const ray &a) {
		R.meet = 0;
		point E1 = A - O, E2 = B - O, T, P = cross(a.dir, E2), Q;
		double det = dot(E1, P);
		if (abs(det) < ERR_EPSILON) return;
		T = a.orig - O;
		double t, u, v;
		u = dot(T, P) / det;
		if (u < 0.0 || u > 1.0) return;
		Q = cross(T, E1);
		v = dot(a.dir, Q) / det;
		if (v < 0.0 || v > 1.0) return;	// For triangles: v < 0.0 || u + v > 1.0
		t = dot(E2, Q) / det;
		if (t < ERR_EPSILON) return;
		R.meet = 1;
		R.dist = t * a.dir.mod();
		R.intrs = (1 - u - v)*O + u * A + v * B;
		point OA = A - R.intrs, OB = B - R.intrs;	// Ps. Sometimes OA∥OB and 0/0=NAN occurs.
		point ON = cross(OA, OB);
		ON *= dot(a.dir, ON) / dot(ON, ON);
		R.reflect = a.dir - 2 * ON;
		return;
	}
	inline point Max() {
		return point(max({ O.x, A.x, B.x, A.x + B.x - O.x }),
			max({ O.y, A.y, B.y, A.y + B.y - O.y }), max({ O.y, A.y, B.y,A.y + B.y - O.y }));
	}
	inline point Min() {
		return point(min({ O.x, A.x, B.x, A.x + B.x - O.x }),
			min({ O.y, A.y, B.y, A.y + B.y - O.y }), min({ O.y, A.y, B.y,A.y + B.y - O.y }));
	}
	inline friend parallelogram operator * (matrix<double> M, parallelogram a) {
		a.O = M * a.O, a.A = M * a.A, a.B = M * a.B;
		return a;
	}
	inline void operator *= (matrix<double> M) {
		A = M * A, B = M * B, O = M * O;
	}
	inline void rotate(const double &rx, const double &ry, const double &rz) {
		matrix<double> M = matrix<double>({ {cos(rz),-sin(rz),0}, {sin(rz),cos(rz),0}, {0,0,1} })
			* matrix<double>({ {cos(ry),0,sin(ry)}, {0,1,0}, {-sin(ry),0,cos(ry)} })
			* matrix<double>({ {1,0,0}, {0,cos(rx),-sin(rx)}, {0,sin(rx),cos(rx)} });
		A = M * A, B = M * B, O = M * O;
	}
	inline void operator *= (const double &t) {
		A *= t, B *= t, O *= t;
	}
	inline friend parallelogram operator + (parallelogram t, const point &a) {
		t.A += a, t.B += a, t.O += a;
		return t;
	}
	inline void operator += (const point &a) {
		A += a, B += a, O += a;
	}
	friend ostream& operator << (ostream& os, const parallelogram &a) {
		/*os << "(" << a.O.x << "," << a.O.y << "," << a.O.z << ") -> ("
			<< a.A.x << "," << a.A.y << "," << a.A.z << "), (" << a.B.x << "," << a.B.y << "," << a.B.z << "), end = ("
			<< (a.A.x + a.B.x - a.O.x) << "," << (a.A.y + a.B.y - a.O.y) << "," << (a.A.z + a.B.z - a.O.z) << ")";*/
			//os << "Polyline(" << a.O << "," << a.A << "," << (a.A + a.B - a.O) << "," << a.B << "," << a.O << ")";
		os << "Surface(" << noshowpos << a.O.x << "*(1-u-v)" << showpos << a.A.x << "*u" << showpos << a.B.x << "*v" << ", "
			<< noshowpos << a.O.y << "*(1-u-v)" << showpos << a.A.y << "*u" << showpos << a.B.y << "*v" << ", "
			<< noshowpos << a.O.z << "*(1-u-v)" << showpos << a.A.z << "*u" << showpos << a.B.z << "*v" << ", "
			<< "u, 0, 1, v, 0, 1)" << noshowpos;
		return os;
	}
	int telltype() const { return Parallelogram_Sign; }
};


/* spacial circle class */
#define Circle_Sign 0x00000003
class circle :public objectSF {
public:
	point C; double r; double rx, ry, rz;	// First rotate x, then rotate z    R = Rz * Rx		// Clockwise??!!
	circle() :r(0), rx(0), ry(0), rz(0) {}
	circle(circle* a) {
		reflect = a->reflect;
		C = a->C, r = a->r; rx = a->rx, ry = a->ry, rz = a->rz;
	}
	circle(const point &C, const double &r) {
		this->C = C, this->r = r, rx = ry = rz = 0;
	}
	circle(const point &C, const double &r, const double &rx, const double &rz) {
		this->C = C, this->r = r, this->rx = rx, this->ry = 0, this->rz = rz;
	}
	circle(const point &C, const double &r, const double &rx, const double &ry, const double &rz) {
		this->C = C, this->r = r, this->rx = rx, this->ry = ry, this->rz = rz;
	}
	inline friend circle operator + (circle c, const point &p) {
		c.C += p; return c;
	}
	inline void operator += (const point &a) {
		C += a;
	}
	/*inline friend circle operator * (const matrix<double> &Orthogonal, circle a) {
		a *= Orthogonal; return a;
	}*/
	/*inline void operator *= (matrix<double> Orthogonal) {
		// Matrix must be the product of a rotation matrix with determinant 1 and a positive constant
		C *= Orthogonal;
		double d = cbrt(det(Orthogonal)); r *= d;
		//cout << (Orthogonal / d) << endl << ((Orthogonal / d) * Vector<double>(1, 1, 1)) << endl << endl;
		//rx += atan2(Orthogonal[2][1], Orthogonal[2][2]); rz += atan2(Orthogonal[1][0], Orthogonal[0][0]);
		rx += acos(Orthogonal[2][2] / d); //rz += acos(Orthogonal[0][0] / d);
		//cout << rx << " " << rz << endl << endl;
		// For Euler's angle: https://stackoverflow.com/questions/15022630/how-to-calculate-the-angle-from-rotation-matrix
	}*/
	inline void rotate(const double &rx, const double &ry, const double &rz) {
		matrix<double> M = matrix<double>({ {cos(rz),-sin(rz),0}, {sin(rz),cos(rz),0}, {0,0,1} })
			* matrix<double>({ {cos(ry),0,sin(ry)}, {0,1,0}, {-sin(ry),0,cos(ry)} })
			* matrix<double>({ {1,0,0}, {0,cos(rx),-sin(rx)}, {0,sin(rx),cos(rx)} });
		C *= M;
		M *= matrix<double>({ {cos(this->rz),-sin(this->rz),0}, {sin(this->rz),cos(this->rz),0}, {0,0,1} })
			* matrix<double>({ {cos(this->ry),0,sin(this->ry)}, {0,1,0}, {-sin(this->ry),0,cos(this->ry)} })
			* matrix<double>({ {1,0,0}, {0,cos(this->rx),-sin(this->rx)}, {0,sin(this->rx),cos(this->rx)} });
		this->rx = atan2(M[2][1], M[2][2]), this->rz = atan2(M[1][0], M[0][0]), this->ry = atan2(-M[2][0], hypot(M[2][1], M[2][2]));
	}
	inline void operator *= (const double &t) {
		C *= t, r *= t;
	}
	inline friend circle operator * (circle B, const double &a) {
		B.C *= a, B.r *= a; return B;
	}
	inline point Max() {
		point R;
		/*R.z = -atan(cos(rx) * tan(rz));
		R.x = r * abs(cos(rz)*cos(R.z) - sin(rz)*cos(rx)*sin(R.z));
		R.z = atan(cos(rx) / tan(rz));
		R.y = r * abs(sin(rz)*cos(R.z) + cos(rx)*cos(rz)*sin(R.z));
		R.z = r * abs(sin(rx));*/
		double M1, M2, Ma;
		M1 = cos(ry)*cos(rz), M2 = -cos(rx)*sin(rz) + sin(rx)*sin(ry)*cos(rz);
		Ma = atan2(M2, M1);
		R.x = r * abs(M1*cos(Ma) + M2 * sin(Ma));
		M1 = cos(ry)*sin(rz), M2 = cos(rx)*cos(rz) + sin(rx)*sin(ry)*sin(rz);
		Ma = atan2(M2, M1);
		R.y = r * abs(M1*cos(Ma) + M2 * sin(Ma));
		M1 = -sin(ry), M2 = sin(rx)*cos(ry);
		Ma = atan2(M2, M1);
		R.z = r * abs(M1*cos(Ma) + M2 * sin(Ma));
		R += C; return R;
	}
	inline point Min() {
		point R;
		/*R.z = -atan(cos(rx) * tan(rz));
		R.x = r * abs(cos(rz)*cos(R.z) - sin(rz)*cos(rx)*sin(R.z));
		R.z = atan(cos(rx) / tan(rz));
		R.y = r * abs(sin(rz)*cos(R.z) + cos(rx)*cos(rz)*sin(R.z));
		R.z = r * abs(sin(rx));*/
		double M1, M2, Ma;
		M1 = cos(ry)*cos(rz), M2 = -cos(rx)*sin(rz) + sin(rx)*sin(ry)*cos(rz);
		Ma = atan2(M2, M1);
		R.x = r * abs(M1*cos(Ma) + M2 * sin(Ma));
		M1 = cos(ry)*sin(rz), M2 = cos(rx)*cos(rz) + sin(rx)*sin(ry)*sin(rz);
		Ma = atan2(M2, M1);
		R.y = r * abs(M1*cos(Ma) + M2 * sin(Ma));
		M1 = -sin(ry), M2 = sin(rx)*cos(ry);
		Ma = atan2(M2, M1);
		R.z = r * abs(M1*cos(Ma) + M2 * sin(Ma));
		return C - R;
	}
	void meet(intersect &R, const ray &a) {
		ray s(a.orig - C, a.dir);
		double sx, cy;
		sx = -sin(rz), cy = cos(rz);
		s.orig = point(cy*s.orig.x - sx * s.orig.y, sx * s.orig.x + cy * s.orig.y, s.orig.z),
			s.dir = point(cy*s.dir.x - sx * s.dir.y, sx * s.dir.x + cy * s.dir.y, s.dir.z);
		sx = -sin(ry), cy = cos(ry);
		s.orig = point(cy*s.orig.x + sx * s.orig.z, s.orig.y, -sx * s.orig.x + cy * s.orig.z),
			s.dir = point(cy*s.dir.x + sx * s.dir.z, s.dir.y, -sx * s.dir.x + cy * s.dir.z);
		sx = -sin(rx), cy = cos(rx);
		s.orig = point(s.orig.x, cy*s.orig.y - sx * s.orig.z, sx * s.orig.y + cy * s.orig.z),
			s.dir = point(s.dir.x, cy*s.dir.y - sx * s.dir.z, sx * s.dir.y + cy * s.dir.z);

		R.meet = 0;
		if ((s.orig.z < -ERR_EPSILON) == (s.dir.z < -ERR_EPSILON)) return;
		if ((s.orig.z > ERR_EPSILON) == (s.dir.z > ERR_EPSILON)) return;
		double t = -s.orig.z / s.dir.z;
		R.intrs.x = s.dir.x*t + s.orig.x, R.intrs.y = s.dir.y*t + s.orig.y;
		R.dist = R.intrs.x * R.intrs.x + R.intrs.y * R.intrs.y;
		if (R.dist > r*r) return;
		R.dist = (R.intrs - s.orig).mod();
		//R.reflect = point(2 * R.intrs.x - ry.orig.x, 2 * R.intrs.y - ry.orig.y, ry.orig.z);
		R.reflect = point(s.dir.x, s.dir.y, -s.dir.z);

		sx = -sx;
		R.reflect = point(R.reflect.x, cy*R.reflect.y - sx * R.reflect.z, sx * R.reflect.y + cy * R.reflect.z),
			R.intrs = point(R.intrs.x, cy*R.intrs.y - sx * R.intrs.z, sx * R.intrs.y + cy * R.intrs.z);
		sx = sin(ry), cy = cos(ry);
		R.reflect = point(cy*R.reflect.x + sx * R.reflect.z, R.reflect.y, -sx * R.reflect.x + cy * R.reflect.z),
			R.intrs = point(cy*R.intrs.x + sx * R.intrs.z, R.intrs.y, -sx * R.intrs.x + cy * R.intrs.z);
		sx = sin(rz), cy = cos(rz);
		R.reflect = point(cy*R.reflect.x - sx * R.reflect.y, sx * R.reflect.x + cy * R.reflect.y, R.reflect.z),
			R.intrs = point(cy*R.intrs.x - sx * R.intrs.y, sx * R.intrs.x + cy * R.intrs.y, R.intrs.z);
		R.intrs += C;
		R.meet = 1;
		return;
	}
	friend ostream& operator << (ostream& os, const circle &a) {
		os << "Rotate(Rotate(Rotate(Surface(";
		os << "v*cos(u)" << showpos << a.C.x << ", ";
		os << "v*sin(u)" << showpos << a.C.y << ", ";
		os << noshowpos << a.C.z << ", ";
		os << "u, 0, 2*pi, v, 0, " << noshowpos << a.r << "), ";
		os << noshowpos << a.rx << ", " << point(a.C) << ", xAxis), "
			<< noshowpos << a.ry << ", " << point(a.C) << ", yAxis), "
			<< noshowpos << a.rz << ", " << point(a.C) << ", zAxis)";
		return os;
	}
	int telltype() const { return Circle_Sign; }
};


/* cylinder surface class, without bottom */
#define Cylinder_Sign 0x00000004
class cylinder :public objectSF {
public:
	point C; double r, h; double rx, ry, rz;
	cylinder() :r(0), h(0), rx(0), ry(0), rz(0) {}
	cylinder(cylinder* a) {
		reflect = a->reflect;
		C = a->C, r = a->r, h = a->h; rx = a->rx, ry = a->ry, rz = a->rz;
	}
	cylinder(const point &C, const double &r, const double &h) {
		this->C = C, this->r = r, this->h = h, rx = ry = rz = 0;
	}
	cylinder(const point &C, const double &r, const double &h, const double &rx, const double &rz) {
		this->C = C, this->r = r, this->h = h, this->rx = rx, this->ry = 0, this->rz = rz;
	}
	cylinder(const point &C, const double &r, const double &h, const double &rx, const double &ry, const double &rz) {
		this->C = C, this->r = r, this->h = h, this->rx = rx, this->ry = ry, this->rz = rz;
	}
	cylinder(const point &L0, const point &L1, const double &r) {
		point b = L1 - L0;
		C = L0, this->r = r, this->h = b.mod();
		b /= h;
		/*point v = cross(point(0, 0, 1), b);
		matrix<double> V = matrix<double>({ {0,-v.z,v.y}, {v.z,0,-v.x}, {-v.y,v.x,0} });
		V += pow(V, 2) / (1 + b.z) + matrix<double>(3, Matrix_type::_M_Identity);
		rx = atan2(V[1][0], V[0][0]);
		rz = atan2(V[2][1], V[2][2]);
		ry = atan2(-V[2][0], hypot(V[2][1], V[2][2]));
		// https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d */
		rx = acos(b.z);
		rz = atan2(b.x, -b.y);
		ry = 0;
	}
	inline void operator += (const point &a) {
		C += a;
	}
	/*inline void operator *= (matrix<double> Orthogonal) {
		C *= Orthogonal;
		double d = cbrt(det(Orthogonal)); r *= d, h *= d;
		cout << (Orthogonal / d) << endl;
		rx += acos(Orthogonal[2][2] / d); //rz += acos(Orthogonal[0][0] / d);
	}*/
	inline void rotate(const double &rx, const double &ry, const double &rz) {
		matrix<double> M = matrix<double>({ {cos(rz),-sin(rz),0}, {sin(rz),cos(rz),0}, {0,0,1} })
			* matrix<double>({ {cos(ry),0,sin(ry)}, {0,1,0}, {-sin(ry),0,cos(ry)} })
			* matrix<double>({ {1,0,0}, {0,cos(rx),-sin(rx)}, {0,sin(rx),cos(rx)} });
		C *= M;
		M *= matrix<double>({ {cos(this->rz),-sin(this->rz),0}, {sin(this->rz),cos(this->rz),0}, {0,0,1} })
			* matrix<double>({ {cos(this->ry),0,sin(this->ry)}, {0,1,0}, {-sin(this->ry),0,cos(this->ry)} })
			* matrix<double>({ {1,0,0}, {0,cos(this->rx),-sin(this->rx)}, {0,sin(this->rx),cos(this->rx)} });
		this->rx = atan2(M[2][1], M[2][2]), this->rz = atan2(M[1][0], M[0][0]), this->ry = atan2(-M[2][0], hypot(M[2][1], M[2][2]));
	}
	inline void operator *= (const double &t) {
		C *= t, r *= t, h *= t;
	}
	inline point Max() {
		// https://www.geogebra.org/m/d6fybucd
		point R;
		double M1, M2, M3, Ma;
		M1 = cos(ry)*cos(rz), M2 = -cos(rx)*sin(rz) + sin(rx)*sin(ry)*cos(rz), M3 = sin(rx)*sin(rz) + cos(rx)*sin(ry)*cos(rz);
		Ma = atan2(M2, M1);
		R.x = M1 * r*cos(Ma) + M2 * r*sin(Ma);
		if (M3 > 0) R.x += M3 * h;
		M1 = cos(ry)*sin(rz), M2 = cos(rx)*cos(rz) + sin(rx)*sin(ry)*sin(rz), M3 = -sin(rx)*cos(rz) + cos(rx)*sin(ry)*sin(rz);
		Ma = atan2(M2, M1);
		R.y = M1 * r*cos(Ma) + M2 * r*sin(Ma);
		if (M3 > 0) R.y += M3 * h;
		M1 = -sin(ry), M2 = sin(rx)*cos(ry), M3 = cos(rx)*cos(ry);
		Ma = atan2(M2, M1);
		R.z = M1 * r*cos(Ma) + M2 * r*sin(Ma);
		if (M3 > 0) R.z += M3 * h;
		return C + R;
	}
	inline point Min() {
		point R;
		double M1, M2, M3, Ma;
		M1 = cos(ry)*cos(rz), M2 = -cos(rx)*sin(rz) + sin(rx)*sin(ry)*cos(rz), M3 = sin(rx)*sin(rz) + cos(rx)*sin(ry)*cos(rz);
		Ma = atan2(M2, M1);
		R.x = -(M1 * r*cos(Ma) + M2 * r*sin(Ma));
		if (M3 < 0) R.x += M3 * h;
		M1 = cos(ry)*sin(rz), M2 = cos(rx)*cos(rz) + sin(rx)*sin(ry)*sin(rz), M3 = -sin(rx)*cos(rz) + cos(rx)*sin(ry)*sin(rz);
		Ma = atan2(M2, M1);
		R.y = -(M1 * r*cos(Ma) + M2 * r*sin(Ma));
		if (M3 < 0) R.y += M3 * h;
		M1 = -sin(ry), M2 = sin(rx)*cos(ry), M3 = cos(rx)*cos(ry);
		Ma = atan2(M2, M1);
		R.z = -(M1 * r*cos(Ma) + M2 * r*sin(Ma));
		if (M3 < 0) R.z += M3 * h;
		return C + R;
	}
	void meet(intersect &R, const ray &a) {
		ray s(a.orig - C, a.dir);
		double sx, cy;
		sx = -sin(rz), cy = cos(rz);
		s.orig = point(cy*s.orig.x - sx * s.orig.y, sx * s.orig.x + cy * s.orig.y, s.orig.z),
			s.dir = point(cy*s.dir.x - sx * s.dir.y, sx * s.dir.x + cy * s.dir.y, s.dir.z);
		sx = -sin(ry), cy = cos(ry);
		s.orig = point(cy*s.orig.x + sx * s.orig.z, s.orig.y, -sx * s.orig.x + cy * s.orig.z),
			s.dir = point(cy*s.dir.x + sx * s.dir.z, s.dir.y, -sx * s.dir.x + cy * s.dir.z);
		sx = -sin(rx), cy = cos(rx);
		s.orig = point(s.orig.x, cy*s.orig.y - sx * s.orig.z, sx * s.orig.y + cy * s.orig.z),
			s.dir = point(s.dir.x, cy*s.dir.y - sx * s.dir.z, sx * s.dir.y + cy * s.dir.z);

		R.meet = 0;
		double a0 = s.dir.x*s.dir.x + s.dir.y*s.dir.y, b0 = s.dir.x * s.orig.x + s.dir.y * s.orig.y;
		double delta = b0 * b0 - a0 * (s.orig.x*s.orig.x + s.orig.y*s.orig.y - r * r);
		if (delta < 0) return;
		delta = sqrt(delta);
		R.dist = (-delta - b0) / a0;
		R.intrs.z = s.dir.z*R.dist + s.orig.z;
		if (R.dist < ERR_EPSILON || (R.intrs.z < 0 || R.intrs.z > h)) {
			R.dist = (delta - b0) / a0;
			R.intrs.z = s.dir.z*R.dist + s.orig.z;
			if (R.dist < ERR_EPSILON || (R.intrs.z < 0 || R.intrs.z > h)) return;
		}
		R.intrs.x = s.dir.x*R.dist + s.orig.x;
		R.intrs.y = s.dir.y*R.dist + s.orig.y;
		R.dist *= s.dir.mod();
		R.reflect.x = R.intrs.x, R.reflect.y = R.intrs.y, R.reflect.z = 0;
		R.reflect /= R.reflect.mod();
		R.reflect *= -2 * dot(s.dir, R.reflect);
		R.reflect += s.dir;

		sx = -sx;
		R.reflect = point(R.reflect.x, cy*R.reflect.y - sx * R.reflect.z, sx * R.reflect.y + cy * R.reflect.z),
			R.intrs = point(R.intrs.x, cy*R.intrs.y - sx * R.intrs.z, sx * R.intrs.y + cy * R.intrs.z);
		sx = sin(ry), cy = cos(ry);
		R.reflect = point(cy*R.reflect.x + sx * R.reflect.z, R.reflect.y, -sx * R.reflect.x + cy * R.reflect.z),
			R.intrs = point(cy*R.intrs.x + sx * R.intrs.z, R.intrs.y, -sx * R.intrs.x + cy * R.intrs.z);
		sx = sin(rz), cy = cos(rz);
		R.reflect = point(cy*R.reflect.x - sx * R.reflect.y, sx * R.reflect.x + cy * R.reflect.y, R.reflect.z),
			R.intrs = point(cy*R.intrs.x - sx * R.intrs.y, sx * R.intrs.x + cy * R.intrs.y, R.intrs.z);
		R.intrs += C;
		R.meet = 1;
		return;
	}
	friend ostream& operator << (ostream& os, const cylinder &a) {
		os << "Rotate(Rotate(Rotate(Surface(";
		os << noshowpos << a.r << "*cos(u)" << showpos << a.C.x << ", ";
		os << noshowpos << a.r << "*sin(u)" << showpos << a.C.y << ", ";
		os << "v" << showpos << a.C.z << ", ";
		os << "u, 0, 2*pi, v, 0, " << noshowpos << a.h << "), ";
		os << noshowpos << a.rx << ", " << point(a.C) << ", xAxis), "
			<< noshowpos << a.ry << ", " << point(a.C) << ", yAxis), "
			<< noshowpos << a.rz << ", " << point(a.C) << ", zAxis)";
		return os;
	}
	int telltype() const { return Cylinder_Sign; }
};


/* sphere surface class */
#define Sphere_Sign 0x00000005
class sphere :public objectSF {
public:
	point C; double r;
	sphere() :r(0) {}
	sphere(const point &C, const double &r) {
		this->C = C, this->r = r;
	}
	sphere(sphere* a) {
		reflect = a->reflect;
		C = a->C, r = a->r;
	}
	sphere(const initializer_list<double> &C, const double &r) {
		this->C = point(C), this->r = r;
	}
	void meet(intersect &R, const ray &a) {
		R.meet = 0;
		point p = C - a.orig;
		if (dot(p, a.dir) < 0) return;
		double sm = a.dir.mod();
		double d = cross(p, a.dir).mod() / sm;
		if (d > r) return;
		d *= d;
		R.dist = sqrt(dot(p, p) - d) - sqrt(r*r - d);
		point s = a.dir * (R.dist / sm), n = s - p;
		R.intrs = s + a.orig;
		R.reflect = s - (2 * dot(s, n) / dot(n, n)) * n;
		R.meet = 1;
		return;
	}
	/*inline friend sphere operator * (matrix<double> Orthogonal, sphere a) {
		a *= Orthogonal; return a;
	}
	inline void operator *= (matrix<double> Orthogonal) {
		// Matrix must be the product of an orthogonal matrix and a constant
		C *= Orthogonal; r *= cbrt(det(Orthogonal));
	}*/
	inline void rotate(const double &rx, const double &ry, const double &rz) {
		matrix<double> M = matrix<double>({ {cos(rz),-sin(rz),0}, {sin(rz),cos(rz),0}, {0,0,1} })
			* matrix<double>({ {cos(ry),0,sin(ry)}, {0,1,0}, {-sin(ry),0,cos(ry)} })
			* matrix<double>({ {1,0,0}, {0,cos(rx),-sin(rx)}, {0,sin(rx),cos(rx)} });
		C *= M;
	}
	inline void operator *= (const double &t) {
		C *= t, r *= t;
	}
	inline friend sphere operator + (sphere B, const point &a) {
		B.C += a; return B;
	}
	inline void operator += (const point &a) {
		C += a;
	}
	inline point Max() {
		return point(C.x + r, C.y + r, C.z + r);
	}
	inline point Min() {
		return point(C.x - r, C.y - r, C.z - r);
	}
	inline friend sphere operator * (sphere B, const double &a) {
		B.C *= a, B.r *= a; return B;
	}

	friend ostream& operator << (ostream& os, const sphere &a) {
		os << "Sphere(" << a.C << "," << a.r << ")";
		return os;
	}
	int telltype() const { return Sphere_Sign; }
};

#endif

#ifndef _INC_object2D_dif
#define _INC_object2D_dif

#include <cstdlib>

#ifndef RANDOM_NORMALDISTRIBUTION
#define RANDOM_NORMALDISTRIBUTION
// not necessary to set random number seeds

// calculate inverse error function
inline double erfinv(double x) {
	//return atanh(x);	// quality gets much worse
	double n = log(1 - x * x);
	double t = 0.5 * n + 2 / (PI*0.147);
	if (signbit(x)) return -sqrt(-t + sqrt(t*t - n / 0.147));
	return sqrt(-t + sqrt(t*t - n / 0.147));
}
// produce a random number with given median and variance
inline double randnor(double median, double variance) {
	return erfinv(2.0 * double(rand()) / double(RAND_MAX + 1) - 1)*sqrt(2)*variance + median;
}
inline double randnor_0(double variance) {
	return erfinv(2.0 * double(rand()) / double(RAND_MAX + 1) - 1) * 1.41421356237309504876 * variance;
}

extern double RAND_LCG_DV = 0.36787944117;
#define RAND_LCG_TMS 13.35717028437795
#define RAND_LCG_ADD 0.841470984807897
inline double erfinv0(const double &x) {
	double n = log(1.0 - x * x);
	double t = 0.5 * n + 4.33074675079987308221;
	if (signbit(x)) return -sqrt(-t + sqrt(t*t - n / 0.147));
	return sqrt(-t + sqrt(t*t - n / 0.147));
}
inline double randnor0(double variance) {
	RAND_LCG_DV = fmod(RAND_LCG_DV * RAND_LCG_TMS + RAND_LCG_ADD, 2.0);
	return 1.41421356237309504876 * erfinv0(RAND_LCG_DV - 1) * variance;
}

#endif

/* Opacity surface with diffuse reflection */
class objectSF_dif : public objectSF {
public:
	double difvar;	// variance of random rays
	objectSF_dif() :difvar(0) {}
	objectSF_dif(const objectSF_dif &a) :difvar(a.difvar) { this->reflect = a.reflect; }
	~objectSF_dif() {}

	void setvar(double d) { difvar = d; }

	// random rotate a vector, for rendering diffuse reflection
	inline void rotate_vec(point &a) {
		double sx, sy, sz, cx, cy, cz;
		//sx = randnor(0, difvar); sy = randnor(0, difvar); sz = randnor(0, difvar);
		sx = randnor_0(difvar); sy = randnor_0(difvar); sz = randnor_0(difvar);
		//sx = randnor0(difvar); sy = randnor0(difvar); sz = randnor0(difvar);
		cx = cos(sx), cy = cos(sy), cz = cos(sz); sx = sin(sx), sy = sin(sy), sz = sin(sz);
		double x = a.x, y = a.y, z = a.z;
		a.x = cy * cz*x + (sx*sy*cz - cx * sz)*y + (sx*sz + cx * sy*cz)*z;
		a.y = cy * sz*x + (cx*cz + sx * sy*sz)*y + (cx*sy*sz - sx * cz)*z;
		a.z = -sy * x + sx * cy*y + cx * cy*z;
	}

};

#define Plane_Dif_Sign 0x00000100
class plane_dif : public objectSF_dif/*, public plane*/ {
public:
	point N; double D;	// Ax+By+Cz=D
	plane_dif() { difvar = 0; D = 0, N.z = 1; }
	// normal N, through origin
	plane_dif(const point &N) {
		this->N = N, D = 0;
		difvar = 0;
	}
	// horizontal plain, z=D
	plane_dif(const double &D) {
		N.z = 1, this->D = D;
		difvar = 0;
	}
	// point P and normal N
	plane_dif(const point &P, const point &N) {
		this->N = N; D = dot(P, N);
		difvar = 0;
	}
	// through three points
	plane_dif(const point &A, const point &B, const point &C) {
		N = cross(B - A, C - A); D = dot(A, N);
		difvar = 0;
	}
	plane_dif(const double &x_int, const double &y_int, const double &z_int) {
		N.x = y_int * z_int, N.y = x_int * z_int, N.z = x_int * y_int;
		D = x_int * y_int * z_int;
		difvar = 0;
	}
	plane_dif(const plane_dif &p) :N(p.N), D(p.D) {
		reflect = p.reflect;
		difvar = 0;
	}
	plane_dif(const plane_dif *p) :N(p->N), D(p->D) {
		reflect = p->reflect;
		difvar = 0;
	}
	~plane_dif() {}

	inline point Max() {
		if (N.x == 0 && N.y == 0 && N.z == 0) return point(0, 0, 0);
		if (N.x == 0 && N.y == 0) return point(INFINITY, INFINITY, D / N.z + 0.01);
		if (N.x == 0 && N.z == 0) return point(INFINITY, D / N.y + 0.01, INFINITY);
		if (N.y == 0 && N.z == 0) return point(D / N.x + 0.01, INFINITY, INFINITY);
		return point(INFINITY, INFINITY, INFINITY);
	}
	inline point Min() {
		if (N.x == 0 && N.y == 0 && N.z == 0) return point(0, 0, 0);
		if (N.x == 0 && N.y == 0) return point(-INFINITY, -INFINITY, D / N.z - 0.01);
		if (N.x == 0 && N.z == 0) return point(-INFINITY, D / N.y - 0.01, -INFINITY);
		if (N.y == 0 && N.z == 0) return point(D / N.x - 0.01, -INFINITY, -INFINITY);
		return point(-INFINITY, -INFINITY, -INFINITY);
	}

	void meet(intersect &R, const ray &a) {
		R.meet = 0;
		double t = (D - dot(N, a.orig)) / dot(N, a.dir);
		if (t < ERR_EPSILON || t > ERR_UPSILON) return;
		R.intrs = t * a.dir + a.orig;
		R.dist = t * a.dir.mod();
		R.reflect = 2 * (-dot(a.dir, N) / dot(N, N)) * N + a.dir;
		R.meet = 1;
		return;
	}

	int telltype() const { return Plane_Dif_Sign; }
};

#endif

#ifndef _INC_object2D_col
#define _INC_object2D_col

// Opacity smooth surface, different part using different colors
class objectSF_col : public objectSF {
public:
	objectSF_col() {}
	~objectSF_col() {}

	// get color with calculated intersection data
	virtual void getcol(const intersect &R, rgblight &c) { cout << "objectSF_col::getcol is called. This function should never be called. \a\n"; }
};


// Plane with grid, usually for debug
#define Plane_Grid_Sign 0x00000200
class plane_grid : public objectSF_col {
protected:
	double z_int;
	double wx, hy;
	rgblight c1, c2;
	/*     ^ y
	  +----+----+
	  | c2 | c1 |
	--+----+----+--> x
	  | c1 | c2 |
	  +----+----+
	*/
public:
	plane_grid() : wx(1), hy(1), c1(LightBlue), c2(Gray), z_int(0) { }
	plane_grid(const double &z_int) : wx(1), hy(1), c1(LightBlue), c2(Gray) { this->z_int = z_int; }
	plane_grid(const double &z_int, const double &side_length) : c1(LightBlue), c2(Gray) { this->z_int = z_int, wx = hy = side_length; }
	plane_grid(const double &z_int, const double &side_length_x, const double &side_length_y) : c1(LightBlue), c2(Gray) { this->z_int = z_int, wx = side_length_x, hy = side_length_y; }
	plane_grid(const double &z_int, const double &side_length, const rgblight &c1, const rgblight &c2) { this->z_int = z_int, wx = hy = side_length, this->c1 = c1, this->c2 = c2; }
	plane_grid(const double &z_int, const double &side_length_x, const double &side_length_y, const rgblight &c1, const rgblight &c2) { this->z_int = z_int, wx = side_length_x, hy = side_length_y, this->c1 = c1, this->c2 = c2; }

	inline point Max() {
		return point(INFINITY, INFINITY, z_int);
	}
	inline point Min() {
		return point(-INFINITY, -INFINITY, z_int);
	}

	void meet(intersect &R, const ray &a) {
		R.meet = 0;
		if (abs(a.orig.z - z_int) < ERR_EPSILON) return;
		if ((a.orig.z > z_int) ^ (a.dir.z < 0)) return;
		R.dist = (z_int - a.orig.z) / a.dir.z;
		R.intrs.x = a.dir.x * R.dist + a.orig.x;
		R.intrs.y = a.dir.y * R.dist + a.orig.y;
		R.intrs.z = z_int;
		R.dist = (R.intrs - a.orig).mod();
		R.reflect.x = a.dir.x, R.reflect.y = a.dir.y, R.reflect.z = -a.dir.z;
		R.meet = 1;
		return;
	}

	void getcol(const intersect &R, rgblight &c) {
		c = (int(floor(R.intrs.x / wx)) ^ int(floor(R.intrs.y / hy))) & 1 ? c2 : c1;
	}

	int telltype() const { return Plane_Grid_Sign; }
};


#define Bitmap_Inc_Sign 0x00000201
namespace insertType {
	enum insert_type {
		normal = 0x00000001, zoom = 0x00000002, proportion = 0x00000004,
		leftbottom = 0x00000000, center = 0x00010000,
	};
}
class bitmap_inc : public objectSF_col {
public:
	point O, A, B;
	bitmap M;

	bitmap_inc(const bitmap &M) {
		this->M = M;
		O = point(0, 0), A = point(M.width(), 0), B = point(M.height(), 0);
	}
	bitmap_inc(const bitmap &M, point O, point X, point Y) {
		// relative coordinate
		this->M = M, this->O = O, this->A = X, this->B = Y;
	}
	bitmap_inc(const bitmap &M, point O, point X, point Y, unsigned t) {
		this->M = M, this->O = O;
		switch (t & 0xFFFF) {
		case insertType::normal: {
			this->A = X, this->B = Y;
			break;
		}
		case insertType::zoom: {
			point V = cross(cross(X, Y), X);
			V *= dot(Y, V) / dot(V, V);
			this->A = X, this->B = V;
			break;
		}
		case insertType::proportion: {
			point V = cross(cross(X, Y), X);
			V *= (M.height()*X.mod()) / (M.width()*V.mod());
			this->A = X, this->B = V;
			break;
		}
		}
		if (t >> 16 == 1)  this->O -= A + B, A *= 2, B *= 2;
	}
	bitmap_inc(const bitmap &M, const point &O, const point &X, const point &Y, const insertType::insert_type &t) {
		bitmap_inc(M, O, X, Y, t);
	}

	inline point Max() {
		return point(max({ O.x, O.x + A.x, O.x + B.x, O.x + A.x + B.x }),
			max({ O.y, O.y + A.y, O.y + B.y, O.y + A.y + B.y }), max({ O.z, O.z + A.z, O.z + B.z, O.z + A.z + B.z }));
	}
	inline point Min() {
		return point(min({ O.x, O.x + A.x, O.x + B.x, O.x + A.x + B.x }),
			min({ O.y, O.y + A.y, O.y + B.y, O.y + A.y + B.y }), min({ O.z, O.z + A.z, O.z + B.z, O.z + A.z + B.z }));
	}


	void meet(intersect &R, const ray &a) {
		R.meet = 0;
		point T, P = cross(a.dir, B), Q;
		double det = dot(A, P);
		if (abs(det) < ERR_EPSILON) return;
		T = a.orig - O;
		R.ut = dot(T, P) / det;
		if (R.ut < 0.0 || R.ut > 1.0) return;
		Q = cross(T, A);
		R.vt = dot(a.dir, Q) / det;
		if (R.vt < 0.0 || R.vt > 1.0) return;
		double t = dot(B, Q) / det;
		if (t < ERR_EPSILON) return;
		R.meet = 1;
		R.dist = t * a.dir.mod();
		R.intrs = (1 - R.ut - R.vt)*O + R.ut * A + R.vt * B;
		point ON = O - R.intrs;
		point OA = ON + A, OB = ON + B;
		ON = cross(OA, OB);
		ON *= dot(a.dir, ON) / dot(ON, ON);
		R.reflect = a.dir - 2 * ON;
		return;
	}

	void getcol(const intersect &R, rgblight &c) {
		unsigned x = R.ut * M.width(), y = R.vt * M.height();
		c = M[y][x];
	}

	int telltype() const { return Bitmap_Inc_Sign; }
};

#endif



#ifndef _INC_object3D
#define _INC_object3D

/* 3D objects, with volumn */
/* DEBUG */
class object3D : public object {
public:
	rgblight attcoe;	// light attenuation coefficient, color doesn't apply
	double ri; // refractive index
	object3D() { ri = 1; }
	object3D(object3D* a) {
		attcoe = a->attcoe;
		ri = a->ri;
	}
	~object3D() {}

	void setAttCof(const double &a) { attcoe.r = attcoe.g = attcoe.b = abs(a); }
	void setAttCof(const double &r, const double &g, const double &b) { attcoe.r = abs(r), attcoe.g = abs(g), attcoe.b = abs(b); }
	void setIndex(const double &c) { ri = c; }

	// Given calculated intersection data and calculate refraction (since this always cost a lot of time)
	// R: calculated intersection data;  a: ray (must match intersection);  mi: refractive index of the other media;
	// rlr: reflect ratio calculated with Fresnel equations; NAN occurs in "refract" means total reflection (also rlr=1)
	virtual void refractData(const intersect &R, const ray &a, const double &mi, point &refract, double &rlr) {
		cout << "\aobject3D::refractData is called. This function should never be called. \n";
	}
	// R.ut=1: refract in (air->obj);  R.ut=0: refract out (obj->air);

	// inside => negative, outside => positive
	virtual double SDF(const point &A) { cout << "\aobject3D::SDF is called. This function should never be called. \n"; return NAN; }
	virtual bool contain(const point &A) { cout << "\aobject3D::inside is called. This function should never be called. \n"; return false; }
};


/* infinite large horizontal plane, with volumn beneath it */
/* DEBUG */
#define WaterSurface_Sign 0x00010000
class WaterSurface : public object3D {
public:
	double z_int;
	WaterSurface() {
		z_int = 0; ri = 1.33;
	}
	WaterSurface(const double &z_int) {
		this->z_int = z_int; ri = 1.33;
	}
	WaterSurface(const double &z_int, const double &c) {
		this->z_int = z_int, this->ri = c > 1 ? c : 1 / c;
	}
	WaterSurface(const double &z_int, const double &c, const double &ac_r, const double &ac_g, const double &ac_b) {
		this->z_int = z_int, this->ri = c > 1 ? c : 1 / c;
		this->attcoe = rgblight(ac_r, ac_g, ac_b);
	}
	inline point Max() {
		return point(INFINITY, INFINITY, z_int);
	}
	inline point Min() {
		return point(-INFINITY, -INFINITY, z_int);
	}

	void meet(intersect &R, const ray &a) {
		R.meet = 0;
		if (abs(a.orig.z - z_int) < ERR_EPSILON) return;
		if ((a.orig.z > z_int) ^ (a.dir.z < 0)) return;
		R.dist = (z_int - a.orig.z) / a.dir.z;
		R.intrs.x = a.dir.x * R.dist + a.orig.x;
		R.intrs.y = a.dir.y * R.dist + a.orig.y;
		R.intrs.z = z_int;
		R.dist = (R.intrs - a.orig).mod();
		R.reflect.x = a.dir.x, R.reflect.y = a.dir.y, R.reflect.z = -a.dir.z;
		R.meet = 1;
		return;
	}
	void refractData(const intersect &R, const ray &a, const double &mi, point &refract, double &rlr) {
		refract.x = a.dir.x, refract.y = a.dir.y;
		if (a.dir.z < 0) {
			refract.z = -sqrt((ri * ri - 1)*(a.dir.x*a.dir.x + a.dir.y*a.dir.y) + ri * ri * a.dir.z*a.dir.z);
		}
		else {
			refract.z = (1 - ri * ri)*(a.dir.x*a.dir.x + a.dir.y*a.dir.y) + a.dir.z*a.dir.z;
			if (refract.z <= 0) {
				refract.x = refract.y = refract.z = NAN;
				rlr = 1;
				return;
			}
			refract.z = sqrt(refract.z) / ri;
		}

		// Fresnel Equations
		double cci = abs(R.reflect.z) / R.reflect.mod(), cco = abs(refract.z) / refract.mod();
		double Rs, Rp;
		if (a.dir.z < 0) {
			Rs = (mi*cci - ri * cco) / (mi*cci + ri * cco); Rs *= Rs;
			Rp = (mi*cco - ri * cci) / (mi*cco + ri * cci); Rp *= Rp;
			rlr = 0.5 * (Rs + Rp);
			const_cast<intersect&>(R).ut = 1;
		}
		else {
			Rs = (ri*cci - mi * cco) / (ri*cci + mi * cco); Rs *= Rs;
			Rp = (ri*cco - mi * cci) / (ri*cco + mi * cci); Rp *= Rp;
			rlr = 0.5 * (Rs + Rp);
			const_cast<intersect&>(R).ut = 0;
		}
		return;
	}

	double SDF(const point &A) {
		return A.z - z_int;
	}
	bool contain(const point &A) {
		return A.z < z_int;
	}

	int telltype() const {
		return WaterSurface_Sign;
	}
	friend ostream& operator << (ostream& os, const WaterSurface &a) {
		os << "z=" << a.z_int << endl;
		return os;
	}
};

/* 3D sphere, "crystal ball" */
/* DEBUG */
#define Sphere3D_Sign 0x00010001
class sphere3D :public object3D {
public:
	point C; double r;
	sphere3D() :r(0) { ri = 1.5; }
	sphere3D(const point &C, const double &r) {
		this->C = C, this->r = r, this->ri = 1.5;
	}
	sphere3D(sphere3D* a) {
		attcoe = a->attcoe, ri = a->ri;
		C = a->C, r = a->r;
	}
	sphere3D(const initializer_list<double> &C, const double &r) {
		this->C = point(C), this->r = r, this->ri = 1.5;
	}

	inline void rotate(const double &rx, const double &ry, const double &rz) {
		matrix<double> M = matrix<double>({ {cos(rz),-sin(rz),0}, {sin(rz),cos(rz),0}, {0,0,1} })
			* matrix<double>({ {cos(ry),0,sin(ry)}, {0,1,0}, {-sin(ry),0,cos(ry)} })
			* matrix<double>({ {1,0,0}, {0,cos(rx),-sin(rx)}, {0,sin(rx),cos(rx)} });
		C *= M;
	}
	inline void operator *= (const double &t) {
		C *= t, r *= t;
	}
	inline friend sphere3D operator + (sphere3D B, const point &a) {
		B.C += a; return B;
	}
	inline void operator += (const point &a) {
		C += a;
	}
	inline friend sphere3D operator * (sphere3D B, const double &a) {
		B.C *= a, B.r *= a; return B;
	}

	inline point Max() {
		return point(C.x + r, C.y + r, C.z + r);
	}
	inline point Min() {
		return point(C.x - r, C.y - r, C.z - r);
	}
	void meet(intersect &R, const ray &a) {
		point p = C - a.orig;
		double pm = p.mod();
		if (abs(pm - r) < ERR_EPSILON) {
			(const_cast<ray&>(a)).orig += ERR_ZETA * a.dir;
			p = C - a.orig;
			pm = p.mod();
		}
		if (pm > r) {
			if (dot(p, a.dir) < 0) return;
			double sm = a.dir.mod();
			double d = cross(p, a.dir).mod() / sm;
			if (d > r) return;
			d *= d;
			R.dist = sqrt(dot(p, p) - d) - sqrt(r*r - d);
			if (R.dist < ERR_EPSILON) return;
			point s = a.dir * (R.dist / sm), n = p - s;
			R.intrs = s + a.orig;
			R.reflect = s - (2 * dot(s, n) / dot(n, n)) * n;
			R.meet = 1;
			R.ut = 1;
			//R.intrs += ERR_ZETA * R.reflect;
			return;
		}
		else {
			double sm = a.dir.mod();
			double d = cross(p, a.dir).mod() / sm;
			d *= d;
			R.dist = sqrt(dot(p, p) - d) + sqrt(r*r - d);
			point s = a.dir * (R.dist / sm), n = s - p;
			R.intrs = s + a.orig;
			R.reflect = s - (2 * dot(s, n) / dot(n, n)) * n;
			R.meet = 1;
			R.ut = 0;
			//R.intrs += ERR_ZETA * R.reflect;
		}
	}
	void refractData(const intersect &R, const ray &a, const double &mi, point &refract, double &rlr) {
		// https://graphics.stanford.edu/courses/cs148-10-summer/docs/2006--degreve--reflection_refraction.pdf
		point n = C - R.intrs;
		double mn = n.mod(), ms = a.dir.mod(), c1 = dot(n, a.dir) / (mn*ms);
		double n1 = mi, n2 = ri;
		if (R.ut == 0) n1 = ri, n2 = mi, c1 = -c1, n = -n; 	// ut = { 1: air->obj;  0: obj->air }
		double c = n1 / n2, c2 = sqrt(1 - c * c * (1 - c1 * c1));
		if (isnan(c2)) {
			refract.x = refract.y = refract.z = NAN; rlr = 1;
			return;
		}
		refract = c * a.dir - ((c*c1 - c2)*ms / mn)*n;
		double Rs = (n2*c1 - n1 * c2) / (n2*c1 + n1 * c2); Rs *= Rs;
		double Rp = (n2*c2 - n1 * c1) / (n2*c2 + n1 * c1); Rp *= Rp;
		rlr = 0.5*(Rs + Rp);
	}

	bool contain(const point &A) {
		return (C - A).mod() < r;
	}
	double SDF(const point &A) {
		return (C - A).mod() - r;
	}

	friend ostream& operator << (ostream& os, const sphere3D &a) {
		os << "Sphere(" << a.C << "," << a.r << ")";
		return os;
	}
	int telltype() const { return Sphere3D_Sign; }
};

#endif



#ifndef _INC_lightsource
#define _INC_lightsource

class lightsource :public object {
public:
	rgblight col;

	lightsource() { col.r = col.g = col.b = 1; }
	lightsource(const lightsource &a) { col = a.col; }
	~lightsource() {}

	void setcolor(const WebSafeColour &c) { col = color(c); }
	void setcolor(const pixel &c) { col = c; }
	void setcolor(const rgblight &c) { col = c; }
	void setcolor(const unsigned &hex) { col = pixel(hex); }

	// Note that in all derivative classes of lightsource, the meet function set "ut" as the cosine of angle of incidence

};

#define SphereBulb_Sign 0x01000000
class spherebulb : public lightsource {
public:
	point C; double r;
	~spherebulb() {}
	spherebulb(const point &C, const double &r) {
		this->C = C, this->r = r, col.r = col.g = col.b = 1;
	}
	spherebulb(const point &C, const double &r, const rgblight &col) {
		this->C = C, this->r = r, this->col = col;
	}
	spherebulb(const spherebulb &a) { C = a.C, r = a.r, col = a.col; }
	spherebulb() :r(0) { col.r = col.g = col.b = 1; }

	void meet(intersect &R, const ray &a) {
		R.meet = 0;
		point p = C - a.orig;
		if (dot(p, a.dir) < 0) return;
		double sm = a.dir.mod();
		double d = cross(p, a.dir).mod() / sm;
		if (d > r) return;
		d *= d;
		R.dist = sqrt(dot(p, p) - d) - sqrt(r*r - d);
		point n = a.dir * (R.dist / sm), t = n - p;
		R.intrs = n + a.orig;
		R.reflect = n + abs(2 * dot(n, t) / dot(t, t)) * t;

		R.ut = abs(dot(t, a.dir)) / (t.mod()*a.dir.mod());
		R.meet = 1;
		return;
	}
	inline point Max() {
		return point(C.x + r, C.y + r, C.z + r);
	}
	inline point Min() {
		return point(C.x - r, C.y - r, C.z - r);
	}
	friend ostream& operator << (ostream& os, const spherebulb &a) {
		os << "Sphere(" << a.C << "," << a.r << ")";
		return os;
	}
	int telltype() const { return SphereBulb_Sign; }

};


#define RectBulb_Sign 0x01000001
class rectbulb : public lightsource {
public:
	point O, A, B;
	rectbulb() {}
	rectbulb(rectbulb* a) {
		col = a->col;
		A = a->A, B = a->B, O = a->O;
	}
	rectbulb(const point &O, const point &A, const point &B) {
		this->O = O, this->A = O + A, this->B = O + B;
	}
	rectbulb(const point &O, const point &A, const point &B, bool absolute) {
		this->O = O, this->A = A, this->B = B;
		if (!absolute) this->A += O, this->B += O;
	}
	rectbulb(const initializer_list<double> &O, const initializer_list<double> &A, const initializer_list<double> &B) {
		this->O = O, this->A = A, this->B = B;
		this->A += O, this->B += O;
	}
	void operator = (const initializer_list<initializer_list<double>> &V) {
		O = *V.begin(), A = *(V.begin() + 1), B = *(V.begin() + 2);
	}
	~rectbulb() {}
	void meet(intersect &R, const ray &a) {
		R.meet = 0;
		point E1 = A - O, E2 = B - O, T, P = cross(a.dir, E2), Q;
		double det = dot(E1, P);
		if (abs(det) < ERR_EPSILON) return;
		T = a.orig - O;
		double t, u, v;
		u = dot(T, P) / det;
		if (u < 0.0 || u > 1.0) return;
		Q = cross(T, E1);
		v = dot(a.dir, Q) / det;
		if (v < 0.0 || v > 1.0) return;
		t = dot(E2, Q) / det;
		if (t < ERR_EPSILON) return;
		R.meet = 1;
		R.dist = t * a.dir.mod();
		R.intrs = (1 - u - v)*O + u * A + v * B;
		point OA = A - R.intrs, OB = B - R.intrs;
		point ON = cross(OA, OB);
		ON *= dot(a.dir, ON) / dot(ON, ON);		// NAN occurs when 0/0
		R.reflect = a.dir - 2 * ON;
		R.ut = dot(a.dir, ON) / (a.dir.mod()*ON.mod());
		return;
	}
	inline point Max() {
		return point(max({ O.x, A.x, B.x, A.x + B.x - O.x }),
			max({ O.y, A.y, B.y, A.y + B.y - O.y }), max({ O.y, A.y, B.y,A.y + B.y - O.y }));
	}
	inline point Min() {
		return point(min({ O.x, A.x, B.x, A.x + B.x - O.x }),
			min({ O.y, A.y, B.y, A.y + B.y - O.y }), min({ O.y, A.y, B.y,A.y + B.y - O.y }));
	}
	inline void operator += (const point &a) {
		A += a, B += a, O += a;
	}
	friend ostream& operator << (ostream& os, const rectbulb &a) {
		os << "Surface(" << noshowpos << a.O.x << "*(1-u-v)" << showpos << a.A.x << "*u" << showpos << a.B.x << "*v" << ", "
			<< noshowpos << a.O.y << "*(1-u-v)" << showpos << a.A.y << "*u" << showpos << a.B.y << "*v" << ", "
			<< noshowpos << a.O.z << "*(1-u-v)" << showpos << a.A.z << "*u" << showpos << a.B.z << "*v" << ", "
			<< "u, 0, 1, v, 0, 1)" << noshowpos;
		return os;
	}
	int telltype() const { return RectBulb_Sign; }
};


#endif




ostream& operator << (ostream& os, const object &a) {
	switch (a.telltype()) {
	case Plane_Sign: {
		cout << *((plane*)(&a));
		break;
	}
	case Triangle_Sign: {
		cout << *((triangle*)(&a));
		break;
	}
	case Parallelogram_Sign: {
		cout << *((parallelogram*)(&a));
		break;
	}
	case Circle_Sign: {
		cout << *((circle*)(&a));
		break;
	}
	case Sphere_Sign: {
		cout << *((sphere*)(&a));
		break;
	}
	case Cylinder_Sign: {
		cout << *((cylinder*)(&a));
		break;
	}
	case WaterSurface_Sign: {
		cout << *((WaterSurface*)(&a));
		break;
	}
	case Sphere3D_Sign: {
		cout << *((sphere3D*)(&a));
		break;
	}
	default: {
		cout << "\aError! a.telltype() return 0x" << hex << a.telltype() << ". ";
	}
	}
	return os;
}

