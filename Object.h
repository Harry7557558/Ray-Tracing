#pragma once

#include "objdef.h"
#include <vector>
#include <initializer_list>

using namespace std;


/* object parent class */
#define Object_Sign 0xFFFFFFFF
class object {
public:

	object() { /*cout << "\aobject::Constructor is called. This function should never be called. \n";*/ }
	object(object* a) {
		//cout << "\aobject::CpyConstructor is called. This function should never be called. \n";
	}
	~object() {}

	virtual object* copy() const {
		return new object;
	}
	virtual void init() {}	// some objects need to be initialized before rendering

	virtual void setcolor(const WebSafeColour &c) {}
	virtual void setcolor(const pixel &c) {}
	virtual void setcolor(const rgblight &c) {}
	virtual void setcolor(const unsigned &c) {}

	/* Intersection Test */
	virtual void meet(intersect &R, const ray &a) const { WARN("object::meet is called. This function should never be called."); return; }
	// a.dir must be a unit vector


	virtual point Max() const { WARN("object::Max is called. This function should never be called."); return point(); }
	virtual point Min() const { WARN("object::Min is called. This function should never be called."); return point(); }

	/* return the type of an object, undefined is -1 */
	virtual unsigned telltype() const { WARN("object::telltype is called. This function should never be called."); return Object_Sign; }
	friend ostream& operator << (ostream& os, const object &a) {
		a.print(os); os << defaultfloat; return os;
	}
	virtual void print(ostream& os) const { os << "Object Parent Class"; }
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

	object* copy() const {
		return new objectSF;
	}

	void setcolor(const WebSafeColour &c) { reflect = color(c); }
	void setcolor(const pixel &c) { reflect = c; }
	void setcolor(const rgblight &c) { reflect = c; }
	void setcolor(const unsigned &hex) { reflect = pixel(hex); }

	void print(ostream& os) const { os << "objectSF parent class"; }
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
	~plane() {}
	object* copy() const {
		return new plane(*this);
	}

	point Max() const {
		if (N.x == 0 && N.y == 0 && N.z == 0) return point(0, 0, 0);
		if (N.x == 0 && N.y == 0) return point(INFINITY, INFINITY, D / N.z + 0.01);
		if (N.x == 0 && N.z == 0) return point(INFINITY, D / N.y + 0.01, INFINITY);
		if (N.y == 0 && N.z == 0) return point(D / N.x + 0.01, INFINITY, INFINITY);
		return point(INFINITY, INFINITY, INFINITY);
	}
	point Min() const {
		if (N.x == 0 && N.y == 0 && N.z == 0) return point(0, 0, 0);
		if (N.x == 0 && N.y == 0) return point(-INFINITY, -INFINITY, D / N.z - 0.01);
		if (N.x == 0 && N.z == 0) return point(-INFINITY, D / N.y - 0.01, -INFINITY);
		if (N.y == 0 && N.z == 0) return point(D / N.x - 0.01, -INFINITY, -INFINITY);
		return point(-INFINITY, -INFINITY, -INFINITY);
	}

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

	void meet(intersect &R, const ray &a) const {
		R.meet = 0;
		double t = (D - dot(N, a.orig)) / dot(N, a.dir);
		if (t < ERR_EPSILON || t > ERR_UPSILON) return;
		R.intrs = t * a.dir + a.orig;
		R.dist = t;
		R.reflect = 2 * (-dot(a.dir, N) / dot(N, N)) * N + a.dir;
		R.meet = 1;
		return;
	}

	void print(ostream& os) const {
		os << "Plane_{" << OS_Pointer << (unsigned)this << "}: ";
		if (abs(N.x) > ERR_EPSILON) os << N.x << "*x";
		if (abs(N.y) > ERR_EPSILON) os << showpos << N.y << "*y";
		if (abs(N.z) > ERR_EPSILON) os << showpos << N.z << "*z";
		os << "=" << noshowpos << D;
	}
	unsigned telltype() const { return Plane_Sign; }
};


/* triangle class */
#define Triangle_Sign 0x00000001
class triangle :public objectSF {
public:
	point A, B, C;
	triangle() {}
	triangle(const triangle& a) {
		reflect = a.reflect;
		A = a.A, B = a.B, C = a.C;
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
	object* copy() const {
		return new triangle(*this);
	}

	void meet(intersect &R, const ray &a) const {
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
		R.dist = t;
		R.intrs = (1 - u - v)*A + u * B + v * C;
		point ON = cross(E1, E2);
		ON *= dot(a.dir, ON) / dot(ON, ON);
		R.reflect = a.dir - 2 * ON;
		return;
	}
	point Max() const { return point(max({ A.x, B.x, C.x }), max({ A.y, B.y, C.y }), max({ A.z, B.z, C.z })); }
	point Min() const { return point(min({ A.x, B.x, C.x }), min({ A.y, B.y, C.y }), min({ A.z, B.z, C.z })); }
	inline friend triangle operator * (matrix<double> M, triangle a) {
		a.A = M * a.A, a.B = M * a.B, a.C = M * a.C;
		return a;
	}
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
	void print(ostream& os) const {
		os << "Triangle_{" << OS_Pointer << (unsigned)this << "}: ";
		os << "Surface(if (u + v < 1, " << noshowpos << A.x << "*(1 - u - v)" << showpos << B.x << "*u" << showpos << C.x << "*v" << "), "
			<< noshowpos << A.y << "*(1-u-v)" << showpos << B.y << "*u" << showpos << C.y << "*v" << ", "
			<< noshowpos << A.z << "*(1-u-v)" << showpos << B.z << "*u" << showpos << C.z << "*v" << ", "
			<< "u, 0, 1, v, 0, 1)" << noshowpos;
	}
	unsigned telltype() const { return Triangle_Sign; }
};

/* parallelogram class */
#define Parallelogram_Sign 0x00000002
class parallelogram :public objectSF {
public:
	point O, A, B;
	parallelogram() {}
	parallelogram(const parallelogram& a) {
		reflect = a.reflect;
		A = a.A, B = a.B, O = a.O;
	}
	parallelogram(const point &O, const point &A, const point &B) {
		this->O = O, this->A = A, this->B = B;
	}
	parallelogram(const point &O, const point &A, const point &B, bool absolute) {
		this->O = O, this->A = A, this->B = B;
		if (absolute) this->A -= O, this->B -= O;
	}
	parallelogram(const initializer_list<double> &O, const initializer_list<double> &A, const initializer_list<double> &B) {
		this->O = O, this->A = A, this->B = B;
	}
	~parallelogram() {}
	object* copy() const {
		return new parallelogram(*this);
	}

	void meet(intersect &R, const ray &a) const {
		R.meet = 0;
		point T, P = cross(a.dir, B), Q;
		double det = dot(A, P);
		if (abs(det) < ERR_EPSILON) return;
		T = a.orig - O;
		double t, u, v;
		u = dot(T, P) / det;
		if (u < 0.0 || u > 1.0) return;
		Q = cross(T, A);
		v = dot(a.dir, Q) / det;
		if (v < 0.0 || v > 1.0) return;	// For triangles: v < 0.0 || u + v > 1.0
		t = dot(B, Q) / det;
		if (t < ERR_EPSILON) return;
		R.meet = 1;
		R.dist = t;
		R.intrs = O + u * A + v * B;
		point ON = cross(A, B);
		ON *= dot(a.dir, ON) / dot(ON, ON);
		R.reflect = a.dir - 2 * ON;
		return;
	}
	point Max() const {
		return point(max({ O.x, O.x + A.x, O.x + B.x, O.x + A.x + B.x }),
			max({ O.y, O.y + A.y, O.y + B.y, O.y + A.y + B.y }), max({ O.z, O.z + A.z, O.z + B.z, O.z + A.z + B.z }));
	}
	point Min() const {
		return point(min({ O.x, O.x + A.x, O.x + B.x, O.x + A.x + B.x }),
			min({ O.y, O.y + A.y, O.y + B.y, O.y + A.y + B.y }), min({ O.z, O.z + A.z, O.z + B.z, O.z + A.z + B.z }));
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
		t.O += a; return t;
	}
	inline void operator += (const point &a) {
		O += a;
	}
	void print(ostream& os) const {
		os << "Parallelogram_{" << OS_Pointer << (unsigned)this << "}: ";
		os << "Surface(" << noshowpos << O.x << showpos << A.x << "*u" << showpos << B.x << "*v" << ", "
			<< noshowpos << O.y << showpos << A.y << "*u" << showpos << B.y << "*v" << ", "
			<< noshowpos << O.z << showpos << A.z << "*u" << showpos << B.z << "*v" << ", "
			<< "u, 0, 1, v, 0, 1)" << noshowpos;
	}
	unsigned telltype() const { return Parallelogram_Sign; }
};


/* spacial circle class */
#define Circle_Sign 0x00000003
class circle :public objectSF {
public:
	point C; double r; double rx, ry, rz;	// First rotate x, then rotate z    R = Rz * Rx		// Clockwise??!!
	circle() :r(0), rx(0), ry(0), rz(0) {}
	circle(const circle& a) {
		reflect = a.reflect;
		C = a.C, r = a.r; rx = a.rx, ry = a.ry, rz = a.rz;
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
	~circle() {}
	object* copy() const {
		return new circle(*this);
	}

	inline friend circle operator + (circle c, const point &p) {
		c.C += p; return c;
	}
	inline void operator += (const point &a) {
		C += a;
	}
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
	point Max() const {
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
	point Min() const {
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
	void meet(intersect &R, const ray &a) const {
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
		R.intrs.x = s.dir.x*t + s.orig.x, R.intrs.y = s.dir.y*t + s.orig.y, R.intrs.z = 0;
		R.dist = R.intrs.x * R.intrs.x + R.intrs.y * R.intrs.y;
		if (R.dist > r*r) return;
		R.dist = (R.intrs - s.orig).mod();
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
	void print(ostream& os) const {
		os << "Circle_{" << OS_Pointer << (unsigned)this << "}: ";
		os << "Rotate(Rotate(Rotate(Surface(";
		os << "v*cos(u)" << showpos << C.x << ", ";
		os << "v*sin(u)" << showpos << C.y << ", ";
		os << noshowpos << C.z << ", ";
		os << "u, 0, 2*pi, v, 0, " << noshowpos << r << "), ";
		os << noshowpos << rx << ", " << point(C) << ", xAxis), "
			<< noshowpos << ry << ", " << point(C) << ", yAxis), "
			<< noshowpos << rz << ", " << point(C) << ", zAxis)";
	}
	unsigned telltype() const { return Circle_Sign; }
};


/* cylinder surface class, without bottom */
#define Cylinder_Sign 0x00000004
class cylinder :public objectSF {
public:
	point C; double r, h; double rx, ry, rz;
	cylinder() :r(0), h(0), rx(0), ry(0), rz(0) {}
	cylinder(const cylinder& a) {
		reflect = a.reflect;
		C = a.C, r = a.r, h = a.h; rx = a.rx, ry = a.ry, rz = a.rz;
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
	~cylinder() {}
	object* copy() const {
		return new cylinder(*this);
	}

	inline void operator += (const point &a) {
		C += a;
	}
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
	point Max() const {
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
	point Min() const {
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
	void meet(intersect &R, const ray &a) const {
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
	void print(ostream& os) const {
		os << "Cylinder_{" << OS_Pointer << (unsigned)this << "}: ";
		os << "Rotate(Rotate(Rotate(Surface(";
		os << noshowpos << r << "*cos(u)" << showpos << C.x << ", ";
		os << noshowpos << r << "*sin(u)" << showpos << C.y << ", ";
		os << "v" << showpos << C.z << ", ";
		os << "u, 0, 2*pi, v, 0, " << noshowpos << h << "), ";
		os << noshowpos << rx << ", " << point(C) << ", xAxis), "
			<< noshowpos << ry << ", " << point(C) << ", yAxis), "
			<< noshowpos << rz << ", " << point(C) << ", zAxis)";
	}
	unsigned telltype() const { return Cylinder_Sign; }
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
	sphere(const sphere& a) {
		reflect = a.reflect;
		C = a.C, r = a.r;
	}
	sphere(const initializer_list<double> &C, const double &r) {
		this->C = point(C), this->r = r;
	}
	~sphere() {}
	object* copy() const {
		return new sphere(*this);
	}

	void meet(intersect &R, const ray &a) const {
		R.meet = 0;
		point P = C - a.orig;
		if (dot(P, a.dir) < 0) return;
		double d = cross(P, a.dir).mod();
		if (d > r) return;
		d *= d;
		R.dist = sqrt(dot(P, P) - d) - sqrt(r*r - d);
		point S = a.dir * R.dist, N = S - P;
		R.intrs = S + a.orig;
		R.reflect = S - (2 * dot(S, N) / dot(N, N)) * N;
		R.meet = 1;
		return;
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
	inline friend sphere operator + (sphere B, const point &a) {
		B.C += a; return B;
	}
	inline void operator += (const point &a) {
		C += a;
	}
	point Max() const {
		return point(C.x + r, C.y + r, C.z + r);
	}
	point Min() const {
		return point(C.x - r, C.y - r, C.z - r);
	}
	inline friend sphere operator * (sphere B, const double &a) {
		B.C *= a, B.r *= a; return B;
	}

	void print(ostream& os) const {
		os << "Sphere_{" << OS_Pointer << (unsigned)this << "}: Sphere(" << C << "," << r << ")";
	}
	unsigned telltype() const { return Sphere_Sign; }
};

/* ring class */
#define Torus_Sign 0x00000006
class torus :public objectSF {
public:
	point C; double R, r; double rx, ry, rz;
	torus() :R(0), r(0), rx(0), ry(0), rz(0) {}
	torus(const torus &a) {
		C = a.C, R = a.R, r = a.r, rx = a.rx, ry = a.ry, rz = a.rz;
	}
	torus(const point &C, const double &R, const double &r) {
		this->C = C, this->R = R, this->r = r, rx = ry = rz = 0;
	}
	torus(const point &C, const double &R, const double &r, const double &rx, const double &rz) {
		this->C = C, this->R = R, this->r = r, this->rx = rx, this->ry = 0, this->rz = rz;
	}
	torus(const point &C, const double &R, const double &r, const double &rx, const double &ry, const double &rz) {
		this->C = C, this->R = R, this->r = r, this->rx = rx, this->ry = ry, this->rz = rz;
	}
	~torus() {}
	object* copy() const {
		return new torus(*this);
	}

	point Max() const {
		return C + point(R + r, R + r, r);
	}
	point Min() const {
		return C - point(R + r, R + r, r);
	}
	void meet(intersect &Res, const ray &a) const {
		Res.meet = false;
		point P = a.orig - C;
		point s = P.mod()*a.dir;	// issues in numerically solving quartic equations
		double x2y2 = P.x*P.x + P.y*P.y, z2 = P.z*P.z, a2b2 = s.x*s.x + s.y*s.y, c2 = s.z*s.z, axby = P.x*s.x + P.y*s.y, cz = P.z*s.z,
			x2y2z2 = x2y2 + z2, a2b2c2 = a2b2 + c2, axbycz = axby + cz, R2pr2 = R * R + r * r, R2mr2 = R * R - r * r;
		double t4 = a2b2c2 * a2b2c2, t3 = 4 * a2b2c2 * axbycz, t2 = 2 * a2b2c2 * x2y2z2 + 4 * axbycz * axbycz, t1 = 4 * x2y2z2 * axbycz, t0 = x2y2z2 * x2y2z2;
		t2 += 2 * (R2mr2*c2 - R2pr2 * a2b2), t1 += 4 * (cz*R2mr2 - axby * R2pr2), t0 += 2 * (R2mr2*z2 - R2pr2 * x2y2) + R2mr2 * R2mr2;
		Res.dist = solveQuartic(t4, t3, t2, t1, t0);
		if (isnan(Res.dist)) return;
		Res.dist *= P.mod();
		Res.intrs = P + Res.dist * a.dir;
		Res.ut = atan2(Res.intrs.y, Res.intrs.x), Res.vt = asin(Res.intrs.z / r); if (isnan(Res.vt)) Res.vt = Res.intrs.z > 0 ? PI / 2 : -PI / 2;
		point N = point(cos(Res.ut)*cos(Res.vt), sin(Res.ut)*cos(Res.vt), sin(Res.vt));
		Res.reflect = a.dir - 2 * dot(a.dir, N)*N;
		Res.intrs += C;
		Res.meet = true;
	}
	void print(ostream& os) const {
		os << "Ring_{" << OS_Pointer << (unsigned)this << "}: ";
		os << "Surface(" << "cos(u)*(" << noshowpos << R << showpos << r << "*cos(v))" << showpos << C.x << ", "
			<< "sin(u)*(" << noshowpos << R << showpos << r << "*cos(v))" << showpos << C.y << "," << noshowpos << r << "*sin(v)" << showpos << C.z << ",u,0,2*pi,v,0,2*pi)";
	}
	unsigned telltype() const { return Cylinder_Sign; }
};

#endif

#ifndef _INC_object2D_dif
#define _INC_object2D_dif


/* Opacity surface with diffuse reflection */
class objectSF_dif : public objectSF {
public:
	objectSF_dif() {}
	objectSF_dif(const objectSF_dif &a) { this->reflect = a.reflect; }
	~objectSF_dif() {}
	object* copy() const {
		return new objectSF_dif(*this);
	}

	void setvar(double d) { }

	void print(ostream& os) const { os << "objectSF_dif parent class"; }

};

#define Plane_Dif_Sign 0x00000100
class plane_dif : public objectSF_dif/*, public plane*/ {
public:
	point N; double D;	// Ax+By+Cz=D
	plane_dif() { D = 0, N.z = 1; }
	// normal N, through origin
	plane_dif(const point &N) {
		this->N = N, D = 0;
	}
	// horizontal plain, z=D
	plane_dif(const double &D) {
		N.z = 1, this->D = D;
	}
	// point P and normal N
	plane_dif(const point &P, const point &N) {
		this->N = N; D = dot(P, N);
	}
	// through three points
	plane_dif(const point &A, const point &B, const point &C) {
		N = cross(B - A, C - A); D = dot(A, N);
	}
	plane_dif(const double &x_int, const double &y_int, const double &z_int) {
		N.x = y_int * z_int, N.y = x_int * z_int, N.z = x_int * y_int;
		D = x_int * y_int * z_int;
	}
	plane_dif(const plane_dif &p) :N(p.N), D(p.D) {
		reflect = p.reflect;
	}
	~plane_dif() {}
	object* copy() const {
		return new plane_dif(*this);
	}

	point Max() const {
		if (N.x == 0 && N.y == 0 && N.z == 0) return point(0, 0, 0);
		if (N.x == 0 && N.y == 0) return point(INFINITY, INFINITY, D / N.z + 0.01);
		if (N.x == 0 && N.z == 0) return point(INFINITY, D / N.y + 0.01, INFINITY);
		if (N.y == 0 && N.z == 0) return point(D / N.x + 0.01, INFINITY, INFINITY);
		return point(INFINITY, INFINITY, INFINITY);
	}
	point Min() const {
		if (N.x == 0 && N.y == 0 && N.z == 0) return point(0, 0, 0);
		if (N.x == 0 && N.y == 0) return point(-INFINITY, -INFINITY, D / N.z - 0.01);
		if (N.x == 0 && N.z == 0) return point(-INFINITY, D / N.y - 0.01, -INFINITY);
		if (N.y == 0 && N.z == 0) return point(D / N.x - 0.01, -INFINITY, -INFINITY);
		return point(-INFINITY, -INFINITY, -INFINITY);
	}

	void meet(intersect &R, const ray &a) const {
		R.meet = 0;
		double t = (D - dot(N, a.orig)) / dot(N, a.dir);
		if (t < ERR_EPSILON || t > ERR_UPSILON) return;
		R.intrs = t * a.dir + a.orig;
		R.dist = t;
		R.reflect = N / N.mod();
		if (dot(a.dir, N) > 0) R.reflect = -R.reflect;
		R.meet = 1;
		return;
	}

	unsigned telltype() const { return Plane_Dif_Sign; }
};

#define Parallelogram_Dif_Sign 0x00000101
class parallelogram_dif :public objectSF_dif {
public:
	point O, A, B;
	parallelogram_dif() {}
	parallelogram_dif(const parallelogram& a) {
		reflect = a.reflect;
		A = a.A, B = a.B, O = a.O;
	}
	parallelogram_dif(const parallelogram_dif& a) {
		reflect = a.reflect;
		A = a.A, B = a.B, O = a.O;
	}
	parallelogram_dif(const point &O, const point &A, const point &B) {
		this->O = O, this->A = A, this->B = B;
	}
	parallelogram_dif(const point &O, const point &A, const point &B, bool absolute) {
		this->O = O, this->A = A, this->B = B;
		if (absolute) this->A -= O, this->B -= O;
	}
	~parallelogram_dif() {}
	object* copy() const {
		return new parallelogram_dif(*this);
	}
	void meet(intersect &R, const ray &a) const {
		R.meet = 0;
		point T, P = cross(a.dir, B), Q;
		double det = dot(A, P);
		if (abs(det) < ERR_EPSILON) return;
		T = a.orig - O;
		double t, u, v;
		u = dot(T, P) / det;
		if (u < 0.0 || u > 1.0) return;
		Q = cross(T, A);
		v = dot(a.dir, Q) / det;
		if (v < 0.0 || v > 1.0) return;
		t = dot(B, Q) / det;
		if (t < ERR_EPSILON) return;
		R.meet = 1;
		R.dist = t;
		R.intrs = O + u * A + v * B;
		R.reflect = cross(A, B);
		if (dot(a.dir, R.reflect) > 0) R.reflect = -R.reflect;
		R.reflect /= R.reflect.mod();
		return;
	}
	point Max() const {
		return point(max({ O.x, O.x + A.x, O.x + B.x, O.x + A.x + B.x }),
			max({ O.y, O.y + A.y, O.y + B.y, O.y + A.y + B.y }), max({ O.z, O.z + A.z, O.z + B.z, O.z + A.z + B.z }));
	}
	point Min() const {
		return point(min({ O.x, O.x + A.x, O.x + B.x, O.x + A.x + B.x }),
			min({ O.y, O.y + A.y, O.y + B.y, O.y + A.y + B.y }), min({ O.z, O.z + A.z, O.z + B.z, O.z + A.z + B.z }));
	}

	void print(ostream& os) const {
		os << "Surface(" << noshowpos << O.x << showpos << A.x << "*u" << showpos << B.x << "*v" << ", "
			<< noshowpos << O.y << showpos << A.y << "*u" << showpos << B.y << "*v" << ", "
			<< noshowpos << O.z << showpos << A.z << "*u" << showpos << B.z << "*v" << ", "
			<< "u, 0, 1, v, 0, 1)" << noshowpos;
	}
	unsigned telltype() const { return Parallelogram_Dif_Sign; }
};

#define Triangle_Dif_Sign 0x00000102
class triangle_dif :public objectSF_dif {
public:
	point A, B, C;
	triangle_dif() {}
	triangle_dif(const triangle& a) {
		reflect = a.reflect;
		A = a.A, B = a.B, C = a.C;
	}
	triangle_dif(const triangle_dif& a) {
		reflect = a.reflect;
		A = a.A, B = a.B, C = a.C;
	}
	triangle_dif(const point &A, const point &B, const point &C) {
		this->A = A, this->B = B, this->C = C;
	}
	triangle_dif(const initializer_list<double> &A, const initializer_list<double> &B, const initializer_list<double> &C) {
		this->A = A, this->B = B, this->C = C;
	}
	~triangle_dif() {}
	object* copy() const {
		return new triangle_dif(*this);
	}
	void meet(intersect &R, const ray &a) const {
		// Algorithm: http://www.graphics.cornell.edu/pubs/1997/MT97.pdf
		R.meet = 0;
		point E1 = B - A, E2 = C - A, T, P = cross(a.dir, E2), Q;
		double det = dot(E1, P);
		if (abs(det) < ERR_EPSILON) return;
		T = a.orig - A;
		double t, u, v;
		u = dot(T, P) / det;
		if (u < 0.0 || u > 1.0) return;
		Q = cross(T, E1);
		v = dot(a.dir, Q) / det;
		if (v < 0.0 || u + v > 1.0) return;
		t = dot(E2, Q) / det;
		if (t < ERR_EPSILON) return;
		R.meet = 1;
		R.dist = t;
		R.intrs = (1 - u - v)*A + u * B + v * C;
		R.reflect = cross(E1, E2);
		if (dot(a.dir, R.reflect) > 0) R.reflect = -R.reflect;
		R.reflect /= R.reflect.mod();
		return;
	}
	point Max() const { return point(max({ A.x, B.x, C.x }), max({ A.y, B.y, C.y }), max({ A.z, B.z, C.z })); }
	point Min() const { return point(min({ A.x, B.x, C.x }), min({ A.y, B.y, C.y }), min({ A.z, B.z, C.z })); }
	inline void operator += (const point &a) {
		A += a, B += a, C += a;
	}
	void print(ostream& os) const {
		os << "Polyline(" << A << "," << B << "," << C << "," << A << ")";
	}
	unsigned telltype() const { return Triangle_Dif_Sign; }
};

#endif

#ifndef _INC_object2D_col
#define _INC_object2D_col

// Opacity smooth surface, different part using different colors
class objectSF_col : public objectSF {
public:
	objectSF_col() {}
	objectSF_col(const objectSF_col &a) {}
	~objectSF_col() {}
	object* copy() const {
		return new objectSF_col(*this);
	}

	// get color with calculated intersection data
	virtual void getcol(const intersect &R, rgblight &c) { WARN("objectSF_col::getcol is called. This function should never be called."); }

	void print(ostream& os) const { os << "objectSF_col class"; }
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
	~plane_grid() {}
	object* copy() const {
		return new plane_grid(*this);
	}

	point Max() const {
		return point(INFINITY, INFINITY, z_int);
	}
	point Min() const {
		return point(-INFINITY, -INFINITY, z_int);
	}

	void meet(intersect &R, const ray &a) const {
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

	unsigned telltype() const { return Plane_Grid_Sign; }
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
	~bitmap_inc() {}
	object* copy() const {
		return new bitmap_inc(*this);
	}

	point Max() const {
		return point(max({ O.x, O.x + A.x, O.x + B.x, O.x + A.x + B.x }),
			max({ O.y, O.y + A.y, O.y + B.y, O.y + A.y + B.y }), max({ O.z, O.z + A.z, O.z + B.z, O.z + A.z + B.z }));
	}
	point Min() const {
		return point(min({ O.x, O.x + A.x, O.x + B.x, O.x + A.x + B.x }),
			min({ O.y, O.y + A.y, O.y + B.y, O.y + A.y + B.y }), min({ O.z, O.z + A.z, O.z + B.z, O.z + A.z + B.z }));
	}


	void meet(intersect &R, const ray &a) const {
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
		R.dist = t;
		R.intrs = O + R.ut * A + R.vt * B;
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

	unsigned telltype() const { return Bitmap_Inc_Sign; }
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
	object3D(const object3D& a) {
		attcoe = a.attcoe;
		ri = a.ri;
	}
	~object3D() {}
	object* copy() const {
		return new object3D(*this);
	}

	void setAttCof(const double &a) { attcoe.r = attcoe.g = attcoe.b = abs(a); }
	void setAttCof(const double &r, const double &g, const double &b) { attcoe.r = abs(r), attcoe.g = abs(g), attcoe.b = abs(b); }
	void setIndex(const double &c) { ri = c; }

	// Given calculated intersection data and calculate refraction (since this always cost a lot of time)
	// R: calculated intersection data;  a: ray (must match intersection);  mi: refractive index of the other media;
	// rlr: reflect ratio calculated with Fresnel equations; NAN occurs in "refract" means total reflection (also rlr=1)
	virtual void refractData(const intersect &R, const ray &a, const double &mi, point &refract, double &rlr) const {
		WARN("\aobject3D::refractData is called. This function should never be called.");
	}
	// R.ut=1: refract in (air->obj);  R.ut=0: refract out (obj->air);

	// inside => negative, outside => positive
	virtual double SDF(const point &A) const { WARN("\aobject3D::SDF is called. This function should never be called."); return NAN; }

	virtual bool contain(const point &A) const { WARN("\aobject3D::inside is called. This function should never be called."); return false; }

	void print(ostream& os) const { os << "object3D parent class"; }
};


/* infinite large horizontal plane, with volumn beneath it */
#define WaterSurface_Sign 0x00010000
class WaterSurface : public object3D {
public:
	double z_int;
	WaterSurface() {
		z_int = 0; ri = 1.33;
	}
	WaterSurface(const WaterSurface* &a) {
		z_int = a->z_int;
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
	point Max() const {
		return point(INFINITY, INFINITY, z_int);
	}
	point Min() const {
		return point(-INFINITY, -INFINITY, z_int);
	}
	object* copy() const {
		return new WaterSurface(*this);
	}
	~WaterSurface() {}

	void meet(intersect &R, const ray &a) const {
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
	void refractData(const intersect &R, const ray &a, const double &mi, point &refract, double &rlr) const {
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
		double cci = abs(R.reflect.z), cco = abs(refract.z);
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

	double SDF(const point &A) const {
		return A.z - z_int;
	}
	bool contain(const point &A) const {
		return A.z < z_int;
	}

	unsigned telltype() const {
		return WaterSurface_Sign;
	}
	void print(ostream& os) const {
		os << "z=" << z_int;
	}
};

/* 3D sphere, "crystal ball" */
#define Sphere3D_Sign 0x00010001
class sphere3D :public object3D {
public:
	point C; double r;
	sphere3D() :r(0) { ri = 1.5; }
	sphere3D(const point &C, const double &r) {
		this->C = C, this->r = r, this->ri = 1.5;
	}
	sphere3D(const sphere3D& a) {
		attcoe = a.attcoe, ri = a.ri;
		C = a.C, r = a.r;
	}
	sphere3D(const initializer_list<double> &C, const double &r) {
		this->C = point(C), this->r = r, this->ri = 1.5;
	}
	object* copy() const {
		return new sphere3D(*this);
	}
	~sphere3D() {}

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

	point Max() const {
		return point(C.x + r, C.y + r, C.z + r);
	}
	point Min() const {
		return point(C.x - r, C.y - r, C.z - r);
	}
	void meet(intersect &R, const ray &a) const {
		point p = C - a.orig;
		double pm = p.mod();
		if (abs(pm - r) < ERR_EPSILON) {
			// This problem mainly occurs on hook surfaces, less common on linear planes
			// about 70% occurs when pm > r
			(const_cast<ray&>(a)).orig += ERR_ZETA * a.dir;
			p = C - a.orig;
			pm = p.mod();
		}
		if (pm > r) {
			if (dot(p, a.dir) < 0) return;
			double d = cross(p, a.dir).mod();
			if (d > r) return;
			d *= d;
			R.dist = sqrt(dot(p, p) - d) - sqrt(r*r - d);
			if (R.dist < ERR_EPSILON) return;
			point s = a.dir * R.dist, n = p - s;
			R.intrs = s + a.orig;
			R.reflect = s - (2 * dot(s, n) / dot(n, n)) * n;
			R.meet = 1;
			R.ut = 1;
		}
		else {
			double d = cross(p, a.dir).mod();
			d *= d;
			R.dist = sqrt(dot(p, p) - d) + sqrt(r*r - d);
			point s = a.dir * R.dist, n = s - p;
			R.intrs = s + a.orig;
			R.reflect = s - (2 * dot(s, n) / dot(n, n)) * n;
			R.meet = 1;
			R.ut = 0;
		}
	}
	void refractData(const intersect &R, const ray &a, const double &mi, point &refract, double &rlr) const {
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

	bool contain(const point &A) const {
		return (C - A).mod() < r;
	}
	double SDF(const point &A) const {
		return (C - A).mod() - r;
	}

	void print(ostream& os) const {
		os << "Sphere(" << C << "," << r << ")";
	}
	unsigned telltype() const { return Sphere3D_Sign; }
};

/* Use to construct polyhedrons, shouldn't be directly added to World class */
#define Triangle_Ref_Sign 0x00010002
class triangle_ref : public object3D {
public:
	point A, B, C, N;	// N is the normal towards air, unit vector
	triangle_ref() {}
	triangle_ref(const triangle_ref& a) {
		A = a.A, B = a.B, C = a.C, N = a.N;
	}
	triangle_ref(const point &A, const point &B, const point &C) {
		this->A = A, this->B = B, this->C = C;
		N = cross(B - A, C - A); double m = N.mod();
		if (abs(m) < ERR_EPSILON) N = cross(B - A, B + C - 2 * A), m = N.mod();
		N /= N.mod();
	}
	~triangle_ref() {}
	object* copy() const {
		return new triangle_ref(*this);
	}
	point Max() const { return point(max({ A.x, B.x, C.x }), max({ A.y, B.y, C.y }), max({ A.z, B.z, C.z })); }
	point Min() const { return point(min({ A.x, B.x, C.x }), min({ A.y, B.y, C.y }), min({ A.z, B.z, C.z })); }

	void operator += (const point &P) {
		A += P, B += P, C += P;
	}

	void meet(intersect &R, const ray &a) const {
		R.meet = 0;
		point E1 = B - A, E2 = C - A, T, P = cross(a.dir, E2), Q;
		double det = dot(E1, P);
		if (abs(det) < ERR_EPSILON) return;
		T = a.orig - A;
		double t, u, v;
		u = dot(T, P) / det;
		if (u < 0.0 || u > 1.0) return;
		Q = cross(T, E1);
		v = dot(a.dir, Q) / det;
		if (v < 0.0 || u + v > 1.0) return;
		t = dot(E2, Q) / det;
		if (t < ERR_EPSILON) return;
		R.meet = 1;
		R.dist = t;
		R.intrs = (1 - u - v)*A + u * B + v * C;
		double ct = 2 * dot(a.dir, N);
		R.ut = ct < 0 ? 1 : 0;
		R.reflect = a.dir - ct * N;
		return;
	}
	void refractData(const intersect &R, const ray &a, const double &mi, point &refract, double &rlr) const {
		double c1 = -dot(N, a.dir);
		double n1 = mi, n2 = ri;
		if (R.ut == 0) n1 = ri, n2 = mi, c1 = -c1;
		double c = n1 / n2, c2 = sqrt(1 - c * c * (1 - c1 * c1));
		if (isnan(c2)) {
			refract.x = refract.y = refract.z = NAN; rlr = 1;
			return;
		}
		refract = c * a.dir - (R.ut == 0 ? (c*c1 - c2) : (c2 - c*c1))*N;
		double Rs = (n2*c1 - n1 * c2) / (n2*c1 + n1 * c2); Rs *= Rs;
		double Rp = (n2*c2 - n1 * c1) / (n2*c2 + n1 * c1); Rp *= Rp;
		rlr = 0.5*(Rs + Rp);
	}

	double SDF(const point &P) const {
		return dot(P - A, N);
	}
	bool contain(const point &P) const {
		return dot(P - A, N) < 0;
	}
	
	void print(ostream& os) const {
		os << "Surface(if(u+v<1," << noshowpos << A.x << "*(1-u-v)" << showpos << B.x << "*u" << showpos << C.x << "*v" << "), "
			<< noshowpos << A.y << "*(1-u-v)" << showpos << B.y << "*u" << showpos << C.y << "*v" << ", "
			<< noshowpos << A.z << "*(1-u-v)" << showpos << B.z << "*u" << showpos << C.z << "*v" << ", "
			<< "u, 0, 1, v, 0, 1)" << noshowpos;
	}
	unsigned telltype() const { return Triangle_Ref_Sign; }
};

#define Parallelogram_Ref_Sign 0x00010003
class parallelogram_ref : public object3D {
public:
	point O, A, B, N;	// N is the normal towards air, unit vector
	parallelogram_ref() {}
	parallelogram_ref(const parallelogram_ref& a) {
		A = a.A, B = a.B, O = a.O, N = a.N;
	}
	parallelogram_ref(const point &O, const point &A, const point &B) {
		this->O = O, this->A = A, this->B = B;
		N = cross(A, B); double m = N.mod();
		if (abs(m) < ERR_EPSILON) N = cross(A, A + B), m = N.mod();
		N /= N.mod();
	}
	parallelogram_ref(const point &O, const point &A, const point &B, bool absolute) {
		this->O = O, this->A = A, this->B = B;
		if (absolute) this->A -= O, this->B -= O;
		N = cross(this->A, this->B); double m = N.mod();
		if (abs(m) < ERR_EPSILON) N = cross(this->A, this->A + this->B), m = N.mod();
		N /= N.mod();
	}
	~parallelogram_ref() {}
	object* copy() const {
		return new parallelogram_ref(*this);
	}
	point Max() const {
		return point(max({ O.x, A.x, B.x, A.x + B.x - O.x }),
			max({ O.y, A.y, B.y, A.y + B.y - O.y }), max({ O.y, A.y, B.y,A.y + B.y - O.y }));
	}
	point Min() const {
		return point(min({ O.x, A.x, B.x, A.x + B.x - O.x }),
			min({ O.y, A.y, B.y, A.y + B.y - O.y }), min({ O.y, A.y, B.y,A.y + B.y - O.y }));
	}

	void operator += (const point &P) {
		A += P, B += P, O += P;
	}

	void meet(intersect &R, const ray &a) const {
		R.meet = 0;
		point T, P = cross(a.dir, B), Q;
		double det = dot(A, P);
		if (abs(det) < ERR_EPSILON) return;
		T = a.orig - O;
		double t, u, v;
		u = dot(T, P) / det;
		if (u < 0.0 || u > 1.0) return;
		Q = cross(T, A);
		v = dot(a.dir, Q) / det;
		if (v < 0.0 || v > 1.0) return;
		t = dot(B, Q) / det;
		if (t < ERR_EPSILON) return;
		R.meet = 1;
		R.dist = t;
		R.intrs = O + u * A + v * B;
		double ct = 2 * dot(a.dir, N);
		R.ut = ct < 0 ? 1 : 0;
		R.reflect = a.dir - ct * N;
		return;
	}
	void refractData(const intersect &R, const ray &a, const double &mi, point &refract, double &rlr) const {
		double c1 = -dot(N, a.dir);
		double n1 = mi, n2 = ri;
		if (R.ut == 0) n1 = ri, n2 = mi, c1 = -c1; // ut=1: air->obj; ut=0: obj->air
		double c = n1 / n2, c2 = sqrt(1 - c * c * (1 - c1 * c1));
		if (isnan(c2)) {
			refract.x = refract.y = refract.z = NAN; rlr = 1;
			return;
		}
		refract = c * a.dir - (R.ut == 0 ? (c*c1 - c2) : (c2 - c*c1))*N;
		double Rs = (n2*c1 - n1 * c2) / (n2*c1 + n1 * c2); Rs *= Rs;
		double Rp = (n2*c2 - n1 * c1) / (n2*c2 + n1 * c1); Rp *= Rp;
		rlr = 0.5*(Rs + Rp);
	}
	double SDF(const point &P) const {
		return dot(P - A, N);
	}
	bool contain(const point &P) const {
		return dot(P - A, N) < 0;
	}
	void print(ostream& os) const {
		os << "Surface(" << noshowpos << O.x << showpos << A.x << "*u" << showpos << B.x << "*v" << ", "
			<< noshowpos << O.y << showpos << A.y << "*u" << showpos << B.y << "*v" << ", "
			<< noshowpos << O.z << showpos << A.z << "*u" << showpos << B.z << "*v" << ", "
			<< "u, 0, 1, v, 0, 1)" << noshowpos;
	}
	unsigned telltype() const { return Parallelogram_Ref_Sign; }
};

#define Polyhedron_Sign 0x00010004
class polyhedron : public object3D {
	bool dynamic_memory;
public:
	vector<const object3D*> tp; // triangles or parallelograms
	borderbox S;
	polyhedron() : dynamic_memory(false) {}
	polyhedron(const polyhedron &a) : dynamic_memory(true) {
		for (int i = 0; i < a.tp.size(); i++) {
			tp.push_back(dynamic_cast<object3D*>(a.tp[i]->copy()));
		}
	}
	polyhedron(const initializer_list<object3D*> objs) : dynamic_memory(false) {
		for (unsigned i = 0, n = objs.size(); i < n; i++) {
			tp.push_back(objs.begin()[i]);
			S.Max = PMax(S.Max, tp.back()->Max());
			S.Min = PMin(S.Min, tp.back()->Min());
		}
	}
	void add(const object3D* obj) {
		tp.push_back(obj);
		S.Max = PMax(S.Max, tp.back()->Max());
		S.Min = PMin(S.Min, tp.back()->Min());
	}
	void add(const initializer_list<object3D*> objs) {
		for (unsigned i = 0, n = objs.size(); i < n; i++) {
			tp.push_back(objs.begin()[i]);
			S.Max = PMax(S.Max, tp.back()->Max());
			S.Min = PMin(S.Min, tp.back()->Min());
		}
	}
	object* copy() const {
		return new polyhedron(*this);
	}
	~polyhedron() {
		if (dynamic_memory) {
			for (int i = 0; i < tp.size(); i++) delete tp[i];
		}
		tp.clear();
	}

	// Before rendering, make sure all added triangles and parallelograms form a legal polyhedron, and the directions of normals are all correct

	point Max() const {
		return S.Max;

		point C(-INFINITY, -INFINITY, -INFINITY);
		point P;
		for (unsigned i = 0, n = tp.size(); i < n; i++) {
			P = tp[i]->Max();
			C.x = max(C.x, P.x), C.y = max(C.y, P.y), C.z = max(C.z, P.z);
		}
		return C;
	}
	point Min() const {
		return S.Min;

		point C(INFINITY, INFINITY, INFINITY);
		point P;
		for (unsigned i = 0, n = tp.size(); i < n; i++) {
			P = tp[i]->Min();
			C.x = min(C.x, P.x), C.y = min(C.y, P.y), C.z = min(C.z, P.z);
		}
		return C;
	}

	void meet(intersect &R, const ray &a) const {
		//cout << a << endl;
		R.meet = 0; if (!S.meet(a)) return;
		unsigned i = 0, n = tp.size();
		for (i = 0; i < n; i++) {
			tp[i]->meet(R, a);
			if (R.meet) {
				R.vt = i; break;
			}
		}
		if (!R.meet) return;
		intersect S;
		for (; i < n; i++) {
			tp[i]->meet(S, a);
			if (S.meet && S.dist < R.dist) R = S, R.vt = i;
		}
		return;
	}
	bool contain(const point &A) const {
		if (!S.contain(A)) return false;
		/*for (unsigned i = 0, n = tp.size(); i < n; i++) {
			if (tp[i]->contain(A)) return true;
		}
		return false;*/
		// This only work for convex polyhedrons

		unsigned N = 0;
		intersect R;
		for (unsigned i = 0, n = tp.size(); i < n; i++) {
			tp[i]->meet(R, ray(A, point(0.324837423289, 0.239844723344, 0.523423543445)));  // using irregular values to prevent special cases
			if (R.meet) N++;
		}
		return N & 1;
		// Work for all legal polyhedrons, but requires a lot of calculating
	}
	/*double SDF(const point &A) const {
		return NAN;
	}*/
	void refractData(const intersect &R, const ray &a, const double &mi, point &refract, double &rlr) const {
		tp[R.vt]->refractData(R, a, mi, refract, rlr);
	}

	unsigned telltype() const {
		return Polyhedron_Sign;
	}

	void print(ostream& os) const {
		for (int i = 0; i < tp.size(); i++) {
			os << "\t"; tp[i]->print(os); os << endl;
		}
	}
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
	object* copy() const {
		return new lightsource(*this);
	}

	void setcolor(const WebSafeColour &c) { col = color(c); }
	void setcolor(const pixel &c) { col = c; }
	void setcolor(const rgblight &c) { col = c; }
	void setcolor(const unsigned &hex) { col = pixel(hex); }

	// Note that in all derivative classes of lightsource, the meet function set "ut" as the cosine of angle of incidence

	void print(ostream& os) const { os << "lightsource"; }
};

#define SphereBulb_Sign 0x01000000
class spherebulb : public lightsource {
public:
	point C; double r;
	spherebulb() :r(0) { col.r = col.g = col.b = 1; }
	spherebulb(const point &C, const double &r) {
		this->C = C, this->r = r, col.r = col.g = col.b = 1;
	}
	spherebulb(const point &C, const double &r, const rgblight &col) {
		this->C = C, this->r = r, this->col = col;
	}
	spherebulb(const spherebulb &a) { C = a.C, r = a.r, col = a.col; }
	object* copy() const {
		return new spherebulb(*this);
	}
	~spherebulb() {}

	void meet(intersect &R, const ray &a) const {
		R.meet = 0;
		point p = C - a.orig;
		if (dot(p, a.dir) < 0) return;
		double d = cross(p, a.dir).mod();
		if (d > r) return;
		d *= d;
		R.dist = sqrt(dot(p, p) - d) - sqrt(r*r - d);
		point n = a.dir * R.dist, t = n - p;
		R.intrs = n + a.orig;
		R.reflect = n + abs(2 * dot(n, t) / dot(t, t)) * t;

		R.ut = abs(dot(t, a.dir)) / t.mod();
		R.meet = 1;
		return;
	}
	point Max() const {
		return point(C.x + r, C.y + r, C.z + r);
	}
	point Min() const {
		return point(C.x - r, C.y - r, C.z - r);
	}
	friend ostream& operator << (ostream& os, const spherebulb &a) {
		os << "Sphere(" << a.C << "," << a.r << ")";
		return os;
	}
	unsigned telltype() const { return SphereBulb_Sign; }

};


#define RectBulb_Sign 0x01000001
class rectbulb : public lightsource {
public:
	point O, A, B;
	rectbulb() {}
	rectbulb(const rectbulb& a) {
		col = a.col;
		A = a.A, B = a.B, O = a.O;
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
	object* copy() const {
		return new rectbulb(*this);
	}
	void meet(intersect &R, const ray &a) const {
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
		R.dist = t;
		R.intrs = (1 - u - v)*O + u * A + v * B;
		point OA = A - R.intrs, OB = B - R.intrs;
		point ON = cross(OA, OB);
		ON *= dot(a.dir, ON) / dot(ON, ON);		// NAN occurs when 0/0
		R.reflect = a.dir - 2 * ON;
		R.ut = dot(a.dir, ON) / ON.mod();
		return;
	}
	point Max() const {
		return point(max({ O.x, A.x, B.x, A.x + B.x - O.x }),
			max({ O.y, A.y, B.y, A.y + B.y - O.y }), max({ O.y, A.y, B.y,A.y + B.y - O.y }));
	}
	point Min() const {
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
	unsigned telltype() const { return RectBulb_Sign; }
};


#endif



