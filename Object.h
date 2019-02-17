#pragma once

#include "Matrix.h"
#include "D:\Explore\Math\Graph\GraphFun\GraphFun\BitMap.h"
#include <vector>
#include <initializer_list>


#include <chrono>
typedef chrono::high_resolution_clock NTime;
typedef chrono::duration<double> fsec;

class point {
	// This class can be used instead of 3D vector
public:
	double x, y, z;
	point() {
		x = y = z = 0;
	}
	point(double X, double Y, double Z) {
		x = X, y = Y, z = Z;
	}
	point(double X, double Y) {
		x = X, y = Y, z = 0;
	}
	point(const point &a) {
		x = a.x, y = a.y, z = a.z;
	}
	point(const initializer_list<double> &a) {
		x = *(a.begin()), y = *(a.begin() + 1), z = *(a.begin() + 2);
	}
	inline void operator = (const point &a) {
		x = a.x, y = a.y, z = a.z;
	}
	inline void operator = (const initializer_list<double> &a) {
		x = *(a.begin()), y = *(a.begin() + 1), z = *(a.begin() + 2);
	}
	~point() {}
	inline friend point Max(const point &A, const point &B) {
		return point(max(A.x, B.x), max(A.y, B.y), max(A.z, B.z));
	}
	inline friend point Min(const point &A, const point &B) {
		return point(min(A.x, B.x), min(A.y, B.y), min(A.z, B.z));
	}
	inline double mod() const { return sqrt(x * x + y * y + z * z); }
	inline friend point operator + (const point &a, const point &b) {
		return point(a.x + b.x, a.y + b.y, a.z + b.z);
	}
	inline void operator += (const point &a) {
		x += a.x, y += a.y, z += a.z;
	}
	inline friend point operator - (const point &a) {
		return point(-a.x, -a.y, -a.z);
	}
	inline friend point operator - (const point &a, const point &b) {
		return point(a.x - b.x, a.y - b.y, a.z - b.z);
	}
	inline void operator -= (const point &a) {
		x -= a.x, y -= a.y, z -= a.z;
	}
	inline friend point operator * (const double &a, const point &b) {
		return point(a*b.x, a*b.y, a*b.z);
	}
	inline friend point operator * (const point &a, const double &b) {
		return point(b*a.x, b*a.y, b*a.z);
	}
	inline void operator *= (const double &a) {
		x *= a, y *= a, z *= a;
	}
	inline friend point operator / (const point &a, const double &b) {
		return point(a.x / b, a.y / b, a.z / b);
	}
	inline void operator /= (const double &a) {
		x /= a, y /= a, z /= a;
	}
	inline friend double operator * (const point &a, const point &b) {
		return a.x*b.x + a.y*b.y + a.z*b.z;		// dot product
	}
	inline friend double dot(const point &a, const point &b) {
		return a.x*b.x + a.y*b.y + a.z*b.z;
	}
	inline friend point cross(const point &a, const point &b) {
		return point(a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x);
	}
	inline friend point operator * (const matrix<double> &M, const point &a) {
		return point(M[0][0] * a.x + M[0][1] * a.y + M[0][2] * a.z,
			M[1][0] * a.x + M[1][1] * a.y + M[1][2] * a.z,
			M[2][0] * a.x + M[2][1] * a.y + M[2][2] * a.z);
	}
	inline void operator *= (const matrix<double> &M) {
		*this = M * (*this);
	}
	friend ostream& operator << (ostream& os, const point &a) {
		os << "(" << a.x << "," << a.y << "," << a.z << ")";
		return os;
	}
};

class ray {
public:
	point orig, dir;
	ray() {}
	ray(const ray &a) {
		orig = a.orig, dir = a.dir;
	}
	ray(double x, double y, double z) {
		orig = { 0, 0, 0 }, dir = { x, y, z };
	}
	ray(point D) {
		orig = { 0,0,0 }, dir = D;
	}
	ray(double Ox, double Oy, double Oz, double Vx, double Vy, double Vz) {
		orig = { Ox, Oy, Oz }, dir = { Vx - Ox, Vy - Oy, Vz - Oz };
	}
	ray(point O, point D) {
		orig = O, dir = D;
	}
	~ray() {}
	friend ostream& operator << (ostream& os, const ray &r) {
		os << "(" << r.orig.x << "," << r.orig.y << "," << r.orig.z << "), d -> (" << r.dir.x << "," << r.dir.y << "," << r.dir.z << ")";
		return os;
	}
};

struct intersect {
	bool meet = 0;
	double dist;
	point intrs;
	point reflect;	// Vector
};
class rgblight {
public:
	double r, g, b;
	rgblight() :r(0), g(0), b(0) {}
	rgblight(const double &r, const double &g, const double &b) {
		this->r = r, this->g = g, this->b = b;
	}
	rgblight(const double &l) {
		r = g = b = l;
	}
	rgblight(const pixel &p) {
		r = p.r / 255.0, g = p.g / 255.0, b = p.b / 255.0;
	}
	inline void operator = (const initializer_list<double> &p) {
		this->r = p.begin()[0], this->g = p.begin()[1], this->b = p.begin()[2];
	}
	inline rgblight& operator = (const pixel &p) {
		r = p.r / 255.0, g = p.g / 255.0, b = p.b / 255.0;
		return *this;
	}
	inline rgblight operator ~() {
		return rgblight(1 - r, 1 - g, 1 - b);
	}
	~rgblight() {}
};
inline pixel rgb(const rgblight &a) {
	return drgb(a.r, a.g, a.b);
}

class object {
#define Object_Sign 0x00000000
public:
	rgblight reflect, refract, absorb;
	object() {}
	object(object* a) {
		reflect = a->reflect, refract = a->refract, absorb = a->absorb;
	}
	~object() {}
	virtual intersect meet(const ray &a) { cout << "\aobject::meet is called. This function should never be called. \n"; return intersect(); }
	virtual inline void operator *= (matrix<double> M) { cout << "\aobject::operator *= is called. This function shouldn't be called. \n"; }
	virtual inline void operator += (const point &a) { cout << "\aobject::operator += is called. This function shouldn't be called. \n"; }
	virtual inline point Max() { cout << "\aobject::Max is called. This function should never be called. \n"; return point(); }
	virtual inline point Min() { cout << "\aobject::Min is called. This function should never be called. \n"; return point(); }
	inline void fixcol(double &reflect, double &refract, double &absorb, double &sum) {
		sum = reflect + refract + absorb;
		if (sum > 0) {
			if (sum < 1 && absorb == 0) absorb = 1 - reflect - refract;
			else if (sum != 1) reflect /= sum, refract /= sum, absorb /= sum;
		}
		else reflect = refract = 0, absorb = 1;
	}
	inline void fix() {
		if (reflect.r < 0) reflect.r = 0;
		if (reflect.g < 0) reflect.g = 0;
		if (reflect.b < 0) reflect.b = 0;
		if (refract.r < 0) refract.r = 0;
		if (refract.g < 0) refract.g = 0;
		if (refract.b < 0) refract.b = 0;
		if (absorb.r < 0) absorb.r = 0;
		if (absorb.g < 0) absorb.g = 0;
		if (absorb.b < 0) absorb.b = 0;
		double sum;
		fixcol(reflect.r, refract.r, absorb.r, sum);
		fixcol(reflect.g, refract.g, absorb.g, sum);
		fixcol(reflect.b, refract.b, absorb.b, sum);

		/*sum = reflect.r + refract.r + absorb.r;
		if (sum != 0) reflect.r /= sum, refract.r /= sum, absorb.r /= sum;
		else reflect.r = refract.r = 0, absorb.r = 1;
		sum = reflect.g + refract.g + absorb.g;
		if (sum != 0) reflect.g /= sum, refract.g /= sum, absorb.g /= sum;
		else reflect.g = refract.g = 0, absorb.g = 1;
		sum = reflect.b + refract.b + absorb.b;
		if (sum != 0) reflect.b /= sum, refract.b /= sum, absorb.b /= sum;
		else reflect.b = refract.b = 0, absorb.b = 1;*/
	}
	virtual int telltype() { return Object_Sign; }
};
class triangle :public object {
#define Triangle_Sign 0x00000001
public:
	point A, B, C;
	triangle() {}
	triangle(triangle* a) {
		reflect = a->reflect, refract = a->refract, absorb = a->absorb;
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
	intersect meet(const ray &a) {
		// Algorithm: http://www.graphics.cornell.edu/pubs/1997/MT97.pdf
		intersect R;
		point E1 = B - A, E2 = C - A, T, P = cross(a.dir, E2), Q;
		double det = dot(E1, P);
		if (abs(det) < 1e-6) return R;
		T = a.orig - A;
		double t, u, v;
		// t: the multiple that is extended when the vector intersects the triangle; 
		// u,v: (1-u-v)·A + u·B + v·C  is where the vector intersects the triangle; 
		// " O + t·D = (1-u-v)·A + u·B + v·C "
		// This three variables may be useful. 
		u = dot(T, P) / det;
		if (u < 0.0 || u > 1.0) return R;
		Q = cross(T, E1);
		v = dot(a.dir, Q) / det;
		if (v < 0.0 || u + v > 1.0) return R;
		t = dot(E2, Q) / det;
		if (t < 1e-6) return R;
		R.meet = 1;
		R.dist = t * a.dir.mod();
		R.intrs = (1 - u - v)*A + u * B + v * C;
		point OB = B - R.intrs, OC = C - R.intrs;
		point ON = cross(OB, OC);
		ON *= dot(a.dir, ON) / dot(ON, ON);
		R.reflect = a.dir - 2 * ON;
		return R;
	}
	inline point Max() { return point(max({ A.x, B.x, C.x }), max({ A.y, B.y, C.y }), max({ A.z, B.z, C.z })); }
	inline point Min() { return point(min({ A.x, B.x, C.x }), min({ A.y, B.y, C.y }), min({ A.z, B.z, C.z })); }
	inline friend triangle operator * (matrix<double> M, triangle a) {
		a.A = M * a.A, a.B = M * a.B, a.C = M * a.C;
		return a;
	}
	inline void operator *= (matrix<double> M) {
		A = M * A, B = M * B, C = M * C;
	}
	inline friend triangle operator + (triangle t, const point &a) {
		t.A += a, t.B += a, t.C += a;
		return t;
	}
	inline void operator += (const point &a) {
		A += a, B += a, C += a;
	}
	friend ostream& operator << (ostream& os, const triangle &a) {
		os << "(" << a.A.x << "," << a.A.y << "," << a.A.z << "), ("
			<< a.B.x << "," << a.B.y << "," << a.B.z << "), (" << a.C.x << "," << a.C.y << "," << a.C.z << ")";
		return os;
	}
	int telltype() { return Triangle_Sign; }
};
class parallelogram :public object {
#define Parallelogram_Sign 0x00000002
public:
	point O, A, B;
	parallelogram() {}
	parallelogram(parallelogram* a) {
		reflect = a->reflect, refract = a->refract, absorb = a->absorb;
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
	intersect meet(const ray &a) {
		// Similar algorithm as that for triangle
		// I have proofed that P = (1-u-v)·O + uA + vB is on the plane formed by O,A,B, but I haven't proofed the condition when P is on triangle/parallelogram AOB. (Try to proof it. )
		intersect R;
		point E1 = A - O, E2 = B - O, T, P = cross(a.dir, E2), Q;
		double det = dot(E1, P);
		if (abs(det) < 1e-6) return R;
		T = a.orig - O;
		double t, u, v;
		u = dot(T, P) / det;
		if (u < 0.0 || u > 1.0) return R;
		Q = cross(T, E1);
		v = dot(a.dir, Q) / det;
		if (v < 0.0 || v > 1.0) return R;	// For triangle: v < 0.0 || u + v > 1.0
		t = dot(E2, Q) / det;
		if (t < 1e-6) return R;
		R.meet = 1;
		R.dist = t * a.dir.mod();
		R.intrs = (1 - u - v)*O + u * A + v * B;
		point OA = A - R.intrs, OB = B - R.intrs;
		point ON = cross(OA, OB);
		ON *= dot(a.dir, ON) / dot(ON, ON);
		R.reflect = a.dir - 2 * ON;
		return R;
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
	inline friend parallelogram operator + (parallelogram t, const point &a) {
		t.A += a, t.B += a, t.O += a;
		return t;
	}
	inline void operator += (const point &a) {
		A += a, B += a, O += a;
	}
	friend ostream& operator << (ostream& os, const parallelogram &a) {
		os << "(" << a.O.x << "," << a.O.y << "," << a.O.z << ") -> ("
			<< a.A.x << "," << a.A.y << "," << a.A.z << "), (" << a.B.x << "," << a.B.y << "," << a.B.z << "), end = ("
			<< (a.A.x + a.B.x - a.O.x) << "," << (a.A.y + a.B.y - a.O.y) << "," << (a.A.z + a.B.z - a.O.z) << ")";
		return os;
	}
	int telltype() { return Parallelogram_Sign; }
};
class circle :public object {
#define Circle_Sign 0x00000003
public:
	point C; double r; double rx, rz;	// First rotate x, than rotate z    R = Rz * Rx
	circle() :r(0), rx(0), rz(0) {}
	circle(circle* a) {
		reflect = a->reflect, refract = a->refract, absorb = a->absorb;
		C = a->C, r = a->r; rx = a->rx, rz = a->rz;
	}
	circle(const point &C, const double &r) {
		this->C = C, this->r = r, rx = rz = 0;
	}
	circle(const point &C, const double &r, const double &rx, const double &rz) {
		this->C = C, this->r = r, this->rx = rx, this->rz = rz;
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
	inline void operator *= (matrix<double> Orthogonal) {
		// Matrix must be the product of a rotation matrix with determinant 1 and a positive constant
		C *= Orthogonal;
		double d = cbrt(det(Orthogonal)); r *= d;
		//cout << (Orthogonal / d) << endl << ((Orthogonal / d) * Vector<double>(1, 1, 1)) << endl << endl;
		//rx += atan2(Orthogonal[2][1], Orthogonal[2][2]); rz += atan2(Orthogonal[1][0], Orthogonal[0][0]);
		rx -= acos(Orthogonal[2][2] / d); //rz += acos(Orthogonal[0][0] / d);
		//cout << rx << " " << rz << endl << endl;
		/*cout << (matrix<double>({ {cos(rz),-sin(rz),0}, {sin(rz),cos(rz),0}, {0,0,1} })
			* matrix<double>({ {1,0,0}, {0,cos(rx),-sin(rx)}, {0,sin(rx),cos(rx)} })) << endl <<
			(matrix<double>({ {cos(rz),-sin(rz),0}, {sin(rz),cos(rz),0}, {0,0,1} })
				* matrix<double>({ {1,0,0}, {0,cos(rx),-sin(rx)}, {0,sin(rx),cos(rx)} }) * Vector<double>(1, 1, 1)) << endl << endl;*/
				// For Euler's angle: https://stackoverflow.com/questions/15022630/how-to-calculate-the-angle-from-rotation-matrix
	}
	inline friend circle operator * (circle B, const double &a) {
		B.C *= a, B.r *= a; return B;
	}
	inline void operator *= (const double &a) {
		C *= a; r *= a;
	}
	inline point Max() {
		point R;
		R.z = -atan(cos(rx) * tan(rz));
		R.x = r * abs(cos(rz)*cos(R.z) - sin(rz)*cos(rx)*sin(R.z));
		R.z = atan(cos(rx) / tan(rz));
		R.y = r * abs(sin(rz)*cos(R.z) + cos(rx)*cos(rz)*sin(R.z));
		R.z = r * abs(sin(rx));
		R += C; return R;
	}
	inline point Min() {
		point R;
		R.z = -atan(cos(rx) * tan(rz));
		R.x = r * abs(cos(rz)*cos(R.z) - sin(rz)*cos(rx)*sin(R.z));
		R.z = atan(cos(rx) / tan(rz));
		R.y = r * abs(sin(rz)*cos(R.z) + cos(rx)*cos(rz)*sin(R.z));
		R.z = r * abs(sin(rx));
		return C - R;
	}
	intersect meet(const ray &a) {
		//cout << a << endl;
		ray ry(a.orig - C, a.dir);
		double sx = -sin(rz), cy = cos(rz);
		ry.orig = point(cy*ry.orig.x - sx * ry.orig.y, sx * ry.orig.x + cy * ry.orig.y, ry.orig.z),
			ry.dir = point(cy*ry.dir.x - sx * ry.dir.y, sx * ry.dir.x + cy * ry.dir.y, ry.dir.z);
		sx = -sin(rx), cy = cos(rx);
		ry.orig = point(ry.orig.x, cy*ry.orig.y - sx * ry.orig.z, sx * ry.orig.y + cy * ry.orig.z),
			ry.dir = point(ry.dir.x, cy*ry.dir.y - sx * ry.dir.z, sx * ry.dir.y + cy * ry.dir.z);
		//cout << ry << endl;
		intersect R;
		if ((ry.orig.z < -1e-6) == (ry.dir.z < -1e-6)) return R;
		if ((ry.orig.z > 1e-6) == (ry.dir.z > 1e-6)) return R;
		double t = -ry.orig.z / ry.dir.z;
		R.intrs.x = ry.dir.x*t + ry.orig.x, R.intrs.y = ry.dir.y*t + ry.orig.y;
		R.dist = R.intrs.x * R.intrs.x + R.intrs.y * R.intrs.y;
		if (R.dist > r*r) return R;
		//cout << R.intrs << " " << R.dist << endl;
		R.dist = (R.intrs - ry.orig).mod();
		R.reflect = point(2 * R.intrs.x - ry.orig.x, 2 * R.intrs.y - ry.orig.y, ry.orig.z);
		//cout << R.reflect << endl;
		R.reflect = point(R.reflect.x, cy*R.reflect.y + sx * R.reflect.z, -sx*R.reflect.y + cy * R.reflect.z),
			R.intrs = point(R.intrs.x, cy*R.intrs.y + sx * R.intrs.z, -sx*R.intrs.y + cy * R.intrs.z);
		sx = -sin(rz), cy = cos(rz);
		R.reflect = point(cy*R.reflect.x + sx * R.reflect.y, -sx*R.reflect.x + cy * R.reflect.y, R.reflect.z),
			R.intrs = point(cy*R.intrs.x + sx * R.intrs.y, -sx*R.intrs.x + cy * R.intrs.y, R.intrs.z);
		R.intrs += C;
		//cout << R.intrs << " " << R.reflect << endl;
		R.meet = 1;
		//cout << endl;
		return R;
	}
	int telltype() { return Circle_Sign; }
};
class sphere :public object {
#define Sphere_Sign 0x00010001
public:
	point C; double r;
	sphere() :r(0) {}
	sphere(const point &C, const double &r) {
		this->C = C, this->r = r;
	}
	sphere(sphere* a) {
		reflect = a->reflect, refract = a->refract, absorb = a->absorb;
		C = a->C, r = a->r;
	}
	sphere(const initializer_list<double> &C, const double &r) {
		this->C = point(C), this->r = r;
	}
	intersect meet(const ray &a) {
		intersect R;
		point p = C - a.orig;
		if (dot(p, a.dir) < 0) return R;
		double sm = a.dir.mod();
		double d = cross(p, a.dir).mod() / sm;
		if (d - r > 1e-6) return R;
		d *= d;
		R.dist = sqrt(dot(p, p) - d) - sqrt(r*r - d);
		point n = a.dir * (R.dist / sm), t = n - p;
		R.intrs = n + a.orig;
		R.reflect = n + abs(2 * dot(n, t) / dot(t, t)) * t;
		R.meet = 1;
		return R;
	}
	inline friend sphere operator * (matrix<double> Orthogonal, sphere a) {
		a *= Orthogonal; return a;
	}
	inline void operator *= (matrix<double> Orthogonal) {
		// Matrix must be the product of an orthogonal matrix and a constant
		C *= Orthogonal; r *= cbrt(det(Orthogonal));
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
	inline void operator *= (const double &a) {
		C *= a; r *= a;
	}


	int telltype() { return Sphere_Sign; }
};

//#define Ray_Tracing_Debug
#define TestX 54
#define TestY 42

#include <thread>
#include <mutex>
#include <Windows.h>

class World {
	vector<object*> Objs;
	vector<World*> GObjs;
	class parallelepiped {
		// Use for border of subgroups of objects
	public:
		point O, A, B, C;
		parallelepiped() {}
		parallelepiped(const point &O, const point &A, const point &B, const point &C) {
			this->O = O, this->A = O + A, this->B = O + B, this->C = O + C;
		}
		parallelepiped(const point &O, const point &A, const point &B, const point &C, bool absolute) {
			this->O = O, this->A = A, this->B = B, this->C = C;
			if (!absolute) this->A += O, this->B += O, this->C += O;
		}
		parallelepiped(const initializer_list<double> &O, const initializer_list<double> &A, const initializer_list<double> &B, const initializer_list<double> &C) {
			this->O = O, this->A = A, this->B = B, this->C = C;
			this->A += O, this->B += O, this->C += O;
		}
		void operator = (const initializer_list<initializer_list<double>> &V) {
			O = *V.begin(), A = *(V.begin() + 1), B = *(V.begin() + 2), C = *(V.begin() + 3);
		}
		inline point Max() {
			/*return point(max({ O.x, O.x + A.x, O.x + B.x, O.x + C.x, O.x + A.x + B.x, O.x + A.x + C.x, O.x + B.x + C.x, O.x + A.x + B.x + C.x }),
				max({ O.y, O.y + A.y, O.y + B.y, O.y + C.y, O.y + A.y + B.y, O.y + A.y + C.y, O.y + B.y + C.y, O.y + A.y + B.y + C.y }),
				max({ O.z, O.z + A.z, O.z + B.z, O.z + C.z, O.z + A.z + B.z, O.z + A.z + C.z, O.z + B.z + C.z, O.z + A.z + B.z + C.z }));*/
			return point(max(O.x, O.x + A.x), max(O.y, O.y + B.y), max(O.z, O.z + C.z));	// This class is only used as a cuboid. 
		}
		inline point Min() {
			/*return point(min({ O.x, O.x + A.x, O.x + B.x, O.x + C.x, O.x + A.x + B.x, O.x + A.x + C.x, O.x + B.x + C.x, O.x + A.x + B.x + C.x }),
				min({ O.y, O.y + A.y, O.y + B.y, O.y + C.y, O.y + A.y + B.y, O.y + A.y + C.y, O.y + B.y + C.y, O.y + A.y + B.y + C.y }),
				min({ O.z, O.z + A.z, O.z + B.z, O.z + C.z, O.z + A.z + B.z, O.z + A.z + C.z, O.z + B.z + C.z, O.z + A.z + B.z + C.z }));*/
			return point(min(O.x, O.x + A.x), min(O.y, O.y + B.y), min(O.z, O.z + C.z));
		}
		intersect meet(ray a) {
			intersect R;
			R = parallelogram(O + C, A, B).meet(a);		// Top
			if (R.meet) return R;
			R = parallelogram(O, C, B).meet(a);		// Front
			if (R.meet) return R;
			R = parallelogram(O + A, C, B).meet(a);		// Right
			if (R.meet) return R;
			R = parallelogram(O, C, A).meet(a);		// Left
			if (R.meet) return R;
			R = parallelogram(O, A, B).meet(a);		// Buttom
			if (R.meet) return R;
			R = parallelogram(O + B, A, C).meet(a);		// Back
			// necessary to include all six faces (since sometimes ray start inside the parallelepiped)
			return R;
		}
		~parallelepiped() {}
		friend ostream& operator << (ostream& os, const parallelepiped &a) {
			os << a.O << " -> " << "{" << a.A.x << ", " << a.B.y << ", " << a.C.z << "}";
			return os;
		}
	};
	class WIntersect {
	public:
		vector<intersect*> Intrs;
		vector<WIntersect*> GIntrs;
		WIntersect() {}
		~WIntersect() {
			while (!Intrs.empty()) {
				delete Intrs.back();
				Intrs.pop_back();
			}
			while (!GIntrs.empty()) {
				GIntrs.back()->~WIntersect(); delete GIntrs.back();
				GIntrs.pop_back();
			}
		}
		void resize(const World* p) {
			if (Intrs.size() != p->Objs.size()) {
				while (!Intrs.empty()) {
					delete Intrs.back();
					Intrs.pop_back();
				}
				Intrs.resize(p->Objs.size());
				for (int i = Intrs.size() - 1; i >= 0; i--) {
					Intrs.at(i) = new intersect;
				}
			}
			while (!GIntrs.empty()) {
				delete GIntrs.back();
				GIntrs.pop_back();
			}
			for (int i = 0; i < p->GObjs.size(); i++) {
				GIntrs.push_back(new WIntersect);
				GIntrs.back()->resize(p->GObjs.at(i));
			}
		}
	};
	WIntersect WIntrs;
	parallelepiped border;
	int RenderingProcess;
	void resize() {
		point maxc, minc;
		if (!Objs.empty()) maxc = Objs.front()->Max(), minc = Objs.front()->Min();
		for (int i = 1; i < Objs.size(); i++) {
			maxc = Max(maxc, Objs.at(i)->Max());		// Don't use "if (xxx > xxx) ...."  that would take more time
			minc = Min(minc, Objs.at(i)->Min());
		}
		for (int i = 0; i < GObjs.size(); i++) GObjs.at(i)->resize();
		if (Objs.size() == 0 && !GObjs.empty()) {
			maxc = GObjs.front()->border.Max(), minc = GObjs.front()->border.Min();
		}
		for (int i = ((Objs.size() == 0 && !GObjs.empty()) ? 1 : 0); i < GObjs.size(); i++) {
			maxc = Max(maxc, GObjs.at(i)->border.Max());
			minc = Min(minc, GObjs.at(i)->border.Min());
		}
		maxc -= minc;
		border.O = minc;
		border.A = point(maxc.x, 0, 0), border.B = point(0, maxc.y, 0), border.C = point(0, 0, maxc.z);
		WIntrs.resize(this);
	}
	friend int main();
public:
	rgblight background;
	World() {}
	World(const World &W) {
		Objs.resize(W.Objs.size());
		for (int i = 0; i < Objs.size(); i++) {
			switch (W.Objs.at(i)->telltype()) {
			case Triangle_Sign: {
				Objs.at(i) = new triangle((triangle*)(W.Objs.at(i)));
				break;
			}
			case Parallelogram_Sign: {
				Objs.at(i) = new parallelogram((parallelogram*)(W.Objs.at(i)));
				break;
			}
			case Circle_Sign: {
				Objs.at(i) = new circle((circle*)(W.Objs.at(i)));
				break;
			}
			case Sphere_Sign: {
				Objs.at(i) = new sphere((sphere*)(W.Objs.at(i)));
				break;
			}
			default: {
				Objs.at(i) = new object(W.Objs.at(i));
			}
			}
		}
		GObjs.resize(W.GObjs.size());
		for (int i = 0; i < GObjs.size(); i++) {
			GObjs.at(i) = new World(*W.GObjs.at(i));
		}
		this->resize();
		RenderingProcess = 0;
		background = W.background;
		return;
	}
	~World() {
		Objs.clear(); 	// DO NOT delete ANY elements there !!!
		GObjs.clear();
	}
	void destruct() {
		// Manually call this function if necessary to delete pointers. 
		for (int i = 0; i < Objs.size(); i++) delete Objs.at(i);
		Objs.clear();
		for (int i = 0; i < GObjs.size(); i++) GObjs.at(i)->destruct(), delete GObjs.at(i);
		GObjs.clear();
		WIntrs.~WIntersect();
	}
	inline void insert(World *a) {
		GObjs.push_back(a);
		this->resize();
	}
	inline void add(object* a) {
		Objs.push_back(a);
	}
	inline void add(object* a, point insert) {
		Objs.push_back(a);
		*Objs.back() += insert;
	}
	inline void add(initializer_list<object*> a) {
		for (int i = 0; i < a.size(); i++) Objs.push_back(a.begin()[i]);
	}
	void operator *= (matrix<double> T) {
		for (int i = 0; i < Objs.size(); i++) {
			*Objs.at(i) *= T;
		}
		for (int i = 0; i < GObjs.size(); i++) {
			*GObjs.at(i) *= T;
		}
		this->resize();
	}
	void operator += (point V) {
		for (int i = 0; i < Objs.size(); i++) {
			*Objs.at(i) += V;
		}
		for (int i = 0; i < GObjs.size(); i++) {
			*GObjs.at(i) += V;
		}
		border.O += V;
	}
	void RayTracing_EnumObjs(const ray &v, intersect* &ni, object* &no, const WIntersect* Ip) const {
		int NB = Objs.size();
		for (int k = 0; k < NB; k++) {
			(*Ip->Intrs.at(k)) = Objs.at(k)->meet(v);
		}	// Check meeting of single objects
		for (int k = 0; k < NB; k++) {
			if ((*Ip->Intrs.at(k)).meet) {
				if (ni == 0 || no == 0) ni = Ip->Intrs.at(k), no = Objs.at(k);
				else if (ni->dist > (*Ip->Intrs.at(k)).dist) ni = Ip->Intrs.at(k), no = Objs.at(k);
			}
		}	// Find the nearest object among single objects

		NB = GObjs.size();
		for (int i = 0; i < NB; i++) {
			if (GObjs.at(i)->border.meet(v).meet) {
				GObjs.at(i)->RayTracing_EnumObjs(v, ni, no, Ip->GIntrs.at(i));
			}
		}
	}
	rgblight RayTracing(const ray &v, int &count, const point &N) const {
		if (count > 100) return rgblight(1, 1, 1);
		int NB = Objs.size();
		//int M = -1, n = -1;
		intersect* ni = 0; object* no = 0;	// keeps the nearest object and it's reflect
		const World* Wp = this;	// World of the nearest object
		const WIntersect* Ip = &WIntrs;	// Intersect of the nearest object

		RayTracing_EnumObjs(v, ni, no, Ip);

		if (ni == 0 || no == 0) {
			//return rgblight(0, 0, 0);
			double ang = dot(N, v.dir);
			if (ang <= 0) return background;
			ang /= N.mod()*v.dir.mod();
			return rgblight(ang);
		}
		else {
			//return rgblight(no->reflect.r, no->reflect.g, no->reflect.b);
			rgblight c;
			c = RayTracing(ray(ni->intrs, ni->reflect), ++count, N);
			c.r *= no->reflect.r, c.g *= no->reflect.g, c.b *= no->reflect.b;
			return c;

		}
	}
	void RenderingProcessCounter(int ps, int &T1, int &T2, int &T3, int &T4) {
		int m, n = 0;
		int sum;
		while (1) {
			sum = T1 + T2 + T3 + T4;
			m = sum * 1000 / ps;
			if (m != n) cout << "\r" << (m / 10) << "." << (m % 10) << "%";
			n = m;
			//if (n == 1000) break;
			if (sum >= ps) break;
			Sleep(50);
		}
	}
	void MultiThread_CC(bitmap &canvas, int begw, int endw, int begh, int endh, ray beg, int cx, int cy, point N) {
		rgblight c;
		int n;
		//*➤*/ auto t0 = NTime::now(); auto t1 = NTime::now(); fsec fs = t1 - t0; pixel* p; double a;
		for (int i = begw; i < endw; i++) {
			beg.dir.x = i - cx;
			for (int j = begh; j < endh; j++) {
				beg.dir.y = j - cy;
				n = 0;
				//*➤*/ t0 = NTime::now();
				c = RayTracing(beg, n, N);
				//*➤*/ t1 = NTime::now(); fs = t1 - t0;
				c.r = sqrt(1 - (c.r - 1)*(c.r - 1)), c.g = sqrt(1 - (c.g - 1)*(c.g - 1)), c.b = sqrt(1 - (c.b - 1)*(c.b - 1));	// Make it lighter
				canvas.dot(i, j, drgb(c.r, c.g, c.b));
				//*➤*/ p = canvas[j] + i; p->r = p->b = 0, p->g /= 2; a = (tanh(log2(fs.count() * 1000)) + 1) / 2; canvas.dot(i, j, drgb(a, 0, 0), 1 - a);
				RenderingProcess++;
			}
		}
		this->destruct();
	}
	void render(bitmap &canvas, point W, double lrx, double lrz, double rx, double ry, double rz, double cx, double cy, double tms) {
		// Algorithm: https://de.wikipedia.org/wiki/Raytracing (Language: German)

		cout << "Initializing...";

		canvas.clear();
		vector<object*> obj = Objs;
		const int NB = obj.size(); int n;
		vector<intersect> P; P.resize(NB);
		matrix<double> T = matrix<double>({ {1,0,0}, {0,cos(rx),-sin(rx)}, {0,sin(rx),cos(rx)} })
			* matrix<double>({ {cos(ry),0,sin(ry)}, {0,1,0}, {-sin(ry),0,cos(ry)} })
			* matrix<double>({ {cos(rz),-sin(rz),0}, {sin(rz),cos(rz),0}, {0,0,1} }) * tms;
		// First rotate around z-axis, then whole object rotate around original y, then x. 
		// A intuitive matrix viewer: https://www.geogebra.org/m/cxqzveav 
		//cout << endl << (T / tms) << endl;
		/*matrix<double> T = matrix<double>({ {cos(rz), sin(rz), 0}, {-sin(rz), cos(rz), 0}, {0, 0, 1} }) *
			matrix<double>({ {1, 0, 0}, {0, cos(ry), sin(ry)}, {0, -sin(ry), cos(ry)} }) *
			matrix<double>({ { cos(rx), sin(rx), 0 }, { -sin(rx), cos(rx), 0 }, { 0, 0, 1 } }) * tms;	// Euler's angle */
		point add(cx, cy, 0);
		for (int i = 0; i < obj.size(); i++) *obj.at(i) *= T, *obj.at(i) += add;
		World* Wp;
		for (int i = 0; i < GObjs.size(); i++) {
			Wp = GObjs.at(i);
			*Wp *= T, *Wp += add;
		}

		ray beg(cx, cy, W.mod()*tms, 0, 0, 0);
		for (int i = 0; i < NB; i++) obj.at(i)->fix();
		point N(cos(lrx), sqrt(1 - cos(lrx)*cos(lrx) - cos(lrz)*cos(lrz)), cos(lrz)); N = T * N;
		rgblight c;
		int percent;

		cout << "\nCloning..." << endl;

		this->resize();
		//cout << *((parallelogram*)(Objs.front())) << endl;
		World This0(*this), This1(*this), This2(*this), This3(*this);

		cout << "\r               \r";
		cout << "Attempting...";

		vector<double> attempt; attempt.resize(canvas.width() / 10);
		auto t0 = NTime::now();
		auto t1 = NTime::now();
		fsec fs = t1 - t0;
		for (int i = 0; i < canvas.width(); i += 10) {
			t0 = NTime::now();
			beg.dir.x = i - cx;
			for (int j = 0; j < canvas.height(); j += 10) {
				beg.dir.y = j - cy;
				n = 0;
				c = RayTracing(beg, n, N);
			}
			t1 = NTime::now();
			fs = t1 - t0;
			attempt.at(i / 10) = fs.count();
		}
		double sum = 0, sumt = 0;
		int B1 = -1, B2 = -1, B3 = -1;
		for (int i = 0; i < attempt.size(); i++) {
			sum += attempt.at(i);
		}
		sum /= 4;
		for (int i = 0; i < attempt.size(); i++) {
			sumt += attempt.at(i);
			if (sumt > sum) {
				if (B1 == -1) B1 = i;
				else if (B2 == -1) B2 = i;
				else if (B3 == -1) B3 = i;
				else break;
				sumt = 0;
			}
		}
		B1 *= 10, B2 *= 10, B3 *= 10;
		if (B3 == -10) B3 = canvas.width();
		if (B2 == -10) B2 = canvas.width();
		if (B1 == -10) B1 = canvas.width();
		sum = floor(sum * 120);

		cout << "\r             \r"; if (sum >= 1) cout << "Estimated Time Required: " << int(sum) << "secs. \n\n";
		cout << "Rendering...\n";

		thread T0([&](World* WC) { WC->MultiThread_CC(canvas, 0, B1, 0, canvas.height(), beg, cx, cy, N); }, &This0);
		thread T1([&](World* WC) { WC->MultiThread_CC(canvas, B1, B2, 0, canvas.height(), beg, cx, cy, N); }, &This1);
		thread T2([&](World* WC) { WC->MultiThread_CC(canvas, B2, B3, 0, canvas.height(), beg, cx, cy, N); }, &This2);
		thread T3([&](World* WC) { WC->MultiThread_CC(canvas, B3, canvas.width(), 0, canvas.height(), beg, cx, cy, N); }, &This3);
		thread Proc([&](World* WC) { WC->RenderingProcessCounter(canvas.height()*canvas.width(),
			This0.RenderingProcess, This1.RenderingProcess, This2.RenderingProcess, This3.RenderingProcess); }, this);
		T0.join(); T1.join(); T2.join(); T3.join(); Proc.join();

		//int debugx = 360, debugy = 190;
		//This0.MultiThread_CC(canvas, debugx, debugx + 1, debugy, debugy + 1, beg, cx, cy, N);
		//canvas.dot(debugx + 1, debugy, rgb(255, 255, 0)); canvas.dot(debugx - 1, debugy, rgb(255, 255, 0)); canvas.dot(debugx, debugy + 1, rgb(255, 255, 0)); canvas.dot(debugx, debugy - 1, rgb(255, 255, 0));

		cout << " Completed. \n\n";

		obj.clear();

	}
};
