#pragma once

/*
	This header file defines all classes, structs, functions "Object.h" needs.
	For convinience in debugging, all console outputs are GeoGebra commands. (This also for "Object.h")
	GeoGebra: https://www.geogebra.org/3d
*/

#include "Matrix.h"
#include "D:\Explore\Math\Graph\GraphFun\GraphFun\BitMap.h"

using namespace std;

#ifndef PI
#define PI 3.1415926535897932384626433832795029L
#endif

/* Time Recorder */
#include <chrono>
typedef chrono::high_resolution_clock NTime;
typedef chrono::duration<double> fsec;

#include <cstdlib>

#define ERR_EPSILON 1e-8	// minimum distance of intersection
#define ERR_EPSILON_SN 1e-9
#define ERR_UPSILON 1e+12
#define ERR_ZETA 1e-6	// minumum distance of SDF test, may cause threads jointing problems when too large
#define ERR_ZETA_SN 1e-6

// For Debugging
extern ofstream fout("IMAGE\\Log.txt");
void WARN(string s) {
	cout << s << "\a\n"; fout << s << endl;
}



/*
	Spacial point class.
	Defined addition, substraction, dot and cross product, etc. Can be used as a spacial vector.
*/
class point {
public:
	double x, y, z;
	point() {
		x = y = z = 0;
	}
	point(const double &X, const double &Y, const double &Z) {
		x = X, y = Y, z = Z;
	}
	point(const double &X, const double &Y) {
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
		os << noshowpos << "(" << a.x << "," << a.y << "," << a.z << ")";
		return os;
	}
};
class point2D {
public:
	double x, y;
	point2D() { x = 0, y = 0; }
	point2D(const double &x, const double &y) { this->x = x, this->y = y; }
	point2D(const point2D &another) { x = another.x, y = another.y; }
	inline void operator = (const point2D &another) { x = another.x, y = another.y; }
	~point2D() {}

	inline point2D operator + (const point2D &a) const {
		return point2D(this->x + a.x, this->y + a.y);
	}
	inline point2D operator - () const {
		return point2D(-this->x, -this->y);
	}
	inline point2D operator - (const point2D &a) const {
		return point2D(this->x - a.x, this->y - a.y);
	}
	inline void operator += (const point2D &a) {
		this->x += a.x, this->y += a.y;
	}
	inline void operator -= (const point2D &a) {
		this->x -= a.x, this->y -= a.y;
	}
	friend inline point2D operator * (const double &k, const point2D &P) {
		return point2D(k*P.x, k*P.y);
	}
	inline point2D operator * (const double &k) const {
		return point2D(this->x*k, this->y*k);
	}
	inline point2D operator / (const double &k) const {
		return point2D(this->x / k, this->y / k);
	}
	inline void operator *= (const double &k) {
		this->x *= k, this->y *= k;
	}
	inline void operator /= (const double &k) {
		this->x /= k, this->y /= k;
	}
	inline double mod() const {
		return hypot(x, y);
	}
	inline friend double dot(const point2D &a, const point2D &b) {
		return a.x*b.x + a.y*b.y;
	}
	inline friend double cross(const point2D &a, const point2D &b) {
		return a.x*b.y - a.y*b.x;
	}
	void rotate(const double &r) {
		*this = point2D(cos(r)*x - sin(r)*y, sin(r)*x + cos(r)*y);
	}
	friend ostream& operator << (ostream& os, const point2D &a) {
		os << "(" << a.x << "," << a.y << ")";
		return os;
	}
};

/* Ray class */
class ray {
public:
	point orig, dir;	// P = O + t·D
	ray() {}
	ray(const ray &a) {
		orig = a.orig, dir = a.dir;
	}
	ray(double x, double y, double z) {
		orig = { 0, 0, 0 }, dir = { x, y, z };
	}
	ray(const point &D) {
		orig = { 0,0,0 }, dir = D;
	}
	ray(double Ox, double Oy, double Oz, double Vx, double Vy, double Vz) {
		orig = { Ox, Oy, Oz }, dir = { Vx - Ox, Vy - Oy, Vz - Oz };
	}
	ray(const point &O, const point &D) {
		orig = O, dir = D;
	}
	~ray() {}
	friend ostream& operator << (ostream& os, const ray &r) {
		os << "Vector(" << r.orig << "," << (r.orig + r.dir) << ")";
		return os;
	}
};
class ray2D {
public:
	point2D orig, dir;
	ray2D() {}
	ray2D(const ray2D &a) { this->orig = a.orig, this->dir = a.dir; }
	ray2D(const double &x, const double &y) { dir.x = x, dir.y = y; }
	ray2D(const point2D &O, const point2D &D) { orig = O, dir = D; }
	~ray2D() {}

	friend ostream& operator << (ostream& os, const ray2D &r) {
		os << "Vector(" << r.orig << "," << (r.orig + r.dir) << ")";
		return os;
	}
};

/* Struct of intersection data */
struct intersect {
	bool meet = 0;	// intersect?
	double dist;	// distance from the intersect to the origin
	point intrs;	// point of intersection
	point reflect;	// reflect vector
	double ut, vt;	// reserve
};
struct intersect2D {
	bool meet = 0;
	double dist;
	point2D intrs;
	point2D reflect;
	double ut;
};

/* Extreme points of an object */
struct maxmin {
	point max;	// maximum x,y,z coordinates
	point min;	// minimum x,y,z coordinates
};

/* Light class, each channel is a float between 0 and 1 */
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
	rgblight(const WebSafeColour &c) {
		unsigned n = c;
		b = (n & 0xFF) / 256.0, n >>= 8;
		g = (n & 0xFF) / 256.0, n >>= 8;
		r = (n & 0xFF) / 256.0;
	}
	inline void operator = (const initializer_list<double> &p) {
		this->r = p.begin()[0], this->g = p.begin()[1], this->b = p.begin()[2];
	}
	inline rgblight& operator = (const pixel &p) {
		r = p.r / 256.0, g = p.g / 256.0, b = p.b / 256.0;
		return *this;
	}
	inline rgblight operator ~() {
		return rgblight(1 - r, 1 - g, 1 - b);
	}
	inline friend rgblight operator + (const rgblight &a, const rgblight &b) {
		return rgblight(a.r + b.r, a.g + b.g, a.b + b.b);
	}
	inline friend void operator += (rgblight &a, const rgblight &b) {
		a.r += b.r, a.g += b.g, a.b += b.b;
	}
	inline friend rgblight operator * (const rgblight &a, const double &b) {
		return rgblight(a.r*b, a.g*b, a.b*b);
	}
	inline friend rgblight operator * (const double &a, const rgblight &b) {
		return rgblight(b.r*a, b.g*a, b.b*a);
	}
	inline friend rgblight operator * (const rgblight &a, const rgblight &b) {
		return rgblight(a.r*b.r, a.g*b.g, a.b*b.b);
	}
	inline friend rgblight operator / (const rgblight &a, const double &b) {
		return rgblight(a.r / b, a.g / b, a.b / b);
	}
	inline double vsl() {
		return 0.299*r + 0.587*g + 0.114*b;
	}
	~rgblight() {}
};
/* operations of rgblight, for debugging */
inline pixel operator + (const pixel &a, const pixel &b) {
	return pixel(a.r + b.r, a.g + b.g, a.b + b.b);
}
inline pixel operator * (const pixel &a, const double &b) {
	return pixel(a.r*b, a.g*b, a.b*b);
}
inline pixel operator / (const pixel &a, const double &b) {
	return pixel(a.r / b, a.g / b, a.b / b);
}
inline pixel rgb(const rgblight &a) {
	return drgb(a.r >= 1 ? 0.9999 : a.r, a.g >= 1 ? 0.9999 : a.g, a.b >= 1 ? 0.9999 : a.b);
}
inline void operator *= (rgblight &a, const double &b) {
	a.r *= b, a.g *= b, a.b *= b;
}
inline void operator /= (rgblight &a, const double &b) {
	a.r /= b, a.g /= b, a.b /= b;
}


class borderbox {
	// Use for border of subgroups of objects
public:
	point B1, B2;
	inline point& Max() { return B2; }
	inline point& Min() { return B1; }
	borderbox() { B1 = point(INFINITY, INFINITY, INFINITY), B2 = point(-INFINITY, -INFINITY, -INFINITY); }
	bool meet(const ray &a) const {
		// http://www.cs.utah.edu/~awilliam/box/box.pdf
		double tmin, tmax, tymin, tymax, tzmin, tzmax;
		tmin = ((a.dir.x < 0 ? B2 : B1).x - a.orig.x) / a.dir.x;
		tmax = ((a.dir.x < 0 ? B1 : B2).x - a.orig.x) / a.dir.x;
		tymin = ((a.dir.y < 0 ? B2 : B1).y - a.orig.y) / a.dir.y;
		tymax = ((a.dir.y < 0 ? B1 : B2).y - a.orig.y) / a.dir.y;
		if ((tmin > tymax) || (tymin > tmax)) return 0;
		if (tymin > tmin) tmin = tymin;
		if (tymax < tmax) tmax = tymax;
		tzmin = ((a.dir.z < 0 ? B2 : B1).z - a.orig.z) / a.dir.z;
		tzmax = ((a.dir.z < 0 ? B1 : B2).z - a.orig.z) / a.dir.z;
		if ((tmin > tzmax) || (tzmin > tmax)) return 0;
		if (tzmin > tmin) tmin = tzmin;
		if (tzmax < tmax) tmax = tzmax;
		return tmax > 0;
	}
	inline bool contain(const point &a) const {
		return (a.x > B1.x && a.x<B2.x && a.y>B1.y && a.y<B2.y && a.z>B1.z && a.z < B2.z);
	}
	~borderbox() {}
	friend ostream& operator << (ostream& os, const borderbox &a) {
		point A = a.B1, B = point(a.B2.x, a.B1.y, a.B1.z), C = point(a.B2.x, a.B2.y, a.B1.z), D = point(a.B1.x, a.B2.y, a.B1.z),
			E = point(a.B1.x, a.B1.y, a.B2.z), F = point(a.B2.x, a.B1.y, a.B2.z), G = a.B2, H = point(a.B1.x, a.B2.y, a.B2.z);
		os << "Polyline(" << A << "," << B << "," << F << "," << G << "," << C << "," << D << "," << H << "," << E << "," <<
			A << "," << E << "," << F << "," << B << "," << C << "," << G << "," << H << "," << D << "," << A << ")";
		return os;
	}
};

inline point PMax(point A, point B) {
	return point(max(A.x, B.x), max(A.y, B.y), max(A.z, B.z));
}
inline point PMin(point A, point B) {
	return point(min(A.x, B.x), min(A.y, B.y), min(A.z, B.z));
}

// Solve equation ax^3+bx^2+cx+d=0 with Cardano formula
// return 1: two or three real roots, r, u, v;
// return 0: one real root r and two complex roots u+vi, u-vi;
bool solveCubic(double a, double b, double c, double d, double &r, double &u, double &v) {
	b /= a, c /= a, d /= a;		// now a=1
	double p = c - b * b / 3, q = (b*b / 13.5 - c / 3) * b + d;		// => t^3+pt+q=0, x=t-b/3
	b /= 3, p /= 3, q /= -2; a = q * q + p * p * p;
	if (a > 0) {
		a = sqrt(a);
		//u = cbrt(q + a), v = cbrt(q - a);
		u = q + a; u = u > 0 ? pow(u, 1. / 3) : -pow(-u, 1. / 3);
		v = q - a; v = v > 0 ? pow(v, 1. / 3) : -pow(-v, 1. / 3);
		r = u + v;
		v = sqrt(0.75)*(u - v);
		u = -0.5 * r - b, r -= b;
		return 0;
	}
	else {
		a = -a; c = pow(q*q + a, 1.0 / 6);
		a = sqrt(a); u = atan2(a, q) / 3;
		d = c * sin(u), c *= cos(u);
		r = 2 * c - b;
		c = -c, d *= sqrt(3);
		u = c - d - b, v = u + 2 * d;
		return 1;
	}
}

// solve ax^4+bx^3+cx^2+dx+e=0, return minimum possitive real root
inline double pick_random(double min, double max);
double solveQuartic(double a, double b, double c, double d, double e) {
	b /= a, c /= a, d /= a, e /= a, a = 1;
	double x = 0, dx;
	double a_ = 4 * a, b_ = 3 * b, c_ = 2 * c, d_ = d;
	double r, u, v;
	if (solveCubic(a_, b_, c_, d_, r, u, v)) {
		double mi = min(min(u, v), r), ma = max(max(u, v), r);	// to minimas
		if ((((a*mi + b)*mi + c)*mi + d)*mi + e > 0 && (((a*ma + b)*ma + c)*ma + d)*ma + e > 0) return NAN;
		// both minimas with values greater than 0 => no real root
	}
	else {
		if ((((a*r + b)*r + c)*r + d)*r + e > 0) return NAN;	// minima with value greate 0
	}
	x = -0.25*b;	// third derivative equal to zero
	unsigned n = 0; do {
		if (++n > 30 && abs(dx) > 1) x = pick_random(-2, 2) - 0.25*b, n -= 30;
		u = (((a*x + b)*x + c)*x + d)*x + e, v = ((a_*x + b_)*x + c_)*x + d_;
		dx = u / v; x -= dx;
	} while (abs(dx) > ERR_EPSILON_SN && ++n < 60);		// finding one root x using Newton's method
	if (n == 60) return NAN;
	c_ = b + x, d_ = x * c_ + c, e = x * d_ + d, d = d_, c = c_, b = 1, a = 0;	// Euclid division
	if (solveCubic(b, c, d, e, r, u, v)) {
		if (x < ERR_EPSILON) x = INFINITY; if (r < ERR_EPSILON) r = INFINITY; if (u < ERR_EPSILON) u = INFINITY; if (v < ERR_EPSILON) v = INFINITY;
		x = min(min(u, v), min(x, r));
		if (x == INFINITY) x = NAN;
	}
	else {
		if (x < ERR_EPSILON) x = r > ERR_EPSILON ? r : NAN;
		else x = r < ERR_EPSILON ? x : min(x, r);
	}
	if (isnan(x)) return x; 
	/*n = 0; do {
		u = (((a*x + b)*x + c)*x + d)*x + e, v = ((a_*x + b_)*x + c_)*x + d_;
		dx = u / v; x -= dx;
	} while (abs(dx) > ERR_EPSILON_SN && ++n < 20);		// dispose error*/
	return x;
}


/* About random, for solving render equations */

// linear congruence method producing random float value
extern double RAND_LCG_DV = 0.36787944117;	// random number seed
#define RAND_LCG_TMS 13.35717028437795
#define RAND_LCG_ADD 0.841470984807897
inline double pick_random(double max) {
	RAND_LCG_DV = fmod(RAND_LCG_DV * RAND_LCG_TMS + RAND_LCG_ADD, max);
	return RAND_LCG_DV;
}
inline double pick_random(double min, double max) {
	RAND_LCG_DV = fmod(RAND_LCG_DV * RAND_LCG_TMS + RAND_LCG_ADD, max - min);
	return RAND_LCG_DV + min;
}

// calculate inverse error function
inline double erfinv(double x) {
	double n = log(1 - x * x);
	double t = 0.5 * n + 2 / (PI*0.147);
	if (signbit(x)) return -sqrt(-t + sqrt(t*t - n / 0.147));
	return sqrt(-t + sqrt(t*t - n / 0.147));
}
// produce random normal-distributed number with given median and variance
inline double randnor(double median, double variance) {
	return erfinv(2.0 * double(rand()) / double(RAND_MAX + 1) - 1)*sqrt(2)*variance + median;
}
inline double randnor_0(double variance) {
	return erfinv(2.0 * double(rand()) / double(RAND_MAX + 1) - 1) * 1.41421356237309504876 * variance;
}
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

// Randomly rotate a vector, for solving render equations, return the cosine of rotate angle
double rotate_normal(point &N) {
	double m = N.mod();
	double x = acos(N.z / m), z = atan2(N.x, -N.y);
	RAND_LCG_DV = fmod(RAND_LCG_DV * RAND_LCG_TMS + RAND_LCG_ADD, PI);
	double rx = RAND_LCG_DV - PI / 2;
	RAND_LCG_DV = fmod(RAND_LCG_DV * RAND_LCG_TMS + RAND_LCG_ADD, 2 * PI);
	double rz = RAND_LCG_DV;
	double nx = m * sin(rx)*sin(rz), ny = -m * sin(rx)*cos(rz), nz = m * cos(rx);
	N.x = cos(z)*nx - cos(x)*sin(z)*ny + sin(x)*sin(z)*nz;
	N.y = sin(z)*nx + cos(x)*cos(z)*ny - sin(x)*cos(z)*nz;
	N.z = sin(x)*ny + cos(x)*nz;
	return cos(rx);
}
