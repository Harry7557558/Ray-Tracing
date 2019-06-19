#pragma once

/*
	This header file defines all classes, structs, functions "Object.h" needs.
	For convinience in debugging, all console outputs are GeoGebra commands. (This also for "Object.h")
	GeoGebra: https://www.geogebra.org/3d
*/

#pragma warning(disable: 4244)	// conversion from 'type1' to 'type2', possible loss of data
#pragma warning(disable: 4305)	// truncation from 'type1' to 'type2'
#pragma warning(disable: 4018)	// signed/unsigned mismatch, sometimes cause problems
#pragma warning(disable: 4010)

#include "BitMap.h"
#include <queue>
#include <stack>
#include <initializer_list>

using namespace std;

#ifndef PI
#define PI 3.1415926535897932384626433832795029L
#endif

/* Time Recorder */
#include <chrono>
typedef chrono::high_resolution_clock NTime;
typedef chrono::duration<double> fsec;

#include <cstdlib>
#include <algorithm>

#define ERR_EPSILON 1e-7	// minimum distance of intersection
#define ERR_UPSILON 1e+8
#define ERR_ZETA 1e-5	// minumum distance of SDF test, may cause threads jointing problems when too large

// For Debugging
extern ofstream fout("IMAGE\\Log.txt");
void WARN(string s) {
	cout << s << "\a\n"; fout << s << endl;
}

#define TeX_space "\\quad"
#define TeX_endl "\\\\[3pt]"

string uint2str(unsigned n, unsigned digits) {
	string s;
	for (unsigned i = 0; i < digits; i++) {
		s = "0" + s; s[0] += n % 10; n /= 10;
	}
	return s;
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
	point(const point &A, const point &B) {
		x = B.x - A.x, y = B.y - A.y, z = B.z - A.z;
	}
	point(const initializer_list<double> &a) {
		x = *(a.begin()), y = *(a.begin() + 1), z = *(a.begin() + 2);
	}
	inline point& operator = (const point &a) {
		x = a.x, y = a.y, z = a.z;
		return *this;
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
	friend inline point unitVec(const point &v) {
		return v / v.mod();
	}
	friend ostream& operator << (ostream& os, const point &a) {
		os << noshowpos << "(" << a.x << "," << a.y << "," << a.z << ")";
		return os;
	}
};
typedef point vec3;
class point2D {
public:
	double x, y;
	point2D() { x = 0, y = 0; }
	point2D(const double &x, const double &y) { this->x = x, this->y = y; }
	point2D(const point2D &another) { x = another.x, y = another.y; }
	inline point2D& operator = (const point2D &another) { x = another.x, y = another.y; return *this; }
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
		//return hypot(x, y);
		return sqrt(x * x + y * y);
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
typedef point2D vec2;


#define NAP (point(NAN,NAN,NAN));

enum matrix_type {
	Zero, Identity,
	Stretching, Rotation, Shearing, Reflection, Projection,
};
class matrix3D {
	double p[3][3];
public:
	matrix3D() {}
	matrix3D(const matrix_type &Transformation, const double &x, const double &y, const double &z) {
		switch (Transformation) {
		case Stretching: {
			p[0][0] = x, p[0][1] = p[0][2] = 0;
			p[1][1] = y, p[1][0] = p[1][2] = 0;
			p[2][2] = z, p[2][0] = p[2][1] = 0;
			return;
		}
		case Rotation: {
			p[0][0] = cos(y)*cos(z), p[0][1] = sin(x)*sin(y)*cos(z) - cos(x)*sin(z), p[0][2] = cos(x)*sin(y)*cos(z) + sin(x)*sin(z);
			p[1][0] = cos(y)*sin(z), p[1][1] = sin(x)*sin(y)*sin(z) + cos(x)*cos(z), p[1][2] = cos(x)*sin(y)*sin(z) - sin(x)*cos(z);
			p[2][0] = -sin(y), p[2][1] = sin(x)*cos(y), p[2][2] = cos(x)*cos(y);
			return;
		}
		}
	}
	matrix3D(const double& _00, const double& _01, const double& _02,
		const double& _10, const double& _11, const double& _12, const double& _20, const double& _21, const double& _22) {
		p[0][0] = _00, p[0][1] = _01, p[0][2] = _02, p[1][0] = _10, p[1][1] = _11, p[1][2] = _12, p[2][0] = _20, p[2][1] = _21, p[2][2] = _22;
	}
	matrix3D(const initializer_list<double> &l) {
		for (int i = 0; i < 9; i++) *(&p[0][0] + i) = *(l.begin() + i);
	}
	matrix3D(const initializer_list<initializer_list<double>> &l) {
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				p[i][j] = *((l.begin() + i)->begin() + j);
			}
		}
	}
	matrix3D(const matrix3D &other) {
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				p[i][j] = other.p[i][j];
			}
		}
	}
	~matrix3D() {}

	inline point operator * (const point &P) const {
		return point(p[0][0] * P.x + p[0][1] * P.y + p[0][2] * P.z,
			p[1][0] * P.x + p[1][1] * P.y + p[1][2] * P.z, p[2][0] * P.x + p[2][1] * P.y + p[2][2] * P.z);
	}
	friend void operator *= (point &P, const matrix3D &A) {
		double x = P.x, y = P.y, z = P.z;
		P.x = A.p[0][0] * x + A.p[0][1] * y + A.p[0][2] * z;
		P.y = A.p[1][0] * x + A.p[1][1] * y + A.p[1][2] * z;
		P.z = A.p[2][0] * x + A.p[2][1] * y + A.p[2][2] * z;
	}
	matrix3D operator * (const matrix3D &A) const {
		matrix3D R;
		for (int m = 0; m < 3; m++) {
			for (int n = 0; n < 3; n++) {
				R.p[m][n] = 0;
				for (int i = 0; i < 3; i++) R.p[m][n] += p[m][i] * A.p[i][n];
			}
		}
		return R;
	}

	friend ostream& operator << (ostream& os, const matrix3D &A) {
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				os << A.p[i][j] << "\t";
			}
			os << endl;
		}
		return os;
	}

	matrix3D invert() {
		matrix3D M(
			p[1][1] * p[2][2] - p[1][2] * p[2][1], p[0][2] * p[2][1] - p[0][1] * p[2][2], p[0][1] * p[1][2] - p[0][2] * p[1][1],
			p[1][2] * p[2][0] - p[1][0] * p[2][2], p[0][0] * p[2][2] - p[0][2] * p[2][0], p[0][2] * p[1][0] - p[0][0] * p[1][2],
			p[1][0] * p[2][1] - p[1][1] * p[2][0], p[0][1] * p[2][0] - p[0][0] * p[2][1], p[0][0] * p[1][1] - p[0][1] * p[1][0]);
		double det = p[0][0] * (p[1][1] * p[2][2] - p[1][2] * p[2][1]) - p[0][1] * (p[1][0] * p[2][2] - p[1][2] * p[2][0]) + p[0][2] * (p[1][0] * p[2][1] - p[1][1] * p[2][0]);
		for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) M.p[i][j] /= det;
		return M;
	}

	friend class matrix3D_affine;
};
class matrix3D_affine {
	double p[4][4];
public:
	matrix3D_affine() {
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				p[i][j] = i == j ? 1 : 0;
			}
		}
	}
	matrix3D_affine(const double& _00, const double& _01, const double& _02,
		const double& _10, const double& _11, const double& _12, const double& _20, const double& _21, const double& _22) {
		p[0][0] = _00, p[0][1] = _01, p[0][2] = _02, p[1][0] = _10, p[1][1] = _11, p[1][2] = _12, p[2][0] = _20, p[2][1] = _21, p[2][2] = _22;
		p[0][3] = p[1][3] = p[2][3] = 0, p[3][0] = p[3][1] = p[3][2] = 0, p[3][3] = 1;
	}
	matrix3D_affine(const double& _00, const double& _01, const double& _02, const double& _03, const double& _10, const double& _11, const double& _12, const double &_13,
		const double& _20, const double& _21, const double& _22, const double& _23, const double& _30, const double& _31, const double& _32, const double& _33) {
		p[0][0] = _00, p[0][1] = _01, p[0][2] = _02, p[0][3] = _03;
		p[1][0] = _10, p[1][1] = _11, p[1][2] = _12, p[1][3] = _13;
		p[2][0] = _20, p[2][1] = _21, p[2][2] = _22, p[2][3] = _23;
		p[3][0] = _30, p[3][1] = _31, p[3][2] = _32, p[3][3] = _33;
	}
	matrix3D_affine(const double *P) {
		for (int i = 0; i < 16; i++) *((&p[0][0]) + i) = *(P + i);
	}
	matrix3D_affine(const matrix3D &other) {
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				this->p[i][j] = other.p[i][j];
			}
			this->p[i][3] = this->p[3][i] = 0;
		}
		this->p[3][3] = 1;
	}
	matrix3D_affine(const matrix3D_affine &other) {
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				this->p[i][j] = other.p[i][j];
			}
		}
	}
	matrix3D_affine& operator = (const matrix3D_affine& other) {
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				this->p[i][j] = other.p[i][j];
			}
		}
		return *this;
	}
	~matrix3D_affine() {}

	inline double* operator [] (const unsigned &n) {
		return &p[n][0];
	}
	bool operator == (const matrix3D_affine &other) {
		double c = p[3][3] / other.p[3][3], dif;
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				dif = other.p[i][j] * c - this->p[i][j];
				if (abs(dif) > ERR_EPSILON) {
					//cout << *this << TeX_space << other << TeX_endl;
					return false;
				}
			}
		}
		return true;
	}
	inline point operator * (const point &P) const {
		double m = p[3][0] * P.x + p[3][1] * P.y + p[3][2] * P.z + p[3][3];
		return point((p[0][0] * P.x + p[0][1] * P.y + p[0][2] * P.z + p[0][3]) / m,
			(p[1][0] * P.x + p[1][1] * P.y + p[1][2] * P.z + p[1][3]) / m, (p[2][0] * P.x + p[2][1] * P.y + p[2][2] * P.z + p[2][3]) / m);
	}
	matrix3D_affine operator * (const matrix3D_affine &A) const {
		matrix3D_affine R;
		for (int m = 0; m < 4; m++) {
			for (int n = 0; n < 4; n++) {
				R.p[m][n] = 0;
				for (int i = 0; i < 4; i++) R.p[m][n] += p[m][i] * A.p[i][n];
			}
		}
		return R;
	}
	double det() const {
		bool sign = false;
		double P[4][4];
		for (int i = 0; i < 4; i++) for (int j = 0; j < 4; j++) P[i][j] = p[i][j];
		for (int i = 0; i < 4; i++) {
			if (abs(P[i][i]) < ERR_EPSILON) {
				P[i][i] = 0;
				for (int j = i + 1; j < 4; j++) {
					if (abs(P[i][j]) > ERR_EPSILON) {
						for (int k = 0; k < 4; k++) swap(P[k][j], P[k][i]);
						sign ^= 1;
						break;
					}
					if (j == 3) return 0;
				}
			}
			for (int j = i + 1; j < 4; j++) {
				if (abs(P[j][i]) > ERR_EPSILON) {
					double c = P[j][i] / P[i][i];
					for (int k = i; k < 4; k++) {
						P[j][k] -= c * P[i][k];
					}
				}
			}
		}

		double c = sign ? -1 : 1;
		for (int i = 0; i < 4; i++) c *= P[i][i];
		return c;
	}
	matrix3D_affine invert() const {
		matrix3D_affine R;
		R.p[0][0] = p[1][1] * p[2][2] * p[3][3] - p[1][1] * p[2][3] * p[3][2] - p[2][1] * p[1][2] * p[3][3] + p[2][1] * p[1][3] * p[3][2] + p[3][1] * p[1][2] * p[2][3] - p[3][1] * p[1][3] * p[2][2];
		R.p[1][0] = -p[1][0] * p[2][2] * p[3][3] + p[1][0] * p[2][3] * p[3][2] + p[2][0] * p[1][2] * p[3][3] - p[2][0] * p[1][3] * p[3][2] - p[3][0] * p[1][2] * p[2][3] + p[3][0] * p[1][3] * p[2][2];
		R.p[2][0] = p[1][0] * p[2][1] * p[3][3] - p[1][0] * p[2][3] * p[3][1] - p[2][0] * p[1][1] * p[3][3] + p[2][0] * p[1][3] * p[3][1] + p[3][0] * p[1][1] * p[2][3] - p[3][0] * p[1][3] * p[2][1];
		R.p[3][0] = -p[1][0] * p[2][1] * p[3][2] + p[1][0] * p[2][2] * p[3][1] + p[2][0] * p[1][1] * p[3][2] - p[2][0] * p[1][2] * p[3][1] - p[3][0] * p[1][1] * p[2][2] + p[3][0] * p[1][2] * p[2][1];
		R.p[0][1] = -p[0][1] * p[2][2] * p[3][3] + p[0][1] * p[2][3] * p[3][2] + p[2][1] * p[0][2] * p[3][3] - p[2][1] * p[0][3] * p[3][2] - p[3][1] * p[0][2] * p[2][3] + p[3][1] * p[0][3] * p[2][2];
		R.p[1][1] = p[0][0] * p[2][2] * p[3][3] - p[0][0] * p[2][3] * p[3][2] - p[2][0] * p[0][2] * p[3][3] + p[2][0] * p[0][3] * p[3][2] + p[3][0] * p[0][2] * p[2][3] - p[3][0] * p[0][3] * p[2][2];
		R.p[2][1] = -p[0][0] * p[2][1] * p[3][3] + p[0][0] * p[2][3] * p[3][1] + p[2][0] * p[0][1] * p[3][3] - p[2][0] * p[0][3] * p[3][1] - p[3][0] * p[0][1] * p[2][3] + p[3][0] * p[0][3] * p[2][1];
		R.p[3][1] = p[0][0] * p[2][1] * p[3][2] - p[0][0] * p[2][2] * p[3][1] - p[2][0] * p[0][1] * p[3][2] + p[2][0] * p[0][2] * p[3][1] + p[3][0] * p[0][1] * p[2][2] - p[3][0] * p[0][2] * p[2][1];
		R.p[0][2] = p[0][1] * p[1][2] * p[3][3] - p[0][1] * p[1][3] * p[3][2] - p[1][1] * p[0][2] * p[3][3] + p[1][1] * p[0][3] * p[3][2] + p[3][1] * p[0][2] * p[1][3] - p[3][1] * p[0][3] * p[1][2];
		R.p[1][2] = -p[0][0] * p[1][2] * p[3][3] + p[0][0] * p[1][3] * p[3][2] + p[1][0] * p[0][2] * p[3][3] - p[1][0] * p[0][3] * p[3][2] - p[3][0] * p[0][2] * p[1][3] + p[3][0] * p[0][3] * p[1][2];
		R.p[2][2] = p[0][0] * p[1][1] * p[3][3] - p[0][0] * p[1][3] * p[3][1] - p[1][0] * p[0][1] * p[3][3] + p[1][0] * p[0][3] * p[3][1] + p[3][0] * p[0][1] * p[1][3] - p[3][0] * p[0][3] * p[1][1];
		R.p[3][2] = -p[0][0] * p[1][1] * p[3][2] + p[0][0] * p[1][2] * p[3][1] + p[1][0] * p[0][1] * p[3][2] - p[1][0] * p[0][2] * p[3][1] - p[3][0] * p[0][1] * p[1][2] + p[3][0] * p[0][2] * p[1][1];
		R.p[0][3] = -p[0][1] * p[1][2] * p[2][3] + p[0][1] * p[1][3] * p[2][2] + p[1][1] * p[0][2] * p[2][3] - p[1][1] * p[0][3] * p[2][2] - p[2][1] * p[0][2] * p[1][3] + p[2][1] * p[0][3] * p[1][2];
		R.p[1][3] = p[0][0] * p[1][2] * p[2][3] - p[0][0] * p[1][3] * p[2][2] - p[1][0] * p[0][2] * p[2][3] + p[1][0] * p[0][3] * p[2][2] + p[2][0] * p[0][2] * p[1][3] - p[2][0] * p[0][3] * p[1][2];
		R.p[2][3] = -p[0][0] * p[1][1] * p[2][3] + p[0][0] * p[1][3] * p[2][1] + p[1][0] * p[0][1] * p[2][3] - p[1][0] * p[0][3] * p[2][1] - p[2][0] * p[0][1] * p[1][3] + p[2][0] * p[0][3] * p[1][1];
		R.p[3][3] = p[0][0] * p[1][1] * p[2][2] - p[0][0] * p[1][2] * p[2][1] - p[1][0] * p[0][1] * p[2][2] + p[1][0] * p[0][2] * p[2][1] + p[2][0] * p[0][1] * p[1][2] - p[2][0] * p[0][2] * p[1][1];
		double det = p[0][0] * R.p[0][0] + p[0][1] * R.p[1][0] + p[0][2] * R.p[2][0] + p[0][3] * R.p[3][0];
		if (det == 0) return false;
		for (int i = 0; i < 4; i++) for (int j = 0; j < 4; j++) R.p[i][j] /= det;
		return R;
	}

	void scale(double x, double y, double z) {
		*this = matrix3D_affine(
			x, 0, 0,
			0, y, 0,
			0, 0, z)*(*this);
	}
	void translate(double x, double y, double z) {
		*this = matrix3D_affine(
			1, 0, 0, x,
			0, 1, 0, y,
			0, 0, 1, z,
			0, 0, 0, 1)*(*this);
	}
	void rotate(double rx, double ry, double rz) {
		*this = matrix3D_affine(
			cos(ry)*cos(rz), sin(rx)*sin(ry)*cos(rz) - cos(rx)*sin(rz), cos(rx)*sin(ry)*cos(rz) + sin(rx)*sin(rz),
			cos(ry)*sin(rz), sin(rx)*sin(ry)*sin(rz) + cos(rx)*cos(rz), cos(rx)*sin(ry)*sin(rz) - sin(rx)*cos(rz),
			-sin(ry), sin(rx)*cos(ry), cos(rx)*cos(ry)) * (*this);
	}
	void perspective(double x, double y, double z) {
		*this = matrix3D_affine(
			1, 0, 0, 0,
			0, 1, 0, 0,
			0, 0, 1, 0,
			x, y, z, 1)*(*this);
	}

	friend ostream& operator << (ostream& os, const matrix3D_affine &A) {
		os << "\\begin{bmatrix}";
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				os << A.p[i][j];
				if (j != 3) os << "&";
			}
			os << TeX_endl;
		}
		os << "\\end{bmatrix}";
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
	double ut, vt, wt;	// reserve
};
struct intersect2D {
	bool meet = 0;
	double dist;
	point2D intrs;
	point2D reflect;
	double ut;
};

/* Extreme points of an object */
class PMaxMin {
public:
	point max;	// maximum x,y,z coordinates
	point min;	// minimum x,y,z coordinates
	PMaxMin(const point &Max, const point &Min) {
		max = Max, min = Min;
	}
	~PMaxMin() {}
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


#define OS_Pointer "0x" << hex << uppercase << setw(8) << setfill('0')
class borderbox {
	// Use for border of subgroups of objects
public:
	point Min, Max;
	borderbox() { Min = point(INFINITY, INFINITY, INFINITY), Max = point(-INFINITY, -INFINITY, -INFINITY); }
	borderbox(const point &Min, const point &Max) {
		this->Min = Min, this->Max = Max;
		fix();
	}
	borderbox(const borderbox &other) {
		this->Min = other.Min, this->Max = other.Max;
	}
	bool meet(const ray &a) const {
		// http://www.cs.utah.edu/~awilliam/box/box.pdf
		double tmin, tmax, tymin, tymax, tzmin, tzmax;
		tmin = ((a.dir.x < 0 ? Max : Min).x - a.orig.x) / a.dir.x;
		tmax = ((a.dir.x < 0 ? Min : Max).x - a.orig.x) / a.dir.x;
		tymin = ((a.dir.y < 0 ? Max : Min).y - a.orig.y) / a.dir.y;
		tymax = ((a.dir.y < 0 ? Min : Max).y - a.orig.y) / a.dir.y;
		if ((tmin > tymax) || (tymin > tmax)) return 0;
		if (tymin > tmin) tmin = tymin;
		if (tymax < tmax) tmax = tymax;
		tzmin = ((a.dir.z < 0 ? Max : Min).z - a.orig.z) / a.dir.z;
		tzmax = ((a.dir.z < 0 ? Min : Max).z - a.orig.z) / a.dir.z;
		if ((tmin > tzmax) || (tzmin > tmax)) return 0;
		if (tzmin > tmin) tmin = tzmin;
		if (tzmax < tmax) tmax = tzmax;
		return tmax > 0;
	}
	inline bool contain(const point &a) const {
		return (a.x > Min.x && a.x<Max.x && a.y>Min.y && a.y<Max.y && a.z>Min.z && a.z < Max.z);
	}
	void fix() {
		if (Max.x < Min.x) swap(Max.x, Min.x);
		if (Max.y < Min.y) swap(Max.y, Min.y);
		if (Max.z < Min.z) swap(Max.z, Min.z);
	}
	~borderbox() {}
	friend ostream& operator << (ostream& os, const borderbox &a) {
		point A = a.Min, B = point(a.Max.x, a.Min.y, a.Min.z), C = point(a.Max.x, a.Max.y, a.Min.z), D = point(a.Min.x, a.Max.y, a.Min.z),
			E = point(a.Min.x, a.Min.y, a.Max.z), F = point(a.Max.x, a.Min.y, a.Max.z), G = a.Max, H = point(a.Min.x, a.Max.y, a.Max.z);
		os << "Borderbox_{" << OS_Pointer << (unsigned)&a << "}: Polyline(" << A << "," << B << "," << F << "," << G << "," << C << "," << D << "," << H << "," << E << "," <<
			A << "," << E << "," << F << "," << B << "," << C << "," << G << "," << H << "," << D << "," << A << ")";
		return os;
	}
	inline point center() const {
		return 0.5*(Min + Max);
	}
};
class borderbox2D {
	// Use for border of subgroups of objects
public:
	point2D Min, Max;
	borderbox2D() { Min = point2D(INFINITY, INFINITY), Max = point2D(-INFINITY, -INFINITY); }
	borderbox2D(const point2D &Min, const point2D &Max) {
		this->Min = Min, this->Max = Max;
		fix();
	}
	borderbox2D(const borderbox2D &other) {
		this->Min = other.Min, this->Max = other.Max;
	}
	void fix() {
		if (Max.x < Min.x) swap(Max.x, Min.x);
		if (Max.y < Min.y) swap(Max.y, Min.y);
	}
	~borderbox2D() {}
};

inline point PMax(point A, point B) {
	return point(max(A.x, B.x), max(A.y, B.y), max(A.z, B.z));
}
inline point PMin(point A, point B) {
	return point(min(A.x, B.x), min(A.y, B.y), min(A.z, B.z));
}
inline point2D PMax(point2D A, point2D B) {
	return point2D(max(A.x, B.x), max(A.y, B.y));
}
inline point2D PMin(point2D A, point2D B) {
	return point2D(min(A.x, B.x), min(A.y, B.y));
}


#define ERR_EPSILON_SN 1e-10L
#define ERR_ZETA_SN 1e-8L

// Solve equation ax^3+bx^2+cx+d=0 with Cardano formula
// return 1: two or three real roots, r, u, v;
// return 0: one real root r and two complex roots u+vi, u-vi;
bool solveCubic(double a, double b, double c, double d, double &r, double &u, double &v) {
	/*if (a == 0) {
		a = b, b = c, c = d;
		double delta = b * b - 4 * a*c; r = NAN;
		if (delta >= 0) {
			delta = sqrt(delta), a *= 2, b = -b;
			u = (b + delta) / a, v = (b - delta) / a; return 1;
		}
		else {
			u = -b, v = sqrt(-delta);
			a *= 2; u /= a, v /= a; return 0;
		}
	}*/
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
namespace MyComplex {
	class complex {
	public:
		double re, im;
		complex() {}
		complex(double r) {
			this->re = r, this->im = 0;
		}
		complex(double re, double im) {
			this->re = re, this->im = im;
		}
		complex(const complex &c) {
			re = c.re, im = c.im;
		}
		complex& operator = (const complex &c) {
			re = c.re, im = c.im; return *this;
		}
		~complex() {}

		inline complex operator + (const double &c) const {
			return complex(re + c, im);
		}
		inline complex operator + (const complex &c) const {
			return complex(re + c.re, im + c.im);
		}
		inline complex operator - () const {
			return complex(-re, -im);
		}
		inline complex operator - (const complex &c) const {
			return complex(re - c.re, im - c.im);
		}
		inline complex operator * (const double &c) const {
			return complex(c*re, c*im);
		}
		inline complex operator / (const double &c) const {
			return complex(re / c, im / c);
		}
		inline friend complex operator / (const double &a, const complex &c) {
			double m = a / (c.re*c.re + c.im * c.im);
			return complex(m*c.re, -m * c.im);
		}

		inline complex sqrt() const {
			double m = std::hypot(re, im);
			return complex(std::sqrt(0.5*(m + re)),
				im > 0 ? std::sqrt(0.5*(m - re)) : -std::sqrt(0.5*(m - re)));
		}
		inline complex cbrt() const {
			double s = std::pow(re*re + im * im, 1. / 6), c = std::atan2(im, re) / 3;
			return complex(s*std::cos(c), s*std::sin(c));
		}

		friend ostream& operator << (ostream &os, const complex &c) {
			os << noshowpos << c.re << showpos << c.im << "i";
			return os;
		}
	};
	inline complex sqrt(double r) {
		if (r >= 0) return complex(std::sqrt(r), 0);
		return complex(0, std::sqrt(-r));
	}
}
double solveQuartic_GeneralFormula(double a, double b, double c, double d, double e) {
	// https://en.wikipedia.org/wiki/Quartic_function#General_formula_for_roots
	// This way occurs a high error (over 1e-4)
	b /= a, c /= a, d /= a, e /= a;
	double b2 = b * b, c2 = c * c;
	double p = c - 0.375 * b2, q = 0.125 * b2*b - 0.5*b*c + d; p *= 2;
	double Delta0 = c2 - 3 * b*d + 12 * e, Delta1 = 2 * c2*c - 9 * b*c*d + 27 * b2 * e + 27 * d*d - 72 * c*e;
	MyComplex::complex Q = ((MyComplex::sqrt(Delta1*Delta1 - 4 * Delta0*Delta0*Delta0) + Delta1)*0.5).cbrt();
	MyComplex::complex S2 = ((Q + Delta0 / Q) - p) / 3, S = S2.sqrt()*0.5;
	MyComplex::complex K;
	K = (-S2 - q / S - p).sqrt()*0.5;
	MyComplex::complex R1 = S + K, R2 = S - K;
	K = (q / S - S2 - p).sqrt()*0.5, S = -S;
	MyComplex::complex R3 = S + K, R4 = S - K;
	b /= 4; R1.re -= b, R2.re -= b, R3.re -= b, R4.re -= b; b *= 4;

	double r, r1 = INFINITY, r2 = INFINITY, r3 = INFINITY, r4 = INFINITY, dx, u, v;
	double a_ = 4, b_ = 3 * b, c_ = 2 * c, d_ = d;
	if (abs(R1.im) < ERR_ZETA_SN && abs(R2.im) < ERR_ZETA_SN) {
		r1 = R1.re, r2 = R2.re;
		if (r1 < ERR_EPSILON) r1 = INFINITY;
		if (r2 < ERR_EPSILON) r2 = INFINITY;
	}
	if (abs(R3.im) < ERR_ZETA_SN && abs(R4.im) < ERR_ZETA_SN) {
		r3 = R3.re, r4 = R4.re;
		if (r3 < ERR_EPSILON) r3 = INFINITY;
		if (r4 < ERR_EPSILON) r4 = INFINITY;
	}
	if (abs(R1.im) < ERR_ZETA_SN && abs(R4.im) < ERR_ZETA_SN) {
		// R1,R4 real, R2,R3 conjugate (I didn't proof it, but it's satisfied in over 20,000,000 random-value tests)
		// Also high error occurs in this situation
		r1 = R1.re, r4 = R2.re;
		if (r1 < ERR_EPSILON) r1 = INFINITY;
		if (r4 < ERR_EPSILON) r4 = INFINITY;
	}
	r = min(min(r1, r2), min(r3, r4));
	if (r == INFINITY) return NAN;
	unsigned n = 0; do {
		u = (((a*r + b)*r + c)*r + d)*r + e, v = ((a_*r + b_)*r + c_)*r + d_;
		dx = u / v; r -= dx;
	} while (abs(dx) > ERR_EPSILON_SN && ++n < 30);
	//if (r < ERR_EPSILON) fout << "W";
	//if (isnan(r)) fout << "M";
	return r;
}
double solveQuartic(double a, double b, double c, double d, double e) {
	b /= a, c /= a, d /= a, e /= a;
	//double B = b, C = c, D = d, E = e;
	double x = 0, dx;
	double a_ = 4, b_ = 3 * b, c_ = 2 * c, d_ = d;
	double r, u, v;
	if (solveCubic(a_, b_, c_, d_, r, u, v)) {
		double mi = min(min(u, v), r), ma = max(max(u, v), r);	// two minimas
		if ((((mi + b)*mi + c)*mi + d)*mi + e > 0 && (((ma + b)*ma + c)*ma + d)*ma + e > 0) return NAN;
		// both minimas with values greater than 0 => no real root
	}
	else {
		if ((((r + b)*r + c)*r + d)*r + e > 0) return NAN;	// minima with value greate 0
	}
	x = -0.25*b;	// third derivative equal to zero
	unsigned n = 0; do {
		if (++n > 30 && abs(dx) > 1) x = pick_random(-2, 2) - 0.25*b, n -= 30;
		u = (((x + b)*x + c)*x + d)*x + e, v = ((a_*x + b_)*x + c_)*x + d_;
		dx = u / v; x -= dx;
	} while (abs(dx) > ERR_EPSILON_SN && ++n < 60);		// finding one root x using Newton's method
	//if (n == 60 && abs(dx) > 1e-3) return NAN;
	c_ = b + x, d_ = x * c_ + c, e = x * d_ + d, d = d_, c = c_, b = 1, a = 0;	// Euclid division
	if (solveCubic(b, c, d, e, r, u, v)) {
		if (x < ERR_EPSILON) x = INFINITY; if (r < ERR_EPSILON) r = INFINITY; if (u < ERR_EPSILON) u = INFINITY; if (v < ERR_EPSILON) v = INFINITY;
		x = min(min(u, v), min(x, r));
		if (x == INFINITY) x = NAN;
	}
	else {
		if (x < ERR_EPSILON) x = r > ERR_EPSILON ? r : NAN;
		else x = r < ERR_EPSILON ? x : (x < r ? x : r);
	}
	/*b = B, c = C, d = D, e = E; a = 1; a_ *= a, b_ *= b, c_ *= d, d_ *= d;
	n = 0; do {
		u = (((a*x + b)*x + c)*x + d)*x + e, v = ((a_*x + b_)*x + c_)*x + d_;
		dx = u / v; x -= dx;
	} while (abs(dx) > ERR_EPSILON_SN && ++n < 20);		// dispose error*/
	return x;
}
void solveQuartic(double a, double b, double c, double d, double e, double &r1, double &r2, double &r3, double &r4) {
	r1 = r2 = r3 = r4 = NAN;
	b /= a, c /= a, d /= a, e /= a;
	double x = 0, dx;
	double a_ = 4, b_ = 3 * b, c_ = 2 * c, d_ = d;
	double r, u, v;
	if (solveCubic(a_, b_, c_, d_, r, u, v)) {
		double mi = min(min(u, v), r), ma = max(max(u, v), r);	// two minimas
		if ((((mi + b)*mi + c)*mi + d)*mi + e > 0 && (((ma + b)*ma + c)*ma + d)*ma + e > 0) return;
	}
	else {
		if ((((r + b)*r + c)*r + d)*r + e > 0) return;	// one minima with value greate 0
	}
	x = -0.25*b;	// fourth derivative equals to zero (Note that this is not always the best point, sometimes trig extreme point
	unsigned n = 0; do {
		if (++n > 30) x = pick_random(-2, 2) - 0.25*b, n -= 30;		// seldom occurs
		u = (((x + b)*x + c)*x + d)*x + e, v = ((a_*x + b_)*x + c_)*x + d_;
		dx = u / v; x -= dx;
	} while (abs(dx) > ERR_EPSILON_SN && ++n < 120);		// finding one root x using Newton's method
	r1 = x;
	c_ = b + x, d_ = x * c_ + c, e = x * d_ + d, d = d_, c = c_, b = 1, a = 0;	// Euclid division
	if (!solveCubic(b, c, d, e, r2, r3, r4)) r3 = r4 = NAN;
	return;
}
void solveQuintic(double a, double b, double c, double d, double e, double f, double &r1, double &r2, double &r3, double &r4, double &r5) {
	//cout << setprecision(15) << noshowpos << a << "*x^5" << showpos << b << "*x^4" << c << "*x^3" << d << "*x^2" << e << "*x" << f << endl;
	r1 = r2 = r3 = r4 = r5 = NAN;
	b /= a, c /= a, d /= a, e /= a, f /= a, a = 1;
	double x = -0.2*b, u, v, dx;
	double a_ = 5, b_ = 4 * b, c_ = 3 * c, d_ = 2 * d, e_ = e;
	unsigned n = 0; do {
		if (++n > 30) x = pick_random(-2, 2) - 0.2*b, n -= 30;
		u = ((((a*x + b)*x + c)*x + d)*x + e)*x + f, v = (((a_*x + b_)*x + c_)*x + d_)*x + e_;
		dx = u / v; x -= dx;
		//const double MN = 10; if (v < 1 / MN) a *= MN, b *= MN, c *= MN, d *= MN, e *= MN, f *= MN, a_ *= MN, b_ *= MN, c_ *= MN, d_ *= MN, e_ *= MN;
	} while (abs(dx) > ERR_EPSILON_SN && ++n < 120);
	r1 = x;
	//if (a != 1) b /= a, c /= a, d /= a, e /= a, f /= a, a = 1;
	c_ = b + x, d_ = x * c_ + c, e_ = x * d_ + d, f = x * e_ + e, e = e_, d = d_, c = c_, b = 1, a = 0;
	solveQuartic(b, c, d, e, f, r2, r3, r4, r5);

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

void mess_vec(vec3 &v) {	// preventing division by zero
	if (v.x == 0) v.x = pick_random(-ERR_EPSILON, ERR_EPSILON);
	if (v.y == 0) v.y = pick_random(-ERR_EPSILON, ERR_EPSILON);
	if (v.z == 0) v.z = pick_random(-ERR_EPSILON, ERR_EPSILON);
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

// Output data of a spacial quadrilateral ABCD with given points, for debug
void print_quadrilateral(ostream& os, point A, point B, point C, point D) {
	os << "Surface(" << "(" << "(1-u)*" << A << "+u*" << B << ")*(1-v)+(" << "(1-u)*" << D << "+u*" << C << ")*v,u,0,1,v,0,1)";
}

// correct unit vector t, make it perpendicular to unit vector s (their plane don't change)
inline void correct_vector(const vec3 &s, vec3 &t) {
	t -= dot(s, t)*s; t /= t.mod();
}
// Get rotation angles, two vectors must be perpendicular unit vector
inline void getRotationAngle_ik(const vec3 &i, const vec3 &k, double &rx, double &ry, double &rz) {
	vec3 j = cross(k, i);
	rx = atan2(j.z, k.z), rz = atan2(i.y, i.x), ry = atan2(-i.z, hypot(i.x, i.y));
}
inline void getRotationAngle_ij(const vec3 &i, const vec3 &j, double &rx, double &ry, double &rz) {
	vec3 k = cross(i, j);
	rx = atan2(j.z, k.z), rz = atan2(i.y, i.x), ry = atan2(-i.z, hypot(i.x, i.y));
}
inline void getRotationAngle_jk(const vec3 &j, const vec3 &k, double &rx, double &ry, double &rz) {
	vec3 i = cross(j, k);
	rx = atan2(j.z, k.z), rz = atan2(i.y, i.x), ry = atan2(-i.z, hypot(j.z, k.z));
}


double sigmoid(double a) {
	return 1 / (exp(-a) + 1);
}
template<typename T> inline T clamp(const T &x, const T &min, const T &max) {
	return x < min ? min : x > max ? max : x;
}
template<typename T> inline T mix(const T &x, const T &y, const double &a) {
	return (1 - a)*x + a * y;
}
template<typename T> inline T Bezier(const T &a, const T &b, const T &c, const double &t) {
	return (1 - t)*(1 - t)*a + 2 * t*(1 - t)*b + t * t*c;
}
template<typename T> T Catmull_Rom(vector<T> p, double t) {
	// Construct a smooth curve between P1 and P2				\
	                |-0.5  1.5 -1.5  0.5 | |P0|					\
	P = [t³ t² t 1] | 1.0 -2.5  2.0 -0.5 | |P1|   0 ≤ t ≤ 1		\
	                |-0.5  0.0  0.5  0.0 | |P2|					\
	                | 0.0  1.0  0.0  0.0 | |P3|					
	p.insert(p.begin(), 2 * p[0] - p[1]); p.push_back(2 * p.back() - p[p.size() - 2]);	// constructing endpoints
	unsigned n = t; t -= n; n++;
	return (((-0.5*t + 1.0)*t - 0.5)*t) * p[n - 1]
		+ (((1.5*t - 2.5)*t)*t + 1.0) * p[n]
		+ (((-1.5*t + 2.0)*t + 0.5)*t) * p[n + 1]
		+ (((0.5*t - 0.5)*t)*t) * p[n + 2];
}

void Catmull_Rom_Test() {
	vector<double> x({ 138.6779485436933, 168.31067737090004, 235.72513545298528, 251.28231808730646, 274.9885011491574, 299.4355024315747, 325.3641401554627, 348.32950499654896, 360.9234147481677, 354.25605076207387, 332.7723223622327, 320.17841261069844, 320.17841261069844, 333.5131405829267, 359.44177830680445, 386.8520524720798, 458.71141987816094, 517.9768775326739, 550.5728792427587, 569.0933347597334, 570.5749712010845, 540.9422423738595, 471.3053296297195, 415.74396307864185, 416.4847812993417, 425.3745999474973, 489.0849669260721, 565.3892436563156, 625.3955195315486, 660.2139759036689, 665.3997034483971, 619.4689737661065, 513.5319682086522, 470.56451140904556, 471.3053296297195, 482.4176029399455, 550.5728792427587, 624.6547013108323, 676.5119767586549, 707.6263420272492, 712.812069571968, 681.6977043033485, 570.5749712010845, 497.233967353591, 445.3766919058968, 440.93178258175135, 507.60542244313876, 545.3871516979514, 595.7627907042681, 627.6179741935744, 614.2832462212773, 573.5382440838049, 526.5693642545089, 451.00590574492605, 370.9975379113052, 290.2205678383551 });
	vector<double> y({ 442.9614502984681, 382.95517442323467, 289.6120786173552, 264.4242591141776, 171.08116330827397, 102.18506878488225, 69.58906707489012, 68.10743063351683, 85.88706792988718, 116.26061497779808, 157.74643533598612, 203.67716501830603, 258.4977133487139, 266.64671377624654, 268.8691684382727, 239.97725783167448, 159.9688899980538, 96.2585230194327, 88.1095225919454, 102.9258870055611, 125.15043362603622, 158.48725355670078, 228.12416630076132, 295.538624382801, 301.4651701482469, 293.316169720755, 235.53234850757514, 174.04443619101832, 120.7055243019117, 130.33616117077855, 165.89543576352006, 212.5669836664461, 303.6876248103258, 334.06117185825997, 342.9509905064583, 354.0632638166745, 319.98562566529296, 280.7222599691738, 247.38544003849177, 259.23853156941715, 289.6120786173552, 319.98562566529296, 387.4000837473489, 428.1450858848215, 480.7431795532375, 505.93099905642555, 508.1534537184577, 500.7452715116415, 517.043272366622, 550.3800922972955, 568.1597295936995, 568.1597295936995, 573.5773096391665, 578.0222189632526, 603.9508566871282, 604.6970658152852 });
	vector<vec2> p; for (int i = 0; i < x.size(); i++) p.push_back(vec2(x[i], 720 - y[i]));
	bitmap canvas(900, 720);
	point2D prev = p[0], cur;
	for (double t = 0; t < p.size(); t += 0.01) {
		cur = Catmull_Rom(p, t);
		canvas.line(prev.x, prev.y, cur.x, cur.y, 1, Black);
		prev = cur;
	}
	canvas.out("IMAGE\\CTMR.bmp");
}

