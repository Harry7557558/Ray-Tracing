#pragma once

#ifndef Object_Sign
#include "Object.h"
#endif

// https://www.iquilezles.org/www/index.htm

namespace XObjs {
	class XObjs_Comp {
	public:
		XObjs_Comp() {}
		~XObjs_Comp() {}
		virtual double SDF(const point &P) const {
			WARN("XObjs::XObjs_Comp.SDF() called!");
			return NAN;
		}
		virtual borderbox MaxMin() const {
			WARN("XObjs::XObjs_Comp.SDF() called!");
			return borderbox(point(NAN, NAN, NAN), point(NAN, NAN, NAN));
		}
	};

	class Sphere : public XObjs_Comp {
		point C; double r;
	public:
		Sphere() :r(0) {}
		Sphere(point C, double r) {
			this->C = C, this->r = r;
		}
		Sphere(const Sphere &other) {
			C = other.C, r = other.r;
		}
		~Sphere() {}

		double SDF(const point &P) const {
			return (P - C).mod() - r;
		}
		borderbox MaxMin() const {
			return borderbox(C - point(r, r, r), C + point(r, r, r));
		}
	};

	class Plane : public XObjs_Comp {
		point N; double D;	// Ax+By+Cz<D, N unit vector
	public:
		Plane() :D(0) { N.z = 1; }
		Plane(const point &N) {
			this->N = N / N.mod(), D = 0;
		}
		Plane(const double &z_int) {
			N.z = 1, this->D = z_int;
		}
		Plane(const point &P, const point &N) {
			this->N = N / N.mod(); D = dot(P, N);
		}
		Plane(const point &A, const point &B, const point &C) {
			N = cross(B - A, C - A); D = dot(A, N); N /= N.mod();
		}
		Plane(const Plane &p) :N(p.N), D(p.D) {}
		~Plane() {}

		double SDF(const point &P) const {
			return dot(N, P) - D;
		}
		borderbox MaxMin() const {
			return borderbox(point(-INFINITY, -INFINITY, -INFINITY), point(INFINITY, INFINITY, INFINITY));
		}
	};

	class Cylinder_std : public XObjs_Comp {
		point C; double r; point D;		// central point, radius, direction vector (unit)
	public:
		Cylinder_std() { r = 0; }
		Cylinder_std(const point &C, const double &r, const  point &d) {
			this->C = C, this->r = r, this->D = d;
			D /= D.mod();
		}
		Cylinder_std(const point &P1, const point &P2, const double &r) {
			this->C = P1, this->r = r, this->D = P2 - P1;
			D /= D.mod();
		}
		Cylinder_std(const Cylinder_std &other) {
			this->C = other.C, this->r = other.r, this->D = other.D;
		}
		~Cylinder_std() {}

		double SDF(const point &P) const {
			return cross(P - C, D).mod() - r;
		}
		borderbox MaxMin() const {
			borderbox A;
			A.Max.x = D.x == 0 ? r : INFINITY;
			A.Max.y = D.y == 0 ? r : INFINITY;
			A.Max.z = D.z == 0 ? r : INFINITY;
			A.Min = -A.Max + C; A.Max += C;
			return A;
		}
	};

	class Cone_std : public XObjs_Comp {
		point C; double alpha; point D;		// 0 < alpha < PI/2, D unit
	public:
		Cone_std() { alpha = 0; }
		Cone_std(const point &C, const point  &P, const double &r) {
			this->C = C, this->D = P - C;
			alpha = atan(r / P.mod());
			D /= D.mod();
		}
		Cone_std(const point &C, const double &alpha, const point &D) {
			this->C = C, this->alpha = alpha, this->D = D / D.mod();
		}
		Cone_std(const Cone_std &other) {
			C = other.C, alpha = other.alpha, D = other.D;
		}
		~Cone_std() {}

		double SDF(const point &P) const {
			point S = P - C; double u = dot(S, D), m = S.mod();
			//if (u < 0) return m;
			return m * sin(acos(u / m) - alpha);
		}
		borderbox MaxMin() const {
			borderbox B;
			B.Max = point(INFINITY, INFINITY, INFINITY); B.Min = -B.Max;	// Optimize
			return B;
		}
	};

	class Torus_xOy : public XObjs_Comp {
		point C; double R, r;
	public:
		Torus_xOy() { R = 0, r = 0; }
		Torus_xOy(double R, double r) {
			this->R = R, this->r = r;
		}
		Torus_xOy(point C, double R, double r) {
			this->C = C, this->R = R, this->r = r;
		}
		Torus_xOy(const Torus_xOy &other) {
			C = other.C, R = other.R, r = other.r;
		}
		~Torus_xOy() {}

		double SDF(const point &P) const {
			//return hypot(hypot(P.x - C.x, P.y - C.y) - R, P.z - C.z) - r;
			point S = P - C; double m = sqrt(S.x*S.x + S.y*S.y) - R; m = sqrt(m*m + S.z*S.z); return m - r;
		}
		borderbox MaxMin() const {
			borderbox B;
			B.Max = point(R + r, R + r, r), B.Min = -B.Max;
			B.Max += C, B.Min += C;
			return B;
		}
	};

	class Box : public XObjs_Comp {
		point O; point A, B, C; // parallelepiped
	public:

	};
}

// Windows null pointer: 0x00000000-0x0000FFFF
// 0xABCC    A: number of operated objects;  B: number of extra parameters;  CC: operation #
#define CSG_UnionOp_Sign 0x2001
#define CSG_IntersectionOp_Sign 0x2002
#define CSG_SubtractionOp_Sign 0x2003
#define CSG_ComplementOp_Sign 0x1004
#define CSG_RoundingOp_Sign 0x1105
#define CSG_OnionOp_Sign 0x1106

#define CSG_UnionOperator ((const XObjs::XObjs_Comp*)CSG_UnionOp_Sign)
#define CSG_IntersectionOperator ((const XObjs::XObjs_Comp*)CSG_IntersectionOp_Sign)
#define CSG_SubtractionOperator ((const XObjs::XObjs_Comp*)CSG_SubtractionOp_Sign)
#define CSG_ComplementOperator ((const XObjs::XObjs_Comp*)CSG_ComplementOp_Sign)
#define CSG_RoundingOperator ((const XObjs::XObjs_Comp*)CSG_RoundingOp_Sign)
#define CSG_OnionOperator ((const XObjs::XObjs_Comp*)CSG_OnionOp_Sign)

#include <thread>
#include <mutex>

#define XSolid_Sign 0x10000001
#define XSolid_Smooth 0x10000010
#define XSolid_Diffuse 0x10000100
#define XSolid_Crystal 0x10001000
class XSolid : public object {
	vector<const XObjs::XObjs_Comp*> objs;	// postfix expression
	vector<void**> objs_tp;
#define XSolid_Max_Parameter 4
	inline void objs_tp_push() {
		objs_tp.push_back(new void*[XSolid_Max_Parameter]);
		for (unsigned m = 0; m < XSolid_Max_Parameter; m++) this->objs_tp.back()[m] = 0;
	}
	void deep_copy(const XSolid &other) {
		// assign objs_dir
		for (unsigned i = 0; i < other.objs_tp.size(); i++) {
			if (other.objs_tp.at(i) == 0) this->objs_tp.push_back(0);
			else {
				this->objs_tp.push_back(new void*[XSolid_Max_Parameter]);
				for (unsigned m = 0; m < XSolid_Max_Parameter; m++) this->objs_tp.back()[m] = 0;
				switch (unsigned(objs.at(i))) {
				case CSG_RoundingOp_Sign:;
				case CSG_OnionOp_Sign: {
					this->objs_tp.back()[0] = new double(*((double*)(other.objs_tp.at(i)[0])));
					break;
				}
				}
			}
		}
	}

	unsigned Max_StackLength; unsigned Stack_Usage[4]; double* *Stack; bool TIed; 	// 4 rendering threads + 1 attempting thread
	borderbox border;
public:
	unsigned type;	// smooth / diffuse / crystal
	rgblight col;	// or damping ratio
	double refractive_index;

	XSolid() {
		Stack = 0;
	}
	XSolid(const XObjs::XObjs_Comp &a) {
		objs.push_back(&a);
		objs_tp.push_back(0);
		border = a.MaxMin();
		Stack = 0;
	}
	XSolid(const XSolid &other) {
		for (int i = 0; i < other.objs.size(); i++) {
			objs.push_back(other.objs.at(i));
		}
		this->deep_copy(other);
		this->border = other.border;
		Stack = 0;
	}
	XSolid& operator = (const XSolid &other) {
		objs.clear(), objs_tp.clear();
		for (int i = 0; i < other.objs.size(); i++) {
			objs.push_back(other.objs.at(i));
		}
		this->deep_copy(other);
		this->border = other.border;
		Stack = 0;
		return *this;
	}
	XSolid& operator = (const XObjs::XObjs_Comp &other) {
		objs.clear(), objs_tp.clear();
		objs.push_back(&other), objs_tp.push_back(0);
		border = other.MaxMin();
		Stack = 0;
		return *this;
	}
	~XSolid() {
		objs.clear();
		for (unsigned i = 0, l = objs_tp.size(); i < l; i++) {
			if (objs_tp.at(i) != 0) {
				for (unsigned j = 0; j < XSolid_Max_Parameter; j++) {
					if (objs_tp.at(i)[j] == 0) break;
					else delete objs_tp.at(i)[j], objs_tp.at(i)[j] = 0;
				}
				delete objs_tp.at(i);
			}
		}
		if (Stack != 0) {
			delete Stack[0], Stack[1], Stack[2], Stack[3];
			delete Stack; Stack = 0;
		}
	}
	void setColor(rgblight c) {
		col = c;
	}
	point Max() const { return border.Max; }
	point Min() const { return border.Min; }
	unsigned telltype() const { return XSolid_Sign; }

	void init() {
		if (objs.empty()) border.Max = border.Min = point(0, 0, 0);
		border.fix();

		if (type != XSolid_Smooth && type != XSolid_Diffuse && type != XSolid_Crystal) type = XSolid_Smooth;
		if (refractive_index > 100 || refractive_index < 0.01) refractive_index = 1.5;

		Max_StackLength = 0;
		for (int i = 0, sz = 0; i < objs.size(); i++) {
			if ((unsigned)(objs.at(i)) >> 16) sz++;
			else sz -= ((unsigned)(objs.at(i)) >> 12) - 1;
			if (Max_StackLength < sz) Max_StackLength = sz;
		}

		if (Stack != 0) delete Stack[0], Stack[1], Stack[2], Stack[3], Stack;
		Stack = new double*[4]; TIed = false;
		for (int i = 0; i < 4; i++) Stack[i] = new double[Max_StackLength], Stack_Usage[i] = false;
	}

	double SDF(const point &P, const unsigned &_thread_No) const {
		double *S = Stack[_thread_No];
		for (unsigned i = 0, l = objs.size(), dir = 0; i < l; i++) {
			switch ((unsigned)(objs.at(i))) {
			case CSG_UnionOp_Sign: {
				dir--;
				if (S[dir - 1] > S[dir]) S[dir - 1] = S[dir];
				break;
			}
			case CSG_IntersectionOp_Sign: {
				dir--;
				if (S[dir - 1] < S[dir]) S[dir - 1] = S[dir];
				break;
			}
			case CSG_SubtractionOp_Sign: {
				dir--;
				if (S[dir - 1] < -S[dir]) S[dir - 1] = -S[dir];
				break;
			}
			case CSG_ComplementOp_Sign: {
				S[dir - 1] = -S[dir - 1];
				break;
			}
			case CSG_RoundingOp_Sign: {
				S[dir - 1] -= *((double*)objs_tp.at(i)[0]);		// doesn't work
				break;
			}
			case CSG_OnionOp_Sign: {
				S[dir - 1] = abs(S[dir - 1]) - *((double*)objs_tp.at(i)[0]);
				break;
			}
			default: {
				S[dir] = objs.at(i)->SDF(P); dir++;
			}
			}
		}
		return S[0];
	}

	// Note that the reflect vector is not calculated in this function
	void meet(intersect &R, const ray &a) const {
		unsigned thread_no = 0, k = *((unsigned*)(&this_thread::get_id()));
#define STACKUSAGE_XSOLID const_cast<unsigned*>(Stack_Usage)
		TIed ? (k == STACKUSAGE_XSOLID[0] ? thread_no = 0 : (k == STACKUSAGE_XSOLID[1] ? thread_no = 1 :
			(k == STACKUSAGE_XSOLID[2] ? thread_no = 2 : (k == STACKUSAGE_XSOLID[3] ? thread_no = 3 : (STACKUSAGE_XSOLID[0] = k) && (thread_no = 0)))))
			: (STACKUSAGE_XSOLID[0] ? (k == STACKUSAGE_XSOLID[0] ? thread_no = 0 : STACKUSAGE_XSOLID[1] ? (k == STACKUSAGE_XSOLID[1] ? thread_no = 0 :
				STACKUSAGE_XSOLID[2] ? (k == STACKUSAGE_XSOLID[2] ? thread_no = 0 : (STACKUSAGE_XSOLID[3] = k) && (thread_no = 3) && (*const_cast<bool*>(&TIed) = true))
				: (STACKUSAGE_XSOLID[2] = k) && (thread_no = 2)) : (STACKUSAGE_XSOLID[1] = k) && (thread_no = 1)) : (STACKUSAGE_XSOLID[0] = k) && (thread_no = 0));
		R.wt = thread_no;
		R.meet = 0;
		if (!border.meet(a)) return;
		double t = ERR_ZETA, sdf; point P;
		P = a.orig + t * a.dir, R.vt = sdf = this->SDF(P, thread_no);
		if (sdf > 0) {
			do {
				t += sdf;
				if (t > ERR_UPSILON) return;
				P = a.orig + t * a.dir, sdf = this->SDF(P, thread_no);
			} while (sdf > ERR_EPSILON);
			R.dist = t, R.intrs = P;
			R.meet = true;
			return;
		}
		else {
			do {
				t -= sdf;
				if (t > ERR_UPSILON) return;
				P = a.orig + t * a.dir, sdf = this->SDF(P, thread_no);
			} while (sdf < -ERR_EPSILON);
			R.dist = t, R.intrs = P;
			R.meet = true;
			return;
		}
	}

	void reflectData(intersect &R, const ray &a) const {
		point N;
		N.x = (this->SDF(point(R.intrs.x + ERR_ZETA, R.intrs.y, R.intrs.z), R.wt) - this->SDF(point(R.intrs.x - ERR_ZETA, R.intrs.y, R.intrs.z), R.wt)) * (0.5 / ERR_ZETA);
		N.y = (this->SDF(point(R.intrs.x, R.intrs.y + ERR_ZETA, R.intrs.z), R.wt) - this->SDF(point(R.intrs.x, R.intrs.y - ERR_ZETA, R.intrs.z), R.wt)) * (0.5 / ERR_ZETA);
		N.z = (this->SDF(point(R.intrs.x, R.intrs.y, R.intrs.z + ERR_ZETA), R.wt) - this->SDF(point(R.intrs.x, R.intrs.y, R.intrs.z - ERR_ZETA), R.wt)) * (0.5 / ERR_ZETA);
		R.reflect = a.dir - 2 * dot(N, a.dir) * N;	// now N is unit vector
		if (type == XSolid_Crystal) {
			// send R.ut as the refractive index of the other media
			// meet = air->obj ? 1 : 0; intrs = refract; reflect = reflect; vt = rate-of-reflection; 
			R.meet = R.vt > 0 ? 1 : 0;
			double c1 = dot(N, a.dir);
			double n1 = R.ut, n2 = refractive_index;
			if (R.vt < 0) n1 = refractive_index, n2 = R.ut;
			else c1 = -c1, N = -N;
			double c = n1 / n2, c2 = sqrt(1 - c * c * (1 - c1 * c1));
			if (isnan(c2)) {
				R.intrs.x = R.intrs.y = R.intrs.z = NAN; R.vt = 1;
				return;
			}
			R.intrs = c * a.dir - (c*c1 - c2)*N;
			double Rs = (n2*c1 - n1 * c2) / (n2*c1 + n1 * c2); Rs *= Rs;
			double Rp = (n2*c2 - n1 * c1) / (n2*c2 + n1 * c1); Rp *= Rp;
			R.vt = 0.5*(Rs + Rp);
		}
	}

	friend XSolid CSG_UnionOp(const XSolid &A, const XSolid &B) {
		XSolid X;
		for (int i = 0; i < A.objs.size(); i++) X.objs.push_back(A.objs.at(i));
		for (int i = 0; i < B.objs.size(); i++) X.objs.push_back(B.objs.at(i));
		X.deep_copy(A), X.deep_copy(B);
		X.objs.push_back(CSG_UnionOperator), X.objs_tp.push_back(0);
		X.border.Max = PMax(A.border.Max, B.border.Max), X.border.Min = PMin(A.border.Min, B.border.Min);
		return X;
	}
	friend XSolid CSG_IntersectionOp(const XSolid &A, const XSolid &B) {
		XSolid X;
		for (int i = 0; i < A.objs.size(); i++) X.objs.push_back(A.objs.at(i));
		for (int i = 0; i < B.objs.size(); i++) X.objs.push_back(B.objs.at(i));
		X.deep_copy(A), X.deep_copy(B);
		X.objs.push_back(CSG_IntersectionOperator), X.objs_tp.push_back(0);
		X.border.Max = PMin(A.border.Max, B.border.Max), X.border.Min = PMax(A.border.Min, B.border.Min);
		return X;
	}
	friend XSolid CSG_SubtractionOp(const XSolid &A, const XSolid &B) {
		XSolid X;
		for (int i = 0; i < A.objs.size(); i++) X.objs.push_back(A.objs.at(i));
		for (int i = 0; i < B.objs.size(); i++) X.objs.push_back(B.objs.at(i));
		X.deep_copy(A), X.deep_copy(B);
		X.objs.push_back(CSG_SubtractionOperator), X.objs_tp.push_back(0);
		X.border = A.border;	// need to optimize
		return X;
	}
	friend XSolid CSG_ComplementOp(const XSolid &A) {
		XSolid X;
		for (int i = 0; i < A.objs.size(); i++) X.objs.push_back(A.objs.at(i));
		X.deep_copy(A);
		X.objs.push_back(CSG_ComplementOperator), X.objs_tp.push_back(0);
		X.border.Max = point(INFINITY, INFINITY, INFINITY), X.border.Min = -X.border.Max;
		return X;
	}
	friend XSolid CSG_RoundingOp(const XSolid &A, double r) {
		XSolid X;
		for (int i = 0; i < A.objs.size(); i++) X.objs.push_back(A.objs.at(i));
		X.deep_copy(A);
		X.objs.push_back(CSG_RoundingOperator), X.objs_tp_push();
		X.objs_tp.back()[0] = new double(r);
		X.border.Max = A.border.Max + point(r, r, r), X.border.Min = A.border.Min - point(r, r, r);
		return X;
	}
	friend XSolid CSG_OnionOp(const XSolid &A, double r) {
		XSolid X;
		for (int i = 0; i < A.objs.size(); i++) X.objs.push_back(A.objs.at(i));
		X.deep_copy(A);
		X.objs.push_back(CSG_OnionOperator), X.objs_tp_push();
		X.objs_tp.back()[0] = new double(r);
		X.border.Max = A.border.Max + point(r, r, r), X.border.Min = A.border.Min - point(r, r, r);
		return X;
	}
};

