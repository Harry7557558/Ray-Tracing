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

	class Box : public XObjs_Comp {
		point O; point A, B, C; // parallelepiped
	public:

	};
}

// Windows null pointer: 0x00000000-0x0000FFFF
#define CSG_UnionOp_Sign 0x0201
#define CSG_IntersectionOp_Sign 0x0202
#define CSG_SubtractionOp_Sign 0x0203

#define CSG_UnionOperator ((const XObjs::XObjs_Comp*)CSG_UnionOp_Sign)
#define CSG_IntersectionOperator ((const XObjs::XObjs_Comp*)CSG_IntersectionOp_Sign)
#define CSG_SubtractionOperator ((const XObjs::XObjs_Comp*)CSG_SubtractionOp_Sign)

#include <thread>
#include <mutex>

#define XSolid_Sign 0x10000001
class XSolid : public object {
	vector<const XObjs::XObjs_Comp*> objs;	// postfix expression
	unsigned Max_StackLength; unsigned Stack_Usage[4]; double* *Stack; bool TIed; mutex _MUTEX;	// 4 rendering threads
	borderbox border;
	XSolid() {}
public:
	rgblight col;
	XSolid(const XObjs::XObjs_Comp &a) {
		objs.push_back(&a);
		border = a.MaxMin();
		Stack = 0;
	}
	XSolid(const XSolid &other) {
		for (int i = 0; i < other.objs.size(); i++) {
			objs.push_back(other.objs.at(i));
		}
		this->border = other.border;
		Stack = 0;
	}
	XSolid& operator = (const XSolid &other) {
		objs.clear();
		for (int i = 0; i < other.objs.size(); i++) {
			objs.push_back(other.objs.at(i));
		}
		this->border = other.border;
		Stack = 0;
		return *this;
	}
	~XSolid() {
		objs.clear();
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
		border.fix();
		int sz = 0; Max_StackLength = 0;
		for (int i = 0; i < objs.size(); i++) {
			if ((unsigned)(objs.at(i)) >> 16) sz++;
			else sz -= ((unsigned)(objs.at(i)) >> 8) - 1;
			if (Max_StackLength < sz) Max_StackLength = sz;
		}
		Stack = new double*[4]; TIed = false;
		for (int i = 0; i < 4; i++) Stack[i] = new double[Max_StackLength], Stack_Usage[i] = false;
	}

	double SDF(const point &P, const unsigned &_thread_No) const {
		double *S = Stack[_thread_No];
		double s;
		for (int i = 0, dir = 0; i < objs.size(); i++) {
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
			default: {
				S[dir] = objs.at(i)->SDF(P); dir++;
			}
			}
		}
		s = S[0];
		return s;
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
		double t = 0, sdf; point P;
		if (this->SDF(a.orig, thread_no) > ERR_EPSILON) {
			P = a.orig, sdf = this->SDF(P, thread_no);
			double n = 0; do {
				t += sdf;
				if (t > ERR_UPSILON) return;
				P = a.orig + t * a.dir, sdf = this->SDF(P, thread_no);
			} while (abs(sdf) > ERR_EPSILON);
			R.dist = t, R.intrs = P;
			R.meet = true;
			//R.reflect = point(0, 0, 1);
			//fout << ray(a.orig, P - a.orig) << endl << R.dist << endl;
			return;
		}
		else return;
	}

	void reflectData(intersect &R, const ray &a) const {
		point N;
		N.x = (this->SDF(point(R.intrs.x + ERR_ZETA, R.intrs.y, R.intrs.z), R.wt) - this->SDF(point(R.intrs.x - ERR_ZETA, R.intrs.y, R.intrs.z), R.wt)) * (0.5 / ERR_ZETA);
		N.y = (this->SDF(point(R.intrs.x, R.intrs.y + ERR_ZETA, R.intrs.z), R.wt) - this->SDF(point(R.intrs.x, R.intrs.y - ERR_ZETA, R.intrs.z), R.wt)) * (0.5 / ERR_ZETA);
		N.z = (this->SDF(point(R.intrs.x, R.intrs.y, R.intrs.z + ERR_ZETA), R.wt) - this->SDF(point(R.intrs.x, R.intrs.y, R.intrs.z - ERR_ZETA), R.wt)) * (0.5 / ERR_ZETA);
		//fout << N << endl;
		R.reflect = a.dir - 2 * dot(N, a.dir) * N;
		//fout << ray(R.intrs, R.reflect) << endl;
	}

	friend XSolid CSG_UnionOp(const XSolid &A, const XSolid &B) {
		XSolid X;
		for (int i = 0; i < A.objs.size(); i++) X.objs.push_back(A.objs.at(i));
		for (int i = 0; i < B.objs.size(); i++) X.objs.push_back(B.objs.at(i));
		X.objs.push_back(CSG_UnionOperator);
		X.border.Max = PMax(A.border.Max, B.border.Max), X.border.Min = PMin(A.border.Min, B.border.Min);
		return X;
	}
	friend XSolid CSG_IntersectionOp(const XSolid &A, const XSolid &B) {
		XSolid X;
		for (int i = 0; i < A.objs.size(); i++) X.objs.push_back(A.objs.at(i));
		for (int i = 0; i < B.objs.size(); i++) X.objs.push_back(B.objs.at(i));
		X.objs.push_back(CSG_IntersectionOperator);
		X.border.Max = PMin(A.border.Max, B.border.Max), X.border.Min = PMax(A.border.Min, B.border.Min);
		return X;
	}
	friend XSolid CSG_SubtractionOp(const XSolid &A, const XSolid &B) {
		XSolid X;
		for (int i = 0; i < A.objs.size(); i++) X.objs.push_back(A.objs.at(i));
		for (int i = 0; i < B.objs.size(); i++) X.objs.push_back(B.objs.at(i));
		X.objs.push_back(CSG_SubtractionOperator);
		X.border = A.border;	// need to optimize
		return X;
	}
};

