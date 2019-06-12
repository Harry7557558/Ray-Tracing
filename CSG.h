#pragma once

#ifndef Object_Sign
#include "Object.h"
#endif

// Reference: https://www.iquilezles.org/www/articles/distfunctions/distfunctions.htm

// All object/SDF/operation/transformations labled "exact" or "not exact" are strictly prooved or disproved
// object exact: legal shape (not an approximation);  SDF exact: legal distance field (issues won't occur when applying operations);
// "not exact" operations/transformations may affect farther operations/transformations

class CatmullRom_2D {
	vector<point2D> points;
	bool close;
public:
#pragma warning(disable: 4010)
	// Create a smooth curve between P1 and P2					\
	                |-0.5  1.5 -1.5  0.5 | |P0|					\
	P = [t³ t² t 1] | 1.0 -2.5  2.0 -0.5 | |P1|   0 ≤ t ≤ 1		\
	                |-0.5  0.0  0.5  0.0 | |P2|					\
	                | 0.0  1.0  0.0  0.0 | |P3|					\
	
	// SDF of P(p0,p1) and {x=x(t),y=y(t)}: derivative p0*x'(t) + p1*y'(t) + x'(t)*x(t) + y'(t)*y(t)

	CatmullRom_2D() : close(false) {}
	CatmullRom_2D(initializer_list<point2D> q) : close(false) { points.assign(q); }
	~CatmullRom_2D() { points.clear(); close = false; }

	void assign(point2D P) { points.push_back(P); close = false; }
	void assign(initializer_list<point2D> q) { points.assign(q); close = false; }
	void construct_vertex() {
		points.insert(points.begin(), 2 * points[0] - points[1]);
		points.push_back(2 * points.back() - points[points.size() - 2]);
	}
	void close_path() {
		points.push_back(points[0]); points.push_back(points[1]);
		close = true;
	}
	point2D& at(unsigned i) { close = false; return points.at(i); }
	void clear() { points.clear(); close = false; }

	point2D eval(double t) {
		unsigned n = t; t -= n;
		return (((-0.5*t + 1.0)*t - 0.5)*t) * points[n - 1]
			+ (((1.5*t - 2.5)*t)*t + 1.0) * points[n]
			+ (((-1.5*t + 2.0)*t + 0.5)*t) * points[n + 1]
			+ (((0.5*t - 0.5)*t)*t) * points[n + 2];
	}

};

double QuadraticBezier_SDF(const point2D &P, const point2D &A, const point2D &B, const point2D &C, const point2D &E, const point2D &F)
{
	double a, b, c, d;
	a = dot(E, E), b = 3 * dot(E, F), c = 2 * dot(F, F) + dot(E, A - P), d = dot(F, A - P);
	double r, u, v;
	auto eval = [](double t, const point2D &A, const point2D &B, const point2D &C) -> point2D {
		return (1 - t)*(1 - t)*A + 2 * t*(1 - t)*B + t * t*C;
	};
	if (solveCubic(a, b, c, d, r, u, v)) {
		if (u < r) swap(u, r); if (r < v) swap(r, v); if (u < r) swap(u, r);	// r is the maxima and u,v are minimas
		if ((u < 0 || u > 1) && (v < 0 || v > 1)) return min((P - A).mod(), (P - C).mod());
		if ((u > 0 && u < 1) && (v > 0 && v < 1)) return min((P - eval(u, A, B, C)).mod(), (P - eval(v, A, B, C)).mod());
		if (u < 0 || u > 1) return min(min((P - A).mod(), (P - C).mod()), (P - eval(v, A, B, C)).mod());
		if (v < 0 || v > 1) return min(min((P - A).mod(), (P - C).mod()), (P - eval(u, A, B, C)).mod());
	}
	else {
		if (r < 0 || r > 1) return min((P - A).mod(), (P - C).mod());
		return (P - eval(r, A, B, C)).mod();
	}
	return NAN;
}
inline double QuadraticBezier_SDF(point2D P, point2D A, point2D B, point2D C)
{
	return QuadraticBezier_SDF(P, A, B, C, C + A - 2 * B, B - A);
}



namespace XObjs {
#define XObjs_Sphere_Sign 0x00000001
#define XObjs_Plane_Sign 0x00000002
#define XObjs_CylinderSTD_Sign 0x00000003
#define XObjs_ConeSTD_Sign 0x00000004
#define XObjs_TorusXOY_Sign 0x00000005
#define XObjs_BoxXOY_Sign 0x00000006
#define XObjs_Cylinder_Sign 0x00000007
#define XObjs_Cone_Sign 0x00000008
#define XObjs_ConeTrunc_Sign 0x00000009
#define XObjs_BoxAffine_Sign 0x00000101
#define XObjs_BezierQuadraticXOY_Sign 0x0000000A

	class XObjs_Comp {
	public:
		XObjs_Comp() {}
		virtual XObjs_Comp* clone() const { return new XObjs_Comp(*this); }		// "virtual constructor"
		virtual ~XObjs_Comp() {}
		virtual double SDF(const point &P) const {
			WARN("XObjs::XObjs_Comp.SDF() called!");
			return NAN;
		}
		virtual borderbox MaxMin() const {
			WARN("XObjs::XObjs_Comp.SDF() called!");
			return borderbox(point(NAN, NAN, NAN), point(NAN, NAN, NAN));
		}
		virtual unsigned telltype() { return 0; }
	};


	class Sphere : public XObjs_Comp {	// exact
		point C; double r;
	public:
		Sphere() :r(0) {}
		Sphere(point C, double r) {
			this->C = C, this->r = r;
		}
		Sphere(const Sphere &other) {
			C = other.C, r = other.r;
		}
		virtual Sphere* clone() const { return new Sphere(*this); }
		~Sphere() {}

		double SDF(const point &P) const {	// exact
			return (P - C).mod() - r;
		}
		borderbox MaxMin() const {
			return borderbox(C - point(r, r, r), C + point(r, r, r));
		}
		virtual unsigned telltype() { return XObjs_Sphere_Sign; }
	};

	class Plane : public XObjs_Comp {	// exact
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
		Plane* clone() const { return new Plane(*this); }
		~Plane() {}

		double SDF(const point &P) const {	// exact
			return dot(N, P) - D;
		}
		borderbox MaxMin() const {
			return borderbox(point(-INFINITY, -INFINITY, -INFINITY), point(INFINITY, INFINITY, INFINITY));
		}
		virtual unsigned telltype() { return XObjs_Plane_Sign; }
	};

	class Cylinder_std : public XObjs_Comp {	// exact
		point C; double r; point D;		// central point, radius, direction vector (unit)
	public:
		Cylinder_std() { r = 0; }
		Cylinder_std(const point &C, const double &r, const  vec3 &d) {
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
		Cylinder_std* clone() const { return new Cylinder_std(*this); }
		~Cylinder_std() {}

		double SDF(const point &P) const {	// exact
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
		virtual unsigned telltype() { return XObjs_CylinderSTD_Sign; }
	};

	class Cone_std : public XObjs_Comp {	// exact
		point C; double alpha; point D;		// 0 < alpha < PI/2, D unit
	public:
		Cone_std() { alpha = 0; }
		Cone_std(const point &C, const point &P, const double &r) {
			this->C = C, this->D = P - C;
			double m = D.mod();
			alpha = atan(r / m);
			D /= m;
		}
		Cone_std(const point &C, const double &alpha, const point &D) {
			this->C = C, this->alpha = alpha, this->D = D / D.mod();
		}
		Cone_std(const Cone_std &other) {
			C = other.C, alpha = other.alpha, D = other.D;
		}
		Cone_std* clone() const { return new Cone_std(*this); }
		~Cone_std() {}

		double SDF(const point &P) const {	// exact
			point S = P - C; double m = S.mod(), theta = acos(dot(S, D) / m) - alpha;
			if (theta < PI / 2) return m * sin(theta);
			else return m;
		}
		borderbox MaxMin() const {
			borderbox B;
			B.Max = point(INFINITY, INFINITY, INFINITY); B.Min = -B.Max;	// Optimize
			return B;
		}
		virtual unsigned telltype() { return XObjs_ConeSTD_Sign; }
	};

	class Torus_xOy : public XObjs_Comp {	// exact
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
		Torus_xOy* clone() const { return new Torus_xOy(*this); }
		~Torus_xOy() {}

		double SDF(const point &P) const {	// exact
			//return hypot(hypot(P.x - C.x, P.y - C.y) - R, P.z - C.z) - r;
			point S = P - C; double m = sqrt(S.x*S.x + S.y*S.y) - R; m = sqrt(m*m + S.z*S.z); return m - r;
		}
		borderbox MaxMin() const {
			borderbox B;
			B.Max = point(R + r, R + r, r), B.Min = -B.Max;
			B.Max += C, B.Min += C;
			return B;
		}
		virtual unsigned telltype() { return XObjs_TorusXOY_Sign; }
	};

	class Box_xOy : public XObjs_Comp {		// exact
		point Min; point Max;
	public:
		Box_xOy() {}
		Box_xOy(const borderbox &B) {
			this->Max = B.Max, this->Min = B.Min;
		}
		Box_xOy(const point &Max, const point &Min) {
			this->Max = Max; this->Min = Min;
		}
		Box_xOy(const point &Ctr, const double &x_rad, const double &y_rad, const double &z_rad) {
			this->Max = point(x_rad, y_rad, z_rad), this->Min = -this->Max;
			this->Max += Ctr, this->Min += Ctr;
		}
		Box_xOy(const Box_xOy &other) {
			this->Max = other.Max, this->Min = other.Min;
		}
		Box_xOy* clone() const { return new Box_xOy(*this); }
		~Box_xOy() {}

		double SDF(const point &P) const {	// exact
			unsigned char xi = (P.x > Min.x && P.x < Max.x), yi = (P.y > Min.y && P.y < Max.y), zi = (P.z > Min.z && P.z < Max.z);
			switch (xi + yi + zi) {
			case 3: {
				return max(max(max(P.x - Max.x, Min.x - P.x), max(P.y - Max.y, Min.y - P.y)), max(P.z - Max.z, Min.z - P.z));
			}
			case 2: {
				if (!xi) return P.x > Max.x ? P.x - Max.x : Min.x - P.x;
				if (!yi) return P.y > Max.y ? P.y - Max.y : Min.y - P.y;
				if (!zi) return P.z > Max.z ? P.z - Max.z : Min.z - P.z;
			}
			case 1: {
				double d1, d2;
				if (xi) { d1 = min(abs(P.y - Max.y), abs(P.y - Min.y)), d2 = min(abs(P.z - Max.z), abs(P.z - Min.z)); return sqrt(d1 * d1 + d2 * d2); }
				if (yi) { d1 = min(abs(P.x - Max.x), abs(P.x - Min.x)), d2 = min(abs(P.z - Max.z), abs(P.z - Min.z)); return sqrt(d1 * d1 + d2 * d2); }
				if (zi) { d1 = min(abs(P.x - Max.x), abs(P.x - Min.x)), d2 = min(abs(P.y - Max.y), abs(P.y - Min.y)); return sqrt(d1 * d1 + d2 * d2); }
			}
			case 0: {
				double d1, d2, d3;
				d1 = min(abs(P.x - Max.x), abs(P.x - Min.x));
				d2 = min(abs(P.y - Max.y), abs(P.y - Min.y));
				d3 = min(abs(P.z - Max.z), abs(P.z - Min.z));
				return sqrt(d1 * d1 + d2 * d2 + d3 * d3);
			}
			}
			return NAN;
		}
		borderbox MaxMin() const {
			return borderbox(Max, Min);
		}
		virtual unsigned telltype() { return XObjs_BoxXOY_Sign; }
	};

	class Cylinder : public XObjs_Comp {	// exact
		point P1, P2; double r, h; vec3 dir;
	public:
		Cylinder() : r(0), h(0) {}
		Cylinder(const point &P1, const point &P2, const double &r) {
			this->P1 = P1, this->P2 = P2, this->r = r;
			dir = P2 - P1; h = dir.mod(), dir /= h;
		}
		Cylinder(const point &P1, const double &r, const vec3 &dir) {
			this->P1 = P1, this->P2 = P1 + dir, this->r = r;
			this->dir = dir; h = dir.mod(), this->dir /= h;
		}
		Cylinder(const point &P1, const double &r, const double &h, const vec3 &dir) {
			this->dir = dir / dir.mod(); this->h = h;
			this->P1 = P1, this->P2 = P1 + h * this->dir, this->r = r;
		}
		Cylinder(const Cylinder &other) {
			this->P1 = other.P1, this->P2 = other.P2, this->r = other.r;
			this->dir = other.dir, this->h = other.h;
		}
		Cylinder* clone() const { return new Cylinder(*this); }
		~Cylinder() {}

		double SDF(const point &P) const {	// exact
			point S = P - P1;
			double n = dot(S, dir);
			double d = cross(S, dir).mod() - r;
			if (n < 0) {
				if (d > 0) return sqrt(d * d + n * n);
				return -n;
			}
			if (n > h) {
				n -= h;
				if (d > 0) return sqrt(d * d + n * n);
				return n;
			}
			if (d > 0) return d;
			return max(max(-n, n - h), d);
		}
		borderbox MaxMin() const {
			vec3 u = P2 - P1; u /= u.mod();
			u.x = sqrt(1 - u.x*u.x), u.y = sqrt(1 - u.y*u.y), u.z = sqrt(1 - u.z*u.z); u *= r;
			borderbox B(P1, P2); B.fix();
			B.Max += u, B.Min -= u;
			return B;
		}
		virtual unsigned telltype() { return XObjs_Cylinder_Sign; }
	};

	class Cone : public XObjs_Comp {	// exact
		point C, O; vec3 dir; double r, h, l, alpha;
	public:
		Cone() { r = h = 0; }
		Cone(const point &C, const point &P, const double &r) {
			this->C = C, this->O = P, this->r = r, this->dir = P - C;
			h = dir.mod(); dir /= h;
			alpha = atan(r / h), l = r / sin(alpha);
		}
		Cone(const point &C, const vec3 &dir, const double &R, const double &theta) {
			this->C = C, this->dir = dir / dir.mod();
			this->r = R * sin(theta), this->h = R * cos(theta), this->l = R;
			this->O = C + dir * this->h;
			alpha = theta;
		}
		Cone(const Cone &other) {
			this->C = other.C, this->O = other.O; this->dir = other.dir;
			this->r = other.r, this->h = other.h, this->l = other.l; this->alpha = other.alpha;
		}
		Cone* clone() const { return new Cone(*this); }
		~Cone() {}

		double SDF(const point &P) const {	// exact
			point S = P - C; double m = S.mod(), u = dot(S, dir), theta = acos(u / m) - alpha;
			if (theta > PI / 2) return m;
			if (u > h) {
				u -= h;
				double v = m * sin(theta + alpha) - r;
				if (v > 0) return sqrt(u * u + v * v);
				return u;
			}
			if (theta > 0) {
				if (m * cos(theta) > l) return sqrt(m*m + l * l - 2 * m*l*cos(theta));
				return m * sin(theta);
			}
			return max(m*sin(theta), u - h);
		}
		borderbox MaxMin() const {
			point ad;
			double rx = acos(dir.z), rz = atan2(dir.x, -dir.y);
			double u = -atan(cos(rx)*tan(rz));
			ad.x = cos(rz)*cos(u) - cos(rx)*sin(rz)*sin(u);
			u = atan(cos(rx) / tan(rz));
			ad.y = sin(rz)*cos(u) + cos(rx)*cos(rz)*sin(u);
			ad.z = sin(rx);
			ad *= r;
			return borderbox(PMax(PMax(O + ad, O - ad), C), PMin(PMin(O + ad, O - ad), C));
		}
		virtual unsigned telltype() { return XObjs_Cone_Sign; }
	};

	class Cone_trunc : public XObjs_Comp {		// exact
		point C, P1, P2; double r1, r2, h1, h2, l1, l2, alpha; vec3 dir;
	public:
		Cone_trunc() { r1 = r2 = h1 = h2 = l1 = l2 = alpha = 0; }
		Cone_trunc(const point &P_1, const point &P_2, const double &r_1, const double &r_2) {
			P1 = P_1, P2 = P_2; r1 = r_1, r2 = r_2; if (r1 > r2) swap(r1, r2), swap(P1, P2);
			dir = P2 - P1; double h = dir.mod(); dir /= h; alpha = atan((r2 - r1) / h);
			l1 = r1 / sin(alpha), l2 = r2 / sin(alpha); h1 = r1 / tan(alpha), h2 = r2 / tan(alpha);
			C = P1 - h1 * dir;
		}
		Cone_trunc(const point &C, const vec3 &dir, const double &h_1, const double &h_2, const double &alpha) {
			this->C = C, this->dir = dir / dir.mod(); h1 = h_1, h2 = h_2, this->alpha = alpha;
			P1 = C + h1 * dir, P2 = C + h2 * dir; r1 = h1 * tan(alpha), r2 = h2 * tan(alpha), l1 = h1 / cos(alpha), l2 = h2 / cos(alpha);
		}
		Cone_trunc(const Cone_trunc &other) {
			C = other.C, P1 = other.P1, P2 = other.P2; dir = other.dir;
			r1 = other.r1, r2 = other.r2, h1 = other.h1, h2 = other.h2, l1 = other.l1, l2 = other.l2, alpha = other.alpha;
		}
		Cone_trunc* clone() const { return new Cone_trunc(*this); }
		~Cone_trunc() {}

		double SDF(const point &P) const {	// exact
			point s = P - C; double m = s.mod(), u = dot(s, dir), theta = acos(u / m), d = m * sin(theta); theta -= alpha;
			if (theta > PI / 2) return d < r1 ? h1 - u : sqrt(m*m + l1 * l1 - 2 * m*l1*cos(theta));
			if (u > h1 && u < h2 && theta < 0) return max(max(h1 - u, u - h2), m*sin(theta));
			if (u > h2) return d < r2 ? u - h2 : sqrt(m*m + l2 * l2 - 2 * m*l2*cos(theta));
			if (u < h1) return d < r1 ? h1 - u : m * cos(theta) > l1 ? m * sin(theta) : sqrt(m*m + l1 * l1 - 2 * m*l1*cos(theta));
			return m * cos(theta) < l2 ? m * sin(theta) : sqrt(m*m + l2 * l2 - 2 * m*l2*cos(theta));
		}
		borderbox MaxMin() const {
			point ad1, ad2;
			double rx = acos(dir.z), rz = atan2(dir.x, -dir.y);
			double u = -atan(cos(rx)*tan(rz));
			ad1.x = cos(rz)*cos(u) - cos(rx)*sin(rz)*sin(u);
			u = atan(cos(rx) / tan(rz));
			ad1.y = sin(rz)*cos(u) + cos(rx)*cos(rz)*sin(u);
			ad1.z = sin(rx);
			ad2 = ad1 * r2, ad1 *= r1;
			return borderbox(PMax(PMax(P1 + ad1, P1 - ad1), PMax(P2 + ad2, P2 - ad2)), PMin(PMin(P1 + ad1, P1 - ad1), PMin(P2 + ad2, P2 - ad2)));
		}
		virtual unsigned telltype() { return XObjs_ConeTrunc_Sign; }
	};

	class Box_affine : public XObjs_Comp {
		matrix3D_affine M;	// applying to an unit cube, the only variable
		matrix3D_affine M_invert;
		point P000, P100, P010, P001, P110, P101, P011, P111;	// points
		vec3 Ex00, Ex01, Ex10, Ex11, Ey00, Ey01, Ey10, Ey11, Ez00, Ez01, Ez10, Ez11;	// unit direction vector of edges
		vec3 _Ex00, _Ex01, _Ex10, _Ex11, _Ey00, _Ey01, _Ey10, _Ey11, _Ez00, _Ez01, _Ez10, _Ez11;	// negative of edges
		double ex00, ex01, ex10, ex11, ey00, ey01, ey10, ey11, ez00, ez01, ez10, ez11;	// length of edges
		vec3 Px0, Px1, Py0, Py1, Pz0, Pz1; double px0, px1, py0, py1, pz0, pz1;		// plane data (including normal) of faces

		vec3 PET_Ex00_n1, PET_Ex00_n2, PET_Ex01_n1, PET_Ex01_n2, PET_Ex10_n1, PET_Ex10_n2, PET_Ex11_n1, PET_Ex11_n2,
			PET_Ey00_n1, PET_Ey00_n2, PET_Ey01_n1, PET_Ey01_n2, PET_Ey10_n1, PET_Ey10_n2, PET_Ey11_n1, PET_Ey11_n2,
			PET_Ez00_n1, PET_Ez00_n2, PET_Ez01_n1, PET_Ez01_n2, PET_Ez10_n1, PET_Ez10_n2, PET_Ez11_n1, PET_Ez11_n2;		// Point-Edge Test

		borderbox border;

		inline bool Point_Vertex_Test(const point &P, const vec3 &n1, const vec3 &n2, const vec3 &n3, double &dist) const {
			if (dot(P, n1) > -ERR_EPSILON && dot(P, n2) > -ERR_EPSILON && dot(P, n3) > -ERR_EPSILON) dist = P.mod();
			else return false;
			return true;
		}
		inline bool Point_Edge_Test(const point &P, const point &V1, const point &V2, const vec3 &dir,
			const vec3 &n1, const vec3 &n2, const vec3 &_n1, const vec3 &_n2, double &dist) const {
			if (dot(P - V1, n1) > -ERR_EPSILON && dot(P - V2, n2) > -ERR_EPSILON &&
				dot(P - V1, _n1) < ERR_EPSILON && dot(P - V2, _n2) < ERR_EPSILON) dist = cross(P - V1, dir).mod();
			else return false;
			return true;
		}
		inline bool Point_Plane_Test(const point &P, const point &A, const point &B, const point &C, const point &D, const vec3 &N, point &I, double &dist) const {
			dist = dot(P - A, N); if (dist < -ERR_EPSILON) return false;
			I = P - dist * N;
			vec3 Cp = cross(B - A, C - A), Cs; Cp /= Cp.mod();
			Cs = cross(B - A, I - A); if (abs(dot(Cp, Cs) / Cs.mod() - 1) > ERR_EPSILON) return false;
			Cs = cross(C - B, I - B); if (abs(dot(Cp, Cs) / Cs.mod() - 1) > ERR_EPSILON) return false;
			Cs = cross(D - C, I - C); if (abs(dot(Cp, Cs) / Cs.mod() - 1) > ERR_EPSILON) return false;
			Cs = cross(A - D, I - D); if (abs(dot(Cp, Cs) / Cs.mod() - 1) > ERR_EPSILON) return false;
			return true;
		}

	public:
		Box_affine() { M = matrix3D_affine(1, 0, 0, 0, 1, 0, 0, 0, 1); }
		Box_affine(const matrix3D_affine &M) {
			this->M = M;
		}
		Box_affine(const Box_affine& other) {
			this->M = other.M;
			this->MaxMin();
			//fout << M << endl;
		}
		Box_affine* clone() const { return new Box_affine(*this); }
		~Box_affine() {}

		inline bool meet_box(const ray &a) {
			return this->border.meet(a);
		}

		double SDF(const point &P) const {	// debugging
			point S = M_invert * P;
			if (S.x > 0 && S.x < 1 && S.y > 0 && S.y < 1 && S.z > 0 && S.z < 1) return -min(min(min(abs(dot(P, Px0) + px0), abs(dot(P, Px1) + px1)),
				min(abs(dot(P, Py0) + py0), abs(dot(P, Py1) + py1))), min(abs(dot(P, Pz0) + pz0), abs(dot(P, Pz1) + pz1)));

			double dist;

			if (Point_Vertex_Test(P - P000, _Ex00, _Ey00, _Ez00, dist)) return dist;
			if (Point_Vertex_Test(P - P001, _Ex01, _Ey01, Ez00, dist)) return dist;
			if (Point_Vertex_Test(P - P010, _Ex10, Ey00, _Ez01, dist)) return dist;
			if (Point_Vertex_Test(P - P011, _Ex11, Ey01, Ez01, dist)) return dist;
			if (Point_Vertex_Test(P - P100, Ex00, _Ey10, _Ez10, dist)) return dist;
			if (Point_Vertex_Test(P - P101, Ex01, _Ey11, Ez10, dist)) return dist;
			if (Point_Vertex_Test(P - P110, Ex10, Ey10, _Ez11, dist)) return dist;
			if (Point_Vertex_Test(P - P111, Ex11, Ey11, Ez11, dist)) return dist;

			if (Point_Edge_Test(P, P000, P100, Ex00, PET_Ex00_n1, PET_Ex00_n2, _Ex00, Ex00, dist)) return dist;
			if (Point_Edge_Test(P, P001, P101, Ex01, PET_Ex01_n1, PET_Ex01_n2, _Ex01, Ex01, dist)) return dist;
			if (Point_Edge_Test(P, P010, P110, Ex10, PET_Ex10_n1, PET_Ex10_n2, _Ex10, Ex10, dist)) return dist;
			if (Point_Edge_Test(P, P011, P111, Ex11, PET_Ex11_n1, PET_Ex11_n2, _Ex11, Ex11, dist)) return dist;
			if (Point_Edge_Test(P, P000, P010, Ey00, PET_Ey00_n1, PET_Ey00_n2, _Ey00, Ey00, dist)) return dist;
			if (Point_Edge_Test(P, P001, P011, Ey01, PET_Ey01_n1, PET_Ey01_n2, _Ey01, Ey01, dist)) return dist;
			if (Point_Edge_Test(P, P100, P110, Ey10, PET_Ey10_n1, PET_Ey10_n2, _Ey10, Ey10, dist)) return dist;
			if (Point_Edge_Test(P, P101, P111, Ey11, PET_Ey11_n1, PET_Ey11_n2, _Ey11, Ey11, dist)) return dist;
			if (Point_Edge_Test(P, P000, P001, Ez00, PET_Ez00_n1, PET_Ez00_n2, _Ez00, Ez00, dist)) return dist;
			if (Point_Edge_Test(P, P010, P011, Ez01, PET_Ez01_n1, PET_Ez01_n2, _Ez01, Ez01, dist)) return dist;
			if (Point_Edge_Test(P, P100, P101, Ez10, PET_Ez10_n1, PET_Ez10_n2, _Ez10, Ez10, dist)) return dist;
			if (Point_Edge_Test(P, P110, P111, Ez11, PET_Ez11_n1, PET_Ez11_n2, _Ez11, Ez11, dist)) return dist;

			if (Point_Plane_Test(P, P000, P001, P011, P010, Px0, S, dist)) return dist;
			if (Point_Plane_Test(P, P100, P101, P111, P110, Px1, S, dist)) return dist;
			if (Point_Plane_Test(P, P000, P100, P101, P001, Py0, S, dist)) return dist;
			if (Point_Plane_Test(P, P010, P110, P111, P011, Py1, S, dist)) return dist;
			if (Point_Plane_Test(P, P000, P100, P110, P010, Pz0, S, dist)) return dist;
			if (Point_Plane_Test(P, P001, P101, P111, P011, Pz1, S, dist)) return dist;

			return ERR_UPSILON;
		}

		borderbox MaxMin() const {
			// Since this function will be called before rendering, it would be okay to invert the matrix there.
			*const_cast<matrix3D_affine*>(&M_invert) = M.invert();
			*const_cast<point*>(&P000) = M * point(0, 0, 0), *const_cast<point*>(&P001) = M * point(0, 0, 1), *const_cast<point*>(&P100) = M * point(1, 0, 0), *const_cast<point*>(&P010) = M * point(0, 1, 0), *const_cast<point*>(&P110) = M * point(1, 1, 0), *const_cast<point*>(&P101) = M * point(1, 0, 1), *const_cast<point*>(&P011) = M * point(0, 1, 1), *const_cast<point*>(&P111) = M * point(1, 1, 1);
			*const_cast<vec3*>(&Ex00) = P100 - P000, *const_cast<vec3*>(&Ex01) = P101 - P001, *const_cast<vec3*>(&Ex10) = P110 - P010, *const_cast<vec3*>(&Ex11) = P111 - P011, *const_cast<vec3*>(&Ey00) = P010 - P000, *const_cast<vec3*>(&Ey01) = P011 - P001, *const_cast<vec3*>(&Ey10) = P110 - P100, *const_cast<vec3*>(&Ey11) = P111 - P101, *const_cast<vec3*>(&Ez00) = P001 - P000, *const_cast<vec3*>(&Ez01) = P011 - P010, *const_cast<vec3*>(&Ez10) = P101 - P100, *const_cast<vec3*>(&Ez11) = P111 - P110;
			*const_cast<double*>(&ex00) = Ex00.mod(), *const_cast<double*>(&ex01) = Ex01.mod(), *const_cast<double*>(&ex10) = Ex10.mod(), *const_cast<double*>(&ex11) = Ex11.mod(), *const_cast<double*>(&ey00) = Ey00.mod(), *const_cast<double*>(&ey01) = Ey01.mod(), *const_cast<double*>(&ey10) = Ey10.mod(), *const_cast<double*>(&ey11) = Ey11.mod(), *const_cast<double*>(&ez00) = Ez00.mod(), *const_cast<double*>(&ez01) = Ez01.mod(), *const_cast<double*>(&ez10) = Ez10.mod(), *const_cast<double*>(&ez11) = Ez11.mod();
			*const_cast<vec3*>(&Ex00) /= ex00, *const_cast<vec3*>(&Ex01) /= ex01, *const_cast<vec3*>(&Ex10) /= ex10, *const_cast<vec3*>(&Ex11) /= ex11, *const_cast<vec3*>(&Ey00) /= ey00, *const_cast<vec3*>(&Ey01) /= ey01, *const_cast<vec3*>(&Ey10) /= ey10, *const_cast<vec3*>(&Ey11) /= ey11, *const_cast<vec3*>(&Ez00) /= ez00, *const_cast<vec3*>(&Ez01) /= ez01, *const_cast<vec3*>(&Ez10) /= ez10, *const_cast<vec3*>(&Ez11) /= ez11;
			*const_cast<vec3*>(&Px0) = cross(Ez00, Ey00), *const_cast<vec3*>(&Px1) = cross(Ey10, Ez10); *const_cast<vec3*>(&Px0) /= Px0.mod(), *const_cast<vec3*>(&Px1) /= Px1.mod(), *const_cast<vec3*>(&Py0) = cross(Ex00, Ez00), *const_cast<vec3*>(&Py1) = cross(Ez01, Ex10); *const_cast<vec3*>(&Py0) /= Py0.mod(), *const_cast<vec3*>(&Py1) /= Py1.mod(), *const_cast<vec3*>(&Pz0) = cross(Ey00, Ex00), *const_cast<vec3*>(&Pz1) = cross(Ex01, Ey01); *const_cast<vec3*>(&Pz0) /= Pz0.mod(), *const_cast<vec3*>(&Pz1) /= Pz1.mod();
			if (M.det() < 0) *const_cast<vec3*>(&Px0) = -Px0, *const_cast<vec3*>(&Px1) = -Px1, *const_cast<vec3*>(&Py0) = -Py0, *const_cast<vec3*>(&Py1) = -Py1, *const_cast<vec3*>(&Pz0) = -Pz0, *const_cast<vec3*>(&Pz1) = -Pz1;
			*const_cast<double*>(&px0) = -dot(Px0, P000), *const_cast<double*>(&px1) = -dot(Px1, P100), *const_cast<double*>(&py0) = -dot(Py0, P000), *const_cast<double*>(&py1) = -dot(Py1, P111), *const_cast<double*>(&pz0) = -dot(Pz0, P000), *const_cast<double*>(&pz1) = -dot(Pz1, P001);
			*const_cast<vec3*>(&_Ex00) = -Ex00, *const_cast<vec3*>(&_Ex01) = -Ex01, *const_cast<vec3*>(&_Ex10) = -Ex10, *const_cast<vec3*>(&_Ex11) = -Ex11, *const_cast<vec3*>(&_Ey00) = -Ey00, *const_cast<vec3*>(&_Ey01) = -Ey01, *const_cast<vec3*>(&_Ey10) = -Ey10, *const_cast<vec3*>(&_Ey11) = -Ey11, *const_cast<vec3*>(&_Ez00) = -Ez00, *const_cast<vec3*>(&_Ez01) = -Ez01, *const_cast<vec3*>(&_Ez10) = -Ez10, *const_cast<vec3*>(&_Ez11) = -Ez11;
			*const_cast<vec3*>(&PET_Ex00_n1) = _Ey00 - dot(Ex00, _Ey00) * Ex00, *const_cast<vec3*>(&PET_Ex01_n1) = _Ey01 - dot(Ex01, _Ey01) * Ex01, *const_cast<vec3*>(&PET_Ex10_n1) = Ey00 - dot(Ex10, Ey00) * Ex10, *const_cast<vec3*>(&PET_Ex11_n1) = Ey01 - dot(Ex11, Ey01) * Ex11, *const_cast<vec3*>(&PET_Ey00_n1) = _Ex00 - dot(Ey00, _Ex00) * Ey00, *const_cast<vec3*>(&PET_Ey01_n1) = _Ex01 - dot(Ey01, _Ex01) * Ey01, *const_cast<vec3*>(&PET_Ey10_n1) = Ex00 - dot(Ey10, Ex00) * Ey10, *const_cast<vec3*>(&PET_Ey11_n1) = Ex01 - dot(Ey11, Ex01) * Ey11, *const_cast<vec3*>(&PET_Ez00_n1) = _Ex00 - dot(Ez00, _Ex00) * Ez00, *const_cast<vec3*>(&PET_Ez01_n1) = _Ex10 - dot(Ez01, _Ex10) * Ez01, *const_cast<vec3*>(&PET_Ez10_n1) = Ex00 - dot(Ez10, Ex00) * Ez10, *const_cast<vec3*>(&PET_Ez11_n1) = Ex10 - dot(Ez11, Ex10) * Ez11;
			*const_cast<vec3*>(&PET_Ex00_n1) /= PET_Ex00_n1.mod(), *const_cast<vec3*>(&PET_Ex01_n1) /= PET_Ex01_n1.mod(), *const_cast<vec3*>(&PET_Ex10_n1) /= PET_Ex10_n1.mod(), *const_cast<vec3*>(&PET_Ex11_n1) /= PET_Ex11_n1.mod(), *const_cast<vec3*>(&PET_Ey00_n1) /= PET_Ey00_n1.mod(), *const_cast<vec3*>(&PET_Ey01_n1) /= PET_Ey01_n1.mod(), *const_cast<vec3*>(&PET_Ey10_n1) /= PET_Ey10_n1.mod(), *const_cast<vec3*>(&PET_Ey11_n1) /= PET_Ey11_n1.mod(), *const_cast<vec3*>(&PET_Ez00_n1) /= PET_Ez00_n1.mod(), *const_cast<vec3*>(&PET_Ez01_n1) /= PET_Ez01_n1.mod(), *const_cast<vec3*>(&PET_Ez10_n1) /= PET_Ez10_n1.mod(), *const_cast<vec3*>(&PET_Ez11_n1) /= PET_Ez11_n1.mod();
			*const_cast<vec3*>(&PET_Ex00_n2) = _Ez00 - dot(Ex00, _Ez00) * Ex00, *const_cast<vec3*>(&PET_Ex01_n2) = Ez00 - dot(Ex01, Ez00) * Ex01, *const_cast<vec3*>(&PET_Ex10_n2) = _Ez01 - dot(Ex10, _Ez01) * Ex10, *const_cast<vec3*>(&PET_Ex11_n2) = Ez01 - dot(Ex11, Ez01) * Ex11, *const_cast<vec3*>(&PET_Ey00_n2) = _Ez00 - dot(Ey00, _Ez00) * Ey00, *const_cast<vec3*>(&PET_Ey01_n2) = Ez00 - dot(Ey01, Ez00) * Ey01, *const_cast<vec3*>(&PET_Ey10_n2) = _Ez10 - dot(Ey10, _Ez10) * Ey10, *const_cast<vec3*>(&PET_Ey11_n2) = Ez10 - dot(Ey11, Ez10) * Ey11, *const_cast<vec3*>(&PET_Ez00_n2) = _Ey00 - dot(Ez00, _Ey00) * Ez00, *const_cast<vec3*>(&PET_Ez01_n2) = Ey00 - dot(Ez01, Ey00) * Ez01, *const_cast<vec3*>(&PET_Ez10_n2) = _Ey10 - dot(Ez10, _Ey10) * Ez10, *const_cast<vec3*>(&PET_Ez11_n2) = Ey10 - dot(Ez11, Ey10) * Ez11;
			*const_cast<vec3*>(&PET_Ex00_n2) /= PET_Ex00_n2.mod(), *const_cast<vec3*>(&PET_Ex01_n2) /= PET_Ex01_n2.mod(), *const_cast<vec3*>(&PET_Ex10_n2) /= PET_Ex10_n2.mod(), *const_cast<vec3*>(&PET_Ex11_n2) /= PET_Ex11_n2.mod(), *const_cast<vec3*>(&PET_Ey00_n2) /= PET_Ey00_n2.mod(), *const_cast<vec3*>(&PET_Ey01_n2) /= PET_Ey01_n2.mod(), *const_cast<vec3*>(&PET_Ey10_n2) /= PET_Ey10_n2.mod(), *const_cast<vec3*>(&PET_Ey11_n2) /= PET_Ey11_n2.mod(), *const_cast<vec3*>(&PET_Ez00_n2) /= PET_Ez00_n2.mod(), *const_cast<vec3*>(&PET_Ez01_n2) /= PET_Ez01_n2.mod(), *const_cast<vec3*>(&PET_Ez10_n2) /= PET_Ez10_n2.mod(), *const_cast<vec3*>(&PET_Ez11_n2) /= PET_Ez11_n2.mod();
			*const_cast<point*>(&border.Min) = PMin(PMin(PMin(P000, P001), PMin(P010, P011)), PMin(PMin(P100, P101), PMin(P110, P111))), *const_cast<point*>(&border.Max) = PMax(PMax(PMax(P000, P001), PMax(P010, P011)), PMax(PMax(P100, P101), PMax(P110, P111)));
			return border;
		}

		friend ostream& operator << (ostream& os, const Box_affine &B) {
			B.MaxMin();
			os << "Polyline(" << B.P000 << "," << B.P100 << "," << B.P101 << "," << B.P111 << "," << B.P110 << "," << B.P010 << "," << B.P011 << "," << B.P001 << "," <<
				B.P000 << "," << B.P010 << "," << B.P011 << "," << B.P111 << "," << B.P110 << "," << B.P100 << "," << B.P101 << "," << B.P001 << ")" << endl;
			os << "P_{x0}: "; print_quadrilateral(os, B.P000, B.P001, B.P011, B.P010); os << endl;
			os << "P_{x1}: "; print_quadrilateral(os, B.P100, B.P101, B.P111, B.P110); os << endl;
			os << "P_{y0}: "; print_quadrilateral(os, B.P000, B.P001, B.P101, B.P100); os << endl;
			os << "P_{y1}: "; print_quadrilateral(os, B.P010, B.P011, B.P111, B.P110); os << endl;
			os << "P_{z0}: "; print_quadrilateral(os, B.P000, B.P010, B.P110, B.P100); os << endl;
			os << "P_{z1}: "; print_quadrilateral(os, B.P001, B.P011, B.P111, B.P101); os << endl;
			return os;
		}
		virtual unsigned telltype() { return XObjs_BoxAffine_Sign; }

		/* This part doesn't require high efficiency */
		point center() {
			return M * point(0.5, 0.5, 0.5);
		}
		void transform(matrix3D M) {
			this->M = matrix3D_affine(M)*this->M;
		}
		void transform(matrix3D_affine M) {
			this->M = M * this->M;
		}
		void scale(double s) {
			M.scale(s, s, s);
		}
		void scale(double x, double y, double z) {
			M.scale(x, y, z);
		}
		void scale(point P, double s) {
			M.translate(-P.x, -P.y, -P.z);
			M.scale(s, s, s);
			M.translate(P.x, P.y, P.z);
		}
		void scale(point P, double x, double y, double z) {
			M.translate(-P.x, -P.y, -P.z);
			M.scale(x, y, z);
			M.translate(P.x, P.y, P.z);
		}
		void translate(double x, double y, double z) {
			M.translate(x, y, z);
		}
		void operator += (point P) {
			M.translate(P.x, P.y, P.z);
		}
		void operator -= (point P) {
			M.translate(-P.x, -P.y, -P.z);
		}
		void rotate(double rz) {
			M.rotate(0, 0, rz);
		}
		void rotate(point P, double rz) {
			M.translate(-P.x, -P.y, -P.z);
			M.rotate(0, 0, rz);
			M.translate(P.x, P.y, P.z);
		}
		void rotate(double rx, double rz) {
			M.rotate(rx, 0, rz);
		}
		void rotate(point P, double rx, double rz) {
			M.translate(-P.x, -P.y, -P.z);
			M.rotate(rx, 0, rz);
			M.translate(P.x, P.y, P.z);
		}
		void rotate(double rx, double ry, double rz) {
			M.rotate(rx, ry, rz);
		}
		void rotate(point P, double rx, double ry, double rz) {
			M.translate(-P.x, -P.y, -P.z);
			M.rotate(rx, ry, rz);
			M.translate(P.x, P.y, P.z);
		}
		void perspective(double x, double y, double z) {
			M.perspective(x, y, z);
			// 1/x, 1/y, 1/z, 1/(x+1), 1/(y+1), 1/(z+1)
		}
	};

	class BezierQuadratic_xOy : public XObjs_Comp {
		point2D A, B, C; double a, b;
	public:
		BezierQuadratic_xOy() { a = b = 0; }
		BezierQuadratic_xOy(const point2D &A, const point2D &B, const point2D &C) {
			this->A = A, this->B = B, this->C = C;
			a = b = 0;
		}
		BezierQuadratic_xOy(const point2D &A, const point2D &B, const point2D &C, double h) {
			this->A = A, this->B = B, this->C = C;
			a = 0, b = h;
			if (h < 0) swap(a, b);
		}
		BezierQuadratic_xOy(const point2D &A, const point2D &B, const point2D &C, double a, double b) {
			this->A = A, this->B = B, this->C = C;
			this->a = a, this->b = b;
			if (a > b) swap(a, b);
		}
		BezierQuadratic_xOy(const BezierQuadratic_xOy &other) {
			this->A = other.A, this->B = other.B, this->C = other.C;
			this->a = other.a, this->b = other.b;
		}
		BezierQuadratic_xOy* clone() const { return new BezierQuadratic_xOy(*this); }
		~BezierQuadratic_xOy() {}
		virtual unsigned telltype() { return XObjs_BezierQuadraticXOY_Sign; }

		double SDF(const point &P) const {
			double sd = QuadraticBezier_SDF(point2D(P.x, P.y), A, B, C);
			if (P.z > b) return max(P.z - b, sd);
			if (P.z < a) return max(a - P.z, sd);
			return sd;
		}
		borderbox MaxMin() const {
			point2D min = PMin(A, C), max = PMax(A, C);
			point2D nm = A - B, dc = A - 2 * B + C;
			double t = nm.x / dc.x;
			if (t > 0 && t < 1) min = PMin(min, ((1 - t)*(1 - t))*A + (2 * t*(1 - t))*B + (t * t)*C),
				max = PMax(max, ((1 - t)*(1 - t))*A + (2 * t*(1 - t))*B + (t * t)*C);
			t = nm.y / dc.y;
			if (t > 0 && t < 1) min = PMin(min, ((1 - t)*(1 - t))*A + (2 * t*(1 - t))*B + (t * t)*C),
				max = PMax(max, ((1 - t)*(1 - t))*A + (2 * t*(1 - t))*B + (t * t)*C);
			return borderbox(point(min.x, min.y, a), point(max.x, max.y, b));
		}
	};
}

// Windows null pointer: 0x00000000-0x0000FFFF
// 0xABCC    A: number of operated objects, 0 means apply to point;  B: number of extra parameters;  CC: operation #
#define CSG_UnionOp_Sign 0x2001
#define CSG_IntersectionOp_Sign 0x2002
#define CSG_SubtractionOp_Sign 0x2003
#define CSG_ComplementOp_Sign 0x1004
#define CSG_RoundingOp_Sign 0x1105
#define CSG_OnionOp_Sign 0x1106
#define CSG_Translation_Sign 0x0107
#define CSG_Rotation_Sign 0x0108
#define CSG_SmoothUnionOp_Sign 0x2109

#define CSG_UnionOperator ((const XObjs::XObjs_Comp*)CSG_UnionOp_Sign)
#define CSG_IntersectionOperator ((const XObjs::XObjs_Comp*)CSG_IntersectionOp_Sign)
#define CSG_SubtractionOperator ((const XObjs::XObjs_Comp*)CSG_SubtractionOp_Sign)
#define CSG_ComplementOperator ((const XObjs::XObjs_Comp*)CSG_ComplementOp_Sign)
#define CSG_RoundingOperator ((const XObjs::XObjs_Comp*)CSG_RoundingOp_Sign)
#define CSG_OnionOperator ((const XObjs::XObjs_Comp*)CSG_OnionOp_Sign)
#define CSG_TranslationOperator ((const XObjs::XObjs_Comp*)CSG_Translation_Sign)
#define CSG_RotationOperator ((const XObjs::XObjs_Comp*)CSG_Rotation_Sign)
#define CSG_SmoothUnionOperator ((const XObjs::XObjs_Comp*)CSG_SmoothUnionOp_Sign)


#include <thread>
#include <mutex>

#define XSolid_Sign 0x10000001
#define XSolid_Smooth 0x10000010
#define XSolid_Diffuse 0x10000100
#define XSolid_Crystal 0x10001000
class XSolid : public object {
	vector<const XObjs::XObjs_Comp*> objs;
	vector<void**> objs_tp;
#define XSolid_Max_Parameter 4
	inline void objs_tp_push() {
		objs_tp.push_back(new void*[XSolid_Max_Parameter]);
		for (unsigned m = 0; m < XSolid_Max_Parameter; m++) this->objs_tp.back()[m] = 0;
	}
	void deep_copy(const XSolid &other) {
		// assign objs_dir, don't copy object
		for (unsigned i = 0; i < other.objs_tp.size(); i++) {
			if (other.objs_tp.at(i) == 0) this->objs_tp.push_back(0);
			else {
				this->objs_tp.push_back(new void*[XSolid_Max_Parameter]);
				for (unsigned m = 0; m < XSolid_Max_Parameter; m++) this->objs_tp.back()[m] = 0;
				switch (unsigned(other.objs.at(i))) {
				case CSG_RoundingOp_Sign:;
				case CSG_OnionOp_Sign:;
				case CSG_SmoothUnionOp_Sign: {
					this->objs_tp.back()[0] = new double(*((double*)(other.objs_tp.at(i)[0])));
					break;
				}
				case CSG_Translation_Sign: {
					this->objs_tp.back()[0] = new point(*((point*)(other.objs_tp.at(i)[0])));
					break;
				}
				case CSG_Rotation_Sign: {
					this->objs_tp.back()[0] = new matrix3D(*((matrix3D*)(other.objs_tp.at(i)[0])));
					break;
				}
				}
			}
		}
	}
	unsigned Max_StackLength; unsigned Stack_Usage[4]; double* *Stack; bool TIed; 	// 4 rendering threads + 1 attempting thread
public:
	borderbox border;

	unsigned type;	// smooth / diffuse / crystal
	rgblight col;	// or damping ratio
	double refractive_index;

	XSolid() {
		Stack = 0;
	}
	XSolid(const XObjs::XObjs_Comp &a) {
		objs.push_back(a.clone());
		objs_tp.push_back(0);
		border = a.MaxMin();
		Stack = 0;
	}
	XSolid(const XSolid &other) {
		for (int i = 0; i < other.objs.size(); i++) {
			if (unsigned(other.objs.at(i)) >> 16) objs.push_back(other.objs.at(i)->clone());
			else objs.push_back(other.objs.at(i));
		}
		this->deep_copy(other);
		this->border = other.border;
		this->type = other.type, this->col = other.col;
		this->TIed = false, Stack = 0;
	}
	XSolid& operator = (const XSolid &other) {
		this->destruct();
		for (int i = 0; i < other.objs.size(); i++) {
			if (unsigned(other.objs.at(i)) >> 16) objs.push_back(other.objs.at(i)->clone());
			else objs.push_back(other.objs.at(i));
		}
		this->deep_copy(other);
		this->border = other.border;
		Stack = 0;
		return *this;
	}
	XSolid& operator = (const XObjs::XObjs_Comp &other) {
		this->destruct();
		objs.push_back(other.clone()), objs_tp.push_back(0);
		border = other.MaxMin();
		Stack = 0;
		return *this;
	}
	object* clone() const {
		return new XSolid(*this);
	}
	void destruct() {
		for (unsigned i = 0, l = objs.size(); i < l; i++) if (unsigned(objs.at(i)) >> 16) delete objs.at(i);
		for (unsigned i = 0, l = objs_tp.size(); i < l; i++) {
			if (objs_tp.at(i) != 0) {
				for (unsigned j = 0; j < XSolid_Max_Parameter; j++) {
					if (objs_tp.at(i)[j] == 0) break;
					else delete objs_tp.at(i)[j], objs_tp.at(i)[j] = 0;
				}
				delete objs_tp.at(i);
			}
		}
		objs.clear(), objs_tp.clear();
		if (Stack != 0) {
			delete Stack[0]; delete Stack[1]; delete Stack[2]; delete Stack[3];
			delete Stack; Stack = 0;
		}
	}
	~XSolid() {
		this->destruct();
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

		if (Stack != 0) {
			delete Stack[0]; delete Stack[1]; delete Stack[2]; delete Stack[3]; delete Stack;
		}
		Stack = new double*[4]; TIed = false;
		for (int i = 0; i < 4; i++) Stack[i] = new double[Max_StackLength], Stack_Usage[i] = false;

	}

	double SDF(point P, const unsigned &_thread_No) const {
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
				S[dir - 1] -= *((double*)objs_tp.at(i)[0]);
				break;
			}
			case CSG_OnionOp_Sign: {
				S[dir - 1] = abs(S[dir - 1]) - *((double*)objs_tp.at(i)[0]);
				break;
			}
			case CSG_Translation_Sign: {
				P += *((point*)objs_tp.at(i)[0]);
				break;
			}
			case CSG_Rotation_Sign: {
				P *= *((matrix3D*)objs_tp.at(i)[0]);
				break;
			}
			case CSG_SmoothUnionOp_Sign: {
				dir--;
				//S[dir - 1] = -*((double*)objs_tp.at(i)[0]) * log(exp(-S[dir] / *((double*)objs_tp.at(i)[0])) + exp(-S[dir - 1] / *((double*)objs_tp.at(i)[0])));
				double h = clamp(0.5 + 0.5*(S[dir] - S[dir - 1]) / *((double*)objs_tp.at(i)[0]), 0.0, 1.0);
				S[dir - 1] = mix(S[dir], S[dir - 1], h) - *((double*)objs_tp.at(i)[0]) * h*(1.0 - h);
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
	void meet(intersect &R, const ray &a, double Max_Dist) const {
#ifdef _4_Threads_Rendering
		unsigned thread_no = 0, k = *((unsigned*)(&this_thread::get_id()));
#define STACKUSAGE_XSOLID const_cast<unsigned*>(Stack_Usage)
		TIed ? (k == STACKUSAGE_XSOLID[0] ? thread_no = 0 : (k == STACKUSAGE_XSOLID[1] ? thread_no = 1 :
			(k == STACKUSAGE_XSOLID[2] ? thread_no = 2 : (k == STACKUSAGE_XSOLID[3] ? thread_no = 3 : (STACKUSAGE_XSOLID[0] = k) && (thread_no = 0)))))
			: (STACKUSAGE_XSOLID[0] ? (k == STACKUSAGE_XSOLID[0] ? thread_no = 0 : STACKUSAGE_XSOLID[1] ? (k == STACKUSAGE_XSOLID[1] ? thread_no = 0 :
				STACKUSAGE_XSOLID[2] ? (k == STACKUSAGE_XSOLID[2] ? thread_no = 0 : (STACKUSAGE_XSOLID[3] = k) && (thread_no = 3) && (*const_cast<bool*>(&TIed) = true))
				: (STACKUSAGE_XSOLID[2] = k) && (thread_no = 2)) : (STACKUSAGE_XSOLID[1] = k) && (thread_no = 1)) : (STACKUSAGE_XSOLID[0] = k) && (thread_no = 0));
		R.wt = thread_no;
#else
		const unsigned thread_no = 0;	// optimized by compiler
		R.wt = 0;
#endif
		R.meet = 0;
		if (!border.meet(a)) return;
		double t = ERR_ZETA, sdf; point P;
		P = a.orig + t * a.dir, R.vt = sdf = this->SDF(P, thread_no);
		if (sdf > 0) {
			do {
				t += sdf;
				if (t > Max_Dist) return;
				P = a.orig + t * a.dir, sdf = this->SDF(P, thread_no);
			} while (sdf > ERR_EPSILON);
			R.dist = t, R.intrs = P;
			R.meet = true;
			return;
		}
		else {
			do {
				t -= sdf;
				if (t > Max_Dist) return;
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

	friend XSolid CSG_UnionOp(const XSolid &A, const XSolid &B) {	// not exact
		XSolid X;
		for (int i = 0; i < A.objs.size(); i++) X.objs.push_back(unsigned(A.objs.at(i)) >> 16 ? A.objs.at(i)->clone() : A.objs.at(i));
		for (int i = 0; i < B.objs.size(); i++) X.objs.push_back(unsigned(B.objs.at(i)) >> 16 ? B.objs.at(i)->clone() : B.objs.at(i));
		X.deep_copy(A), X.deep_copy(B);
		X.objs.push_back(CSG_UnionOperator), X.objs_tp.push_back(0);
		X.border.Max = PMax(A.border.Max, B.border.Max), X.border.Min = PMin(A.border.Min, B.border.Min);
		return X;
	}
	friend XSolid CSG_IntersectionOp(const XSolid &A, const XSolid &B) {	// not exact
		XSolid X;
		for (int i = 0; i < A.objs.size(); i++) X.objs.push_back(unsigned(A.objs.at(i)) >> 16 ? A.objs.at(i)->clone() : A.objs.at(i));
		for (int i = 0; i < B.objs.size(); i++) X.objs.push_back(unsigned(B.objs.at(i)) >> 16 ? B.objs.at(i)->clone() : B.objs.at(i));
		X.deep_copy(A), X.deep_copy(B);
		X.objs.push_back(CSG_IntersectionOperator), X.objs_tp.push_back(0);
		X.border.Max = PMin(A.border.Max, B.border.Max), X.border.Min = PMax(A.border.Min, B.border.Min);
		return X;
	}
	friend XSolid CSG_SubtractionOp(const XSolid &A, const XSolid &B) {		// not exact
		XSolid X;
		for (int i = 0; i < A.objs.size(); i++) X.objs.push_back(unsigned(A.objs.at(i)) >> 16 ? A.objs.at(i)->clone() : A.objs.at(i));
		for (int i = 0; i < B.objs.size(); i++) X.objs.push_back(unsigned(B.objs.at(i)) >> 16 ? B.objs.at(i)->clone() : B.objs.at(i));
		X.deep_copy(A), X.deep_copy(B);
		X.objs.push_back(CSG_SubtractionOperator), X.objs_tp.push_back(0);
		X.border = A.border;	// need to optimize
		return X;
	}
	friend XSolid CSG_ComplementOp(const XSolid &A) {	// exact
		XSolid X;
		for (int i = 0; i < A.objs.size(); i++) X.objs.push_back(unsigned(A.objs.at(i)) >> 16 ? A.objs.at(i)->clone() : A.objs.at(i));
		X.deep_copy(A);
		X.objs.push_back(CSG_ComplementOperator), X.objs_tp.push_back(0);
		X.border.Max = point(INFINITY, INFINITY, INFINITY), X.border.Min = -X.border.Max;
		return X;
	}
	friend XSolid CSG_RoundingOp(const XSolid &A, double r) {	// makes the object larger, doesn't apply for Union/Intersection/Subtraction objects
		XSolid X;
		for (int i = 0; i < A.objs.size(); i++) X.objs.push_back(unsigned(A.objs.at(i)) >> 16 ? A.objs.at(i)->clone() : A.objs.at(i));
		X.deep_copy(A);
		X.objs.push_back(CSG_RoundingOperator), X.objs_tp_push();
		X.objs_tp.back()[0] = new double(r);
		X.border.Max = A.border.Max + point(r, r, r), X.border.Min = A.border.Min - point(r, r, r);
		//cout << A.border << endl << X.border << endl;
		return X;
	}
	friend XSolid CSG_OnionOp(const XSolid &A, double r) {	// make larger, problems may occur to concave side
		XSolid X;
		for (int i = 0; i < A.objs.size(); i++) X.objs.push_back(unsigned(A.objs.at(i)) >> 16 ? A.objs.at(i)->clone() : A.objs.at(i));
		X.deep_copy(A);
		X.objs.push_back(CSG_OnionOperator), X.objs_tp_push();
		X.objs_tp.back()[0] = new double(r);
		X.border.Max = A.border.Max + point(r, r, r), X.border.Min = A.border.Min - point(r, r, r);
		return X;
	}
	friend XSolid CSG_Translation(const XSolid &A, const vec3 &P) {		// exact
		XSolid X;
		X.objs.push_back(CSG_TranslationOperator), X.objs_tp_push();
		X.objs_tp.back()[0] = new point(-P);
		for (int i = 0; i < A.objs.size(); i++) X.objs.push_back(unsigned(A.objs.at(i)) >> 16 ? A.objs.at(i)->clone() : A.objs.at(i));
		X.deep_copy(A);
		X.objs.push_back(CSG_TranslationOperator), X.objs_tp_push();
		X.objs_tp.back()[0] = new point(P);
		X.border.Max = A.border.Max + P, X.border.Min = A.border.Min + P;
		//cout << A.border << endl << X.border << endl;
		return X;
	}
	friend XSolid CSG_Rotation(const XSolid &A, double rx, double ry, double rz) {	// exact, rotate about origin
		matrix3D M(Rotation, rx, ry, rz);
		XSolid X;
		X.objs.push_back(CSG_RotationOperator), X.objs_tp_push();
		X.objs_tp.back()[0] = new matrix3D(M.invert());
		for (int i = 0; i < A.objs.size(); i++) X.objs.push_back(unsigned(A.objs.at(i)) >> 16 ? A.objs.at(i)->clone() : A.objs.at(i));
		X.deep_copy(A);
		X.objs.push_back(CSG_RotationOperator), X.objs_tp_push();
		X.objs_tp.back()[0] = new matrix3D(M);

		X.border.Max = PMax(PMax(PMax(M*point(A.border.Min.x, A.border.Min.y, A.border.Min.z), M*point(A.border.Min.x, A.border.Min.y, A.border.Max.z)),
			PMax(M*point(A.border.Min.x, A.border.Max.y, A.border.Min.z), M*point(A.border.Max.x, A.border.Min.y, A.border.Min.z))),
			PMax(PMax(M*point(A.border.Max.x, A.border.Max.y, A.border.Max.z), M*point(A.border.Max.x, A.border.Max.y, A.border.Min.z)),
				PMax(M*point(A.border.Max.x, A.border.Min.y, A.border.Max.z), M*point(A.border.Min.x, A.border.Max.y, A.border.Max.z))));
		X.border.Min = PMin(PMin(PMin(M*point(A.border.Min.x, A.border.Min.y, A.border.Min.z), M*point(A.border.Min.x, A.border.Min.y, A.border.Max.z)),
			PMin(M*point(A.border.Min.x, A.border.Max.y, A.border.Min.z), M*point(A.border.Max.x, A.border.Min.y, A.border.Min.z))),
			PMin(PMin(M*point(A.border.Max.x, A.border.Max.y, A.border.Max.z), M*point(A.border.Max.x, A.border.Max.y, A.border.Min.z)),
				PMin(M*point(A.border.Max.x, A.border.Min.y, A.border.Max.z), M*point(A.border.Min.x, A.border.Max.y, A.border.Max.z))));
		//cout << A.border << endl << X.border << endl;
		return X;
	}
	friend XSolid CSG_UnionOp_Smooth(const XSolid &A, const XSolid &B, double r) {
		XSolid X;
		for (int i = 0; i < A.objs.size(); i++) X.objs.push_back(unsigned(A.objs.at(i)) >> 16 ? A.objs.at(i)->clone() : A.objs.at(i));
		for (int i = 0; i < B.objs.size(); i++) X.objs.push_back(unsigned(B.objs.at(i)) >> 16 ? B.objs.at(i)->clone() : B.objs.at(i));
		X.deep_copy(A), X.deep_copy(B);
		X.objs.push_back(CSG_SmoothUnionOperator), X.objs_tp_push();
		X.objs_tp.back()[0] = new double(r);
		X.border.Max = PMax(A.border.Max, B.border.Max), X.border.Min = PMin(A.border.Min, B.border.Min);
		X.border.Max += point(r, r, r), X.border.Min -= point(r, r, r);
		return X;
	}
};



#define GXSolid_Sign 0x10000002
class GXSolid : public object {
	vector<XSolid*> V;
	borderbox border;
public:
	GXSolid() {}
	GXSolid(initializer_list<XSolid> A) {
		for (int i = 0; i < A.size(); i++) V.push_back(new XSolid(*(A.begin() + i)));
	}
	GXSolid(const GXSolid &other) {
		for (int i = 0; i < other.V.size(); i++) V.push_back(new XSolid(*other.V.at(i)));
	}
	object* clone() const {
		return new GXSolid(*this);
	}
	~GXSolid() {
		for (int i = 0; i < V.size(); i++) delete V.at(i);
		V.clear();
	}

	void clear() {
		for (int i = 0; i < V.size(); i++) delete V.at(i);
		V.clear();
	}
	void push(XSolid A) {
		V.push_back(new XSolid(A));
	}
	void init() {
		if (V.empty()) {
			border.Min = border.Max = point(0, 0, 0); return;
		}
		V.at(0)->init(); border = V.at(0)->border;
		for (int i = 0; i < V.size(); i++) {
			V.at(i)->init();
			border.Min = PMin(border.Min, V.at(i)->border.Min);
			border.Max = PMin(border.Max, V.at(i)->border.Max);
		}
	}

	void meet(intersect &R, const ray &a, double Max_Dist) const {
		for (int i = 0; i < V.size(); i++) {

		}
	}

};


rgblight SolarDeepsea(double t) {	// color function, -1 < t < 1
	double r, g, b;
	if (t >= 0) r = ((0.880183 * t - 2.27602) * t + 1.94992) * t + 0.442456, g = ((-1.15153 * t + 1.85887) * t + 0.0899384) * t + 0.0120354, b = ((-0.204874 * t + 0.464309) * t - 0.151561) * t + 0.0178089;
	else t = -t,
		r = (((((-7.31545 * t + 13.31) * t - 3.75123) * t - 1.83143) * t - 0.560149) * t + 0.779769) * t + 0.156216,
		g = (((((-3.53306 * t + 14.265) * t - 21.3179) * t + 13.3875) * t - 2.16183) * t + 0.286874) * t - 9.48835e-05,
		b = (((((5.62199 * t - 16.1938) * t + 17.7046) * t - 9.38477) * t + 1.85801) * t + 1.08822) * t + 0.29789;
	if (r > 1) r = 1; if (g > 1) g = 1; if (b > 1) b = 1;
	if (r < 0) r = 0; if (g < 0) g = 0; if (b < 0) b = 0;
	return rgblight(r, g, b);
}
void VisualizeSDF_core(XSolid &X, parallelogram &p, bitmap &canvas) {
	double u, v, sdf, x, y, mag, arg; point P;
	unsigned w = canvas.width() / 2, h = canvas.height() / 2;

	for (unsigned i = 0; i < h; i++) {
		for (unsigned j = 0; j < w; j++) {
			u = double(j) / w, v = double(i) / h;
			P = p.O + u * p.A + v * p.B;
			sdf = X.SDF(P, 0);
			x = X.SDF(p.O + (u + ERR_ZETA) * p.A + v * p.B, 0) - sdf;
			y = X.SDF(p.O + u * p.A + (v + ERR_ZETA) * p.B, 0) - sdf;
			x /= ERR_ZETA, y /= ERR_ZETA;
			sdf *= 10, x *= 10, y *= 10;
			mag = abs(x) + abs(y);
			arg = atan2(y, x - y);
			if (isnan(arg)) mag = arg = 0;

			// General Visualization of Magnitude
			double t = tanh(0.5 * sdf);
			double r, g, b;
			canvas[i + h][j] = rgb(SolarDeepsea(t));

			// Magnitude with Contour, difference = 0.1
			t = pow(cos(10 * PI * sdf), 12); r = g = b = t;
			if (abs(sdf) < 0.05) r = b = 0;
			canvas[i + h][j + w] = drgb(r, g, b);

			// Derivative of SDF
			canvas[i][j] = fromHSL(arg / (2 * PI), 0.8, 1 - pow(0.7, log(log(mag + 1) + 1.05)));

			// Derivative with Contour
			t = cos(10 * PI * arg)/* + cos((2 * PI * mag))*/;
			t = pow(abs(t) > 1 ? 1 : t, 20);
			canvas[i][j + w] = drgb(t, t, t);

			canvas[i][j + w] = rgb(0.5*(rgblight(canvas[i][j]) + rgblight(canvas[i + h][j])));
		}
	}

}
void VisualizeSDF(XSolid X, parallelogram p, double S) {
	X.init();
	double w = p.A.mod(), h = p.B.mod(), t = h / w;
	w = sqrt(S / t), h = w * t;
	bitmap canvas(2 * ceil(w), 2 * ceil(h));
	VisualizeSDF_core(X, p, canvas);
	canvas.out("IMAGE\\SDF.bmp");
}
void VisualizeSDF(XSolid X, parallelogram p, double S, double u, double v) {
	X.init();
	double w = p.A.mod(), h = p.B.mod(), t = h / w;
	w = sqrt(S / t), h = w * t;

	X.SDF(p.O + u * p.A + v * p.B, 0);

	bitmap canvas(2 * ceil(w), 2 * ceil(h));
	VisualizeSDF_core(X, p, canvas);

	canvas.dot(u*w - 1, v*h + h, Green); canvas.dot(u*w, v*h + h - 1, Green); canvas.dot(u*w + 1, v*h + h, Green); canvas.dot(u*w, v*h + h + 1, Green);


	canvas.out("IMAGE\\SDF.bmp");
}

void ScanXSolid(XSolid X, unsigned Nx, unsigned Ny, unsigned Nz, double S) {
	X.init();
	point C = 0.5*(X.border.Max + X.border.Min);
	double Dx = X.border.Max.x - X.border.Min.x, Dy = X.border.Max.y - X.border.Min.y, Dz = X.border.Max.z - X.border.Min.z;
	system("md IMAGE\\SDF");
	system("del /f /s /q IMAGE\\SDF\\*.bmp /s");

	string filename;

	double w = Dx, h = Dy, t = h / w; w = sqrt(S / t), h = w * t;
	parallelogram P(X.border.Min, point(Dx, 0, 0), point(0, Dy, 0));
	bitmap canvas(2 * ceil(w), 2 * ceil(h));
	unsigned digits = log10(Nz) + 1;
	for (unsigned i = 1; i <= Nz; i++) {
		VisualizeSDF_core(X, P, canvas);
		canvas.out(&(("IMAGE\\SDF\\Z" + uint2str(i, digits) + ".bmp")[0]));
		P.O.z += Dz / (Nz - 1);
	}

	w = Dx, h = Dz, t = h / w; w = sqrt(S / t), h = w * t;
	P = parallelogram(X.border.Min, point(Dx, 0, 0), point(0, 0, Dz));
	canvas = bitmap(2 * ceil(w), 2 * ceil(h));
	digits = log10(Ny) + 1;
	for (unsigned i = 1; i <= Ny; i++) {
		VisualizeSDF_core(X, P, canvas);
		canvas.out(&(("IMAGE\\SDF\\Y" + uint2str(i, digits) + ".bmp")[0]));
		P.O.y += Dy / (Ny - 1);
	}

	w = Dy, h = Dz, t = h / w; w = sqrt(S / t), h = w * t;
	P = parallelogram(X.border.Min, point(0, Dy, 0), point(0, 0, Dz));
	canvas = bitmap(2 * ceil(w), 2 * ceil(h));
	digits = log10(Nx) + 1;
	for (unsigned i = 1; i <= Nx; i++) {
		VisualizeSDF_core(X, P, canvas);
		canvas.out(&(("IMAGE\\SDF\\X" + uint2str(i, digits) + ".bmp")[0]));
		P.O.x += Dx / (Nx - 1);
	}
}
