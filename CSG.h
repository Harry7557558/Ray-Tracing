#pragma once

#ifndef Object_Sign
#include "Object.h"
#endif

// Reference: https://www.iquilezles.org/www/articles/distfunctions/distfunctions.htm

// All object/SDF/operation/transformations labled "exact" or "not exact" are strictly prooved or disproved
// object exact: legal shape (not an approximation);  SDF exact: legal distance field (issues won't occur when applying operations); 
// "not exact" operations/transformations may affect farther operations/transformations


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
		~Sphere() {}

		double SDF(const point &P) const {	// exact
			return (P - C).mod() - r;
		}
		borderbox MaxMin() const {
			return borderbox(C - point(r, r, r), C + point(r, r, r));
		}
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
		~Plane() {}

		double SDF(const point &P) const {	// exact
			return dot(N, P) - D;
		}
		borderbox MaxMin() const {
			return borderbox(point(-INFINITY, -INFINITY, -INFINITY), point(INFINITY, INFINITY, INFINITY));
		}
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
			point U(abs(P1.x - P2.x), abs(P1.y - P2.y), abs(P1.z - P2.z));
			U *= r / h;
			U.x = sqrt(r*r - U.x*U.x), U.y = sqrt(r*r - U.y*U.y), U.z = sqrt(r*r - U.z*U.z);
			borderbox B(P1, P2); B.fix();
			B.Max += U, B.Min -= U;
			return B;
		}
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
			this->r = r, this->h = h, this->l = l; alpha = atan(r / h);
		}
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
	};

	class Cone_capped : public XObjs_Comp {		// exact
		point C, P1, P2; double r1, r2, h1, h2, l1, l2, alpha; vec3 dir;
	public:
		Cone_capped() { r1 = r2 = h1 = h2 = l1 = l2 = alpha = 0; }
		Cone_capped(const point &P_1, const point &P_2, const double &r_1, const double &r_2) {
			P1 = P_1, P2 = P_2; r1 = r_1, r2 = r_2; if (r1 > r2) swap(r1, r2), swap(P1, P2);
			dir = P2 - P1; double h = dir.mod(); dir /= h; alpha = atan((r2 - r1) / h);
			l1 = r1 / sin(alpha), l2 = r2 / sin(alpha); h1 = r1 / tan(alpha), h2 = r2 / tan(alpha);
			C = P1 - h1 * dir;
		}
		Cone_capped(const point &C, const vec3 &dir, const double &h_1, const double &h_2, const double &alpha) {
			this->C = C, this->dir = dir / dir.mod(); h1 = h_1, h2 = h_2, this->alpha = alpha;
			P1 = C + h1 * dir, P2 = C + h2 * dir; r1 = h1 * tan(alpha), r2 = h2 * tan(alpha), l1 = h1 / cos(alpha), l2 = h2 / cos(alpha);
		}
		Cone_capped(const Cone_capped &other) {
			C = other.C, P1 = other.P1, P2 = other.P2; dir = other.dir;
			r1 = other.r1, r2 = other.r2, h1 = other.h1, h2 = other.h2, l1 = other.l1, l2 = other.l2, alpha = other.alpha;
		}
		~Cone_capped() {}

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
	};

	class Box_affine : public XObjs_Comp {
		matrix3D_affine M;	// applying to an unit cube, the only variable
		matrix3D_affine M_invert;
		point P000, P100, P010, P001, P110, P101, P011, P111;
		vec3 Ex00, Ex01, Ex10, Ex11, Ey00, Ey01, Ey10, Ey11, Ez00, Ez01, Ez10, Ez11;
		vec3 Px0, Px1, Py0, Py1, Pz0, Pz1; double px0, px1, py0, py1, pz0, pz1;
	public:
		Box_affine() { M = matrix3D_affine(1, 0, 0, 0, 1, 0, 0, 0, 1); }
		Box_affine(const matrix3D_affine &M) {
			this->M = M;
		}
		~Box_affine() {}

		double SDF(const point &P) const {

			point S = M_invert * P;
			bool xi = S.x > 0 && S.x < 1, yi = S.y > 0 && S.y < 1, zi = S.z > 0 && S.z < 1;

			if (xi && yi && zi) return -min(min(min(abs(dot(P, Px0) + px0), abs(dot(P, Px1) + px1)),
				min(abs(dot(P, Py0) + py0), abs(dot(P, Py1) + py1))), min(abs(dot(P, Pz0) + pz0), abs(dot(P, Pz1) + pz1)));

			if (!xi && !yi && !zi) return min(min(min((P - P000).mod(), (P - P100).mod()), min((P - P010).mod(), (P - P001).mod())),
				min(min((P - P011).mod(), (P - P101).mod()), min((P - P110).mod(), (P - P111).mod())));

			if (xi + yi + zi == 2) {
				if (!xi) return min(abs(dot(P, Px0) + px0), abs(dot(P, Px1) + px1));
				if (!yi) return min(abs(dot(P, Py0) + py0), abs(dot(P, Py1) + py1));
				if (!zi) return min(abs(dot(P, Pz0) + pz0), abs(dot(P, Pz1) + pz1));
			}

			if (xi + yi + zi == 1) {
				if (xi) return min(min(abs(cross(P - P000, Ex00).mod()), abs(cross(P - P001, Ex01).mod())),
					min(abs(cross(P - P010, Ex10).mod()), abs(cross(P - P011, Ex11).mod())));
				if (yi) return min(min(abs(cross(P - P000, Ey00).mod()), abs(cross(P - P001, Ey01).mod())),
					min(abs(cross(P - P100, Ey10).mod()), abs(cross(P - P101, Ey11).mod())));
				if (zi) return min(min(abs(cross(P - P000, Ez00).mod()), abs(cross(P - P010, Ez01).mod())),
					min(abs(cross(P - P100, Ez10).mod()), abs(cross(P - P110, Ez11).mod())));
			}

			return 0;
		}

		borderbox MaxMin() const {
			point P1 = M * point(0, 0, 0), P2 = M * point(1, 0, 0), P3 = M * point(0, 1, 0), P4 = M * point(0, 0, 1),
				P5 = M * point(0, 1, 1), P6 = M * point(1, 0, 1), P7 = M * point(1, 1, 0), P8 = M * point(1, 1, 1);

			// Since this function will be called before rendering, it would be okay to invert the matrix there.
			*const_cast<matrix3D_affine*>(&M_invert) = M.invert();
			*const_cast<point*>(&P000) = M * point(0, 0, 0), *const_cast<point*>(&P001) = M * point(0, 0, 1);
			*const_cast<point*>(&P100) = M * point(1, 0, 0), *const_cast<point*>(&P010) = M * point(0, 1, 0);
			*const_cast<point*>(&P110) = M * point(1, 1, 0), *const_cast<point*>(&P101) = M * point(1, 0, 1);
			*const_cast<point*>(&P011) = M * point(0, 1, 1), *const_cast<point*>(&P111) = M * point(1, 1, 1);
			*const_cast<vec3*>(&Ex00) = P100 - P000, *const_cast<vec3*>(&Ex01) = P101 - P001, *const_cast<vec3*>(&Ex10) = P110 - P010, *const_cast<vec3*>(&Ex11) = P111 - P011;
			*const_cast<vec3*>(&Ey00) = P010 - P000, *const_cast<vec3*>(&Ey01) = P011 - P001, *const_cast<vec3*>(&Ey10) = P110 - P100, *const_cast<vec3*>(&Ey11) = P111 - P101;
			*const_cast<vec3*>(&Ez00) = P001 - P000, *const_cast<vec3*>(&Ez01) = P011 - P010, *const_cast<vec3*>(&Ez10) = P101 - P100, *const_cast<vec3*>(&Ez11) = P111 - P110;
			*const_cast<vec3*>(&Px0) = cross(Ez00, Ey00), *const_cast<vec3*>(&Px1) = cross(Ey10, Ez10); *const_cast<vec3*>(&Px0) /= Px0.mod(), *const_cast<vec3*>(&Px1) /= Px1.mod();
			*const_cast<vec3*>(&Py0) = cross(Ex00, Ez00), *const_cast<vec3*>(&Py1) = cross(Ez11, Ex11); *const_cast<vec3*>(&Py0) /= Py0.mod(), *const_cast<vec3*>(&Py1) /= Py1.mod();
			*const_cast<vec3*>(&Pz0) = cross(Ey00, Ex00), *const_cast<vec3*>(&Pz1) = cross(Ex01, Ey01); *const_cast<vec3*>(&Pz0) /= Pz0.mod(), *const_cast<vec3*>(&Pz1) /= Pz1.mod();
			*const_cast<vec3*>(&Ex00) /= Ex00.mod(), *const_cast<vec3*>(&Ex01) /= Ex01.mod(), *const_cast<vec3*>(&Ex10) /= Ex10.mod(), *const_cast<vec3*>(&Ex11) /= Ex11.mod();
			*const_cast<vec3*>(&Ey00) /= Ey00.mod(), *const_cast<vec3*>(&Ey01) /= Ey01.mod(), *const_cast<vec3*>(&Ey10) /= Ey10.mod(), *const_cast<vec3*>(&Ey11) /= Ey11.mod();
			*const_cast<vec3*>(&Ez00) /= Ez00.mod(), *const_cast<vec3*>(&Ez01) /= Ez01.mod(), *const_cast<vec3*>(&Ez10) /= Ez10.mod(), *const_cast<vec3*>(&Ez11) /= Ez11.mod();
			*const_cast<double*>(&px0) = -dot(Px0, P000), *const_cast<double*>(&px1) = -dot(Px1, P100);
			*const_cast<double*>(&py0) = -dot(Py0, P000), *const_cast<double*>(&py1) = -dot(Py1, P111);
			*const_cast<double*>(&pz0) = -dot(Pz0, P000), *const_cast<double*>(&pz1) = -dot(Pz1, P001);

			return borderbox(PMin(PMin(PMin(P1, P2), PMin(P3, P4)), PMin(PMin(P5, P6), PMin(P7, P8))),
				PMax(PMax(PMax(P1, P2), PMax(P3, P4)), PMax(PMax(P5, P6), PMax(P7, P8))));
		}
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

	friend XSolid CSG_UnionOp(const XSolid &A, const XSolid &B) {	// not exact
		XSolid X;
		for (int i = 0; i < A.objs.size(); i++) X.objs.push_back(A.objs.at(i));
		for (int i = 0; i < B.objs.size(); i++) X.objs.push_back(B.objs.at(i));
		X.deep_copy(A), X.deep_copy(B);
		X.objs.push_back(CSG_UnionOperator), X.objs_tp.push_back(0);
		X.border.Max = PMax(A.border.Max, B.border.Max), X.border.Min = PMin(A.border.Min, B.border.Min);
		return X;
	}
	friend XSolid CSG_IntersectionOp(const XSolid &A, const XSolid &B) {	// not exact
		XSolid X;
		for (int i = 0; i < A.objs.size(); i++) X.objs.push_back(A.objs.at(i));
		for (int i = 0; i < B.objs.size(); i++) X.objs.push_back(B.objs.at(i));
		X.deep_copy(A), X.deep_copy(B);
		X.objs.push_back(CSG_IntersectionOperator), X.objs_tp.push_back(0);
		X.border.Max = PMin(A.border.Max, B.border.Max), X.border.Min = PMax(A.border.Min, B.border.Min);
		return X;
	}
	friend XSolid CSG_SubtractionOp(const XSolid &A, const XSolid &B) {		// not exact
		XSolid X;
		for (int i = 0; i < A.objs.size(); i++) X.objs.push_back(A.objs.at(i));
		for (int i = 0; i < B.objs.size(); i++) X.objs.push_back(B.objs.at(i));
		X.deep_copy(A), X.deep_copy(B);
		X.objs.push_back(CSG_SubtractionOperator), X.objs_tp.push_back(0);
		X.border = A.border;	// need to optimize
		return X;
	}
	friend XSolid CSG_ComplementOp(const XSolid &A) {	// exact
		XSolid X;
		for (int i = 0; i < A.objs.size(); i++) X.objs.push_back(A.objs.at(i));
		X.deep_copy(A);
		X.objs.push_back(CSG_ComplementOperator), X.objs_tp.push_back(0);
		X.border.Max = point(INFINITY, INFINITY, INFINITY), X.border.Min = -X.border.Max;
		return X;
	}
	friend XSolid CSG_RoundingOp(const XSolid &A, double r) {	// makes the object larger, doesn't apply for Union/Intersection/Subtraction objects
		XSolid X;
		for (int i = 0; i < A.objs.size(); i++) X.objs.push_back(A.objs.at(i));
		X.deep_copy(A);
		X.objs.push_back(CSG_RoundingOperator), X.objs_tp_push();
		X.objs_tp.back()[0] = new double(r);
		X.border.Max = A.border.Max + point(r, r, r), X.border.Min = A.border.Min - point(r, r, r);
		return X;
	}
	friend XSolid CSG_OnionOp(const XSolid &A, double r) {	// make larger, problems may occur to concave side
		XSolid X;
		for (int i = 0; i < A.objs.size(); i++) X.objs.push_back(A.objs.at(i));
		X.deep_copy(A);
		X.objs.push_back(CSG_OnionOperator), X.objs_tp_push();
		X.objs_tp.back()[0] = new double(r);
		X.border.Max = A.border.Max + point(r, r, r), X.border.Min = A.border.Min - point(r, r, r);
		return X;
	}
};



void VisualizeSDF(XSolid &X, parallelogram p, double S) {
	X.init();
	double w = p.A.mod(), h = p.B.mod(), t = h / w;
	w = sqrt(S / t), h = w * t;

	bitmap canvas(2 * ceil(w), 2 * ceil(h));
	double u, v, sdf, x, y, mag, arg; point P;

	for (unsigned i = 0; i < h; i++) {
		for (unsigned j = 0; j < w; j++) {
			u = j / w, v = i / h;
			P = p.O + u * p.A + v * p.B;
			sdf = X.SDF(P, 0);
			x = X.SDF(p.O + (u + ERR_ZETA) * p.A + v * p.B, 0) - sdf;
			y = X.SDF(p.O + u * p.A + (v + ERR_ZETA) * p.B, 0) - sdf;
			x /= ERR_ZETA, y /= ERR_ZETA;
			mag = 1 / sqrt(x * x + y * y);
			arg = atan2(y, x - y);

			// General Visualization of Magnitude
			double t = tanh(0.5 * sdf);
			double r, g, b;
			if (t >= 0) r = ((0.880183 * t - 2.27602) * t + 1.94992) * t + 0.442456, g = ((-1.15153 * t + 1.85887) * t + 0.0899384) * t + 0.0120354, b = ((-0.204874 * t + 0.464309) * t - 0.151561) * t + 0.0178089;
			else t = -t,
				r = (((((-7.31545 * t + 13.31) * t - 3.75123) * t - 1.83143) * t - 0.560149) * t + 0.779769) * t + 0.156216,
				g = (((((-3.53306 * t + 14.265) * t - 21.3179) * t + 13.3875) * t - 2.16183) * t + 0.286874) * t - 9.48835e-05,
				b = (((((5.62199 * t - 16.1938) * t + 17.7046) * t - 9.38477) * t + 1.85801) * t + 1.08822) * t + 0.29789;
			if (r > 1) r = 1; if (g > 1) g = 1; if (b > 1) b = 1;
			if (r < 0) r = 0; if (g < 0) g = 0; if (b < 0) b = 0;
			canvas[i + h][j] = drgb(r, g, b);

			// Magnitude with Contour, difference = 0.2
			t = pow(cos(10 * PI * sdf), 12); r = g = b = t;
			canvas[i + h][j + int(w)] = drgb(r, g, b);

			// Derivative of SDF
			canvas[i][j] = fromHSL(arg / (2 * PI), 1, 1 - pow(0.4, log(log(mag + 1) + 1.05)));

			// Magnitude + Change
			canvas[i][j + int(w)] = rgb(0.5*(rgblight(canvas[i][j]) + rgblight(canvas[i + h][j])));
		}
	}

	canvas.out("IMAGE\\SDF.bmp");
}

