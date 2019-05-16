#pragma once

#include "Object.h"

#include <thread>
#include <mutex>
#include <Windows.h>

extern ofstream fout("D:\\Coding\\AboutPhysics\\RayTracing\\RayTracing\\IMAGE\\Log.txt");

class World {
public:
	vector<object*> Objs;	// series of objects
	vector<World*> GObjs;	// series of sub-worlds
	point N;	// direction of global light source (ideal sky, suppose it's fine, white and cloudness)
	rgblight background;
	class borderbox {
		// Use for border of subgroups of objects
	public:
		point B1, B2;
		borderbox() {}
		inline point Max() {
			return B2;
		}
		inline point Min() {
			return B1;
		}
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
	borderbox border;
	void resize() {
		point maxc, minc;
		if (!Objs.empty()) maxc = Objs.front()->Max(), minc = Objs.front()->Min();
		for (int i = 1; i < Objs.size(); i++) {
			maxc = Max(maxc, Objs.at(i)->Max());
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
		border.B1 = minc, border.B2 = maxc;
	}
	friend int main();
	friend void Render_GTest01();
	friend void Render_GTest02();
	friend void Render_GTest03();
public:
	World() {}
	World(const World &W) {
		Objs.resize(W.Objs.size());
		for (int i = 0; i < Objs.size(); i++) {
			switch (W.Objs.at(i)->telltype()) {
			case Plane_Sign: {
				Objs.at(i) = new plane((plane*)(W.Objs.at(i)));
				break;
			}
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
			case Cylinder_Sign: {
				Objs.at(i) = new cylinder((cylinder*)(W.Objs.at(i)));
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
		background = W.background;
		N = W.N;
		return;
	}
	~World() {
		Objs.clear(); 	// DO NOT delete ANY elements there !!!
		GObjs.clear();
	}
	inline void insert(World *a) {
		GObjs.push_back(a);
		this->resize();
	}
	inline void insert(World *a, double x, double y, double z) {
		*a += point(x, y, z);
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
	void setGlobalLightSource(double lrx, double lrz) {
		N = point(cos(lrx), sqrt(1 - cos(lrx)*cos(lrx) - cos(lrz)*cos(lrz)), cos(lrz));
	}
	void setGlobalLightSource(double x, double y, double z) {
		N = point(x, y, z);
	}
	void operator += (point V) {
		for (int i = 0; i < Objs.size(); i++) {
			*Objs.at(i) += V;
		}
		for (int i = 0; i < GObjs.size(); i++) {
			*GObjs.at(i) += V;
		}
		border.B1 += V, border.B2 += V;
	}
	void rotate_mult(const double &rx, const double &ry, const double &rz, const double &tms) {
		for (int i = 0; i < Objs.size(); i++) {
			Objs.at(i)->rotate(rx, ry, rz); *Objs.at(i) *= tms;
		}
		for (int i = 0; i < GObjs.size(); i++) {
			GObjs.at(i)->rotate_mult(rx, ry, rz, tms);
		}
		this->resize();
	}

	// Find the nearest object
	void RayTracing_EnumObjs(const ray &v, intersect &ni, object* &no, intersect& nt) const {
		int NB = Objs.size();
		for (int k = 0; k < NB; k++) {
			Objs.at(k)->meet(nt, v);
			if (nt.meet) {
				if (no == 0 || (ni.dist > nt.dist)) {
					ni = nt, no = Objs.at(k);
				}
			}
		}	// Find the nearest object among single objects

		NB = GObjs.size();
		for (int i = 0; i < NB; i++) {
			if (GObjs.at(i)->border.meet(v)) {
				GObjs.at(i)->RayTracing_EnumObjs(v, ni, no, nt);
			}
		}
	}

	// Check if the viewing point is inside a 3D object (eg.under the water). If yes, return its adress; or return 0;
	object3D* Render_CheckContaining(const point &W) const {
		int NB = Objs.size();
		for (int k = 0; k < NB; k++) {
			if (Objs.at(k)->telltype() >> 16 == 1) {
				if (((object3D*)(Objs.at(k)))->contain(W)) {
					return (object3D*)(Objs.at(k));
				}
			}
		}
		object3D* p;
		NB = GObjs.size();
		for (int i = 0; i < NB; i++) {
			if (GObjs.at(i)->border.contain(W)) {
				p = GObjs.at(i)->Render_CheckContaining(W);
				if (p != 0) return p;
			}
		}
		return 0;
	}

	// Calculate the result rgb color with a giving ray
	double RayTracing(const ray &v, rgblight &c, const double count, const int n) const {
		if (count < 0.001 || n > 30) { c = rgblight(1, 1, 1); return INFINITY; }
		if (isnan(count)) { c = rgblight(NAN, NAN, NAN); return NAN; }
		//for (int i = 0; i < n; i++) fout << " ";

		//for (int i = 0; i < n; i++) fout << "    "; fout << fixed << setprecision(3) << count << "\t" << defaultfloat << v << endl;

		int NB = Objs.size();
		intersect ni; object* no = 0;	// keeps the nearest object and its intersect data
		const World* Wp = this;	// World of the nearest object
		intersect nt;

		RayTracing_EnumObjs(v, ni, no, nt);

		if (no == 0) {
			//fout << v << endl;
			double ang = dot(N, v.dir);
			if (ang <= 0) {
				c = background; return INFINITY;
			}
			ang /= N.mod()*v.dir.mod();
			c = rgblight(ang); return INFINITY;
		}
		else if (isnan(ni.reflect.x)) {
			// must exit function when nan occur; or may occur endless recursive
			c = rgblight(NAN, NAN, NAN); return NAN;
		}
		else {
			//fout << ray(v.orig, ni.intrs - v.orig) << endl;
			// Meet
			switch (no->telltype() >> 16) {
			case 0: {	// surface
				switch ((no->telltype() >> 8) & 0b1111) {
				case 0: {	// smooth opacity surface
					RayTracing(ray(ni.intrs, ni.reflect), c, count * ((objectSF*)no)->reflect.vsl(), n + 1);
					c.r *= ((objectSF*)no)->reflect.r, c.g *= ((objectSF*)no)->reflect.g, c.b *= ((objectSF*)no)->reflect.b;	// should use Fresnel's formula
					break;

					/*// brightness of refraction: proper to exp(1-sec(θ))
					// (1-r)*L*exp(1-1/cos(t))+r*L, which L is the brightness of the light source, t is the angle of incidence, and r is the ratio of reflection calculated with Fresnel equations
					// when 0<2x<π, cos(x)=sqrt((cos(2x)+1)/2), and cos(2x) can be calculated with the angle between two vectors
					const double n1 = 1.0, n2 = 1.5;
					double ct = dot(v.dir, ni.reflect) / (v.dir.mod()*ni.reflect.mod()); ct = sqrt(0.5*(ct + 1));	// cosine of angle of incidence
					double ct2 = sqrt(1 - (n1*n1) / (n2*n2)*(1 - ct * ct));		// cosine of angle of refraction
					double refl = 0.5 * (pow((n1*ct - n2 * ct2) / (n1*ct + n2 * ct2), 2) + pow((n1*ct2 - n2 * ct) / (n1*ct2 + n2 * ct), 2));	// ratio of reflection
					c = ((1 - refl)*exp(1 - 1 / ct)*((objectSF*)no)->reflect) * c + refl * c;	// probably bug

					break;*/
				}
				case 1: {	// diffuse reflection, opacity
					point Ni = ni.reflect;
					//((objectSF_dif*)no)->rotate_vec(Ni);
					rotate_normal(Ni);
					RayTracing(ray(ni.intrs, Ni), c, count * ((objectSF_dif*)no)->reflect.vsl(), n + 1);
					c.r *= ((objectSF_dif*)no)->reflect.r, c.g *= ((objectSF_dif*)no)->reflect.g, c.b *= ((objectSF_dif*)no)->reflect.b;
					break;
				}
				case 2: {	// surface with color
					rgblight d; ((objectSF_col*)no)->getcol(ni, d);
					RayTracing(ray(ni.intrs, ni.reflect), c, count * d.vsl(), n + 1);
					c.r *= d.r, c.g *= d.g, c.b *= d.b;
					break;
				}
				}
				break;
			}
			case 1: {	// 3D objects with volumn, debugging
				point refract; double rlr;
				((object3D*)no)->refractData(ni, v, 1, refract, rlr);
				//fout << count << " " << rlr << " " << (refract.z > 0 ? "U" : "D") << endl;

				double d = RayTracing(ray(ni.intrs, ni.reflect), c, count * rlr, n + 1);
				if (ni.ut == 0) { // obj -> air
					c.r *= exp(-((object3D*)no)->attcoe.r * d), c.g *= exp(-((object3D*)no)->attcoe.g * d), c.b *= exp(-((object3D*)no)->attcoe.b * d);
				}
				if (!isnan(refract.z)) {
					rgblight c1; double d1 = RayTracing(ray(ni.intrs, refract), c1, count * (1 - rlr), n + 1);
					if (ni.ut == 1) {	// air -> obj
						c1.r *= exp(-((object3D*)no)->attcoe.r * d1), c1.g *= exp(-((object3D*)no)->attcoe.g * d1), c1.b *= exp(-((object3D*)no)->attcoe.b * d1);
					}
					if (!isnan(c1.b)) c = c * rlr + c1 * (1 - rlr);
				}
				//else fout << "nan" << endl;
				break;
			}
			case 0x100: {	// light source
				// Light sources don't look well when global lightsource defined

				RayTracing(ray(ni.intrs, ni.reflect), c, count * (1 - ni.ut), n + 1);
				c *= ni.ut; c += ((lightsource*)no)->col * ni.ut;
				//c = ((lightsource*)no)->col * ni.ut;
				break;
			}
			}
			return ni.dist;
		}
	}

	// oi: object that contains ray end, must match
	rgblight CalcRGB(const ray &v, const object3D* oi) {
		rgblight c;
		double d = -RayTracing(v, c, 1.0, 0);
		if (oi != 0) {
			c.r *= exp(oi->attcoe.r * d), c.g *= exp(oi->attcoe.g * d), c.b *= exp(oi->attcoe.b * d);
		}
		return c;
	}

	// Time Recorder
	void RenderingProcessCounter(int ps, int &T1, int &T2, int &T3, int &T4) {
		int m, n = 0;
		int sum;
		while (1) {
			sum = T1 + T2 + T3 + T4;
			m = (sum * 1000.0) / ps;	// using float operation, sometimes overflow occurs when sum is large enough
			if (m != n) cout << "\r" << (m / 10) << "." << (m % 10) << "%";
			n = m;
			if (sum >= ps) break;
			Sleep(50);
		}
	}

	// Multithread entrance, note that "old" and "bad" functions are discarded
	void MultiThread_CC(bitmap &canvas, int begw, int endw, int begh, int endh, const point W, const parallelogram &sc, const object3D* oi, int &RenderingProcess) {
		ray beg;
		beg.orig = W;
		double u, v;
		rgblight c, s;
		//*➤*/ auto t0 = NTime::now(); auto t1 = NTime::now(); fsec fs = t1 - t0; pixel* p; double a;

		unsigned NP;
		for (unsigned i = begw; i < endw; i++) {
			for (unsigned j = begh; j < endh; j++) {
				NP = 0;
				s.r = s.g = s.b = 0;
				for (unsigned m = 0; m < Render_Sampling; m++) {
					for (unsigned n = 0; n < Render_Sampling; n++) {
						u = (i + double(m) / Render_Sampling) / canvas.width(), v = (j + double(n) / Render_Sampling) / canvas.height();
						beg.orig = W;	// necessary
						beg.dir = (1 - u - v)*sc.O + u * sc.A + v * sc.B - beg.orig;
						//*➤*/ t0 = NTime::now();
						c = CalcRGB(beg, oi);
						//*➤*/ t1 = NTime::now(); fs = t1 - t0;
						if (!(isnan(c.r) || isnan(c.g) || isnan(c.b))) {
							s += c;
							NP++;
						}
					}
				}
				if (NP != 0) c = s / NP;
				else {
					if (j != begh && i != begw && j + 1 != endh) canvas[j][i] = rgb((rgblight(canvas[j - 1][i]) + rgblight(canvas[j][i - 1])) / 3
						+ (rgblight(canvas[j - 1][i - 1]) + rgblight(canvas[j + 1][i - 1])) / 6);
					else if (j == begh && j + 1 != endh) canvas[j][i] = rgb(rgblight(canvas[j][i - 1])*0.666667 + rgblight(canvas[j + 1][i - 1]) / 3);
					else if (i == begw && j != begh) canvas[j][i] = canvas[j - 1][i];
					else;
				}

				if (c.r >= 1) c.r = 0.9999; if (c.g >= 1) c.g = 0.9999; if (c.b >= 1) c.b = 0.9999;
				//c.r = sqrt(1 - (c.r - 1)*(c.r - 1)), c.g = sqrt(1 - (c.g - 1)*(c.g - 1)), c.b = sqrt(1 - (c.b - 1)*(c.b - 1));	// Make it lighter
				canvas.dot(i, j, drgb(c.r, c.g, c.b));
				//c *= 2; canvas.dot(i, j, drgb(1 - exp(-c.r), 1 - exp(-c.g), 1 - exp(-c.b)));
				//canvas.dot(i, j, drgb(tanh(c.r), tanh(c.g), tanh(c.b)));

				//*➤*/ p = canvas[j] + i; p->r = p->b = 0, p->g /= 2; a = (tanh(log2(fs.count() * 10000)) + 1) / 2; canvas.dot(i, j, drgb(a, 0, 0), 1 - a);
				//*➤*/ if (i == begw || i == endw - 1 || j == begh || j == begh - 1) canvas.dot(i, j, color(Yellow));
				RenderingProcess++;
			}
		}
	}



	// class won't be dagamed after rendering
	/* Parameters (all (solid) angles are in radians):
		canvas      canvas for rendering
		C           center of camera
		O           point appears in the center of the canvas
		rt          rotation of screen of camera
		sr          solid angle which the camera watches the canvas
	*/
	static unsigned Render_Sampling;
	inline double rf_SR(double w, double r, double t) const {
		// calculate solid angle with given parameters, for somehow use
		double h = t * w;
		double aw = acos((2 * r*r - w * w) / (2 * r*r));
		double ah = acos((2 * r*r - h * h) / (2 * r*r));
		double ac = acos((2 * r*r - w * w - h * h) / (2 * r*r));
		double am = (aw + ah + ac) / 2;
		double SR = 8 * atan(sqrt(tan(am / 2) * tan((am - aw) / 2) * tan((am - ah) / 2) * tan((am - ac) / 2)));
		return SR;
	}
	void render(bitmap &canvas, point C, point O, double rt, double sr) {

		cout << "Initializing...";
		auto Time_Beg = NTime::now();
		canvas.clear();

		const point OC = O - C;
		const double r = OC.mod();

		// numerical solve the width and height of canvas in world coordinate
		double t = double(canvas.height()) / double(canvas.width());
		double x = 1, _x, xc, y;
		if (rf_SR(x, r, t) > sr) {
			while (rf_SR(x, r, t) > sr) _x = x, x /= 2;
			swap(_x, x);
		}
		else if (rf_SR(x, r, t) < sr) {
			while (rf_SR(x, r, t) < sr) _x = x, x *= 2;
		}
		if (rf_SR(x, r, t) != 0) {
			for (int i = 0; i < 60; i++) {
				xc = (x + _x) / 2;
				y = rf_SR(xc, r, t);
				if (y == sr) break;
				else if (y > sr) x = xc;
				else if (y < sr) _x = xc;
			}
		}
		const double w = x, h = t * x;
		const double rs = 0.5 * sqrt(4 * r*r - w * w - h * h);

		// calculate rotation
		double orx = acos(OC.z / r), orz = atan2(OC.x, -OC.y);
		const double rx = atan2(sin(orx)*cos(rt), cos(orx)), ry = atan2(-sin(orx)*sin(rt), hypot(sin(orx)*cos(rt), cos(orx))),
			rz = atan2(sin(orz)*cos(rt) + sin(rt)*cos(orx)*cos(orz), cos(rt)*cos(orz) - sin(rt)*cos(orx)*sin(orz));		// first x, then y, finally z
		Matrix R = Matrix({ {cos(rz),-sin(rz),0}, {sin(rz),cos(rz),0}, {0,0,1} }) *
			Matrix({ { cos(ry),0,sin(ry) }, {0,1,0}, {-sin(ry),0,cos(ry)} }) *
			Matrix({ {1,0,0}, {0,cos(rx),-sin(rx)}, {0,sin(rx),cos(rx)} });
		parallelogram CV(point(w / 2, -h / 2, rs), point(-w / 2, -h / 2, rs), point(w / 2, h / 2, rs), true);
		CV *= R;
		CV += C;


		this->resize();
		object3D* oi = Render_CheckContaining(C);


		cout << "\nAttempting...";

		// calculate pixels allocates to each thread
		ray beg; beg.orig = C;
		rgblight c, s;
		fsec fs;
#ifndef DEBUG
		const int STEP = 8;
		vector<double> attempt; attempt.resize(canvas.width() / STEP);
		double u, v;
		auto t0 = NTime::now();
		auto t1 = NTime::now();
		fs = t1 - t0;
		for (int i = 0; i < canvas.width(); i += STEP) {
			t0 = NTime::now();
			for (int j = 0; j < canvas.height(); j += STEP) {
				unsigned NP = 0; s.r = s.g = s.b = 0;
				for (unsigned m = 0; m < Render_Sampling; m++) {
					for (unsigned n = 0; n < Render_Sampling; n++) {
						u = (i + double(m) / Render_Sampling) / canvas.width(), v = (j + double(n) / Render_Sampling) / canvas.height();
						beg.orig = C;
						beg.dir = (1 - u - v) * CV.O + u * CV.A + v * CV.B - beg.orig;
						c = CalcRGB(beg, oi);
						if (!(isnan(c.r) || isnan(c.g) || isnan(c.b))) {
							s += c;
							NP++;
						}
					}
				}
				canvas.dot(i, j, rgb(c));
				//cout << "* ";
			}
			//cout << endl;
			t1 = NTime::now();
			fs = t1 - t0;
			attempt.at(i / STEP) = fs.count();
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
		B1 *= STEP, B2 *= STEP, B3 *= STEP;
		if (B3 == -STEP) B3 = canvas.width();
		if (B2 == -STEP) B2 = canvas.width();
		if (B1 == -STEP) B1 = canvas.width();
		sum = floor(sum * STEP * STEP * 1.2);

		cout << "\r             \r"; if (sum >= 1) cout << "Estimated Time Required: " << int(sum) << "secs. \n\n";
		cout << "Rendering...\n";

		// multithread rendering
		int RenderingProcess0 = 0, RenderingProcess1 = 0, RenderingProcess2 = 0, RenderingProcess3 = 0;
		thread T0([&](World* WC) { WC->MultiThread_CC(canvas, 0, B1, 0, canvas.height(), C, CV, oi, RenderingProcess0); }, this);
		thread T1([&](World* WC) { WC->MultiThread_CC(canvas, B1, B2, 0, canvas.height(), C, CV, oi, RenderingProcess1); }, this);
		thread T2([&](World* WC) { WC->MultiThread_CC(canvas, B2, B3, 0, canvas.height(), C, CV, oi, RenderingProcess2); }, this);
		thread T3([&](World* WC) { WC->MultiThread_CC(canvas, B3, canvas.width(), 0, canvas.height(), C, CV, oi, RenderingProcess3); }, this);
		thread Proc([&](World* WC) { WC->RenderingProcessCounter(canvas.height()*canvas.width(),
			RenderingProcess0, RenderingProcess1, RenderingProcess2, RenderingProcess3); }, this);
		T0.join(); T1.join(); T2.join(); T3.join(); Proc.join();

#endif

		cout << " Completed. \n\n";
		auto Time_End = NTime::now();
		fs = Time_End - Time_Beg;
		cout << "Elapsed Time: " << setprecision(3) << fs.count() << "s. \n\n\n";



		// Debug single pixel
		int debugx = 245, debugy = 222, RP = 0;
		this->MultiThread_CC(canvas, debugx, debugx + 1, debugy, debugy + 1, C, CV, oi, RP);
		//canvas.dot(debugx + 1, debugy, Red); canvas.dot(debugx - 1, debugy, Red); canvas.dot(debugx, debugy + 1, Red); canvas.dot(debugx, debugy - 1, Red);

	}

	friend ostream& operator << (ostream& os, const World &W) {
		for (unsigned i = 0; i < W.Objs.size(); i++) os << *(W.Objs.at(i)) << endl;
		for (unsigned i = 0; i < W.GObjs.size(); i++) os << *(W.GObjs.at(i)) << endl;
		return os;
	}
};

unsigned World::Render_Sampling = 1;

// https://zhuanlan.zhihu.com/p/41269520
