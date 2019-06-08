#pragma once

#include "Object.h"
#include "CSG.h"

#include <thread>
#include <mutex>
#include <Windows.h>

using namespace std;

#ifndef _INC_WORLD_H
#define _INC_WORLD_H

class World {
public:
	vector<const object*> Objs;	// series of objects
	vector<World*> GObjs;	// series of sub-worlds
	point N;	// direction of global light source (ideal sky, suppose it's fine, white and cloudness)
	rgblight background;
	borderbox border;
	void resize() {
		point maxc, minc;
		if (!Objs.empty()) maxc = Objs.front()->Max(), minc = Objs.front()->Min();
		const_cast<object*>(Objs.at(0))->init();
		for (int i = 1; i < Objs.size(); i++) {
			maxc = Max(maxc, Objs.at(i)->Max());
			minc = Min(minc, Objs.at(i)->Min());
			const_cast<object*>(Objs.at(i))->init();
		}
		for (int i = 0; i < Objs.size(); i++) {
			/*if (Objs.at(i)->telltype() == XSolid_Sign) {
				Objs.push_back(Objs.at(i));
				Objs.erase(Objs.begin() + i);
			}*/
		}
		for (int i = 0; i < GObjs.size(); i++) GObjs.at(i)->resize();
		if (Objs.size() == 0 && !GObjs.empty()) {
			maxc = GObjs.front()->border.Max, minc = GObjs.front()->border.Min;
		}
		for (int i = ((Objs.size() == 0 && !GObjs.empty()) ? 1 : 0); i < GObjs.size(); i++) {
			maxc = Max(maxc, GObjs.at(i)->border.Max);
			minc = Min(minc, GObjs.at(i)->border.Min);
		}
		border.Min = minc, border.Max = maxc;
		border.fix();
	}
	bool dynamic_memory;	// indicates whether the destructor should delete elements or not
#ifndef _4_Threads_Rendering
	mutex ML;
#endif
	friend int main();
public:
	World() :dynamic_memory(false) {
		background = rgblight(1, 1, 1);
	}
	World(const World &W) : dynamic_memory(true) {
		for (int i = 0; i < W.Objs.size(); i++) {
			Objs.push_back(W.Objs[i]->clone());
		}
		for (int i = 0; i < W.GObjs.size(); i++) {
			GObjs.push_back(new World(*(W.GObjs.at(i))));
		}
		dynamic_memory = true;
		N = W.N, background = W.background, border = W.border;
	}
	~World() {
		if (dynamic_memory) {
			for (int i = 0; i < Objs.size(); i++) delete Objs.at(i);
			for (int i = 0; i < GObjs.size(); i++) delete GObjs.at(i);
		}
		Objs.clear(); GObjs.clear();
	}
	inline void clear() {
		if (dynamic_memory) {
			for (int i = 0; i < Objs.size(); i++) delete Objs.at(i);
			for (int i = 0; i < GObjs.size(); i++) delete GObjs.at(i);
		}
		Objs.clear(); GObjs.clear();
	}
	inline void insert(World *a) {
		GObjs.push_back(a);
		this->resize();
	}
	inline void add(object* a) {
		Objs.push_back(a);
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

	// Find the closest object
	void RayTracing_EnumObjs(const ray &v, intersect &ni, const object* &no, intersect& nt) const {
		int NB = Objs.size();
		for (int k = 0; k < NB; k++) {
			if (Objs.at(k)->telltype() == XSolid_Sign) {
				((XSolid*)(Objs.at(k)))->meet(nt, v, ni.meet ? ni.dist : ERR_UPSILON);
			}
			else Objs.at(k)->meet(nt, v);
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
		// Note: "object under water" "glass in water" "camera under water" situations are still debugging

		if (count < 0.004 || n > 20) {
			double ang = dot(N, v.dir) / N.mod()*v.dir.mod();
			c = ang > 0 ? ang * background : rgblight(0, 0, 0);
			return INFINITY;
		}
		if (isnan(count)) { c = rgblight(NAN, NAN, NAN); return NAN; }
		//for (int i = 0; i < n; i++) fout << " ";

		//for (int i = 0; i < n; i++) fout << "    "; fout << fixed << setprecision(3) << count << "\t" << defaultfloat << v << endl;

		int NB = Objs.size();
		intersect ni; const object* no = 0;	// keeps the nearest object and its intersect data
		const World* Wp = this;	// World of the nearest object
		intersect nt;

		RayTracing_EnumObjs(v, ni, no, nt);

		if (no == 0) {
			//fout << v << endl;
			double ang = dot(N, v.dir);
			if (ang <= 0) {
				c = rgblight(0, 0, 0); return INFINITY;
			}
			ang /= N.mod()*v.dir.mod();
			c = background * ang; return INFINITY;
		}
		else if (isnan(ni.reflect.x)) {
			// must exit function when nan occur; or may occur endless recursive
			c = rgblight(NAN, NAN, NAN); return NAN;
		}
		else {
			//fout << ray(v.orig, ni.intrs - v.orig) << endl; fout.flush();
			// Meet
			ni.reflect /= ni.reflect.mod();
			switch (no->telltype() >> 16) {
			case 0: {	// surface
				switch ((no->telltype() >> 8) & 0b1111) {
				case 0: {	// smooth opacity surface
					RayTracing(ray(ni.intrs, ni.reflect), c, count * ((objectSF*)no)->reflect.vsl(), n + 1);
					c.r *= ((objectSF*)no)->reflect.r, c.g *= ((objectSF*)no)->reflect.g, c.b *= ((objectSF*)no)->reflect.b;	// should use Fresnel's formula
					if (isnan(c.r) || isnan(c.b) || isnan(c.g)) {
						return ni.dist;
					}
					break;
				}
				case 1: {	// diffuse reflection, opacity
					point Ni = ni.reflect;
					double fr = rotate_normal(Ni);
					RayTracing(ray(ni.intrs, Ni), c, count * ((objectSF_dif*)no)->reflect.vsl() * fr, n + 1);
					c.r *= ((objectSF_dif*)no)->reflect.r * fr, c.g *= ((objectSF_dif*)no)->reflect.g * fr, c.b *= ((objectSF_dif*)no)->reflect.b * fr;
					if (isnan(c.r) || isnan(c.b) || isnan(c.g)) {
						return ni.dist;
					}
					break;
					// I probably found a hidden bug: as the modulus of reflect ray getting larger, the OA effect gets weaker. 
				}
				case 2: {	// surface with color
					rgblight d; ((objectSF_col*)no)->getcol(ni, d);
					RayTracing(ray(ni.intrs, ni.reflect), c, count * d.vsl(), n + 1);
					c.r *= d.r, c.g *= d.g, c.b *= d.b;
					if (isnan(c.r) || isnan(c.b) || isnan(c.g)) {
						return ni.dist;
					}
					break;
				}
				}
				break;
			}
			case 1: {	// 3D objects with volumn, debugging
				point refract; double rlr;
				((object3D*)no)->refractData(ni, v, 1, refract, rlr);
				refract /= refract.mod();
				//fout << count << " " << rlr << " " << (refract.z > 0 ? "U" : "D") << endl;

				double d = RayTracing(ray(ni.intrs, ni.reflect), c, count * rlr, n + 1);
				if (ni.ut == 0) { // obj -> air
					if (((object3D*)no)->attcoe.r != 0) c.r *= exp(-((object3D*)no)->attcoe.r * d);
					if (((object3D*)no)->attcoe.g != 0) c.g *= exp(-((object3D*)no)->attcoe.g * d);
					if (((object3D*)no)->attcoe.b != 0) c.b *= exp(-((object3D*)no)->attcoe.b * d);
				}
				if (!isnan(refract.z)) {
					rgblight c1; double d1 = RayTracing(ray(ni.intrs, refract), c1, count * (1 - rlr), n + 1);
					if (ni.ut == 1) {	// air -> obj
						c1.r *= exp(-((object3D*)no)->attcoe.r * d1), c1.g *= exp(-((object3D*)no)->attcoe.g * d1), c1.b *= exp(-((object3D*)no)->attcoe.b * d1);
					}
					if (!isnan(c1.b)) c = c * rlr + c1 * (1 - rlr);
				}
				if (isnan(c.r) || isnan(c.b) || isnan(c.g)) {
					return ni.dist;
				}
				break;
			}
			case 0x100: {	// light source
				// Light sources don't look well when global lightsource defined

				RayTracing(ray(ni.intrs, ni.reflect), c, count * (1 - ni.ut), n + 1);
				c *= ni.ut; c += ((lightsource*)no)->col * ni.ut;
				//c = ((lightsource*)no)->col * ni.ut;
				if (isnan(c.r) || isnan(c.b) || isnan(c.g)) {
					return ni.dist;
				}
				break;
			}
			case 0x1000: {	// XSolid
				switch (((XSolid*)no)->type) {
				case XSolid_Crystal: {
					// send R.ut as the refractive index of the other media
					// meet = air->obj ? 1 : 0; intrs = refract; reflect = reflect; vt = rate-of-reflection; 
					ni.ut = 1;
					point P = ni.intrs;
					((XSolid*)no)->reflectData(ni, v);
					double d = RayTracing(ray(P, ni.reflect), c, count * ni.vt, n + 1);
					//fout << ni.vt << endl;
					if (!ni.meet) {
						// 0*INF => NAN
						if (((XSolid*)no)->col.r != 0) c.r *= exp(-((XSolid*)no)->col.r * d);
						if (((XSolid*)no)->col.g != 0) c.g *= exp(-((XSolid*)no)->col.g * d);
						if (((XSolid*)no)->col.b != 0) c.b *= exp(-((XSolid*)no)->col.b * d);
					}
					if (!isnan(ni.intrs.z)) {
						rgblight c1; double d1 = RayTracing(ray(P, ni.intrs), c1, count * (1 - ni.vt), n + 1);
						if (ni.meet) {	// air -> obj
							if (((XSolid*)no)->col.r != 0) c1.r *= exp(-((XSolid*)no)->col.r * d1);
							if (((XSolid*)no)->col.g != 0) c1.g *= exp(-((XSolid*)no)->col.g * d1);
							if (((XSolid*)no)->col.b != 0) c1.b *= exp(-((XSolid*)no)->col.b * d1);
						}
						if (!isnan(c1.b)) c = c * ni.vt + c1 * (1 - ni.vt);
					}
					break;
				}
				case XSolid_Diffuse: {

					break;
				}
				default: {	// smooth surface
					((XSolid*)no)->reflectData(ni, v);
					RayTracing(ray(ni.intrs, ni.reflect), c, count * ((XSolid*)no)->col.vsl(), n + 1);
					c.r *= ((XSolid*)no)->col.r, c.g *= ((XSolid*)no)->col.g, c.b *= ((XSolid*)no)->col.b;
					break;
				}
				}

			}
			}
			return ni.dist;
		}
	}

	// oi: object that contains ray end, must match
	rgblight CalcRGB(ray &v, const object3D* oi) {
		v.dir /= v.dir.mod();
		rgblight c;
		double d = -RayTracing(v, c, 1.0, 0);
		if (oi != 0) {
			c.r *= exp(oi->attcoe.r * d), c.g *= exp(oi->attcoe.g * d), c.b *= exp(oi->attcoe.b * d);
		}
		return c;
	}

	// Time Recorder
#ifdef _4_Threads_Rendering
	void RenderingProcessCounter(int ps, int &T1, int &T2, int &T3, int &T4) {
		int m, n = 0;
		int sum;
		while (1) {
			sum = T1 + T2 + T3 + T4;
			m = (sum * 1000.0) / ps;	// using float operation, sometimes overflow occurs when sum is large enough
			if (m != n) cout << "\r" << (m / 10) << "." << (m % 10) << "%";
			n = m;
			if (sum >= ps) break;
			this_thread::sleep_for(chrono::milliseconds(50));
		}
	}
#endif

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
				//*➤*/ t0 = NTime::now();
				for (unsigned m = 0; m < Render_Sampling; m++) {
					for (unsigned n = 0; n < Render_Sampling; n++) {
						u = (i + double(m) / Render_Sampling) / canvas.width(), v = (j + double(n) / Render_Sampling) / canvas.height();
						beg.orig = W;	// necessary
						beg.dir = sc.O + u * sc.A + v * sc.B - beg.orig;
						c = CalcRGB(beg, oi);
						if (!(isnan(c.r) || isnan(c.g) || isnan(c.b))) {
							s += c;
							NP++;
						}
					}
				}
				//*➤*/ t1 = NTime::now(); fs = t1 - t0;
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



	/* Parameters (all (solid) angles are in radians):
		canvas      canvas for rendering
		C           center of camera
		O           point appears in the center of the canvas
		rt          rotation of screen of camera
		sr          solid angle which the camera watches the canvas
	*/
	static unsigned Render_Sampling;
	void render(bitmap &canvas, point C, point O, double rt, double sr) {

#ifndef FoldUp


#ifdef _4_Threads_Rendering
		cout << "Initializing...";
#else
		//ML.lock();
#endif
		auto Time_Beg = NTime::now();
		canvas.clear();

		const point OC = O - C;
		const double r = OC.mod();

		// numerical solve the width and height of canvas in world coordinate
		auto rf_SR = [](double w, double r, double t) -> double {
			// calculate solid angle with given parameters
			double h = t * w;
			double aw = acos((2 * r*r - w * w) / (2 * r*r));
			double ah = acos((2 * r*r - h * h) / (2 * r*r));
			double ac = acos((2 * r*r - w * w - h * h) / (2 * r*r));
			double am = (aw + ah + ac) / 2;
			double SR = 8 * atan(sqrt(tan(am / 2) * tan((am - aw) / 2) * tan((am - ah) / 2) * tan((am - ac) / 2)));
			return SR;
		};
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
		double orx = acos(OC.z / r), orz = atan2(OC.x, -OC.y); if (isnan(orz)) orz = 0;
		const double rx = atan2(sin(orx)*cos(rt), cos(orx)), ry = atan2(-sin(orx)*sin(rt), hypot(sin(orx)*cos(rt), cos(orx))),
			rz = atan2(sin(orz)*cos(rt) + sin(rt)*cos(orx)*cos(orz), cos(rt)*cos(orz) - sin(rt)*cos(orx)*sin(orz));		// first x, then y, finally z
		
		parallelogram CV(point(w / 2, -h / 2, rs), point(-w, 0), point(0, h));
		CV = matrix3D(Rotation, rx, ry, rz) * CV + C;


		this->resize();
		object3D* oi = Render_CheckContaining(C);

#ifdef _4_Threads_Rendering
		cout << "\nAttempting...";
#endif

#endif

		// calculate pixels allocates to each thread
		ray beg; beg.orig = C;
		rgblight c, s;
		fsec fs;
#ifndef DEBUG
#ifdef _4_Threads_Rendering
		const int STEP = 8;
		vector<double> attempt; attempt.resize(canvas.width() / STEP + 1);
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
						beg.dir = CV.O + u * CV.A + v * CV.B - beg.orig;
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
#else
		int NONE;
		this->MultiThread_CC(canvas, 0, canvas.width(), 0, canvas.height(), C, CV, oi, NONE);
		//ML.unlock();
#endif


#endif

		cout << " Completed. \n\n";
		auto Time_End = NTime::now();
		fs = Time_End - Time_Beg;
		cout << "Elapsed Time: " << defaultfloat << setprecision(3) << fs.count() << "s. \n\n\n";

		fout << "Elapsed Time: " << defaultfloat << setprecision(3) << fs.count() << "s. \n\n\n";


		// Debug single pixel
		int debugx = 320, debugy = 120, RP = 0;
		this->MultiThread_CC(canvas, debugx, debugx + 1, debugy, debugy + 1, C, CV, oi, RP);
		//pixel col = Red; canvas.dot(debugx + 1, debugy, col); canvas.dot(debugx - 1, debugy, col); canvas.dot(debugx, debugy + 1, col); canvas.dot(debugx, debugy - 1, col);

	}

	void print(ostream &os, unsigned n) {
		string space; for (unsigned i = 0; i < n; i++) space += " ";
		os << space << border << endl; space += " ";
		for (unsigned i = 0; i < Objs.size(); i++) os << space << *(Objs.at(i)) << endl;
		for (unsigned i = 0; i < GObjs.size(); i++) GObjs.at(i)->print(os, n + 1);
	}
	friend ostream& operator << (ostream& os, const World &W) {
		os << W.border << endl;
		for (unsigned i = 0; i < W.Objs.size(); i++) os << " " << *(W.Objs.at(i)) << endl;
		for (unsigned i = 0; i < W.GObjs.size(); i++) W.GObjs.at(i)->print(os, 1);
		os << defaultfloat;
		return os;
	}
};

unsigned World::Render_Sampling = 1;

// https://zhuanlan.zhihu.com/p/41269520

#endif
