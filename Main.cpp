using namespace std;

#include <iostream>

#include "Matrix.h"
#include "D:\Explore\Math\Graph\GraphFun\GraphFun\BitMap.h"


template<typename T>
void sort(T array[], int length) {
	int h = 1;
	while (h < length / 3) {
		h = 3 * h + 1;
	}
	while (h >= 1) {
		for (int i = h; i < length; i++) {
			for (int j = i; j >= h && array[j] < array[j - h]; j -= h) {
				std::swap(array[j], array[j - h]);
			}
		}
		h = h / 3;
	}
}
void elimination(matrix<double> &M) {
	const double Float_Limit = 1e-9;
	const unsigned h = *((unsigned*)&M);
	const unsigned w = *((unsigned*)&M + 1);
	double *data = &M[0][0];
	if (h == 0 || w == 0 || h == 1) return;
	if (w == 1) {
		for (int i = 0; i < h; i++) {
			*data = 0, data++;
		}
		data -= h;
		*data = 1;
		return;
	}
	double *p, *q, c;
	for (int i = 0, j = 0; i < w; i++) {
		while (abs(M[j][i]) < Float_Limit) {	// begin with 0
			if (abs(M[j][i]) < Float_Limit) M[j][i] = round(M[j][i]);
			for (int k = j + 1; k < h; k++) {
				if (abs(M[k][i]) > Float_Limit) {
					// Interchange(j, k)
					p = data + j * w + i, q = data + k * w + i;
					for (int n = i; n < w; n++) {
						c = *p, *p = *q, *q = c;
						p++, q++;
					}
					p = q = 0;
					//cout << "Interchange the " << (j + 1) << "th and the " << (k + 1) << "th row. \n" << M << endl;
					break;
				}
				else M[k][i] = round(M[k][i]);
				if (abs(M[k][i] - floor(M[k][i])) < Float_Limit) M[k][i] = round(M[k][i]);
				if (abs(M[k][i] - ceil(M[k][i])) < Float_Limit) M[k][i] = round(M[k][i]);
				if (k + 1 >= h) {	// whole colume is 0
					i++; break;
				}
			}
			if (i + 1 >= w) break;
			if (j + 1 >= h) break;
			if (abs(M[j][i]) < Float_Limit) M[j][i] = round(M[j][i]);
			if (abs(M[j][i] - floor(M[j][i])) < Float_Limit) M[j][i] = round(M[j][i]);
			if (abs(M[j][i] - ceil(M[j][i])) < Float_Limit) M[j][i] = round(M[j][i]);
		}
		for (int i = 0; i < h; i++) {
			for (int j = 0; j < w; j++) {
				if (abs(M[i][j]) < Float_Limit) M[i][j] = round(M[i][j]);
				if (abs(M[i][j] - floor(M[i][j])) < Float_Limit) M[i][j] = round(M[i][j]);
				if (abs(M[i][j] - ceil(M[i][j])) < Float_Limit) M[i][j] = round(M[i][j]);
			}
		}
		for (int k = j + 1; k < h; k++) {
			if (abs(M[k][i]) > Float_Limit) {
				// Multiple j and Add to k
				c = -M[k][i] / M[j][i];
				p = data + j * w + i, q = data + k * w + i;
				for (int n = i; n < w; n++) {
					*q += c * (*p);
					q++, p++;
				}
				p = q = 0;
				//cout << "Add " << c << " times the " << (j + 1) << "th row to the " << (k + 1) << "th row. \n" << M << endl;
			}
			else M[k][i] = round(M[k][i]);
			if (abs(M[k][i] - floor(M[k][i])) < Float_Limit) M[k][i] = round(M[k][i]);
			if (abs(M[k][i] - ceil(M[k][i])) < Float_Limit) M[k][i] = round(M[k][i]);
		}
		j++;
		if (j >= h) break;
	}
	for (int i = 0; i < h; i++) {
		for (int j = 0; j < w; j++) {
			if (abs(M[i][j]) < Float_Limit) M[i][j] = round(M[i][j]);
			if (abs(M[i][j] - floor(M[i][j])) < Float_Limit) M[i][j] = round(M[i][j]);
			if (abs(M[i][j] - ceil(M[i][j])) < Float_Limit) M[i][j] = round(M[i][j]);
		}
	}
	for (int i = h - 1; i > 0; i--) {
		p = data + i * w;
		int k;
		for (k = 0; k < w; k++) {
			if (abs(*p) > Float_Limit) break;
			else *p = 0;
			p++;
		}
		if (k != w) {
			for (int j = i - 1; j >= 0; j--) {
				if (abs(M[j][k]) > Float_Limit) {
					c = -(M[j][k] / M[i][k]);
					p = data + i * w + k, q = data + j * w + k;
					for (int n = k; n < w; n++) {
						*q += c * (*p);
						q++, p++;
					}
					p = q = 0;
					//cout << "Add " << c << " times the " << (i + 1) << "th row to the " << (j + 1) << "th row. \n" << M << endl;
				}
				else M[j][k] = round(M[j][k]);
				if (abs(M[j][k] - floor(M[j][k])) < Float_Limit) M[j][k] = round(M[j][k]);
				if (abs(M[j][k] - ceil(M[j][k])) < Float_Limit) M[j][k] = round(M[j][k]);
			}
		}
	}
	for (int i = 0; i < h; i++) {
		for (int j = 0; j < w; j++) {
			if (abs(M[i][j]) < Float_Limit) M[i][j] = round(M[i][j]);
			if (abs(M[i][j] - floor(M[i][j])) < Float_Limit) M[i][j] = round(M[i][j]);
			if (abs(M[i][j] - ceil(M[i][j])) < Float_Limit) M[i][j] = round(M[i][j]);
		}
	}
	p = data;
	for (int i = 0; i < h; i++) {
		c = 0;
		for (int j = 0; j < w; j++) {
			if (abs(*p) > Float_Limit) {
				if (abs(c) < Float_Limit) c = 0;
				if (abs(c) < Float_Limit) c = 1 / (*p), *p = 1;
				else (*p) *= c;
			}
			else *p = 0;
			p++;
		}
	}
	p = 0;
	for (int i = 0; i < h; i++) {
		for (int j = 0; j < w; j++) {
			if (abs(M[i][j]) < Float_Limit) M[i][j] = round(M[i][j]);
			if (abs(M[i][j] - floor(M[i][j])) < Float_Limit) M[i][j] = round(M[i][j]);
			if (abs(M[i][j] - ceil(M[i][j])) < Float_Limit) M[i][j] = round(M[i][j]);
		}
	}
}

inline bool baddouble(double a) {
	return ((unsigned long long)((*((unsigned long long*)&a)) & 0b0111111111110000000000000000000000000000000000000000000000000000)
		== 0b0111111111110000000000000000000000000000000000000000000000000000);
	// 0 11111111111 1000000000000000000000000000000000000000000000000000		nan
	// 1 11111111111 1000000000000000000000000000000000000000000000000000		-nan(ind)
	// 0 11111111111 0000000000000000000000000000000000000000000000000000		inf
	// 1 11111111111 0000000000000000000000000000000000000000000000000000		-inf
}

inline fraction randFrac() {
	return fraction(floor(abs(randnor(0, 3))), ceil(abs(randnor(0, 0.8)) + 0.000001));
}
matrix<fraction> randMatrix() {
	if (rand() % 3) {
		matrix<fraction> M(round(abs(randnor(4, 0.8))));
		for (int i = 0; i < M.height(); i++) {
			for (int j = 0; j < M.width(); j++) {
				M[i][j] = randFrac();
			}
		}
		return M;
	}
	else {
		matrix<fraction> M(round(abs(randnor(4, 0.7))), round(abs(randnor(4, 0.7))));
		for (int i = 0; i < M.height(); i++) {
			for (int j = 0; j < M.width(); j++) {
				M[i][j] = randFrac();
			}
		}
		return M;
	}
}

#include <fstream>
#include <cstdlib>
#include <ctime>

struct line {
	Vector<double> beg, end;
	pixel col;
};
inline line operator * (matrix<double> &M, line l) {
	l.beg = M * l.beg, l.end = M * l.end;
	return l;
}
inline void operator *= (line &l, double a) {
	l.beg *= a, l.end *= a;
}
inline line operator * (line l, double a) {
	l.beg *= a, l.end *= a;
	return l;
}
inline void operator /= (line &l, double a) {
	l.beg /= a, l.end /= a;
}
inline void operator += (line &l1, line &l2) {
	l1.beg += l2.beg, l1.end += l2.end;
}
inline line operator + (line l1, line &l2) {
	l1.beg += l2.beg, l1.end += l2.end;
	return l1;
}
inline line operator + (line l, Vector<double> v) {
	l.beg += v, l.end += v; return l;
}
inline void cvb(line &l) {
	l.beg.~Vector();
	l.end.~Vector();
	return;
}
inline bool nal(line &l) {
	return l.beg.dimension() == 0 || l.end.dimension() == 0;
}

void draw(bitmap &img, double d) {

}
inline void rotate(Vector<double> &V) {
	V *= 50;
	//V = Tr * V;
	V.x() += 100, V.y() += 100;
}

void SortVectorByZ(line *V, const int size) {
	if (size == 2) {
		if (V->end.z() > (V + 1)->end.z()) swap(*V, V[1]);
		return;
	}
	if (size == 3) {
		if (V->end.z() > (V + 1)->end.z()) swap(*V, V[1]);
		if ((V + 1)->end.z() > (V + 2)->end.z()) swap(V[1], V[2]);
		if (V->end.z() > (V + 1)->end.z()) swap(*V, V[1]);
		return;
	}
	int sp = size / 2, sq = size - sp;
	line *P, *Q;
	P = new line[sp];
	for (int i = 0; i < sp; i++) {
		*P = *V; P++, V++;
	}
	P -= sp;
	SortVectorByZ(P, sp);
	Q = new line[sq];
	for (int i = 0; i < sq; i++) {
		*Q = *V; Q++, V++;
	}
	Q -= sq;
	SortVectorByZ(Q, sq);
	V -= size;
	int mp = sp, mq = sq;
	while (mp > 0 && mq > 0) {
		if (P->end.z() <= Q->end.z()) *V = *P, P++, mp--;
		else *V = *Q, Q++, mq--;
		V++;
	}
	while (mp > 0) {
		*V = *P;
		V++, P++, mp--;
	}
	while (mq > 0) {
		*V = *Q;
		V++, Q++, mq--;
	}
	V -= size;
	P -= sp, Q -= sq;
	delete[] P, Q;
}

#include "Object.h"

inline double fx(double u, double v) {
	return cos(u)*(0.8 + 0.3*sin(v));
}
inline double fy(double u, double v) {
	return sin(u)*(0.8 + 0.3*sin(v));
}
inline double fz(double u, double v) {
	return 0.3*cos(v) + 0.3;
}

#include <crtdbg.h>
int main()
{
	//_CrtSetBreakAlloc(6718);

	/*ray ra;
	intersect Int;
	while (1) {
		cout << "Enter 3 vertexes of the triangle: " << endl;
		cout << "    A: "; cin >> tri.A.x >> tri.A.y >> tri.A.z;
		cout << "    B: "; cin >> tri.B.x >> tri.B.y >> tri.B.z;
		cout << "    C: "; cin >> tri.C.x >> tri.C.y >> tri.C.z;
		cout << "Enter the origin point of the ray: "; cin >> ra.orig.x >> ra.orig.y >> ra.orig.z;
		cout << "Enter another point that the ray passes through: "; cin >> ra.dir.x >> ra.dir.y >> ra.dir.z;
		ra.dir -= ra.orig;
		Int = tri.meet(ra);
		if (Int.meet) {
			cout << "The ray intersects the triangle at " << Int.intrs << " . " << endl;
			cout << "And the distance between the origin point and the intersection point is " << Int.dist << ". " << endl;
		}
		else {
			cout << "The ray doesn't intersect the triangle. " << endl;
		}
		cout << endl << endl;
	}*/

#define dat double

	bitmap img(600, 400);
	clock_t tm = clock();
	/*img.dot(100, 100, rgb(0, 255, 255));
	img.rect(120, 120, 120, 80, rgb(255, 255, 0));
	img.rect(500, 300, 400, 200, rgb(0, 128, 255));
	img.rect(160, 160, 100, 60, rgb(0, 184, 0), 0.4);
	img.rect(200, 150, -300, -200, rgb(128, 0, 64), 0.2);
	img.rect(-200, 60, 280, -200, rgb(255, 0, 128));
	img.line(140, 20, 140, 240, 3, rgb(255, 0, 0));
	img.line(30, 30, 240, 300, 3, rgb(0, 255, 0));
	img.circle(300, 200, 30, rgb(255, 184, 128));*/
	//img.out("D:\\Coding\\AboutMath\\SuperCalculator\\Image.bmp");
	/*bitmap image(600, 600, rgb(255, 255, 255));
	for (int j = 0; j < 5; j++) {
		double r1 = 50 * (j + 0.5);
		double r2 = 50 * (j + 1.5);
		double t = j * PI / 64.0;
		for (int i = 1; i <= 64; i++, t += PI / 32) {
			double ct = cos(t), st = sin(t);
			image.line(300 + r1 * ct, 300 - r1 * st, 300 + r2 * ct, 300 - r2 * st, j + 1, rgb(0, 0, 0));
		}
	}*/
	//image.out("D:\\Coding\\AboutMath\\SuperCalculator\\LineTest.bmp");
	img.clear();
	int NP = 5;
	//Vector<double> V[NP] = { {0,0,0}, {0,0,1}, {0,1,0}, {0,1,1}, {1,0,0}, {1,0,1}, {1,1,0}, {1,1,1} };	// Cube
	//Vector<double> V[NP] = { {0,0,0}, {0,1.4,0}, {1.4,0,0}, {1.4,1.4,0}, {0.7,0.7,1.1} };		// Pyramid
	/*NP = 540;
	Vector<double> *V;
	V = new Vector<double>[NP];
	double tnf;
	for (int i = 0; i <= 12; i++) {
		for (int j = -18; j < 18; j++) {
			tnf = _const_pi / 18 * j;
			*V = { (cos(tnf) + 1) * 0.7, (sin(tnf) + 1) * 0.7, i * 0.1 }, V++;
		}
	}
	for (int j = -18; j < 18; j++) {
		tnf = _const_pi / 18 * j;
		*V = { (cos(tnf) + 1) * 0.7, (sin(tnf) + 1) * 0.7, 0 }, V++;
		*V = { (cos(tnf) + 1) * 0.7, (sin(tnf) + 1) * 0.7, 1.2 }, V++;
	}
	V -= NP;*/
	// Cylinder
	/*NP = 578;
	Vector<double> *V;
	V = new Vector<double>[NP];
	cout << V << endl;
	double dir;
	*V = { 0.75,0.75,0 }, V++;
	*V = { 0.75,0.75,1.5 }, V++;
	for (int i = -18; i < 18; i++) {
		dir = _const_pi / 18 * i;
		*V = { (cos(dir) + 1) * 0.75, (sin(dir) + 1) * 0.75, 0 }, V++;
	}
	for (int i = 0; i < 15; i++) {
		for (int j = -18; j < 18; j++) {
			dir = _const_pi / 18 * j;
			*V = { 0.75 * (cos(dir) * (1 - i / 15.0) + 1), 0.75 * (sin(dir) * (1 - i / 15.0) + 1), i / 10.0 }, V++;
		}
	}
	V -= NP;*/
	// Cone
	/*NP = 614;
	Vector<double> *V = new Vector<double>[NP];
	cout << V << endl;
	double dir, ra, ht;
	*V = { 0.9, 0.9, 0 }, V++;
	*V = { 0.9, 0.9, 1.8 }, V++;
	for (int i = 1; i < 18; i++) {
		ht = sin((i - 9.0) / 18 * _const_pi) * 0.9 + 0.9;
		ra = 2 * sqrt(0.81 - pow((ht - 0.9), 2));
		for (int j = -18; j < 18; j++) {
			dir = _const_pi / 18 * j;
			*V = { 0.5 * (cos(dir) * ra + 1.8), 0.5 * (sin(dir) * ra + 1.8), ht }, V++;
		}
	}
	V -= NP;
	cout << V << endl;*/
	// Sphere
	/*NP = 2600;
	Vector<double> *V = new Vector<double>[NP];
	double px, py;
	for (int i = -5; i < 45; i++) {
		for (int j = -5; j < 45; j++) {
			px = i / 20.0, py = j / 20.0;
			*V = { px--, py--, sqrt(1 - px * px - py * py - abs(px)*py) };
			V++;
		}
	}
	for (int i = -5; i < 45; i++) {
		px = i / 20.0;
		*V = { px--, sqrt(1 - 0.75*px*px) - abs(px) / 2 + 1, 0 }, V++, px++;
		*V = { px--, -sqrt(1 - 0.75*px*px) - abs(px) / 2 + 1, 0 }, V++;
	}
	V -= NP;*/
	// Heart*
	/*NP = 3721;
	Vector<double> *V = new Vector<double>[NP];
	double kx, ky;
	for (int i = 0; i <= 60; i++) {
		kx = (i - 30.0) / 8;
		for (int j = 0; j <= 60; j++) {
			ky = (j - 30.0) / 8;
			*V = { kx, ky, 3 * (1 - kx)*(1 - kx) * exp(-kx * kx - (ky + 1)*(ky + 1))
				- 10 * (kx / 5 - kx * kx * kx - pow(ky, 5)) * exp(-kx * kx - ky * ky) - 1 / 3 * exp(-(kx + 1)*(kx + 1) - ky * ky) };
			V->z() /= 4; *V /= 1.9; V++;
		}
	}
	V -= NP;*/
	// MATLAB Peaks Function, rotate=(-1.3,0,-2.1)
	/*NP = 12322;
	line *V = new line[NP];
	double kx, ky;
	complex C;
	for (int i = -35; i <= 65; i++) {
		kx = (i - 30.0) / 8;
		for (int j = 0; j <= 60; j++) {
			ky = (j - 30.0) / 8;
			C = tgamma(complex(kx, ky));
			V->end = { kx, ky, abs(C) };
			V->end.z() /= 2;
			*V /= 2;
			V->col = fromHSL(arg(C) * 0.15915494309189533577, 1, 1 - pow(0.5, abs(C)));
			if (j != 0) V->beg = (V - 1)->end;
			else V->beg = V->end;
			V++;
		}
	}
	for (int i = 0; i <= 60; i++) {
		ky = (i - 30.0) / 8;
		for (int j = -35; j <= 65; j++) {
			kx = (j - 30.0) / 8;
			C = tgamma(complex(kx, ky));
			V->end = { kx, ky, abs(C) };
			V->end.z() /= 2;
			*V /= 2;
			V->col = fromHSL(arg(C) * 0.15915494309189533577, 1, 1 - pow(0.5, abs(C)));
			if (j != -35) V->beg = (V - 1)->end;
			else V->beg = V->end;
			V++;
		}
	}
	V -= NP;*/
	// Gamma Function
	/*NP = 1296;
	line *V = new line[NP];
	double u, v;
	for (int i = 0; i < 36; i++) {
		u = i * PI / 18;
		for (int j = -9; j <= 9; j++) {
			v = j * PI / 18;
			V->end = { 0.8*(sin(u)*cos(v) + 1), cos(u)*cos(v) + 1, 0.6*(sin(v) + 1) };
			if (j != -9) V->beg = (V - 1)->end;
			else V->beg = V->end;
			V->col = rgb(255, 255, 255);
			V++;
		}
	}
	for (int j = -8; j < 9; j++) {
		v = j * PI / 18;
		for (int i = 0; i < 36; i++) {
			u = i * PI / 18;
			V->end = { 0.8*(sin(u)*cos(v) + 1), cos(u)*cos(v) + 1, 0.6*(sin(v) + 1) };
			if (i != 0) V->beg = (V - 1)->end;
			if (i == 35) (V - 35)->beg = V->end;
			V->col = rgb(255, 255, 255);
			V++;
		}
	}
	V -= NP;*/
	// Ellipsoid
	/*NP = 1296;
	line *V = new line[NP];
	double u, v;
	for (int i = 0; i < 36; i++) {
		u = i * PI / 18;
		for (int j = 0; j < 18; j++) {
			v = j * PI / 9;
			V->end = { cos(u)*(0.8 + 0.3*sin(v)) + 0.8, sin(u)*(0.8 + 0.3*sin(v)) + 0.8, 0.3*cos(v) + 0.3 };
			if (j != 0) V->beg = (V - 1)->end;
			if (j == 17) (V - 17)->beg = V->end;
			V->col = rgb(255, 255, 255);
			V++;
		}
	}
	for (int j = 0; j < 18; j++) {
		v = j * PI / 9;
		for (int i = 0; i < 36; i++) {
			u = i * PI / 18;
			V->end = { cos(u)*(0.8 + 0.3*sin(v)) + 0.8, sin(u)*(0.8 + 0.3*sin(v)) + 0.8, 0.3*cos(v) + 0.3 };
			if (i != 0) V->beg = (V - 1)->end;
			if (i == 35) (V - 35)->beg = V->end;
			V->col = rgb(255, 255, 255);
			V++;
		}
	}
	V -= NP;*/
	// Ring
	/*NP = 30;
	line *V = new line[NP];
	Vector<double> P[12];
	const double p = (sqrt(5) - 1) / 2;
	P[0] = { 1,p,0 }, P[1] = { 1,-p,0 }, P[2] = { -1,p,0 }, P[3] = { -1,-p,0 };
	P[4] = { p,0,1 }, P[5] = { p,0,-1 }, P[6] = { -p,0,1 }, P[7] = { -p,0,-1 };
	P[8] = { 0,1,p }, P[9] = { 0,1,-p }, P[10] = { 0,-1,p }, P[11] = { 0,-1,-p };
	for (int i = 0; i < 12; i++) {
		for (int j = i + 1; j < 12; j++) {
			if (abs(mod(P[j] - P[i]) - 2 * p) < 0.01) {
				V->beg = P[i], V->end = P[j];
				V->col = rgb(255, 255, 255); V++;
			}
		}
	}
	V -= NP;
	double rt = - PI / 10;
	matrix<double> Xr = { {1,0,0}, {0,cos(rt),-sin(rt)}, {0,sin(rt),cos(rt)} };
	rt = PI / 6;
	matrix<double> Zr = { {cos(rt),-sin(rt),0}, {sin(rt),cos(rt),0}, {0,0,1} };
	for (int i = 0; i < NP; i++) V[i] = Xr * V[i];
	for (int i = 0; i < NP; i++) V[i] = Zr * V[i];
	for (int i = 0; i < NP; i++) V[i].beg = V[i].beg * 0.7 + Vector<double>({1, 1, 1}), V[i].end = V[i].end * 0.7 + Vector<double>({ 1, 1, 1 });*/
	// Icosahedron*

	NP = 12;
	/*line *V = new line[NP];
	V[0].beg = { 0,0,0 }, V[0].end = { 0,1,0 };
	V[1].beg = { 0,1,0 }, V[1].end = { 1,1,0 };
	V[2].beg = { 1,1,0 }, V[2].end = { 1,0,0 };
	V[3].beg = { 1,0,0 }, V[3].end = { 0,0,0 };
	for (int i = 4; i < 8; i++) {
		V[i] = V[i - 4]; V[i].beg.z() = V[i].end.z() = 1;
	}
	for (int i = 8; i < 12; i++) {
		V[i].beg = V[i - 8].beg, V[i].end = V[i - 4].beg;
	}
	for (int i = 0; i < 12; i++) V[i].col = rgb(255, 255, 255);*/


#ifdef RealN
	Vector<double> X(3, 0, 0), Y(0, 3, 0), Z(0, 0, 2);
	double lx = 0.1, ly = 0.1, lz = 0.1;
	double mols = sqrt(lx * lx + ly * ly + lz * lz);
	//double rx = -acos(lx / mols), ry = -acos(ly / mols), rz = -acos(lz / mols);
	//double rx = -0.645771823, ry = 0, rz = -0.523598776;
	double rx = -1.1, ry = 0, rz = -2.1;
	//double rx = 0, ry = 0, rz = 0;
	matrix<double> Rx = { {1,0,0}, {0,cos(rx),-sin(rx)}, {0,sin(rx),cos(rx)} };
	matrix<double> Ry = { {cos(ry),0,sin(ry)}, {0,1,0}, {-sin(ry),0,cos(ry)} };
	matrix<double> Rz = { {cos(rz),-sin(rz),0}, {sin(rz),cos(rz),0}, {0,0,1} };
	matrix<double> Tr = Rx * Ry * Rz;
	for (int i = 0; i < NP; i++) V[i] = Tr * V[i];
	X = Tr * X, Y = Tr * Y, Z = Tr * Z;
	double tms = 100;
	for (int i = 0; i < NP; i++) V[i] *= tms;
	X *= tms, Y *= tms, Z *= tms;
	Vector<double> add(300, 200, 0);
	for (int i = 0; i < NP; i++) V[i] += add;
	X += add, Y += add, Z += add;
	img.line(add.x(), add.y(), X.x(), X.y(), 3, rgb(255, 0, 0));
	img.line(add.x(), add.y(), Y.x(), Y.y(), 3, rgb(0, 255, 0));
	img.line(add.x(), add.y(), Z.x(), Z.y(), 3, rgb(0, 0, 255));
#endif

#ifdef Parallel_Proj
	line X, Y, Z;
	X.beg = Y.beg = Z.beg = { 0, 0, 0 };
	X.end = { 3, 0, 0 }, Y.end = { 0, 3, 0 }, Z.end = { 0, 0, 2 };
	double lx = 0.1, ly = 0.1, lz = 0.1;
	double mols = sqrt(lx * lx + ly * ly + lz * lz);
	//double rx = -acos(lx / mols), ry = -acos(ly / mols), rz = -acos(lz / mols);
	//double rx = -0.645771823, ry = 0, rz = -0.523598776;
	double rx = -1.1, ry = 0, rz = -2.1;
	//double rx = 0, ry = 0, rz = 0;
	matrix<double> Rx = { {1,0,0}, {0,cos(rx),-sin(rx)}, {0,sin(rx),cos(rx)} };
	matrix<double> Ry = { {cos(ry),0,sin(ry)}, {0,1,0}, {-sin(ry),0,cos(ry)} };
	matrix<double> Rz = { {cos(rz),-sin(rz),0}, {sin(rz),cos(rz),0}, {0,0,1} };
	matrix<double> Tr = Rx * Ry * Rz;
	for (int i = 0; i < NP; i++) V[i] = Tr * V[i];
	X = Tr * X, Y = Tr * Y, Z = Tr * Z;
	double tms = 100;
	for (int i = 0; i < NP; i++) V[i] *= tms;
	X *= tms, Y *= tms, Z *= tms;
	line add;
	add.beg = add.end = { 300, 200, 0 };
	for (int i = 0; i < NP; i++) V[i] += add;
	X += add, Y += add, Z += add;
	img.line(add.end.x(), add.end.y(), X.end.x(), X.end.y(), 3, rgb(255, 0, 0));
	img.line(add.end.x(), add.end.y(), Y.end.x(), Y.end.y(), 3, rgb(0, 255, 0));
	img.line(add.end.x(), add.end.y(), Z.end.x(), Z.end.y(), 3, rgb(0, 0, 255));
	//img.clear();
	SortVectorByZ(V, NP);
	for (int i = 0; i < NP; i++) {
		if (!baddouble(V->beg.x()) && !baddouble(V->beg.y()) && !baddouble(V->beg.z())
			&& !baddouble(V->end.x()) && !baddouble(V->end.y()) && !baddouble(V->end.z())) {
			img.line(V->beg.x(), V->beg.y(), V->end.x(), V->end.y(), 3, V->col);
		}
		V++;
	}
	V -= NP;
	delete[] V;
#endif

	//Vector<double> X, Y, Z, O;
	//O = { 0, 0, 0, 1 }, X = { 3, 0, 0, 1 }, Y = { 0, 3, 0, 1 }, Z = { 0, 0, 2, 1 };
	double wx = 30, wy = 60, wz = 90;
	double rx = -1.1, rz = -2.1;
	double cx = 300, cy = 200, tms = 100;
	matrix<double> T = {
		{cos(rz), -sin(rz), 0, cx / tms},
		{cos(rx)*sin(rz), cos(rx)*cos(rz), -sin(rx), cy / tms},
		{sin(rx)*sin(rz), sin(rx)*cos(rz), cos(rx), 0},
		{0, 0, 0, 1 / tms}
	};
	/*T *= matrix<double>({
		{1, 0, 0, 0},
		{0, 1, 0, 0},
		{0, 0, 1, 0},
		{-1 / wx, -1 / wy, -1 / wz, 1}
		});
	Vector<double> D4;
	for (int i = 0; i < NP; i++) {
		D4 = { V[i].beg.x(), V[i].beg.y(), V[i].beg.z(), 1 };
		D4 = T * D4;
		D4 /= D4.at(4);
		V[i].beg = { D4.x(), D4.y(), D4.z() };
		D4 = { V[i].end.x(), V[i].end.y(), V[i].end.z(), 1 };
		D4 = T * D4;
		D4 /= D4.at(4);
		V[i].end = { D4.x(), D4.y(), D4.z() };
	}
	X = T * X, Y = T * Y, Z = T * Z, O = T * O;
	X /= X.at(4), Y /= Y.at(4), Z /= Z.at(4), O /= O.at(4);
	img.line(O.x(), O.y(), X.x(), X.y(), 3, rgb(255, 0, 0));
	img.line(O.x(), O.y(), Y.x(), Y.y(), 3, rgb(0, 255, 0));
	img.line(O.x(), O.y(), Z.x(), Z.y(), 3, rgb(0, 0, 255));
	SortVectorByZ(V, NP);
	for (int i = 0; i < NP; i++) {
		if (!baddouble(V->beg.x()) && !baddouble(V->beg.y()) && !baddouble(V->beg.z())
			&& !baddouble(V->end.x()) && !baddouble(V->end.y()) && !baddouble(V->end.z())) {
			img.line(V->beg.x(), V->beg.y(), V->end.x(), V->end.y(), 3, V->col);
		}
		V++;
	}
	V -= NP;
	delete[] V;*/

	// Cube
	/*img.line(V[0].x(), V[0].y(), V[1].x(), V[1].y(), 3, rgb(255, 255, 255));
	img.line(V[0].x(), V[0].y(), V[2].x(), V[2].y(), 3, rgb(255, 255, 255));
	img.line(V[1].x(), V[1].y(), V[3].x(), V[3].y(), 3, rgb(255, 255, 255));
	img.line(V[2].x(), V[2].y(), V[3].x(), V[3].y(), 3, rgb(255, 255, 255));
	img.line(V[0].x(), V[0].y(), V[4].x(), V[4].y(), 3, rgb(255, 255, 255));
	img.line(V[1].x(), V[1].y(), V[5].x(), V[5].y(), 3, rgb(255, 255, 255));
	img.line(V[2].x(), V[2].y(), V[6].x(), V[6].y(), 3, rgb(255, 255, 255));
	img.line(V[3].x(), V[3].y(), V[7].x(), V[7].y(), 3, rgb(255, 255, 255));
	img.line(V[4].x(), V[4].y(), V[5].x(), V[5].y(), 3, rgb(255, 255, 255));
	img.line(V[4].x(), V[4].y(), V[6].x(), V[6].y(), 3, rgb(255, 255, 255));
	img.line(V[5].x(), V[5].y(), V[7].x(), V[7].y(), 3, rgb(255, 255, 255));
	img.line(V[6].x(), V[6].y(), V[7].x(), V[7].y(), 3, rgb(255, 255, 255));*/
	// Pyramid
	/*img.line(V[0].x(), V[0].y(), V[1].x(), V[1].y(), 3, rgb(255, 255, 255));
	img.line(V[0].x(), V[0].y(), V[2].x(), V[2].y(), 3, rgb(255, 255, 255));
	img.line(V[1].x(), V[1].y(), V[3].x(), V[3].y(), 3, rgb(255, 255, 255));
	img.line(V[2].x(), V[2].y(), V[3].x(), V[3].y(), 3, rgb(255, 255, 255));
	for (int i = 0; i < 4; i++) {
		img.line(V[i].x(), V[i].y(), V[4].x(), V[4].y(), 3, rgb(255, 255, 255));
	}*/
	// Cylinder
	/*for (int i = 0; i <= 12; i++) {
		for (int j = 0; j < 35; j++) {
			img.line(V->x(), V->y(), (V + 1)->x(), (V + 1)->y(), 3, rgb(255, 255, 255)), V++;
		}
		img.line(V->x(), V->y(), (V - 35)->x(), (V - 35)->y(), 3, rgb(255, 255, 255)), V++;
	}
	for (int j = 0; j < 36; j++) {
		img.line(V->x(), V->y(), (V + 1)->x(), (V + 1)->y(), 3, rgb(255, 255, 255)), V += 2;
	}
	V -= NP;
	delete[] V;*/
	// Cone
	/*Vector<double> VK = *V; V++;
	Vector<double> VL = *V; V++;
	for (int i = 0; i < 36; i++) {
		img.line(V->x(), V->y(), VK.x(), VK.y(), 3, rgb(255, 255, 255));
		img.line(V->x(), V->y(), VL.x(), VL.y(), 3, rgb(255, 255, 255));
		V++;
	}
	for (int i = 0; i < 15; i++) {
		for (int i = 0; i < 35; i++) {
			img.line(V->x(), V->y(), (V + 1)->x(), (V + 1)->y(), 3, rgb(255, 255, 255)), V++;
		}
		img.line(V->x(), V->y(), (V - 35)->x(), (V - 35)->y(), 3, rgb(255, 255, 255)), V++;
	}
	V -= NP;
	cout << V << endl;
	delete[] V;*/
	// Sphere
	/*Vector<double> sp = *V; V++;
	Vector<double> np = *V; V++;
	for (int i = 0; i < 17; i++) {
		for (int j = 0; j < 35; j++) {
			img.line(V->x(), V->y(), (V + 1)->x(), (V + 1)->y(), 3, rgb(255, 255, 255)), V++;
		}
		img.line(V->x(), V->y(), (V - 35)->x(), (V - 35)->y(), 3, rgb(255, 255, 255)), V++;
	}
	V -= NP - 2;
	for (int i = 0; i < 36; i++) {
		img.line(sp.x(), sp.y(), V->x(), V->y(), 3, rgb(255, 255, 255));
		for (int i = 0; i < 16; i++) {
			img.line(V->x(), V->y(), (V + 36)->x(), (V + 36)->y(), 3, rgb(255, 255, 255)); V += 36;
		}
		img.line(np.x(), np.y(), V->x(), V->y(), 3, rgb(255, 255, 255));
		V -= 575;
	}
	V -= 38;
	cout << V << endl;
	delete[] V;*/
	// Heart(Infinished)
	/*cout << V << endl;
	for (int i = 0; i < 50; i++) {
		for (int j = 0; j < 49; j++) {
			if (!baddouble(V->x()) && !baddouble(V->y()) && (V->x() >= 0 && V->x() < 600 && V->y() >= 0 && V->y() < 400)
				&& ((V + 1)->x() >= 0 && (V + 1)->x() < 600 && (V + 1)->y() >= 0 && (V + 1)->y() < 400)) {
				img.line(V->x(), V->y(), (V + 1)->x(), (V + 1)->y(), 3, rgb(255, 255, 255));
				//img.dot(V->x(), V->y(), rgb(255, 255, 255));
			}
			V++;
		}
		V++;
	}
	V -= NP - 100;
	for (int i = 0; i < 50; i++) {
		for (int j = 0; j < 49; j++) {
			if (!baddouble(V->x()) && !baddouble(V->y()) && (V->x() >= 0 && V->x() < 600 && V->y() >= 0 && V->y() < 400)
				&& ((V + 50)->x() >= 0 && (V + 50)->x() < 600 && (V + 50)->y() >= 0 && (V + 50)->y() < 400)) {
				img.line(V->x(), V->y(), (V + 50)->x(), (V + 50)->y(), 3, rgb(255, 255, 255));
				//img.dot(V->x(), V->y(), rgb(255, 255, 255));
			}
			V++;
		}
		V++;
	}
	for (int i = 0; i < 49; i++) {
		if (!baddouble(V->x()) && !baddouble(V->y()) && (V->x() >= 0 && V->x() < 600 && V->y() >= 0 && V->y() < 400)
			&& ((V + 2)->x() >= 0 && (V + 2)->x() < 600 && (V + 2)->y() >= 0 && (V + 2)->y() < 400)) {
			img.line(V->x(), V->y(), (V + 2)->x(), (V + 2)->y(), 3, rgb(255, 255, 255));
		}
		V += 2;
	}
	V -= 97;
	for (int i = 0; i < 49; i++) {
		if (!baddouble(V->x()) && !baddouble(V->y()) && (V->x() >= 0 && V->x() < 600 && V->y() >= 0 && V->y() < 400)
			&& ((V + 2)->x() >= 0 && (V + 2)->x() < 600 && (V + 2)->y() >= 0 && (V + 2)->y() < 400)) {
			img.line(V->x(), V->y(), (V + 2)->x(), (V + 2)->y(), 3, rgb(255, 255, 255));
		}
		V += 2;
	}
	V -= NP - 1;
	cout << V << endl;
	delete[] V;*/
	// MATLAB Peaks Function
	/*img.clear();
	for (int i = 0; i <= 60; i++) {
		for (int j = 0; j < 60; j++) {
			img.line(V->x(), V->y(), (V + 1)->x(), (V + 1)->y(), 3, rgb(0, 255, 255)), V++;
		}
		V++;
	}
	V -= NP;
	for (int i = 0; i < 60; i++) {
		for (int j = 0; j <= 60; j++) {
			img.line(V->x(), V->y(), (V + 61)->x(), (V + 61)->y(), 3, rgb(0, 255, 255)), V++;
		}
	}
	V -= NP - 61;
	delete[] V;*/

	//img.out("D:\\Coding\\AboutMath\\SuperCalculator\\Image+.bmp");

#define WORLD 6

	T = matrix<double>({ {1,0,0}, {0,cos(rx),-sin(rx)}, {0,sin(rx),cos(rx)} })
		* matrix<double>({ {cos(rz),-sin(rz),0}, {sin(rz),cos(rz),0}, {0,0,1} }) * tms;
	img.clear();

	tm = clock();

#if WORLD == 1

	triangle G1_tri{ {1,2,-0.5}, {-1,2,1}, {-0.8,1.4,0.7} };
	triangle G1_triXOY({ 0,0,0 }, { 3,0,0 }, { 0,3,0 });
	triangle G1_triXOZ({ 0,0,0 }, { 3,0,0 }, { 0,0,2 });
	triangle G1_triYOZ({ 0,0,0 }, { 0,3,0 }, { 0,0,2 });
	G1_tri.absorb = { 1,0.5,0 }, G1_tri.reflect = { 0,0.5,1 };
	G1_triXOY.absorb = { 0.5,0.5,0.5 }, G1_triXOY.reflect = { 0.5,0.5,0.5 };
	G1_triXOZ.absorb = { 0.8,0.8,0.8 }, G1_triXOZ.reflect = { 0.2,0.2,0.2 };
	G1_triYOZ.absorb = { 0.7,0.7,0.7 }, G1_triYOZ.reflect = { 0.3,0.3,0.3 };
	triangle G1_trix({ 0,-1,0 }, { 1.5,0.2,0.3 }, { 0.5,1,1 });
	G1_trix.absorb = { 0.5,1,0.75 }, G1_trix.reflect = { 0.5,0,0.25 };
	World G1; G1.add({ &G1_tri, &G1_trix, &G1_triXOY, &G1_triXOZ, &G1_triYOZ }); G1.background = rgblight(0.05, 0.08, 0.1);
	G1.render(img, point(30, 60, 90), -1.1, -2.1, 300, 200, 100);

#elif WORLD == 2

	point G2_V1(0, 0, 0.9), G2_V2(1, 0, 0.9), G2_V3(1.0 / 2, sqrt(3) / 2, 0.9), G2_V0(1.0 / 2, sqrt(3) / 6, 0.9 - sqrt(6) / 3);
	triangle G2_tri0(G2_V1, G2_V2, G2_V3);
	triangle G2_tri1(G2_V1, G2_V2, G2_V0), G2_tri2(G2_V2, G2_V3, G2_V0), G2_tri3(G2_V3, G2_V1, G2_V0);
	triangle G2_trib(point(0.0, sqrt(3) / 3), point(1.0, sqrt(3) / 3), point(1.0 / 2, -sqrt(3) / 6));
	parallelogram G2_prlb(point(-0.4, -0.4, 0), point(0, 1.4, 0), point(1.8, 0, 0));
	G2_tri1.absorb = { 0.5,0.5,0.5 }, G2_tri1.reflect = { 0.5,0.5,0.5 };
	G2_tri2.absorb = { 0.8,0.8,0.8 }, G2_tri2.reflect = { 0.2,0.2,0.2 };
	G2_tri3.absorb = { 0.7,0.7,0.7 }, G2_tri3.reflect = { 0.3,0.3,0.3 };
	G2_tri0.absorb = { 0.4,0.4,0.4 }, G2_tri0.reflect = { 0.6,0.6,0.6 };
	G2_prlb.absorb = { 0.6,0.5,0.4 }, G2_prlb.reflect = { 0.4,0.5,0.6 };
	World G2; G2.add({ &G2_tri0, &G2_tri1, &G2_tri2, &G2_tri3, &G2_prlb }); G2.background = rgblight(0.05, 0.08, 0.1);
	G2.render(img, point(30, 60, 90), -1.1, -2.1, 300, 200, 200);

#elif WORLD == 3

	World G3;
	parallelogram G3_U(point(0, 0, 0), point(1, 0, 0), point(0, 1, 0), 1);
	parallelogram G3_T(point(0, 0, 1), point(1, 0, 1), point(0, 1, 1), 1);
	parallelogram G3_L(point(0, 0, 0), point(0, 0, 1), point(1, 0, 0), 1);
	parallelogram G3_R(point(0, 1, 0), point(0, 1, 1), point(1, 1, 0), 1);
	parallelogram G3_B(point(0, 0, 0), point(0, 0, 1), point(0, 1, 0), 1);
	parallelogram G3_F(point(1, 0, 0), point(1, 0, 1), point(1, 1, 0), 1);
	G3_U.reflect = G3_T.reflect = G3_L.reflect = G3_R.reflect = G3_B.reflect = G3_F.reflect = rgblight(0.6, 0.6, 0.6);
	G3_U.absorb = G3_T.absorb = G3_L.absorb = G3_R.absorb = G3_B.absorb = G3_F.absorb = rgblight(0.4, 0.4, 0.4);
	G3.add({ &G3_U, &G3_T, &G3_L, &G3_R, &G3_B, &G3_F });
	triangle G3_xOy(point(0, 0, 0), point(0, 3, 0), point(3, 0, 0));
	triangle G3_xOz(point(0, 0, 0), point(3, 0, 0), point(0, 0, 2));
	triangle G3_yOz(point(0, 0, 0), point(0, 3, 0), point(0, 0, 2));
	G3_xOy.reflect = G3_xOz.reflect = G3_yOz.reflect = rgblight(0.4, 0.4, 0.7);
	G3_xOy.absorb = G3_xOz.absorb = G3_yOz.absorb = rgblight(0.6, 0.6, 0.3);
	//G3.add({ &G3_xOy, &G3_xOz, &G3_yOz });
	parallelogram G3_LD(point(-0.5, -0.5, 0), point(2, -0.5, 0), point(-0.5, 2, 0), 1);
	G3_LD.reflect = rgblight(0.4, 0.4, 0.7), G3_LD.absorb = rgblight(0.6, 0.6, 0.3);
	triangle G3_LDP(point(0.9, 1.4), point(1.5, 0.4), point(1.6, 1.9)); G3_LDP += point(0, 0, 1e-4);
	G3_LDP.absorb = { 0.5, 1, 0.75 }, G3_LDP.reflect = { 0.5, 0, 0.25 };
	G3.add({ &G3_LD, &G3_LDP });
	G3.background = rgblight(0.05, 0.08, 0.1);
	G3_U += point(0, 0, 0.1), G3_T += point(0, 0, 0.1), G3_L += point(0, 0, 0.1), G3_R += point(0, 0, 0.1), G3_B += point(0, 0, 0.1), G3_F += point(0, 0, 0.1);
	G3.render(img, point(30, 60, 90), -1.1, -2.1, 300, 200, 120);

#elif WORLD == 4

#define DPS 3
#define DIF1 8
#define DIF2 10
#define DIF3 32
	World G4, G4_R;
#if DPS==3
	vector<vector<vector<World*>>> G4_RS; G4_RS.resize(DIF1);
	for (int i = 0; i < DIF1; i++) {
		G4_RS.at(i).resize(DIF2);
		for (int j = 0; j < DIF2; j++) {
			G4_RS.at(i).at(j).resize(DIF3);
			for (int k = 0; k < DIF3; k++) {
				G4_RS.at(i).at(j).at(k) = new World;
			}
		}
	}
#elif DPS==2
	vector<vector<World*>> G4_RS; G4_RS.resize(DIF1);
	for (int i = 0; i < DIF1; i++) {
		G4_RS.at(i).resize(DIF2);
		for (int j = 0; j < DIF2; j++) {
			G4_RS.at(i).at(j) = new World;
		}
	}
#elif DPS==1
	vector<World*> G4_RS; G4_RS.resize(DIF1);
	for (int i = 0; i < DIF1; i++) {
		G4_RS.at(i) = new World;
	}
#endif
	triangle *W; NP = 409600;
	W = new triangle[NP];
	point C00, C01, C10, C11; double u, v, un, vn;
	for (int i = 0; i < 640; i++) {
		u = i * PI / 320;
		un = (i + 1) * PI / 320;
		for (int j = 0; j < 320; j++) {
			v = j * PI / 160;
			vn = (j + 1) * PI / 160;
			C00 = { fx(u,v),fy(u,v),fz(u,v) }, C01 = { fx(u,vn),fy(u,vn),fz(u,vn) },
				C10 = { fx(un,v),fy(un,v),fz(un,v) }, C11 = { fx(un,vn),fy(un,vn),fz(un,vn) };
			/*
			Parametrics:
			inline double fx(double u, double v) {
				return cos(u)*(0.8 + 0.3*sin(v));
			}
			inline double fy(double u, double v) {
				return sin(u)*(0.8 + 0.3*sin(v));
			}
			inline double fz(double u, double v) {
				return 0.3*cos(v) + 0.3;
			}
			Result would be a ring.
			*/
			*W = { C00, C10, C11 }, W++;
			*W = { C00, C01, C11 }, W++;
		}
	}
	W -= NP;
#if DPS==3
	int DF = NP / (DIF1*DIF2*DIF3);
	for (int i = 0; i < DIF1; i++) {
		G4_R.insert(new World);
		for (int j = 0; j < DIF2; j++) {
			G4_R.GObjs.back()->insert(new World);
			for (int k = 0; k < DIF3; k++) {
				for (int l = 0; l < DF; l++) {
					G4_RS.at(i).at(j).at(k)->add(W);
					W++;
				}
				G4_R.GObjs.back()->GObjs.back()->insert(G4_RS.at(i).at(j).at(k));
			}
		}
	}
#elif DPS==2
	int DF = NP / (DIF1*DIF2);
	for (int i = 0; i < DIF1; i++) {
		G4_R.insert(new World);
		for (int j = 0; j < DIF2; j++) {
			G4_R.GObjs.back()->insert(new World);
			for (int l = 0; l < DF; l++) {
				G4_RS.at(i).at(j)->add(W);
				W++;
			}
			G4_R.GObjs.back()->insert(G4_RS.at(i).at(j));
		}
	}
#elif DPS==1
	int DF = NP / (DIF1);
	for (int i = 0; i < DIF1; i++) {
		G4_R.insert(new World);
		for (int l = 0; l < DF; l++) {
			G4_R.GObjs.back()->add(W);
			W++;
		}
	}
#endif
	W -= NP;
	parallelogram G4_LD({ -2,-2,0 }, { 2,-2,0 }, { -2,2,0 }, 1);
	//G4_LD.reflect = rgblight(0.4, 0.4, 0.7), G4_LD.absorb = rgblight(0.6, 0.6, 0.3);	// gray-blue
	//G4_LD.reflect = rgb(240, 230, 140), G4_LD.absorb = rgb(15, 115, 25);	// pink
	G4_LD.reflect = rgb(255, 255, 100), G4_LD.absorb = ~G4_LD.reflect;	// bright yellow
	//G4_LD.reflect = rgb(178, 102, 255), G4_LD.absorb = ~G4_LD.reflect;	// mauve
	for (int i = 0; i < NP; i++) {
		//W->reflect = rgblight(0.8, 0.8, 0.8), W->absorb = rgblight(0.2, 0.2, 0.2);	// silver
		//W->reflect = rgblight(rgb(218, 165, 32)), W->absorb = rgb(15, 115, 25);	// gold
		W->reflect = rgblight(rgb(50, 255, 153)), W->absorb = ~W->reflect;	// cyan
		//W->reflect = rgblight(rgb(51, 153, 255)), W->absorb = ~W->reflect;	// blue
		W++;
	}
	W -= NP;
	G4.add(&G4_LD);
	G4.insert(&G4_R);

	G4.render(img, point(30, 60, 90), PI / 2, 0, -1.1, -2.1, 300, 200, 100);
	//G4.render(img, point(30, 60, 90), 1.7, 0.8, -1.1, -2.1, 300, 200, 100);


#if DPS==3
	for (int i = 0; i < DIF1; i++) {
		for (int j = 0; j < DIF2; j++) {
			for (int k = 0; k < DIF3; k++) {
				delete G4_RS.at(i).at(j).at(k);
			}
			//G4_RS.at(j).clear();
		}
		G4_RS.at(i).clear();
	}
#elif DPS==2
	for (int i = 0; i < DIF1; i++) {
		for (int j = 0; j < DIF2; j++) {
			delete G4_RS.at(i).at(j);
		}
		G4_RS.at(i).clear();
	}
#elif DPS==1
	for (int i = 0; i < DIF1; i++) {
		delete G4_RS.at(i);
	}

#endif
	delete[] W;

#elif WORLD == 5
#define DPS 3
#define DIF1 8
#define DIF2 16
	//#define DIF1 4
	//#define DIF2 8
#define DIF3 16
	World G5, G5_R;
#if DPS==3
	vector<vector<vector<World*>>> G5_RS; G5_RS.resize(DIF1);
	for (int i = 0; i < DIF1; i++) {
		G5_RS.at(i).resize(DIF2);
		for (int j = 0; j < DIF2; j++) {
			G5_RS.at(i).at(j).resize(DIF3);
			for (int k = 0; k < DIF3; k++) {
				G5_RS.at(i).at(j).at(k) = new World;
			}
		}
	}
#elif DPS==2
	vector<vector<World*>> G5_RS; G5_RS.resize(DIF1);
	for (int i = 0; i < DIF1; i++) {
		G5_RS.at(i).resize(DIF2);
		for (int j = 0; j < DIF2; j++) {
			G5_RS.at(i).at(j) = new World;
		}
	}
#elif DPS==1
	vector<World*> G5_RS; G5_RS.resize(DIF1);
	for (int i = 0; i < DIF1; i++) {
		G5_RS.at(i) = new World;
	}
#endif

#define Min_X -8.125
#define Max_X 4.375
#define Min_Y -3.75
#define Max_Y 3.75
#define X_Dif 192
#define Y_Dif 192
#define XC_Dif 4
#define YC_Dif 4
	//#define X_Dif 12
	//#define Y_Dif 8
	//#define XC_Dif 1
	//#define YC_Dif 1

	NP = X_Dif * Y_Dif;
	if ((DPS == 1 && (NP % (DIF1) != 0)) || (DPS == 2 && (NP % (DIF1 * DIF2) != 0)) || (DPS == 3 && (NP % (DIF1 * DIF2 * DIF3) != 0)))
		cout << "\aWarning: Number of triangles is not divisible by the number of blocks! \n\n", exit(-1);
	/*if ((DPS == 1 && (XC_Dif * YC_Dif * DIF1 != NP)) || (DPS == 2 && (XC_Dif * YC_Dif * DIF1 * DIF2 != NP)) || (DPS == 3 && (XC_Dif * YC_Dif * DIF1 * DIF2 * DIF3 != NP)))
		cout << "\aWarning: Product of Sub-Differentials is not equal to the number of blocks! \n\n", exit(-1);*/
	NP *= 2 * XC_Dif * YC_Dif;

	triangle *W; W = new triangle[NP];
	complex C;
	point C00, C01, C10, C11; double u, v, un, vn;

	cout << "Calculating Gamma Function..." << endl;

	double sx = (Max_X - Min_X) / (X_Dif * XC_Dif), sy = (Max_Y - Min_Y) / (Y_Dif * YC_Dif);
	for (int i = 0; i < X_Dif; i++) {
		//u = sx * i + Min_X, un = u + sx;
		for (int j = 0; j < Y_Dif; j++) {
			//v = sy * j + Min_Y, vn = v + sy;
			for (int id = 0; id < XC_Dif; id++) {
				u = sx * (i*XC_Dif + id) + Min_X, un = u + sx;
				for (int jd = 0; jd < YC_Dif; jd++) {
					v = sy * (j*YC_Dif + jd) + Min_Y, vn = v + sy;

					C00 = { u, v, abs(tgamma(complex(u, v))) };
					C01 = { u, vn, abs(tgamma(complex(u, vn))) };
					C10 = { un, v, abs(tgamma(complex(un, v))) };
					C11 = { un, vn, abs(tgamma(complex(un, vn))) };
					if (baddouble(C00.z) || abs(C00.z) > 100) C00.z = 0;
					if (baddouble(C01.z) || abs(C01.z) > 100) C01.z = 0;
					if (baddouble(C10.z) || abs(C10.z) > 100) C10.z = 0;
					if (baddouble(C11.z) || abs(C11.z) > 100) C11.z = 0;
					*W = triangle(C00, C01, C10), W++;
					*W = triangle(C11, C10, C01);
					C = tgamma(complex(u, v));
					if (baddouble(C.rel) || abs(C.rel) > 100) C.rel = 0;
					if (baddouble(C.ima) || abs(C.ima) > 100) C.rel = 0;
					W->reflect = (W - 1)->reflect = fromHSL(arg(C) / (2 * PI), 1, 1 - pow(0.5, abs(C)));
					W->absorb = (W - 1)->absorb = ~W->reflect;
					W++;

				}
			}

		}
	}
	W -= NP;

	cout << "Sending Data for Rendering..." << endl;

#if DPS==3
	int DF = NP / (DIF1*DIF2*DIF3);
	for (int i = 0; i < DIF1; i++) {
		G5_R.insert(new World);
		for (int j = 0; j < DIF2; j++) {
			G5_R.GObjs.back()->insert(new World);
			for (int k = 0; k < DIF3; k++) {
				for (int l = 0; l < DF; l++) {
					G5_RS.at(i).at(j).at(k)->add(W);
					W++;
				}
				G5_R.GObjs.back()->GObjs.back()->insert(G5_RS.at(i).at(j).at(k));
			}
		}
	}
#elif DPS==2
	int DF = NP / (DIF1*DIF2);
	for (int i = 0; i < DIF1; i++) {
		G5_R.insert(new World);
		for (int j = 0; j < DIF2; j++) {
			for (int l = 0; l < DF; l++) {
				G5_RS.at(i).at(j)->add(W);
				W++;
			}
			G5_R.GObjs.back()->insert(G5_RS.at(i).at(j));
		}
	}
#elif DPS==1
	int DF = NP / (DIF1);
	for (int i = 0; i < DIF1; i++) {
		G5_R.insert(new World);
		for (int l = 0; l < DF; l++) {
			G5_R.GObjs.back()->add(W);
			W++;
		}
	}
#endif
	W -= NP;


	G5.insert(&G5_R);
	img.clear(150, 100, rgb(0, 0, 255));
	img.clear(600, 600);
	G5.render(img, point(30, 60, 90), 2, 1, -0.9, 0.5, 320, 200, 40/*80, 50, 10*/);
#if DPS==3
	for (int i = 0; i < DIF1; i++) {
		for (int j = 0; j < DIF2; j++) {
			for (int k = 0; k < DIF3; k++) {
				delete G5_RS.at(i).at(j).at(k);
			}
		}
		G5_RS.at(i).clear();
	}
#elif DPS==2
	for (int i = 0; i < DIF1; i++) {
		for (int j = 0; j < DIF2; j++) {
			delete G5_RS.at(i).at(j);
		}
		G5_RS.at(i).clear();
	}
#elif DPS==1
	for (int i = 0; i < DIF1; i++) {
		delete G5_RS.at(i);
	}

#endif
	delete[] W;

#elif WORLD == 6
	World G6, G6_I;
	triangle G6_T1(point(1, -1), point(1, 1), point(0, -1)); G6_T1 += point(1.6, -0.4, 1e-4); G6_T1.reflect = rgb(255, 0, 0);
	sphere G6_B1(point(-1, -2, 0.8), 0.8); G6_B1.reflect = rgb(255, 128, 255);
	sphere G6_B2(point(0, 0, 0.5), 0.5); G6_B2.reflect = rgb(128, 255, 255);
	sphere G6_B3(point(0.8, 2, 0.6), 0.6); G6_B3.reflect = rgb(255, 255, 128);
	parallelogram G6_G(point(-2, -2, 0), point(-2, 2, 0), point(2, -2, 0), 1); G6_G.reflect = rgb(128, 192, 255); G6.add(&G6_G);
	//G6.add({ &G6_T1, &G6_B1, &G6_B2, &G6_B3 });
	srand(clock() ^ unsigned(&tm));
	vector<sphere*> vs;
#ifdef RandBall
	double bcx, bcy, bcz, r, Cr, Cg, Cb;
	for (int i = 0; i < 100; i++) {
		bcx = randnor(0, 1), bcy = randnor(0, 1); bcz = randnor(0, 0.2);
		//r = randnor(0.5, 0.1);
		r = randnor(0.2, 0.05);
		vs.push_back(new sphere(point(bcx, bcy, abs(bcz) + r), r));
		bool ns;
		do {
			ns = 0;
			for (int j = 0; j < int(vs.size()) - 2; j++) {
				if ((vs.back()->C - vs.at(j)->C).mod() < (vs.back()->r + vs.at(j)->r)) {
					bcx = randnor(0, 1), bcy = randnor(0, 1); bcz = randnor(0, 0.2);
					delete vs.back(); vs.back() = new sphere(point(bcx, bcy, abs(bcz) + r), r);
					ns = 1; break;
				}
			}
		} while (ns);
		do {
			Cr = randnor(0.7, 0.2);
		} while (Cr < 0 || Cr > 1);
		do {
			Cg = randnor(0.7, 0.2);
		} while (Cr < 0 || Cr > 1);
		do {
			Cb = randnor(0.7, 0.2);
		} while (Cr < 0 || Cr > 1);
		vs.back()->reflect = drgb(randnor(0.7, 0.2), randnor(0.7, 0.2), randnor(0.7, 0.2));
		G6_I.add(vs.back());
	}
	G6.insert(&G6_I);
	img.clear(1920, 1080); G6.render(img, point(30, 60, 90), 1.6, 0.2, -1.1, -2.1, 960, 540, 200);
#endif
	const int N = 24;
	const double R = 1.2;
	double r = R / 2 * sqrt(2 * (1 - cos(2 * PI / N)));
	for (int i = 0; i < N; i++) {
		vs.push_back(new sphere(point(R*cos(2 * PI * i / N), R*sin(2 * PI * i / N), r), r));
		vs.back()->reflect = pixel(SeaShell);
		G6_I.add(vs.back());
	}
	//G6.insert(&G6_I);
	
	circle G6_TC(point(0, 0, 0.1), 1, 0.2, 0.5); G6_TC.reflect = rgb(255, 0, 0);
	G6.add(&G6_TC);

	G6.render(img,	// canvas
		point(30, 60, 90),	// point to view (debugging)
		1.6, 0.2,	// normal of light source (debugging)
		-1.1, 0, -2.1,	// rotation, Euler's angle (debugging)
		300, 200,	// center of image
		100		// multiply
	);
#endif

	cout << "Elapsed Time: " << ((clock() - tm) / 1000.0) << "s. \n";

	img.out("D:\\Coding\\AboutMath\\SuperCalculator\\comment\\Image++.bmp");



	_CrtDumpMemoryLeaks();
	return 0;




}


// Inserted inside main(), usually for debug, DO NOT remove

#ifdef IntegerClassDebug

srand((unsigned)time(0));

integer a, b, c;
unsigned A, B, C;
unsigned long long k;
int i = 0;
while (0) {
	cin >> a >> B;
	cout << (a + B) << endl;
}
while (0) {
	a = RandomInteger_NearNode();
	b = RandomInteger_NearNode();
	a++, b++;
	c = a + b; c--;
	a = a - 1, b = b - 1;
	a--, b--;
	c -= a + b;
	c -= 3;
	if (c != 0) {
		cout << "Error! " << a << " " << b << " " << c << endl;
	}
	cout << "Test " << ++i << " Succeed.  " << /*a << " " << b << " " << c << */endl;
}

while (0) {
	//a = RandomInteger(), b = RandomInteger();
	a = RandomInteger_NearNode(), b = RandomInteger_NearNode();
	a = a * b, b = b * a, c = a + b;
	a++, ++b, ++c, c++;
	//cout << a << " " << b << " " << c << endl;
	c -= a, c += b, c -= b - a;
	a++, ++b, ++c, c++;
	c += -b, c += -a;
	if (!c.empty()) {
		cout << "Error! " << a << " " << b << " " << c << endl;
	}
	i++;
	cout << "Test " << i << " Succeed.  " << /*a << " " << b << */ c << endl;
	//cout << RandomInteger() << endl;
}

#endif

#ifdef SPLINE

#define Polyline

#ifndef PolySurface

#ifdef Polyline
matrix<dat> *M, *N;
int m, n, h;
dat a, b, c;
vector<dat> x, y, s;
vector<int> d;

while (1) {
	cout << "Enter the number of points the curve pass through: ";
	cin >> n;
	x.clear(), y.clear(), s.clear(), d.clear();
	cout << "Enter the x-coordinate and y-coordinate of points: \n";
	for (int i = 0; i < n; i++) {
		cout << "  Point " << ((char)(65 + i)) << ":  ";
		cin >> a >> b;
		x.push_back(a), y.push_back(b);
	}
	cout << endl;

	M = new Vandermonde<dat>(x, Matrix_type::_M_Row);
	//cout << "Matrix: \n" << *M << endl;
	s = M->solve(y);
	if (s.empty()) cout << "No Solution Found. Make sure each points have different x-coordinates. \n";
	else {
		/*for (int i = s.size() / 2; i >= 0; i--) {
			swap(s.at(i), s.at(s.size() - i - 1));
		}*/
		for (int i = 0; i < n; i++) {
			if (s.at(i) != 0) d.push_back(i);
		}
		for (int i = 0; i < s.size(); i++) {
			if (s.at(i) == 0) s.erase(s.begin() + i), i--;
		}
		cout << "Equation of Curve: y = ";
		for (int i = s.size() - 1; i >= 0; i--) {
			if (s.at(i) == 1);
			else if (s.at(i) == -1) cout << "-";
			else cout << s.at(i);
			if (d.at(i) == 1) cout << "x";
			else if (d.at(i) == 0) {
				if (s.at(i) == -1) cout << "1";
			}
			else cout << "x^" << d.at(i);
			if (i != 0 && s.at(i - 1) > 0) cout << "+";

		}


	}
	cout << "  \n";
	delete M;

	cout << endl << endl;
}

#endif

#ifdef SpLine
while (1) {

	cout << "Enter the number of points the curve pass through: ";
	cin >> n;
	x = new dat[n], y = new dat[n];
	cout << "Enter the x-coordinate and y-coordinate of points: \n";
	for (int i = 0; i < n; i++) {
		cout << "  Point " << ((char)(65 + i)) << ":  ";
		cin >> x[i] >> y[i];
	}
	cout << endl;


	for (int i = 1; i < n; i++) {
		if (x[i] == x[i - 1] && y[i] == y[i - 1]) {
			n--;
			for (int j = i; j < n; j++) {
				x[j] = x[j + 1];
			}
			i--;
		}
	}

	M = new matrix<dat>(n, n + 1);
	N = new matrix<dat>(n, n + 1);

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			(*M)[i][j] = (*N)[i][j] = pow(dat(i), n - j - 1);
		}
		(*M)[i][n] = x[i], (*N)[i][n] = y[i];
	}

	cout << *M << endl << *N << endl;

	elimination(*M), elimination(*N);
	//M->elimination(), N->elimination();
	cout << *M << endl << *N << endl;
	cout << "The equation of the curve is: \n";
	cout << "  x = ";
	for (int i = 0; i < n; i++) {
		if (i == 0) {
			while ((*M)[i][n] == 0) i++;
		}
		if ((*M)[i][n] != 0 && (*M)[i][n] != 1) {
			cout << (*M)[i][n];
			cout << "t^";
			if (n - i > 2) cout << (n - i - 1);
			if (n - i <= 2) cout << "\b";
			if (n - i <= 1) cout << "\b";
			if ((*M)[i + 1][n] > 0 && i + 1 != n) cout << "+";
		}
		if ((*M)[i][n] == 1) {
			cout << "t^";
			if (n - i > 2) cout << (n - i - 1);
			if (n - i <= 2) cout << "\b";
			if (n - i <= 1) cout << "\b1 ";
			if ((*M)[i + 1][n] > 0 && i + 1 != n) cout << "+";
		}
		if ((*M)[i][n] == 0) {
			if ((*M)[i + 1][n] > 0 && i + 1 != n) cout << "+";
		}
	}
	cout << "  \n  y = ";
	for (int i = 0; i < n; i++) {
		if (i == 0) {
			while ((*N)[i][n] == 0) i++;
		}
		if ((*N)[i][n] != 0 && (*N)[i][n] != 1) {
			cout << (*N)[i][n];
			cout << "t^";
			if (n - i > 2) cout << (n - i - 1);
			if (n - i <= 2) cout << "\b";
			if (n - i <= 1) cout << "\b";
			if ((*N)[i + 1][n] > 0 && i + 1 != n) cout << "+";
		}
		if ((*N)[i][n] == 1) {
			cout << "t^";
			if (n - i > 2) cout << (n - i - 1);
			if (n - i <= 2) cout << "\b";
			if (n - i <= 1) cout << "\b1 ";
			if ((*N)[i + 1][n] > 0 && i + 1 != n) cout << "+";
		}
		if ((*N)[i][n] == 0) {
			if ((*N)[i + 1][n] > 0 && i + 1 != n) cout << "+";
		}
	}
	cout << "  \n  t in [0," << (n - 1) << "]\n";

	cout << endl << endl;

	delete M, N, x, y;
}
#endif

#else

matrix<dat> *M;
int m, n, h, k; bool s;
double *x, *y, *z, c;

while (1) {
	cout << "Enter the degree of the polynomial surface: ";
	cin >> m;
	n = (m + 2)*(m + 1) / 2;

	x = new double[n], y = new double[n], z = new double[n];

	cout << "Enter the x, y and z-coordinate of " << n << " points: \n";
	for (int i = 0; i < n; i++) {
		cout << "  Point #" << (i + 1) << ":\t\b\b\b";
		cin >> x[i] >> y[i] >> z[i];
	}
	cout << endl;


	h = 1;
	while (h < n / 3) {
		h = 3 * h + 1;
	}
	while (h >= 1) {
		for (int i = h; i < n; i++) {
			for (int j = i; j >= h && y[j] < y[j - h]; j -= h) {
				swap(x[j], x[j - h]), swap(y[j], y[j - h]), swap(z[j], z[j - h]);
			}
		}
		h /= 3;
	}
	h = 1;
	while (h < n / 3) {
		h = 3 * h + 1;
	}
	while (h >= 1) {
		for (int i = h; i < n; i++) {
			for (int j = i; j >= h && x[j] < x[j - h]; j -= h) {
				swap(x[j], x[j - h]), swap(y[j], y[j - h]), swap(z[j], z[j - h]);
			}
		}
		h /= 3;
	}
	// Shell Sorting Alogorithm

	s = 1;
	for (int i = 1; i < n; i++) {
		if (x[i] == x[i - 1] && y[i] == y[i - 1]) {
			cout << "x and y-coordinate cannot have the same value. \n";
			delete x, y, z;
			s = 0;
			break;
		}
	}

	if (s) {
		M = new matrix<dat>(n, n + 1);
		for (int l = 0; l < n; l++) {
			k = n;
			for (int i = 0; i <= m; i++) {
				for (int j = 0; j <= i; j++) {
					k--;
					(*M)[l][k] = pow(x[l], j)*pow(y[l], i - j);
					//cout << "M[" << l << "][" << k << "] = (" << x[l] << ")^" << j << "*(" << y[l] << ")^" << (i-j) << " = " << ((*M)[l][k]) << endl;
				}
			}
			(*M)[l][n] = z[l];
		}
		//cout << "Matrix: \n" << *M << endl;
		do {
			//ofstream coutf("D:\\Coding\\AboutMath\\SuperCalculator\\Matrix.txt");
			//coutf << *M << endl;
			//coutf.close();

			//M->elimination();
			elimination(*M);
			//cout << "Eliminationed Matrix: \n" << *M << endl;
			s = 0;
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < n; j++) {
					if (i != j && (*M)[i][j] != 0) s = 1;
					if (i == j && (*M)[i][j] == 0) s = 1;
				}
			}
			if (s == 1) {
				cout << "The final solution may not be correct due to the restriction of float point bit depth. \n";
				cout << *M << endl;
				//pause;
			}
		} while (s);
		cout << "Equation of Surface: z = ";
		k = 0;
		for (int i = m; i >= 0; i--) {
			for (int j = i; j >= 0; j--) {
				if ((*M)[k][n] == 0);
				else if (j != 0 && j != i) {
					cout << (*M)[k][n];
					if (j != 1) cout << "*x^" << j;
					else cout << "*x";
					if (i - j != 1) cout << "*y^" << (i - j);
					else cout << "*y";
				}
				else if (j == 0 && j != i) {
					cout << (*M)[k][n];
					if (i != 1) cout << "*y^" << i;
					else cout << "*y";
				}
				else if (j != 0 && j == i) {
					cout << (*M)[k][n];
					if (j != 1) cout << "*x^" << j;
					else cout << "*x";
				}
				else if (j == 0 && i == 0) cout << (*M)[k][n];
				k++;
				if (k != n && (*M)[k][n] > 0) cout << "+";
			}
		}
		cout << endl << endl;
		delete M, x, y;
	}

	cout << endl << endl;
}

#endif

#endif

#ifdef MatrixTest

matrix<dat> M1({
		{2, 1, -1, 8},
		{-3, -1, 2, -11},
		{-2, 1, 2, -3},
	});
matrix<dat> M2({
		{2, 0, -1, },
		{1, 0, 0, -1, },
		{3, 0, 0, -2, -1, },
		{0, 1, 0, 0, -2, },
		{0, 1, -1, 0, 0, },
		{1, 0, 0, 0, 0, 1, },
	});
matrix<dat> M3(3, 4, {
		1, 1, 1, 1,
		1, 2, -5, 2,
		2, 3, -4, 5,
	});
matrix<dat> M4(4, 6, {
		1, -1, -1, 0, 3, -1,
		2, -2, -1, 2, 4, -2,
		3, -3, -1, 4, 5, -3,
		1, -1, 1, 1, 8, 2,
	});
matrix<dat> M5(2, {
		1, 2,
		2, 1
	});
matrix<dat> M6(3, {
		1, 2, 3,
		4, 5, 7,
		6, 8, 9
	});
matrix<dat> M7(4, {
		1, 4, -1, 4,
		2, 1, 4, 3,
		4, 2, 3, 11,
		3, 0, 9, 2,
	});
matrix<dat> M8(3, {
		1, 2, 3,
		4, 5, 6,
		7, 8, 9
	});

matrix<dat> O3(3);
matrix<dat> I3(3, Matrix_type::_M_Identity);
matrix<dat> O3_4(3, 4);
Vandermonde<dat> V3({ 0, 1, 2 }, Matrix_type::M_type(2));

srand((unsigned)clock() ^ (unsigned)&M1);

#define use V3

matrix<dat> MA, MB, MC, MD, ME, MF;
dat a, b, c, k;
unsigned m, n, u, v;
#ifdef A_B_times_C_times__equal__A_B_C_times_times
for (int i = 0; i < 1000; i++) {
	MA = randMatrix();
	while (MA.width() != MB.height()) MB = randMatrix();
	while (MC.height() != MB.width()) MC = randMatrix();
	ME = (MA*MB)*MC;
	MF = MA * (MB*MC);
	if (ME != MF) {
		cout << "MA: \n" << MA << endl;
		cout << "MB: \n" << MB << endl;
		cout << "MC: \n" << MC << endl;
		cout << "(MA*MB)*MC: \n" << ME << endl;
		cout << "MA*(MB*MC): \n" << MF << endl << endl;
		cout << endl << "===============================================" << endl << endl;
		pause;
	}
}
#endif // (AB)*C = A*(BC)
#ifdef c_A_B_times_times__equal__A_c_times_B_times
for (int i = 0; i < 1000; i++) {
	MA = randMatrix();
	while (MA.width() != MB.height()) MB = randMatrix();
	c = randFrac();
	ME = (c*MA) * MB;
	MF = c * (MA*MB);
	if (ME != MF) {
		cout << "MA: \n" << MA << endl;
		cout << "MB: \n" << MB << endl;
		cout << "c: " << c << endl;
		cout << "(c*MA)*MB: \n" << ME << endl;
		cout << "c*(MA*MB): \n" << MF << endl << endl;
		cout << endl << "===============================================" << endl << endl;
		pause;
	}
}
#endif // c*(AB) = (A*c)*B
#ifdef A_B_plus_C_times__equal__A_C_times_B_C_times_plus
for (int i = 0; i < 1000; i++) {
	MA = randMatrix();
	while (MB.width() != MA.width() || MB.height() != MA.height()) MB = randMatrix();
	while (MA.width() != MC.height()) MC = randMatrix();
	//ME = (MA + MB)*MC;
	ME = MA + MB; ME *= MC;
	MF = MA * MC + MB * MC;
	if (ME != MF) {
		cout << "MA: \n" << MA << endl;
		cout << "MB: \n" << MB << endl;
		cout << "MC: \n" << MC << endl;
		cout << "(MA+MB)*MC: \n" << ME << endl;
		cout << "MA*MC+MB*MC: \n" << MF << endl << endl;
		cout << endl << "===============================================" << endl << endl;
		pause;
	}
}
#endif // (A+B)*C = AC+BC
#ifdef A_B_C_plus_times__equal__A_B_times_A_C_times_plus
for (int i = 0; i < 1000; i++) {
	MA = randMatrix();
	while (MA.width() != MB.height()) MB = randMatrix();
	while (MB.width() != MC.width() || MB.height() != MC.height()) MC = randMatrix();
	ME = MA * (MB + MC);
	MF = MA * MB + MA * MC;
	if (ME != MF) {
		cout << "MA: \n" << MA << endl;
		cout << "MB: \n" << MB << endl;
		cout << "MC: \n" << MC << endl;
		cout << "(MA+MB)*MC: \n" << ME << endl;
		cout << "MA*MC+MB*MC: \n" << MF << endl << endl;
		cout << endl << "===============================================" << endl << endl;
		pause;
	}
}
#endif // A*(B+C) = AB+AC
#ifdef A_transpose_transpose__equal__A
for (int i = 0; i < 1000; i++) {
	MA = randMatrix();
	MB = MA.transpose().transpose();
	if (MB != MA) {
		cout << "A: \n" << MA << endl;
		cout << "(A^T)^T: \n" << MB << endl << endl;
		cout << endl << "===============================================" << endl << endl;
		pause;
	}
}
#endif // A^T^T = A
#ifdef A_B_plus_transpose__equal__A_transpose_B_transpose_plus
for (int i = 0; i < 1000; i++) {
	MA = randMatrix();
	while (MA.width() != MB.width() || MA.height() != MB.height()) MB = randMatrix();
	ME = (MA + MB).transpose();
	MF = MA.transpose() + MB.transpose();
	if (ME != MF) {
		cout << "A: \n" << MA << endl;
		cout << "B: \n" << MB << endl;
		cout << "(A+B)^T: \n" << ME << endl;
		cout << "A^T+B^T: \n" << MF << endl << endl;
		cout << endl << "===============================================" << endl << endl;
		pause;
	}
}
#endif // (A+B)^T = A^T+B^T
#ifdef c_A_times_transpose__equal__c_A_transpose_times
for (int i = 0; i < 1000; i++) {
	MA = randMatrix();
	c = randFrac();
	ME = (c*MA).transpose();
	MF = MA.transpose()*c;
	if (ME != MF) {
		cout << "A: \n" << MA << endl;
		cout << "c: " << c << endl;
		cout << "(c*A)^T: \n" << ME << endl;
		cout << "c*A^T: \n" << MF << endl << endl;
		cout << endl << "===============================================" << endl << endl;
		pause;
	}
}
#endif // (c*A)^T = c*(A^T)
#ifdef A_B_times_transpose__equal__B_transpose_A_transpose_times
for (int i = 0; i < 1000; i++) {
	MA = randMatrix();
	while (MA.width() != MB.height()) MB = randMatrix();
	ME = (MA*MB).transpose();
	MF = MB.transpose()*MA.transpose();
	if (ME != MF) {
		cout << "A: \n" << MA << endl;
		cout << "B: \n" << MB << endl;
		cout << "(A*B)^T: \n" << ME << endl;
		cout << "B^T*A^T: \n" << MF << endl << endl;
		cout << endl << "===============================================" << endl << endl;
		pause;
	}
}
#endif // (A*B)^T = B^T*A^T
#ifdef A_B_times_C_times_transpose__equal__C_transpose_B_transpose_times_A_transpose_times
for (int i = 0; i < 1000; i++) {
	MA = randMatrix();
	while (MA.width() != MB.height()) MB = randMatrix();
	while (MB.width() != MC.height()) MC = randMatrix();
	ME = (MA*MB*MC).transpose();
	MF = MC.transpose()*MB.transpose()*MA.transpose();
	if (ME != MF) {
		cout << "A: \n" << MA << endl;
		cout << "B: \n" << MB << endl;
		cout << "C: \n" << MC << endl;
		cout << "(A*B*C)^T: \n" << ME << endl;
		cout << "C^T*B^T*A^T: \n" << MF << endl << endl;
		cout << endl << "===============================================" << endl << endl;
		pause;
	}
}
#endif // (A*B*C)^T = C^T*B^T*A^T
#ifdef A_inverse_inverse__equal__A
/*while (1) {
	while (!MA.square()) MA = randMatrix();
	MB = MA;
	MA.invertible();
	cout << MB << endl << MA << endl;
	if (MB != MA) cout << "\a";
	*(unsigned*)&MA = 0;
	pause;
}*/
for (int i = 0; i < 1000; i++) {
	while (!MA.invertible()) MA = randMatrix();
	MB = MA.inverse();
	MC = MB.inverse();
	if (MC != MA) {
		cout << "Test " << (i + 1) << endl;
		cout << "A: \n" << MA << endl;
		cout << "A^-1: \n" << MB << endl;
		cout << "A^-1^-1: \n" << MC << endl << endl;
		cout << endl << "===============================================" << endl << endl;
		pause;
	}
	MA.clear();
}
#endif // A^-1^-1 = A
#ifdef c_A_times_inverse__equal__A_inverse_c_divide
for (int i = 0; i < 1000; i++) {
	while (!MA.invertible()) MA = randMatrix();
	c = 0;
	while (c == 0) c = randFrac();
	MB = (MA*c).inverse();
	MC = MA.inverse() / c;
	if (MB != MC) {
		cout << "Test " << (i + 1) << endl;
		cout << "A: \n" << MA << endl;
		cout << "c: \n" << c << endl;
		cout << "(c*A)^-1: \n" << MB << endl;
		cout << "(A^-1)/c: \n" << MC << endl << endl;
		cout << endl << "===============================================" << endl << endl;
		pause;
	}
	MA.clear();
}
#endif // (c*A)^-1 = (A^-1)/c
#ifdef A_B_times_inverse__equal__B_inverse_A_inverse_times
for (int i = 0; i < 1000; i++) {
	MA = randMatrix(), MB = randMatrix(); c = det(MA*MB);
	while ((!MA.invertible() || MA.height() > 6 || (c.denominator() + c.numerator()) > 256)
		|| (!(MA.width() == MB.width() && MB.invertible()))) {
		MA = randMatrix(), MB = randMatrix(), c = det(MA*MB);
	}
	ME = (MA*MB).inverse();
	MF = MB.inverse()*MA.inverse();
	if (ME != MF) {
		cout << "Test " << (i + 1) << endl;
		cout << "A: \n" << MA << endl;
		cout << "B: \n" << MB << endl;
		cout << "(A*B)^-1: \n" << ME << endl;
		cout << "(B^-1)*(A^-1): \n" << MF << endl << endl;
		cout << endl << "===============================================" << endl << endl;
		pause;
	}
}
#endif // (A*B)^-1 = (B^-1)*(A^-1)
#ifdef A_inverse_transpose__equal__A_transpose_inverse
for (int i = 0; i < 1000; i++) {
	while (!MA.invertible()) MA = randMatrix();
	MB = MA.inverse().transpose();
	MC = MA.transpose().inverse();
	if (MB != MC) {
		cout << "Test " << (i + 1) << endl;
		cout << "A: \n" << MA << endl;
		cout << "A^-1^T: \n" << MB << endl;
		cout << "A^T^-1: \n" << MC << endl << endl;
		cout << endl << "===============================================" << endl << endl;
		pause;
	}
	MA.clear();
}
#endif // A^-1^T = A^T^-1
#ifdef A_k1_power_A_k2_power_times__equal__A_k1_k2_plus_power
int e1, e2;
for (int i = 0; i < 1000; i++) {
	e1 = round(randnor(1.5, 2)), e2 = round(randnor(1.5, 2));
	while (abs(e1) + abs(e2) > 6) {
		e1 = round(randnor(1.5, 2)), e2 = round(randnor(1.5, 2));
	}
	MA = randMatrix(); c = det(MA) ^ (abs(e1) + abs(e2));
	while (!MA.square() || MA.height() > 6 || ((c.denominator() + c.numerator()) > 256) || ((e1 < 0 || e2 < 0) && !MA.invertible())) {
		MA = randMatrix(), c = det(MA) ^ (abs(e1) + abs(e2));
	}
	ME = (MA^e1)*(MA^e2);
	MF = MA ^ (e1 + e2);
	if (ME != MF) {
		cout << "Test " << (i + 1) << endl;
		cout << "A: \n" << MA << endl;
		cout << "k1: " << e1 << endl << endl;
		cout << "k2: " << e2 << endl << endl;
		cout << "(A^k1)*(A^k2): \n" << ME << endl;
		cout << "A^(k1+k2): \n" << MF << endl << endl;
		cout << endl << "===============================================" << endl << endl;
		pause;
	}
}
#endif // (A^k1)*(A^k2) = A^(k1+k2)
#ifdef A_k1_power_k2_power__equal__A_k1_k2_times_power
int e1, e2;
for (int i = 0; i < 1000; i++) {
	e1 = round(randnor(1.5, 2)), e2 = round(randnor(1.5, 2));
	while (abs(e1) * abs(e2) + (e1 < 0) + (e2 < 0) > 5) {
		e1 = round(randnor(1.5, 2)), e2 = round(randnor(1.5, 2));
	}
	MA = randMatrix(); c = det(MA) ^ (abs(e1) * abs(e2) + (e1 < 0) + (e2 < 0));
	while (!MA.square() || MA.height() > 5 || ((c.denominator() + c.numerator()) > 256) || ((e1 < 0 || e2 < 0) && !MA.invertible())) {
		MA = randMatrix(), c = det(MA) ^ (abs(e1) * abs(e2) + (e1 < 0) + (e2 < 0));
	}
	ME = (MA^e1) ^ e2;
	MF = MA ^ (e1 * e2);
	if (ME != MF) {
		cout << "Test " << (i + 1) << endl;
		cout << "A: \n" << MA << endl;
		cout << "k1: " << e1 << endl << endl;
		cout << "k2: " << e2 << endl << endl;
		cout << "(A^k1)^k2: \n" << ME << endl;
		cout << "A^(k1*k2): \n" << MF << endl << endl;
		cout << endl << "===============================================" << endl << endl;
		pause;
	}
}
#endif // (A^k1)^k2 = A^(k1*k2)
#ifdef A_det__equal__A_determinant
for (int i = 0; i < 10000; i++) {
	while (!MA.square()) MA = randMatrix();
	a = det(MA);
	b = determinant(MA);
	if (a != b) {
		cout << "A: \n" << MA << endl;
		cout << "det(A): " << a << endl;
		cout << "|A|: " << b << endl << endl;
		cout << endl << "===============================================" << endl << endl;
		pause;
	}
	*((unsigned*)&MA) = 0;
}
#endif // det(A) = |A|
#ifdef A_det__equal__A_transpose_det
for (int i = 0; i < 10000; i++) {
	while (!MA.square()) MA = randMatrix();
	MB = MA.transpose();
	a = det(MA); b = det(MB);
	if (a != b) {
		cout << "A: \n" << MA << endl;
		cout << "A^T: \n" << MB << endl;
		cout << "det(A): " << a << endl;
		cout << "det(A^T): " << b << endl << endl;
		cout << endl << "===============================================" << endl << endl;
		pause;
	}
	*((unsigned*)&MA) = 0;
}
#endif // det(A) = det(A^T)
#ifdef A_det_reciprocate__equal__A_inverse_det
for (int i = 0; i < 1000; i++) {
	while (!MA.invertible()) MA = randMatrix();
	MB = MA.inverse();
	a = det(MA); b = det(MB);
	if (a * b != 1) {
		cout << "Test " << (i + 1) << endl;
		cout << "A: \n" << MA << endl;
		cout << "A^-1: \n" << MB << endl;
		cout << "det(A): " << a << endl;
		cout << "det(A^-1): " << b << endl << endl;
		cout << endl << "===============================================" << endl << endl;
		pause;
	}
	*((unsigned*)&MA) = 0;
}
#endif // det(A)^-1 = det(A^-1)
#ifdef A_B_times_det__equal__A_det_B_det_times
for (int i = 0; i < 1000; i++) {
	MA = randMatrix(), MB = randMatrix(); c = det(MA*MB);
	while (!MA.square() || !MB.square() || MA.height() > 6 || MA.width() != MB.width() || (c.denominator() + c.numerator()) > 256) {
		MA = randMatrix(), MB = randMatrix(), c = det(MA*MB);
	}
	a = det(MA*MB);
	b = det(MA)*det(MB);
	if (a != b) {
		cout << "Test " << (i + 1) << endl;
		cout << "A: \n" << MA << endl;
		cout << "B: \n" << MB << endl;
		cout << "det(A): " << det(MA) << endl;
		cout << "det(B): " << det(MB) << endl;
		cout << "det(A*B): " << a << endl;
		cout << "det(A)*det(B): " << b << endl << endl;
		cout << endl << "===============================================" << endl << endl;
		pause;
	}
}
#endif // det(A*B) = det(A)*det(B)
#ifdef A_interchange_det__equal__A_det_negative
int l1, l2;
for (int i = 0; i < 10000; i++) {
	do { MA = randMatrix(); } while (!MA.sqr() || MA.height() == 1);
	MB = MA;
	l1 = rand() % MA.height() + 1;
	do {
		l2 = rand() % MA.height() + 1;
	} while (l2 == l1);
	if (rand() & 1) MB.interchange(l1, l2);
	else MB.col_interchange(l1, l2);
	a = det(MA), b = det(MB);
	if (a + b != 0) {
		cout << "Test " << (i + 1) << endl;
		cout << "A: \n" << MA << endl;
		cout << "Interchange: \n" << MB << endl;
		cout << "det(A): " << a << endl;
		cout << "det(B): " << b << endl << endl;
		cout << endl << "===============================================" << endl << endl;
		pause;
	}
}
#endif // interchange rows/cols, determinant change sign
#ifdef A_multiply_det__equal__A_det_multiply
int l1;
for (int i = 0; i < 10000; i++) {
	do { MA = randMatrix(); } while (!MA.sqr());
	c = randFrac();
	MB = MA;
	l1 = rand() % MA.height() + 1;
	if (rand() & 1) MB.multiply(l1, c);
	else MB.col_multiply(l1, c);
	a = det(MA), b = det(MB);
	if (b != a * c) {
		cout << "Test " << (i + 1) << endl;
		cout << "A: \n" << MA << endl;
		cout << "Multiply: \n" << MB << endl;
		cout << "det(A): " << a << endl;
		cout << "det(B): " << b << endl << endl;
		cout << endl << "===============================================" << endl << endl;
		pause;
	}
}
#endif // multiply rows/cols, determinant multiply
#ifdef A_multiply_det__equal__A_det_multiply
int l1, l2;
for (int i = 0; i < 10000; i++) {
	do { MA = randMatrix(); } while (!MA.sqr() || MA.size() == 1);
	c = randFrac();
	MB = MA;
	l1 = rand() % MA.height() + 1;
	do {
		l2 = rand() % MA.height() + 1;
	} while (l2 == l1);
	if (rand() & 1) MB.multiply(l1, c, l2);
	else MB.col_multiply(l1, c, l2);
	a = det(MA), b = det(MB);
	if (b != a) {
		cout << "Test " << (i + 1) << endl;
		cout << "A: \n" << MA << endl;
		cout << "Multiply: \n" << MB << endl;
		cout << "det(A): " << a << endl;
		cout << "det(B): " << b << endl << endl;
		cout << endl << "===============================================" << endl << endl;
		pause;
	}
}
#endif // multiply rows/cols, determinant multiply
#ifdef A_transpose_adjugate__equal__A_adjugate_transpose
for (int i = 0; i < 1000; i++) {
	do { MA = randMatrix(); } while (!MA.sqr());
	ME = MA.transpose().adjugate();
	MF = MA.adjugate().transpose();
	if (ME != MF) {
		cout << "Test " << (i + 1) << endl;
		cout << "A: \n" << MA << endl;
		cout << "(A^T)*: \n" << ME << endl;
		cout << "(A*)^T: \n" << MF << endl << endl;
		cout << endl << "===============================================" << endl << endl;
		pause;
	}
}
#endif // (A^T)* = (A*)^T
#ifdef A_inverse_adjugate__equal__A_adjugate_inverse
for (int i = 0; i < 1000; i++) {
	do { MA = randMatrix(); c = det(MA); } while (!MA.invertible() || c.numerator() + c.denominator() >= 256);
	ME = MA.inverse().adjugate();
	MF = MA.adjugate().inverse();
	if (ME != MF) {
		cout << "Test " << (i + 1) << endl;
		cout << "A: \n" << MA << endl;
		cout << "(A^-1)*: \n" << ME << endl;
		cout << "(A*)^-1: \n" << MF << endl << endl;
		cout << endl << "===============================================" << endl << endl;
		pause;
	}
}
#endif // (A^-1)* = (A*)^-1
#ifdef A_B_times_adjugate__equal__B_adjugate_A_adjugate_times
for (int i = 0; i < 1000; i++) {
	do {
		MA = randMatrix(); MB = randMatrix(); c = det(MA*MB) ^ MA.height();
	} while (!MA.square() || !MB.square() || MA.height() != MB.height() || c.denominator() + c.numerator() > 256);
	//MA = { {0, fraction(3,2), fraction(3,2)}, {1, 3, 2}, {2, 0, 2} };
	//MB = { {0, 0, 5}, {0, 0, 4}, {1, 1, fraction(3,2)} };
	ME = adj(MA*MB);
	MF = adj(MB)*adj(MA);
	if (ME != MF) {
		cout << "Test " << (i + 1) << endl;
		cout << "A: \n" << MA << endl;
		cout << "B: \n" << MB << endl;
		cout << "AB: \n" << (MA*MB) << endl;
		cout << "(AB)*: \n" << ME << endl;
		cout << "B*: \n" << adj(MB) << endl;
		cout << "A*: \n" << adj(MA) << endl;
		cout << "B*A*: \n" << MF << endl << endl;
		cout << endl << "===============================================" << endl << endl;
		pause;
	}
}
#endif // (A×B)* = (B*)×(A*)
#ifdef A_adjugate__equal__A_det_A_inverse_times
for (int i = 0; i < 1000; i++) {
	do {
		MA = randMatrix();
	} while (!MA.invertible());
	MD = adj(MA) / det(MA);
	ME = MD * MA;
	MF = MA * MD;
	if (!(ME.identity() || MF.identity())) {
		cout << "Test " << (i + 1) << endl;
		cout << "A: \n" << MA << endl;
		cout << "A*/|A|(X): \n" << MD << endl;
		cout << "X*A: \n" << ME << endl << endl;
		cout << "A*X: \n" << MF << endl << endl;
		cout << endl << "===============================================" << endl << endl;
		pause;
	}
}
#endif // A* = det(A)*A^-1
#ifdef A_transpose_rank__equal__A_rank
for (int i = 0; i < 10000; i++) {
	MA = randMatrix();
	MB = MA.transpose();
	m = R(MA); n = R(MB);
	if (m != n) {
		cout << "A: \n" << MA << endl;
		cout << "A^T: \n" << MB << endl;
		cout << "R(A): " << m << endl;
		cout << "R(A^T): " << n << endl << endl;
		cout << endl << "===============================================" << endl << endl;
		pause;
	}
	*((unsigned*)&MA) = 0;
}
#endif // R(A) = R(A^T)
#ifdef A_fullrank__equivalent__A_invertible
for (int i = 0; i < 1000; i++) {
	do { MA = randMatrix(); } while (!MA.sqr());
	if (MA.invertible() ^ (R(MA) == MA.height())) {
		cout << "A: \n" << MA << endl;
		cout << "A^-1: \n" << MA.inverse() << endl;
		cout << "R(A): " << R(MA) << endl << endl;
		cout << endl << "===============================================" << endl << endl;
		pause;
	}
}
#endif // R(A)=n <==> A invertible
#ifdef A_transpose_A_times_rank__equal__A_rank
for (int i = 0; i < 1000; i++) {
	MA = randMatrix();
	//while (MA.size() > 20) MA = randMatrix();
	/*MA = {
		{6,2,fraction(1,2),3,2},
		{3,0,0,3,fraction(1,2)},
		{fraction(1,2),3,0,fraction(1,2),2},
		{1,0,7,0,0}
	};*/
	MB = MA.transpose()*MA;
	m = R(MA); n = R(MB);
	if (m != n) {
		cout << "Test " << (i + 1) << endl;
		cout << "A: \n" << MA << endl;
		cout << "A^T: \n" << MA.transpose() << endl;
		cout << "A^T*A: \n" << MB << endl;
		cout << "R(A): " << m << endl;
		cout << "R(A^T*A): " << n << endl << endl;
		cout << endl << "===============================================" << endl << endl;
		pause;
	}
}
#endif // R(A) = R(A^T*A)
#ifdef k_A_times_rank__equal__A_rank
for (int i = 0; i < 1000; i++) {
	MA = randMatrix();
	do { k = randFrac(); } while (k == 0);
	MB = k * MA;
	m = R(MA); n = R(MB);
	if (m != n) {
		cout << "Test " << (i + 1) << endl;
		cout << "A: \n" << MA << endl;
		cout << "k: " << k << endl;
		cout << "k*A: \n" << MB << endl;
		cout << "R(A): " << m << endl;
		cout << "R(k*A): " << n << endl << endl;
		cout << endl << "===============================================" << endl << endl;
		pause;
	}
}
#endif // R(k*A) = R(A)
#ifndef rank_of_adjugate_matrix
for (int i = 0; i < 1000; i++) {
	do { MA = randMatrix(); } while (!MA.square() || (MA.invertible() && rand() % 4) || MA.size() <= 1);
	unsigned h = MA.height();
	MB = adj(MA);
	m = R(MA); n = R(MB);
	if ((m == h && n != h) || (m == h - 1 && n != 1) || (m <= h - 2 && n != 0)) {
		cout << "A: \n" << MA << endl;
		cout << "A*: \n" << MB << endl;
		cout << "R(A): " << m << endl;
		cout << "R(A*): " << n << endl << endl;
		cout << endl << "===============================================" << endl << endl;
		pause;
	}
}
#endif
return 0;
#endif

#ifdef GIFTest

VLCLZWEnco(3); return 0;
//LZWEnco(); LZWDeco(); return 0;

//bitmap BMPIn("D:\\Coding\\AboutMath\\SuperCalculator\\Uncompressed_GIF.bmp");
bitmap BMPIn("D:\\Coding\\AboutMath\\SuperCalculator\\Image.bmp");
GIF GIFOut(BMPIn.width(), BMPIn.height());
GIFOut = BMPIn;
//GIFOut.out("D:\\Coding\\AboutMath\\SuperCalculator\\Uncompressed_GIF_My.gif");
GIFOut.out("D:\\Coding\\AboutMath\\SuperCalculator\\Image.gif");
return 0;

#endif

