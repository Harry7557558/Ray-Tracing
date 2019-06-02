#include "World.h"

using namespace std;

#define ADD_AXIS(W, r) { \
	cylinder xAxis(point(0, 0, 0), point(ERR_UPSILON, 0, 0), r); xAxis.setcolor(Red);    \
	cylinder yAxis(point(0, 0, 0), point(0, ERR_UPSILON, 0), r); yAxis.setcolor(Green);	 \
	cylinder zAxis(point(0, 0, 0), point(0, 0, ERR_UPSILON), r); zAxis.setcolor(Blue);	 \
	W.add({ &xAxis, &yAxis, &zAxis }); }


void sizetest() {
	cout << "point          " << sizeof(point) << endl;
	cout << "ray            " << sizeof(ray) << endl;
	cout << "intersect      " << sizeof(intersect) << endl;
	cout << "rgblight       " << sizeof(rgblight) << endl;
	cout << endl;

	cout << "object         " << sizeof(object) << endl;
	cout << "object2D       " << sizeof(objectSF) << endl;
	cout << "plane          " << sizeof(plane) << endl;
	cout << "triangle       " << sizeof(triangle) << endl;
	cout << "parallelogram  " << sizeof(parallelogram) << endl;
	cout << "circle         " << sizeof(circle) << endl;
	cout << "cylinder       " << sizeof(cylinder) << endl;
	cout << "sphere         " << sizeof(sphere) << endl;
	cout << "object3D       " << sizeof(object3D) << endl;
	cout << endl;

	cout << "World          " << sizeof(World) << endl;
	cout << endl;
}

void MonteCarlo_Test() {
	const unsigned MaxSamp = 1024;
	bitmap res(MaxSamp, 100);
	double sum = 0;
	for (int j = 0; j < 100; j++) {
		sum = 0;
		for (int i = 0; i < MaxSamp; i++) {
			sum += pick_random(1);
			res[j][i] = rgb(rgblight(sum / (i + 1)));
		}
	}
	for (int i = 0; i < MaxSamp; i++) {
		for (int j = 0; j < 100; j++) {
			swap(res[rand() % 100][i], res[rand() % 100][i]);
		}
		int s = sqrt(i + 1);
		if (s * s == i + 1) {
			res[0][i] = res[99][i] = rgb(255, 0, 0);
			if (s % 4 == 0) res[1][i] = res[98][i] = rgb(255, 255, 0);
		}
	}
	res.out("IMAGE\\MonteCarlo_Test.bmp");
}


// light-blue sphere and plane
void Render_Test01() {
	World W;
	sphere S(point(0, 0, 1), 1);
	S.setcolor(White);
	parallelogram P(point(-2, -2, 0), point(4, 0, 0), point(0, 4, 0));
	P.setcolor(White);
	plane Pi(point(-2.4, 0, 0), point(1, 0, 0));
	Pi.setcolor(LightBlue);
	W.add({ &S, &P, &Pi });
	W.setGlobalLightSource(1, 1);
	bitmap canvas(600, 400);
	W.render(canvas, point(40, -40, 40), point(0, 0, 0.4), 0, 0.006);
	canvas.out("IMAGE\\RT1.bmp");
	W.render(canvas, point(50, 20, 30), point(0, 0, 0.4), 0, 0.004);
	canvas.out("IMAGE\\RT2.bmp");
}

void Render_Test02() {
	World W;
	sphere S(point(0, 0, 1), 1);
	S.setcolor(White);
	parallelogram P(point(-2, -2, 0), point(4, 0, 0), point(0, 4, 0));
	P.setcolor(White);
	plane Pi(point(-2.4, 0, 0), point(1, 0, 0));
	Pi.setcolor(LightBlue);
	W.add({ &S, &P, &Pi });
	W.setGlobalLightSource(1, 1);
	bitmap canvas(600, 400);
	W.render(canvas, point(50, 20, 30), point(0, 0, 0.5), 0, 0.004);
	canvas.out("IMAGE\\RT.bmp");
}

// triangle-based pyramid
void Render_GTest01() {
	// This function has serious memory leaks, only for debug purpose
	World G6, G6_I, G6_IS, G6_IST;
	parallelogram G6_G(point(-200, -200, 0), point(-200, 200, 0), point(200, -200, 0), 1); G6_G.setcolor(LightBlue);
	double r = 1;
	double Cr = r / 2, Ch = 6, Cho2 = Ch / 2;
	int stair = 3;
	int substair = 2;	// A bug occurs when this is greater than 2
	int allstair = stair * substair;
	const double rt2 = sqrt(2);
	const double rt3 = sqrt(3);
	const double rt3o2 = sqrt(3) / 2;
	const double rt3o3 = sqrt(3) / 3;
	const double rt3o6 = sqrt(3) / 6;
	const double rt6o3 = sqrt(6) / 3;
	const double rt6t2o3 = sqrt(6) * 2 / 3;
	for (int i = 0; i < stair; i++) {
		for (int j = 0; j < stair - i; j++) {
			for (int k = 0; k < stair - i - j; k++) {
				for (int is = i * substair; (is < (i + 1)*substair && is < allstair) || (i + 1 == stair && is == allstair); is++) {
					for (int js = j * substair; (js < (j + 1)*substair && is + js < allstair) || (i + j + 1 == stair && is + js == allstair); js++) {
						for (int ks = k * substair; (ks < (k + 1)*substair && is + js + ks < allstair) || (i + j + k + 1 == stair && is + js + ks == allstair); ks++) {
							G6_IST.add(new sphere(point(Ch*is + Cho2 * js + Cho2 * ks, rt3o2*Ch*js + rt3o6 * Ch*ks, rt6o3*Ch*ks + r), r));
							const_cast<object*>(G6_IST.Objs.back())->setcolor(Silver);
							if (is + js + ks < allstair) {
								G6_IST.add(new cylinder(point(Ch*is + Cho2 * js + Cho2 * ks, rt3o2*Ch*js + rt3o6 * Ch*ks, rt6o3*Ch*ks + r),
									point(Ch*(is + 1) + Cho2 * js + Cho2 * ks, rt3o2*Ch*js + rt3o6 * Ch*ks, rt6o3*Ch*ks + r), Cr));
								const_cast<object*>(G6_IST.Objs.back())->setcolor(LightSkyBlue);
								G6_IST.add(new cylinder(point(Ch*is + Cho2 * js + Cho2 * ks, rt3o2*Ch*js + rt3o6 * Ch*ks, rt6o3*Ch*ks + r),
									point(Ch*is + Cho2 * (js + 1) + Cho2 * ks, rt3o2*Ch*(js + 1) + rt3o6 * Ch*ks, rt6o3*Ch*ks + r), Cr));
								const_cast<object*>(G6_IST.Objs.back())->setcolor(LightSkyBlue);
								G6_IST.add(new cylinder(point(Ch*is + Cho2 * js + Cho2 * ks, rt3o2*Ch*js + rt3o6 * Ch*ks, rt6o3*Ch*ks + r),
									point(Ch*is + Cho2 * js + Cho2 * (ks + 1), rt3o2*Ch*js + rt3o6 * Ch*(ks + 1), rt6o3*Ch*(ks + 1) + r), Cr));
								const_cast<object*>(G6_IST.Objs.back())->setcolor(LightSkyBlue);
							}
							if (js != 0) {
								G6_IST.add(new cylinder(point(Ch*is + Cho2 * js + Cho2 * ks, rt3o2*Ch*js + rt3o6 * Ch*ks, rt6o3*Ch*ks + r),
									point(Ch*is + Cho2 * (js - 1) + Cho2 * (ks + 1), rt3o2*Ch*(js - 1) + rt3o6 * Ch*(ks + 1), rt6o3*Ch*(ks + 1) + r), Cr));
								const_cast<object*>(G6_IST.Objs.back())->setcolor(LightGreen);
								G6_IST.add(new cylinder(point(Ch*is + Cho2 * js + Cho2 * ks, rt3o2*Ch*js + rt3o6 * Ch*ks, rt6o3*Ch*ks + r),
									point(Ch*(is + 1) + Cho2 * (js - 1) + Cho2 * ks, rt3o2*Ch*(js - 1) + rt3o6 * Ch*ks, rt6o3*Ch*ks + r), Cr));
								const_cast<object*>(G6_IST.Objs.back())->setcolor(LightGreen);
							}
							if (is != 0) {
								G6_IST.add(new cylinder(point(Ch*is + Cho2 * js + Cho2 * ks, rt3o2*Ch*js + rt3o6 * Ch*ks, rt6o3*Ch*ks + r),
									point(Ch*(is - 1) + Cho2 * js + Cho2 * (ks + 1), rt3o2*Ch*js + rt3o6 * Ch*(ks + 1), rt6o3*Ch*(ks + 1) + r), Cr));
								const_cast<object*>(G6_IST.Objs.back())->setcolor(Pink);
							}
							G6_IS.insert(new World(G6_IST));
							while (!G6_IST.Objs.empty()) delete G6_IST.Objs.back(), G6_IST.Objs.pop_back();
						}
					}
				}
				G6_I.insert(new World(G6_IS));
				for (int i = 0; i < G6_IS.GObjs.size(); i++) delete G6_IS.GObjs.at(i);
				G6_IS.GObjs.clear();
			}
		}
	}
	G6.insert(&G6_I);
	G6.add(&G6_G);
	G6.setGlobalLightSource(1.6, 0.2);


	WaterSurface water(20); water.setAttCof(0.54, 0.05, 0.02); water.setIndex(1.2); G6.add(&water);
	G6.Render_Sampling = 3;

	bitmap canvas(1200, 800);
	G6.render(canvas, point(500, 200, 300), point(0, 0, 0), 0, 0.01);
	canvas.out("IMAGE\\RT.bmp");
}

// ring, memory test
void Render_GTest02() {
	/*
		Note: This object goes through the most wonderful optimization I have ever made.
		CPU:     44min => 5.4s
		Memory:  1.26GB => 56MB
		Although it's just for debug purpose, with some memory leaks.
	*/

#define fx(u,v) (cos(u)*(0.8 + 0.3*sin(v)))
#define fy(u,v) (sin(u)*(0.8 + 0.3*sin(v)))
#define fz(u,v) (0.3*cos(v) + 0.3)

	bitmap img(600, 400);
	unsigned NP;

#define DPS 3
#define DIF1 8
#define DIF2 10
#define DIF3 32

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
	World G4, G4_R;
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
	//G4_LD.reflect = rgb(255, 255, 100);	// bright yellow
	G4_LD.reflect = rgb(178, 102, 255);
	for (int i = 0; i < NP; i++) {
		//W->reflect = rgblight(0.8, 0.8, 0.8), W->absorb = rgblight(0.2, 0.2, 0.2);	// silver
		//W->reflect = rgblight(rgb(218, 165, 32)), W->absorb = rgb(15, 115, 25);	// gold
		//W->reflect = rgblight(rgb(50, 255, 153));	// cyan
		W->reflect = rgblight(rgb(51, 153, 255));	// blue
		W++;
	}
	W -= NP;
	G4.add(&G4_LD);
	G4.insert(&G4_R);
	G4.setGlobalLightSource(PI / 2, 0);

	G4.render(img, point(30, 60, 40), point(0, 0, 0), 0, 0.004);


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

	img.out("IMAGE\\RT.bmp");

}

// complex Γ function
#ifdef MY_COMPLEX_INC

#pragma warning(push, 0)
#include "D:\Coding\AboutMath\SuperCalculator\SuperCalculator\Matrix.h"
inline bool baddouble(double a) {
	return isnan(a) || 1 / a == 0;
}
void Render_GTest03() {
	bitmap img(600, 600);
	unsigned NP;

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
#define Min_Y -3.8
#define Max_Y 3.8
#define X_Dif 192
#define Y_Dif 192
#define XC_Dif 8
#define YC_Dif 8
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
					W++;

				}
			}

		}
	}
	W -= NP;

	cout << "Packaging Data for Rendering..." << endl;

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
	G5.setGlobalLightSource(-1, -2, 1);
	G5.render(img, point(-30, 50, 55), point(-0.6, 0, 3.6), 0, 0.04);
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

	img.out("IMAGE\\RT.bmp");

}
#pragma warning(pop)

#endif

void Render_GTest04() {
	World W;
	//plane_grid PB(-3); W.add(&PB);
	plane_dif PB(-3); PB.setcolor(LightBlue); W.add(&PB);
	//WaterSurface WB(-1); WB.setAttCof(0.54, 0.05, 0.02); WB.setIndex(1.33); W.add(&WB);

	World W1;
	parallelogram G({ -2,-2,0 }, { 2,-2,0 }, { -2,2,0 }, 1); G.setcolor(Silver); W1.add(&G);
	cylinder C1(point(+1.5, +1.5, 0.1), point(+1.5, +1.5, -3), 0.3); C1.setcolor(Red); W1.add(&C1);
	cylinder C2(point(-1.5, +1.5, 0.1), point(-1.5, +1.5, -3), 0.3); C2.setcolor(Red); W1.add(&C2);
	cylinder C3(point(+1.5, -1.5, 0.1), point(+1.5, -1.5, -3), 0.3); C3.setcolor(Red); W1.add(&C3);
	cylinder C4(point(-1.5, -1.5, 0.1), point(-1.5, -1.5, -3), 0.3); C4.setcolor(Red); W1.add(&C4);
	cylinder C0(point(0, 0, 0.2), point(0, 0, -1.2), 0.5); C0.setcolor(Brown); W1.add(&C0);
	circle CC1(point(0, 0, 0.2), 0.5); CC1.setcolor(Brown); W1.add(&CC1);
	circle CC2(point(0, 0, -1.2), 0.5); CC2.setcolor(Brown); W1.add(&CC2);
	sphere B1(point(0, 0, 0.3), 0.1); B1.setcolor(SeaShell); W1.add(&B1);
	torus Rg(point(0, 0, 0.12), 1.2, 0.12); Rg.setcolor(Gold); W1.add(&Rg);
	W.insert(&W1);

	World W2;
	const double Tl = 4, Th = 0.8;
	parallelogram_ref T1(point(-0.5 * Tl, -0.5 * Tl, 0), point(0, Tl), point(Tl, 0));
	parallelogram_ref T2(point(-0.5 * Tl, -0.5 * Tl, Th), point(Tl, 0), point(0, Tl));
	parallelogram_ref T3(point(-0.5 * Tl, -0.5 * Tl, 0), point(Tl, 0), point(0, 0, Th));
	parallelogram_ref T4(point(-0.5 * Tl, -0.5 * Tl, Th), point(0, Tl), point(0, 0, -Th));
	parallelogram_ref T5(point(0.5 * Tl, 0.5 * Tl, 0), point(0, 0, Th), point(0, -Tl));
	parallelogram_ref T6(point(0.5 * Tl, 0.5 * Tl, 0), point(-Tl, 0), point(0, 0, Th));
	polyhedron T({ &T1, &T2, &T3, &T4, &T5, &T6 }); T.setAttCof(0.5, 0.2, 0.1); T.setIndex(1.5); W2.add(&T);
	//W.insert(&W2);

	fout << W << endl;

	bitmap img(3000, 2000);
	W.setGlobalLightSource(-1, -1, 2);
	W.background = rgblight(1);
	W.Render_Sampling = 8;
	W.render(img, point(30, 60, 40), point(0, 0, -1), 0, 0.01);
	img.out("IMAGE\\RT.bmp");
}

// water and two "pillars"
void Render_CTest01() {
	World W;
	plane B(-200); B.setcolor(rgb(255, 153, 0)); B.setcolor(Gray); W.add(&B);
	//plane_dif B(-80); B.setcolor(rgb(255, 153, 0)); B.setcolor(Gray); W.add(&B); W.Render_Sampling = 4;
	parallelogram P1(point(0, -2, -1), point(3, -2, -1), point(0, 0, -1), 1);
	P1.setcolor(LightGreen);
	parallelogram P2(point(0, -2, -200), point(3, -2, -200), point(0, -2, -1), 1);
	P2.setcolor(Gray);
	parallelogram P3 = P2 + point(0, 2, 0);
	parallelogram P4(point(0, -2, -200), point(0, -2, -1), point(0, 0, -200), 1);
	P4.setcolor(Gray);
	parallelogram P5 = P4 + point(3, 0, 0);
	W.add({ &P1, &P2, &P3, &P4, &P5 });
	parallelogram L1(point(-18, 8, 0.4), point(8, 8, 0.4), point(-18, 14, 0.4), 1);
	L1.setcolor(White);
	parallelogram L2(point(-18, 8, 0.4), point(8, 8, 0.4), point(-18, 8, -200), 1);
	L2.setcolor(Gray);
	parallelogram L3 = L2 + point(0, 6, 0);
	W.add({ &L1, &L2, &L3 });
	WaterSurface water(0); water.setAttCof(0.054, 0.005, 0.002); water.setIndex(1.33); W.add(&water);
	//for (int i = 0; i < W.Objs.size(); i++) cout << *W.Objs.at(i) << endl;

	sphere SO(point(0, 0, 0), 1); SO.setcolor(Silver);
	sphere SX(point(10, 0, 0), 1); SX.setcolor(Red);
	sphere SY(point(0, 10, 0), 1); SY.setcolor(Green);
	sphere SZ(point(0, 0, 10), 1); SZ.setcolor(Blue);
	//W.add({ &SO, &SX, &SY, &SZ });

	bitmap canvas(600, 400);
	//W.setGlobalLightSource(0, 0, 1); W.render(canvas, point(200, -200, -40), point(0, 0, 0), 0, 0.06);
	W.setGlobalLightSource(0, 2, 1); W.render(canvas, point(200, -200, 100), point(0, 0, 0), 0, 0.006);
	canvas.out("IMAGE\\RT.bmp");

	return;
	for (string name = "000"; name[0] < '2'; name[2]++) {
		if (name[2] > '9') name[1]++, name[2] = '0';
		if (name[1] > '9') name[0]++, name[1] = '0';
		W.render(canvas, point(200, -200, (name[0] - '0') * 100 + (name[1] - '0') * 10 + name[2] - '0'), point(0, 0, 0), 0, 0.006);
		canvas.out("IMAGE\\T\\RT_z." + name + ".bmp");
	}
}

// sphere-water test
void Render_CTest02() {
	World W;
	plane P(-1); P.setcolor(LightSkyBlue); W.add(&P);
	WaterSurface PW(-ERR_ZETA); PW.setAttCof(0.54, 0.05, 0.02); PW.setIndex(1.33); W.add(&PW);
	//sphere A(point(0, 0, 1), 1); A.setcolor(Red); W.add(&A);

	sphere3D A3(point(0, 0, 2), 2); A3.setAttCof(0.5, 0.1, 0.2); A3.setIndex(1.5); W.add(&A3);

	vector<parallelogram> Pr;
	parallelogram Pr0(point(-4, -4, 0), point(1, 0, 0), point(0, 1, 0));
	for (int i = 0; i < 8; i++) {
		for (int j = 0; j < 8; j++) {
			Pr0.setcolor((i + j) & 1 ? LightBlue : Gray);
			Pr.push_back(Pr0);
			Pr0 += point(0, 1, 0);
		}
		Pr0 += point(1, -8, 0);
	}
	World W1;
	for (int i = 0; i < Pr.size(); i++) W1.add(&Pr.at(i));
	W.insert(&W1);

	bitmap img(600, 400);
	W.Render_Sampling = 2;
	W.setGlobalLightSource(0, 0, 1);
	W.render(img, point(100, 30, 50), point(0, 0, 0), 0, 0.02);
	img.out("IMAGE\\RT.bmp");
}

// sphere-word test and Moana scene
void Render_CTest03() {
	World W1, W2;
	plane_grid P1(0); W1.add(&P1);

	bitmap_inc WT(bitmap("D:\\Coding\\AboutPhysics\\RayTracing\\RayTracing\\IMAGE\\WT.bmp"), point(0, 0, ERR_ZETA), point(4, 0), point(0, 1), insertType::center | insertType::proportion);
	W1.add(&WT);
	// screenshot: <i style="margin:50px;font-family:impact;font-size:128px;color:black;">Ray-Tracing</i>

	sphere3D S1(point(0, 0, 20), 2); S1.setAttCof(0.5, 0.1, 0.2); S1.setIndex(1.1); W1.add(&S1);

	bitmap img(600, 400);
	W1.setGlobalLightSource(0, 0, 1);
	for (string i = "00"; i[0] < '3'; i[1]++) {
		S1.C.z = (i[0] - '0') * 10 + i[1] - '0';
		//W1.render(img, point(0, 0, 100), point(0, 0, 0), 0, 0.01);
		//img.out("IMAGE\\T\\RT" + i + ".bmp");
		if (i[1] == '9') i[0]++, i[1] = '0' - 1;
	}
	//return;

	bitmap_inc BM(bitmap("D:\\Coding\\AboutImage\\ColorAnalysisTest\\2DFitting\\Origin\\SK01.bmp"), point(0, 0, 0), point(19.1, 0, 0), point(0, 0, 8)); W2.add(&BM);
	// capture from 3D animation film Moana, 16:20, a girl standing on the rock and staring at the sea

	S1.C = point(0, 0, 0), S1.r = 2; S1.setIndex(1.5); W2.add(&S1);
	sphere3D S2(point(8, -5, 4), 2); S2.setAttCof(0.5, 0.1, 0.2); S2.setIndex(1.5); W2.add(&S2);
	plane_grid P2(-2, 2); W2.add(&P2);
	WaterSurface PW(0, 1.33); PW.setAttCof(0.54, 0.05, 0.02); W2.add(&PW);
	img = bitmap(8000, 3500);
	W2.setGlobalLightSource(0, -1, 1);
	W2.render(img, point(0, -100, 40), point(8, 0, 2.4), 0, 0.04); img.out("IMAGE\\RT.bmp");	// it would be a very wonderful picture

	return;

	img = bitmap(600, 400);
	S1.setIndex(1.1);
	for (string name = "000"; name[1] < '2'; name[2]++) {
		if (name[2] > '9') name[1]++, name[2] = '0';
		if (name[1] > '9') name[0]++, name[1] = '0';
		W1.render(img, point(0, 0, 100), point(0, 0, 0), 0, 0.02);
		img.out("IMAGE\\T\\RT_z." + name + ".bmp");
		S1.C.z++;
	}
}

// triangular prism
void Render_CTest04() {
	World W;
	triangle_ref T1(point(0, 0, 0), point(0, 1, 0), point(1, 0, 0));
	parallelogram_ref T2(point(0, 0, 0), point(1, 0, 0), point(0, 0, 1), 1);
	parallelogram_ref T3(point(0, 0, 0), point(0, 0, 1), point(0, 1, 0), 1);
	parallelogram_ref T4(point(1, 0, 0), point(0, 1, 0), point(1, 0, 1), 1);
	triangle_ref T5(point(0, 0, 1), point(1, 0, 1), point(0, 1, 1));
	double cn = 1.3; T1.setIndex(cn), T2.setIndex(cn), T3.setIndex(cn), T4.setIndex(cn), T5.setIndex(cn);
	point P(0, 0, ERR_ZETA); T1 += P, T2 += P, T3 += P, T4 += P, T5 += P;
	polyhedron T({ &T1, &T2, &T3, &T4, &T5 }); sphere3D S(point(0.5, 0.2, 0.5), 0.3);
	rgblight cf(1.08, 0.10, 0.04); T.attcoe = S.attcoe = cf;
	W.add(&T); W.add(&S);
	plane_grid G(0.0); W.add(&G);
	//plane_dif G(0.0); G.setcolor(DarkSeaGreen); W.add(&G);
	parallelogram G1(point(-1, 0, ERR_ZETA), point(cos(2.8), sin(2.8)), 2 * point(cos(2.8 - PI / 2), sin(2.8 - PI / 2))); G1.setcolor(Brown); W.add(&G1);
	fout << W << endl;

	W.Render_Sampling = 2;
	W.setGlobalLightSource(0, 0, 1);
	bitmap img(600, 400);
	W.render(img, point(10, -6, 5), point(0, 0, 0), 0, 0.06);
	img.out("IMAGE\\RT.bmp");
}

// "chessboard"
void Render_LTest01() {
	World W;
	spherebulb L(point(0, 0, 3), 0.8, rgblight(1, 1, 2)); W.add(&L);
	//sphere L(point(0, 0, 3), 0.8); L.setcolor(rgblight(1, 1, 2)); W.add(&L);
	sphere B(point(1, 1.2, 1), 1); B.setcolor(Red); W.add(&B);
	plane_dif Pb(point(-2, 0, 0), point(1, -0.8, 0)); Pb.setcolor(LightYellow); Pb.setvar(0.1); W.add(&Pb);
	plane_dif P(-1); P.setcolor(LightBlue); P.setvar(0.1); W.add(&P);
	rectbulb Pr0(point(-4, -4, 0), point(1, 0, 0), point(0, 1, 0)); Pr0.setcolor(LightBlue); //W.add(&Pr0);
	vector<rectbulb> Pr;
	W.Render_Sampling = 8;
	for (int i = 0; i < 8; i++) {
		for (int j = 0; j < 8; j++) {
			Pr0.setcolor((i + j) & 1 ? LightBlue : Gray);
			Pr.push_back(Pr0);
			Pr0 += point(0, 1, 0);
		}
		Pr0 += point(1, -8, 0);
	}
	for (int i = 0; i < Pr.size(); i++) W.add(&Pr.at(i));
	//for (int i = 0; i < 2; i++) W.add(&Pr.at(i));
	//W.setGlobalLightSource(PI / 2, 0);

	bitmap canvas(600, 400);
	W.render(canvas, point(200, -200, 200), point(0, 0, 0), 0, 0.001);
	canvas.out("IMAGE\\RT.bmp");

}

// two spheres inside a room
void Render_LTest02() {
	World W;
	plane_dif B(0.0); B.setcolor(Gray); B.setvar(0.1);
	plane_dif T(4.0); T.setcolor(Gray); T.setvar(0.1);
	plane_dif L(point(0, -2, 0), point(0, 1, 0)); L.setcolor(Orange); L.setvar(0.1);
	plane_dif R(point(0, 2, 0), point(0, -1, 0)); R.setcolor(SkyBlue); R.setvar(0.1);
	plane_dif Bk(point(-6, 0, 0), point(1, 0, 0)); Bk.setcolor(LightGray); Bk.setvar(0.1);
	sphere S1(point(0, 2, 1), 1); S1.setcolor(Gray);
	sphere S2(point(4, -1, 0.6), 0.6); S2.setcolor(Gray);
	spherebulb Lb(point(0, 0, 4.2), 1); Lb.setcolor(White); Lb.setcolor(rgblight(10));
	W.add({ &B, &T, &L, &R, &Bk, &S1, &S2, &Lb });
	//W.setGlobalLightSource(PI / 2, 0);

	bitmap canvas(600, 400);
	W.Render_Sampling = 12;
	W.render(canvas, point(12, 0, 1), point(0, 0, 1), 0, 0.6);
	canvas.out("IMAGE\\RT.bmp");
}

// a house on the grassland near a pool
void Render_LTest03() {
	World W;

	point P1[7] = { point(0,0), point(0,4), point(1.2,4), point(2,3.2), point(2,0), point(0,0), point(0,4) };
	parallelogram_dif P[5], B[5];
	for (int i = 0; i < 5; i++) {
		P[i] = parallelogram(P1[i], P1[i + 1] - P1[i], point(0, 0, -1));
		B[i] = parallelogram(P1[i + 1], ERR_UPSILON*(P1[i + 2] - P1[i + 1]), ERR_UPSILON*(P1[i + 1] - P1[i]));
		P[i].setcolor(SeaShell), B[i].setcolor(SeaGreen);
		W.add(&P[i]); W.add(&B[i]);
	}

	point KP[5] = { point(2.8,1), point(2.8,2.8), point(3.2,2.8), point(3.2,1), point(2.8,1) };
	parallelogram_dif K[4], K1, K2;
	for (int i = 0; i < 4; i++) {
		K[i] = parallelogram(KP[i], KP[i + 1] - KP[i], point(0, 0, 0.24));
		K[i].setcolor(Brown);
		W.add(&K[i]);
	}
	K1 = parallelogram(KP[0], KP[1], 0.5*(KP[0] + KP[3]) + point(0, 0, 0.08), 1) + point(0, 0, 0.24); K1.setcolor(Brown); W.add(&K1);
	K2 = parallelogram(KP[3], KP[2], 0.5*(KP[0] + KP[3]) + point(0, 0, 0.08), 1) + point(0, 0, 0.24); K2.setcolor(Brown); W.add(&K2);
	triangle_dif T1(0.5*(KP[0] + KP[3]) + point(0, 0, 0.08), KP[0], KP[3]); T1 += point(0, 0, 0.24); T1.setcolor(Brown); W.add(&T1);
	triangle_dif T2(0.5*(KP[1] + KP[2]) + point(0, 0, 0.08), KP[1], KP[2]); T2 += point(0, 0, 0.24); T2.setcolor(Brown); W.add(&T2);

	plane_dif F(-1); F.setcolor(LightGray); W.add(&F);
	WaterSurface WF(-0.1); WF.setAttCof(0.54, 0.05, 0.02); WF.setIndex(1.33); W.add(&WF);

	//sphere S1(0.5*(P1[2] + P1[3]) + point(0, 0, -1), 0.08); S1.setcolor(Red); W.add(&S1);

	W.setGlobalLightSource(0.2, 0.2, 1);
	W.Render_Sampling = 16;


	bitmap canvas(600, 400);
	W.render(canvas, point(-1, -3, 4), point(1.8, 1.8, 0), 0, 0.6);
	canvas.out("IMAGE\\RT.bmp");
}


#include "CSG.h"

void Render_XTest00() {
	World W;
	plane_grid P(0); W.add(&P);
	XObjs::Sphere XS1(point(0, 0, 1), 1);
	XObjs::Sphere XS2(point(0, 1, 1), 1);
	XObjs::Plane XP1(point(0, 0.5, 1), point(-0.4, 1, 1));
	XObjs::Plane XP2(point(0, 0, 1.2), point(0, 0, 1));
	XObjs::Cylinder_std XCSx(point(0, 0, 1), 1, point(1, 0, 0));
	XObjs::Cylinder_std XCSy(point(0, 0, 1), 1, point(0, 1, 0));
	XObjs::Cylinder_std XCSz(point(0, 0, 1), 1, point(0, 0, 1));
	XObjs::Cylinder XC1(point(0, 0, 0.5), point(0, 1, 1.5), 0.5);
	XObjs::Cone_std XCNS1(point(0, 0, 0), point(0, 0, 1), 0.5);
	XObjs::Cone_std XCNS2(point(0, 0, 1), point(0, 0, 0), 0.5);
	XObjs::Cone XCN1(point(1, 0, 0.1), point(-4, 0, 1), 2, atan(0.25));
	XObjs::Cone_capped XCC1(point(0, 0, 0.1001), point(0, 0, 1), 0.6, 0.2);
	XObjs::Torus_xOy XT1(point(0, 0, 1), 3, 1);
	XObjs::Box_xOy XB1(point(1, 1, 1.2), point(-2, 0, 0.5));
	XObjs::Box_affine XX1(matrix3D_affine(
		2, 1, 0, 0,
		-1, 1, 0, 0,
		0, 0, 1, 1,
		0.6, 0, 0, 1));
	XSolid X;
	//X = CSG_SubtractionOp(CSG_IntersectionOp(XSolid(XS1), XSolid(XS2)), XSolid(XP1));
	//X = CSG_IntersectionOp(CSG_IntersectionOp(XSolid(XCSx), XSolid(XCSy)), XSolid(XCSz));
	//X = CSG_IntersectionOp(XSolid(XS1), XSolid(XCNS1));
	//X = CSG_IntersectionOp(CSG_OnionOp(XSolid(XS1), 0.1), XSolid(XP2));
	X = XX1;
	X = CSG_RoundingOp(X, 0.2);
	//X = CSG_OnionOp(X, 0.1);
	X.setColor(White); W.add(&X);

	X.type = XSolid_Crystal; X.setColor(rgblight(0, 0, 0));
	X.col = rgblight(0.5, 0.2, 0.1);

	sphere3D S(point(-1.5, 1.5, 1), 1); S.C = point(0, 0, 1);
	S.setIndex(1.5); S.setAttCof(0.8, 0.4, 0.2);
	//W.add(&S);

	parallelogram Pr(point(3, 0.5, -1), point(-6, 0, 0), point(0, 0, 4));
	Pr = parallelogram(point(-2.5, -1, -1), point(6, 0, 0), point(0, 3, 3));
	//Pr.setcolor(White); W.add(&Pr);


	//ADD_AXIS(W, 0.05);
	VisualizeSDF(X, Pr, 600 * 400);

	bitmap img(600, 400);
	W.setGlobalLightSource(0, 0, 1);
	//W.Render_Sampling = 2;
	W.render(img, point(10, -10, 10), point(pick_random(-0.001, 0.001), 0, 0), 0, 0.1);
	img.out("IMAGE\\RT.bmp");
}
