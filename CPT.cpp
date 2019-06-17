#pragma once

#include "World.h"
#include "GlassMan.h"
#include "Test.cpp"

using namespace std;

#define CANVAS_WIDTH 192
#define CANVAS_HEIGHT 108
#define RENDER_SAMPLING 1


void CPT_T1() {
	World W;
	plane_grid P(0); W.add(&P);

	GlassMan G;
	//G.top = point(0, 0, 1), G.dir = point(-1, 1, 0);

	G.Heel_l = point(0, 1, 0.3), G.Heel_r = point(0, -1, 0.5);
	G.foot_dir_l = point(2, 0.2, -0.3), G.foot_dir_r = point(2, -0.1, -0.5); G.foot_side_l = point(0, 1, 0), G.foot_side_r = point(0, -1, 0);
	G.width_of_heel = 1.0, G.width_of_tiptoe = 0.8, G.height_of_heel = 0.6; G.foot_rounding = 0.15;
	G.length_of_foot = 3;

	G.Knee_l = point(0.3, 1, 4.8), G.Knee_r = point(0.5, -0.9, 4.8);
	G.upper_radius_of_shank = 0.4, G.lower_radius_of_shank = 0.24, G.shank_rounding = 0.2;

	G.Buttock_l = point(0, 0.7, 7.8), G.Buttock_r = point(0, -0.6, 7.8);
	G.upper_radius_of_thigh = 0.6, G.lower_radius_of_thigh = 0.4, G.thigh_rounding = 0.35;

	G.Waist = point(0, 0, 10); G.lower_width_of_waist = 2.8, G.width_of_waist = 2.2, G.lower_depth_of_waist = 1.4; G.waist_rounding = 0.6;

	G.Shoulder_l = point(0, 1.7, 13), G.Shoulder_r = point(0, -1.7, 13); G.depth_of_chest = 1.4; G.chest_rounding = 0.5; G.upper_radius_of_upperarm = 0.4;

	G.Elbow_l = point(0.2, 1.8, 9.5), G.Elbow_r = point(0.6, -1.7, 9.5); G.lower_radius_of_upperarm = G.upperarm_rounding = 0.3;

	G.Hand_l = point(0.7, 1.6, 7.5), G.Hand_r = point(0.9, -1.2, 7.8); G.upper_radius_of_forearm = 0.3, G.lower_radius_of_forearm = 0.25, G.forearm_rounding = 0.25;

	G.Head_lower = point(0, 0, 13.8), G.Head_top = point(0, 0, 16.5); G.Head_side = vec3(0, 1, 0);
	G.width_of_forehead = 2.4, G.width_of_chin = 1.8, G.depth_of_forehead = 1.4; G.head_rounding = 0.699;
	G.neck_radius = 0.25, G.neck_rounding = 0.1;


	//sphere SH1(G.Heel_l, 0.05); SH1.setcolor(Red); W.add(&SH1);
	//sphere SH2(G.Heel_r, 0.05); SH2.setcolor(Red); W.add(&SH2);
	//sphere SK1(G.Knee_l, 0.05); SK1.setcolor(Red); W.add(&SK1);
	//sphere SK2(G.Knee_r, 0.05); SK2.setcolor(Red); W.add(&SK2);
	//sphere SB1(G.Buttock_l, 0.05); SB1.setcolor(Red); W.add(&SB1);
	//sphere SB2(G.Buttock_r, 0.05); SB2.setcolor(Red); W.add(&SB2);
	//sphere SW(G.Waist, 0.05); SW.setcolor(Red); W.add(&SW);
	//sphere SS1(G.Shoulder_l, 0.05); SS1.setcolor(Red); W.add(&SS1);
	//sphere SS2(G.Shoulder_r, 0.05); SS2.setcolor(Red); W.add(&SS2);
	//sphere SE1(G.Elbow_l, 0.05); SE1.setcolor(Red); W.add(&SE1);
	//sphere SE2(G.Elbow_r, 0.05); SE2.setcolor(Red); W.add(&SE2);
	//sphere SHd1(G.Elbow_l, 0.05); SHd1.setcolor(Red); W.add(&SHd1);
	//sphere SHd2(G.Elbow_r, 0.05); SHd2.setcolor(Red); W.add(&SHd2);

	G.construct(); G.push(W);

	parallelogram Pr(point(-6, -4, 1), point(12, 0, 0), point(0, 8, 0)); //Pr.setcolor(White), W.add(&Pr);
	//VisualizeSDF(G.construct(), Pr, 600 * 400);
	//ScanXSolid(G.construct(), 20, 20, 0, 240000);

	bitmap img(400, 600);
	//ADD_AXIS(W, 0.1);
	W.setGlobalLightSource(0, 0, 1);
	W.Render_Sampling = 1;
	//W.render(img, point(200, -100, 100), point(0, 1, 0), 0, 0.0001);
	W.render(img, point(200, -100, 100), point(0, 0, 8), 0, 0.004);
	img.out("IMAGE\\CPT.bmp");
}

void CPT_T2() {
	World W;
	plane_grid P(0); W.add(&P);

	GlassMan_std G;

	G.Pos = point(0, 0, 0), G.top = vec3(0, 0, 1), G.dir = vec3(1, 0, 0);
	G.v_head = vec3(0, 0, 1), G.v_head_side = vec3(-sin(0.2), cos(0.2), 0), G.v_neck = vec3(0, 0, 1);
	G.v_chest = vec3(0, 0, 1), G.v_chest_side = vec3(-sin(0.1), cos(0.1), 0), G.v_waist = vec3(0, 0, -1), G.v_waist_side = vec3(0, 1, 0);
	G.v_upperarm_l = vec3(0, 0.1, -1), G.v_forearm_l = vec3(0.2, -0.2, -1); G.v_upperarm_r = vec3(-0.05, -0.15, -1), G.v_forearm_r = vec3(0.3, 0.2, -1);
	G.v_thigh_l = vec3(0.05, 0, -1), G.v_shank_l = vec3(0, 0, -1); G.v_thigh_r = vec3(0.15, -0.15, -1), G.v_shank_r = vec3(0.1, -0.1, -1);
	G.v_foot_l = vec3(cos(0.1), sin(0.1), -0.1), G.v_foot_l_side = vec3(0, 1, 0); G.v_foot_r = vec3(cos(0.7), -sin(0.7), 0), G.v_foot_r_side = vec3(-sin(0.7), -cos(0.7), 0);
	G.auto_fit = true;

	//ScanXSolid(G.construct(), 10, 10, 10, 240000);

	World W_ = W;
	G.push(W_);
	W.setGlobalLightSource(0, 0, 1);
	W.Render_Sampling = RENDER_SAMPLING;
	bitmap img(CANVAS_WIDTH, CANVAS_HEIGHT);
	W.render(img, point(200, -100, 100), point(0, 0, 0.8), 0, 0.00005);
	img.out("IMAGE\\CPT.bmp");

}

void CPT_Animation_T1() {
	const string HEADER = "Animation\\CPTAT1_";
	const unsigned N = 2;
	const int FRAMES = 25;
	const int THREADS = thread::hardware_concurrency();

	World W;
	plane_grid P(0); W.add(&P);
	GlassMan_std G;

	W.setGlobalLightSource(0, 0, 1);
	W.Render_Sampling = RENDER_SAMPLING;

	double t = 0;
	vector<World*> Ws; vector<point> Cs, Vs; vector<double> Va;
	vector<thread*> Ts; vector<char> Rs;	// Rs: 1 finished  0 not finished  -1 finished

	for (unsigned i = 0; i < THREADS; i++) {
		Ws.push_back(0), Cs.push_back(point()), Vs.push_back(point()), Va.push_back(0);
		Ts.push_back(0), Rs.push_back(1);
	}

	int finished = 0;
	while (1) {
		for (int i = 0; i < THREADS; i++) {
			if (Rs.at(i) == 1) {
				if (Ws.at(i) != 0) delete Ws.at(i), Ws.at(i) = 0;
				if (Ts.at(i) != 0) delete Ts.at(i), Ts.at(i) = 0; Rs.at(i) = 0;

				if (finished == -1) Rs.at(i) = -1;
				else if (++finished > FRAMES) Rs.at(i) = -1, finished = -1;	// all finished
				else {
					t = double(finished) / FRAMES;

					G.Pos = point(0, 0, 0), G.top = vec3(0, 0, 1), G.dir = vec3(cos(t), -sin(t), 0);
					G.v_head = vec3(0, 0, 1), G.v_head_side = vec3(-sin(0.2), cos(0.2), 0), G.v_neck = vec3(0, 0, 1);
					G.v_chest = vec3(0, 0, 1), G.v_chest_side = vec3(-sin(0.1), cos(0.1), 0), G.v_waist = vec3(0, 0, -1), G.v_waist_side = vec3(0, 1, 0);
					G.v_upperarm_l = vec3(0, 0.1, -1), G.v_forearm_l = vec3(0.2, -0.2, -1); G.v_upperarm_r = vec3(-0.05, -0.15, -1), G.v_forearm_r = vec3(0.3, 0.2, -1);
					G.v_thigh_l = vec3(0.05, 0, -1), G.v_shank_l = vec3(0, 0, -1); G.v_thigh_r = vec3(0.15, -0.15, -1), G.v_shank_r = vec3(0.1, -0.1, -1);
					G.v_foot_l = vec3(cos(0.1), sin(0.1), -0.1), G.v_foot_l_side = vec3(0, 1, 0); G.v_foot_r = vec3(cos(0.7), -sin(0.7), 0), G.v_foot_r_side = vec3(-sin(0.7), -cos(0.7), 0);
					G.auto_fit = true;

					W.clear(); W.add(&P); G.push(W);

					Ws.at(i) = new World(W);
					Cs.at(i) = point(200, -100, 120 * t + 30), Vs.at(i) = point(0, 0, 0.8), Va.at(i) = 0.00005;

					Ts.at(i) = new thread([&](World *_W, unsigned th, unsigned no, unsigned N) {
						bitmap img(CANVAS_WIDTH, CANVAS_HEIGHT);
						_W->render(img, Cs.at(th), Vs.at(th), 0, Va.at(th));
						img.out(&(HEADER + uint2str(no, N) + ".bmp")[0]);
						Sleep(50); Rs.at(th) = 1;
					}, Ws.at(i), i, finished, N);
					Ts.at(i)->detach();
				}
			}
		}
		if (finished == -1) {
			for (int i = 0; i < THREADS; i++) {
				if (Rs.at(i) != -1) break;
				if (i + 1 == THREADS) goto Finish;
			}
		}
		Sleep(50);
	}

Finish:;

	for (int i = 0; i < THREADS; i++) {
		if (Ws.at(i) != 0) delete Ws.at(i), Ws.at(i) = 0;
		if (Ts.at(i) != 0) {
			delete Ts.at(i), Ts.at(i) = 0;
		}
	}
}

void CPT_Animation_T2() {
	const string HEADER = "Animation\\CPTAT2_";
	const unsigned N = 2;
	const int FRAMES = 25;
	const int THREADS = thread::hardware_concurrency();

	World W;
	parallelogram_grid P(parallelogram(point(-2, -2, 0), point(4, 0), point(0, 4))); W.add(P);
	GlassMan_std G;

	W.setGlobalLightSource(0, 0, 1);
	W.Render_Sampling = RENDER_SAMPLING;

	double t = 0;
	vector<World*> Ws; vector<point> Cs, Vs; vector<double> Va;
	vector<thread*> Ts; vector<char> Rs;	// Rs: 1 finished  0 not finished  -1 finished

	for (unsigned i = 0; i < THREADS; i++) {
		Ws.push_back(0), Cs.push_back(point()), Vs.push_back(point()), Va.push_back(0);
		Ts.push_back(0), Rs.push_back(1);
	}

	int finished = 0;
	while (1) {
		for (int i = 0; i < THREADS; i++) {
			if (Rs.at(i) == 1) {
				if (Ws.at(i) != 0) delete Ws.at(i), Ws.at(i) = 0;
				if (Ts.at(i) != 0) delete Ts.at(i), Ts.at(i) = 0; Rs.at(i) = 0;

				if (finished == -1) Rs.at(i) = -1;
				else if (++finished > FRAMES) Rs.at(i) = -1, finished = -1;	// all finished
				else {
					t = double(finished) / FRAMES;

					G.Pos = point(0, 0, 0), G.top = vec3(0, 0, 1), G.dir = vec3(cos(t), -sin(t), 0);
					G.v_head = vec3(0, 0, 1), G.v_head_side = vec3(-sin(0.2), cos(0.2), 0), G.v_neck = vec3(0, 0, 1);
					G.v_chest = vec3(0, 0, 1), G.v_chest_side = vec3(-sin(0.1), cos(0.1), 0), G.v_waist = vec3(0, 0, -1), G.v_waist_side = vec3(0, 1, 0);
					G.v_upperarm_l = vec3(0, 0.1, -1), G.v_forearm_l = vec3(0.2, -0.2, -1); G.v_upperarm_r = vec3(-0.05, -0.15, -1), G.v_forearm_r = vec3(0.3, 0.2, -1);
					G.v_thigh_l = vec3(0.05, 0, -1), G.v_shank_l = vec3(0, 0, -1); G.v_thigh_r = vec3(0.15, -0.15, -1), G.v_shank_r = vec3(0.1, -0.1, -1);
					G.v_foot_l = vec3(cos(0.1), sin(0.1), -0.1), G.v_foot_l_side = vec3(0, 1, 0); G.v_foot_r = vec3(cos(0.7), -sin(0.7), 0), G.v_foot_r_side = vec3(-sin(0.7), -cos(0.7), 0);
					G.auto_fit = true;

					W.clear(); W.add(&P); G.push(W);

					Ws.at(i) = new World(W);
					Cs.at(i) = point(200, -100, 120 * t + 30), Vs.at(i) = point(0, 0, 0.8), Va.at(i) = 0.00005;

					Ts.at(i) = new thread([&](World *_W, unsigned th, unsigned no, unsigned N) {
						bitmap img(CANVAS_WIDTH, CANVAS_HEIGHT);
						_W->render(img, Cs.at(th), Vs.at(th), 0, Va.at(th));
						img.out(&(HEADER + uint2str(no, N) + ".bmp")[0]);
						Sleep(50); Rs.at(th) = 1;
					}, Ws.at(i), i, finished, N);
					Ts.at(i)->detach();
				}
			}
		}
		if (finished == -1) {
			for (int i = 0; i < THREADS; i++) {
				if (Rs.at(i) != -1) break;
				if (i + 1 == THREADS) goto Finish;
			}
		}
		Sleep(50);
	}

Finish:;

	for (int i = 0; i < THREADS; i++) {
		if (Ws.at(i) != 0) delete Ws.at(i), Ws.at(i) = 0;
		if (Ts.at(i) != 0) {
			delete Ts.at(i), Ts.at(i) = 0;
		}
	}
}


void CPT_Animation(void(*s_World)(double, World*), point(*s_Camera)(double), point(*s_ViewPoint)(double), double(*s_SolidAngle)(double),
	double t0, double t1, double canvas_width, double canvas_height, unsigned fps, unsigned render_sampling, string filename_header, unsigned filename_reserve) {

	const int FRAMES = fps * (t1 - t0);
	const int THREADS = thread::hardware_concurrency();

	World W;
	W.setGlobalLightSource(0, 0, 1);
	W.Render_Sampling = render_sampling;

	double t = 0;
	vector<World*> Ws; vector<point> Cs, Vs; vector<double> Va;
	vector<thread*> Ts; vector<char> Rs;	// Rs: 1 finished  0 not finished  -1 finished

	for (unsigned i = 0; i < THREADS; i++) {
		Ws.push_back(0), Cs.push_back(point()), Vs.push_back(point()), Va.push_back(0);
		Ts.push_back(0), Rs.push_back(1);
	}

	int finished = 0;
	while (1) {
		for (int i = 0; i < THREADS; i++) {
			if (Rs.at(i) == 1) {
				if (Ws.at(i) != 0) delete Ws.at(i), Ws.at(i) = 0;
				if (Ts.at(i) != 0) delete Ts.at(i), Ts.at(i) = 0; Rs.at(i) = 0;

				if (finished == -1) Rs.at(i) = -1;
				else if (++finished > FRAMES) Rs.at(i) = -1, finished = -1;	// all finished
				else {
					t = mix(t0, t1, double(finished) / FRAMES);

					W.clear(); s_World(t, &W);

					Ws.at(i) = new World(W);
					Cs.at(i) = s_Camera(t), Vs.at(i) = s_ViewPoint(t), Va.at(i) = s_SolidAngle(t);

					Ts.at(i) = new thread([&](World *_W, unsigned th, unsigned no, unsigned N) {
						bitmap img(canvas_width, canvas_height);
						_W->render(img, Cs.at(th), Vs.at(th), 0, Va.at(th));
						img.out(&(filename_header + uint2str(no, filename_reserve) + ".bmp")[0]);
						Rs.at(th) = 1;
					}, Ws.at(i), i, finished, filename_reserve);
					Ts.at(i)->detach();
				}
			}
		}
		if (finished == -1) {
			for (int i = 0; i < THREADS; i++) {
				if (Rs.at(i) != -1) break;
				if (i + 1 == THREADS) goto Finish;
			}
		}
		this_thread::sleep_for(chrono::milliseconds(50));
	}

Finish:;

	for (int i = 0; i < THREADS; i++) {
		if (Ws.at(i) != 0) delete Ws.at(i), Ws.at(i) = 0;
		if (Ts.at(i) != 0) {
			delete Ts.at(i), Ts.at(i) = 0;
		}
	}
}

void CPT_Animation_T3() {
	CPT_Animation([](double t, World* W) {	// World
		parallelogram_grid P(parallelogram(point(-3, -ERR_UPSILON, 0), point(7, 0), point(0, 2 * ERR_UPSILON))); W->add(P);
		GlassMan_std G;
		G.Pos = point(0, 0, 0), G.top = vec3(0, 0, 1), G.dir = vec3(cos(t), -sin(t), 0);
		G.v_head = vec3(0, 0, 1), G.v_head_side = vec3(-sin(0.2), cos(0.2), 0), G.v_neck = vec3(0, 0, 1);
		G.v_chest = vec3(0, 0, 1), G.v_chest_side = vec3(-sin(0.1), cos(0.1), 0), G.v_waist = vec3(0, 0, -1), G.v_waist_side = vec3(0, 1, 0);
		G.v_upperarm_l = vec3(0, 0.1, -1), G.v_forearm_l = vec3(0.2, -0.2, -1); G.v_upperarm_r = vec3(-0.05, -0.15, -1), G.v_forearm_r = vec3(0.3, 0.2, -1);
		G.v_thigh_l = vec3(0.05, 0, -1), G.v_shank_l = vec3(0, 0, -1); G.v_thigh_r = vec3(0.15, -0.15, -1), G.v_shank_r = vec3(0.1, -0.1, -1);
		G.v_foot_l = vec3(cos(0.1), sin(0.1), -0.1), G.v_foot_l_side = vec3(0, 1, 0); G.v_foot_r = vec3(cos(0.7), -sin(0.7), 0), G.v_foot_r_side = vec3(-sin(0.7), -cos(0.7), 0);
		G.auto_fit = true;
		G.push(*W);
	}, [](double t) -> point {	// camera
		return mix(point(10, -20, 15), point(-4, -30, 2), t);
	}, [](double t) -> point {	// view point
		return point(0, 0, 0.8 - 0.8 * t);
	}, [](double t) -> double {		// solid angle
		return 0.2 / (2 * t + 1);
	}, 0, 1, 1920, 1080, 60, 2, "Animation\\CPTAT3_", 3);
}

#include "HumanWalkingData.h"

void CPT_Animation_T4() {
	CPT_Animation([](double t, World* W) {	// World
		double Pos_m, Pos_b; WalkingMan::LinearRegression(WalkingMan::waist, 0, Pos_m, Pos_b);
		const unsigned N = 6;
		vec3 thigh_l_a[N], thigh_l_b[N]; WalkingMan::FourierSeries(WalkingMan::v_thigh_l, N, thigh_l_a, thigh_l_b);
		vec3 thigh_r_a[N], thigh_r_b[N]; WalkingMan::FourierSeries(WalkingMan::v_thigh_r, N, thigh_r_a, thigh_r_b);
		vec3 shank_l_a[N], shank_l_b[N]; WalkingMan::FourierSeries(WalkingMan::v_shank_l, N, shank_l_a, shank_l_b);
		vec3 shank_r_a[N], shank_r_b[N]; WalkingMan::FourierSeries(WalkingMan::v_shank_r, N, shank_r_a, shank_r_b);
		vec3 foot_l_a[N], foot_l_b[N]; WalkingMan::FourierSeries(WalkingMan::v_foot_l, N, foot_l_a, foot_l_b);
		vec3 foot_r_a[N], foot_r_b[N]; WalkingMan::FourierSeries(WalkingMan::v_foot_r, N, foot_r_a, foot_r_b);
		vec3 upperarm_r_a[N], upperarm_r_b[N]; WalkingMan::FourierSeries(WalkingMan::v_upperarm_r, N, upperarm_r_a, upperarm_r_b);
		vec3 forearm_r_a[N], forearm_r_b[N]; WalkingMan::FourierSeries(WalkingMan::v_forearm_r, N, forearm_r_a, forearm_r_b);
		vec3 chest_a[N], chest_b[N]; WalkingMan::FourierSeries(WalkingMan::v_chest, N, chest_a, chest_b);
		vec3 waist_a[N], waist_b[N]; WalkingMan::FourierSeries(WalkingMan::v_waist, N, waist_a, waist_b);
		vec3 head_a[N], head_b[N]; WalkingMan::FourierSeries(WalkingMan::v_head, N, head_a, head_b);
		vec3 neck_a[N], neck_b[N]; WalkingMan::FourierSeries(WalkingMan::v_neck, N, neck_a, neck_b);
		auto FourierEval = [](double t, vec3 *a, vec3 *b, int n) -> vec3 {
			const double omega = 2 * PI / 1.2; double nwt;
			vec3 r; r.x = r.y = r.z = 0;
			for (int m = 0; m < n; m++) {
				nwt = m * omega*t;
				r.x += a[m].x*cos(nwt) + b[m].x*sin(nwt);
				r.y += a[m].y*cos(nwt) + b[m].y*sin(nwt);
				r.z += a[m].z*cos(nwt) + b[m].z*sin(nwt);
				if (m == 0) r /= 2;
			}
			return r;
		};
		parallelogram_grid P(parallelogram(point(-3, -ERR_UPSILON, 0), point(7, 0), point(0, 2 * ERR_UPSILON))); W->add(P);
		GlassMan_std G;
		G.Pos = point(0, Pos_m*t, 0), G.top = vec3(0, 0, 1), G.dir = vec3(0, 1, 0);
		G.v_head = FourierEval(t, head_a, head_b, N), G.v_head_side = vec3(0, 1, 0), G.v_neck = FourierEval(t, neck_a, neck_b, N);
		G.v_chest = FourierEval(t, chest_a, chest_b, N), G.v_chest_side = vec3(0, 1, 0);
		G.v_waist = FourierEval(t, waist_a, waist_b, N), G.v_waist_side = vec3(0, 1, 0);
		G.v_upperarm_l = FourierEval(t + 0.6, upperarm_r_a, upperarm_r_b, N), G.v_forearm_l = FourierEval(t + 0.6, forearm_r_a, forearm_r_b, N);
		G.v_upperarm_r = FourierEval(t, upperarm_r_a, upperarm_r_b, N), G.v_forearm_r = FourierEval(t, forearm_r_a, forearm_r_b, N);
		G.v_thigh_l = FourierEval(t, thigh_l_a, thigh_l_b, N), G.v_shank_l = FourierEval(t, shank_l_a, shank_l_b, N);
		G.v_thigh_r = FourierEval(t, thigh_r_a, thigh_r_b, N), G.v_shank_r = FourierEval(t, shank_r_a, shank_r_b, N);
		G.v_foot_l = FourierEval(t, foot_l_a, foot_l_b, N), G.v_foot_l_side = vec3(0, 1, 0);
		G.v_foot_r = FourierEval(t, foot_r_a, foot_r_b, N), G.v_foot_r_side = vec3(0, -1, 0);
		G.auto_fit = true;
		G.push(*W);
		//ADD_AXIS(*W, 0.1);
	}, [](double t) -> point {	// camera
		return point(-10, 5, 10);
	}, [](double t) -> point {	// view point
		double Pos_m, Pos_b; WalkingMan::LinearRegression(WalkingMan::waist, 0, Pos_m, Pos_b);
		return point(0, Pos_m*t, 0.8);
	}, [](double t) -> double {		// solid angle
		return 0.05;
	}, -10, 10, 600, 400, 25, 2, "Animation\\CPTAT4_", 3);
}


void SetWalkingAttitudes(double t, GlassMan_std &G) {
	double Pos_m, Pos_b; WalkingMan::LinearRegression(WalkingMan::waist, 0, Pos_m, Pos_b);
	const unsigned N = 6;
	vec3 thigh_l_a[N], thigh_l_b[N]; WalkingMan::FourierSeries(WalkingMan::v_thigh_l, N, thigh_l_a, thigh_l_b);
	vec3 thigh_r_a[N], thigh_r_b[N]; WalkingMan::FourierSeries(WalkingMan::v_thigh_r, N, thigh_r_a, thigh_r_b);
	vec3 shank_l_a[N], shank_l_b[N]; WalkingMan::FourierSeries(WalkingMan::v_shank_l, N, shank_l_a, shank_l_b);
	vec3 shank_r_a[N], shank_r_b[N]; WalkingMan::FourierSeries(WalkingMan::v_shank_r, N, shank_r_a, shank_r_b);
	vec3 foot_l_a[N], foot_l_b[N]; WalkingMan::FourierSeries(WalkingMan::v_foot_l, N, foot_l_a, foot_l_b);
	vec3 foot_r_a[N], foot_r_b[N]; WalkingMan::FourierSeries(WalkingMan::v_foot_r, N, foot_r_a, foot_r_b);
	vec3 upperarm_r_a[N], upperarm_r_b[N]; WalkingMan::FourierSeries(WalkingMan::v_upperarm_r, N, upperarm_r_a, upperarm_r_b);
	vec3 forearm_r_a[N], forearm_r_b[N]; WalkingMan::FourierSeries(WalkingMan::v_forearm_r, N, forearm_r_a, forearm_r_b);
	vec3 chest_a[N], chest_b[N]; WalkingMan::FourierSeries(WalkingMan::v_chest, N, chest_a, chest_b);
	vec3 waist_a[N], waist_b[N]; WalkingMan::FourierSeries(WalkingMan::v_waist, N, waist_a, waist_b);
	vec3 head_a[N], head_b[N]; WalkingMan::FourierSeries(WalkingMan::v_head, N, head_a, head_b);
	vec3 neck_a[N], neck_b[N]; WalkingMan::FourierSeries(WalkingMan::v_neck, N, neck_a, neck_b);
	auto FourierEval = [](double t, vec3 *a, vec3 *b, int n) -> vec3 {
		const double omega = 2 * PI / 1.2; double nwt;
		vec3 r; r.x = r.y = r.z = 0;
		for (int m = 0; m < n; m++) {
			nwt = m * omega*t;
			r.x += a[m].x*cos(nwt) + b[m].x*sin(nwt);
			r.y += a[m].y*cos(nwt) + b[m].y*sin(nwt);
			r.z += a[m].z*cos(nwt) + b[m].z*sin(nwt);
			if (m == 0) r /= 2;
		}
		return r;
	};
	G.Pos = point(0, Pos_m*t, 0), G.top = vec3(0, 0, 1), G.dir = vec3(0, 1, 0);
	G.v_head = FourierEval(t, head_a, head_b, N), G.v_head_side = vec3(0, 1, 0), G.v_neck = FourierEval(t, neck_a, neck_b, N);
	G.v_chest = FourierEval(t, chest_a, chest_b, N), G.v_chest_side = vec3(0, 1, 0);
	G.v_waist = FourierEval(t, waist_a, waist_b, N), G.v_waist_side = vec3(0, 1, 0);
	G.v_upperarm_l = FourierEval(t + 0.6, upperarm_r_a, upperarm_r_b, N), G.v_forearm_l = FourierEval(t + 0.6, forearm_r_a, forearm_r_b, N);
	G.v_upperarm_r = FourierEval(t, upperarm_r_a, upperarm_r_b, N), G.v_forearm_r = FourierEval(t, forearm_r_a, forearm_r_b, N);
	G.v_thigh_l = FourierEval(t, thigh_l_a, thigh_l_b, N), G.v_shank_l = FourierEval(t, shank_l_a, shank_l_b, N);
	G.v_thigh_r = FourierEval(t, thigh_r_a, thigh_r_b, N), G.v_shank_r = FourierEval(t, shank_r_a, shank_r_b, N);
	G.v_foot_l = FourierEval(t, foot_l_a, foot_l_b, N), G.v_foot_l_side = vec3(0, 1, 0);
	G.v_foot_r = FourierEval(t, foot_r_a, foot_r_b, N), G.v_foot_r_side = vec3(0, -1, 0);
	G.auto_fit = true;
};
void CPT_Animation_R() {
	const unsigned width = 1920, height = 1080, sampling = 2, fps = 25;		// standard
	//const unsigned width = 768, height = 432, sampling = 1, fps = 10;
	const unsigned R = 0;

#pragma region R0
	if (R == -1 || R == 0) CPT_Animation([](double t, World* W) {	// World
		parallelogram_grid P(parallelogram(point(6, 0, 0), vec3(-ERR_UPSILON, 0), vec3(0, 3)), 1, 1, Gray, LightBlue); W->add(P);
		parallelogram_grid Ps(parallelogram(point(8, 0, 0), vec3(3, 0), vec3(0, 3)), 1, 1, Gray, LightBlue); W->add(Ps);
		World X = GoldenCoin_Constructor(point(0, 0, 0), 0);
		*((XSolid*)(X.Objs[0])) = CSG_Translation(CSG_Scale(*((XSolid*)(X.Objs[0])), 1.8), point(9.5, 1.5, 0.72));
		*((XSolid*)(X.Objs[1])) = CSG_Translation(CSG_Scale(*((XSolid*)(X.Objs[1])), 1.8), point(9.5, 1.5, 0.72));
		W->insert(X);
		GlassMan_std G; SetWalkingAttitudes(1.25 * (t + 1.4), G); G.Pos.x = G.Pos.y, G.Pos.y = 1.5; G.dir = vec3(1, 0, 0);
		if (t > 4.0) G = mix(G, ManStandingStraight_Constructor(vec3(7 - 0.5 * GlassMan_std::length_of_feet, G.Pos.y), G.dir), clamp(2 * (t - 4.3), 0., 1.));
		if (t > 7) G.Pos.z -= 4.9 * pow(t - 7, 2);
		G.push(*W);

		Figure2D TheGlass("M322.44766 130.79099l-5.296875 24.890625l-19.203125 0l-21.140625 99.484375l-32.343765 0l21.14064 -99.484375l-19.125015 0l5.296875 -24.890625l70.67189 0zM375.0911 164.52536l-19.265625 90.640625l-23.5625 0l8.09375 -38.0625l-7.0625 0l-8.09375 38.0625l-23.5625 0l19.265625 -90.640625l23.5625 0l-6.890625 32.421875l7.0625 0l6.890625 -32.421875l23.5625 0zm9.466858 0l39.296875 0l-3.859375 18.140625l-15.734375 0l-3.65625 17.1875l14.734375 0l-3.65625 17.25l-14.734375 0l-4.234375 19.921875l17.3125 0l-3.859375 18.140625l-40.875 0l19.265625 -90.640625zM525.80176 176.57224l-32.34372 0l2.40625 -11.28125q2.265625 -10.6875 1.90625 -13.375q-0.34375 -2.6875 -3.796875 -2.6875q-3.0 0 -4.5625 2.3125q-1.5625 2.296875 -3.578125 11.828125l-12.6875 59.6875q-1.78125 8.375 -1.28125 11.03125q0.515625 2.640625 3.75 2.640625q3.53125 0 5.421875 -2.984375q1.90625 -3.0 3.75 -11.6875l3.140625 -14.75l-6.53125 0l4.015625 -18.890625l37.87497 0l-14.1874695 66.75l-20.359375 0l-1.09375 -8.90625q-4.53125 5.75 -10.171875 8.640625q-5.640625 2.875 -12.46875 2.875q-8.140625 0 -14.421875 -3.953125q-6.265625 -3.953125 -8.71875 -9.796875q-2.4375 -5.84375 -2.0 -12.25q0.453125 -6.421875 3.171875 -19.25l7.859375 -36.9375q3.78125 -17.828125 7.421875 -25.890625q3.640625 -8.078125 14.15625 -14.796875q10.53125 -6.71875 24.96875 -6.71875q14.21875 0 22.359344 5.84375q8.140625 5.828125 9.265625 13.859375q1.125 8.03125 -2.125 23.3125l-1.140625 5.375zM562.70386 164.52536l-15.40625 72.5l14.34375 0l-3.859375 18.140625l-37.90625 0l19.265625 -90.640625l23.5625 0zm60.752563 0l-5.78125 90.640625l-24.125 0l2.296875 -16.296875l-8.453125 0l-4.875 16.296875l-24.40625 0l31.25 -90.640625l34.09375 0zm-24.875 58.28125q1.484375 -15.390625 4.484375 -38.015625q-9.09375 25.984375 -12.546875 38.015625l8.0625 0zm86.50513 -30.84375l-21.890625 0l1.421875 -6.71875q1.0 -4.703125 0.4375 -5.984375q-0.5625 -1.296875 -2.53125 -1.296875q-2.125 0 -3.59375 1.734375q-1.453125 1.734375 -2.203125 5.265625q-0.96875 4.53125 -0.21875 6.828125q0.6875 2.296875 5.484375 5.546875q13.75 9.34375 16.5625 15.328125q2.8125 6.0 -0.03125 19.328125q-2.046875 9.671875 -5.296875 14.265625q-3.234375 4.59375 -10.390625 7.703125q-7.15625 3.109375 -15.78125 3.109375q-9.453125 0 -15.390625 -3.578125q-5.921875 -3.59375 -6.828125 -9.125q-0.890625 -5.546875 1.28125 -15.734375l1.265625 -5.9375l21.890625 0l-2.34375 11.03125q-1.09375 5.09375 -0.484375 6.546875q0.625 1.453125 2.984375 1.453125q2.34375 0 3.875 -1.84375q1.546875 -1.84375 2.3125 -5.484375q1.71875 -8.015625 0.046875 -10.46875q-1.71875 -2.46875 -9.265625 -8.234375q-7.5625 -5.828125 -9.859375 -8.453125q-2.28125 -2.625 -3.171875 -7.265625q-0.890625 -4.65625 0.640625 -11.875q2.21875 -10.421875 5.890625 -15.234375q3.6875 -4.8125 10.203125 -7.53125q6.515625 -2.71875 14.90625 -2.71875q9.1875 0 15.015625 2.96875q5.828125 2.96875 6.96875 7.484375q1.15625 4.5 -1.140625 15.296875l-0.765625 3.59375zm59.287537 0l-21.890625 0l1.421875 -6.71875q1.0 -4.703125 0.4375 -5.984375q-0.5625 -1.296875 -2.53125 -1.296875q-2.125 0 -3.59375 1.734375q-1.453125 1.734375 -2.203125 5.265625q-0.96875 4.53125 -0.21875 6.828125q0.6875 2.296875 5.484375 5.546875q13.75 9.34375 16.5625 15.328125q2.8125 6.0 -0.03125 19.328125q-2.046875 9.671875 -5.296875 14.265625q-3.234375 4.59375 -10.390625 7.703125q-7.15625 3.109375 -15.78125 3.109375q-9.453125 0 -15.390625 -3.578125q-5.921875 -3.59375 -6.828125 -9.125q-0.890625 -5.546875 1.28125 -15.734375l1.265625 -5.9375l21.890625 0l-2.34375 11.03125q-1.09375 5.09375 -0.484375 6.546875q0.625 1.453125 2.984375 1.453125q2.34375 0 3.875 -1.84375q1.546875 -1.84375 2.3125 -5.484375q1.71875 -8.015625 0.046875 -10.46875q-1.71875 -2.46875 -9.265625 -8.234375q-7.5625 -5.828125 -9.859375 -8.453125q-2.28125 -2.625 -3.171875 -7.265625q-0.890625 -4.65625 0.640625 -11.875q2.21875 -10.421875 5.890625 -15.234375q3.6875 -4.8125 10.203125 -7.53125q6.515625 -2.71875 14.90625 -2.71875q9.1875 0 15.015625 2.96875q5.828125 2.96875 6.96875 7.484375q1.15625 4.5 -1.140625 15.296875l-0.765625 3.59375z");
		Figure2D TheGlass_T("M322.44766 130.79099l-5.296875 24.890625l-19.203125 0l-21.140625 99.484375l-32.343765 0l21.14064 -99.484375l-19.125015 0l5.296875 -24.890625l70.67189 0z");
		Figure2D TheGlass_H("M375.0911,164.52536L355.825475,255.165985L332.262975,255.165985L340.356725,217.103485L333.294225,217.103485L325.200475,255.165985L301.637975,255.165985L320.9036,164.52536L344.4661,164.52536L337.575475,196.947235L344.637975,196.947235L351.5286,164.52536L375.0911,164.52536Z");
		Figure2D TheGlass_E("M384.557958,164.52536L423.854833,164.52536L419.995458,182.665985L404.261083,182.665985L400.604833,199.853485L415.339208,199.853485L411.682958,217.103485L396.948583,217.103485L392.714208,237.02536L410.026708,237.02536L406.167333,255.165985L365.292333,255.165985L384.557958,164.52536Z");
		Figure2D TheGlass_G("M525.80176 176.57224l-32.34372 0l2.40625 -11.28125q2.265625 -10.6875 1.90625 -13.375q-0.34375 -2.6875 -3.796875 -2.6875q-3.0 0 -4.5625 2.3125q-1.5625 2.296875 -3.578125 11.828125l-12.6875 59.6875q-1.78125 8.375 -1.28125 11.03125q0.515625 2.640625 3.75 2.640625q3.53125 0 5.421875 -2.984375q1.90625 -3.0 3.75 -11.6875l3.140625 -14.75l-6.53125 0l4.015625 -18.890625l37.87497 0l-14.1874695 66.75l-20.359375 0l-1.09375 -8.90625q-4.53125 5.75 -10.171875 8.640625q-5.640625 2.875 -12.46875 2.875q-8.140625 0 -14.421875 -3.953125q-6.265625 -3.953125 -8.71875 -9.796875q-2.4375 -5.84375 -2.0 -12.25q0.453125 -6.421875 3.171875 -19.25l7.859375 -36.9375q3.78125 -17.828125 7.421875 -25.890625q3.640625 -8.078125 14.15625 -14.796875q10.53125 -6.71875 24.96875 -6.71875q14.21875 0 22.359344 5.84375q8.140625 5.828125 9.265625 13.859375q1.125 8.03125 -2.125 23.3125l-1.140625 5.375z");
		Figure2D TheGlass_L("M562.70386,164.52536L547.29761,237.02536L561.64136,237.02536L557.781985,255.165985L519.875735,255.165985L539.14136,164.52536L562.70386,164.52536Z");
		Figure2D TheGlass_A("M623.456423,164.52536L617.675173,255.165985L593.550173,255.165985L595.847048,238.86911L587.393923,238.86911L582.518923,255.165985L558.112673,255.165985L589.362673,164.52536L623.456423,164.52536ZM598.581423,222.80661Q600.065798,207.415985,603.065798,184.790985Q593.972048,210.77536,590.518923,222.80661L598.581423,222.80661Z");
		Figure2D TheGlass_S1("M685.086553,191.96286L663.195928,191.96286L664.617803,185.24411Q665.617803,180.540985,665.055303,179.259735Q664.492803,177.96286,662.524053,177.96286Q660.399053,177.96286,658.930303,179.697235Q657.477178,181.43161,656.727178,184.96286Q655.758428,189.49411,656.508428,191.790985Q657.195928,194.08786,661.992803,197.33786Q675.742803,206.68161,678.555303,212.665985Q681.367803,218.665985,678.524053,231.99411Q676.477178,241.665985,673.227178,246.259735Q669.992803,250.853485,662.836553,253.96286Q655.680303,257.072235,647.055303,257.072235Q637.602178,257.072235,631.664678,253.49410999999998Q625.742803,249.90035999999998,624.836553,244.36910999999998Q623.945928,238.82223499999998,626.117803,228.63473499999998L627.383428,222.69723499999998L649.274053,222.69723499999998L646.930303,233.72848499999998Q645.836553,238.82223499999998,646.445928,240.27535999999998Q647.070928,241.72848499999998,649.430303,241.72848499999998Q651.774053,241.72848499999998,653.305303,239.88473499999998Q654.852178,238.04098499999998,655.617803,234.40035999999998Q657.336553,226.38473499999998,655.664678,223.93160999999998Q653.945928,221.46285999999998,646.399053,215.69723499999998Q638.836553,209.86910999999998,636.539678,207.24410999999998Q634.258428,204.61910999999998,633.367803,199.97848499999998Q632.477178,195.32223499999998,634.008428,188.10348499999998Q636.227178,177.68160999999998,639.899053,172.86910999999998Q643.586553,168.05660999999998,650.102178,165.33785999999998Q656.617803,162.61910999999998,665.008428,162.61910999999998Q674.195928,162.61910999999998,680.024053,165.58785999999998Q685.852178,168.55660999999998,686.992803,173.07223499999998Q688.149053,177.57223499999998,685.852178,188.36910999999998L685.086553,191.96285999999998Z");
		Figure2D TheGlass_S2("M744.37409,191.96286L722.483465,191.96286L723.90534,185.24411Q724.90534,180.540985,724.34284,179.259735Q723.78034,177.96286,721.81159,177.96286Q719.68659,177.96286,718.21784,179.697235Q716.764715,181.43161,716.014715,184.96286Q715.045965,189.49411,715.795965,191.790985Q716.483465,194.08786,721.28034,197.33786Q735.03034,206.68161,737.84284,212.665985Q740.65534,218.665985,737.81159,231.99411Q735.764715,241.665985,732.514715,246.259735Q729.28034,250.853485,722.12409,253.96286Q714.96784,257.072235,706.34284,257.072235Q696.889715,257.072235,690.952215,253.49410999999998Q685.03034,249.90035999999998,684.12409,244.36910999999998Q683.233465,238.82223499999998,685.40534,228.63473499999998L686.670965,222.69723499999998L708.56159,222.69723499999998L706.21784,233.72848499999998Q705.12409,238.82223499999998,705.733465,240.27535999999998Q706.358465,241.72848499999998,708.71784,241.72848499999998Q711.06159,241.72848499999998,712.59284,239.88473499999998Q714.139715,238.04098499999998,714.90534,234.40035999999998Q716.62409,226.38473499999998,714.952215,223.93160999999998Q713.233465,221.46285999999998,705.68659,215.69723499999998Q698.12409,209.86910999999998,695.827215,207.24410999999998Q693.545965,204.61910999999998,692.65534,199.97848499999998Q691.764715,195.32223499999998,693.295965,188.10348499999998Q695.514715,177.68160999999998,699.18659,172.86910999999998Q702.87409,168.05660999999998,709.389715,165.33785999999998Q715.90534,162.61910999999998,724.295965,162.61910999999998Q733.483465,162.61910999999998,739.31159,165.58785999999998Q745.139715,168.55660999999998,746.28034,173.07223499999998Q747.43659,177.57223499999998,745.139715,188.36910999999998L744.37409,191.96285999999998Z");

		Figure2D ByNames1("m345.71906 308.19638q0 0.765625 -0.21875 1.453125q-0.203125 0.6875 -0.625 1.25q-0.40625 0.5625 -0.984375 1.0q-0.578125 0.4375 -1.296875 0.703125q0.578125 0.09375 1.0625 0.359375q0.484375 0.265625 0.828125 0.671875q0.34375 0.390625 0.53125 0.953125q0.1875 0.546875 0.1875 1.234375q0 0.625 -0.171875 1.296875q-0.171875 0.671875 -0.53125 1.3125q-0.359375 0.640625 -0.953125 1.21875q-0.578125 0.578125 -1.40625 1.015625q-0.828125 0.4375 -1.921875 0.6875q-1.078125 0.25 -2.5 0.25l-4.671875 0q-0.53125 0 -0.78125 -0.25q-0.234375 -0.265625 -0.125 -0.828125l2.953125 -14.6875q0.09375 -0.5625 0.453125 -0.8125q0.359375 -0.265625 0.828125 -0.265625l4.046875 0q1.34375 0 2.3125 0.21875q0.984375 0.203125 1.640625 0.625q0.671875 0.421875 1.0 1.0625q0.34375 0.640625 0.34375 1.53125zm-3.53125 0.640625q0 -0.375 -0.125 -0.640625q-0.125 -0.265625 -0.390625 -0.453125q-0.265625 -0.203125 -0.65625 -0.3125q-0.390625 -0.109375 -1.015625 -0.109375l-1.859375 0l-0.875 4.375l1.921875 0q0.84375 0 1.40625 -0.265625q0.5625 -0.265625 0.921875 -0.6875q0.359375 -0.421875 0.515625 -0.921875q0.15625 -0.515625 0.15625 -0.984375zm-0.5625 7.15625q0 -0.828125 -0.625 -1.3125q-0.609375 -0.484375 -1.96875 -0.484375l-2.28125 0l-0.953125 4.796875l2.421875 0q0.921875 0 1.546875 -0.265625q0.640625 -0.265625 1.046875 -0.671875q0.421875 -0.421875 0.609375 -0.953125q0.203125 -0.546875 0.203125 -1.109375zm17.613953 -6.484375q0 0.171875 -0.046875 0.390625q-0.03125 0.203125 -0.109375 0.421875q-0.4375 1.1875 -0.953125 2.46875q-0.515625 1.28125 -1.1875 2.640625q-0.65625 1.34375 -1.5 2.796875q-0.84375 1.4375 -1.9375 2.984375l-3.171875 4.5q-0.125 0.171875 -0.296875 0.28125q-0.171875 0.109375 -0.4375 0.171875q-0.265625 0.078125 -0.65625 0.109375q-0.375 0.046875 -0.859375 0.046875q-0.5625 0 -0.90625 -0.046875q-0.34375 -0.03125 -0.5 -0.140625q-0.15625 -0.109375 -0.15625 -0.265625q0.015625 -0.15625 0.1875 -0.375l2.953125 -3.890625l-1.90625 -11.171875q-0.03125 -0.171875 -0.0625 -0.34375q-0.015625 -0.1875 -0.015625 -0.3125q0 -0.234375 0.078125 -0.375q0.09375 -0.15625 0.28125 -0.234375q0.203125 -0.078125 0.546875 -0.109375q0.34375 -0.03125 0.828125 -0.03125q0.53125 0 0.828125 0.03125q0.3125 0.015625 0.46875 0.078125q0.15625 0.0625 0.203125 0.203125q0.0625 0.125 0.09375 0.3125l1.203125 8.140625l0.015625 0q0.65625 -0.984375 1.1875 -1.984375q0.53125 -1.0 0.96875 -2.03125q0.453125 -1.03125 0.8125 -2.09375q0.375 -1.0625 0.6875 -2.140625q0.03125 -0.140625 0.125 -0.234375q0.09375 -0.09375 0.296875 -0.15625q0.21875 -0.0625 0.5625 -0.09375q0.359375 -0.03125 0.921875 -0.03125q0.828125 0 1.140625 0.109375q0.3125 0.109375 0.3125 0.375zm19.41983 11.640625q-0.015625 0.125 -0.125 0.234375q-0.109375 0.109375 -0.3125 0.171875q-0.203125 0.0625 -0.546875 0.09375q-0.34375 0.03125 -0.828125 0.03125q-0.5 0 -0.828125 -0.03125q-0.328125 -0.03125 -0.5 -0.09375q-0.171875 -0.0625 -0.234375 -0.171875q-0.0625 -0.109375 -0.03125 -0.234375l1.359375 -6.828125l-6.296875 0l-1.359375 6.828125q-0.03125 0.125 -0.140625 0.234375q-0.09375 0.109375 -0.3125 0.171875q-0.21875 0.0625 -0.546875 0.09375q-0.328125 0.03125 -0.828125 0.03125q-0.5 0 -0.828125 -0.03125q-0.3125 -0.03125 -0.5 -0.09375q-0.171875 -0.0625 -0.25 -0.171875q-0.0625 -0.109375 -0.015625 -0.234375l3.1875 -15.9375q0.015625 -0.125 0.109375 -0.21875q0.109375 -0.109375 0.3125 -0.171875q0.21875 -0.078125 0.546875 -0.109375q0.34375 -0.03125 0.84375 -0.03125q0.484375 0 0.8125 0.03125q0.328125 0.03125 0.5 0.109375q0.1875 0.0625 0.25 0.171875q0.0625 0.09375 0.03125 0.21875l-1.234375 6.203125l6.3125 0l1.21875 -6.203125q0.03125 -0.125 0.125 -0.21875q0.109375 -0.109375 0.3125 -0.171875q0.21875 -0.078125 0.546875 -0.109375q0.34375 -0.03125 0.84375 -0.03125q0.5 0 0.8125 0.03125q0.328125 0.03125 0.515625 0.109375q0.1875 0.0625 0.234375 0.171875q0.0625 0.09375 0.03125 0.21875l-3.1875 15.9375zm14.26947 0.03125q-0.046875 0.265625 -0.40625 0.390625q-0.359375 0.109375 -1.15625 0.109375q-0.40625 0 -0.65625 -0.03125q-0.25 -0.015625 -0.40625 -0.078125q-0.15625 -0.0625 -0.203125 -0.15625q-0.046875 -0.09375 -0.015625 -0.234375l0.359375 -1.875q-0.203125 0.421875 -0.640625 0.890625q-0.421875 0.453125 -1.0 0.84375q-0.578125 0.375 -1.296875 0.609375q-0.703125 0.25 -1.46875 0.25q-1.03125 0 -1.734375 -0.359375q-0.6875 -0.359375 -1.109375 -0.984375q-0.40625 -0.625 -0.578125 -1.40625q-0.171875 -0.796875 -0.171875 -1.671875q0 -0.828125 0.15625 -1.796875q0.15625 -0.96875 0.46875 -1.921875q0.328125 -0.96875 0.828125 -1.859375q0.515625 -0.890625 1.21875 -1.578125q0.71875 -0.6875 1.640625 -1.09375q0.921875 -0.421875 2.078125 -0.421875q1.078125 0 1.90625 0.453125q0.828125 0.453125 1.453125 1.3125l0.21875 -1.03125q0.046875 -0.28125 0.40625 -0.390625q0.359375 -0.125 1.15625 -0.125q0.390625 0 0.640625 0.03125q0.265625 0.03125 0.40625 0.09375q0.15625 0.0625 0.203125 0.15625q0.0625 0.09375 0.03125 0.234375l-2.328125 11.640625zm-1.609375 -8.15625q-0.453125 -0.71875 -1.015625 -1.078125q-0.546875 -0.359375 -1.28125 -0.359375q-0.546875 0 -1.0 0.265625q-0.453125 0.265625 -0.8125 0.703125q-0.34375 0.4375 -0.609375 1.015625q-0.265625 0.5625 -0.453125 1.1875q-0.171875 0.609375 -0.265625 1.25q-0.078125 0.625 -0.078125 1.15625q0 0.390625 0.0625 0.75q0.0625 0.34375 0.21875 0.625q0.15625 0.265625 0.40625 0.4375q0.265625 0.15625 0.671875 0.15625q0.578125 0 1.171875 -0.34375q0.59375 -0.34375 1.109375 -0.9375q0.53125 -0.59375 0.921875 -1.375q0.40625 -0.796875 0.578125 -1.6875l0.375 -1.765625zm14.538757 -3.65625q0 0.09375 -0.03125 0.328125q-0.03125 0.234375 -0.078125 0.546875q-0.046875 0.3125 -0.125 0.640625q-0.078125 0.328125 -0.1875 0.625q-0.09375 0.28125 -0.21875 0.46875q-0.125 0.171875 -0.265625 0.171875q-0.109375 0 -0.234375 -0.046875q-0.109375 -0.046875 -0.265625 -0.078125q-0.140625 -0.046875 -0.3125 -0.09375q-0.15625 -0.046875 -0.375 -0.046875q-0.453125 0 -0.953125 0.34375q-0.5 0.328125 -0.953125 0.890625q-0.453125 0.546875 -0.8125 1.28125q-0.359375 0.734375 -0.515625 1.546875l-1.046875 5.21875q-0.015625 0.125 -0.125 0.234375q-0.09375 0.09375 -0.296875 0.15625q-0.203125 0.0625 -0.53125 0.09375q-0.3125 0.03125 -0.78125 0.03125q-0.484375 0 -0.796875 -0.03125q-0.296875 -0.03125 -0.484375 -0.09375q-0.171875 -0.0625 -0.21875 -0.15625q-0.046875 -0.109375 -0.03125 -0.234375l2.328125 -11.625q0.015625 -0.140625 0.109375 -0.234375q0.09375 -0.09375 0.265625 -0.15625q0.171875 -0.0625 0.453125 -0.09375q0.28125 -0.03125 0.671875 -0.03125q0.40625 0 0.65625 0.03125q0.265625 0.03125 0.40625 0.09375q0.140625 0.0625 0.1875 0.15625q0.046875 0.09375 0.015625 0.234375l-0.359375 1.84375q0.3125 -0.5625 0.71875 -1.03125q0.40625 -0.46875 0.859375 -0.8125q0.46875 -0.34375 0.96875 -0.53125q0.5 -0.203125 0.96875 -0.203125q0.203125 0 0.40625 0.03125q0.203125 0.015625 0.375 0.0625q0.171875 0.03125 0.3125 0.09375q0.140625 0.0625 0.21875 0.125q0.078125 0.078125 0.078125 0.25zm9.384369 0q0 0.09375 -0.03125 0.328125q-0.03125 0.234375 -0.078125 0.546875q-0.046875 0.3125 -0.125 0.640625q-0.078125 0.328125 -0.1875 0.625q-0.09375 0.28125 -0.21875 0.46875q-0.125 0.171875 -0.265625 0.171875q-0.109375 0 -0.234375 -0.046875q-0.109375 -0.046875 -0.265625 -0.078125q-0.140625 -0.046875 -0.3125 -0.09375q-0.15625 -0.046875 -0.375 -0.046875q-0.453125 0 -0.953125 0.34375q-0.5 0.328125 -0.953125 0.890625q-0.453125 0.546875 -0.8125 1.28125q-0.359375 0.734375 -0.515625 1.546875l-1.046875 5.21875q-0.015625 0.125 -0.125 0.234375q-0.09375 0.09375 -0.296875 0.15625q-0.203125 0.0625 -0.53125 0.09375q-0.3125 0.03125 -0.78125 0.03125q-0.484375 0 -0.796875 -0.03125q-0.296875 -0.03125 -0.484375 -0.09375q-0.171875 -0.0625 -0.21875 -0.15625q-0.046875 -0.109375 -0.03125 -0.234375l2.328125 -11.625q0.015625 -0.140625 0.109375 -0.234375q0.09375 -0.09375 0.265625 -0.15625q0.171875 -0.0625 0.453125 -0.09375q0.28125 -0.03125 0.671875 -0.03125q0.40625 0 0.65625 0.03125q0.265625 0.03125 0.40625 0.09375q0.140625 0.0625 0.1875 0.15625q0.046875 0.09375 0.015625 0.234375l-0.359375 1.84375q0.3125 -0.5625 0.71875 -1.03125q0.40625 -0.46875 0.859375 -0.8125q0.46875 -0.34375 0.96875 -0.53125q0.5 -0.203125 0.96875 -0.203125q0.203125 0 0.40625 0.03125q0.203125 0.015625 0.375 0.0625q0.171875 0.03125 0.3125 0.09375q0.140625 0.0625 0.21875 0.125q0.078125 0.078125 0.078125 0.25zm12.212463 0.140625q0 0.171875 -0.046875 0.390625q-0.03125 0.203125 -0.109375 0.421875q-0.4375 1.1875 -0.953125 2.46875q-0.515625 1.28125 -1.1875 2.640625q-0.65625 1.34375 -1.5 2.796875q-0.84375 1.4375 -1.9375 2.984375l-3.171875 4.5q-0.125 0.171875 -0.296875 0.28125q-0.171875 0.109375 -0.4375 0.171875q-0.265625 0.078125 -0.65625 0.109375q-0.375 0.046875 -0.859375 0.046875q-0.5625 0 -0.90625 -0.046875q-0.34375 -0.03125 -0.5 -0.140625q-0.15625 -0.109375 -0.15625 -0.265625q0.015625 -0.15625 0.1875 -0.375l2.953125 -3.890625l-1.90625 -11.171875q-0.03125 -0.171875 -0.0625 -0.34375q-0.015625 -0.1875 -0.015625 -0.3125q0 -0.234375 0.078125 -0.375q0.09375 -0.15625 0.28125 -0.234375q0.203125 -0.078125 0.546875 -0.109375q0.34375 -0.03125 0.828125 -0.03125q0.53125 0 0.828125 0.03125q0.3125 0.015625 0.46875 0.078125q0.15625 0.0625 0.203125 0.203125q0.0625 0.125 0.09375 0.3125l1.203125 8.140625l0.015625 0q0.65625 -0.984375 1.1875 -1.984375q0.53125 -1.0 0.96875 -2.03125q0.453125 -1.03125 0.8125 -2.09375q0.375 -1.0625 0.6875 -2.140625q0.03125 -0.140625 0.125 -0.234375q0.09375 -0.09375 0.296875 -0.15625q0.21875 -0.0625 0.5625 -0.09375q0.359375 -0.03125 0.921875 -0.03125q0.828125 0 1.140625 0.109375q0.3125 0.109375 0.3125 0.375zm20.48233 -3.125q0 0.375 -0.125 0.984375q-0.125 0.609375 -0.3125 0.984375q-0.1875 0.359375 -0.421875 0.359375q-0.171875 0 -0.40625 -0.203125q-0.234375 -0.21875 -0.609375 -0.46875q-0.375 -0.25 -0.96875 -0.453125q-0.578125 -0.21875 -1.453125 -0.21875q-0.96875 0 -1.765625 0.390625q-0.78125 0.375 -1.40625 1.015625q-0.625 0.625 -1.09375 1.453125q-0.46875 0.8125 -0.78125 1.703125q-0.296875 0.875 -0.453125 1.75q-0.15625 0.875 -0.15625 1.640625q0 0.90625 0.234375 1.59375q0.234375 0.6875 0.65625 1.15625q0.4375 0.453125 1.046875 0.6875q0.625 0.21875 1.40625 0.21875q0.9375 0 1.59375 -0.1875q0.65625 -0.203125 1.109375 -0.4375q0.46875 -0.234375 0.765625 -0.421875q0.3125 -0.203125 0.5625 -0.203125q0.171875 0 0.21875 0.125q0.0625 0.125 0.0625 0.34375q0 0.09375 -0.03125 0.28125q-0.015625 0.171875 -0.0625 0.40625q-0.03125 0.21875 -0.09375 0.484375q-0.046875 0.265625 -0.109375 0.5q-0.0625 0.21875 -0.15625 0.421875q-0.078125 0.1875 -0.21875 0.328125q-0.140625 0.140625 -0.5625 0.375q-0.421875 0.21875 -1.015625 0.421875q-0.59375 0.1875 -1.359375 0.3125q-0.75 0.140625 -1.609375 0.140625q-1.421875 0 -2.5625 -0.375q-1.125 -0.375 -1.921875 -1.140625q-0.796875 -0.765625 -1.234375 -1.921875q-0.421875 -1.15625 -0.421875 -2.71875q0 -1.203125 0.25 -2.5q0.25 -1.3125 0.75 -2.5625q0.5 -1.25 1.265625 -2.375q0.765625 -1.125 1.8125 -1.953125q1.046875 -0.84375 2.375 -1.328125q1.34375 -0.5 2.96875 -0.5q0.90625 0 1.71875 0.1875q0.828125 0.1875 1.421875 0.5q0.609375 0.3125 0.84375 0.5625q0.25 0.25 0.25 0.640625zm10.275848 14.78125q-0.015625 0.125 -0.125 0.234375q-0.09375 0.09375 -0.296875 0.15625q-0.203125 0.0625 -0.53125 0.09375q-0.3125 0.03125 -0.796875 0.03125q-0.46875 0 -0.78125 -0.03125q-0.296875 -0.03125 -0.46875 -0.09375q-0.171875 -0.0625 -0.234375 -0.15625q-0.0625 -0.109375 -0.03125 -0.234375l1.390625 -6.96875q0.0625 -0.265625 0.09375 -0.609375q0.046875 -0.359375 0.046875 -0.625q0 -0.28125 -0.0625 -0.515625q-0.0625 -0.25 -0.203125 -0.421875q-0.125 -0.1875 -0.328125 -0.28125q-0.203125 -0.09375 -0.5 -0.09375q-0.53125 0 -1.109375 0.34375q-0.578125 0.328125 -1.09375 0.921875q-0.515625 0.578125 -0.921875 1.359375q-0.390625 0.78125 -0.5625 1.65625l-1.0625 5.234375q-0.015625 0.125 -0.125 0.234375q-0.09375 0.09375 -0.296875 0.15625q-0.203125 0.0625 -0.53125 0.09375q-0.3125 0.03125 -0.78125 0.03125q-0.484375 0 -0.796875 -0.03125q-0.296875 -0.03125 -0.484375 -0.09375q-0.171875 -0.0625 -0.21875 -0.15625q-0.046875 -0.109375 -0.03125 -0.234375l3.421875 -17.125q0.03125 -0.125 0.125 -0.21875q0.09375 -0.109375 0.296875 -0.1875q0.203125 -0.078125 0.515625 -0.109375q0.328125 -0.046875 0.8125 -0.046875q0.484375 0 0.78125 0.046875q0.3125 0.03125 0.484375 0.109375q0.171875 0.078125 0.21875 0.1875q0.046875 0.09375 0.03125 0.21875l-0.921875 4.71875q-0.078125 0.28125 -0.15625 0.609375q-0.078125 0.328125 -0.171875 0.65625q-0.09375 0.3125 -0.203125 0.609375q-0.09375 0.296875 -0.171875 0.53125q0.234375 -0.421875 0.65625 -0.828125q0.4375 -0.421875 0.984375 -0.75q0.5625 -0.34375 1.203125 -0.5625q0.65625 -0.21875 1.328125 -0.21875q0.890625 0 1.5 0.25q0.609375 0.25 1.0 0.6875q0.390625 0.4375 0.5625 1.0625q0.171875 0.609375 0.171875 1.328125q0 0.46875 -0.046875 0.9375q-0.046875 0.46875 -0.140625 0.921875l-1.4375 7.171875zm15.307007 -9.0625q0 0.9375 -0.421875 1.6875q-0.421875 0.75 -1.3125 1.296875q-0.890625 0.53125 -2.265625 0.8125q-1.375 0.28125 -3.265625 0.28125l-1.203125 0q-0.046875 0.3125 -0.078125 0.609375q-0.03125 0.28125 -0.03125 0.53125q0 1.046875 0.546875 1.59375q0.546875 0.53125 1.765625 0.53125q0.859375 0 1.53125 -0.125q0.6875 -0.125 1.203125 -0.265625q0.515625 -0.15625 0.84375 -0.28125q0.34375 -0.125 0.5 -0.125q0.15625 0 0.21875 0.09375q0.078125 0.09375 0.078125 0.28125q0 0.21875 -0.046875 0.515625q-0.03125 0.28125 -0.109375 0.5625q-0.078125 0.28125 -0.1875 0.53125q-0.109375 0.234375 -0.25 0.375q-0.15625 0.15625 -0.578125 0.3125q-0.40625 0.15625 -1.0 0.28125q-0.578125 0.125 -1.3125 0.203125q-0.71875 0.09375 -1.484375 0.09375q-1.25 0 -2.1875 -0.265625q-0.921875 -0.25 -1.546875 -0.78125q-0.625 -0.53125 -0.9375 -1.359375q-0.3125 -0.828125 -0.3125 -1.96875q0 -0.875 0.15625 -1.859375q0.171875 -1.0 0.53125 -1.96875q0.375 -0.96875 0.9375 -1.84375q0.578125 -0.890625 1.390625 -1.5625q0.8125 -0.671875 1.875 -1.078125q1.078125 -0.40625 2.4375 -0.40625q1.1875 0 2.03125 0.28125q0.859375 0.28125 1.40625 0.75q0.546875 0.453125 0.8125 1.046875q0.265625 0.59375 0.265625 1.21875zm-3.21875 0.25q0 -0.53125 -0.40625 -0.875q-0.40625 -0.34375 -1.1875 -0.34375q-0.65625 0 -1.171875 0.234375q-0.5 0.234375 -0.890625 0.640625q-0.390625 0.390625 -0.671875 0.9375q-0.28125 0.53125 -0.453125 1.140625l1.09375 0q1.03125 0 1.734375 -0.125q0.703125 -0.140625 1.125 -0.375q0.4375 -0.234375 0.625 -0.546875q0.203125 -0.328125 0.203125 -0.6875zm16.703217 -0.21875q0 0.46875 -0.0625 0.921875q-0.046875 0.453125 -0.140625 0.9375l-1.453125 7.171875q-0.015625 0.125 -0.125 0.234375q-0.09375 0.09375 -0.296875 0.15625q-0.203125 0.0625 -0.515625 0.09375q-0.3125 0.03125 -0.796875 0.03125q-0.484375 0 -0.796875 -0.03125q-0.296875 -0.03125 -0.46875 -0.09375q-0.171875 -0.0625 -0.234375 -0.15625q-0.046875 -0.109375 -0.03125 -0.234375l1.40625 -6.96875q0.0625 -0.3125 0.09375 -0.640625q0.046875 -0.328125 0.046875 -0.5625q0 -0.59375 -0.25 -0.96875q-0.25 -0.375 -0.84375 -0.375q-0.53125 0 -1.109375 0.34375q-0.578125 0.328125 -1.09375 0.90625q-0.515625 0.578125 -0.921875 1.359375q-0.390625 0.765625 -0.5625 1.640625l-1.0625 5.265625q-0.015625 0.125 -0.125 0.234375q-0.09375 0.09375 -0.296875 0.15625q-0.203125 0.0625 -0.53125 0.09375q-0.3125 0.03125 -0.78125 0.03125q-0.484375 0 -0.796875 -0.03125q-0.296875 -0.03125 -0.484375 -0.09375q-0.171875 -0.0625 -0.21875 -0.15625q-0.046875 -0.109375 -0.03125 -0.234375l2.328125 -11.625q0.03125 -0.140625 0.109375 -0.234375q0.09375 -0.09375 0.265625 -0.15625q0.171875 -0.0625 0.453125 -0.09375q0.28125 -0.03125 0.671875 -0.03125q0.40625 0 0.65625 0.03125q0.265625 0.03125 0.40625 0.09375q0.140625 0.0625 0.1875 0.15625q0.046875 0.09375 0.015625 0.234375l-0.375 1.90625q0.265625 -0.453125 0.71875 -0.921875q0.46875 -0.484375 1.046875 -0.859375q0.59375 -0.375 1.296875 -0.609375q0.703125 -0.25 1.453125 -0.25q0.90625 0 1.515625 0.265625q0.625 0.265625 1.015625 0.734375q0.390625 0.453125 0.546875 1.0625q0.171875 0.59375 0.171875 1.265625zm5.853882 6.59375q0 0.96875 -0.390625 1.921875q-0.390625 0.9375 -1.15625 1.765625l-2.234375 2.5q-0.234375 0.265625 -0.609375 0.375q-0.375 0.109375 -1.015625 0.109375q-0.375 0 -0.59375 -0.046875q-0.234375 -0.03125 -0.328125 -0.109375q-0.09375 -0.0625 -0.0625 -0.203125q0.015625 -0.125 0.125 -0.3125l2.28125 -3.59375l0.390625 -1.921875q0.09375 -0.453125 0.234375 -0.71875q0.15625 -0.28125 0.375 -0.421875q0.234375 -0.140625 0.5625 -0.1875q0.34375 -0.046875 0.84375 -0.046875q0.46875 0 0.78125 0.046875q0.3125 0.046875 0.484375 0.15625q0.171875 0.109375 0.234375 0.28125q0.078125 0.171875 0.078125 0.40625z");
		Figure2D ByNames2("m23.007965 -12.34375q0 0.375 -0.125 0.984375q-0.125 0.609375 -0.3125 0.984375q-0.1875 0.359375 -0.421875 0.359375q-0.171875 0 -0.40625 -0.203125q-0.234375 -0.21875 -0.609375 -0.46875q-0.375 -0.25 -0.96875 -0.453125q-0.578125 -0.21875 -1.4530945 -0.21875q-0.96875 0 -1.765625 0.390625q-0.78125 0.375 -1.40625 1.015625q-0.625 0.625 -1.09375 1.453125q-0.46875 0.8125 -0.78125 1.703125q-0.296875 0.875 -0.453125 1.75q-0.15625 0.875 -0.15625 1.640625q0 0.90625 0.234375 1.59375q0.234375 0.6875 0.65625 1.15625q0.4375 0.453125 1.046875 0.6875q0.625 0.21875 1.40625 0.21875q0.9375 0 1.59375 -0.1875q0.65625 -0.203125 1.109375 -0.4375q0.46871948 -0.234375 0.7655945 -0.421875q0.3125 -0.203125 0.5625 -0.203125q0.171875 0 0.21875 0.125q0.0625 0.125 0.0625 0.34375q0 0.09375 -0.03125 0.28125q-0.015625 0.171875 -0.0625 0.40625q-0.03125 0.21875 -0.09375 0.484375q-0.046875 0.265625 -0.109375 0.5q-0.0625 0.21875 -0.15625 0.421875q-0.078125 0.1875 -0.21875 0.328125q-0.140625 0.140625 -0.5625 0.375q-0.42184448 0.21875 -1.0155945 0.421875q-0.59375 0.1875 -1.359375 0.3125q-0.75 0.140625 -1.609375 0.140625q-1.421875 0 -2.5625 -0.375q-1.125 -0.375 -1.921875 -1.140625q-0.796875 -0.765625 -1.234375 -1.921875q-0.421875 -1.15625 -0.421875 -2.71875q0 -1.203125 0.25 -2.5q0.25 -1.3125 0.75 -2.5625q0.5 -1.25 1.265625 -2.375q0.765625 -1.125 1.8125 -1.953125q1.046875 -0.84375 2.375 -1.328125q1.34375 -0.5 2.96875 -0.5q0.9062195 0 1.7187195 0.1875q0.828125 0.1875 1.421875 0.5q0.609375 0.3125 0.84375 0.5625q0.25 0.25 0.25 0.640625zm10.213379 14.796875q-0.046875 0.265625 -0.40625 0.390625q-0.359375 0.109375 -1.15625 0.109375q-0.40625 0 -0.65625 -0.03125q-0.25 -0.015625 -0.40625 -0.078125q-0.15625 -0.0625 -0.203125 -0.15625q-0.046875 -0.09375 -0.015625 -0.234375l0.359375 -1.875q-0.203125 0.421875 -0.640625 0.890625q-0.421875 0.453125 -1.0 0.84375q-0.578125 0.375 -1.296875 0.609375q-0.703125 0.25 -1.46875 0.25q-1.03125 0 -1.734375 -0.359375q-0.6875 -0.359375 -1.109375 -0.984375q-0.40625 -0.625 -0.578125 -1.40625q-0.171875 -0.796875 -0.171875 -1.671875q0 -0.828125 0.15625 -1.796875q0.15625 -0.96875 0.46875 -1.921875q0.328125 -0.96875 0.828125 -1.859375q0.515625 -0.890625 1.21875 -1.578125q0.71875 -0.6875 1.640625 -1.09375q0.921875 -0.421875 2.078125 -0.421875q1.078125 0 1.90625 0.453125q0.828125 0.453125 1.453125 1.3125l0.21875 -1.03125q0.046875 -0.28125 0.40625 -0.390625q0.359375 -0.125 1.15625 -0.125q0.390625 0 0.640625 0.03125q0.265625 0.03125 0.40625 0.09375q0.15625 0.0625 0.203125 0.15625q0.0625 0.09375 0.03125 0.234375l-2.328125 11.640625zm-1.609375 -8.15625q-0.453125 -0.71875 -1.015625 -1.078125q-0.546875 -0.359375 -1.28125 -0.359375q-0.546875 0 -1.0 0.265625q-0.453125 0.265625 -0.8125 0.703125q-0.34375 0.4375 -0.609375 1.015625q-0.265625 0.5625 -0.453125 1.1875q-0.171875 0.609375 -0.265625 1.25q-0.078125 0.625 -0.078125 1.15625q0 0.390625 0.0625 0.75q0.0625 0.34375 0.21875 0.625q0.15625 0.265625 0.40625 0.4375q0.265625 0.15625 0.671875 0.15625q0.578125 0 1.171875 -0.34375q0.59375 -0.34375 1.109375 -0.9375q0.53125 -0.59375 0.921875 -1.375q0.40625 -0.796875 0.578125 -1.6875l0.375 -1.765625zm14.538757 -3.65625q0 0.09375 -0.03125 0.328125q-0.03125 0.234375 -0.078125 0.546875q-0.046875 0.3125 -0.125 0.640625q-0.078125 0.328125 -0.1875 0.625q-0.09375 0.28125 -0.21875 0.46875q-0.125 0.171875 -0.265625 0.171875q-0.109375 0 -0.234375 -0.046875q-0.109375 -0.046875 -0.265625 -0.078125q-0.140625 -0.046875 -0.3125 -0.09375q-0.15625 -0.046875 -0.375 -0.046875q-0.453125 0 -0.953125 0.34375q-0.5 0.328125 -0.953125 0.890625q-0.453125 0.546875 -0.8125 1.28125q-0.359375 0.734375 -0.515625 1.546875l-1.046875 5.21875q-0.015625 0.125 -0.125 0.234375q-0.09375 0.09375 -0.296875 0.15625q-0.203125 0.0625 -0.53125 0.09375q-0.3125 0.03125 -0.78125 0.03125q-0.484375 0 -0.796875 -0.03125q-0.296875 -0.03125 -0.484375 -0.09375q-0.171875 -0.0625 -0.21875 -0.15625q-0.046875 -0.109375 -0.03125 -0.234375l2.328125 -11.625q0.015625 -0.140625 0.109375 -0.234375q0.09375 -0.09375 0.265625 -0.15625q0.171875 -0.0625 0.453125 -0.09375q0.28125 -0.03125 0.671875 -0.03125q0.40625 0 0.65625 0.03125q0.265625 0.03125 0.40625 0.09375q0.140625 0.0625 0.1875 0.15625q0.046875 0.09375 0.015625 0.234375l-0.359375 1.84375q0.3125 -0.5625 0.71875 -1.03125q0.40625 -0.46875 0.859375 -0.8125q0.46875 -0.34375 0.96875 -0.53125q0.5 -0.203125 0.96875 -0.203125q0.203125 0 0.40625 0.03125q0.203125 0.015625 0.375 0.0625q0.171875 0.03125 0.3125 0.09375q0.140625 0.0625 0.21875 0.125q0.078125 0.078125 0.078125 0.25zm9.384338 0q0 0.09375 -0.03125 0.328125q-0.03125 0.234375 -0.078125 0.546875q-0.046875 0.3125 -0.125 0.640625q-0.078125 0.328125 -0.1875 0.625q-0.09375 0.28125 -0.21875 0.46875q-0.125 0.171875 -0.265625 0.171875q-0.109375 0 -0.234375 -0.046875q-0.109375 -0.046875 -0.265625 -0.078125q-0.140625 -0.046875 -0.3125 -0.09375q-0.15625 -0.046875 -0.375 -0.046875q-0.453125 0 -0.953125 0.34375q-0.5 0.328125 -0.953125 0.890625q-0.453125 0.546875 -0.8125 1.28125q-0.359375 0.734375 -0.515625 1.546875l-1.046875 5.21875q-0.015625 0.125 -0.125 0.234375q-0.09375 0.09375 -0.296875 0.15625q-0.203125 0.0625 -0.53125 0.09375q-0.3125 0.03125 -0.78125 0.03125q-0.484375 0 -0.796875 -0.03125q-0.296875 -0.03125 -0.484375 -0.09375q-0.171875 -0.0625 -0.21875 -0.15625q-0.046875 -0.109375 -0.03125 -0.234375l2.328125 -11.625q0.015625 -0.140625 0.109375 -0.234375q0.09375 -0.09375 0.265625 -0.15625q0.171875 -0.0625 0.453125 -0.09375q0.28125 -0.03125 0.671875 -0.03125q0.40625 0 0.65625 0.03125q0.265625 0.03125 0.40625 0.09375q0.140625 0.0625 0.1875 0.15625q0.046875 0.09375 0.015625 0.234375l-0.359375 1.84375q0.3125 -0.5625 0.71875 -1.03125q0.40625 -0.46875 0.859375 -0.8125q0.46875 -0.34375 0.96875 -0.53125q0.5 -0.203125 0.96875 -0.203125q0.203125 0 0.40625 0.03125q0.203125 0.015625 0.375 0.0625q0.171875 0.03125 0.3125 0.09375q0.140625 0.0625 0.21875 0.125q0.078125 0.078125 0.078125 0.25zm6.3687134 -3.71875q-0.09375 0.484375 -0.25 0.8125q-0.15625 0.328125 -0.421875 0.53125q-0.25 0.203125 -0.625 0.296875q-0.375 0.078125 -0.9375 0.078125q-0.546875 0 -0.890625 -0.078125q-0.34375 -0.09375 -0.53125 -0.296875q-0.171875 -0.203125 -0.203125 -0.53125q-0.015625 -0.328125 0.09375 -0.8125q0.09375 -0.46875 0.25 -0.796875q0.15625 -0.328125 0.40625 -0.53125q0.265625 -0.21875 0.640625 -0.296875q0.375 -0.09375 0.9375 -0.09375q0.546875 0 0.890625 0.09375q0.34375 0.078125 0.515625 0.296875q0.1875 0.203125 0.203125 0.53125q0.015625 0.328125 -0.078125 0.796875zm-3.359375 15.515625q-0.015625 0.125 -0.125 0.234375q-0.09375 0.09375 -0.296875 0.15625q-0.203125 0.0625 -0.53125 0.09375q-0.3125 0.03125 -0.78125 0.03125q-0.484375 0 -0.796875 -0.03125q-0.296875 -0.03125 -0.484375 -0.09375q-0.171875 -0.0625 -0.21875 -0.15625q-0.046875 -0.109375 -0.03125 -0.234375l2.3125 -11.59375q0.015625 -0.125 0.109375 -0.21875q0.109375 -0.109375 0.3125 -0.171875q0.203125 -0.078125 0.515625 -0.109375q0.328125 -0.046875 0.796875 -0.046875q0.484375 0 0.78125 0.046875q0.3125 0.03125 0.484375 0.109375q0.1875 0.0625 0.234375 0.171875q0.0625 0.09375 0.03125 0.21875l-2.3125 11.59375zm15.375061 -9.0625q0 0.9375 -0.421875 1.6875q-0.421875 0.75 -1.3125 1.296875q-0.890625 0.53125 -2.265625 0.8125q-1.375 0.28125 -3.265625 0.28125l-1.203125 0q-0.046875 0.3125 -0.078125 0.609375q-0.03125 0.28125 -0.03125 0.53125q0 1.046875 0.546875 1.59375q0.546875 0.53125 1.765625 0.53125q0.859375 0 1.53125 -0.125q0.6875 -0.125 1.203125 -0.265625q0.515625 -0.15625 0.84375 -0.28125q0.34375 -0.125 0.5 -0.125q0.15625 0 0.21875 0.09375q0.078125 0.09375 0.078125 0.28125q0 0.21875 -0.046875 0.515625q-0.03125 0.28125 -0.109375 0.5625q-0.078125 0.28125 -0.1875 0.53125q-0.109375 0.234375 -0.25 0.375q-0.15625 0.15625 -0.578125 0.3125q-0.40625 0.15625 -1.0 0.28125q-0.578125 0.125 -1.3125 0.203125q-0.71875 0.09375 -1.484375 0.09375q-1.25 0 -2.1875 -0.265625q-0.921875 -0.25 -1.546875 -0.78125q-0.625 -0.53125 -0.9375 -1.359375q-0.3125 -0.828125 -0.3125 -1.96875q0 -0.875 0.15625 -1.859375q0.171875 -1.0 0.53125 -1.96875q0.375 -0.96875 0.9375 -1.84375q0.578125 -0.890625 1.390625 -1.5625q0.8125 -0.671875 1.875 -1.078125q1.078125 -0.40625 2.4375 -0.40625q1.1875 0 2.03125 0.28125q0.859375 0.28125 1.40625 0.75q0.546875 0.453125 0.8125 1.046875q0.265625 0.59375 0.265625 1.21875zm-3.21875 0.25q0 -0.53125 -0.40625 -0.875q-0.40625 -0.34375 -1.1875 -0.34375q-0.65625 0 -1.171875 0.234375q-0.5 0.234375 -0.890625 0.640625q-0.390625 0.390625 -0.671875 0.9375q-0.28125 0.53125 -0.453125 1.140625l1.09375 0q1.03125 0 1.734375 -0.125q0.703125 -0.140625 1.125 -0.375q0.4375 -0.234375 0.625 -0.546875q0.203125 -0.328125 0.203125 -0.6875zm19.979492 6.953125q0 0.078125 -0.015625 0.28125q-0.015625 0.203125 -0.0625 0.46875q-0.046875 0.265625 -0.125 0.546875q-0.0625 0.265625 -0.171875 0.5q-0.09375 0.21875 -0.21875 0.359375q-0.125 0.140625 -0.296875 0.140625l-7.71875 0q-0.203125 0 -0.390625 -0.046875q-0.171875 -0.0625 -0.28125 -0.203125q-0.109375 -0.140625 -0.15625 -0.34375q-0.03125 -0.203125 0.015625 -0.484375l3.078125 -15.3125q0.015625 -0.125 0.109375 -0.21875q0.109375 -0.109375 0.328125 -0.171875q0.21875 -0.078125 0.546875 -0.109375q0.328125 -0.03125 0.84375 -0.03125q0.484375 0 0.8125 0.03125q0.328125 0.03125 0.5 0.109375q0.1875 0.0625 0.25 0.171875q0.0625 0.09375 0.03125 0.21875l-2.703125 13.5625l5.3125 0q0.171875 0 0.234375 0.15625q0.078125 0.140625 0.078125 0.375zm12.3185425 1.875q-0.046875 0.265625 -0.40625 0.390625q-0.359375 0.109375 -1.15625 0.109375q-0.40625 0 -0.65625 -0.03125q-0.25 -0.015625 -0.40625 -0.078125q-0.15625 -0.0625 -0.203125 -0.15625q-0.046875 -0.09375 -0.015625 -0.234375l0.359375 -1.875q-0.203125 0.421875 -0.640625 0.890625q-0.421875 0.453125 -1.0 0.84375q-0.578125 0.375 -1.296875 0.609375q-0.703125 0.25 -1.46875 0.25q-1.03125 0 -1.734375 -0.359375q-0.6875 -0.359375 -1.109375 -0.984375q-0.40625 -0.625 -0.578125 -1.40625q-0.171875 -0.796875 -0.171875 -1.671875q0 -0.828125 0.15625 -1.796875q0.15625 -0.96875 0.46875 -1.921875q0.328125 -0.96875 0.828125 -1.859375q0.515625 -0.890625 1.21875 -1.578125q0.71875 -0.6875 1.640625 -1.09375q0.921875 -0.421875 2.078125 -0.421875q1.078125 0 1.90625 0.453125q0.828125 0.453125 1.453125 1.3125l0.21875 -1.03125q0.046875 -0.28125 0.40625 -0.390625q0.359375 -0.125 1.15625 -0.125q0.390625 0 0.640625 0.03125q0.265625 0.03125 0.40625 0.09375q0.15625 0.0625 0.203125 0.15625q0.0625 0.09375 0.03125 0.234375l-2.328125 11.640625zm-1.609375 -8.15625q-0.453125 -0.71875 -1.015625 -1.078125q-0.546875 -0.359375 -1.28125 -0.359375q-0.546875 0 -1.0 0.265625q-0.453125 0.265625 -0.8125 0.703125q-0.34375 0.4375 -0.609375 1.015625q-0.265625 0.5625 -0.453125 1.1875q-0.171875 0.609375 -0.265625 1.25q-0.078125 0.625 -0.078125 1.15625q0 0.390625 0.0625 0.75q0.0625 0.34375 0.21875 0.625q0.15625 0.265625 0.40625 0.4375q0.265625 0.15625 0.671875 0.15625q0.578125 0 1.171875 -0.34375q0.59375 -0.34375 1.109375 -0.9375q0.53125 -0.59375 0.921875 -1.375q0.40625 -0.796875 0.578125 -1.6875l0.375 -1.765625zm15.663757 8.15625q-0.015625 0.125 -0.109375 0.21875q-0.078125 0.09375 -0.265625 0.15625q-0.171875 0.0625 -0.453125 0.09375q-0.265625 0.03125 -0.671875 0.03125q-0.390625 0 -0.65625 -0.03125q-0.25 -0.03125 -0.390625 -0.09375q-0.140625 -0.0625 -0.1875 -0.15625q-0.046875 -0.09375 -0.015625 -0.21875l0.359375 -1.890625q-0.265625 0.453125 -0.71875 0.921875q-0.453125 0.453125 -1.046875 0.84375q-0.59375 0.375 -1.28125 0.609375q-0.6875 0.234375 -1.421875 0.234375q-0.90625 0 -1.546875 -0.25q-0.625 -0.25 -1.015625 -0.6875q-0.375 -0.4375 -0.5625 -1.046875q-0.171875 -0.609375 -0.171875 -1.328125q0 -0.453125 0.0625 -0.90625q0.0625 -0.46875 0.15625 -0.953125l1.421875 -7.1875q0.03125 -0.140625 0.125 -0.234375q0.109375 -0.09375 0.296875 -0.15625q0.203125 -0.0625 0.53125 -0.09375q0.328125 -0.03125 0.796875 -0.03125q0.484375 0 0.78125 0.03125q0.3125 0.03125 0.484375 0.09375q0.171875 0.0625 0.21875 0.15625q0.0625 0.09375 0.046875 0.234375l-1.390625 6.96875q-0.0625 0.296875 -0.09375 0.609375q-0.03125 0.3125 -0.03125 0.59375q0 0.3125 0.046875 0.5625q0.0625 0.25 0.1875 0.4375q0.140625 0.1875 0.34375 0.28125q0.203125 0.078125 0.484375 0.078125q0.53125 0 1.109375 -0.34375q0.578125 -0.34375 1.09375 -0.9375q0.53125 -0.609375 0.9375 -1.4375q0.421875 -0.828125 0.625 -1.796875l0.984375 -5.015625q0.015625 -0.140625 0.109375 -0.234375q0.109375 -0.09375 0.3125 -0.15625q0.203125 -0.0625 0.53125 -0.09375q0.328125 -0.03125 0.796875 -0.03125q0.484375 0 0.78125 0.03125q0.3125 0.03125 0.484375 0.09375q0.171875 0.0625 0.21875 0.15625q0.0625 0.09375 0.03125 0.234375l-2.328125 11.640625z");

		TheGlass.mirrors(400);
		borderbox TGB = XSolid(XObjs::Extrusion_xOy(TheGlass, 0.01)).border;
		auto TheGlass_Transform = [](Figure2D &F, borderbox A, double t) -> XSolid {
			F.mirrors(400);
			F += -vec2(A.center().x, A.center().y); F *= 6 / (A.Max.x - A.Min.x); F += vec2(5.8, 4.2);
			XSolid R(XObjs::Extrusion_xOy(F, 1.5, 2));
			R = CSG_Rotation(R, PI / 2, 0, 0);
			R = CSG_Translation(R, mix(vec3(0, -15, 5), vec3(0, 0, 0), sigmoid(8 * (t - 10.5))));
			R.type = XSolid_LightingCrystal; R.col = drgb(0.02, 0.025, 0.03);
			return R;
		};
		World THE_GLASS; THE_GLASS.add({
			&TheGlass_Transform(TheGlass_T, TGB, t), &TheGlass_Transform(TheGlass_H, TGB, t), &TheGlass_Transform(TheGlass_E, TGB, t),
			&TheGlass_Transform(TheGlass_G, TGB, t), &TheGlass_Transform(TheGlass_L, TGB, t), &TheGlass_Transform(TheGlass_A, TGB, t),
			&TheGlass_Transform(TheGlass_S1, TGB, t), &TheGlass_Transform(TheGlass_S2, TGB, t)
		});
		World BY_NAMES; BY_NAMES.add({ &TheGlass_Transform(ByNames1, TGB, t), &TheGlass_Transform(ByNames2, TGB, t) });

		if (t > 10) W->insert(&THE_GLASS)/*, W->insert(&BY_NAMES)*/;

	}, [](double t) -> point {	// camera
		t = mix(-2.5, -PI / 2, sigmoid(20 * (t - 6)));
		point P(20 * cos(t) + 6, 20 * sin(t), 8);
		if (t < 9) return P;
		return mix(P, point(6, -20, 2.0), sigmoid(10 * (t - 10.5)));
	}, [](double t) -> point {	// view point
		return mix(point(6, 1.5, 0), point(6, 1.5, 2.0), sigmoid(10 * (t - 10)));
	}, [](double t) -> double {		// solid angle
		return 0.15;
	}, 9, 12, width, height, fps, sampling, "Animation\\CPTAR0_", 3);	// beginning

#pragma endregion

#pragma region R1
	if (R == -1 || R == 1) CPT_Animation([](double t, World* W) {	// World
		parallelogram_grid P(parallelogram(point(-3, -ERR_UPSILON, 0), point(7, 0), point(0, 2 * ERR_UPSILON))); W->add(P);
		GlassMan_std G0; SetWalkingAttitudes(t, G0); G0.Pos += point(-0.4, 0); G0.push(*W);
		GlassMan_std G1; SetWalkingAttitudes(t + 0.4, G1); G1.Pos += point(3.2, 12); G1.push(*W);
		GlassMan_std G; SetWalkingAttitudes(t - 1, G); G.Pos += point(1.8, -8); G.push(*W);		// Protagonist, t = 12 => (1.77834,1.62188,0.94572)
		GlassMan_std G3; SetWalkingAttitudes(t - 0.7, G3); G3.Pos += point(-1.2, 32); G3.push(*W);
		GlassMan_std G4; SetWalkingAttitudes(t + 1.1, G4); G4.Pos += point(3.6, 80); G4.push(*W);
		GlassMan_std G5; SetWalkingAttitudes(t + 0.2, G5); G5.Pos += point(-0.6, -22); G5.push(*W);
		//ADD_AXIS(*W, 0.1);
	}, [](double t) -> point {	// camera
		if (t > 0 && t <= 4) return mix(point(10, -20, 15), point(-4, -30, 2), tanh(t));
		if (t > 4 && t <= 8) return Bezier(point(-4, -30, 2), point(-30, 10, 8), point(-2, 20, 2), tanh(t - 4) * (1 - exp(-(t - 4)*(t - 4))));
		if (t > 8 && t <= 12) return mix(point(-2, 20, 2), mix(point(-2, 20, 2), ([](double t)->point {
			GlassMan_std G; SetWalkingAttitudes(t - 1, G); G.Pos += point(1.8, -8); return G.construct(); })(t), 0.5), tanh(t - 8) * (1 - exp(-(t - 8)*(t - 8))));
		return NAP;
	}, [](double t) -> point {	// view point
		if (t > 0 && t <= 4) return point(0, 0, 0.8 - 0.8 * tanh(t));
		if (t > 4 && t <= 8) return Bezier(point(0, 0, 0), point(0, -5, 3), point(0, -10, 0), tanh(t - 4));
		if (t > 8 && t <= 12) return mix(point(0, -10, 0), ([](double t)->point {
			GlassMan_std G; SetWalkingAttitudes(t - 1, G); G.Pos += point(1.8, -8); return G.construct(); })(t), tanh(t - 8) * (1 - exp(-(t - 8)*(t - 8))));
		return NAP;
	}, [](double t) -> double {		// solid angle
		if (t > 0 && t <= 4) return 0.2 / (2 * tanh(t) + 1);
		if (t > 4 && t <= 8) return Bezier(0.2 / 3, 1.2, 0.2, tanh(t - 4) * (1 - exp(-(t - 4)*(t - 4))));
		if (t > 8 && t <= 12) return mix(0.2, 0.08, sigmoid(4 * ((t - 8) - 2.5)));
		return NAN;
	}, 0, 12, width, height, fps, sampling, "Animation\\CPTAR1_", 3);	// 0s-12s
#pragma endregion The glass man walks on a path with several identical glass mans.

#pragma region R2
	if (R == -1 || R == 2 || R == 201) CPT_Animation([](double t, World* W) {	// World
		parallelogram_grid P(parallelogram(point(-3, -ERR_UPSILON, 0), point(7, 0), point(0, 2 * ERR_UPSILON))); W->add(P);
		GlassMan_std G0; SetWalkingAttitudes(t, G0); G0.Pos += point(-0.4, 0); G0.push(*W);
		GlassMan_std G1; SetWalkingAttitudes(t + 0.4, G1); G1.Pos += point(3.2, 12); G1.push(*W);
		GlassMan_std G3; SetWalkingAttitudes(t - 0.7, G3); G3.Pos += point(-1.2, 32); G3.push(*W);
		GlassMan_std G4; SetWalkingAttitudes(t + 1.1, G4); G4.Pos += point(3.6, 80); G4.push(*W);
		GlassMan_std G5; SetWalkingAttitudes(t + 0.2, G5); G5.Pos += point(-0.6, -22); G5.push(*W);
		if (t > 12.68) t = 12.68;
		GlassMan_std G; SetWalkingAttitudes(t - 1, G); G.Pos += point(1.8, -8);		// Protagonist
		mix(G, ManStaringAtTheGoldenCoin_Constructor(point(1.8, 2.3), point(6.2, 7.2, -0.1)), clamp(2 * (t - 12.3), 0.0, 1.0)).push(*W);
		parallelogram_grid Ps(parallelogram(point(5, 6, -0.5), point(2.4, 0), point(0, 2.4)), 0.8, 0.8); W->add(Ps);
		W->insert(&GoldenCoin_Constructor(point(6.2, 7.2, -0.1), 0.2));
		//ADD_AXIS(*W, 0.1);
	}, [](double t) -> point {	// camera
		if (t < 14) return mix(point(-2, 20, 2), ([](double t)->point {
			GlassMan_std G; SetWalkingAttitudes(t - 1, G); G.Pos += point(1.8, -8);
			G = mix(G, ManStaringAtTheGoldenCoin_Constructor(point(1.8, 2.3), point(6.2, 7.2, -0.1)), clamp(2 * (t - 12.3), 0.0, 1.0));
			return G.construct(); })(t), 0.5);	// max t: (-0.0963206,11.1622,1.4699)
		return Catmull_Rom(vector<point>({ point(-0.0963206,11.1622,1.4699), point(-6.30194,16.31809,0.93712),
			point(-12.26837,4.85286,0.86076), point(-10.92227,-7.02197,2), point(-4, -6, 4) }), 4 * sigmoid(12 * ((t - 14) - 0.5)));
		return point(-4, -6, 4);
	}, [](double t) -> point {	// view point
		if (t < 14) return ([](double t)->point {
			GlassMan_std G; SetWalkingAttitudes(t - 1, G); G.Pos += point(1.8, -8);
			G = mix(G, ManStaringAtTheGoldenCoin_Constructor(point(1.8, 2.3), point(6.2, 7.2, -0.1)), clamp(2 * (t - 12.3), 0.0, 1.0));
			return G.construct(); })(t);	// max t: (1.80736,2.32443,0.939794)
		else return mix(point(1.80736, 2.32443, 0.939794), point(6.2, 7.2, -0.4), sigmoid(20 * ((t - 14) - 0.6)));
		return point(1.77982, 3.37131, 0.946132);
	}, [](double t) -> double {		// solid angle
		return 0.08;
	}, 12, 15, width, height, fps, sampling, "Animation\\CPTAR2.1_", 3);	// 12s-16s, still the last theme for 1s

	if (R == -1 || R == 2 || R == 202) CPT_Animation([](double t, World* W) {	// World
		parallelogram_grid P(parallelogram(point(-3, -ERR_UPSILON, 0), point(7, 0), point(0, 2 * ERR_UPSILON))); W->add(P);
		GlassMan_std G0; SetWalkingAttitudes(t, G0); G0.Pos += point(-0.4, 0); G0.push(*W);
		GlassMan_std G1; SetWalkingAttitudes(t + 0.4, G1); G1.Pos += point(3.2, 12); G1.push(*W);
		GlassMan_std G3; SetWalkingAttitudes(t - 0.7, G3); G3.Pos += point(-1.2, 32); G3.push(*W);
		GlassMan_std G4; SetWalkingAttitudes(t + 1.1, G4); G4.Pos += point(3.6, 80); G4.push(*W);
		GlassMan_std G5; SetWalkingAttitudes(t + 0.2, G5); G5.Pos += point(-0.6, -22); G5.push(*W);
		ManStaringAtTheGoldenCoin_Constructor(point(1.8, 2.3), point(6.2, 7.2, -0.1)).push(*W);		// Protagonist
		auto vertical = [](double t) -> double {
			t -= 16;
			return (t < PI * 2. / 3) ? (0.05 * pow(sin(3 * t), 2)) : 0;
		};
		parallelogram_grid Ps(parallelogram(point(5, 6, vertical(t) - 0.5), point(2.4, 0), point(0, 2.4)), 0.8, 0.8); W->add(Ps);
		W->insert(&GoldenCoin_Constructor(point(6.2, 7.2, vertical(t) - 0.1), 0.2 + ((t > 18.5 && t < 19.5) ? 2 * (t - 18.5)*PI : 0)));
	}, [](double t) -> point {	// camera
		return point(0, 3, 2.8);
	}, [](double t) -> point {	// view point
		return point(6.2, 7.2, -0.4);
	}, [](double t) -> double {		// solid angle
		return 0.15;
	}, 16, 20, width, height, fps, sampling, "Animation\\CPTAR2.2_", 3);	// 12s-16s, still the last theme for 1s

#pragma endregion A small lower platform with same material floats lower under the path. As the glass man passes by, a golden coin shines on it.

#pragma region R3
	if (R == -1 || R == 3) CPT_Animation([](double t, World* W) {	// World
		parallelogram_grid P(parallelogram(point(-3, -ERR_UPSILON, 0), point(7, 0), point(0, 2 * ERR_UPSILON))); W->add(P);
		GlassMan_std G0; SetWalkingAttitudes(t, G0); G0.Pos += point(-0.4, 0); G0.push(*W);
		GlassMan_std G1; SetWalkingAttitudes(t + 0.4, G1); G1.Pos += point(3.2, 12); G1.push(*W);
		GlassMan_std G3; SetWalkingAttitudes(t - 0.7, G3); G3.Pos += point(-1.2, 32); G3.push(*W);
		GlassMan_std G4; SetWalkingAttitudes(t + 1.1, G4); G4.Pos += point(3.6, 80); G4.push(*W);
		GlassMan_std G5; SetWalkingAttitudes(t + 0.2, G5); G5.Pos += point(-0.6, -22); G5.push(*W);
		ManStaringAtTheGoldenCoin_Constructor(point(1.8, 2.3), point(6.2, 7.2, -0.1)).push(*W);		// Protagonist
		auto vertical = [](double t) -> double {
			t -= 16;
			return (t < PI * 2. / 3) ? (0.05 * pow(sin(3 * t), 2)) : 0;
		};
		parallelogram_grid Ps(parallelogram(point(5, 6, vertical(t) - 0.5), point(2.4, 0), point(0, 2.4)), 0.8, 0.8); W->add(Ps);
		W->insert(&GoldenCoin_Constructor(point(6.2, 7.2, vertical(t) - 0.1), 0.2 + ((t > 18.5 && t < 19.5) ? 2 * (t - 18.5)*PI : 0)));
	}, [](double t) -> point {	// camera
		return point(0, 3, 2.8);
	}, [](double t) -> point {	// view point
		return point(6.2, 7.2, -0.4);
	}, [](double t) -> double {		// solid angle
		return 0.15;
	}, 16, 20, width, height, fps, sampling, "Animation\\CPTAR2.2_", 3);	// 20s-##s
#pragma endregion The glass man saw it and stopped, then jumped onto the platform to get the golden coin.

}
