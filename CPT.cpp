
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
