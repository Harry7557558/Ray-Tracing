
#include "World.h"
#include "GlassMan.h"
#include "Test.cpp"

using namespace std;

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

	G.push(W);
	W.setGlobalLightSource(0, 0, 1);
	W.Render_Sampling = 1;
	bitmap img(400, 600);
	W.render(img, point(200, -100, 100), point(0, 0, 0.8), 0, 0.00005);
	img.out("IMAGE\\CPT.bmp");

}

void CPT_Animation_T1() {
	const string HEADER = "Animation\\CPTAT1_";
	const unsigned N = 3;
	const unsigned FRAMES = 32;

	World W;
	plane_grid P(0); W.add(&P);
	GlassMan_std G;

	W.setGlobalLightSource(0, 0, 1);
	W.Render_Sampling = 2;


	double t = 0;
	vector<World> Ws; vector<point> Cs, Vs; vector<double> Va;
	vector<thread> Ts; vector<byte> Rs;
	for (unsigned i = 0; i < FRAMES; i++) {
		t = 0.1*i;

		G.Pos = point(0, 0, 0), G.top = vec3(0, 0, 1), G.dir = vec3(1, 0, 0);
		G.v_head = vec3(0, 0, 1), G.v_head_side = vec3(-sin(0.2), cos(0.2), 0), G.v_neck = vec3(0, 0, 1);
		G.v_chest = vec3(0, 0, 1), G.v_chest_side = vec3(-sin(0.1), cos(0.1), 0), G.v_waist = vec3(0, 0, -1), G.v_waist_side = vec3(0, 1, 0);
		G.v_upperarm_l = vec3(0, 0.1, -1), G.v_forearm_l = vec3(0.2, -0.2, -1); G.v_upperarm_r = vec3(-0.05, -0.15, -1), G.v_forearm_r = vec3(0.3, 0.2, -1);
		G.v_thigh_l = vec3(0.05, 0, -1), G.v_shank_l = vec3(0, 0, -1); G.v_thigh_r = vec3(0.15, -0.15, -1), G.v_shank_r = vec3(0.1, -0.1, -1);
		G.v_foot_l = vec3(cos(0.1), sin(0.1), -0.1), G.v_foot_l_side = vec3(0, 1, 0); G.v_foot_r = vec3(cos(0.7), -sin(0.7), 0), G.v_foot_r_side = vec3(-sin(0.7), -cos(0.7), 0);
		G.auto_fit = true;

		W.clear(); W.add(&P); G.push(W);
		Ws.push_back(W); Cs.push_back(point(200, -100, 15*(0.25*i+2))), Vs.push_back(point(0, 0, 0.8)), Va.push_back(0.00005);
		
		Rs.push_back(0);
		Ts.push_back(thread([&](World _W, unsigned no, unsigned N) {
			bitmap img(800, 1200);
			_W.render(img, Cs.at(no), Vs.at(no), 0, Va.at(no));
			img.out(&(HEADER + uint2str(no, N) + ".bmp")[0]);
			Rs.at(no) = 1;
		}, Ws.at(i), i, N));
		Ts.at(i).detach();
	}

	while (1) {
		for (int i = 0; i < FRAMES; i++) {
			if (Rs.at(i) == 0) break;
			if (i + 1 == FRAMES) return;
		}
		Sleep(50);
	}
}
