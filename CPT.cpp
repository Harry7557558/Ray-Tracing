
#include "World.h"
#include "GlassMan.h"
#include "Test.cpp"

using namespace std;

void CPT_T1() {
	World W;
	plane_grid P(0); W.add(&P);

	GlassMan G;
	//G.top = point(1, 0, 0), G.dir = point(0, -1, 1);

	G.Heel_l = point(0, 1, 0.3), G.Heel_r = point(0, -1, 0.5);
	G.foot_dir_l = point(2, 0.2, -0.3), G.foot_dir_r = point(2, -0.1, -0.5);
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
	G.hand_radius = 0.3;

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

	bitmap img(1080, 1920);
	//ADD_AXIS(W, 0.1);
	W.setGlobalLightSource(0, 0, 1);
	W.Render_Sampling = 2;
	//W.render(img, point(200, -100, 100), point(0, 1, 0), 0, 0.0001);
	W.render(img, point(200, -100, 100), point(0, 0, 8), 0, 0.004);
	img.out("IMAGE\\CPT.bmp");
}
