
#include "World.h"
#include "GlassMan.h"
#include "Test.cpp";

using namespace std;

void CPT_T1() {
	World W;
	plane_grid P(0); W.add(&P);

	GlassMan G;

	G.Heel_l = point(0, 1, 0.3), G.Heel_r = point(0, -1, 0.5);
	G.foot_dir_l = point(2, 0.2, -0.3), G.foot_dir_r = point(2, -0.1, -0.5);
	G.width_of_heel = 1.0, G.width_of_tiptoe = 0.8, G.height_of_heel = 0.6; G.foot_rounding = 0.15;
	G.length_of_foot = 3;

	G.Knee_l = point(0.3, 1, 4.8), G.Knee_r = point(0.5, -0.9, 4.8);
	G.upper_radius_of_shank = 0.4, G.lower_radius_of_shank = 0.24, G.shank_rounding = 0.2;

	G.Buttock_l = point(0, 0.7, 7.8), G.Buttock_r = point(0, -0.6, 7.8);
	G.upper_radius_of_thigh = 0.6, G.lower_radius_of_thigh = 0.4, G.thigh_rounding = 0.3;

	//sphere S1(G.Knee_r, 0.05); S1.setcolor(Red); W.add(&S1);
	//sphere S2(G.Buttock_r, 0.05); S2.setcolor(Red); W.add(&S2);

	G.construct(); G.push(W);

	parallelogram Pr(point(-6, -4, 8), point(12, 0, 0), point(0, 8, 0)); //Pr.setcolor(White), W.add(&Pr);
	VisualizeSDF(G.construct(), Pr, 600 * 400);

	bitmap img(600, 400);
	//ADD_AXIS(W, 0.1);
	W.setGlobalLightSource(0, 0, 1);
	W.Render_Sampling = 1;
	W.render(img, point(200, -100, 100), point(0, 0, 4), 0, 0.003);
	img.out("IMAGE\\CPT.bmp");
}
