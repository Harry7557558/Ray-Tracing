#include "bitmap.h"

// Manual statistic data, for interpolation/fitting
namespace WalkingMan {
	enum comps { head_t, head_l, shdlr_c, elbow_l, elbow_r, hand_l, hand_r, waist, butt, knee_l, knee_r, heel_l, heel_r, ttoe_l, ttoe_r };		// total 15 points
	const double _t[25][15][2] = {
		//{ {71,522}, {74,455}, {60,420}, {NAN,NAN}, {65,336}, {NAN,NAN}, {109,265}, {51,313}, {56,254}, {92,157}, {37,169}, {102,49}, {NAN,NAN}, {140,30}, {NAN,NAN} },
		//{ {83,522}, {96,455}, {81,418}, {NAN,NAN}, {84,333}, {NAN,NAN}, {117,263}, {67,313}, {73,254}, {99,154}, {88,166}, {102,46}, {NAN,NAN}, {143,28}, {NAN,NAN} },
		//{ {109,523}, {120,456}, {100,416}, {NAN,NAN}, {101,338}, {NAN,NAN}, {126,260}, {98,313}, {108,245}, {112,161}, {127,169}, {102,38}, {29,87}, {150,27}, {50,36} },
		//{ {141,526}, {153,457}, {128,415}, {NAN,NAN}, {122,336}, {180,266}, {135,261}, {131,304}, {142,243}, {129,165}, {180,179}, {101,39}, {127,54}, {145,28}, {176,28} },
		//{ {164,525}, {174,454}, {147,411}, {NAN,NAN}, {135,334}, {219,270}, {146,257}, {157,317}, {160,242}, {136,154}, {202,176}, {101,36}, {208,40}, {154,24}, {256,34} },
		{ {201,516}, {212,444}, {184,401}, {222,325}, {166,332}, {267,272}, {168,248}, {188,309}, {185,240}, {153,154}, {232,160}, {104,41}, {287,38}, {153,24}, {338,55} },	// reserve
		{ {222,514}, {237,446}, {210,398}, {242,325}, {194,327}, {288,268}, {191,243}, {220,294}, {216,233}, {183,157}, {256,159}, {110,48}, {297,35}, {153,26}, {354,36} },
		{ {276,522}, {293,451}, {270,407}, {NAN,NAN}, {262,331}, {NAN,NAN}, {263,245}, {280,305}, {287,241}, {296,157}, {296,157}, {187,88}, {302,35}, {191,29}, {355,23} },
		{ {317,526}, {331,458}, {310,410}, {NAN,NAN}, {317,331}, {NAN,NAN}, {336,251}, {309,305}, {322,243}, {364,173}, {315,161}, {289,77}, {303,35}, {312,44}, {354,23} },
		{ {341,527}, {356,455}, {341,412}, {NAN,NAN}, {354,336}, {NAN,NAN}, {386,257}, {333,313}, {346,252}, {399,177}, {322,156}, {372,50}, {305,38}, {418,26}, {358,23} },
		{ {368,522}, {382,450}, {370,410}, {NAN,NAN}, {389,333}, {NAN,NAN}, {432,260}, {360,307}, {368,241}, {420,172}, {335,153}, {448,43}, {306,38}, {502,44}, {355,24} },
		{ {396,515}, {409,445}, {400,404}, {NAN,NAN}, {424,329}, {NAN,NAN}, {472,260}, {386,300}, {392,233}, {436,160}, {345,151}, {493,41}, {308,42}, {544,57}, {356,22} },
		{ {423,511}, {436,442}, {429,400}, {NAN,NAN}, {452,325}, {NAN,NAN}, {502,261}, {415,304}, {421,238}, {466,153}, {374,150}, {493,37}, {310,46}, {553,34}, {360,22} },
		{ {462,516}, {479,449}, {470,407}, {NAN,NAN}, {488,329}, {NAN,NAN}, {525,259}, {461,305}, {477,241}, {500,146}, {462,162}, {503,32}, {354,80}, {560,22}, {365,26} },
		{ {492,518}, {508,453}, {491,409}, {NAN,NAN}, {504,331}, {NAN,NAN}, {534,258}, {489,299}, {502,240}, {509,151}, {512,158}, {505,32}, {409,89}, {558,22}, {416,36} },
		{ {533,522}, {545,455}, {525,411}, {NAN,NAN}, {521,334}, {NAN,NAN}, {539,256}, {525,311}, {530,249}, {526,167}, {573,173}, {505,37}, {505,62}, {549,25}, {552,25} },
		{ {562,522}, {568,451}, {548,410}, {NAN,NAN}, {533,332}, {605,268}, {543,254}, {551,318}, {559,240}, {536,156}, {601,172}, {507,35}, {589,39}, {555,22}, {642,26} },
		{ {589,518}, {595,446}, {569,403}, {NAN,NAN}, {553,327}, {646,272}, {550,247}, {573,304}, {577,240}, {548,153}, {622,160}, {506,35}, {666,36}, {559,21}, {716,48} },
		{ {655,517}, {664,444}, {637,398}, {NAN,NAN}, {621,326}, {703,265}, {610,245}, {655,312}, {655,240}, {632,150}, {683,150}, {537,64}, {706,33}, {564,18}, {758,22} },
		{ {681,519}, {692,451}, {669,405}, {NAN,NAN}, {656,329}, {NAN,NAN}, {653,246}, {678,311}, {687,254}, {699,153}, {697,153}, {584,82}, {706,34}, {586,26}, {756,21} },
		{ {710,521}, {718,455}, {697,407}, {NAN,NAN}, {696,332}, {NAN,NAN}, {706,248}, {706,319}, {712,261}, {741,161}, {716,162}, {646,82}, {707,33}, {662,28}, {756,21} },
		{ {735,524}, {743,457}, {727,409}, {NAN,NAN}, {736,336}, {NAN,NAN}, {761,255}, {727,315}, {734,258}, {784,174}, {726,163}, {724,50}, {705,33}, {751,26}, {760,21} },
		{ {764,518}, {781,449}, {777,407}, {NAN,NAN}, {792,337}, {NAN,NAN}, {840,263}, {764,315}, {768,250}, {820,165}, {735,149}, {841,37}, {709,36}, {891,32}, {758,22} },
		{ {797,511}, {809,444}, {805,403}, {NAN,NAN}, {830,331}, {NAN,NAN}, {884,266}, {790,313}, {792,240}, {839,151}, {750,147}, {895,36}, {710,40}, {941,51}, {762,20} },
		{ {818,508}, {836,441}, {829,401}, {NAN,NAN}, {862,330}, {NAN,NAN}, {915,265}, {817,309}, {821,242}, {866,147}, {784,162}, {902,30}, {716,48}, {954,30}, {762,21} },
		{ {860,511}, {878,446}, {876,407}, {NAN,NAN}, {894,336}, {NAN,NAN}, {943,263}, {866,311}, {875,251}, {899,150}, {863,156}, {907,29}, {759,77}, {961,17}, {772,25} },
		{ {881,513}, {904,448}, {899,409}, {NAN,NAN}, {914,330}, {NAN,NAN}, {950,260}, {888,320}, {905,254}, {905,158}, {914,162}, {907,30}, {815,89}, {959,17}, {822,37} },
		{ {907,522}, {923,452}, {918,417}, {NAN,NAN}, {928,334}, {NAN,NAN}, {954,257}, {914,320}, {917,260}, {917,161}, {960,167}, {912,29}, {877,81}, {959,18}, {899,30} },
		{ {971,516}, {987,446}, {973,406}, {NAN,NAN}, {962,327}, {1045,277}, {965,245}, {969,311}, {976,250}, {946,153}, {1023,156}, {911,29}, {1062,33}, {959,17}, {1111,37} },
		{ {1000,513}, {1013,440}, {996,397}, {NAN,NAN}, {984,321}, {1074,279}, {977,239}, {996,308}, {998,246}, {964,150}, {1044,147}, {914,37}, {1090,34}, {961,16}, {1115,34} }
		// two steps
	};
	const unsigned width = 1108, height = 542;

	void drawMan(unsigned n) {
		bitmap canvas(width, height, rgb(255, 255, 255));
		canvas.line(_t[n][head_t][0], _t[n][head_t][1], _t[n][head_l][0], _t[n][head_l][1], 1, color(Black));	// head
		canvas.line(_t[n][head_l][0], _t[n][head_l][1], _t[n][shdlr_c][0], _t[n][shdlr_c][1], 1, color(Black));		// neck??
		canvas.line(_t[n][shdlr_c][0], _t[n][shdlr_c][1], _t[n][elbow_l][0], _t[n][elbow_l][1], 1, color(Green));	// upper arm left
		canvas.line(_t[n][shdlr_c][0], _t[n][shdlr_c][1], _t[n][elbow_r][0], _t[n][elbow_r][1], 1, color(Blue));	// upper arm right
		canvas.line(_t[n][elbow_l][0], _t[n][elbow_l][1], _t[n][hand_l][0], _t[n][hand_l][1], 1, color(Green));	// lower arm left
		canvas.line(_t[n][elbow_r][0], _t[n][elbow_r][1], _t[n][hand_r][0], _t[n][hand_r][1], 1, color(Blue));	// lower arm right
		canvas.line(_t[n][shdlr_c][0], _t[n][shdlr_c][1], _t[n][waist][0], _t[n][waist][1], 1, color(Black));	// chest??
		canvas.line(_t[n][waist][0], _t[n][waist][1], _t[n][butt][0], _t[n][butt][1], 1, color(Black));	// waist??
		canvas.line(_t[n][butt][0], _t[n][butt][1], _t[n][knee_l][0], _t[n][knee_l][1], 1, color(Green));	// thigh left
		canvas.line(_t[n][butt][0], _t[n][butt][1], _t[n][knee_r][0], _t[n][knee_r][1], 1, color(Blue));	// thigh right
		canvas.line(_t[n][heel_l][0], _t[n][heel_l][1], _t[n][knee_l][0], _t[n][knee_l][1], 1, color(Green));	// shank left
		canvas.line(_t[n][heel_r][0], _t[n][heel_r][1], _t[n][knee_r][0], _t[n][knee_r][1], 1, color(Blue));	// shank right
		canvas.line(_t[n][heel_l][0], _t[n][heel_l][1], _t[n][ttoe_l][0], _t[n][ttoe_l][1], 1, color(Green));	// foot left
		canvas.line(_t[n][heel_r][0], _t[n][heel_r][1], _t[n][ttoe_r][0], _t[n][ttoe_r][1], 1, color(Blue));	// foot right
		canvas.out("Animation\\WalkingMan_" + uint2str(n, 2) + ".bmp");
	}

	void drawTrace() {
		bitmap canvas(width, height, rgb(255, 255, 255));
		pixel col = Black;
		for (int i = 0; i < 15; i++) {
			if (i == hand_r) col = Red;
			else col = Black;
			for (int j = 0; j < 24; j++) {
				if (_t[j][i][0] != NAN && _t[j + 1][i][0] != NAN) canvas.line(_t[j][i][0], _t[j][i][1], _t[j + 1][i][0], _t[j + 1][i][1], 1, col);
				canvas.dot(_t[j][i][0], _t[j][i][1], color(Red));
			}
		}
		canvas.out("Animation\\WalkingMan.bmp");
	}


	enum vecs {
		v_head, v_neck, v_chest, v_waist, v_upperarm_l, v_upperarm_r, v_forearm_l, v_forearm_r,
		v_thigh_l, v_thigh_r, v_shank_l, v_shank_r, v_foot_l, v_foot_r,
	};	// total 14 vecs
	vec3 _v[25][14];

	void initVec() {
		for (int i = 0; i < 25; i++) {
			_v[i][v_head] = vec3(_t[i][head_t][0] - _t[i][head_l][0], 0, _t[i][head_t][1] - _t[i][head_l][1]);
			_v[i][v_neck] = vec3(_t[i][head_l][0] - _t[i][shdlr_c][0], 0, _t[i][head_l][1] - _t[i][shdlr_c][1]);
			_v[i][v_chest] = vec3(_t[i][shdlr_c][0] - _t[i][waist][0], 0, _t[i][shdlr_c][1] - _t[i][waist][1]);
			_v[i][v_waist] = vec3(_t[i][butt][0] - _t[i][waist][0], 0, _t[i][butt][1] - _t[i][waist][1]);
			_v[i][v_upperarm_l] = vec3(_t[i][elbow_l][0] - _t[i][shdlr_c][0], 0, _t[i][elbow_l][1] - _t[i][shdlr_c][1]);
			_v[i][v_upperarm_r] = vec3(_t[i][elbow_r][0] - _t[i][shdlr_c][0], 0, _t[i][elbow_r][1] - _t[i][shdlr_c][1]);
			_v[i][v_forearm_l] = vec3(_t[i][hand_l][0] - _t[i][elbow_l][0], 0, _t[i][hand_l][1] - _t[i][elbow_l][1]);
			_v[i][v_forearm_r] = vec3(_t[i][hand_r][0] - _t[i][elbow_r][0], 0, _t[i][hand_r][1] - _t[i][elbow_r][1]);
			_v[i][v_thigh_l] = vec3(_t[i][knee_l][0] - _t[i][butt][0], 0, _t[i][knee_l][1] - _t[i][butt][1]);
			_v[i][v_thigh_r] = vec3(_t[i][knee_r][0] - _t[i][butt][0], 0, _t[i][knee_r][1] - _t[i][butt][1]);
			_v[i][v_shank_l] = vec3(_t[i][heel_l][0] - _t[i][knee_l][0], 0, _t[i][heel_l][1] - _t[i][knee_l][1]);
			_v[i][v_shank_r] = vec3(_t[i][heel_r][0] - _t[i][knee_r][0], 0, _t[i][heel_r][1] - _t[i][knee_r][1]);
			_v[i][v_foot_l] = vec3(_t[i][ttoe_l][0] - _t[i][heel_l][0], 0, _t[i][ttoe_l][1] - _t[i][heel_l][1]);
			_v[i][v_foot_r] = vec3(_t[i][ttoe_r][0] - _t[i][heel_r][0], 0, _t[i][ttoe_r][1] - _t[i][heel_r][1]);
			for (int j = 0; j < 14; j++) _v[i][j] /= _v[i][j].mod();
		}
	}

	void FourierSeries(vecs c, int n, vec3 *a, vec3 *b) {
		const double T = 1.2, omega = 2 * PI / T;
		double t, dt, nwt; vec3 ft;

		for (int m = 0; m < n; m++) {
			a[m] = b[m] = vec3(0, 0, 0);
			for (int i = 1; i <= 24; i++) {
				t = 0.1 * (i - 0.5), nwt = m * omega * t, dt = 0.1;
				ft = 0.5 * (_v[i - 1][c] + _v[i][c]);
				a[m] += ft * (cos(nwt)*dt);
				b[m] += ft * (sin(nwt)*dt);
			}
			a[m] /= T, b[m] /= T;
		}
	}

	void VisualizeFourierSeries(vecs c, int n) {
		bitmap canvas(240, 100, White);
		auto toCanvas = [](double i, double j, unsigned &x, unsigned &y) {	// 0 < i < 24, -1 < j < 1
			x = i * 10, y = 50 * (j + 1);
		};
		initVec();

		unsigned x, y, x_, y_;
		for (int i = 1; i <= 24; i++) {
			toCanvas(i - 1, _v[i - 1][c].x, x_, y_), toCanvas(i, _v[i][c].x, x, y);
			canvas.line(x_, y_, x, y, 1, LightPink);
			toCanvas(i - 1, _v[i - 1][c].y, x_, y_), toCanvas(i, _v[i][c].y, x, y);
			canvas.line(x_, y_, x, y, 1, LightGreen);
			toCanvas(i - 1, _v[i - 1][c].z, x_, y_), toCanvas(i, _v[i][c].z, x, y);
			canvas.line(x_, y_, x, y, 1, LightBlue);
		}

		vec3 *a = new vec3[n], *b = new vec3[n];
		FourierSeries(c, n, a, b);
		auto FourierEvaluate = [](double t, vec3 *a, vec3 *b, int n) -> vec3 {
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
		for (int i = 0; i <= 240; i++) {
			vec3 r = FourierEvaluate(0.01*i, a, b, n);
			toCanvas(0.1*i, r.x, x, y);
			canvas.dot(i, y, Red);
			toCanvas(0.1*i, r.y, x, y);
			canvas.dot(i, y, Green);
			toCanvas(0.1*i, r.z, x, y);
			canvas.dot(i, y, Blue);
		}


		canvas.out("Animation\\VisualizeFourierSeries.bmp");

		delete[] a, b;
	}
}
