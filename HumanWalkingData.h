#pragma once
#include "bitmap.h"
#include "GlassMan.h"

// Manual statistic data, for interpolation/fitting
namespace WalkingMan {
	// height of this man: 540 units
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
			if (i == waist) col = Red;
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

	void LinearRegression(comps c, unsigned j, double &m, double &b) {
		double sumx2, sumx, n, sumxy, sumy;
		double x, y;
		for (n = 0, sumx2 = sumx = sumxy = sumy = 0; n < 25; n++) {
			x = 0.1*n, y = _t[unsigned(n)][c][j] * (GlassMan_std::height() / 540);
			sumx += x, sumx2 += x * x;
			sumxy += x * y, sumy += y;
		}
		double det = sumx2 * n - sumx * sumx;
		m = (sumxy*n - sumx * sumy) / det, b = (sumx2*sumy - sumxy * sumx) / det;
	}

	void FourierSeries(vecs c, int n, vec3 *a, vec3 *b) {
		initVec();
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
		for (int i = 0; i <= 240; i++) {
			vec3 r = FourierEval(0.01*i, a, b, n);
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

namespace TuringMan {
	// height of this man: 580 units
	enum comps { head_t, head_l, shdlr_l, elbow_l, hand_l, waist, butt_l, knee_l, knee_r, heel_l, heel_r, ttoe_l, ttoe_r };		// total 13 points
	const double _t[20][13][2] = {
		//{ {622,1089}, {619,1020}, {633,973}, {610,888}, {556,825}, {621,868}, {640,812}, {658,693}, {570,696}, {732,563}, {511,571}, {692,543}, {475,585} },
		{ {595,1088}, {594,1023}, {607,984}, {597,891}, {547,814}, {584,862}, {615,819}, {637,697}, {551,692}, {728,570}, {505,559}, {695,544}, {462,562} },
		{ {559,1090}, {558,1030}, {579,983}, {576,897}, {548,808}, {550,886}, {576,819}, {575,684}, {533,697}, {687,604}, {502,558}, {679,576}, {469,560} },
		{ {532,1094}, {536,1031}, {558,985}, {565,901}, {547,813}, {536,889}, {546,821}, {517,688}, {528,697}, {640,614}, {503,560}, {630,589}, {470,561} },
		{ {500,1097}, {504,1035}, {530,991}, {548,911}, {538,820}, {509,891}, {519,844}, {453,700}, {511,703}, {550,599}, {503,558}, {529,575}, {469,557} },
		{ {469,1098}, {473,1034}, {499,993}, {523,918}, {521,823}, {488,893}, {484,844}, {405,711}, {496,701}, {446,570}, {505,558}, {415,556}, {471,559} },
		{ {449,1096}, {456,1034}, {478,992}, {499,912}, {497,817}, {463,890}, {456,835}, {382,708}, {482,698}, {379,557}, {505,554}, {342,557}, {468,559} },
		{ {389,1088}, {390,1022}, {391,977}, {408,889}, {390,798}, {377,875}, {376,825}, {315,680}, {431,685}, {268,541}, {500,565}, {234,541}, {467,558} },
		{ {357,1039}, {352,1029}, {342,985}, {353,888}, {322,799}, {329,880}, {329,817}, {278,684}, {380,683}, {265,537}, {479,583}, {232,539}, {462,563} },
		{ {337,1103}, {331,1027}, {317,991}, {315,891}, {274,807}, {324,881}, {317,825}, {268,681}, {340,688}, {267,539}, {444,592}, {230,536}, {432,571} },
		{ {306,1108}, {305,1035}, {284,998}, {272,902}, {221,825}, {292,888}, {289,835}, {264,688}, {295,701}, {266,537}, {378,592}, {232,538}, {360,570} },
		{ {291,1110}, {283,1037}, {258,1005}, {240,904}, {183,833}, {304,890}, {276,831}, {261,689}, {269,696}, {263,538}, {309,578}, {233,537}, {309,578} },
		{ {281,1112}, {270,1037}, {244,996}, {226,908}, {163,835}, {286,888}, {266,836}, {256,682}, {232,702}, {266,536}, {265,574}, {236,537}, {246,571} },
		{ {270,1113}, {262,1036}, {236,1001}, {222,900}, {161,832}, {281,888}, {256,831}, {254,678}, {223,703}, {266,539}, {232,571}, {229,536}, {211,564} },
		{ {269,1113}, {263,1038}, {230,997}, {224,901}, {166,824}, {278,887}, {254,813}, {249,682}, {211,703}, {269,538}, {233,536}, {215,570}, {192,558} },
		{ {265,1112}, {251,1035}, {227,998}, {223,898}, {181,810}, {273,885}, {252,829}, {242,676}, {210,696}, {265,537}, {214,563}, {233,538}, {193,559} },
		{ {257,1112}, {254,1034}, {224,999}, {230,901}, {194,808}, {266,887}, {205,833}, {239,682}, {207,695}, {269,538}, {215,563}, {228,538}, {191,556} },
		{ {255,1115}, {251,1036}, {218,997}, {230,900}, {197,810}, {265,886}, {244,834}, {248,679}, {208,694}, {267,536}, {215,565}, {230,536}, {190,557} },
		{ {244,1112}, {249,1039}, {212,998}, {225,902}, {192,807}, {262,886}, {250,838}, {242,678}, {208,692}, {266,539}, {216,567}, {237,540}, {193,560} },
		{ {240,1113}, {247,1037}, {213,995}, {222,906}, {179,809}, {256,888}, {242,831}, {246,680}, {206,693}, {265,540}, {215,565}, {229,539}, {193,558} },
		{ {241,1113}, {244,1038}, {216,994}, {220,903}, {172,811}, {252,887}, {237,822}, {247,681}, {209,691}, {265,537}, {213,564}, {231,537}, {193,557} }
	};
	const unsigned width = 1200, height = 657;

}

namespace RunningAttitude {
	// These interpolating data comes from a girl of height 120 units, not a man. ^_^
	// This is not a very good set of data.
	enum comps { head_t, head_l, shoulder_l, shoulder_r, elbow_l, elbow_r, hand_l, hand_r, waist, buttocks, knee_l, knee_r, heel_l, heel_r, tiptoe_l, tiptoe_r };		// total 16 points
	enum vecs { v_head, v_chest, v_chest_side, v_waist, v_upperarm_l, v_upperarm_r, v_forearm_l, v_forearm_r, v_thigh_l, v_thigh_r, v_shank_l, v_shank_r, v_foot_l, v_foot_r };		// total 14 vectors
	const unsigned width = 1200, height = 382;
	vec3 _t[43][16] = {	// each 0.05s
		{ vec3(121,166.5), vec3(118,146.5), vec3(112,137.5), vec3(113,137.5), vec3(110,121.5), vec3(104,122.5), vec3(107,105.5), vec3(122,108.5), vec3(106,115.5), vec3(97,106.5), vec3(120,83.5), vec3(93,82.5), vec3(95,84.5), vec3(77,51.5), vec3(89,72.5), vec3(89,48.5) },
		{ vec3(157,172.5), vec3(157,149.5), vec3(148,145.5), vec3(157,146.5), vec3(133,139.5), vec3(166,130.5), vec3(127,128.5), vec3(186,134.5), vec3(147,124.5), vec3(138,114.5), vec3(172,106.5), vec3(120,82.5), vec3(157,76.5), vec3(77,60.5), vec3(168,72.5), vec3(82,47.5) },
		{ vec3(158,171.5), vec3(157,151.5), vec3(148,146.5), vec3(157,145.5), vec3(133,138.5), vec3(165,129.5), vec3(126,127.5), vec3(184,133.5), vec3(145,123.5), vec3(136,114.5), vec3(172,105.5), vec3(114,82.5), vec3(156,76.5), vec3(77,60.5), vec3(168,72.5), vec3(82,47.5) },
		{ vec3(193,172.5), vec3(191,154.5), vec3(178,147.5), vec3(194,147.5), vec3(169,137.5), vec3(203,131.5), vec3(165,125.5), vec3(221,134.5), vec3(184,123.5), vec3(174,110.5), vec3(192,86.5), vec3(146,88.5), vec3(201,53.5), vec3(114,94.5), vec3(213,55.5), vec3(108,89.5) },
		{ vec3(209,170.5), vec3(210,149.5), vec3(201,143.5), vec3(210,143.5), vec3(193,130.5), vec3(212,126.5), vec3(195,114.5), vec3(228,121.5), vec3(203,118.5), vec3(195,109.5), vec3(209,85.5), vec3(191,78.5), vec3(203,51.5), vec3(169,95.5), vec3(215,50.5), vec3(158,90.5) },
		{ vec3(228,169.5), vec3(229,146.5), vec3(224,141.5), vec3(224,141.5), vec3(217,126.5), vec3(220,123.5), vec3(221,108.5), vec3(226,107.5), vec3(221,117.5), vec3(214,108.5), vec3(215,80.5), vec3(230,80.5), vec3(201,51.5), vec3(199,84.5), vec3(213,49.5), vec3(196,74.5) },
		{ vec3(249,168.5), vec3(247,146.5), vec3(252,141.5), vec3(239,141.5), vec3(256,122.5), vec3(227,128.5), vec3(268,120.5), vec3(225,111.5), vec3(241,119.5), vec3(232,110.5), vec3(222,80.5), vec3(260,98.5), vec3(203,58.5), vec3(235,81.5), vec3(213,49.5), vec3(240,68.5) },
		{ vec3(271,171.5), vec3(265,149.5), vec3(273,145.5), vec3(255,144.5), vec3(278,131.5), vec3(241,139.5), vec3(287,133.5), vec3(235,123.5), vec3(260,123.5), vec3(251,117.5), vec3(227,84.5), vec3(281,104.5), vec3(211,68.5), vec3(264,74.5), vec3(211,56.5), vec3(278,70.5) },
		{ vec3(290,173.5), vec3(285,151.5), vec3(293,148.5), vec3(276,148.5), vec3(301,137.5), vec3(259,143.5), vec3(308,140.5), vec3(249,132.5), vec3(277,127.5), vec3(271,118.5), vec3(251,89.5), vec3(301,105.5), vec3(228,77.5), vec3(295,73.5), vec3(223,66.5), vec3(312,73.5) },
		{ vec3(307,176.5), vec3(303,150.5), vec3(311,150.5), vec3(292,150.5), vec3(323,141.5), vec3(277,145.5), vec3(329,143.5), vec3(265,136.5), vec3(298,125.5), vec3(287,117.5), vec3(271,88.5), vec3(314,98.5), vec3(248,84.5), vec3(324,66.5), vec3(241,74.5), vec3(338,68.5) },
		{ vec3(347,171.5), vec3(342,149.5), vec3(349,146.5), vec3(331,147.5), vec3(358,138.5), vec3(317,140.5), vec3(365,135.5), vec3(310,133.5), vec3(335,123.5), vec3(326,112.5), vec3(326,77.5), vec3(339,82.5), vec3(296,87.5), vec3(345,50.5), vec3(288,80.5), vec3(359,51.5) },
		{ vec3(396,157.5), vec3(387,140.5), vec3(381,134.5), vec3(382,134.5), vec3(365,125.5), vec3(386,116.5), vec3(368,108.5), vec3(404,115.5), vec3(373,115.5), vec3(364,106.5), vec3(393,88.5), vec3(362,81.5), vec3(371,65.5), vec3(343,51.5), vec3(383,56.5), vec3(355,46.5) },
		{ vec3(417,155.5), vec3(409,138.5), vec3(399,136.5), vec3(408,133.5), vec3(381,128.5), vec3(419,123.5), vec3(373,115.5), vec3(434,130.5), vec3(393,115.5), vec3(385,107.5), vec3(415,89.5), vec3(371,80.5), vec3(416,59.5), vec3(344,55.5), vec3(432,58.5), vec3(354,47.5) },
		{ vec3(440,152.5), vec3(433,138.5), vec3(422,138.5), vec3(430,135.5), vec3(401,134.5), vec3(444,125.5), vec3(386,125.5), vec3(457,133.5), vec3(413,113.5), vec3(406,107.5), vec3(428,86.5), vec3(389,78.5), vec3(436,51.5), vec3(357,67.5), vec3(449,54.5), vec3(356,58.5) },
		{ vec3(461,150.5), vec3(452,134.5), vec3(437,132.5), vec3(449,130.5), vec3(420,129.5), vec3(450,112.5), vec3(406,124.5), vec3(465,118.5), vec3(434,110.5), vec3(423,103.5), vec3(445,81.5), vec3(411,74.5), vec3(442,46.5), vec3(384,81.5), vec3(455,48.5), vec3(373,72.5) },
		{ vec3(479,150.5), vec3(469,132.5), vec3(462,129.5), vec3(462,129.5), vec3(444,118.5), vec3(456,111.5), vec3(434,109.5), vec3(465,98.5), vec3(450,106.5), vec3(442,100.5), vec3(461,75.5), vec3(434,67.5), vec3(442,47.5), vec3(406,83.5), vec3(454,43.5), vec3(396,78.5) },
		{ vec3(510,158.5), vec3(504,135.5), vec3(509,127.5), vec3(493,132.5), vec3(517,114.5), vec3(477,122.5), vec3(531,114.5), vec3(465,108.5), vec3(493,110.5), vec3(482,105.5), vec3(466,80.5), vec3(514,94.5), vec3(443,52.5), vec3(484,78.5), vec3(452,42.5), vec3(492,70.5) },
		{ vec3(530,165.5), vec3(523,144.5), vec3(531,139.5), vec3(512,141.5), vec3(540,132.5), vec3(494,138.5), vec3(555,137.5), vec3(482,130.5), vec3(511,121.5), vec3(505,115.5), vec3(477,85.5), vec3(534,109.5), vec3(457,65.5), vec3(520,78.5), vec3(457,54.5), vec3(533,74.5) },
		{ vec3(550,173.5), vec3(542,151.5), vec3(551,148.5), vec3(532,149.5), vec3(563,147.5), vec3(514,149.5), vec3(576,152.5), vec3(504,145.5), vec3(535,128.5), vec3(525,121.5), vec3(504,96), vec3(555,116.5), vec3(475,76.5), vec3(549,80.5), vec3(472,63.5), vec3(564,80.5) },
		{ vec3(567,174.5), vec3(563,155.5), vec3(571,153.5), vec3(551,154.5), vec3(584,154.5), vec3(535,153.5), vec3(596,161.5), vec3(524,152.5), vec3(551,129.5), vec3(543,122.5), vec3(528,100.5), vec3(577,115.5), vec3(492,85.5), vec3(577,80.5), vec3(488,73.5), vec3(593,80.5) },
		{ vec3(589,174.5), vec3(582,154.5), vec3(591,155.5), vec3(571,156.5), vec3(606,161.5), vec3(558,155.5), vec3(616,165.5), vec3(547,152.5), vec3(573,130.5), vec3(565,123.5), vec3(549,97.5), vec3(598,112.5), vec3(512,91.5), vec3(602,78.5), vec3(505,80.5), vec3(617,77.5) },
		{ vec3(610,174.5), vec3(603,152.5), vec3(612,154.5), vec3(590,155.5), vec3(624,162.5), vec3(578,155.5), vec3(639,168.5), vec3(567,150.5), vec3(592,128.5), vec3(581,122.5), vec3(566,95.5), vec3(616,106.5), vec3(532,93.5), vec3(621,72.5), vec3(523,86.5), vec3(636,71.5) },
		{ vec3(628,172.5), vec3(621,151.5), vec3(629,156.5), vec3(612,152.5), vec3(639,162.5), vec3(595,152.5), vec3(655,167.5), vec3(584,143.5), vec3(612,128.5), vec3(601,116.5), vec3(585,92.5), vec3(632,101.5), vec3(554,98.5), vec3(641,67.5), vec3(543,91.5), vec3(657,65.5) },
		{ vec3(667,167.5), vec3(659,147.5), vec3(666,150.5), vec3(648,146.5), vec3(673,153.5), vec3(636,140.5), vec3(680,152.5), vec3(627,124.5), vec3(651,121.5), vec3(643,112.5), vec3(626,86.5), vec3(665,87.5), vec3(598,102.5), vec3(676,52.5), vec3(586,95.5), vec3(686,51.5) },
		{ vec3(686,164.5), vec3(678,143.5), vec3(682,141.5), vec3(669,141.5), vec3(685,139.5), vec3(659,130.5), vec3(691,138.5), vec3(655,113.5), vec3(668,118.5), vec3(659,107.5), vec3(648,82.5), vec3(682,83.5), vec3(619,97.5), vec3(678,60.5), vec3(608,90.5), vec3(684,49.5) },
		{ vec3(699.5,161.5), vec3(695.5,140.5), vec3(688.5,136.5), vec3(688.5,136.5), vec3(696.5,127.5), vec3(681.5,122.5), vec3(700.5,123.5), vec3(682.5,101.5), vec3(687.5,114.5), vec3(679.5,107.5), vec3(676.5,76.5), vec3(695.5,77.5), vec3(647.5,90.5), vec3(681.5,45.5), vec3(638.5,82.5), vec3(688.5,44.5) },
		{ vec3(722.5,161.5), vec3(719.5,138.5), vec3(711.5,136.5), vec3(715.5,133.5), vec3(700.5,129.5), vec3(712.5,114.5), vec3(708.5,117.5), vec3(724.5,101.5), vec3(706.5,110.5), vec3(698.5,105.5), vec3(711.5,70.5), vec3(701.5,71.5), vec3(676.5,79.5), vec3(679.5,47.5), vec3(669.5,73.5), vec3(689.5,45.5) },
		{ vec3(746.5,163.5), vec3(740.5,141.5), vec3(730.5,139.5), vec3(741.5,135.5), vec3(715.5,131.5), vec3(745.5,116.5), vec3(709.5,119.5), vec3(761.5,116.5), vec3(727.5,116.5), vec3(717.5,110.5), vec3(741.5,82.5), vec3(741.5,82.5), vec3(712.5,79.5), vec3(702.5,74.5), vec3(707.5,68.5), vec3(684.5,45.5) },
		{ vec3(768.5,169.5), vec3(763.5,146.5), vec3(752.5,144.5), vec3(765.5,141.5), vec3(731.5,137.5), vec3(772.5,126.5), vec3(722.5,127.5), vec3(789.5,129.5), vec3(748.5,122.5), vec3(738.5,114.5), vec3(770.5,106.5), vec3(719.5,85.5), vec3(752.5,75.5), vec3(692.5,65.5), vec3(763.5,69.5), vec3(694.5,55.5) },
		{ vec3(785.5,170.5), vec3(784.5,150.5), vec3(770.5,148.5), vec3(785,145), vec3(749.5,140.5), vec3(795.5,134.5), vec3(739.5,133.5), vec3(812.5,140.5), vec3(769.5,126.5), vec3(763.5,118.5), vec3(791.5,106.5), vec3(743.5,88.5), vec3(789.5,72.5), vec3(710.5,79.5), vec3(800.5,74.5), vec3(703.5,67.5) },
		{ vec3(819.5,170.5), vec3(820.5,148.5), vec3(805.5,149.5), vec3(822.5,146.5), vec3(780.5,131.5), vec3(843.5,142.5), vec3(793.5,129.5), vec3(857.5,148.5), vec3(807.5,122.5), vec3(799.5,115.5), vec3(823.5,92.5), vec3(785.5,82.5), vec3(836.5,58.5), vec3(765.5,99.5), vec3(853.5,54.5), vec3(743.5,91.5) },
		{ vec3(838.5,167.5), vec3(840.5,146.5), vec3(827.5,145.5), vec3(841.5,145.5), vec3(810.5,138.5), vec3(859.5,137.5), vec3(819.5,126.5), vec3(878.5,142.5), vec3(827.5,118.5), vec3(820.5,107.5), vec3(841.5,82.5), vec3(804.5,79.5), vec3(848.5,45.5), vec3(777.5,97.5), vec3(856.5,43.5), vec3(763.5,93.5) },
		{ vec3(883.5,162.5), vec3(879.5,140.5), vec3(876.5,137.5), vec3(878.5,134.5), vec3(862.5,125.5), vec3(875.5,117.5), vec3(868.5,117.5), vec3(886.5,106.5), vec3(869.5,112.5), vec3(863.5,105.5), vec3(866.5,73.5), vec3(879.5,70.5), vec3(848.5,44.5), vec3(845.5,81.5), vec3(862.5,42.5), vec3(834.5,81.5) },
		{ vec3(909.5,166.5), vec3(901.5,142.5), vec3(906.5,138.5), vec3(891.5,138.5), vec3(912.5,125.5), vec3(882.5,125.5), vec3(924.5,116.5), vec3(882.5,110.5), vec3(890.5,113.5), vec3(877.5,104.5), vec3(870.5,79.5), vec3(915.5,98.5), vec3(850.5,54.5), vec3(888.5,78.5), vec3(862.5,43.5), vec3(897.5,65.5) },
		{ vec3(933.5,169.5), vec3(926.5,146.5), vec3(931.5,143.5), vec3(910.5,146.5), vec3(943.5,137.5), vec3(894.5,138.5), vec3(956.5,138.5), vec3(886.5,126.5), vec3(909.5,121.5), vec3(902.5,114.5), vec3(881.5,86.5), vec3(937.5,105.5), vec3(863.5,63.5), vec3(923.5,75.5), vec3(863.5,51.5), vec3(940.5,68.5) },
		{ vec3(954.5,169.5), vec3(944.5,148.5), vec3(953.5,146.5), vec3(932.5,148.5), vec3(966.5,142.5), vec3(917.5,146.5), vec3(976.5,145.5), vec3(903.5,137.5), vec3(933.5,125.5), vec3(922.5,116.5), vec3(907.5,90.5), vec3(959.5,106.5), vec3(878.5,76.5), vec3(951.5,75.5), vec3(876.5,63.5), vec3(968.5,74.5) },
		{ vec3(994.5,168.5), vec3(990.5,145.5), vec3(995.5,142.5), vec3(977.5,149.5), vec3(1002.5,131.5), vec3(959.5,144.5), vec3(1013.5,137.5), vec3(949.5,135.5), vec3(972.5,124.5), vec3(963.5,117.5), vec3(952.5,88.5), vec3(983.5,89.5), vec3(921.5,91.5), vec3(985.5,52.5), vec3(909.5,83.5), vec3(1002.5,51.5) },
		{ vec3(1021.5,165.5), vec3(1011.5,144.5), vec3(1010.5,138.5), vec3(1001.5,140.5), vec3(1007.5,119.5), vec3(986.5,131.5), vec3(1023.5,117.5), vec3(986.5,115.5), vec3(992.5,118.5), vec3(986.5,112.5), vec3(963.5,80.5), vec3(1001.5,82.5), vec3(954.5,97.5), vec3(986.5,50.5), vec3(939.5,90.5), vec3(1003.5,48.5) },
		{ vec3(1039.5,165.5), vec3(1031.5,142.5), vec3(1028.5,138.5), vec3(1028.5,137.5), vec3(1010.5,127.5), vec3(1022.5,121.5), vec3(1018.5,116.5), vec3(1034.5,106.5), vec3(1017.5,116.5), vec3(1009.5,110.5), vec3(1019.5,78.5), vec3(1003.5,80.5), vec3(997.5,86.5), vec3(986.5,52.5), vec3(981.5,80.5), vec3(999.5,48.5) },
		{ vec3(1058.5,166.5), vec3(1054.5,144.5), vec3(1044.5,143.5), vec3(1057.5,139.5), vec3(1024.5,133.5), vec3(1061.5,119.5), vec3(1015.5,117.5), vec3(1078.5,119.5), vec3(1039.5,120.5), vec3(1029.5,112.5), vec3(1058.5,101.5), vec3(1014.5,87.5), vec3(1026.5,87.5), vec3(1009.5,83.5), vec3(1035.5,76.5), vec3(1000.5,49.5) },
		{ vec3(1080.5,169.5), vec3(1080.5,148.5), vec3(1063.5,147.5), vec3(1082.5,143.5), vec3(1045.5,139.5), vec3(1091.5,127.5), vec3(1032.5,128.5), vec3(1104.5,131.5), vec3(1063.5,126.5), vec3(1055.5,118.5), vec3(1085.5,109.5), vec3(1034.5,95.5), vec3(1069.5,76.5), vec3(999.5,68.5), vec3(1085.5,72.5), vec3(1001.5,56.5) },
		{ vec3(1104.5,171.5), vec3(1102.5,148.5), vec3(1085.5,147.5), vec3(1105.5,144.5), vec3(1068.5,138.5), vec3(1116.5,130.5), vec3(1059.5,123.5), vec3(1128.5,135.5), vec3(1087.5,123.5), vec3(1078.5,116.5), vec3(1104.5,101.5), vec3(1030.5,90.5), vec3(1104.5,66.5), vec3(1027.5,85.5), vec3(1122.5,64.5), vec3(1020.5,75.5) },
		{ vec3(1128.5,169.5), vec3(1125.5,148.5), vec3(1112.5,144.5), vec3(1127.5,143.5), vec3(1095.5,133.5), vec3(1131.5,125.5), vec3(1093.5,116.5), vec3(1146.5,124.5), vec3(1109.5,119.5), vec3(1099.5,111.5), vec3(1116.5,86.5), vec3(1091.5,80.5), vec3(1122.5,48.5), vec3(1059.5,94.5), vec3(1136.5,51.5), vec3(1049.5,86.5) },
	};
	const double speed = (_t[42][waist].x - _t[0][waist].x) / (0.05 * 43);
	vec3 _v[43][14];
	enum lengths { upperarm_l, upperarm_r, forearm_l, forearm_r, shank_l, shank_r };
	double _l[44][6], *l = _l[43]; // 0-42: data; 43:max;

	void visualize() {
		bitmap img(width, height, color(White));
		for (int i = 0; i < 16; i++) {
			for (int j = 1; j <= 42; j++) {
				img.line(_t[j - 1][i].x, _t[j - 1][i].y, _t[j][i].x, _t[j][i].y, 1, rgb(0, 0, 0));
			}
		}
		img.out("IMAGE\\VSL.bmp");
	}


	// Fourier: 0-27 forms a period
	const unsigned N = 8;
	const double T = 27 * 0.05, omega = 2 * PI / T;
	vec3 FourierCoffs_a[16][N], FourierCoffs_b[16][N];
	void calcFourier(vecs v) {
		vec3 a, b;
		double t, dt, nwt; vec3 ft;
		for (int n = 0; n < N; n++) {
			a = b = vec3(0, 0, 0);
			for (int i = 1; i <= 27; i++) {
				t = 0.05 * (i - 0.5), nwt = n * omega * t, dt = 0.05;
				ft = 0.5 * (_v[i - 1][v] + _v[i][v]);
				a += ft * (cos(nwt)*dt);
				b += ft * (sin(nwt)*dt);
			}
			a /= T, b /= T;
			if (n != 0) a *= 2, b *= 2;
			FourierCoffs_a[v][n] = a, FourierCoffs_b[v][n] = b;
		}
	}
	vec3 FourierEval(vecs v, double t) {
		vec3 r(0, 0, 0);
		for (int n = 0; n < N; n++) {
			double nwt = n * omega*t;
			r += FourierCoffs_a[v][n] * cos(nwt) + FourierCoffs_b[v][n] * sin(nwt);
		}
		return r / r.mod();
	}

	point WaistFourierCoffs_a[N], WaistFourierCoffs_b[N];	// derivative of waist
	void calcFourier() {
		point a, b;
		double t, dt, nwt; point ft;
		for (int n = 0; n < N; n++) {
			a = b = point(0, 0, 0);
			for (int i = 1; i <= 27; i++) {
				t = 0.05 * (i - 0.5), nwt = n * omega * t, dt = 0.05;
				ft = 0.5 * ((_t[i][waist] - _t[i - 1][waist]) / 0.05 + (_t[i + 1][waist] - _t[i][waist]) / 0.05);
				a += ft * (cos(nwt)*dt);
				b += ft * (sin(nwt)*dt);
			}
			a /= T, b /= T;
			if (n != 0) a *= 2, b *= 2;
			WaistFourierCoffs_a[n] = a, WaistFourierCoffs_b[n] = b;
		}
		ft = point(0, 0, 0); for (int i = 0; i < 27; i++) ft += _t[i][waist]; ft /= 27;
		WaistFourierCoffs_b[0] = ft;
	}
	point FourierEval(double t) {	// integral
		//point r(_t[0][waist]);
		point r(0, 0, WaistFourierCoffs_b[0].z - 45);
		r += WaistFourierCoffs_a[0] * t;
		for (int n = 1; n < N; n++) {
			double nw = n * omega, nwt = nw * t;
			r += (WaistFourierCoffs_a[n] * sin(nwt) - WaistFourierCoffs_b[n] * cos(nwt)) / nw;
		}
		return r;
	}

	void init() {
		if (_t[0][0].z == 0) for (int i = 0; i < 43; i++) {
			for (int j = 0; j < 16; j++) {
				_t[i][j].z = _t[i][j].y, _t[i][j].y = _t[i][j].x, _t[i][j].x = 0;
			}
		}

		for (int i = 0; i < 43; i++) {
			_l[i][upperarm_l] = (_t[i][shoulder_l] - _t[i][elbow_l]).mod(), _l[i][upperarm_r] = (_t[i][shoulder_r] - _t[i][elbow_r]).mod();
			_l[i][forearm_l] = (_t[i][hand_l] - _t[i][elbow_l]).mod(), _l[i][forearm_r] = (_t[i][hand_r] - _t[i][elbow_r]).mod();
			_l[i][shank_l] = (_t[i][heel_l] - _t[i][knee_l]).mod(), _l[i][shank_r] = (_t[i][heel_r] - _t[i][knee_r]).mod();
			//for (int j = 0; j < 6; j++) _l[i][j] *= 0.6;
		}
		for (int i = 0; i < 6; i++) l[i] = 0;
		for (int i = 0; i < 43; i++) {
			for (int j = 0; j < 6; j++) {
				if (_l[i][j] > l[j]) l[j] = _l[i][j];
			}
		}

		for (int i = 0; i < 43; i++) {
			vec3 *v = _v[i], *t = _t[i];
			v[v_head] = t[head_t] - t[head_l];
			v[v_chest] = 0.5*(t[shoulder_l] + t[shoulder_r]) - t[waist];
			v[v_chest_side] = t[shoulder_l] * 1.000001 - t[shoulder_r]; v[v_chest_side].x = -50;
			v[v_waist] = t[buttocks] - t[waist];
			v[v_upperarm_l] = t[elbow_l] - t[shoulder_l], v[v_forearm_l] = t[hand_l] - t[elbow_l];
			v[v_upperarm_r] = t[elbow_r] - t[shoulder_r], v[v_forearm_r] = t[hand_r] - t[elbow_r];
			v[v_thigh_l] = t[knee_l] - t[buttocks], v[v_shank_l] = t[heel_l] - t[knee_l], v[v_foot_l] = t[tiptoe_l] - t[heel_l];
			v[v_thigh_r] = t[knee_r] - t[buttocks], v[v_shank_r] = t[heel_r] - t[knee_r], v[v_foot_r] = t[tiptoe_r] - t[heel_r];

			const double tm = 0.6;
			v[v_upperarm_l].x = -tm * sqrt(l[upperarm_l] * l[upperarm_l] - dot(v[v_upperarm_l], v[v_upperarm_l]));
			v[v_upperarm_r].x = tm * sqrt(l[upperarm_r] * l[upperarm_r] - dot(v[v_upperarm_r], v[v_upperarm_r]));
			v[v_forearm_l].x = -tm * sqrt(l[forearm_l] * l[forearm_l] - dot(v[v_forearm_l], v[v_forearm_l]));	// not alway negative
			v[v_forearm_r].x = tm * sqrt(l[forearm_r] * l[forearm_r] - dot(v[v_forearm_r], v[v_forearm_r]));
			//v[v_shank_l].x = -tm * sqrt(l[shank_l] * l[shank_l] - dot(v[v_shank_l], v[v_shank_l]));
			//v[v_shank_r].x = tm * sqrt(l[shank_r] * l[shank_r] - dot(v[v_shank_r], v[v_shank_r]));
			if (isnan(v[v_upperarm_l].x)) v[v_upperarm_l].x = -1; if (isnan(v[v_upperarm_r].x)) v[v_upperarm_r].x = 1;
			if (isnan(v[v_forearm_l].x)) v[v_forearm_l].x = -1; if (isnan(v[v_forearm_r].x)) v[v_forearm_r].x = 1;
			//if (isnan(v[v_shank_l].x)) v[v_shank_l].x = -1; if (isnan(v[v_shank_r].x)) v[v_shank_r].x = 1;

			for (int j = 0; j < 14; j++) v[j] /= v[j].mod();
		}

		for (int i = 0; i < 14; i++) calcFourier(vecs(i));
		calcFourier();
	}

	void VisualizeFourier(vecs v) {
		init();
		bitmap img(1000, 200);
		img.line(0, 100, 1000, 100, 1, color(Gray));
		for (int i = 1; i < 43; i++) {
			img.line(20 * (i - 1), 100 * (_v[i - 1][v].x + 1), 20 * i, 100 * (_v[i][v].x + 1), 1, color(LightPink));
			img.line(20 * (i - 1), 100 * (_v[i - 1][v].y + 1), 20 * i, 100 * (_v[i][v].y + 1), 1, color(LightGreen));
			img.line(20 * (i - 1), 100 * (_v[i - 1][v].z + 1), 20 * i, 100 * (_v[i][v].z + 1), 1, color(LightBlue));
		}

		calcFourier(v);
		vec3 r, _r; double t;
		for (int i = 0; i < 1000; i++) {
			t = i * 0.05 / 20;
			_r = r, r = FourierEval(v, t); r.x++, r.y++, r.z++;
			img.line(i - 1, 100 * _r.x, i, 100 * r.x, 1, color(Red));
			img.line(i - 1, 100 * _r.y, i, 100 * r.y, 1, color(Green));
			img.line(i - 1, 100 * _r.z, i, 100 * r.z, 1, color(Blue));
		}

		img.out("IMAGE\\VSP.bmp");
	}
	void VisualizeFourier() {
		init();
		bitmap img(width, height, color(White));

		for (int j = 1; j <= 42; j++) {
			img.line(_t[j - 1][waist].y, _t[j - 1][waist].z, _t[j][waist].y, _t[j][waist].z, 1, color(LightGray));
		}

		point r, _r; double t;
		for (int i = -100; i < 200; i++) {
			t = (i / 100.)*(27 * 0.05);
			_r = r, r = FourierEval(t);
			img.line(_r.y, _r.z, r.y, r.z, 1, color(Black));
		}

		img.out("IMAGE\\VSW.bmp");
	}

	void visualize(comps t) {
		init();
		bitmap img(width, height, color(White));
		for (int j = 1; j < 43; j++) {
			img.line(_t[j - 1][t].y, _t[j - 1][t].z, _t[j][t].y, _t[j][t].z, 1, rgb(0, 0, 0));
		}
		img.out("IMAGE\\VSL.bmp");
	}
	void visualize(comps t1, comps t2) {
		init();
		bitmap img(width, height, color(White));
		for (int j = 1; j < 43; j++) {
			img.line(_t[j - 1][t1].y, _t[j - 1][t1].z, _t[j][t1].y, _t[j][t1].z, 1, color(Red));
			img.line(_t[j - 1][t2].y, _t[j - 1][t2].z, _t[j][t2].y, _t[j][t2].z, 1, color(Blue));
			img.line(_t[j][t1].y, _t[j][t1].z, _t[j][t2].y, _t[j][t2].z, 1, color(Gray));
		}
		img.out("IMAGE\\VSL.bmp");
	}
	void visualize(lengths t) {
		init();
		bitmap img(860, 200);
		for (int i = 1; i < 43; i++) {
			img.line(20 * (i - 1), 200 * (_l[i - 1][t] / l[t]), 20 * i, 200 * (_l[i][t] / l[t]), 1, color(Black));
		}
		img.out("IMAGE\\VSL.bmp");
	}
	void visualize(vecs v) {
		init();
		bitmap img(860, 200);
		img.line(0, 100, 840, 100, 1, color(Gray));
		for (int i = 1; i < 43; i++) {
			img.line(20 * (i - 1), 100 * (_v[i - 1][v].x + 1), 20 * i, 100 * (_v[i][v].x + 1), 1, color(Red));
			img.line(20 * (i - 1), 100 * (_v[i - 1][v].y + 1), 20 * i, 100 * (_v[i][v].y + 1), 1, color(Green));
			img.line(20 * (i - 1), 100 * (_v[i - 1][v].z + 1), 20 * i, 100 * (_v[i][v].z + 1), 1, color(Blue));
		}
		img.out("IMAGE\\VSP.bmp");
	}


}
