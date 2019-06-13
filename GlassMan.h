#pragma once

#include "World.h"

class GlassMan {

	XSolid Body;

	XObjs::Box_affine _HEAD, _CHEST, _WAIST;
	XObjs::Cylinder _NECK;
	XObjs::Box_affine _SHOULDER_L, _SHOULDER_R;
	XObjs::Cone_trunc _UPPER_ARM_L, _UPPER_ARM_R, _FOREARM_L, _FOREARM_R;
	XObjs::Sphere _ELBOW_L, _ELBOW_R, _HAND_L, _HAND_R;
	XObjs::Box_affine _THIGH_L, _THIGH_R;
	XObjs::Cone_trunc _SHANK_L, _SHANK_R;
	XObjs::Box_affine _FOOT_L, _FOOT_R;

	friend class GlassMan_std;

public:
	double damping_r, damping_g, damping_b;
	double refractive_index;

	/* Fill all of these parameters legally before constructing */
	point Pos; 	// position the glass man standing/sitting at
	vec3 top;	// direction the glass man's body towards (up)
	vec3 dir;	// direction the glass man facing

	// The following points/vectors are relative to Pos (default 0,0,0), top (default 0,0,1) and dir (default 1,0,0)
	point Heel_l, Heel_r;		// position of heel/ankle
	vec3 foot_dir_l, foot_dir_r, foot_side_l, foot_side_r;	// direction the foot pointing and a vector to side of foot
	double width_of_heel, height_of_heel, width_of_tiptoe, length_of_foot, foot_rounding;	// width and height of heels and tiptoes; length of feet; rounding radius
	point Knee_l, Knee_r;	// position of knees
	double upper_radius_of_shank, lower_radius_of_shank, shank_rounding;		// thickness of upper and lower end of shanks; rounding radius
	point Buttock_l, Buttock_r;
	double upper_radius_of_thigh, lower_radius_of_thigh, thigh_rounding;
	point Waist;
	double lower_width_of_waist, lower_depth_of_waist, width_of_waist, waist_rounding;
	point Shoulder_l, Shoulder_r;
	double depth_of_chest, chest_rounding;
	point Elbow_l, Elbow_r;
	double upper_radius_of_upperarm, lower_radius_of_upperarm, upperarm_rounding;
	point Hand_l, Hand_r;
	double upper_radius_of_forearm, lower_radius_of_forearm, forearm_rounding;
	point Head_lower, Head_top; vec3 Head_side;		// extreme low and top of head (not chin and forehead); any vector to the side of head
	double width_of_forehead, width_of_chin, depth_of_forehead, head_rounding; double neck_radius, neck_rounding;


	GlassMan() {
		damping_r = damping_g = damping_b = 0;
		refractive_index = 1.5;
		Pos = point(0, 0, 0), top = point(0, 0, 1), dir = point(1, 0, 0);
	}
	GlassMan(const GlassMan &other) {
		const unsigned n = sizeof(GlassMan);
		byte* p = (byte*)this;
		const byte* q = (byte*)&other;
		for (unsigned i = 0; i < n; i++) {
			*p = *q, p++, q++;
		}
		Body = XSolid();
	}
	~GlassMan() {}

	XSolid construct() {
		top /= top.mod(), dir /= dir.mod();
		correct_vector(top, dir);
		double rx, ry, rz; getRotationAngle_ik(dir, top, rx, ry, rz);

		Body = XSolid();

		// constructing feet
		XSolid tmp_FEET; {
			_FOOT_L = XObjs::Box_affine();
			foot_dir_l /= foot_dir_l.mod(), foot_dir_r /= foot_dir_r.mod(),
				foot_side_l /= foot_side_l.mod(), foot_side_r /= foot_side_r.mod();
			_FOOT_L.translate(0, -0.5, 0);
			double foot_h_w = width_of_heel - 2 * foot_rounding, foot_h_h = height_of_heel - 2 * foot_rounding,
				foot_t_w = width_of_tiptoe - 2 * foot_rounding, foot_l = length_of_foot - 2 * foot_rounding;
			_FOOT_L.scale(foot_l * foot_h_w / foot_t_w, foot_h_w, foot_h_h);
			_FOOT_L.perspective(foot_h_w / foot_t_w - 1, 0, 0);
			_FOOT_R = _FOOT_L;
			double foot_rx, foot_ry, foot_rz;
			correct_vector(foot_dir_l, foot_side_l); getRotationAngle_ij(foot_dir_l, foot_side_l, foot_rx, foot_ry, foot_rz);
			_FOOT_L.rotate(foot_rx, foot_ry, foot_rz);
			correct_vector(foot_dir_r, foot_side_r); getRotationAngle_ij(foot_dir_r, -foot_side_r, foot_rx, foot_ry, foot_rz);
			_FOOT_R.rotate(foot_rx, foot_ry, foot_rz);
			_FOOT_L += Heel_l, _FOOT_R += Heel_r;
			XSolid tmp_FOOT_L = _FOOT_L, tmp_FOOT_R = _FOOT_R;
			tmp_FOOT_L = CSG_RoundingOp(tmp_FOOT_L, foot_rounding), tmp_FOOT_R = CSG_RoundingOp(tmp_FOOT_R, foot_rounding);
			tmp_FEET = CSG_UnionOp(tmp_FOOT_L, tmp_FOOT_R);
		}

		// constructing shanks
		XSolid tmp_SHANKS; {
			double shank_up_r = upper_radius_of_shank - shank_rounding, shank_lw_r = lower_radius_of_shank - shank_rounding;
			point shank_root = Heel_l + shank_lw_r * foot_dir_l, shank_dir = shank_root - Knee_l; shank_dir /= shank_dir.mod();
			_SHANK_L = XObjs::Cone_trunc(Knee_l + shank_rounding * shank_dir, shank_root - (shank_rounding + 0.3 * height_of_heel) * shank_dir, shank_up_r, shank_lw_r);
			shank_root = Heel_r + shank_lw_r * foot_dir_r; shank_dir = shank_root - Knee_r, shank_dir /= shank_dir.mod();
			_SHANK_R = XObjs::Cone_trunc(Knee_r + shank_rounding * shank_dir, shank_root - (shank_rounding + 0.3 * height_of_heel) * shank_dir, shank_up_r, shank_lw_r);
			XSolid tmp_SHANK_L = _SHANK_L, tmp_SHANK_R = _SHANK_R;
			tmp_SHANK_L = CSG_RoundingOp(tmp_SHANK_L, shank_rounding), tmp_SHANK_R = CSG_RoundingOp(tmp_SHANK_R, shank_rounding);
			tmp_SHANKS = CSG_UnionOp(tmp_SHANK_L, tmp_SHANK_R);
		}

		// constructing thighs
		XSolid tmp_THIGHS; {
			double thigh_R = 2 * (upper_radius_of_thigh - thigh_rounding), thigh_r = 2 * (lower_radius_of_thigh - thigh_rounding);
			vec3 thigh_dir = Knee_l - Buttock_l, thigh_outer = Buttock_l - Buttock_r; thigh_outer /= thigh_outer.mod();
			double thigh_length = thigh_dir.mod(); thigh_dir /= thigh_length; thigh_length -= 2 * thigh_rounding;
			thigh_outer -= dot(thigh_outer, thigh_dir) * thigh_dir; thigh_outer /= thigh_outer.mod();
			_THIGH_L = XObjs::Box_affine();
			_THIGH_L.transform(matrix3D_affine(thigh_length*thigh_R / thigh_r, 0, 0, 0, 0, thigh_R, 0, 0, 0, 0, thigh_R, 0, thigh_R / thigh_r - 1, 0, 0, 1));
			_THIGH_L.scale(1, -1, 1);
			double thigh_rz = atan2(thigh_dir.y, thigh_dir.x), thigh_ry = atan2(-thigh_dir.z, hypot(thigh_dir.x, thigh_dir.y));
			vec3 thigh_front = cross(thigh_dir, thigh_outer); thigh_front /= thigh_front.mod(); double thigh_rx = atan2(-thigh_front.z, thigh_outer.z);
			_THIGH_L.rotate(thigh_rx, thigh_ry, thigh_rz); _THIGH_L += Buttock_l + thigh_rounding * thigh_dir - 0.5*thigh_r*thigh_outer - 0.5*thigh_r*thigh_front;
			XSolid tmp_THIGH_L = _THIGH_L; tmp_THIGH_L = CSG_RoundingOp(tmp_THIGH_L, thigh_rounding);
			thigh_dir = Knee_r - Buttock_r, thigh_length = thigh_dir.mod(), thigh_dir /= thigh_length; thigh_length -= 2 * thigh_rounding;
			thigh_outer = Buttock_r - Buttock_l; thigh_outer /= thigh_outer.mod(); thigh_outer -= dot(thigh_outer, thigh_dir) * thigh_dir; thigh_outer /= thigh_outer.mod();
			_THIGH_R = XObjs::Box_affine();
			_THIGH_R.transform(matrix3D_affine(thigh_length*thigh_R / thigh_r, 0, 0, 0, 0, thigh_R, 0, 0, 0, 0, thigh_R, 0, thigh_R / thigh_r - 1, 0, 0, 1));
			thigh_rz = atan2(thigh_dir.y, thigh_dir.x), thigh_ry = atan2(-thigh_dir.z, hypot(thigh_dir.x, thigh_dir.y));
			thigh_front = cross(thigh_outer, thigh_dir); thigh_front /= thigh_front.mod(); thigh_rx = atan2(thigh_front.z, thigh_outer.z);
			_THIGH_R.rotate(thigh_rx, thigh_ry, thigh_rz); _THIGH_R += Buttock_r + thigh_rounding * thigh_dir - 0.5*thigh_r*thigh_outer - 0.5*thigh_r*thigh_front;
			XSolid tmp_THIGH_R = _THIGH_R; tmp_THIGH_R = CSG_RoundingOp(tmp_THIGH_R, thigh_rounding);
			tmp_THIGHS = CSG_UnionOp(tmp_THIGH_L, tmp_THIGH_R);
		}

		// constructing waist
		XSolid tmp_WAIST; {
			vec3 waist_dir = Waist - 0.5*(Buttock_l + Buttock_r);
			double waist_h = waist_dir.mod(); waist_dir /= waist_h; waist_h -= waist_rounding/* + waist_rounding*/;
			double waist_lw = lower_width_of_waist - 2 * waist_rounding, waist_ld = lower_depth_of_waist - 2 * waist_rounding, waist_uw = width_of_waist - 2 * waist_rounding;
			_WAIST = XObjs::Box_affine(); _WAIST += point(-0.5, -0.5);
			_WAIST.perspective(0, 0, waist_lw / waist_uw - 1);
			_WAIST.scale(waist_ld, waist_lw, waist_h * waist_lw / waist_uw);
			vec3 waist_p = Buttock_l - Buttock_r; waist_p /= waist_p.mod();
			double waist_rx, waist_ry, waist_rz; correct_vector(waist_dir, waist_p); getRotationAngle_jk(waist_p, waist_dir, waist_rx, waist_ry, waist_rz);
			_WAIST.rotate(_WAIST.center(), waist_rx, waist_ry, waist_rz);
			_WAIST += Waist - (waist_h/* + waist_rounding*/)* waist_dir;
			tmp_WAIST = _WAIST;
			tmp_WAIST = CSG_RoundingOp(tmp_WAIST, waist_rounding);
		}

		// constructing chest
		XSolid tmp_CHEST; {
			point Chest = 0.5 * (Shoulder_l + Shoulder_r);
			vec3 chest_dir = Chest - Waist;
			double chest_h = chest_dir.mod(); chest_dir /= chest_h; chest_h -= chest_rounding;
			vec3 chest_p = Shoulder_l - Shoulder_r; double chest_w = chest_p.mod(); chest_p /= chest_w; chest_w -= 2 * (chest_rounding + upper_radius_of_upperarm);
			double chest_d = depth_of_chest - 2 * chest_rounding, chest_lw = width_of_waist - 2 * chest_rounding;
			_CHEST = XObjs::Box_affine(); _CHEST += point(-0.5, -0.5);
			_CHEST.perspective(0, 0, chest_lw / chest_w - 1);
			_CHEST.scale(chest_d * chest_lw / chest_w, chest_lw, chest_h * chest_lw / chest_w);
			double chest_rx, chest_ry, chest_rz; correct_vector(chest_dir, chest_p); getRotationAngle_jk(chest_p, chest_dir, chest_rx, chest_ry, chest_rz);
			_CHEST.rotate(_CHEST.center(), chest_rx, chest_ry, chest_rz);
			_CHEST += Chest - chest_h * chest_dir;
			tmp_CHEST = _CHEST;
			tmp_CHEST = CSG_RoundingOp(tmp_CHEST, chest_rounding);
		}
		//XSolid tmp_TRUNC = CSG_UnionOp_Smooth(tmp_WAIST, tmp_CHEST, 0.5 * (chest_rounding + waist_rounding));
		XSolid tmp_TRUNC = CSG_UnionOp(tmp_WAIST, tmp_CHEST);

		// constructing arms
		XSolid tmp_UPPERARM; {
			vec3 dir = Elbow_l - Shoulder_l; dir /= dir.mod();
			_UPPER_ARM_L = XObjs::Cone_trunc(Shoulder_l, Elbow_l - upperarm_rounding * dir, upper_radius_of_upperarm - upperarm_rounding, lower_radius_of_upperarm - upperarm_rounding);
			XSolid tmp_UPPERARM_L = CSG_RoundingOp(XSolid(_UPPER_ARM_L), upperarm_rounding);
			dir = Elbow_r - Shoulder_r; dir /= dir.mod();
			_UPPER_ARM_R = XObjs::Cone_trunc(Shoulder_r, Elbow_r - upperarm_rounding * dir, upper_radius_of_upperarm - upperarm_rounding, lower_radius_of_upperarm - upperarm_rounding);
			XSolid tmp_UPPERARM_R = CSG_RoundingOp(XSolid(_UPPER_ARM_R), upperarm_rounding);
			tmp_UPPERARM = CSG_UnionOp(tmp_UPPERARM_L, tmp_UPPERARM_R);
		}
		XSolid tmp_FOREARM; {
			vec3 dir = Hand_l - Elbow_l; dir /= dir.mod();
			_FOREARM_L = XObjs::Cone_trunc(Elbow_l, Hand_l - forearm_rounding * dir, upper_radius_of_forearm - forearm_rounding, lower_radius_of_forearm - forearm_rounding);
			XSolid tmp_FOREARM_L = CSG_RoundingOp(XSolid(_FOREARM_L), forearm_rounding);
			dir = Hand_r - Elbow_r; dir /= dir.mod();
			_FOREARM_R = XObjs::Cone_trunc(Elbow_r, Hand_r - forearm_rounding * dir, upper_radius_of_forearm - forearm_rounding, lower_radius_of_forearm - forearm_rounding);
			XSolid tmp_FOREARM_R = CSG_RoundingOp(XSolid(_FOREARM_R), forearm_rounding);
			tmp_FOREARM = CSG_UnionOp(tmp_FOREARM_L, tmp_FOREARM_R);
		}

		// constructing head
		XSolid tmp_HEAD; {
			vec3 head_dir = Head_top - Head_lower; double head_h = head_dir.mod(); head_dir /= head_h; head_h -= 2 * head_rounding;
			Head_side /= Head_side.mod(); correct_vector(head_dir, Head_side);
			double head_rx, head_ry, head_rz; getRotationAngle_jk(Head_side, head_dir, head_rx, head_ry, head_rz);
			double head_w = width_of_forehead - 2 * head_rounding, head_d = depth_of_forehead - 2 * head_rounding, head_lw = width_of_chin - 2 * head_rounding;
			_HEAD = XObjs::Box_affine(); _HEAD -= point(0.5, 0.5, 0);
			_HEAD.perspective(0, 0, head_lw / head_w - 1);
			_HEAD.scale(head_d * head_lw / head_w, head_lw, head_h * head_lw / head_w);
			_HEAD.rotate(head_rx, head_ry, head_rz);
			_HEAD += Head_lower + head_rounding * head_dir;
			tmp_HEAD = CSG_RoundingOp(XSolid(_HEAD), head_rounding);
		}
		XSolid tmp_NECK; {
			vec3 neck_dir = 0.5*(Shoulder_l + Shoulder_r) - Head_lower; double neck_h = neck_dir.mod(); neck_dir /= neck_h; neck_h -= 2 * neck_rounding + chest_rounding;
			_NECK = XObjs::Cylinder(Head_lower + neck_rounding * neck_dir, neck_radius, neck_h, neck_dir);
			tmp_NECK = CSG_RoundingOp(XSolid(_NECK), neck_rounding);
		}


		Body = CSG_UnionOp(tmp_FEET, tmp_SHANKS);
		Body = CSG_UnionOp(Body, tmp_THIGHS);
		Body = CSG_UnionOp(Body, tmp_TRUNC);
		Body = CSG_UnionOp(Body, CSG_UnionOp(tmp_UPPERARM, tmp_FOREARM));
		Body = CSG_UnionOp(Body, tmp_HEAD);
		Body = CSG_UnionOp(Body, tmp_NECK);

		Body = CSG_Rotation(Body, rx, ry, rz);
		Body = CSG_Translation(Body, Pos);

		Body.col = rgb(damping_r, damping_g, damping_b);
		Body.type = XSolid_Crystal;
		Body.refractive_index = this->refractive_index;

		//Body.col = rgblight(White); Body.type = XSolid_Smooth;

		Body.init();

		return Body;
	}

	void push(World &W) {
		W.add(&Body);
	}
};


class GlassMan_std {
	GlassMan G;
	const static double length_of_shanks, length_of_thighs, length_of_upperarms, length_of_forearms,
		length_of_head, length_of_chest, length_of_waist, distance_between_buttocks, distance_between_shoulders;
	const static double width_of_forehead, width_of_chin, depth_of_forehead, radius_of_neck, length_of_neck,
		depth_of_chest, width_of_waist, lower_width_of_waist, lower_depth_of_waist, height_of_waist;
	const static double upper_radius_of_shanks, lower_radius_of_shanks, upper_radius_of_thighs, lower_radius_of_thighs,
		upper_radius_of_upperarms, lower_radius_of_upperarms, upper_radius_of_forearms, lower_radius_of_forearms;
	const static double width_of_heels, height_of_heels, width_of_tiptoes, length_of_feet;
	const static double head_rounding, neck_rounding, chest_rounding, waist_rounding,
		upperarms_rounding, forearms_rounding, thighs_rounding, shanks_rounding, feet_rounding;

public:
	GlassMan_std() {}
	GlassMan_std(const GlassMan_std &other) {
		this->Pos = other.Pos, this->top = other.top, this->dir = other.dir, this->v_head = other.v_head, this->v_neck = other.v_neck, this->v_chest = other.v_chest, this->v_waist = other.v_waist, this->v_upperarm_l = other.v_upperarm_l, this->v_upperarm_r = other.v_upperarm_r, this->v_forearm_l = other.v_forearm_l, this->v_forearm_r = other.v_forearm_r, this->v_thigh_l = other.v_thigh_l, this->v_thigh_r = other.v_thigh_r, this->v_shank_l = other.v_shank_l, this->v_shank_r = other.v_shank_r, this->v_foot_l = other.v_foot_l, this->v_foot_r = other.v_foot_r, this->v_head_side = other.v_head_side, this->v_chest_side = other.v_chest_side, this->v_waist_side = other.v_waist_side, this->v_foot_l_side = other.v_foot_l_side, this->v_foot_r_side = other.v_foot_r_side, this->auto_fit = other.auto_fit;
	}
	~GlassMan_std() {}

	point Pos; vec3 top, dir;
	vec3 v_head, v_neck, v_chest;	// fore vectors
	vec3 v_waist, v_upperarm_l, v_upperarm_r, v_forearm_l, v_forearm_r,
		v_thigh_l, v_thigh_r, v_shank_l, v_shank_r, v_foot_l, v_foot_r;
	vec3 v_head_side, v_chest_side, v_waist_side, v_foot_l_side, v_foot_r_side;		// to left side, beside foot

	bool auto_fit;

	static double height() {
		return height_of_heels + length_of_shanks + length_of_thighs + length_of_waist + length_of_chest + length_of_neck + length_of_neck;
	}

	point construct();

	void push(World &W) {
		this->construct();
		G.push(W);
	}
};

#pragma region GlassMan_Initialize
const double GlassMan_std::length_of_shanks = 0.45;
const double GlassMan_std::length_of_thighs = 0.30;
const double GlassMan_std::length_of_upperarms = 0.28;
const double GlassMan_std::length_of_forearms = 0.28;
const double GlassMan_std::length_of_head = 0.33;
const double GlassMan_std::length_of_chest = 0.33;
const double GlassMan_std::length_of_waist = 0.22;
const double GlassMan_std::distance_between_buttocks = 0.13;
const double GlassMan_std::distance_between_shoulders = 0.34;
const double GlassMan_std::width_of_forehead = 0.24;
const double GlassMan_std::width_of_chin = 0.18;
const double GlassMan_std::depth_of_forehead = 0.14;
const double GlassMan_std::radius_of_neck = 0.024;
const double GlassMan_std::length_of_neck = 0.03;
const double GlassMan_std::depth_of_chest = 0.14;
const double GlassMan_std::width_of_waist = 0.22;
const double GlassMan_std::lower_width_of_waist = 0.28;
const double GlassMan_std::lower_depth_of_waist = 0.14;
const double GlassMan_std::height_of_waist = 0.20;
const double GlassMan_std::upper_radius_of_shanks = 0.04;
const double GlassMan_std::lower_radius_of_shanks = 0.024;
const double GlassMan_std::upper_radius_of_thighs = 0.06;
const double GlassMan_std::lower_radius_of_thighs = 0.04;
const double GlassMan_std::upper_radius_of_upperarms = 0.04;
const double GlassMan_std::lower_radius_of_upperarms = 0.03;
const double GlassMan_std::upper_radius_of_forearms = 0.03;
const double GlassMan_std::lower_radius_of_forearms = 0.025;
const double GlassMan_std::width_of_heels = 0.1;
const double GlassMan_std::height_of_heels = 0.06;
const double GlassMan_std::width_of_tiptoes = 0.06;
const double GlassMan_std::length_of_feet = 0.13;
const double GlassMan_std::head_rounding = 0.0699;
const double GlassMan_std::neck_rounding = 0.01;
const double GlassMan_std::chest_rounding = 0.05;
const double GlassMan_std::waist_rounding = 0.06;
const double GlassMan_std::upperarms_rounding = 0.03;
const double GlassMan_std::forearms_rounding = 0.025;
const double GlassMan_std::thighs_rounding = 0.035;
const double GlassMan_std::shanks_rounding = 0.02;
const double GlassMan_std::feet_rounding = 0.015;
#pragma endregion

point GlassMan_std::construct() {
	top /= top.mod(), dir /= dir.mod();
	G.Pos = Pos, G.top = top, G.dir = dir;
	v_head /= v_head.mod(), v_neck /= v_neck.mod(), v_chest /= v_chest.mod(),
		v_waist /= v_waist.mod(), v_upperarm_l /= v_upperarm_l.mod(), v_upperarm_r /= v_upperarm_r.mod(), v_forearm_l /= v_forearm_l.mod(), v_forearm_r /= v_forearm_r.mod(),
		v_thigh_l /= v_thigh_l.mod(), v_thigh_r /= v_thigh_r.mod(), v_shank_l /= v_shank_l.mod(), v_shank_r /= v_shank_r.mod(), v_foot_l /= v_foot_l.mod(), v_foot_r /= v_foot_r.mod(),
		v_head_side /= v_head_side.mod(), v_chest_side /= v_chest_side.mod(), v_waist_side /= v_waist_side.mod(), v_foot_l_side /= v_foot_l_side.mod(), v_foot_r_side /= v_foot_r_side.mod();
	correct_vector(v_head, v_head_side), correct_vector(v_chest, v_chest_side), correct_vector(v_waist, v_waist_side), correct_vector(v_foot_l, v_foot_l_side), correct_vector(v_foot_r, v_foot_r_side);

	G.width_of_waist = width_of_waist, G.lower_width_of_waist = lower_width_of_waist, G.lower_depth_of_waist = lower_depth_of_waist, G.waist_rounding = waist_rounding;
	G.Waist = point(0, 0, 0);
	point P = G.Waist + height_of_waist * v_waist;
	G.Buttock_l = P + 0.5*distance_between_buttocks * v_waist_side, G.Buttock_r = G.Buttock_l - distance_between_buttocks * v_waist_side;
	G.Knee_l = G.Buttock_l + length_of_thighs * v_thigh_l, G.Knee_r = G.Buttock_r + length_of_thighs * v_thigh_r;
	G.Heel_l = G.Knee_l + length_of_shanks * v_shank_l, G.Heel_r = G.Knee_r + length_of_shanks * v_shank_r;
	G.upper_radius_of_thigh = upper_radius_of_thighs, G.lower_radius_of_thigh = lower_radius_of_thighs, G.thigh_rounding = thighs_rounding;
	G.upper_radius_of_shank = upper_radius_of_shanks, G.lower_radius_of_shank = lower_radius_of_shanks, G.shank_rounding = shanks_rounding;
	G.foot_dir_l = v_foot_l, G.foot_dir_r = v_foot_r, G.foot_side_l = v_foot_l_side, G.foot_side_r = v_foot_r_side, G.length_of_foot = length_of_feet, G.foot_rounding = feet_rounding;
	G.width_of_heel = width_of_heels, G.height_of_heel = height_of_heels, G.width_of_tiptoe = width_of_tiptoes;
	point tiptoe_l = G.Heel_l + (length_of_feet * feet_rounding) * v_foot_l, tiptoe_r = G.Heel_r + (length_of_feet * feet_rounding) * v_foot_r;

	G.depth_of_chest = depth_of_chest, G.chest_rounding = chest_rounding;
	P = G.Waist + length_of_chest * v_chest;
	G.Shoulder_l = P + 0.5*distance_between_shoulders * v_chest_side, G.Shoulder_r = G.Shoulder_l - distance_between_shoulders * v_chest_side;
	G.Elbow_l = G.Shoulder_l + length_of_upperarms * v_upperarm_l, G.Elbow_r = G.Shoulder_r + length_of_upperarms * v_upperarm_r;
	G.Hand_l = G.Elbow_l + length_of_forearms * v_forearm_l, G.Hand_r = G.Elbow_r + length_of_forearms * v_forearm_r;
	G.upper_radius_of_upperarm = upper_radius_of_upperarms, G.lower_radius_of_upperarm = lower_radius_of_upperarms, G.upperarm_rounding = upperarms_rounding;
	G.upper_radius_of_forearm = upper_radius_of_forearms, G.lower_radius_of_forearm = lower_radius_of_forearms, G.forearm_rounding = forearms_rounding;

	G.Head_lower = P + (length_of_neck + chest_rounding) * v_chest, G.Head_top = G.Head_lower + length_of_head * v_head; G.Head_side = v_head_side;
	G.width_of_forehead = width_of_forehead, G.width_of_chin = width_of_chin, G.depth_of_forehead = depth_of_forehead, G.head_rounding = head_rounding;
	G.neck_radius = radius_of_neck, G.neck_rounding = neck_rounding;

	if (auto_fit) {
		P = 0.5*(G.Buttock_l + G.Buttock_r);
		P.z = min(tiptoe_l.z, tiptoe_r.z);
		double ht = width_of_tiptoes / width_of_heels * height_of_heels;
		ht *= 0.5*(tiptoe_l.z < tiptoe_r.z ? sqrt(1 - v_foot_l.z*v_foot_l.z) : sqrt(1 - v_foot_r.z*v_foot_r.z));
		//P.z += ht;	// bug
		P = Pos - P;
	}
	else P = Pos;
	//G.Heel_l += P, G.Heel_r += P, G.Knee_l += P, G.Knee_r += P, G.Buttock_l += P, G.Buttock_r += P,
	//	G.Waist += P, G.Shoulder_l += P, G.Shoulder_r += P, G.Elbow_l += P, G.Elbow_r += P, G.Hand_l += P, G.Hand_r += P, G.Head_lower += P, G.Head_top += P, G.Head_side += P;
	G.Pos = P;

	//return G.construct();
	G.construct(); return G.Waist + P;
}


XSolid GoldenCoin_Constructor(point Pos, double rz) {
	const double R = 0.4, r = 0.35, d = 0.05, dc = 0.01;
	point P = 0.5*point(d*cos(rz), d*sin(rz));
	XObjs::Cylinder M(Pos + P, Pos - P, R);

	XObjs::Cylinder M_(Pos - P, Pos - (1 - (dc / d))*P, r);
	XSolid X = CSG_SubtractionOp(XSolid(M), XSolid(M_));
	//X = M_;

	Figure2D FBase; {
		point2D Dr(R, R);
		point2D A(0.42692, 0.63541), A_(0.55867, 0.6258), B(0.59236, 0.54138), C(0.48977, 0.53258), C_(0.47414, 0.57165),
			D(0.42682, 0.58243), E(0.42688, 0.44652), E_(0.61255, 0.4178), F(0.61291, 0.32093), F_(0.61265, 0.20646),
			G(0.42677, 0.1839), H(0.42744, 0.12651), I(0.37215, 0.12664), J(0.37213, 0.18316), J_(0.21878, 0.19258),
			K(0.1877, 0.30306), L(0.29337, 0.31029), L_(0.31275, 0.25804), M(0.37233, 0.24025), N(0.42669, 0.23671),
			N_(0.51266, 0.24762), O(0.51281, 0.31097), O_(0.513, 0.3594), P(0.42674, 0.37457), Q(0.37189, 0.38538),
			Q_(0.2, 0.42), R(0.20327, 0.5187), R_(0.20661, 0.61941), S(0.37225, 0.63516), T(0.37237, 0.58324),
			U(0.37237, 0.45838), U_(0.30528, 0.47743), V(0.30466, 0.51988), V_(0.30425, 0.56594), W(0.37215, 0.66494), X(0.42691, 0.66522);
		A -= Dr, A_ -= Dr, B -= Dr, C -= Dr, C_ -= Dr, D -= Dr, E -= Dr, E_ -= Dr, F -= Dr, F_ -= Dr, G -= Dr, H -= Dr, I -= Dr, J -= Dr, J_ -= Dr, K -= Dr, L -= Dr, L_ -= Dr,
			M -= Dr, N -= Dr, N_ -= Dr, O -= Dr, O_ -= Dr, P -= Dr, Q -= Dr, Q_ -= Dr, R -= Dr, R_ -= Dr, S -= Dr, T -= Dr, U -= Dr, U_ -= Dr, V -= Dr, V_ -= Dr, W -= Dr, X -= Dr;

		QuadraticBezier2D c_01(A, A_, B), c_03(C, C_, D), c_05(E, E_, F), c_06(F, F_, G), c_10(J, J_, K),
			c_12(L, L_, M), c_13(N, N_, O), c_14(O, O_, P), c_17(Q, Q_, R), c_18(R, R_, S), c_20(U, U_, V), c_21(V, V_, T);
		Segment2D c_02(B, C), c_04(D, E), c_07(G, H), c_08(H, I), c_09(I, J), c_11(K, L), c_15(P, N), c_16(M, Q), c_19(T, U), c_22(S, W), c_23(W, X), c_24(X, A);

		FBase.assign({ &c_01, &c_02, &c_03, &c_04, &c_05, &c_06, &c_07, &c_08, &c_09, &c_10, &c_11, &c_12, &c_13, &c_14, &c_15, &c_16, &c_17, &c_18, &c_19, &c_20, &c_21, &c_22, &c_23, &c_24 });
	}

	XSolid Pt(XObjs::Extrusion_xOy(FBase, -0.5*d, 0.5*d));
	Pt = CSG_Rotation(Pt, PI / 2, 0, rz - PI / 2);
	Pt = CSG_Translation(Pt, Pos);

	//ScanXSolid(Pt, 10, 10, 10, 120000);

	X = CSG_UnionOp(X, Pt);
	//X = Pt;

	X.type = XSolid_Smooth; X.setColor(Gold);
	X.type = XSolid_LightSource; X.setColor(GoldenRod);
	return X;
}
