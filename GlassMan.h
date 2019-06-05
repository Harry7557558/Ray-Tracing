#pragma once

#include "World.h"

class GlassMan {
	// Control Points
	point Head, Head_f, Head_le, Head_re, Head_ch, Neck; 	// center, forehead, left eye, right eye, chin
	point Chest_u, Chest_m, Chest_l, Waist, Belly;	// upper, median, lower
	point Shoulder_l, shoulder_r, Elbow_l, Elbow_r, Hand_l, Hand_r;	// left and right
	point /*Butt_l, Butt_r, Knee_l, Knee_r, Ankle_l, Ankle_r, */Foot_l, Foot_r;

	XSolid Body;

	XObjs::Box_affine _HEAD, _CHEST, _WAIST, _BELLY;
	XObjs::Cylinder _NECK;
	XObjs::Box_affine _SHOULDER_L, _SHOULDER_R;
	XObjs::Cone_capped _UPPER_ARM_L, _UPPER_ARM_R, _FOREARM_L, _FOREARM_R;
	XObjs::Sphere _ELBOW_L, _ELBOW_R, _HAND_L, _HAND_R;
	XObjs::Box_affine _BUTTOCK_L, _BUTTOCK_R;
	XObjs::Box_affine _THIGH_L, _THIGH_R;
	XObjs::Box_affine _KNEE_L, _KNEE_R;
	XObjs::Cone_capped _SHANK_L, _SHANK_R;
	XObjs::Box_affine _FOOT_L, _FOOT_R;

public:
	double damping_r, damping_g, damping_b;
	double refractive_index;

	/* Fill all of these parameters legally before constructing */
	point Pos; 	// position the glass man standing/sitting at
	//vec3 top;	// direction the glass man's body towards (up)
	//vec3 dir;	// direction the glass man facing

	// The following points/vectors are relative to Pos (default 0,0,0), top (default 0,0,1) and dir (default 1,0,0)
	point Heel_l, Heel_r;		// position of heel/ankle
	vec3 foot_dir_l, foot_dir_r;	// direction the foot pointing
	double width_of_heel, height_of_heel, width_of_tiptoe, length_of_foot, foot_rounding;	// width and height of heels and tiptoes; length of feet; rounding radius
	point Knee_l, Knee_r;	// position of knees
	double upper_radius_of_shank, lower_radius_of_shank, shank_rounding;		// thickness of upper and lower end of shanks; rounding radius
	point Buttock_l, Buttock_r;
	double upper_radius_of_thigh, lower_radius_of_thigh, thigh_rounding;


	GlassMan() {
		damping_r = damping_g = damping_b = 0;
		refractive_index = 1.5;
		Pos = point(0, 0, 0)/*, top = point(0, 0, 1), dir = point(1, 0, 0)*/;
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
		//top /= top.mod(), dir /= dir.mod();
		//dir -= dot(dir, top) * top; dir /= dir.mod();	// now dir is perpendicular to top
		//vec3 js = cross(dir, top); js /= js.mod();
		//double rz = atan2(dir.y, dir.x), ry = atan2(-dir.z, hypot(dir.x, dir.y)), rx = atan2(js.z, top.z);
		/*if (dir.y == 0 && dir.x == 0) {
			double a = cos(rx)*sin(ry), b = sin(rx)*sin(ry), c = top.x;
			rz = 2 * atan((b + sqrt(a*a + b * b - c * c)) / (a - c));
		}*/

		Body = XSolid();

		// constructing feet
		XSolid tmp_FEET; {
			_FOOT_L = XObjs::Box_affine();
			foot_dir_l /= foot_dir_l.mod(), foot_dir_r /= foot_dir_r.mod();
			_FOOT_L.translate(0, -0.5, 0);
			double foot_h_w = width_of_heel - 2 * foot_rounding, foot_h_h = height_of_heel - 2 * foot_rounding,
				foot_t_w = width_of_tiptoe - 2 * foot_rounding, foot_l = length_of_foot - 2 * foot_rounding;
			_FOOT_L.scale(foot_l * foot_h_w / foot_t_w, foot_h_w, foot_h_h);
			_FOOT_L.perspective(foot_h_w / foot_t_w - 1, 0, 0);
			_FOOT_R = _FOOT_L;
			double foot_ry = atan2(-foot_dir_l.z, hypot(foot_dir_l.x, foot_dir_l.y)), foot_rz = atan2(foot_dir_l.y, foot_dir_l.x);
			_FOOT_L.rotate(0, foot_ry, foot_rz);
			foot_ry = atan2(-foot_dir_r.z, hypot(foot_dir_r.x, foot_dir_r.y)), foot_rz = atan2(foot_dir_r.y, foot_dir_r.x);
			_FOOT_R.rotate(0, foot_ry, foot_rz);
			_FOOT_L += Heel_l, _FOOT_R += Heel_r;
			XSolid tmp_FOOT_L = _FOOT_L, tmp_FOOT_R = _FOOT_R;
			tmp_FOOT_L = CSG_RoundingOp(tmp_FOOT_L, foot_rounding), tmp_FOOT_R = CSG_RoundingOp(tmp_FOOT_R, foot_rounding);
			tmp_FEET = CSG_UnionOp(tmp_FOOT_L, tmp_FOOT_R);
		}

		// constructing shanks
		XSolid tmp_SHANKS; {
			double shank_up_r = upper_radius_of_shank - shank_rounding, shank_lw_r = lower_radius_of_shank - shank_rounding;
			point shank_root = Heel_l + shank_lw_r * foot_dir_l, shank_dir = shank_root - Knee_l; shank_dir /= shank_dir.mod();
			_SHANK_L = XObjs::Cone_capped(Knee_l + shank_rounding * shank_dir, shank_root - (shank_rounding + 0.3 * height_of_heel) * shank_dir, shank_up_r, shank_lw_r);
			shank_root = Heel_r + shank_lw_r * foot_dir_r; shank_dir = shank_root - Knee_r, shank_dir /= shank_dir.mod();
			_SHANK_R = XObjs::Cone_capped(Knee_r + shank_rounding * shank_dir, shank_root - (shank_rounding + 0.3 * height_of_heel) * shank_dir, shank_up_r, shank_lw_r);
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
		



		Body = CSG_UnionOp(tmp_FEET, tmp_SHANKS);
		Body = CSG_UnionOp(Body, tmp_THIGHS);

		XObjs::Cone Arrow(point(0, 0, 1), point(2, 0, 1), 0.5);
		//Body = Arrow;

		//Body = CSG_Rotation(Body, rx, ry, rz);
		Body = CSG_Translation(Body, Pos);

		Body.col = rgb(damping_r, damping_g, damping_b);
		Body.type = XSolid_Crystal;
		Body.refractive_index = this->refractive_index;

		Body.init();

		return Body;
	}

	void push(World &W) {
		W.add(&Body);
	}
};


