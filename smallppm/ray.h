#pragma once

#include "linagl.h"

NAMESPACE_BEGIN

struct Ray {
	Ray() {
		tMax = Inf;
	};
	Ray(Vec3 o_, Vec3 d_, real tmax_ = Inf) : o(o_), d(d_), tMax(tmax_) {}
	Vec3 o, d;
	real tMax;
};

NAMESPACE_END