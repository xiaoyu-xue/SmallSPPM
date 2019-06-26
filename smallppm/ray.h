#pragma once

#include "linagl.h"

NAMESPACE_BEGIN

struct Ray {
	Ray() {
		tMax = Inf;
	};
	Ray(Vec o_, Vec d_, real tmax_ = Inf) : o(o_), d(d_), tMax(tmax_) {}
	Vec o, d;
	real tMax;
};

NAMESPACE_END