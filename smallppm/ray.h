#pragma once

#include "linagl.h"

struct Ray {
	Ray() {
		tMax = Inf;
	};
	Ray(Vec o_, Vec d_, real tmax_ = Inf) : o(o_), d(d_), tMax(tmax_) {}
	Vec o, d;
	real tMax;
};
