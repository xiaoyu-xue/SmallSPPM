#pragma once

#include "linagl.h"
#include "utils.h"

NAMESPACE_BEGIN

struct Ray {
	Ray() {
		tMax = Inf;
	};
	Ray(Vec3 o_, Vec3 d_, real tmax_ = Inf) : o(o_ + d_ * rayeps), d(d_), tMax(tmax_ - rayeps) {}
	Vec3 o, d;
	real tMax;
};

NAMESPACE_END