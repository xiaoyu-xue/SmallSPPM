#pragma once

#include "linagl.h"
#include "utils.h"

NAMESPACE_BEGIN

struct Ray {
	Ray() {
		tMax = Inf;
	};
	Ray(Vec3 o_, Vec3 d_, real tmax_ = Inf, real tmin_ = 0.f) : 
		o(o_ + d_ * rayeps), d(d_), tMax(tmax_), tMin(tmin_) {}
	Vec3 operator()(real t) const {
		return o + d * t;
	}
	Vec3 o, d;
	mutable real tMin, tMax;
};

std::ostream& operator<<(std::ostream &os, const Ray &ray);

NAMESPACE_END