#pragma once

#include "linagl.h"
#include "utils.h"

NAMESPACE_BEGIN

class Medium;
struct Ray {
	Ray() : medium(nullptr) {
		tMax = Inf;
		tMin = 0.f;
	};
	Ray(Vec3 o_, Vec3 d_, real tmax_ = Inf, real tmin_ = 0.f, const Medium* medium = nullptr) : 
		o(o_ + d_ * rayeps), d(d_), tMax(tmax_), tMin(tmin_), medium(medium) {}
	Vec3 operator()(real t) const {
		return o + d * t;
	}
	Vec3 o, d;
	mutable real tMin, tMax;
	const Medium* medium;
};

std::ostream& operator<<(std::ostream &os, const Ray &ray);

NAMESPACE_END