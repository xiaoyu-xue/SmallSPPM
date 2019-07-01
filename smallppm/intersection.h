#pragma once

#include "linagl.h"
#include "ray.h"

NAMESPACE_BEGIN

struct Intersection {
	Intersection() {}
	Vec3 hit, n, nl, wo;
	Ray SpawnTo(const Intersection &isect) const {
		Vec3 dir = isect.hit - (hit + nl * nEps);
		real d = dir.Length();
		dir.Normalize();
		return Ray(hit + dir * rayeps + nl * rayeps, dir, d - shadowRayEps);
	}
};

NAMESPACE_END