#pragma once

#include "linagl.h"
#include "ray.h"

NAMESPACE_BEGIN

struct Intersection {
	Intersection() { rayEps = 1e-3f; }
	Vec3 hit, n, nl, wo;
	real rayEps;
	Ray SpawnTo(const Intersection &isect) const {
		Vec3 dir = isect.hit - (hit + nl * nEps);
		real d = dir.Length();
		dir.Normalize();
		//return Ray(hit + dir * rayeps + nl * rayeps, dir, d - shadowRayEps);
		return Ray(hit, dir, d * (1.f - isect.rayEps), rayEps);
	}
};

NAMESPACE_END