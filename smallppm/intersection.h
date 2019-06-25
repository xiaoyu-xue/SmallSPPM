#pragma once

#include "linagl.h"
#include "ray.h"

struct Intersection {
	Intersection() {}
	Vec hit, n, nl, wo;
	Ray SpawnTo(const Intersection &isect) const {
		Vec dir = isect.hit - (hit + nl * nEps);
		real d = dir.length();
		dir.normalize();
		return Ray(hit + dir * rayeps + nl * rayeps, dir, d - shadowRayEps);
	}
};

