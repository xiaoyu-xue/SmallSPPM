#pragma once

#include "linagl.h"
#include "ray.h"
#include "geometry_util.h"

NAMESPACE_BEGIN

struct Intersection {
	Intersection() { rayEps = 1e-3f; }
	Intersection(const Vec3 &hit, const Vec3 &n, const Vec3 &nl, const Vec3 &wo, const Vec3 &pError):
		hit(hit), n(n), nl(nl), wo(wo), pError(pError){}
	Vec3 hit, n, nl, wo;
	real rayEps;
	Vec3 pError;
	//Ray SpawnTo(const Intersection &isect) const {
	//	Vec3 dir = isect.hit - (hit + nl * nEps);
	//	real d = dir.Length();
	//	dir.Normalize();
	//	//return Ray(hit + dir * rayeps + nl * rayeps, dir, d - shadowRayEps);
	//	return Ray(hit, dir, d * (1.f - isect.rayEps), rayEps);
	//}

	Ray SpawnTo(const Intersection &it) const {
		Vec3 origin = OffsetRayOrigin(hit, pError, n, it.hit - hit);
		Vec3 target = OffsetRayOrigin(it.hit, it.pError, it.n, origin - it.hit);
		Vec3 d = target - origin;
		return Ray(origin, d, 1 - shadowRayEps, 0.f);
	}

	Ray SpawnRay(const Vec3 &d) const {
		Vec3 o = OffsetRayOrigin(hit, pError, n, d);
		return Ray(o, d, Infinity, 0.f);
	}

};

NAMESPACE_END