#pragma once

#include "linagl.h"
#include "ray.h"
#include "geometry_util.h"
#include "material.h"

NAMESPACE_BEGIN

class BSDF;
class Primitive;

class Intersection {
public:
	Vec3 hit, n, nl, wo;
	Vec3 dpdu, dpdv;
	Vec2 uv;
	real rayEps;
	Vec3 pError;

	const Primitive *primitive;
	std::shared_ptr<BSDF> bsdf;

	Intersection() { rayEps = 1e-3f; }

	Intersection(const Vec3 &hit, const Vec3 &n, const Vec3 &nl, const Vec3 &wo, const Vec3 &pError):
		hit(hit), n(n), nl(nl), wo(wo), pError(pError), rayEps(1e-3), primitive(nullptr){
		bsdf = nullptr;
	}

	Intersection(const Vec3& hit, const Vec3& n, const Vec3& nl, 
		const Vec3& dpdu, const Vec3& dpdv, const Vec3& wo, const Vec3& pError) :
		hit(hit), dpdu(dpdu), dpdv(dpdv), n(n), nl(nl), wo(wo), pError(pError), rayEps(1e-3), primitive(nullptr) {
		bsdf = nullptr;
	}

	void ComputeScatteringFunction(TransportMode mode = TransportMode::Radiance);

	Ray SpawnTo(const Intersection &it) const {
		Vec3 origin = OffsetRayOrigin(hit, pError, nl, it.hit - hit);
		Vec3 target = OffsetRayOrigin(it.hit, it.pError, it.nl, origin - it.hit);
		Vec3 d = target - origin;
		return Ray(origin, d, 1 - shadowRayEps, 0.f);
	}

	Ray SpawnRay(const Vec3 &d) const {
		Vec3 o = OffsetRayOrigin(hit, pError, nl, d);
		return Ray(o, d, Infinity, 0.f);
	}

	//Ray SpawnTo(const Intersection &isect) const {
	//	Vec3 dir = isect.hit - (hit + nl * nEps);
	//	real d = dir.Length();
	//	dir.Normalize();
	//	//return Ray(hit + dir * rayeps + nl * rayeps, dir, d - shadowRayEps);
	//	return Ray(hit, dir, d * (1.f - isect.rayEps), rayEps);
	//}
};

NAMESPACE_END