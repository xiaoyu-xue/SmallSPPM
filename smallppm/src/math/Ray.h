#pragma once

#include "Linagl.h"
#include "common/Core.h"

NAMESPACE_BEGIN

class Medium;
struct Ray {
	Ray() 
		: mpMedium(nullptr)
	{
		tMax = Inf;
		tMin = 0.f;
	}

	Ray(Vec3 orig, Vec3 dir, real tmax = Inf, real tmin = 0.f, const Medium* medium = nullptr) 
		: mOrig(orig + dir * RayEps), mDir(dir), tMax(tmax), tMin(tmin), mpMedium(medium) 
	{}

	Vec3 operator()(real t) const 
	{
		return mOrig + mDir * t;
	}

	Vec3 mOrig, mDir;
	mutable real tMin, tMax;
	const Medium* mpMedium;
};

std::ostream& operator<<(std::ostream &os, const Ray &ray);

NAMESPACE_END