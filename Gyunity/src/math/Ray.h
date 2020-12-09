#pragma once

#include "Linagl.h"
#include "common/Core.h"

GYT_NAMESPACE_BEGIN

class Medium;
struct Ray {
	Ray() 
		: mpMedium(nullptr)
	{
		m_tMax = Inf;
		m_tMin = 0.f;
	}

	Ray(Vec3 orig, Vec3 dir, real tmax = Inf, real tmin = 0.f, const Medium* medium = nullptr) 
		: mOrig(orig + dir * RayEps), mDir(dir), m_tMax(tmax), m_tMin(tmin), mpMedium(medium) 
	{}

	Vec3 operator()(real t) const 
	{
		return mOrig + mDir * t;
	}

	Vec3 mOrig, mDir;
	mutable real m_tMin, m_tMax;
	const Medium* mpMedium;
};

std::ostream& operator<<(std::ostream &os, const Ray &ray);

GYT_NAMESPACE_END