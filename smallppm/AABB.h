#pragma once

#include "linagl.h"
#include <algorithm>

NAMESPACE_BEGIN

struct AABB {
	Vec3 minPoint, maxPoint; // axis aligned bounding box
	inline void Fit(const Vec3 &p)
	{
		if (p.x < minPoint.x) minPoint.x = p.x; // min
		if (p.y < minPoint.y) minPoint.y = p.y; // min
		if (p.z < minPoint.z) minPoint.z = p.z; // min
		maxPoint.x = std::max(p.x, maxPoint.x);
		maxPoint.y = std::max(p.y, maxPoint.y);
		maxPoint.z = std::max(p.z, maxPoint.z);
	}
	inline void Reset() {
		minPoint = Vec3(1e20, 1e20, 1e20);
		maxPoint = Vec3(-1e20, -1e20, -1e20);
	}
};

NAMESPACE_END