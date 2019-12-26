#include "AABB.h"

AABB Union(const AABB& a, const AABB& b) {

	Vec3 minPoint(
		std::min(a.minPoint.x, b.minPoint.x),
		std::min(a.minPoint.y, b.minPoint.y),
		std::min(a.minPoint.z, b.minPoint.z)
	);

	Vec3 maxPoint(
		std::max(a.maxPoint.x, b.maxPoint.x),
		std::max(a.maxPoint.y, b.maxPoint.y),
		std::max(a.maxPoint.z, b.maxPoint.z)
	);

	return AABB(minPoint, maxPoint);
}

AABB Union(const AABB& a, const Vec3& p) {
	AABB ret = a;
	ret.minPoint = Vec3(
		std::min(ret.minPoint.x, p.x),
		std::min(ret.minPoint.y, p.y),
		std::min(ret.minPoint.z, p.z)
	);
	ret.maxPoint = Vec3(
		std::max(ret.maxPoint.x, p.x),
		std::max(ret.maxPoint.y, p.y),
		std::max(ret.maxPoint.z, p.z)
	);
	return ret;
}