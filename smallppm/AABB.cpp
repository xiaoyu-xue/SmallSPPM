#include "AABB.h"

AABB AABB::Union(const AABB& a, const AABB& b) {

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