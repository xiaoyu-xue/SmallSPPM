#pragma once

#include "linagl.h"
#include <algorithm>
#include "utils.h"

NAMESPACE_BEGIN

FORCE_INLINE void CoordinateSystem(const Vec3 &v1, Vec3 *v2, Vec3 *v3) {
	if (std::abs(v1.x) > std::abs(v1.y))
		*v2 = Vec3(-v1.z, 0, v1.x) / std::sqrt(v1.x * v1.x + v1.z * v1.z);
	else
		*v2 = Vec3(0, v1.z, -v1.y) / std::sqrt(v1.y * v1.y + v1.z * v1.z);
	*v3 = v1.Cross(*v2);
}


Vec3 Reflect(const Vec3 &inDir, const Vec3 &n) {
	return 2.f * inDir.Dot(n) * n - inDir;
}

NAMESPACE_END