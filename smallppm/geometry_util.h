#pragma once

#include "linagl.h"
#include <algorithm>
#include "utils.h"

FORCE_INLINE void CoordinateSystem(const Vec &v1, Vec *v2, Vec *v3) {
	if (std::abs(v1.x) > std::abs(v1.y))
		*v2 = Vec(-v1.z, 0, v1.x) / std::sqrt(v1.x * v1.x + v1.z * v1.z);
	else
		*v2 = Vec(0, v1.z, -v1.y) / std::sqrt(v1.y * v1.y + v1.z * v1.z);
	*v3 = v1 % (*v2);
}


Vec Reflect(const Vec &inDir, const Vec &n) {
	return 2.0 * inDir.dot(n) * n - inDir;
}