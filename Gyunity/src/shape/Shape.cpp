#include "Shape.h"
#include "core/Intersection.h"

GYT_NAMESPACE_BEGIN

Intersection Shape::Sample(const Intersection& isect, real* pdf, const Vec2& u) const {
	Intersection it = Sample(pdf, u);
	Vec3 wi = isect.mPos - it.mPos;
	if (wi.Length() == 0) {
		*pdf = 0;
	}
	else {
		wi.Normalize();
		*pdf *= Distance2(isect.mPos, it.mPos) / std::abs(Dot(it.mNormal, -wi));
		if (std::isinf(*pdf)) *pdf = 0;
	}
	it.mShapeId = mShapeId;
	return it;
}

GYT_NAMESPACE_END