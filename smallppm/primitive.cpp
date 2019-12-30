#include "primitive.h"

NAMESPACE_BEGIN

bool GeometryPrimitive::Intersect(const Ray& r, Intersection* isect, real* t) const {

	if (shape->Intersect(r, isect, t)) {
		r.tMax = *t;
		isect->primitive = this;
		return true;
	}
	return false;
}

bool GeometryPrimitive::Intersect(const Ray& r) const {
	return shape->Intersect(r);
}

void GeometryPrimitive::ComputeScatteringFunction(Intersection* isect,
	TransportMode mode) const {
	if (material) {
		material->ComputeScatteringFunction(isect, mode);
	}
}


NAMESPACE_END