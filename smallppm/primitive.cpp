#include "primitive.h"

NAMESPACE_BEGIN

bool GeometryPrimitive::Intersect(const Ray& r, Intersection* isect, real* t) const {

	return shape->Intersect(r, isect, t);
}

bool GeometryPrimitive::Intersect(const Ray& r) const {
	return shape->Intersect(r);
}

void GeometryPrimitive::ComputeScatteringFunction(Intersection* isect,
	TransportMode mode) {
	if (material) {
		material->ComputeScatteringFunction(isect, mode);
	}
}


NAMESPACE_END