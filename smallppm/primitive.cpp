#include "primitive.h"

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