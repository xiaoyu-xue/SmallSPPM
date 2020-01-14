#include "primitive.h"

NAMESPACE_BEGIN

bool GeometryPrimitive::Intersect(const Ray& r, Intersection* isect, real* t) const {

	if (shape->Intersect(r, isect, t)) {
		r.tMax = *t;
		isect->primitive = this;
		isect->shapeId = this->GetShape()->shapeId;
		//DEBUG_PIXEL_IF() {
		//	std::cout << shape->shapeId << ": " << "t " << *t << " ray.t: " << r.tMax << std::endl;
		//}
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