#include "primitive.h"

NAMESPACE_BEGIN

bool Primitive::Intersect(const Ray& r, Intersection* isect) const {
	real t;
	if (shape->Intersect(r, isect, &t)) {
		r.tMax = t;
		isect->primitive = this;
		isect->primId = this->primId;
		if (mediumInterface.IsMediumTransition())
			isect->mediumInterface = mediumInterface;
		else
			isect->mediumInterface = MediumInterface(r.medium);
		return true;
	}
	return false;
}

bool Primitive::Intersect(const Ray& r) const {
	return shape->Intersect(r);
}

void Primitive::ComputeScatteringFunction(Intersection* isect, MemoryArena& arena,
	TransportMode mode) const {
	if (material) {
		material->ComputeScatteringFunction(isect, arena, mode);
	}
}

void Primitive::QueryIntersectionInfo(const Ray& ray, Intersection* isect) const {
	shape->QueryIntersectionInfo(ray, isect);
}

NAMESPACE_END