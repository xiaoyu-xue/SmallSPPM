#include "Primitive.h"

NAMESPACE_BEGIN

bool Primitive::Intersect(const Ray& r, Intersection* isect) const {
	real t;
	if (mpShape->Intersect(r, isect, &t)) {
		r.m_tMax = t;
		isect->primitive = this;
		isect->primId = this->mPrimId;
		if (mMediumInterface.IsMediumTransition())
			isect->mediumInterface = mMediumInterface;
		else
			isect->mediumInterface = MediumInterface(r.mpMedium);
		return true;
	}
	return false;
}

bool Primitive::Intersect(const Ray& r) const {
	return mpShape->Intersect(r);
}

void Primitive::ComputeScatteringFunction(Intersection* isect, MemoryPool& arena,
	TransportMode mode) const {
	if (mpMaterial) {
		mpMaterial->ComputeScatteringFunction(isect, arena, mode);
	}
}

void Primitive::QueryIntersectionInfo(const Ray& ray, Intersection* isect) const {
	mpShape->QueryIntersectionInfo(ray, isect);
}

NAMESPACE_END