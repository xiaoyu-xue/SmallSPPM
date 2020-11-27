#pragma once

#include "visual/Medium.h"

NAMESPACE_BEGIN

//class MediumIntersection;

class HomogeneousMedium : public Medium {
public:
	HomogeneousMedium(const Vec3 &sigmaA, const Vec3 &sigmaS, real g):
		sigma_a(sigmaA), sigma_s(sigmaS), sigma_t(sigma_s + sigma_a), g(g) {}

	Vec3 Tr(const Ray& ray, StateSequence& rand) const;

	Vec3 Sample(const Ray& ray, StateSequence& rand, MemoryPool& arena, MediumIntersection* mi) const;

	Vec3 EquiAngularSampling(const Ray& ray, StateSequence& rand, MemoryPool& arena,
		const Intersection& lightPoint, MediumIntersection* mi) const override;

private:
	Vec3 sigma_a, sigma_s, sigma_t;
	real g;
};

NAMESPACE_END
