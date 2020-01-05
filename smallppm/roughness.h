#pragma once

#include "def.h"
#include "material.h"
#include "texture.h"
#include "bsdf.h"
#include "microfacet.h"
NAMESPACE_BEGIN

class RoughnessMaterial : public Material {
public:
	RoughnessMaterial(const std::shared_ptr<Texture>& reflectence, real alpha) : reflectence(reflectence), alpha(alpha) {

	}
	void ComputeScatteringFunction(Intersection* isect,
		TransportMode mode = TransportMode::Radiance) const {
		MicrofacetDistribution *distribution = new GGXDistribution(alpha);
		isect->bsdf = std::shared_ptr<BSDF>(new MicrofacetReflectionBSDF(*isect, distribution, reflectence->Sample(*isect)));
	}
private:
	std::shared_ptr<Texture> reflectence;
	real alpha;
};

NAMESPACE_END